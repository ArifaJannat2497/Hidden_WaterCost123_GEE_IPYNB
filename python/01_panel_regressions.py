import pandas as pd, numpy as np, statsmodels.formula.api as smf, pymannkendall as mk, itertools, json
from linearmodels.panel import PanelOLS
d=pd.read_csv('../data/Boro_district_season_panel_2001_2025.csv').sort_values(['District','BoroYear']).reset_index(drop=True)
d['ET']=d.BoroET_ERA5Land; d['Rain']=d.BoroRain_mm; d['LST']=d.Mean_LST-273.15
d['PET']=d.BoroPET_ERA5Land; d['ETPET']=d.ET_PET_ratio; d['BI']=d.Boro_intensity_pct
d['WPress']=(d.ET-d.Rain)/d.ET
d['NIRv']=d.Max_NIRv; d['iNDVI']=d.Integrated_NDVI
d['WP']=d.NIRv/d.ET*1e4   # NIRv per 10,000 mm  -> report as x10^-4 mm^-1
d['WPi']=d.iNDVI/d.ET*1e3
d['post2013']=(d.BoroYear>=2013).astype(int)
d.to_csv('analysis_panel.csv',index=False)
out={}
def cen(df,cols):
    for c in cols: df[c+'_c']=df[c]-df[c].mean()
    return df
d=cen(d,['BI','LST','Rain','WPress','ETPET'])
# ---------- wild cluster bootstrap (Webb 6-point, full enumeration) ----------
W6=np.array([-np.sqrt(1.5),-1,-np.sqrt(.5),np.sqrt(.5),1,np.sqrt(1.5)])
def wcb(formula,df,var,cluster='District'):
    import statsmodels.api as sm, patsy
    y,X=patsy.dmatrices(formula,df,return_type='dataframe')
    cl=df.loc[y.index,cluster].values; G=np.unique(cl)
    j=list(X.columns).index(var)
    def tstat(yv):
        b=np.linalg.lstsq(X.values,yv,rcond=None)[0]; e=yv-X.values@b
        XtXi=np.linalg.pinv(X.values.T@X.values); meat=np.zeros((X.shape[1],)*2)
        for g in G:
            s=X.values[cl==g].T@e[cl==g]; meat+=np.outer(s,s)
        n,k=X.shape; c=len(G)/(len(G)-1)*(n-1)/(n-k)
        V=c*XtXi@meat@XtXi; return b[j]/np.sqrt(V[j,j])
    t0=tstat(y.values.ravel())
    Xr=X.drop(columns=var).values; br=np.linalg.lstsq(Xr,y.values.ravel(),rcond=None)[0]
    fit=Xr@br; er=y.values.ravel()-fit
    ts=[]
    for w in itertools.product(W6,repeat=len(G)):
        wm=dict(zip(G,w)); ystar=fit+er*np.array([wm[c] for c in cl]); ts.append(tstat(ystar))
    ts=np.array(ts); return float(np.mean(np.abs(ts)>=abs(t0)))
def dk(formula_rhs,dep,df,entity_fe=True,time_fe=False):
    p=df.set_index(['District','BoroYear'])
    rhs=formula_rhs+(' + EntityEffects' if entity_fe else '')+(' + TimeEffects' if time_fe else '')
    m=PanelOLS.from_formula(f'{dep} ~ 1 + {rhs}',p).fit(cov_type='kernel',kernel='bartlett',bandwidth=3)
    return m
rows=[]
def record(model_name,spec,var,coef,se,p,extra=None):
    rows.append(dict(model=model_name,spec=spec,term=var,coef=coef,se=se,p=p,**(extra or {})))
# ---------------- LEVEL MODEL ----------------
X='BI_c + LST_c + Rain_c + WPress_c'
specs={'A_districtFE_HC1':f'WP ~ {X} + C(District)','B_twowayFE_HC1':f'WP ~ {X} + C(District) + C(BoroYear)',
       'C_districtFE_ETPET':'WP ~ BI_c + LST_c + Rain_c + ETPET_c + C(District)','D_twowayFE_ETPET':'WP ~ BI_c + LST_c + Rain_c + ETPET_c + C(District) + C(BoroYear)',
       'E_districtFE_trend':f'WP ~ {X} + C(District) + C(District):BoroYear','F_iNDVI_twowayFE':f'WPi ~ {X} + C(District) + C(BoroYear)'}
lev={}
for k,f in specs.items():
    m=smf.ols(f,d).fit(cov_type='HC1'); lev[k]=m
    for v in ['BI_c','LST_c','Rain_c','WPress_c','ETPET_c']:
        if v in m.params: record('level',k,v,m.params[v],m.bse[v],m.pvalues[v],dict(n=int(m.nobs),r2=m.rsquared))
# DK + wild bootstrap for key
for k,(rhs,te) in {'A_districtFE_DK':(X,False),'B_twowayFE_DK':(X,True)}.items():
    m=dk(rhs,'WP',d,True,te)
    for v in ['BI_c','WPress_c']: record('level',k,v,m.params[v],m.std_errors[v],m.pvalues[v],dict(n=int(m.nobs)))
for k,f in [('A_districtFE_WCB',specs['A_districtFE_HC1']),('B_twowayFE_WCB',specs['B_twowayFE_HC1'])]:
    for v in ['BI_c','WPress_c']:
        m=lev[k[0]+k[1:].replace('_WCB','_HC1')] if False else None
        base=lev['A_districtFE_HC1' if k.startswith('A') else 'B_twowayFE_HC1']
        record('level',k,v,base.params[v],np.nan,wcb(f,d,v))
# leave-one-district-out (two-way FE)
for dist in d.District.unique():
    m=smf.ols(specs['B_twowayFE_HC1'],d[d.District!=dist]).fit(cov_type='HC1')
    record('level','LODO_drop_'+dist,'BI_c',m.params['BI_c'],m.bse['BI_c'],m.pvalues['BI_c'],dict(n=int(m.nobs)))
# ---------------- CHANGE MODEL ----------------
g=d.groupby('District')
for c in ['WP','WPi','ET','LST','Rain','NIRv','iNDVI']: d['d_'+c]=g[c].diff()
d['dln_NIRv']=g['NIRv'].transform(lambda s: np.log(s).diff()); d['dln_ET']=g['ET'].transform(lambda s: np.log(s).diff())
d['dln_iNDVI']=g['iNDVI'].transform(lambda s: np.log(s).diff()); d['dln_WP']=d.dln_NIRv-d.dln_ET
c=d.dropna(subset=['d_WP']).copy()
c['BIc']=c.BI-c.BI.mean(); c['dETc']=c.d_ET-c.d_ET.mean(); c['dET_x_BI']=c.dETc*c.BIc
c['dlnETc']=c.dln_ET-c.dln_ET.mean(); c['dlnET_x_BI']=c.dlnETc*c.BIc
chs={'G_FD_districtFE':'d_WP ~ d_ET + d_LST + d_Rain + BIc + dET_x_BI + C(District)',
     'H_FD_twowayFE':'d_WP ~ d_ET + d_LST + d_Rain + BIc + dET_x_BI + C(District) + C(BoroYear)',
     'I_FD_iNDVI':'d_WPi ~ d_ET + d_LST + d_Rain + BIc + dET_x_BI + C(District) + C(BoroYear)',
     'J_elasticity_NIRv':'dln_NIRv ~ dln_ET + d_LST + d_Rain + BIc + dlnET_x_BI + C(District)',
     'K_elasticity_NIRv_twFE':'dln_NIRv ~ dln_ET + d_LST + d_Rain + BIc + dlnET_x_BI + C(District) + C(BoroYear)',
     'L_elasticity_iNDVI':'dln_iNDVI ~ dln_ET + d_LST + d_Rain + BIc + dlnET_x_BI + C(District)',
     'M_elasticity_iNDVI_twFE':'dln_iNDVI ~ dln_ET + d_LST + d_Rain + BIc + dlnET_x_BI + C(District) + C(BoroYear)'}
chm={}
for k,f in chs.items():
    m=smf.ols(f,c).fit(cov_type='HC1'); chm[k]=m
    for v in ['d_ET','dln_ET','d_LST','d_Rain','BIc','dET_x_BI','dlnET_x_BI']:
        if v in m.params: record('change',k,v,m.params[v],m.bse[v],m.pvalues[v],dict(n=int(m.nobs),r2=m.rsquared))
    if 'dln_ET' in m.params:
        t=(m.params['dln_ET']-1)/m.bse['dln_ET']; from scipy import stats
        record('change',k,'H0: elasticity=1',m.params['dln_ET']-1,m.bse['dln_ET'],2*stats.norm.sf(abs(t)))
for k in ['G_FD_districtFE','H_FD_twowayFE','J_elasticity_NIRv','K_elasticity_NIRv_twFE']:
    v='d_ET' if 'FD' in k else 'dln_ET'; f=chs[k]
    record('change',k.replace('_',' ',0)+'_WCB',v,chm[k].params[v],np.nan,wcb(f,c,v))
    iv='dET_x_BI' if 'FD' in k else 'dlnET_x_BI'
    record('change',k+'_WCB',iv,chm[k].params[iv],np.nan,wcb(f,c,iv))
cp=c.copy()
m=dk('d_ET + d_LST + d_Rain + BIc + dET_x_BI','d_WP',cp,True,False)
for v in ['d_ET','dET_x_BI','BIc']: record('change','G_FD_districtFE_DK',v,m.params[v],m.std_errors[v],m.pvalues[v])
m=dk('dln_ET + d_LST + d_Rain + BIc + dlnET_x_BI','dln_NIRv',cp,True,False)
for v in ['dln_ET','dlnET_x_BI']: record('change','J_elasticity_NIRv_DK',v,m.params[v],m.std_errors[v],m.pvalues[v])
for dist in d.District.unique():
    m=smf.ols(chs['J_elasticity_NIRv'],c[c.District!=dist]).fit(cov_type='HC1')
    record('change','LODO_drop_'+dist,'dln_ET',m.params['dln_ET'],m.bse['dln_ET'],m.pvalues['dln_ET'])
# per-district elasticity (simple)
for dist,s in c.groupby('District'):
    m=smf.ols('dln_NIRv ~ dln_ET + d_Rain + d_LST',s).fit(cov_type='HC1')
    record('change','district_'+dist,'dln_ET',m.params['dln_ET'],m.bse['dln_ET'],m.pvalues['dln_ET'],dict(n=int(m.nobs)))
res=pd.DataFrame(rows); res.to_csv('regression_results.csv',index=False)
pd.set_option('display.width',250); print(res.to_string())
c.to_csv('change_panel.csv',index=False)
