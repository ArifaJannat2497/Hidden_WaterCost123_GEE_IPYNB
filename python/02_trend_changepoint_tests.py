import pandas as pd, numpy as np, pymannkendall as mk
d=pd.read_csv('analysis_panel.csv')
def pettitt(x):
    x=np.asarray(x); n=len(x); U=[np.sum(np.sign(x[t+1:,None]-x[None,:t+1]).sum()) for t in range(n-1)]
    # U_t = sum_{i<=t} sum_{j>t} sgn(x_j - x_i)
    U=np.array([sum(np.sign(x[j]-x[i]) for i in range(t+1) for j in range(t+1,n)) for t in range(n-1)])
    K=np.max(np.abs(U)); t=int(np.argmax(np.abs(U))); p=2*np.exp(-6*K**2/(n**3+n**2)); return t,K,min(p,1)
rows=[]
vars_={'BI':'Boro intensity (%)','ET':'Boro-season ET (mm)','Rain':'Boro-season rainfall (mm)','PET':'Boro-season PET (mm)','ETPET':'ET/PET','LST':'Mean LST (°C)','NIRv':'Max NIRv','iNDVI':'Integrated NDVI','WP':'NIRv-WP (×10⁻⁴ mm⁻¹)','WPress':'Water pressure (ET−P)/ET'}
for dist,s in d.groupby('District'):
    s=s.sort_values('BoroYear')
    for v,lab in vars_.items():
        r=mk.original_test(s[v].values); m2=mk.hamed_rao_modification_test(s[v].values)
        t,K,pp=pettitt(s[v].values)
        e=s[s.BoroYear<=2012][v].mean(); l=s[s.BoroYear>=2013][v].mean()
        rows.append(dict(District=dist,var=v,label=lab,mean=s[v].mean(),sd=s[v].std(),sen_slope=r.slope,tau=r.Tau,p_MK=r.p,p_MK_HR=m2.p,pettitt_year=int(s.BoroYear.iloc[t+1]) ,pettitt_p=pp,mean_2001_12=e,mean_2013_25=l,pct_change=(l-e)/e*100))
T=pd.DataFrame(rows); T.to_csv('trend_results.csv',index=False)
pd.set_option('display.width',250)
print(T[['District','var','mean','sen_slope','p_MK','p_MK_HR','pettitt_year','pettitt_p','mean_2001_12','mean_2013_25','pct_change']].round(4).to_string())
