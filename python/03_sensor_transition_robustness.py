"""
Supplementary Table 12 | Robustness of the main estimates to the Landsat 7–8 sensor transition
------------------------------------------------------------------------------------------------
Input : ../data/Boro_district_season_panel_2001_2025.csv  (district-season panel, 4 districts x 25 Boro seasons, 2001–2025)
        required columns: District, BoroYear, BoroET_ERA5Land, BoroRain_mm, Mean_LST (K),
                          Boro_intensity_pct, Max_NIRv
Output: SuppTable12_sensor_robustness.csv  (+ printed table)

Tests
  (i)   exclude first differences that span the transition (ending 2013; 2012–2014; 2012–2015)
  (ii)  district-specific sensor offsets: district x 2013 indicators (change models),
        district x post-2012 indicators (level model)
  (iii) elasticity estimated within each sensor era (2002–2012; 2014–2025)
All models: OLS with HC1 standard errors; controls and fixed effects as in Table 1.
Requires: pandas, numpy, statsmodels, scipy   (pip install pandas numpy statsmodels scipy)
"""
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats

# ------------------------------------------------------------------ 1. panel
d = pd.read_csv('../data/Boro_district_season_panel_2001_2025.csv').sort_values(['District', 'BoroYear']).reset_index(drop=True)
d['ET'] = d.BoroET_ERA5Land                  # mm, Nov–May
d['Rain'] = d.BoroRain_mm                    # mm, Nov–May (CHIRPS)
d['LST'] = d.Mean_LST - 273.15               # °C
d['BI'] = d.Boro_intensity_pct               # %
d['NIRv'] = d.Max_NIRv
d['WPress'] = (d.ET - d.Rain) / d.ET         # water pressure
d['WP'] = d.NIRv / d.ET * 1e4                # NIRv-WP, x10^-4 mm^-1
d['post2013'] = (d.BoroYear >= 2013).astype(int)   # Landsat 8 era (seasons 2013–2025)
for v in ['BI', 'LST', 'Rain', 'WPress']:    # mean-centred level covariates
    d[v + '_c'] = d[v] - d[v].mean()

# ------------------------------------------------------------------ 2. first differences
g = d.groupby('District')
for v in ['WP', 'ET', 'LST', 'Rain']:
    d['d_' + v] = g[v].diff()
d['dln_NIRv'] = g['NIRv'].transform(lambda s: np.log(s).diff())
d['dln_ET'] = g['ET'].transform(lambda s: np.log(s).diff())
c = d.dropna(subset=['d_WP']).copy()         # 96 obs (differences end in 2002–2025)
c['BIc'] = c.BI - c.BI.mean()
c['dET_x_BI'] = (c.d_ET - c.d_ET.mean()) * c.BIc
c['dlnET_x_BI'] = (c.dln_ET - c.dln_ET.mean()) * c.BIc

# ------------------------------------------------------------------ 3. model formulas (Table 1)
LEVEL = 'WP ~ BI_c + LST_c + Rain_c + WPress_c + C(District) + C(BoroYear)'                     # col. 2
CHANGE = 'd_WP ~ d_ET + d_LST + d_Rain + BIc + dET_x_BI + C(District) + C(BoroYear)'            # col. 3
ELAST = 'dln_NIRv ~ dln_ET + d_LST + d_Rain + BIc + dlnET_x_BI + C(District) + C(BoroYear)'     # col. 4
ELAST_DFE = 'dln_NIRv ~ dln_ET + d_LST + d_Rain + BIc + dlnET_x_BI + C(District)'               # district FE only
OFF13 = ' + C(District):I(BoroYear == 2013)'   # district-specific sensor offset in differences
OFFPOST = ' + C(District):post2013'            # district-specific sensor offset in levels

rows = []
def fit(panel, spec, formula, data, term, test_eps1=False):
    m = smf.ols(formula, data).fit(cov_type='HC1')
    b, se = m.params[term], m.bse[term]
    p_eps1 = 2 * stats.norm.sf(abs((b - 1) / se)) if test_eps1 else np.nan
    rows.append(dict(Panel=panel, Specification=spec, Estimate=b, SE=se,
                     CI_low=b - 1.96 * se, CI_high=b + 1.96 * se,
                     P=m.pvalues[term], P_eps_eq_1=p_eps1, n=int(m.nobs)))

drop = lambda yrs: c[~c.BoroYear.isin(yrs)]

# ------------------------------------------------------------------ 4. estimates
P1 = 'Elasticity of NIRv with respect to ET (two-way FE unless stated)'
fit(P1, 'Main estimate (Table 1, column 4)', ELAST, c, 'dln_ET', True)
fit(P1, 'Excluding first differences ending in 2013', ELAST, drop([2013]), 'dln_ET', True)
fit(P1, 'Excluding first differences ending in 2012–2014', ELAST, drop([2012, 2013, 2014]), 'dln_ET', True)
fit(P1, 'Excluding first differences ending in 2012–2015', ELAST, drop([2012, 2013, 2014, 2015]), 'dln_ET', True)
fit(P1, 'District x 2013 sensor offsets', ELAST + OFF13, c, 'dln_ET', True)
fit(P1, 'District fixed effects only, excluding 2012–2014', ELAST_DFE, drop([2012, 2013, 2014]), 'dln_ET', True)
fit(P1, 'Landsat 5/7 era only (2002–2012)', ELAST, c[c.BoroYear <= 2012], 'dln_ET', True)
fit(P1, 'Landsat 7/8/9 era only (2014–2025)', ELAST, c[c.BoroYear >= 2014], 'dln_ET', True)

P2 = 'Effect of dET (mm) on dNIRv-WP (x10^-4 mm^-1)'
fit(P2, 'Main estimate (Table 1, column 3)', CHANGE, c, 'd_ET')
fit(P2, 'Excluding first differences ending in 2012–2014', CHANGE, drop([2012, 2013, 2014]), 'd_ET')
fit(P2, 'District x 2013 sensor offsets', CHANGE + OFF13, c, 'd_ET')

P3 = 'Level effect of Boro intensity (pp) on NIRv-WP (x10^-4 mm^-1)'
fit(P3, 'Main estimate (Table 1, column 2)', LEVEL, d, 'BI_c')
fit(P3, 'District x post-2012 sensor offsets', LEVEL + OFFPOST, d, 'BI_c')

# ------------------------------------------------------------------ 5. output
t = pd.DataFrame(rows)
t.to_csv('SuppTable12_sensor_robustness.csv', index=False)
pd.set_option('display.width', 200, 'display.max_colwidth', 50)
for panel, s in t.groupby('Panel', sort=False):
    print('\n' + panel)
    print(s.drop(columns='Panel').to_string(index=False, float_format=lambda x: f'{x:.4g}'))
