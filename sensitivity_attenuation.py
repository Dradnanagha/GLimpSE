"""
One-way sensitivity analysis of the GLiMPSE real-world adjustment layer.
Uses the archived simulator code (synthetic data, seed 42; models random_state=1) unchanged.
Run from the repository root:  python sensitivity_attenuation.py
Real-world estimate = trial-style estimate x factor x (1 - k x discontinuation probability)
Base values: HbA1c factor 0.55 (k = 0.5); weight factor 0.72 (k = 0.4).
"""
import numpy as np, pandas as pd
import app as g   # the GLiMPSE simulator (app.py in this repository)

AGENTS = g.AGENT_KEYS
REF = dict(age=55, sex=0.0, bmi=34, weight=95, hba1c=8.2, duration=6,
           insulin=0.0, su=0.0, metformin=1.0, adherence=0.85)   # app default profile
GRID = [round(x, 2) for x in np.arange(0.40, 1.0001, 0.05)] + [0.72]
GRID = sorted(set(GRID))

def predict_rows(X):
    return (g.H1_M.predict(X), g.WL_M.predict(X),
            g.GI_M.predict_proba(X)[:, 1], g.DC_M.predict_proba(X)[:, 1])

rows = []
# reference profile
for ag in AGENTS:
    x = np.array([g.featurize(dict(REF, agent=ag))])
    h1, wl, gi, dc = [float(v[0]) for v in predict_rows(x)]
    for f in GRID:
        rows.append(dict(level='reference profile', agent=ag, factor=f,
                         h1_trial=h1, wl_trial=wl, gi_prob=gi, disc_prob=dc,
                         h1_adj=h1 * f * (1 - 0.5 * dc), wl_adj=wl * f * (1 - 0.4 * dc)))
# synthetic cohort (n = 2,500), means by agent
X = g._DF[g.FEATURES].values
h1v, wlv, giv, dcv = predict_rows(X)
onehot = g._DF[[c for c in g.FEATURES if c.startswith('ag_')]].values
for j, ag in enumerate(AGENTS):
    m = onehot[:, j] == 1
    for f in GRID:
        rows.append(dict(level='synthetic cohort mean', agent=ag, factor=f, n=int(m.sum()),
                         h1_trial=h1v[m].mean(), wl_trial=wlv[m].mean(),
                         gi_prob=giv[m].mean(), disc_prob=dcv[m].mean(),
                         h1_adj=(h1v[m] * f * (1 - 0.5 * dcv[m])).mean(),
                         wl_adj=(wlv[m] * f * (1 - 0.4 * dcv[m])).mean()))
df = pd.DataFrame(rows)
df['h1_retained'] = df.h1_adj / df.h1_trial
df['wl_retained'] = df.wl_adj / df.wl_trial
df.to_csv('GLiMPSE_attenuation_sensitivity_grid.csv', index=False, float_format='%.4f')

# discontinuation-drag coefficients (reference profile, base factors)
drag = []
for ag in AGENTS:
    x = np.array([g.featurize(dict(REF, agent=ag))])
    h1, wl, gi, dc = [float(v[0]) for v in predict_rows(x)]
    for kh, kw in [(0.0, 0.0), (0.25, 0.2), (0.5, 0.4), (0.75, 0.6), (1.0, 0.8)]:
        drag.append(dict(agent=ag, k_h1=kh, k_wl=kw, disc_prob=dc,
                         h1_adj=h1 * 0.55 * (1 - kh * dc), wl_adj=wl * 0.72 * (1 - kw * dc)))
pd.DataFrame(drag).to_csv('GLiMPSE_drag_sensitivity.csv', index=False, float_format='%.4f')
print('written', len(df), 'grid rows and', len(drag), 'drag rows')
