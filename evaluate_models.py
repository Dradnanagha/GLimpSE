"""
Internal performance of the simulator's four models on a held-out 25% partition.

Uses the synthetic cohort (n = 2,500, seed 42) and the model settings in app.py
(gradient boosting; regressors 180 trees, classifiers 160 trees, depth 3, random_state=1).
Split: sklearn train_test_split(test_size=0.25, random_state=42).
Run from the repository root:  python evaluate_models.py
"""
import app as g
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score, roc_auc_score
from sklearn.ensemble import GradientBoostingRegressor, GradientBoostingClassifier

df = g._DF
X = df[g.FEATURES].values
print("Model                         Metric   Value")
for label, col, kind in [("Weight-loss regressor", "weight_loss", "R2"), ("HbA1c regressor", "h1_red", "R2"),
                         ("GI-event classifier", "gi", "AUC"), ("Discontinuation classifier", "disc", "AUC")]:
    Xtr, Xte, ytr, yte = train_test_split(X, df[col].values, test_size=0.25, random_state=42)
    if kind == "R2":
        m = GradientBoostingRegressor(n_estimators=180, max_depth=3, random_state=1).fit(Xtr, ytr)
        v = r2_score(yte, m.predict(Xte))
    else:
        m = GradientBoostingClassifier(n_estimators=160, max_depth=3, random_state=1).fit(Xtr, ytr)
        v = roc_auc_score(yte, m.predict_proba(Xte)[:, 1])
    print(f"{label:29s} {kind:6s} {v:.2f}")
