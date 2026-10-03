"""
Recomputes the active-comparator (restricted-background) reporting odds ratios
from the deposited two-by-two counts in data/sensitivity_active_comparator.csv.

Counts were retrieved from openFDA on 15 July 2026 (class: GLP-1 and dual GIP/GLP-1
receptor agonists; comparator: other glucose-lowering drugs, as listed in the paper's
Methods). ROR = (a/b)/(c/d); 95% CI by the Woolf (log) method.
Run:  python active_comparator_ror.py
"""
import csv, math, os
path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "sensitivity_active_comparator.csv")
rows = list(csv.DictReader(open(path, newline="", encoding="utf-8-sig")))
all_match = True
print(f"{'event':32s} {'ROR':>6s} {'95% CI':>15s}   matches deposited")
for r in rows:
    a, b, c, d = (int(r[k]) for k in ("a_glp1_event", "b_glp1_noevent", "c_comp_event", "d_comp_noevent"))
    ror = (a / b) / (c / d)
    se = math.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    lo, hi = math.exp(math.log(ror) - 1.96 * se), math.exp(math.log(ror) + 1.96 * se)
    ok = (f"{ror:.2f}", f"{lo:.2f}", f"{hi:.2f}") == (f"{float(r['ROR_restricted']):.2f}", f"{float(r['CI_low']):.2f}", f"{float(r['CI_high']):.2f}")
    all_match &= ok
    print(f"{r['event']:32s} {ror:6.2f} {lo:7.2f}-{hi:<7.2f}   {ok}")
print("\nAll recomputed values match the deposited file:", all_match)
