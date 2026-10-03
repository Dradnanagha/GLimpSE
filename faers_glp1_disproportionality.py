#!/usr/bin/env python3
"""
FAERS GLP-1 class disproportionality engine (+ geographic stratification)
=========================================================================
Pulls counts from the public openFDA drug-event API and computes four standard
pharmacovigilance disproportionality metrics for the full GLP-1 / dual GIP-GLP-1
receptor agonist class, per agent and pooled at class level, across all reported
MedDRA Preferred Terms plus a curated signal list.

GEOGRAPHIC MODULE: also stratifies the class-level analysis by reporting country
(openFDA `occurcountry`) to contrast UAE / Gulf-Asian vs Western vs Rest-of-World.
Read the caveats: country-stratified FAERS is confounded by reporting behaviour and
is hypothesis-generating only (see GEOGRAPHIC_CAVEATS.txt and the on-screen banner).

Metrics: ROR, PRR (+chi-square), IC (BCPNN), EBGM (empirical-Bayes Gamma-Poisson).
Study author: A. Agha. Code written with AI assistance and verified by the author.
No patient-identifiable data are handled; openFDA returns aggregate counts only.

USAGE
-----
    pip install -r requirements.txt
    export OPENFDA_API_KEY="your_free_key"     # optional but recommended
    python faers_glp1_disproportionality.py
"""

from __future__ import annotations
import os, sys, time, json, math, urllib.parse
import requests
import numpy as np
import pandas as pd
from scipy import stats
from scipy.special import digamma
from scipy.optimize import minimize

# --------------------------------------------------------------------------- #
# 1. CONFIGURATION
# --------------------------------------------------------------------------- #
API = "https://api.fda.gov/drug/event.json"
API_KEY = os.environ.get("OPENFDA_API_KEY", "")
DATE_FROM, DATE_TO = "20040101", "20251231"
TOP_N_EVENTS = 150
MIN_A_FOR_SIGNAL = 3
OUTDIR = "out"
CACHE_PATH = os.path.join(OUTDIR, "contingency_cache.json")

RUN_GEOGRAPHIC = True          # set False to skip the country module
GEO_MIN_DRUG = 200            # below this many class reports a stratum is flagged 'unstable'

GLP1_CLASS = {
    "SEMAGLUTIDE":  ["Ozempic", "Wegovy", "Rybelsus"],
    "TIRZEPATIDE":  ["Mounjaro", "Zepbound"],
    "DULAGLUTIDE":  ["Trulicity"],
    "LIRAGLUTIDE":  ["Victoza", "Saxenda"],
    "EXENATIDE":    ["Byetta", "Bydureon"],
    "LIXISENATIDE": ["Adlyxin", "Lyxumia"],
    "ALBIGLUTIDE":  ["Tanzeum"],
    "ORFORGLIPRON": ["Foundayo"],
}

# Geographic strata (openFDA occurcountry uses ISO-2 codes)
UAE_CODE = "AE"
GULF_ASIAN = {"AE": "UAE", "SA": "Saudi Arabia", "IN": "India", "CN": "China", "JP": "Japan"}
WESTERN    = {"US": "USA", "GB": "UK", "DE": "Germany", "FR": "France", "CA": "Canada"}

# Clinical signal concepts -> keyword fragments matched against FAERS's OWN reaction
# vocabulary (stored in CAPITALS), so the exact spelling is always correct even when
# FAERS files a concept under a synonym (e.g. gastroparesis -> IMPAIRED GASTRIC EMPTYING).
SIGNAL_CONCEPTS = {
    "GI tolerability":            ["NAUSEA", "VOMITING", "DIARRHOEA", "CONSTIPATION",
                                   "ABDOMINAL PAIN", "DYSPEPSIA", "ERUCTATION",
                                   "GASTROOESOPHAGEAL REFLUX"],
    "Gastroparesis/hypomotility": ["GASTROPARESIS", "IMPAIRED GASTRIC EMPTYING",
                                   "GASTROINTESTINAL HYPOMOTILITY"],
    "Ileus/obstruction":          ["ILEUS", "INTESTINAL OBSTRUCTION", "BOWEL OBSTRUCTION",
                                   "GASTROINTESTINAL OBSTRUCTION"],
    "Pancreatitis":               ["PANCREATITIS"],
    "Biliary":                    ["CHOLELITHIASIS", "CHOLECYSTITIS", "BILIARY",
                                   "GALLBLADDER", "CHOLANGITIS"],
    "Ocular incl. NAION":         ["DIABETIC RETINOPATHY", "OPTIC ISCHAEMIC NEUROPATHY",
                                   "ISCHAEMIC OPTIC", "BLINDNESS", "VISUAL IMPAIRMENT"],
    "Aspiration":                 ["ASPIRATION"],
    "Renal/volume":               ["ACUTE KIDNEY INJURY", "RENAL IMPAIRMENT",
                                   "RENAL FAILURE", "DEHYDRATION"],
    "Psychiatric":                ["SUICIDAL", "DEPRESSION"],
    "Thyroid neoplasm":           ["THYROID CANCER", "THYROID NEOPLASM",
                                   "MEDULLARY THYROID", "THYROID CARCINOMA"],
    "Nutrition/muscle":           ["MALNUTRITION", "DECREASED APPETITE", "SARCOPENIA"],
    "Thrombosis/VTE":             ["THROMBOSIS", "EMBOLISM", "THROMBOTIC", "THROMBUS",
                                   "THROMBOPHLEBITIS"],
}

# Specific terms always tested directly, even if too rare to reach the top-1000
# vocabulary. Covers the headline safety signals and your VTE case-report lane.
MUST_TEST_TERMS = [
    "OPTIC ISCHAEMIC NEUROPATHY", "GASTROPARESIS", "ILEUS",
    "DEEP VEIN THROMBOSIS", "PULMONARY EMBOLISM", "VENOUS THROMBOSIS",
    "PORTAL VEIN THROMBOSIS", "THYROID CANCER", "PULMONARY ASPIRATION",
]

# --------------------------------------------------------------------------- #
# 2. HTTP LAYER (retry + on-disk cache)
# --------------------------------------------------------------------------- #
_session = requests.Session()
_cache: dict = {}

def _load_cache():
    global _cache
    if os.path.exists(CACHE_PATH):
        with open(CACHE_PATH) as f:
            _cache = json.load(f)

def _save_cache():
    os.makedirs(OUTDIR, exist_ok=True)
    with open(CACHE_PATH, "w") as f:
        json.dump(_cache, f)

def _get(params: dict, retries: int = 4) -> dict:
    if API_KEY:
        params = {**params, "api_key": API_KEY}
    key = urllib.parse.urlencode(sorted(params.items()))
    if key in _cache:
        return _cache[key]
    for attempt in range(retries):
        try:
            r = _session.get(API, params=params, timeout=30)
            if r.status_code == 404:
                out = {"meta": {"results": {"total": 0}}, "results": []}
                _cache[key] = out
                return out
            r.raise_for_status()
            out = r.json()
            _cache[key] = out
            time.sleep(0.34 if API_KEY else 0.6)
            return out
        except Exception as e:
            wait = 2 ** attempt
            sys.stderr.write(f"  retry {attempt+1}/{retries} after {wait}s ({e})\n")
            time.sleep(wait)
    raise RuntimeError(f"openFDA request failed: {params}")

def _date_clause() -> str:
    return f"receivedate:[{DATE_FROM} TO {DATE_TO}]"

def count(search: str, country_clause: str | None = None) -> int:
    """Total reports matching a Lucene search, within the date window and optional
    country clause. Empty search = denominator (all reports)."""
    parts = []
    if search:
        parts.append(f"({search})")
    parts.append(_date_clause())
    if country_clause:
        parts.append(country_clause)
    full = " AND ".join(parts)
    js = _get({"search": full, "limit": 1})
    return int(js["meta"]["results"]["total"])

# --------------------------------------------------------------------------- #
# 3. QUERY BUILDERS
# --------------------------------------------------------------------------- #
def drug_clause(substances: list[str]) -> str:
    terms = [f'patient.drug.openfda.substance_name:"{s}"' for s in substances]
    return "(" + " OR ".join(terms) + ")"

def event_clause(pt: str) -> str:
    return f'patient.reaction.reactionmeddrapt.exact:"{pt}"'

def country_clause(codes: list[str]) -> str:
    return "(" + " OR ".join(f"occurcountry:{c}" for c in codes) + ")"

def top_events(drug_search: str, n: int) -> list[str]:
    full = f"({drug_search}) AND {_date_clause()}"
    js = _get({"search": full, "count": "patient.reaction.reactionmeddrapt.exact"})
    return [row["term"] for row in js.get("results", [])[:n]]

def resolve_signal_terms(drug_search: str) -> list[str]:
    """Match each clinical concept to FAERS's actual reaction vocabulary (correct
    spelling guaranteed), then add the must-test terms. Prints a clear concept-by-
    concept report instead of a 'not found' warning."""
    full = f"({drug_search}) AND {_date_clause()}"
    js = _get({"search": full, "count": "patient.reaction.reactionmeddrapt.exact"})
    vocab = [row["term"] for row in js.get("results", [])]      # FAERS exact terms (CAPS)
    chosen, report = [], []
    for concept, keys in SIGNAL_CONCEPTS.items():
        hits = sorted({t for t in vocab if any(k in t for k in keys)})
        chosen += hits
        report.append(f"    {concept}: {len(hits)} term(s)"
                      + (f", e.g. {hits[0]}" if hits else " (none in this vocabulary)"))
    chosen += [t for t in MUST_TEST_TERMS if t not in chosen]
    sys.stderr.write("  signal concepts matched to FAERS vocabulary:\n"
                     + "\n".join(report) + "\n")
    return list(dict.fromkeys(chosen))

# --------------------------------------------------------------------------- #
# 4. DISPROPORTIONALITY METRICS
# --------------------------------------------------------------------------- #
def metrics_from_2x2(a, b, c, d):
    a_, b_, c_, d_ = (max(a, 0.5), max(b, 0.5), max(c, 0.5), max(d, 0.5))
    ror = (a_ * d_) / (b_ * c_)
    se = math.sqrt(1/a_ + 1/b_ + 1/c_ + 1/d_)
    ror_lo = math.exp(math.log(ror) - 1.96 * se)
    ror_hi = math.exp(math.log(ror) + 1.96 * se)
    prr = (a_ / (a_ + b_)) / (c_ / (c_ + d_))
    n = a + b + c + d
    chi2 = 0.0
    if n and (a + b) > 0 and (a + c) > 0:
        try:
            chi2 = stats.chi2_contingency(np.array([[a, b], [c, d]], float), correction=True)[0]
        except Exception:
            chi2 = 0.0
    expected = (a + b) * (a + c) / n if n else 0.5
    ic = math.log2((a + 0.5) / (expected + 0.5))
    ic_lo = ic - 3.3 * (a + 0.5) ** -0.5 - 2.0 * (a + 0.5) ** -1.5
    ic_hi = ic + 2.4 * (a + 0.5) ** -0.5 - 0.5 * (a + 0.5) ** -1.5
    return dict(a=a, b=b, c=c, d=d, expected=expected,
                ROR=ror, ROR025=ror_lo, ROR975=ror_hi, PRR=prr, chi2=chi2,
                IC=ic, IC025=ic_lo, IC975=ic_hi)

def add_ebgm(df: pd.DataFrame) -> pd.DataFrame:
    a = df["a"].to_numpy(float); E = df["expected"].to_numpy(float)
    E = np.where(E <= 0, 1e-6, E)
    def neg_ll(p):
        alpha, beta = np.exp(p)
        prob = beta / (beta + E)
        return -(stats.nbinom.logpmf(a, alpha, prob)).sum()
    res = minimize(neg_ll, x0=np.log([1.0, 1.0]), method="Nelder-Mead")
    alpha, beta = np.exp(res.x)
    post_shape, post_rate = alpha + a, beta + E
    df = df.copy()
    df["EBGM"] = np.exp(digamma(post_shape) - np.log(post_rate))
    df["EB05"] = stats.gamma.ppf(0.05, a=post_shape, scale=1.0/post_rate)
    df["EB95"] = stats.gamma.ppf(0.95, a=post_shape, scale=1.0/post_rate)
    return df

def flag_signals(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["signal_ROR"] = (df["ROR025"] > 1) & (df["a"] >= MIN_A_FOR_SIGNAL)
    df["signal_PRR"] = (df["PRR"] >= 2) & (df["chi2"] >= 4) & (df["a"] >= MIN_A_FOR_SIGNAL)
    df["signal_IC"]  = (df["IC025"] > 0) & (df["a"] >= MIN_A_FOR_SIGNAL)
    if "EB05" in df:
        df["signal_EBGM"] = (df["EB05"] > 2) & (df["a"] >= MIN_A_FOR_SIGNAL)
        df["signal_count"] = df[["signal_ROR","signal_PRR","signal_IC","signal_EBGM"]].sum(axis=1)
    else:
        df["signal_count"] = df[["signal_ROR","signal_PRR","signal_IC"]].sum(axis=1)
    return df

# --------------------------------------------------------------------------- #
# 5. MAIN (non-geographic) DRIVER
# --------------------------------------------------------------------------- #
def analyse(label: str, substances: list[str], N: int) -> pd.DataFrame:
    dsearch = drug_clause(substances)
    n_drug = count(dsearch)
    print(f"[{label}] reports for drug query: {n_drug:,}")
    if n_drug == 0:
        return pd.DataFrame()
    events = list(dict.fromkeys(top_events(dsearch, TOP_N_EVENTS) + resolve_signal_terms(dsearch)))
    print(f"[{label}] evaluating {len(events)} events")
    rows = []
    for i, pt in enumerate(events, 1):
        n_event = count(event_clause(pt))
        if n_event == 0:                      # term not present anywhere in FAERS
            continue
        a = count(f"{dsearch} AND {event_clause(pt)}")
        b, c = n_drug - a, n_event - a
        d = N - a - b - c
        if min(b, c, d) < 0:
            continue
        m = metrics_from_2x2(a, b, c, d)
        m.update(analysis=label, event=pt, n_drug=n_drug, n_event=n_event, N=N)
        rows.append(m)
        if i % 25 == 0:
            print(f"  ...{i}/{len(events)}")
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    df = flag_signals(add_ebgm(df))
    cols = ["analysis","event","a","b","c","d","expected","ROR","ROR025","ROR975",
            "PRR","chi2","IC","IC025","IC975","EBGM","EB05","EB95",
            "signal_ROR","signal_PRR","signal_IC","signal_EBGM","signal_count"]
    return df[cols].sort_values("signal_count", ascending=False)

# --------------------------------------------------------------------------- #
# 6. GEOGRAPHIC MODULE
# --------------------------------------------------------------------------- #
def geo_cells(substances, events, cclause):
    """Within-stratum raw cells for the class vs the rest of the database in that stratum."""
    dsearch = drug_clause(substances)
    N_s = count("", cclause)
    n_drug_s = count(dsearch, cclause)
    cells = {}
    for pt in events:
        n_event_s = count(event_clause(pt), cclause)
        a_s = count(f"{dsearch} AND {event_clause(pt)}", cclause)
        cells[pt] = dict(a=a_s, n_drug=n_drug_s, n_event=n_event_s, N=N_s)
    return cells, N_s, n_drug_s

def cells_to_rows(cells, stratum, n_drug_s, N_s):
    rows = []
    for pt, c in cells.items():
        a, n_drug, n_event, N = c["a"], c["n_drug"], c["n_event"], c["N"]
        b, cc, d = n_drug - a, n_event - a, N - n_drug - n_event + a
        if min(a, b, cc, d) < 0:
            continue
        m = metrics_from_2x2(a, b, cc, d)
        m.update(stratum=stratum, event=pt, n_class_reports=n_drug_s, N_stratum=N_s,
                 unstable=(n_drug_s < GEO_MIN_DRUG or a < 10))
        rows.append(m)
    return rows

def run_geographic(events):
    print("\n" + "=" * 70)
    print("GEOGRAPHIC STRATIFICATION  (hypothesis-generating only -- see caveats)")
    print("=" * 70)
    subs = list(GLP1_CLASS.keys())

    global_cells, _, ndg = geo_cells(subs, events, None);            _save_cache()
    uae_cells,    _, ndu = geo_cells(subs, events, country_clause([UAE_CODE])); _save_cache()
    gulf_cells, Ng, ndgulf = geo_cells(subs, events, country_clause(list(GULF_ASIAN))); _save_cache()
    west_cells, Nw, ndwest = geo_cells(subs, events, country_clause(list(WESTERN)));    _save_cache()

    # Rest-of-world by subtraction (global minus UAE); includes unspecified-country reports
    row_cells = {pt: dict(a=global_cells[pt]["a"] - uae_cells[pt]["a"],
                          n_drug=global_cells[pt]["n_drug"] - uae_cells[pt]["n_drug"],
                          n_event=global_cells[pt]["n_event"] - uae_cells[pt]["n_event"],
                          N=global_cells[pt]["N"] - uae_cells[pt]["N"]) for pt in events}
    ndrow = ndg - ndu; Nrow = global_cells[events[0]]["N"] - uae_cells[events[0]]["N"]

    rows = []
    rows += cells_to_rows(global_cells, "Global", ndg, global_cells[events[0]]["N"])
    rows += cells_to_rows(uae_cells, "UAE", ndu, uae_cells[events[0]]["N"])
    rows += cells_to_rows(row_cells, "Rest_of_World", ndrow, Nrow)
    rows += cells_to_rows(gulf_cells, "Gulf_Asian", ndgulf, Ng)
    rows += cells_to_rows(west_cells, "Western", ndwest, Nw)

    df = pd.DataFrame(rows)
    cols = ["stratum","event","n_class_reports","N_stratum","a","ROR","ROR025","ROR975",
            "PRR","IC","IC025","IC975","unstable"]
    df = df[cols].sort_values(["event","stratum"])
    df.to_csv(os.path.join(OUTDIR, "geographic_comparison.csv"), index=False)

    print(f"  class reports by stratum: Global {ndg:,} | UAE {ndu:,} | "
          f"Gulf/Asian {ndgulf:,} | Western {ndwest:,}")
    unstable_n = int(df["unstable"].sum())
    print(f"  geographic_comparison.csv written ({len(df)} rows; {unstable_n} flagged unstable)")
    _write_geo_caveats()

def _write_geo_caveats():
    txt = """GEOGRAPHIC ANALYSIS -- CAVEATS (read before interpreting or citing)

Country-stratified FAERS disproportionality is HYPOTHESIS-GENERATING ONLY. Differences
between strata reflect REPORTING behaviour at least as much as any clinical difference:

1. The United States dominates FAERS. Non-US reports are a small, non-random subset.
2. Reporting culture, regulatory reporting requirements, and pharmacovigilance maturity
   differ sharply by country, so a higher or lower reporting ratio is not a higher or
   lower clinical risk.
3. `occurcountry` is frequently missing; 'Rest of World' here is computed as Global minus
   UAE and therefore includes reports with no specified country.
4. UAE and other Gulf strata have small report counts, giving unstable estimates with wide
   confidence intervals (rows are flagged 'unstable' where class reports < %d or a < %d).
5. FAERS captures no exposure denominator, so none of these are incidence rates.

Conclusion: present geographic contrasts as a descriptive, caveated sub-analysis that
generates regionally specific hypotheses. They do NOT establish ethnic or physiological
differences in benefit or harm. Testing those requires real regional outcome data under a
THREC/IRB-approved real-world-evidence study.
""" % (GEO_MIN_DRUG, MIN_A_FOR_SIGNAL)
    with open(os.path.join(OUTDIR, "GEOGRAPHIC_CAVEATS.txt"), "w") as f:
        f.write(txt)

# --------------------------------------------------------------------------- #
# 7. MAIN
# --------------------------------------------------------------------------- #
def main():
    os.makedirs(OUTDIR, exist_ok=True)
    _load_cache()
    N = count("")
    print(f"Total FAERS reports in window {DATE_FROM}-{DATE_TO}: {N:,}\n")

    per_agent = []
    for sub in GLP1_CLASS:
        per_agent.append(analyse(sub, [sub], N)); _save_cache()
    per_agent_df = pd.concat([d for d in per_agent if not d.empty], ignore_index=True)
    class_df = analyse("GLP1_CLASS_POOLED", list(GLP1_CLASS.keys()), N); _save_cache()

    per_agent_df.to_csv(os.path.join(OUTDIR, "per_agent_disproportionality.csv"), index=False)
    class_df.to_csv(os.path.join(OUTDIR, "class_level_disproportionality.csv"), index=False)
    sig = pd.concat([per_agent_df, class_df], ignore_index=True)
    sig = sig[sig["signal_count"] >= 1].sort_values(
        ["analysis","signal_count","ROR"], ascending=[True, False, False])
    sig.to_csv(os.path.join(OUTDIR, "signals_summary.csv"), index=False)
    print(f"\nNon-geographic: {len(sig)} drug-event signals flagged by >=1 metric.")

    if RUN_GEOGRAPHIC:
        geo_events = resolve_signal_terms(drug_clause(list(GLP1_CLASS.keys())))
        run_geographic(geo_events)

    print("\nDone. Outputs in ./out/  "
          "(signals_summary, per_agent, class_level, geographic_comparison, GEOGRAPHIC_CAVEATS)")

if __name__ == "__main__":
    main()
