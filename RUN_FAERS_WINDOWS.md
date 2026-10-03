# FAERS GLP-1 Disproportionality Engine

Reproducible pharmacovigilance signal detection for the full GLP-1 / dual GIP-GLP-1
receptor agonist class, using the public US FDA Adverse Event Reporting System (FAERS).
Computes four standard disproportionality metrics (ROR, PRR, IC, EBGM) per agent and
pooled at class level, across all reported MedDRA Preferred Terms plus a curated signal list.

No patient-identifiable data are handled. openFDA returns aggregate counts only.

================================================================================
## QUICK START FOR WINDOWS (PowerShell) -- step by step
================================================================================

You do NOT need to be a programmer. Follow these in order.

### Step 0. Check that Python is installed
In PowerShell, type:

    py --version

- If you see something like `Python 3.12.4`, you are ready. Go to Step 1.
- If you see "py is not recognized", install Python first:
  - Go to https://www.python.org/downloads/ , download the latest Windows installer,
    run it, and on the FIRST screen tick the box **"Add python.exe to PATH"**, then install.
  - Close and reopen PowerShell, and run `py --version` again.

### Step 1. Put the three files together in one folder
Create a folder, for example `C:\Users\<you>\faers-glp1`, and save all three files into it:
- `faers_glp1_disproportionality.py`
- `requirements.txt`
- `README.md`  (this file)

### Step 2. Move PowerShell into that folder
Use `cd` (change directory) with your real path:

    cd C:\Users\<you>\faers-glp1

Tip: type `cd ` then drag the folder from File Explorer into the PowerShell window; it
pastes the path for you. Press Enter.

### Step 3. (Recommended) Create a clean, isolated environment
This keeps these packages separate from the rest of your system.

    py -m venv venv
    .\venv\Scripts\Activate.ps1

You should now see `(venv)` at the start of the prompt line.

- If you get a red error about "running scripts is disabled on this system", run this once,
  then repeat the Activate line above:

    Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass

(This relaxes the policy for the current window only; it resets when you close PowerShell.)

### Step 4. Install the required packages
    py -m pip install -r requirements.txt

This downloads requests, pandas, numpy, and scipy. It runs once; afterwards it is cached.

### Step 5. Set your openFDA API key
You already did this correctly. The key must be in quotes:

    $env:OPENFDA_API_KEY="PASTE_YOUR_KEY_HERE"

- This lasts only for the CURRENT PowerShell window. If you close it, set it again.
- To check it is set:  `echo $env:OPENFDA_API_KEY`  (it should print your key).
- To make it permanent across all future windows (optional), run once, then open a NEW window:

    setx OPENFDA_API_KEY "PASTE_YOUR_KEY_HERE"

A free key is requested at https://open.fda.gov/apis/authentication/ (email only, instant).
Without a key the script still runs but is limited to ~1,000 requests/day instead of 120,000.

### Step 6. Run the analysis
    py faers_glp1_disproportionality.py

What you will see, in order:
- `Total FAERS reports in window ...: 18,xxx,xxx`
- For each agent: `[SEMAGLUTIDE] reports for drug query: ...` then `evaluating N events`
  and progress like `...25/130`.
- It works through every agent, then the pooled class, then prints
  `Done. NN drug-event signals flagged by >=1 metric.`

The FIRST run pulls live from the FDA API and can take several minutes (it makes many small
requests). Every response is saved to `out\contingency_cache.json`, so if you run it again it
finishes almost instantly and reproduces the exact same numbers.

### Step 7. Open your results
A new `out\` folder appears in the same directory. Open the CSV files in Excel:

| File | What it contains |
|---|---|
| `signals_summary.csv` | The headline: only drug-event pairs flagged by >=1 metric, strongest first |
| `per_agent_disproportionality.csv` | Every agent x event pair with all four metrics + confidence intervals |
| `class_level_disproportionality.csv` | The whole class pooled vs the rest of the database |
| `contingency_cache.json` | Saved raw counts (do not edit; gives instant, reproducible re-runs) |

--------------------------------------------------------------------------------
### Windows troubleshooting
--------------------------------------------------------------------------------
- **"py is not recognized"** -> Python is not installed or not on PATH. Reinstall Python and
  tick "Add python.exe to PATH" (Step 0). You can also try `python` instead of `py`.
- **"pip is not recognized"** -> use `py -m pip ...` instead of `pip ...` (as shown above).
- **Activate.ps1 is blocked** -> run the `Set-ExecutionPolicy -Scope Process ...` line in Step 3.
- **It seems stuck / slow** -> normal on the first run; it is making many API calls. Leave it.
  The progress counter (e.g. `...50/130`) shows it is working. Re-runs are fast (cache).
- **"429 Too Many Requests"** -> you hit the daily limit. Make sure your key is set (Step 5);
  partial progress is saved in the cache, so just run it again later to continue.
- **0 reports for an agent** -> expected for very new drugs (orforglipron) or unusual spellings.
- **Start completely fresh** -> delete the whole `out\` folder, then run again.

### macOS / Linux equivalents
Same steps; only these lines differ:
- Activate venv:  `python3 -m venv venv && source venv/bin/activate`
- Set key:        `export OPENFDA_API_KEY="PASTE_YOUR_KEY_HERE"`
- Run:            `python3 faers_glp1_disproportionality.py`

================================================================================
## HOW THE DATA IS OBTAINED FROM FDA -- two routes
================================================================================

### Route A -- openFDA API (used by this script; fastest, fully reproducible)
The openFDA `drug/event` endpoint exposes the entire FAERS database as a query. The script
never downloads files; it asks the API for the aggregate counts it needs to build each 2x2
table. The "download" is just setting your API key (Step 5). Fields used:
- Drug exposure: `patient.drug.openfda.substance_name` (active ingredient, normalised)
- Reaction: `patient.reaction.reactionmeddrapt.exact` (MedDRA Preferred Term)
- Country (for the geographic version): `occurcountry` (two-letter code, e.g. AE, IN, US)
- Denominator window: `receivedate:[20040101 TO 20251231]` (edit DATE_FROM / DATE_TO at top)

### Route B -- Raw FAERS Quarterly Data Extract Files (publication-grade)
For the final manuscript table, many reviewers prefer raw FAERS with explicit de-duplication,
which the API does not expose.
1. Search the FDA site for "FAERS Quarterly Data Extract Files" (data run from 2004 Q4).
2. Download the ASCII ZIP per quarter. Each has seven `$`-delimited tables: DEMO, DRUG, REAC,
   OUTC, RPSR, THER, INDI, joined on `PRIMARYID`.
3. De-duplicate: keep only the latest version of each case by `CASEID` + `CASEVERSION`, and
   drop the quarter's deleted cases.
4. Build the same 2x2 tables from the de-duplicated DRUG x REAC join and feed them to the
   `metrics_from_2x2`, `add_ebgm`, and `flag_signals` functions (they are data-source-agnostic).

Recommendation: Route A for the fast signal screen and live demo; Route B for the final
publication table if reviewers request de-duplicated counts.

================================================================================
## GEOGRAPHIC ANALYSIS (country stratification)
================================================================================
With `RUN_GEOGRAPHIC = True` (default), the script also stratifies the pooled-class
analysis by reporting country (`occurcountry`) over the curated signal set, and writes:

| File | What it contains |
|---|---|
| `geographic_comparison.csv` | Class-vs-rest disproportionality computed WITHIN each stratum |
| `GEOGRAPHIC_CAVEATS.txt` | The interpretation caveats (read before citing) |

Strata produced: Global, UAE, Rest-of-World (Global minus UAE), Gulf/Asian (AE, SA, IN,
CN, JP) and Western (US, GB, DE, FR, CA). Each row carries `n_class_reports`, `N_stratum`,
and an `unstable` flag (set where class reports < 50 or a < 3).

READ THE CAVEATS. Country-stratified FAERS is confounded by reporting behaviour (the US
dominates the database; reporting culture and requirements differ by country; `occurcountry`
is often missing; Gulf strata are small and unstable). These contrasts are HYPOTHESIS-
GENERATING ONLY and do not establish ethnic or physiological differences; testing those
needs real regional outcome data under a THREC/IRB-approved study.

================================================================================
## METHODS NOTES
================================================================================
- **2x2 table.** a = drug AND event; b = drug AND not-event; c = not-drug AND event; d =
  neither. Built from openFDA `meta.results.total` counts. 0.5 continuity correction where a
  cell is zero.
- **ROR** = ad/bc; 95% CI from the log-odds SE. Signal: ROR025 > 1 and a >= 3.
- **PRR** = [a/(a+b)] / [c/(c+d)] with Yates chi-square. Signal: PRR >= 2, chi2 >= 4, a >= 3.
- **IC** (BCPNN) = log2(observed/expected); Noren 2006 heuristic CI. Signal: IC025 > 0.
- **EBGM** (Empirical-Bayes Gamma-Poisson, single-component GPS); EBGM is the posterior
  geometric mean, EB05/EB95 its 5th/95th percentiles. Signal: EB05 > 2. The five-parameter
  MGPS (DuMouchel 1999) can be swapped in via R `openEBGM`/`PhViD` if a journal requires it.
- **Curated PTs** are matched to the live FAERS vocabulary so spelling is exact; any not found
  verbatim are printed for MedDRA cross-check (e.g., the precise PT for NAION).

================================================================================
## LIMITATIONS (state these in the manuscript)
================================================================================
Disproportionality quantifies REPORTING association, not incidence or causation. FAERS has no
exposure denominator, and is subject to reporting and notoriety bias (strong for litigated or
media-prominent signals), duplicates, and missing data. Signals are hypotheses for further
study, not confirmed risks. This analysis supports an external-consistency check of the
educational simulator's safety flags; it is not a clinical validation.
