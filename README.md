# GLimpSE: A GLP-1RA Reality Check
**A pharmacovigilance-aware educational simulator for GLP-1 receptor agonist prescribing**
*Short name: GLiMPSE. Author: Dr Adnan Agha, Consultant Endocrinologist (FRCP London & Glasgow).*

> ⚠️ **Educational tool — not for clinical use.** This simulator is trained on **synthetic** data generated from published trial averages. It has **no clinical predictive validity** and must not be used to guide the care of any individual patient. No patient-identifiable data are used anywhere in this project.

---

## What this is

Trial numbers describe the *promise* of a drug; the real world describes the *lesson*. This project pairs a small prescribing **simulator** with an independent **FDA safety-signal analysis** to make one teaching point concrete: a model can look impressive on efficacy and still be unsafe to rely on if it ignores tolerability, discontinuation, and the gap between trial and real-world performance.

Two parts:
1. **The simulator** — estimates expected weight loss and HbA1c change for a chosen GLP-1 agent, with a toggle between *trial-style* and *attenuated real-world* outputs, plus a tolerability/safety flag.
2. **The FAERS consistency check** — a reproducible disproportionality analysis of the openFDA Adverse Event Reporting System across the whole GLP-1 / dual GIP–GLP-1 class, used to show that the simulator's safety dimension is *directionally* consistent with real-world reporting (face validity, not clinical validation).

---

## What's in this repository

| File | Purpose |
|------|---------|
| `index.html` | The public web simulator (single file, no build step). Open locally or via the live link below. |
| `app.py` | Interactive Gradio version of the simulator (for Hugging Face Spaces or local use). |
| `requirements.txt` | Python dependencies for `app.py` and the FAERS script. |
| `faers_glp1_disproportionality.py` | Reproducible FAERS disproportionality engine (ROR, PRR, IC, EBGM) per agent and pooled. |
| `RUN_FAERS_WINDOWS.md` | Detailed step-by-step guide for running the FAERS script on Windows (PowerShell). |
| `LICENSE` | MIT licence. |
| `data/` *(optional)* | The analysis outputs (`class_level_disproportionality.csv`, `per_agent_disproportionality.csv`, `signals_summary.csv`, `geographic_comparison.csv`) if you choose to archive them alongside the code. |

---

## Three ways to use it

### 1. Live web app (no install)
Open `index.html` in any modern browser, or visit the hosted version:

    https://glp-sim.netlify.app/

### 2. Interactive Gradio app (local)
    pip install -r requirements.txt
    python app.py

Then open the local URL it prints (usually http://127.0.0.1:7860). Hosted version:

    https://huggingface.co/spaces/naturally-intuitive/glimpse

### 3. Reproduce the FAERS analysis
The script reads an optional openFDA API key from an environment variable (a key raises the rate limit; the script also runs without one). **Never commit your key** — it is read from the environment, not the code.

macOS / Linux:

    export OPENFDA_API_KEY="your_free_key"     # optional; get one at https://open.fda.gov/apis/authentication/
    pip install -r requirements.txt
    python faers_glp1_disproportionality.py

Windows (PowerShell):

    $env:OPENFDA_API_KEY="your_free_key"
    pip install -r requirements.txt
    python faers_glp1_disproportionality.py

The script writes the per-agent, class-level, signal-summary, and geographic CSVs. Full Windows walkthrough is in `RUN_FAERS_WINDOWS.md`.

---

## Agents modelled (simulator)

The simulator distinguishes dose and route rather than drug name alone:

| Option | Trial weight anchor | HbA1c anchor |
|--------|--------------------|--------------|
| Semaglutide 1 mg injectable (Ozempic) | ≈ 6% (SURPASS-2) | −1.86% (SURPASS-2) |
| Semaglutide 2.4 mg injectable (Wegovy) | 13.7% (SURMOUNT-5) | −1.86% (SURPASS-2) |
| Oral semaglutide 14 mg (Rybelsus) | ≈ 3% (PIONEER 1) | −1.4% (PIONEER 1) |
| Tirzepatide 10 mg | 19.0% (SURMOUNT-5) | −2.24% (SURPASS-2) |
| Tirzepatide 15 mg | 20.2% (SURMOUNT-5) | −2.30% (SURPASS-2) |

All magnitudes are illustrative trial-derived anchors. The real-world toggle attenuates them to reflect the documented trial-to-practice gap.

---

## Data and sources

- **Safety data:** US FDA Adverse Event Reporting System via the public openFDA interface (aggregate counts only).
- **Reporting standard:** disproportionality results follow the READUS-PV statement.
- **Efficacy anchors:** SURMOUNT-5, SURPASS-2, and PIONEER 1.
- The FAERS analysis is **documented for transparency and reproducibility** (it is exploratory and hypothesis-generating, not a pre-registered confirmatory study).

---

## How to cite

- **Manuscript:** Agha A. A pharmacovigilance-aware educational simulator for GLP-1 receptor agonist prescribing: development and consistency with FDA safety signals. *Submitted.* 
- **Code/data archive:** [10.5281/zenodo.20772220]

---

## Licence

MIT — see `LICENSE`.

## Author

Dr Adnan Agha, Consultant Endocrinologist; College of Medicine and Health Sciences, United Arab Emirates University. Director Clinical Trials Unit CMHS. Chair CardioMetabolic Research Priority Group. Contact via LinkedIn: https://www.linkedin.com/in/adnan-agha-1354b6187/
