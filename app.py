"""
GLP-1 Response & Safety Stress-Test
A Clinical-AI Validation Demonstration by Dr Adnan Agha
(Consultant Endocrinologist, FRCP London & Glasgow | Clinical-AI Validation Lead)

PURPOSE
-------
This Space is a DEMONSTRATION of clinical-AI validation methodology, not a clinical tool.
It shows how a tabular ML model can estimate GLP-1 efficacy AND how a *validation lens*
exposes where an efficacy-only model overstates real-world benefit by ignoring adverse
events and the trial-to-real-world gap.

DATA
----
The model is trained on SYNTHETIC data generated from published clinical-trial averages
(SURMOUNT-5, SURPASS-2 and related). No real patient data is used. Outputs are illustrative.

NOT A MEDICAL DEVICE. NOT FOR CLINICAL USE.
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import gradio as gr
from sklearn.ensemble import GradientBoostingRegressor, GradientBoostingClassifier

RNG = np.random.default_rng(42)

# ----------------------------------------------------------------------------
# Agent parameters grounded in published averages (illustrative magnitudes)
#   weight-loss % (magnitude), HbA1c reduction (% points), baseline GI AE rate
#   Sources: SURMOUNT-5 (sema -13.7% / tirz -20.2% weight); SURPASS-2 HbA1c
#   reductions (-1.86% sema 1mg; -2.01/-2.24/-2.30% tirz 5/10/15mg).
# ----------------------------------------------------------------------------
AGENTS = {
    "Semaglutide 1 mg injectable (Ozempic)":   dict(wl=6.0,  h1=1.86, gi=0.40),
    "Semaglutide 2.4 mg injectable (Wegovy)":  dict(wl=13.7, h1=1.86, gi=0.45),
    "Oral semaglutide 14 mg (Rybelsus)":       dict(wl=3.0,  h1=1.4,  gi=0.40),
    "Tirzepatide 10 mg":  dict(wl=19.0, h1=2.24, gi=0.38),
    "Tirzepatide 15 mg":  dict(wl=20.2, h1=2.30, gi=0.42),
}
AGENT_KEYS = list(AGENTS.keys())

FEATURES = ["age", "sex", "bmi", "weight", "hba1c", "duration",
            "insulin", "su", "metformin", "adherence",
            "ag_sema1", "ag_sema", "ag_oral", "ag_t10", "ag_t15"]

# Real-world shrinkage vs trial (real-world HbA1c reductions ~ half of trial;
# weight loss somewhat attenuated). Illustrative, for the validation lens.
RW_H1 = 0.55
RW_WL = 0.72


def _agent_onehot(agent):
    return [1.0 if agent == k else 0.0 for k in AGENT_KEYS]


def featurize(r):
    return [r["age"], r["sex"], r["bmi"], r["weight"], r["hba1c"], r["duration"],
            r["insulin"], r["su"], r["metformin"], r["adherence"],
            *_agent_onehot(r["agent"])]


def generate_synthetic_data(n=2500):
    rows, wl, h1, gi, disc = [], [], [], [], []
    for _ in range(n):
        agent = RNG.choice(AGENT_KEYS)
        p = AGENTS[agent]
        age = float(RNG.normal(54, 11)); age = float(np.clip(age, 25, 80))
        sex = float(RNG.integers(0, 2))
        bmi = float(np.clip(RNG.normal(34, 5), 25, 55))
        weight = float(np.clip(bmi * RNG.normal(2.7, 0.2), 60, 200))
        hba1c = float(np.clip(RNG.normal(8.0, 1.2), 5.8, 12.0))
        duration = float(np.clip(RNG.exponential(7), 0, 30))
        insulin = float(RNG.random() < 0.25)
        su = float(RNG.random() < 0.30)
        metformin = float(RNG.random() < 0.75)
        adherence = float(np.clip(RNG.normal(0.85, 0.12), 0.4, 1.0))

        bmi_factor = float(np.clip(1 + 0.012 * (bmi - 32), 0.80, 1.25))
        age_factor = float(np.clip(1 - 0.003 * (age - 50), 0.85, 1.10))
        weight_loss = -(p["wl"] * bmi_factor * age_factor * adherence + RNG.normal(0, 2.5))
        weight_loss = float(np.clip(weight_loss, -35, 0))

        h1_red = (p["h1"] * (0.6 + 0.12 * (hba1c - 7.0))
                  * (1 - 0.018 * duration)
                  * (0.85 if insulin else 1.0)
                  * (0.95 if su else 1.0)
                  * adherence + RNG.normal(0, 0.30))
        h1_red = float(np.clip(h1_red, 0, 4.0))

        gi_prob = float(np.clip(p["gi"] + 0.15 * (1 - adherence)
                                - 0.002 * (age - 50) + RNG.normal(0, 0.05), 0.05, 0.90))
        gi_event = float(RNG.random() < gi_prob)

        disc_prob = float(np.clip(0.10 + 0.30 * (gi_prob - 0.40)
                                  + 0.20 * (1 - adherence) + RNG.normal(0, 0.03), 0.02, 0.50))
        disc_event = float(RNG.random() < disc_prob)

        r = dict(age=age, sex=sex, bmi=bmi, weight=weight, hba1c=hba1c,
                 duration=duration, insulin=insulin, su=su, metformin=metformin,
                 adherence=adherence, agent=agent)
        rows.append(featurize(r)); wl.append(weight_loss); h1.append(h1_red)
        gi.append(gi_event); disc.append(disc_event)

    df = pd.DataFrame(rows, columns=FEATURES)
    df["weight_loss"] = wl; df["h1_red"] = h1; df["gi"] = gi; df["disc"] = disc
    return df


def train_models(df):
    X = df[FEATURES].values
    wl_m = GradientBoostingRegressor(n_estimators=180, max_depth=3, random_state=1).fit(X, df["weight_loss"])
    h1_m = GradientBoostingRegressor(n_estimators=180, max_depth=3, random_state=1).fit(X, df["h1_red"])
    gi_m = GradientBoostingClassifier(n_estimators=160, max_depth=3, random_state=1).fit(X, df["gi"])
    dc_m = GradientBoostingClassifier(n_estimators=160, max_depth=3, random_state=1).fit(X, df["disc"])
    return wl_m, h1_m, gi_m, dc_m


# Train once at startup (synthetic data -> fast, reproducible)
_DF = generate_synthetic_data()
WL_M, H1_M, GI_M, DC_M = train_models(_DF)


def _trajectory_fig(wl_trial, wl_rw):
    weeks = np.array([0, 8, 16, 24, 36, 52])
    tau = 16.0
    trial = wl_trial * (1 - np.exp(-weeks / tau))
    rw = wl_rw * (1 - np.exp(-weeks / tau))
    fig, ax = plt.subplots(figsize=(5.2, 3.2), dpi=120)
    ax.plot(weeks, trial, marker="o", linewidth=2.2, label="Trial-style estimate")
    ax.plot(weeks, rw, marker="s", linewidth=2.2, linestyle="--", label="Real-world-adjusted")
    ax.fill_between(weeks, trial, rw, alpha=0.12)
    ax.set_xlabel("Weeks since initiation"); ax.set_ylabel("Weight change (%)")
    ax.set_title("Illustrative weight-loss trajectory (synthetic)", fontsize=10)
    ax.axhline(0, color="#999", linewidth=0.8)
    ax.legend(fontsize=8, loc="lower left"); ax.grid(alpha=0.25)
    fig.tight_layout()
    return fig


def run_prediction(age, sex, bmi, weight, hba1c, duration, insulin, su, metformin, adherence, agent):
    sex_v = 1.0 if sex == "Male" else 0.0
    r = dict(age=age, sex=sex_v, bmi=bmi, weight=weight, hba1c=hba1c, duration=duration,
             insulin=float(bool(insulin)), su=float(bool(su)), metformin=float(bool(metformin)),
             adherence=adherence / 100.0, agent=agent)
    x = np.array([featurize(r)])

    wl_trial = float(WL_M.predict(x)[0])
    h1_trial = float(H1_M.predict(x)[0])
    gi_prob = float(GI_M.predict_proba(x)[0][1])
    disc_prob = float(DC_M.predict_proba(x)[0][1])

    # Validation lens: real-world adjustment + discontinuation drag
    h1_rw = h1_trial * RW_H1 * (1 - 0.5 * disc_prob)
    wl_rw = wl_trial * RW_WL * (1 - 0.4 * disc_prob)
    h1_gap = h1_trial - h1_rw

    flags = []
    if gi_prob > 0.50:
        flags.append("**Elevated GI adverse-event risk** — tolerability and titration strategy materially affect real-world results.")
    if disc_prob > 0.20:
        flags.append("**Elevated discontinuation risk** — an efficacy-only model overstates *population* benefit because a meaningful fraction stop therapy.")
    if h1_gap > 0.7:
        flags.append(f"**Large trial-to-real-world gap (~{h1_gap:.1f}% HbA1c)** — headline efficacy should be discounted for this profile.")
    if not flags:
        flags.append("No high-risk flags for this profile, but real-world attenuation still applies.")

    md = f"""
### ⚠️ Demonstration only — not for clinical use
*Synthetic model. Illustrative outputs. See **Methodology** below.*

#### Efficacy — trial-style estimate
- **HbA1c reduction:** −{h1_trial:.1f}% points
- **Weight change (52 wk):** {wl_trial:.1f}%

#### Safety — pharmacovigilance signals
- **GI adverse-event risk:** {gi_prob*100:.0f}%
- **Treatment-discontinuation risk:** {disc_prob*100:.0f}%

#### 🔍 Validation lens — what an efficacy-only model misses
| | Trial-style | Real-world-adjusted |
|---|---|---|
| HbA1c reduction | −{h1_trial:.1f}% | **−{h1_rw:.1f}%** |
| Weight change | {wl_trial:.1f}% | **{wl_rw:.1f}%** |

{chr(10).join('- ' + f for f in flags)}

> The point of this demo: a model can be highly *accurate* on efficacy and still be *unsafe to deploy* if it ignores adverse events, discontinuation, and the trial-to-real-world gap. Validating for that is the work.
"""
    return md, _trajectory_fig(wl_trial, wl_rw)


DISCLAIMER = """
# 🩺 GLP-1 Response & Safety Stress-Test
### A clinical-AI validation demonstration — Dr Adnan Agha, Consultant Endocrinologist (FRCP London & Glasgow)

> **⚠️ DEMONSTRATION ONLY — NOT FOR CLINICAL USE.** This is an educational demonstration of clinical-AI
> *validation methodology*. It is **not** a medical device, is **not** validated for patient care, and must
> **not** be used to make treatment decisions. The model is trained on **synthetic data generated from
> published trial averages** — not on real patient records. All outputs are illustrative.
"""

METHODOLOGY = """
### Methodology & data provenance
- **Synthetic data only.** Training data is generated from published clinical-trial average effect sizes
  (e.g., SURMOUNT-5 weight loss: tirzepatide −20.2% vs semaglutide −13.7%; SURPASS-2 HbA1c reductions
  −1.86% to −2.30%), with plausible predictor relationships and added noise. **No real patient data is used.**
- **Model.** Gradient-boosted trees estimate efficacy (HbA1c reduction, weight change) and adverse-event /
  discontinuation risk — the tabular approach well suited to structured clinical prediction.
- **The validation lens.** A real-world shrinkage factor and a discontinuation drag are applied to the
  trial-style estimate to show how efficacy-only models overstate population benefit. This is the core
  message: validation must stress-test for safety and real-world behaviour, not just headline accuracy.
- **Predictors used (illustrative):** baseline BMI and weight (strongest for weight loss); baseline HbA1c,
  diabetes duration, and concurrent insulin/sulfonylurea (for glycemic response); adherence; agent.

*For a real validated tool, this would require IRB/ethics-approved real-world data, prospective validation,
and regulatory assessment — services my practice provides under the DoH/THREC framework.*
"""

ABOUT = """
### About this demonstration
Built by **Dr Adnan Agha** — Consultant Endocrinologist (FRCP London & Glasgow; dual UK CCT) and clinical-AI
developer. I build and independently **validate** clinical-AI systems, with a focus on **GLP-1 / cardiometabolic
real-world evidence** and **medical-AI safety**. My fine-tuned model *Pentabrid* ranked **#3 globally on the
MedXpertQA leaderboard (ICML 2025)**.

**This tool exists to make one point tangible:** most predictive models in this space optimise for efficacy and
quietly ignore adverse events and the trial-to-real-world gap. Closing that gap — validating clinical AI for
*safety*, not just accuracy — is what I do for health-tech teams and medical-affairs groups.

**Independent clinical-AI validation · GLP-1 / cardiometabolic RWE · advisory & evaluation.**
Contact: adnanagha@uaeu.ac.ae · [LinkedIn URL] · Al Ain / Abu Dhabi, UAE
*(External engagements undertaken in line with institutional approval.)*
"""

with gr.Blocks(title="GLP-1 Response & Safety Stress-Test") as demo:
    gr.Markdown(DISCLAIMER)
    with gr.Row():
        with gr.Column(scale=1):
            gr.Markdown("#### Patient parameters *(hypothetical)*")
            agent = gr.Dropdown(AGENT_KEYS, value="Tirzepatide 15 mg", label="GLP-1 agent")
            age = gr.Slider(25, 80, value=55, step=1, label="Age (years)")
            sex = gr.Dropdown(["Female", "Male"], value="Female", label="Sex")
            bmi = gr.Slider(25, 55, value=34, step=0.5, label="Baseline BMI (kg/m²)")
            weight = gr.Slider(60, 200, value=95, step=1, label="Baseline weight (kg)")
            hba1c = gr.Slider(5.8, 12.0, value=8.2, step=0.1, label="Baseline HbA1c (%)")
            duration = gr.Slider(0, 30, value=6, step=1, label="Diabetes duration (years)")
            adherence = gr.Slider(40, 100, value=85, step=1, label="Expected adherence (%)")
            with gr.Row():
                insulin = gr.Checkbox(label="On insulin")
                su = gr.Checkbox(label="On sulfonylurea")
                metformin = gr.Checkbox(value=True, label="On metformin")
            btn = gr.Button("Run stress-test", variant="primary")
        with gr.Column(scale=1):
            out_md = gr.Markdown()
            out_plot = gr.Plot()

    with gr.Accordion("Methodology & data provenance", open=False):
        gr.Markdown(METHODOLOGY)
    with gr.Accordion("About / contact", open=False):
        gr.Markdown(ABOUT)

    btn.click(run_prediction,
              [age, sex, bmi, weight, hba1c, duration, insulin, su, metformin, adherence, agent],
              [out_md, out_plot])
    demo.load(run_prediction,
              [age, sex, bmi, weight, hba1c, duration, insulin, su, metformin, adherence, agent],
              [out_md, out_plot])

if __name__ == "__main__":
    demo.launch()
