 Modeling the Circadian Clock: Deterministic and Stochastic Simulation

This project reproduces key results from the influential article  
**"Mechanisms of noise-resistance in genetic oscillators"**  
by José M. G. Vilar et al. (PNAS, 2002), using both deterministic and stochastic mathematical models.

🧬 Project Overview

Circadian clocks are molecular mechanisms used by organisms to synchronize biological processes with daily environmental cycles. In this mini-project, we model a simplified biological system involving two proteins — an **activator** and a **repressor** — and investigate their dynamic interactions over time using both:

- **Ordinary Differential Equations (ODEs)** — deterministic model  
- **Stochastic Simulation Algorithm (SSA)** — stochastic model

The aim is to explore how noise affects biological rhythms, and compare the behavior of deterministic and stochastic formulations of the same system.

⚙️ Methods

### 1. Deterministic Model (ODE)
- Implemented the system of stiff ODEs describing protein interactions over 400 hours.
- Solved using `scipy.integrate.solve_ivp` with the **BDF method** (for stiff problems).
- Reproduced **Figure 2a & 2b** from the article, showing concentration curves of proteins A and R.

### 2. Stochastic Model (SSA)
- Modeled 16 reactions using either:
  - A custom SSA implementation, or
  - The `Gillespy2` Python library in the **stochSS** environment.
- Reproduced **Figure 2c & 2d** via repeated SSA runs to observe variability across trials.

### 3. Comparison & Analysis
- Simulated and compared both models for long-term behavior.
- Reproduced **Figure 5**, showing qualitative divergence due to stochastic noise.
- Reflected on when stochastic modeling is more appropriate than deterministic.

📊 Files

- `circadian_ode.py` – ODE implementation of the deterministic model
- `circadian_ssa.py` – Stochastic implementation (SSA or Gillespy2)
- `figures/` – Plots replicating results (e.g., Fig 2a–2d, Fig 5)
- `report.pdf` – Final report with explanation, figures, and analysis
- `README.md` – Project overview (this file)

🔬 Biological Background

The model tracks:
- Activator gene \( D_A \), Repressor gene \( D_R \)
- Transcription, binding, and degradation of mRNA and proteins A and R
- Feedback loops where A activates, and R represses gene expression

Despite its biological complexity, the system is reduced to a manageable number of ODEs and reactions.

📦 Technologies Used

- Python 3  
- NumPy, SciPy (`solve_ivp`)  
- Matplotlib for plotting  
- Gillespy2 / StochSS (optional for SSA)

📈 Results Summary

- **Deterministic model**: Predictable oscillations in A and R with stable amplitude and frequency.
- **Stochastic model**: Similar average behavior, but exhibits variability between runs due to molecular noise.
- **Insight**: Noise-resistance mechanisms are crucial in biological systems; stochastic models are necessary for capturing cell-scale variability.

📚 Source Article

Vilar, J.M.G., Kueh, H.Y., Barkai, N., & Leibler, S. (2002).  
*Mechanisms of noise-resistance in genetic oscillators*.  
PNAS, 99(9), 5988–5992.  

---

Feel free to explore, adapt, or cite this code for educational or research purposes.
