# Results (Genetic Oscillator)


## Part 1 — Deterministic ODE (BDF)
Numerical solution of the deterministic model shows sustained oscillations in activator **A** and repressor **R**.

![Activator A — deterministic](assets/img/(a)A.png)
![Repressor R — deterministic](assets/img/(a)R.png)

---

## Part 2 — Stochastic SSA (Gillespie)
Stochastic simulations reproduce the oscillatory behavior with run-to-run variability characteristic of intrinsic noise.

![Activator A — stochastic](assets/img/(b)A.png)
![Repressor R — stochastic](assets/img/(b)R.png)

---

## Part 3 — Parameter change
With $\delta_R = 0.05\,\mathrm{h}^{-1}$, we compare deterministic vs. stochastic dynamics.

![Repressor R — deterministic, $\delta_R=0.05$](assets/img/(c)deterministic.png)
![Repressor R — stochastic, $\delta_R=0.05$](assets/img/(c)stochastic.png)
