"""
Assignment C – Deterministic vs. Stochastic Comparison
------------------------------------------------------
Reproduce Fig. 5 from Vilar et al. (PNAS, 2002), showing how random
noise in the stochastic model can qualitatively change behavior compared
to the deterministic model.

Goal:
- Simulate the deterministic ODE system with adjusted parameters.
- Compare the repressor protein R over 400 hours.
- Highlight differences vs. the stochastic solution (see Assignment B).
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


# -------------------------------------------------------------------
# Parameters (note: delta_r adjusted compared to Assignment A)
# -------------------------------------------------------------------
alpha_a, alpha_ap, alpha_r, alpha_rp = 50, 500, 0.01, 50    # transcription
beta_a, beta_r = 50, 5                                      # translation
gamma_a, gamma_r, gamma_c = 1, 1, 2                         # binding
delta_a, delta_r, delta_ma, delta_mr = 1, 0.05, 10, 0.5     # degradation
phi_a, phi_r = 50, 100                                      # unbinding

# Initial conditions: [DA, DR, DAp, DRp, MA, MR, A, R, C]
Initial = [1, 1, 0, 0, 0, 0, 0, 0, 0]

FinalTime = 400  # hours


# -------------------------------------------------------------------
# ODE system (Eq. [1], modified parameters)
# -------------------------------------------------------------------
def molecules(t, y):
    """Right-hand side of the circadian clock ODE system (Assignment C)."""
    DA, DR, DAp, DRp, MA, MR, A, R, C = y
    yprime = np.zeros(9)

    # Gene promoter occupancy
    yprime[0] = phi_a * DAp - gamma_a * DA * A
    yprime[1] = phi_r * DRp - gamma_r * DR * A
    yprime[2] = gamma_a * DA * A - phi_a * DAp
    yprime[3] = gamma_r * DR * A - phi_r * DRp

    # mRNA
    yprime[4] = alpha_ap * DAp + alpha_a * DA - delta_ma * MA
    yprime[5] = alpha_rp * DRp + alpha_r * DR - delta_mr * MR

    # Proteins and complex
    yprime[6] = (
        beta_a * MA + phi_a * DAp + phi_r * DRp
        - A * (gamma_a * DA + gamma_r * DR + gamma_c * R + delta_a)
    )  # dA/dt
    yprime[7] = beta_r * MR - gamma_c * A * R + delta_a * C - delta_r * R  # dR/dt
    yprime[8] = gamma_c * A * R - delta_a * C                              # dC/dt

    return yprime


# -------------------------------------------------------------------
# Solve system
# -------------------------------------------------------------------
teval = np.linspace(0, FinalTime, 1000)
sol = solve_ivp(molecules, [0, FinalTime], Initial, method="BDF", t_eval=teval)


# -------------------------------------------------------------------
# Plot result
# -------------------------------------------------------------------
plt.figure(figsize=(8, 4))
plt.plot(sol.t, sol.y[7], linestyle="solid", color="red", label="R (repressor)")
plt.xlabel("Time (hours)")
plt.ylabel("Number of molecules")
plt.title("Deterministic simulation – Repressor R (Assignment C)")
plt.legend(loc="upper right")
plt.tight_layout()
plt.show()
