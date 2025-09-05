
"""
Assignment A – Deterministic Model of Circadian Clock
-----------------------------------------------------
Reproduce Fig. 2a and Fig. 2b from Vilar et al. (PNAS, 2002) by solving
the stiff ODE model (Eq. [1]) over 400 hours.

Goal:
- Plot protein A (activator) and protein R (repressor) over time.
- Use SciPy's solve_ivp with the BDF method (suitable for stiff ODEs).
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

# -------------------------------------------------------------------
# Parameters (from article, Fig. 1 caption)
# -------------------------------------------------------------------
alpha_a, alpha_ap, alpha_r, alpha_rp = 50, 500, 0.01, 50   # transcription rates
beta_a, beta_r = 50, 5                                     # translation rates
gamma_a, gamma_r, gamma_c = 1, 1, 2                        # binding rates
delta_a, delta_r, delta_ma, delta_mr = 1, 0.2, 10, 0.5     # degradation rates
phi_a, phi_r = 50, 100                                     # unbinding rates

# Initial conditions: [DA, DR, DAp, DRp, MA, MR, A, R, C]
Initial = [1, 1, 0, 0, 0, 0, 0, 0, 0]

FinalTime = 400  # total simulation time (hours)


# -------------------------------------------------------------------
# ODE system (Eq. [1])
# -------------------------------------------------------------------
def molecules(t, y):
    """
    Right-hand side of the circadian clock ODE system.
    State vector y = [DA, DR, DAp, DRp, MA, MR, A, R, C]
    """
    DA, DR, DAp, DRp, MA, MR, A, R, C = y
    yprime = np.zeros(9)

    # Gene dynamics
    yprime[0] = phi_a * DAp - gamma_a * DA * A                # d(DA)/dt
    yprime[1] = phi_r * DRp - gamma_r * DR * A                # d(DR)/dt
    yprime[2] = gamma_a * DA * A - phi_a * DAp                # d(DAp)/dt
    yprime[3] = gamma_r * DR * A - phi_r * DRp                # d(DRp)/dt

    # mRNA dynamics
    yprime[4] = alpha_ap * DAp + alpha_a * DA - delta_ma * MA # d(MA)/dt
    yprime[5] = alpha_rp * DRp + alpha_r * DR - delta_mr * MR # d(MR)/dt

    # Protein dynamics
    yprime[6] = (
        beta_a * MA + phi_a * DAp + phi_r * DRp
        - A * (gamma_a * DA + gamma_r * DR + gamma_c * R + delta_a)
    )  # dA/dt

    yprime[7] = beta_r * MR - gamma_c * A * R + delta_a * C - delta_r * R  # dR/dt

    # Complex dynamics
    yprime[8] = gamma_c * A * R - delta_a * C                             # dC/dt

    return yprime


# -------------------------------------------------------------------
# Solve the ODE system
# -------------------------------------------------------------------
teval = np.linspace(0, FinalTime, 1000)  # evaluation time points
sol = solve_ivp(molecules, [0, FinalTime], Initial, method="BDF", t_eval=teval)


# -------------------------------------------------------------------
# Plot results
# -------------------------------------------------------------------
# Activator A (Fig. 2a)
plt.figure(figsize=(8, 4))
plt.plot(sol.t, sol.y[6], linestyle="solid", color="blue", label="A (activator)")
plt.xlabel("Time (hours)")
plt.ylabel("Number of molecules")
plt.title("Deterministic simulation – Activator A")
plt.legend(loc="upper right")
plt.show()

# Repressor R (Fig. 2b)
plt.figure(figsize=(8, 4))
plt.plot(sol.t, sol.y[7], linestyle="solid", color="red", label="R (repressor)")
plt.xlabel("Time (hours)")
plt.ylabel("Number of molecules")
plt.title("Deterministic simulation – Repressor R")
plt.legend(loc="upper right")
plt.show()



