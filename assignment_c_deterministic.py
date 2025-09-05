"""
Assignment C – Deterministic vs. Stochastic Comparison
------------------------------------------------------
Reproduce Fig. 5 from Vilar et al. (PNAS, 2002).

Goal:
- Simulate the circadian clock using both the deterministic ODE model and
  the stochastic SSA model (with modified parameters).
- Compare the dynamics of the repressor protein R over 400 hours.
- Highlight qualitative differences caused by noise in the stochastic model.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


# -------------------------------------------------------------------
# Shared parameters (same for deterministic and stochastic models)
# -------------------------------------------------------------------
alpha_a, alpha_ap, alpha_r, alpha_rp = 50, 500, 0.01, 50
beta_a, beta_r = 50, 5
gamma_a, gamma_r, gamma_c = 1, 1, 2
delta_a, delta_r, delta_ma, delta_mr = 1, 0.05, 10, 0.5   # note: delta_r differs from Assignment A
phi_a, phi_r = 50, 100
theta_a, theta_r = 50, 100

FinalTime = 400  # hours


# -------------------------------------------------------------------
# Deterministic ODE model
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


# Initial conditions (deterministic)
Initial_det = [1, 1, 0, 0, 0, 0, 0, 0, 0]

# Solve ODE
teval = np.linspace(0, FinalTime, 1000)
sol = solve_ivp(molecules, [0, FinalTime], Initial_det, method="BDF", t_eval=teval)


# -------------------------------------------------------------------
# Stochastic SSA model
# -------------------------------------------------------------------
def RandExp(lam, N):
    """Generate N random numbers from an exponential distribution with rate lam."""
    U = np.random.rand(N)
    return -1 / lam * np.log(1 - U)


def RandDisct(x, p, N):
    """Sample N random values from a discrete distribution with states x and probabilities p."""
    cdf = np.cumsum(p)
    U = np.random.rand(N)
    idx = np.searchsorted(cdf, U)
    return x[idx]


def SSA(Initial, StateChangeMat, FinalTime):
    """
    Run SSA simulation.

    Parameters
    ----------
    Initial : list or array
        Initial state vector.
    StateChangeMat : ndarray
        State-change matrix (reactions x states).
    FinalTime : float
        End time of simulation.

    Returns
    -------
    AllTimes : dict
        Event times.
    AllStates : dict
        State vectors at each event.
    """
    m, n = StateChangeMat.shape
    ReactNum = np.arange(m)

    AllTimes, AllStates = {}, {}
    AllStates[0] = Initial
    AllTimes[0] = [0]

    k, t, State = 0, 0, np.array(Initial)

    while True:
        w = PropensityFunc(State, m)        # propensities
        a = np.sum(w)                       # total rate
        tau = RandExp(a, 1)                 # time to next reaction
        t += tau
        if t > FinalTime:
            break

        which = RandDisct(ReactNum, w / a, 1)       # which reaction occurs
        State = State + StateChangeMat[which.item(), :]  # update state
        k += 1
        AllTimes[k] = t
        AllStates[k] = State

    return AllTimes, AllStates


# State-change matrix
StateChangeMat = np.array([
    [-1, -1,  1,  0,  0,  0,  0,  0,  0],
    [-1,  0,  0,  0,  0,  0,  0,  0,  0],
    [ 0,  1, -1,  0,  0,  0,  0,  0,  0],
    [ 0, -1,  0,  0,  0,  0,  0,  0,  0],
    [-1,  0,  0, -1,  1,  0,  0,  0,  0],
    [-1,  0,  0,  0,  0, -1,  1,  0,  0],
    [ 1,  0,  0,  1, -1,  0,  0,  0,  0],
    [ 0,  0,  0,  0,  0,  0,  0,  1,  0],
    [ 0,  0,  0,  0,  0,  0,  0,  1,  0],
    [ 0,  0,  0,  0,  0,  0,  0, -1,  0],
    [ 1,  0,  0,  0,  0,  0,  0,  0,  0],
    [ 1,  0,  0,  0,  0,  1, -1,  0,  0],
    [ 0,  0,  0,  0,  0,  0,  0,  0,  1],
    [ 0,  0,  0,  0,  0,  0,  0,  0,  1],
    [ 0,  0,  0,  0,  0,  0,  0,  0, -1],
    [ 0,  1,  0,  0,  0,  0,  0,  0,  0]
])


def PropensityFunc(State, ReactNo):
    """Compute propensities for all reactions."""
    w = np.zeros(ReactNo)
    A, R, C, D_A, D_Ap, D_R, D_Rp, M_A, M_R = State

    w[0]  = gamma_c * A * R
    w[1]  = delta_a * A
    w[2]  = delta_a * C
    w[3]  = delta_r * R
    w[4]  = gamma_a * D_A * A
    w[5]  = gamma_r * D_R * A
    w[6]  = theta_r * D_Ap
    w[7]  = alpha_a * D_A
    w[8]  = alpha_ap * D_Ap
    w[9]  = delta_ma * M_A
    w[10] = beta_a * M_A
    w[11] = theta_r * D_Rp
    w[12] = alpha_r * D_R
    w[13] = alpha_rp * D_Rp
    w[14] = delta_mr * M_R
    w[15] = beta_r * M_R

    return w


# Initial conditions (stochastic)
Initial_sto = [0, 0, 0, 1, 0, 1, 0, 0, 0]

# Run SSA
Times, States = SSA(Initial_sto, StateChangeMat, FinalTime)

# Convert dicts to lists
n = len(Times)
t_sto = [Times[i][0] if isinstance(Times[i], (list, np.ndarray)) else Times[i] for i in range(n)]
R_vals = [States[i][1] for i in range(n)]


# -------------------------------------------------------------------
# Plot comparison
# -------------------------------------------------------------------
plt.figure(figsize=(10, 5))

# Deterministic R
plt.plot(sol.t, sol.y[7], linestyle="-", color="black", label="Deterministic R")

# Stochastic R
plt.plot(t_sto, R_vals, linestyle="--", color="red", label="Stochastic R")

plt.xlabel("Time (hours)")
plt.ylabel("Number of molecules")
plt.title("Assignment C – Deterministic vs Stochastic (Repressor R)")
plt.legend(loc="upper right")
plt.tight_layout()
plt.show()
