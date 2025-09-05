"""
Assignment B – Stochastic Model of Circadian Clock
--------------------------------------------------
Reproduce Fig. 2c and Fig. 2d from Vilar et al. (PNAS, 2002) by simulating
the circadian clock using the Stochastic Simulation Algorithm (SSA).

Goal:
- Simulate the dynamics of A (activator) and R (repressor) over 400 hours.
- Use the same parameters as in the deterministic model (Assignment A).
- Observe variability across multiple runs and compare with the ODE solution.
"""

import numpy as np
import matplotlib.pyplot as plt


# -------------------------------------------------------------------
# Random number generators
# -------------------------------------------------------------------
def RandExp(lam, N):
    """Generate N random numbers from an exponential distribution with rate lam."""
    U = np.random.rand(N)
    return -1 / lam * np.log(1 - U)


def RandDisct(x, p, N):
    """Sample N random values from discrete distribution defined by states x and probabilities p."""
    cdf = np.cumsum(p)
    U = np.random.rand(N)
    idx = np.searchsorted(cdf, U)
    return x[idx]


# -------------------------------------------------------------------
# Stochastic Simulation Algorithm (SSA)
# -------------------------------------------------------------------
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
        Time points for each event.
    AllStates : dict
        States at each event time.
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

        which = RandDisct(ReactNum, w / a, 1)     # which reaction occurs
        State = State + StateChangeMat[which.item(), :]  # update state
        k += 1
        AllTimes[k] = t
        AllStates[k] = State

    return AllTimes, AllStates


# -------------------------------------------------------------------
# Parameters (same as Assignment A)
# -------------------------------------------------------------------
alpha_A, alpha_Ap, alpha_R, alpha_Rp = 50, 500, 0.01, 50
beta_A, beta_R = 50, 5
gamma_A, gamma_R, gamma_C = 1, 1, 2
delta_A, delta_R, delta_MA, delta_MR = 1, 0.2, 10, 0.5
theta_A, theta_R = 50, 100

# Initial conditions
Initial = [0, 0, 0, 1, 0, 1, 0, 0, 0]

FinalTime = 400  # hours


# -------------------------------------------------------------------
# State-change matrix
# -------------------------------------------------------------------
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


# -------------------------------------------------------------------
# Propensity function
# -------------------------------------------------------------------
def PropensityFunc(State, ReactNo):
    """Compute propensities for all reactions."""
    w = np.zeros(ReactNo)
    A, R, C, D_A, D_Ap, D_R, D_Rp, M_A, M_R = State

    w[0]  = gamma_C * A * R
    w[1]  = delta_A * A
    w[2]  = delta_A * C
    w[3]  = delta_R * R
    w[4]  = gamma_A * D_A * A
    w[5]  = gamma_R * D_R * A
    w[6]  = theta_R * D_Ap
    w[7]  = alpha_A * D_A
    w[8]  = alpha_Ap * D_Ap
    w[9]  = delta_MA * M_A
    w[10] = beta_A * M_A
    w[11] = theta_R * D_Rp
    w[12] = alpha_R * D_R
    w[13] = alpha_Rp * D_Rp
    w[14] = delta_MR * M_R
    w[15] = beta_R * M_R

    return w


# -------------------------------------------------------------------
# Run SSA simulation
# -------------------------------------------------------------------
Times, States = SSA(Initial, StateChangeMat, FinalTime)

# Convert dicts to lists for plotting
n = len(Times)
t = [Times[i][0] if isinstance(Times[i], (list, np.ndarray)) else Times[i] for i in range(n)]
A_vals = [States[i][0] for i in range(n)]
R_vals = [States[i][1] for i in range(n)]


# -------------------------------------------------------------------
# Plot results
# -------------------------------------------------------------------
plt.figure(figsize=(8, 4))
plt.plot(t, A_vals, linestyle="-", color="blue", label="A (activator)")
plt.xlabel("Time (hours)")
plt.ylabel("Number of molecules")
plt.title("Stochastic simulation – Activator A")
plt.legend(loc="upper right")
plt.show()

plt.figure(figsize=(8, 4))
plt.plot(t, R_vals, linestyle="-", color="red", label="R (repressor)")
plt.xlabel("Time (hours)")
plt.ylabel("Number of molecules")
plt.title("Stochastic simulation – Repressor R")
plt.legend(loc="upper right")
plt.show()
