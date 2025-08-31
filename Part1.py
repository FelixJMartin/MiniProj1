"""

Deterministic) Reproduce the result in Fig 2a, and Fig 2b in the article by solving
the ODE model given in page 1 (Eq. [1]) to see how both proteins A (activator) and R (repressor)
vary over 400 hours. The model parameters and initial values to be used
are given in the caption of Fig. 1.

"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

alpha_a, alpha_ap, alpha_r, alpha_rp = 50, 500, 0.01, 50 # constant rates
beta_a, beta_r = 50, 5                                   # constant rates
gamma_a, gamma_r, gamma_c = 1, 1, 2                      # constant rates
delta_a, delta_r, delta_ma, delta_mr = 1, 0.2, 10, 0.5   # constant rates
phi_a, phi_r = 50, 100   

Initial = [1, 1, 0, 0, 0, 0, 0, 0, 0]   #DA, DR, DAp, DRp, MA, MR, A, R, C
FinalTime = 400                         #Final time of simulation

def molecules(t, y):
    DA, DR, DAp, DRp, MA, MR, A, R, C = y
    
    yprime = np.zeros(9);
    yprime[0] = phi_a*DAp - gamma_a*DA*A #DA
    yprime[1] = phi_r*DRp - gamma_r*DR*A  #DR
    yprime[2] = gamma_a*DA*A - phi_a*DAp  #DAp
    yprime[3] = gamma_r*DR*A - phi_r*DRp #DRp
    yprime[4] = alpha_ap * DAp + alpha_a*DA - delta_ma * MA #MA
    yprime[5] = alpha_rp * DRp + alpha_r*DR - delta_mr * MR #MR
    yprime[6] = beta_a*MA + phi_a*DAp + phi_r*DRp  - A * ( gamma_a*DA + gamma_r * DR + gamma_c * R + delta_a) # dA/dt
    yprime[7] = beta_r*MR - gamma_c*A*R + delta_a*C - delta_r*R  # dR/dt
    yprime[8] = gamma_c*A*R - delta_a*C


    return yprime

teval = np.linspace(0, FinalTime,1000)    
sol = solve_ivp(molecules, [0,FinalTime], Initial, method = 'BDF', t_eval = teval)


# deterministic
plt.figure(figsize = (8, 4))
plt.plot(sol.t,sol.y[6],linestyle = 'solid', color='blue', label = 'A')
plt.xlabel('Time'); 
plt.ylabel('Number of molecules')
plt.legend(loc='upper right')
plt.show()

plt.figure(figsize = (8, 4))
plt.plot(sol.t,sol.y[7],linestyle = 'solid', color='blue', label = 'R')
plt.xlabel('Time'); 
plt.ylabel('Number of molecules')
plt.legend(loc='upper right')
plt.show()


Stochastic process. 

def Predator_Prey(t, y):
    F,R = y
    yprime = np.zeros(2);
    yprime[0] = beta*F*R - gamma*F
    yprime[1] = alpha*R - beta*F*R
    return yprime

teval = np.linspace(0, FinalTime,1000)      # fine evaluation time samples
sol = solve_ivp(Predator_Prey, [0,FinalTime], Initial, method = 'BDF', t_eval = teval)

plt.figure(figsize = (6, 3))
plt.plot(sol.t,sol.y[0],linestyle = 'solid', color='blue', label = 'Rabbits')
plt.plot(sol.t,sol.y[1],linestyle = 'solid', color='red', label = 'Foxes')
plt.xlabel('Time'); plt.ylabel('Number of animals')
plt.title('Deterministic solution using BDF')
plt.legend(loc='upper right')
plt.show()
