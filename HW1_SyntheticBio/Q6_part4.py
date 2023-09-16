def f_1(u, t):
    E, S1, I, EI, ES1, P1, F, S2, FS2, P2 = u
    k1 = 2 * (10**(-3))
    k1n = 1 * (10**(-3))
    k2 = 5 * (10**(-3))
    k3 = 3 * (10**(-3))
    k3n = 1.5 * (10**(-3))
    k4 = 5 * (10**(-3))
    ki = 1 * (10**(-3))
    kin = 7 * (10**(-4))

    dE = -ki*E*I + kin*EI + k1n*ES1 - k1*E*S1 + k4*FS2
    dS1 = k1n*ES1 - k1*E*S1
    dI = -ki*E*I + kin*EI
    dEI = ki*E*I - kin*EI
    dES1 = -k1n*ES1 + k1*E*S1 - k2*ES1
    dP1 = k2*ES1
    dF = -k3*F*S2 + k3n*FS2 + k2*ES1
    dS2 = -k3*F*S2 + k3n*FS2 + k2*ES1
    dFS2 = k3*F*S2 - k3n*FS2 - k4*FS2
    dP2 = k4*FS2

    return [dE,dS1,dI,dEI,dES1,dP1,dF,dS2,dFS2,dP2]

def f_2(u, t):
    E, S1, I, FI, ES1, P1, F, S2, FS2, P2 = u
    k1 = 2 * (10**(-3))
    k1n = 1 * (10**(-3))
    k2 = 5 * (10**(-3))
    k3 = 3 * (10**(-3))
    k3n = 1.5 * (10**(-3))
    k4 = 5 * (10**(-3))
    ki = 1 * (10**(-3))
    kin = 7 * (10**(-4))

    dE = k1n*ES1 - k1*E*S1 + k4*FS2
    dS1 = k1n*ES1 - k1*E*S1
    dI = -ki*F*I + kin*FI
    dFI = ki*F*I - kin*FI
    dES1 = -k1n*ES1 + k1*E*S1 - k2*ES1
    dP1 = k2*ES1
    dF = -k3*F*S2 + k3n*FS2 + k2*ES1 - ki*F*I + kin*FI
    dS2 = -k3*F*S2 + k3n*FS2 + k2*ES1
    dFS2 = k3*F*S2 - k3n*FS2 - k4*FS2
    dP2 = k4*FS2

    return [dE,dS1,dI,dFI,dES1,dP1,dF,dS2,dFS2,dP2]

def f_3(u, t):
    E, S1, I, EI, FI, ES1, P1, F, S2, FS2, P2 = u
    k1 = 2 * (10**(-3))
    k1n = 1 * (10**(-3))
    k2 = 5 * (10**(-3))
    k3 = 3 * (10**(-3))
    k3n = 1.5 * (10**(-3))
    k4 = 5 * (10**(-3))
    ki = 1 * (10**(-3))
    kin = 7 * (10**(-4))
    kj = 1 * (10**(-3))
    kjn = 7 * (10**(-4))

    dE = -ki*E*I + kin*EI + k1n*ES1 - k1*E*S1 + k4*FS2
    dS1 = k1n*ES1 - k1*E*S1
    dI = -kj*F*I + kjn*FI - ki*E*I + kin*EI
    dEI = ki*E*I - kin*EI
    dFI = kj*F*I - kjn*FI
    dES1 = -k1n*ES1 + k1*E*S1 - k2*ES1
    dP1 = k2*ES1
    dF = -k3*F*S2 + k3n*FS2 + k2*ES1 - kj*F*I + kjn*FI
    dS2 = -k3*F*S2 + k3n*FS2 + k2*ES1
    dFS2 = k3*F*S2 - k3n*FS2 - k4*FS2
    dP2 = k4*FS2

    return [dE,dS1,dI,dEI,dFI,dES1,dP1,dF,dS2,dFS2,dP2]


from ODESolver import ForwardEuler
import numpy as np
import matplotlib.pyplot as plt


# Initial condition, start at the origin:
E0 = 10
S10 = 300
I0 = 30
EI0 = 0
FI0 = 0
ES10 =0
P10 = 0
F0 = 0
S20 = 300
FS20 = 0
P20 = 0 

U0_1 = [E0,S10,I0,EI0,ES10,P10,F0,S20,FS20,P20]
U0_2 = [E0,S10,I0,FI0,ES10,P10,F0,S20,FS20,P20]
U0_3 = [E0,S10,I0,EI0,FI0,ES10,P10,F0,S20,FS20,P20]


def plotter(f,U0,No):
    solver = ForwardEuler(f)
    solver.set_initial_condition(U0)
    time_points = np.linspace(0, 30000, 300000)
    u, t = solver.solve(time_points)
    
    p_1 = u[:,5]; x = time_points
    p_2 = u[:,9]
    
    plt.subplot(1, 2, 1)
    plt.plot(x, p_1, color="r")
    plt.xlabel("Time")
    plt.ylabel("P1")
    
    plt.subplot(1, 2, 2)
    plt.plot(x, p_2, color='c')
    plt.xlabel("Time")
    plt.ylabel("P2")
    plot_name = "Reaction NO" + str(No)
    path_name = "D:\Master\Current-Term\Synthetic-biology\HW03\\"
    plt.suptitle(plot_name)
    plt.savefig(path_name + plot_name + ".png")
    plt.close()
    return p_1, p_2

def plotter_comparison(r1,r2,r3):
    time_points = np.linspace(0, 30000, 300000)
    x = time_points
    r1_p_1 = r1[0]; r1_p_2 = r1[1]
    r2_p_1 = r2[0]; r2_p_2 = r2[1]
    r3_p_1 = r3[0]; r3_p_2 = r3[1]
    
    plt.subplot(1, 2, 1)
    plt.plot(x, r1_p_1, color="g")
    plt.plot(x, r2_p_1, color="b")
    plt.plot(x, r3_p_1, color="y")
    plt.xlabel("Time")
    plt.ylabel("P1")
    
    plt.subplot(1, 2, 2)
    plt.plot(x, r1_p_2, color="g")
    plt.plot(x, r2_p_2, color="b")
    plt.plot(x, r3_p_2, color="y")
    plt.xlabel("Time")
    plt.ylabel("P2")
    plot_name = "Comparison"
    path_name = "D:\Master\Current-Term\Synthetic-biology\HW03\\"
    plt.suptitle(plot_name)
    plt.savefig(path_name + plot_name + ".png")
    plt.close()

# Reaction No.1
r1 = plotter(f_1,U0_1,1)

# Reaction No.2
r2 = plotter(f_2,U0_2,2)

# Reaction No.3
r3 = plotter(f_3,U0_3,3)

# Comparison
plotter_comparison(r1,r2,r3)