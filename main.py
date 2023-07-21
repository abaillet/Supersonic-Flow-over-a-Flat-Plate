##################################################################
# Navier-Stokes solution of supersonic flow over flat plate.
# Following J.D. Anderson, Computational Fluid Dynamics: The Basics with Applications, Chapter 10
##################################################################

from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt

import setup as setup
from MacCormack import mac_step
import flow_viz as viz

# 3D array of primitive values prims: [[u], [v], [p], [rho], [e], [mu], [k]]

def mass_flow_check(u0, rho0, un, rhon, dy):
    mass_flow0 = np.trapz(u0 * rho0, dx=dy)
    mass_flow_out = np.trapz(un * rhon, dx=dy)

    delta_mass_flow = abs(mass_flow_out - mass_flow0)/mass_flow0 * 100
    print()
    print(f"Mass flow in: {mass_flow0:.5f}")
    print(f"Mass flow out: {mass_flow_out:.5f}")
    return delta_mass_flow

def main(timeSteps, MAX_TOL=1e-8):

    Re_L = setup.Re_L

    print(f"Reynolds number: {Re_L:.3e}")
    print(f"M: {setup.M_inf}")

    prims = setup.prims

    ADIABATIC = setup.ADIABATIC

    print("Adiabatic boundary condition is", ADIABATIC)
    print()
    for i in tqdm(range(timeSteps), desc="Calculating", colour="blue"):
        rho_prev = prims[3]
        prims = mac_step(prims, setup.dx, setup.dy, ADIABATIC)
        rho = prims[3]
        check_arr = np.abs(rho - rho_prev)
        conv_condition = np.all(check_arr < MAX_TOL)
        if conv_condition:
            print(f"Converged in {i} iterations.")
            break
    if i == (timeSteps - 1):
        print(f"Max difference = {np.max(check_arr):.4e}")
       
    u, v, pressure, rho, e, mu, k = prims
    u0 = u[:,0]
    rho0 = rho[:,0]
    un = u[:,-1]
    rhon = rho[:,-1]

    T = e/setup.cv

    sound_sp = np.sqrt(setup.gamma*setup.R*T)
    V = np.sqrt(u**2 + v**2)
    M = V/sound_sp

    X = setup.X
    Y = setup.Y

    dy = setup.dy

    d_mdot = mass_flow_check(u0, rho0, un, rhon, dy) # check conservation of mass

    print(f"Percentage difference between mass flow in and mass flow out: {d_mdot:.3f}%")
    ylim = int(len(Y)//2) + 10

    viz.plot_xprofile([pressure[0]/setup.p_inf], list(range(len(X[0]))), r"$P/P_\infty$", "Pressure along surface of plate, M=4, Adiabatic Wall")
    viz.plot_map(pressure, X, Y, "Pressure [Pa]")
    viz.plot_map(M, X, Y, "M")
    viz.plot_map(u, X, Y, "u [m/s]")
    viz.plot_map(T, X, Y, "T [K]")
    viz.plot_yprofile(T[:,-1]/setup.T_inf, Y[:,0], r"$T/T_\infty$", "Temperature distribution at trailing edge, M=4, Adiabatic Wall")
    # viz.plot_3D(pressure, X, Y, "Pressure [Pa]")
    # viz.plot_3D(M[:ylim], X[:ylim], Y[:ylim], "M")
    # slice_ = 5
    # viz.plot_field(u[::slice_,::slice_], v[::slice_,::slice_], X[::slice_, ::slice_], Y[::slice_, ::slice_], "Velocity field")
    plt.show()

if __name__ == '__main__':

    print("\n##################################################################\n")
    print("Navier-Stokes solution of supersonic flow over flat plate.")
    print("Following J.D. Anderson, Computational Fluid Dynamics, chapter 10\n")
    print("##################################################################\n")

    timeSteps = 10_000 # max number of time steps

    TOL = 1e-6

    main(timeSteps, TOL)

