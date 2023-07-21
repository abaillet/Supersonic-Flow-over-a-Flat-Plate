import numpy as np

##############################################
# Atmospheric properties
##############################################
R = 287 # J/(kg K)
gamma = 1.4
cv = R/(gamma - 1)
cp = gamma * cv
Pr = 0.71

##############################################
# Freestream conditions
##############################################
p_inf = 101_325 # N/m^2
T_inf = 288.16 # K
e_inf = T_inf*cv
a_inf = np.sqrt(gamma*R*T_inf) # m/s
M_inf = 4 # Mach number
u_inf = M_inf * a_inf # m/s
v_inf = 0 # m/s
Tw = T_inf * 1 # temperature of the wall, not technically a freestream condition
rho_inf = p_inf/(R * T_inf) # Density

ADIABATIC = True # Temperture condition on the plate

T0 = 288.16 # Reference values, change if altitude changes
mu0 = 1.7894e-5  # Reference values, change if altitude changes

def get_mu(T):
    mu = mu0*(T/T0)**1.5 * (T0 + 110)/(T + 110)
    return mu

mu_inf = get_mu(T_inf)
k_inf = mu_inf*cp/Pr

##############################################
# Domain sizing
##############################################

LHORI = 1e-5 # m, length of plate

xSteps = 70
ySteps = 70

dx = LHORI/(xSteps - 1)

x_arr = np.linspace(0, LHORI, xSteps)

# Verical height should be 5 times the boundary layer at trailing edge (delta)

Re_L = rho_inf*u_inf*LHORI/mu_inf

delta = 5 * LHORI/(np.sqrt(Re_L))

LVERT = 5 * delta

dy = LVERT/(ySteps - 1)

y_arr = np.linspace(0, LVERT, ySteps)

X, Y = np.meshgrid(x_arr, y_arr)

##############################################
# Initial condition
##############################################

#3D array of primitive values prims: [[u], [v], [p], [rho], [e], [mu], [k]]

# Initial condition
u = np.ones((xSteps, ySteps)) * u_inf
v = np.ones((xSteps, ySteps)) * v_inf
p = np.ones((xSteps, ySteps)) * p_inf
rho = np.ones((xSteps, ySteps)) * rho_inf
e = np.ones((xSteps, ySteps)) * e_inf
mu = np.ones((xSteps, ySteps)) * mu_inf
k = np.ones((xSteps, ySteps)) * k_inf

u[0] = 0 # no slip B.C.
v[0] = 0 # no slip B.C.

prims = np.array([u, v, p, rho, e, mu, k]) # initial array

##############################################