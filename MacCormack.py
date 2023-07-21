import numpy as np
from setup import *

#3D array of primitive values prims: [[u], [v], [p], [rho], [e], [mu], [k]]

############################################
# Differences
############################################

def df_dx_f(f, dx):
    # f is a 2d array
    res = (f[:,1:] - f[:,:-1])/dx 
    res = np.pad(res, ((0, 0), (0, 1))) # forward difference
    res[:,-1] = (f[:,-1] - f[:,-2])/dx # back diff at the rightmost boundary

    return res

def df_dx_b(f, dx):
    # f is a 2d array
    res = (f[:,1:] - f[:,:-1])/dx 
    res = np.pad(res, ((0, 0), (1, 0))) # back difference
    res[:,0] = (f[:,1] - f[:,0])/dx # for diff at the left boundary

    return res

def df_dy_f(f, dy):
    # f is a 2d array
    res = (f[1:] - f[:-1])/dy 
    res = np.pad(res, ((0, 1), (0, 0))) # f difference
    res[-1] = (f[-1] - f[-2])/dy # b diff at the upper boundary

    return res

def df_dy_b(f, dy):
    # f is a 2d array
    res = (f[1:] - f[:-1])/dy 
    res = np.pad(res, ((1, 0), (0, 0))) # b difference
    res[0] = (f[1] - f[0])/dy # f diff at the lower boundary

    return res

def df_dx_c(f, dx):
    res = (f[:,2:] - f[:,:-2])/(2 * dx)
    res = np.pad(res, ((0, 0), (1, 1)))
    res[:,0] = (f[:,1] - f[:,0])/dx
    res[:,-1] = (f[:,-1] - f[:,-2])/dx
    return res

def df_dy_c(f, dy):
    res = (f[2:] - f[:-2])/(2 * dy)
    res = np.pad(res, ((1, 1), (0, 0)))
    res[0] = (f[1] - f[0])/dy
    res[-1] = (f[-1] - f[-2])/dy
    return res

############################################
def get_U(prims):
    u, v, p, rho, e, mu, k = prims

    V_sq = u**2 + v**2
    Et = rho*(e + V_sq/2)

    U1 = rho
    U2 = rho*u
    U3 = rho*v
    U4 = Et

    U = np.array([U1, U2, U3, U4])

    return U

def get_prims(U):
    rho = U[0]
    u = U[1]/U[0]
    v = U[2]/U[0]
    e = U[3]/U[0] - (u**2 + v**2)/2
    T = e/cv

    p = rho*R*T
    mu = get_mu(T)
    k = mu*cp/Pr  

    prims = np.array([u, v, p, rho, e, mu, k])

    return prims

def get_shearStress(prims, step_case, vec_case, dx, dy):
    # step_case is 0 if predictor, 1 if corrector
    # vec_case is "E" if vector E and "F" if vector "F"
    # Only for internal points so far
    u, v = prims[:2]
    mu = prims[5]
    lmbda = -2/3 * mu
    if step_case == 0 and vec_case == "E":
        du_dx = df_dx_b(u, dx) # back
        dv_dy = df_dy_c(v, dy) # central
        du_dy = df_dy_c(u, dy) # central
        dv_dx = df_dx_b(v, dx) # back

        tau_xx = lmbda*(du_dx + dv_dy) + 2*mu*du_dx
        tau_xy = mu*(du_dy + dv_dx)

        return tau_xx, tau_xy

    elif step_case == 0 and vec_case == "F":
        du_dx = df_dx_c(u, dx)
        dv_dy = df_dy_b(v, dy)
        du_dy = df_dy_b(u, dy)
        dv_dx = df_dx_c(v, dx)

        tau_yy = lmbda*(du_dx + dv_dy) + 2*mu*dv_dy
        tau_xy = mu*(du_dy + dv_dx)

        return tau_yy, tau_xy
    
    elif step_case == 1 and vec_case == "E":
        du_dx = df_dx_f(u, dx)
        dv_dy = df_dy_c(v, dy)
        du_dy = df_dy_c(u, dy)
        dv_dx = df_dx_f(v, dx)

        tau_xx = lmbda*(du_dx + dv_dy) + 2*mu*du_dx
        tau_xy = mu*(du_dy + dv_dx)

        return tau_xx, tau_xy

    elif step_case == 1 and vec_case == "F":
        du_dx = df_dx_c(u, dx)
        dv_dy = df_dy_f(v, dy)
        du_dy = df_dy_f(u, dy)
        dv_dx = df_dx_c(v, dx)

        tau_yy = lmbda*(du_dx + dv_dy) + 2*mu*dv_dy
        tau_xy = mu*(du_dy + dv_dx)

        return tau_yy, tau_xy
    
    else:
        raise ValueError("Shear stress function recieved incorrect step_case or vec_case")

def get_heatFlux(prims, step_case, vec_case, dx, dy, ADIABATIC):
    e, mu, k = prims[4:]
    T = e/cv
    if step_case == 0 and vec_case == "E":
        dT_dx = df_dx_b(T, dx)
        qx = -k * dT_dx
        return qx
    elif step_case == 0 and vec_case == "F":
        dT_dy = df_dy_b(T, dy)
        if ADIABATIC:
            dT_dy[0] *= 0
        qy = -k * dT_dy
        return qy
    elif step_case == 1 and vec_case == "E":
        dT_dx = df_dx_f(T, dx)
        qx = -k * dT_dx
        return qx
    elif step_case == 1 and vec_case == "F":
        dT_dy = df_dy_f(T, dy)
        if ADIABATIC:
            dT_dy[0] *= 0
        qy = -k * dT_dy
        return qy
    
def get_E_arr(prims, step_case, dx, dy): 
    # takes in 3D array of primitive values prims: [[u], [v], [p], [rho], [e], [mu], [k]]
    u, v, p, rho, e, mu, k = prims

    V_sq = u**2 + v**2
    Et = rho*(e + V_sq/2)

    tau_xx, tau_xy = get_shearStress(prims, step_case, "E", dx, dy)
    qx = get_heatFlux(prims, step_case, "E", dx, dy, False)

    E1 = rho * u
    E2 = rho * u**2 + p - tau_xx
    E3 = rho*u*v - tau_xy
    E4 = (Et + p)*u - u*tau_xx - v*tau_xy + qx

    E_arr = np.array([E1, E2, E3, E4])

    return E_arr

def get_F_arr(prims, step_case, dx, dy, ADIABATIC): 
    # takes in 3D array of primitive values prims: [[u], [v], [p], [rho], [e], [mu], [k]]
    u, v, p, rho, e, mu, k = prims

    V_sq = u**2 + v**2
    Et = rho*(e + V_sq/2)

    tau_yy, tau_xy = get_shearStress(prims, step_case, "F", dx, dy)
    qy = get_heatFlux(prims, step_case, "F", dx, dy, ADIABATIC)

    F1 = rho * v
    F2 = rho * u * v - tau_xy
    F3 = rho * v**2 + p - tau_yy
    F4 = (Et + p)*v - u*tau_xy - v*tau_yy + qy

    F_arr = np.array([F1, F2, F3, F4])

    return F_arr

def get_dt(prims, dx, dy, K=0.8):
    u, v, p, rho, e, mu, k = prims
    T = e/cv
    a = np.sqrt(gamma * R * T)
    v_bar = np.max([4/3*mu, (gamma*mu/Pr)])/rho
    dt_cfl = 1/(abs(u)/dx + abs(v)/dy + a*np.sqrt(1/(dx**2) + 1/(dy**2)) + 2*v_bar*(1/(dx**2) + 1/(dy**2)))

    dt = np.min(K * dt_cfl)

    return dt

def get_boundaries(prims, ADIABATIC):
    u, v, p, rho, e, mu, k = prims

    # Leading edge
    u[0, 0] = 0
    v[0, 0] = 0
    p[0, 0] = p_inf
    e[0, 0] = T_inf * cv

    # Inflow
    u[1:,0] = u_inf
    v[:,0] = 0
    p[:,0] = p_inf
    e[:,0] = T_inf * cv

    # Upper boundary
    u[-1] = u_inf
    v[-1] = 0
    p[-1] = p_inf
    e[-1] = T_inf * cv


    # Surface
    u[0] = 0
    v[0] = 0
    p[0,1:] = 2*p[1, 1:] - p[2, 1:]
    e[0, 1:] = Tw * cv
    if ADIABATIC:
        e[0, 1:] = e[1, 1:]

    # Outflow
    u[1:-1, -1] = 2*u[1:-1, -2] - u[1:-1, -3]
    v[1:-1, -1] = 2*v[1:-1, -2] - v[1:-1, -3]
    p[1:-1, -1] = 2*p[1:-1, -2] - p[1:-1, -3]
    e[1:-1, -1] = 2*e[1:-1, -2] - e[1:-1, -3]
    
    T = e/cv
    rho = p/(R*T)
    mu = get_mu(T)
    k = mu*cp/Pr

    new_prims = np.array([u, v, p, rho, e, mu, k])

    return new_prims

def mac_step(prims, dx, dy, ADIABATIC=False):
    U = get_U(prims)
    U_copy = np.copy(U)

    dt = get_dt(prims, dx, dy)

    E_p = get_E_arr(prims, 0, dx, dy)
    F_p = get_F_arr(prims, 0, dx, dy, ADIABATIC)

    dE_p = (E_p[:,:,1:] - E_p[:,:,:-1])/dx
    dF_p = (F_p[:,1:] - F_p[:,:-1])/dy

    dE_p = np.pad(dE_p, ((0, 0), (0, 0), (0, 1))) # pad is ok here, points are not used
    dF_p = np.pad(dF_p, ((0, 0), (0, 1), (0, 0)))

    dU_dt_p = -(dE_p + dF_p)

    U_copy += dU_dt_p * dt
    new_prims = get_prims(U_copy)

    E_c = get_E_arr(new_prims, 1, dx, dy)
    F_c = get_F_arr(new_prims, 1, dx, dy, ADIABATIC)

    dE_c = (E_c[:,:,1:] - E_c[:,:,:-1])/dx
    dF_c = (F_c[:,1:] - F_c[:,:-1])/dy

    dE_c = np.pad(dE_c, ((0, 0), (0, 0), (1, 0)))
    dF_c = np.pad(dF_c, ((0, 0), (1, 0), (0, 0)))

    dU_dt_c = -(dE_c + dF_c)

    dU_dt = 0.5*(dU_dt_p + dU_dt_c)

    U += dU_dt * dt 

    prims = get_prims(U)

    prims = get_boundaries(prims, ADIABATIC)

    return prims
