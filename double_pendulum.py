import numpy as np
import scipy.integrate
import init

def double_pendulum(t, y, M0 = init.M0, L0 = init.L0, I01 = init.I01, I02 = init.I02):
    """Define equations governing the double pendulum system"""
    phi1, phi2, phi3, phi4 = y
    dydt = [phi2, -(np.cos(phi1-phi3)*(np.sin(phi3)-L0*np.sin(phi1-phi3)*phi2**2)/L0 - (I02 + 1/L0) * (M0*np.sin(phi1)+np.sin(phi1-phi3)*phi4**2))/(-(I02/L0)*(I01+L0*M0)+np.cos(phi1-phi3)**2), phi4, -(L0*M0*np.sin(2*phi1-phi3)+(2*I01+3*L0*M0)*np.sin(phi3) - 2*L0*(I01 + L0*M0)*np.sin(phi1 - phi3)*phi2**2 + L0*np.sin(2*(phi1 - phi3))*phi4**2)/(L0 + 2*(1 + I02*L0)*(I01 + L0*M0)+L0*np.cos(2*(phi1 - phi3)))]
    return dydt

def solve_double_pendulum(init_cond, t_min, t_max, eval_points):
    """Numerically solve the ode"""
    sol = scipy.integrate.solve_ivp(fun = double_pendulum, t_span = (t_min, t_max), y0 = init_cond, t_eval = eval_points)
    return (sol.t, sol.y[0], sol.y[1], sol.y[2], sol.y[3])
