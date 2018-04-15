import numpy as np
from scipy.integrate import solve_ivp
import double_pendulum

def lyapunov_single(dt, d0, p_orig, p_offset):
    '''calculate the lyapunov exponent for 1 time step'''
    diff = p_orig-p_offset
    d = np.sqrt(np.dot(diff, diff))
    return np.log(d/d0)/dt

def lyapunov_exp(p0, t_step = 0.01, length = 1e-8, max_t = 100):
    '''set up system to calculate lyapunov exponent'''
    t = 0
    p0d_not_scale = np.random.rand(4)
    p0d = length/np.sqrt(np.dot(p0d_not_scale, p0d_not_scale))*p0d_not_scale
    lyapunovs = []
    l_mean = 0
    lyapunov_means = []
    while t < max_t:
        # solve ode for the original line
        _, phi1_orig, phi2_orig, phi3_orig, phi4_orig = double_pendulum.solve_double_pendulum(init_cond = p0, t_min= t, t_max = t+t_step, eval_points = [t + t_step])
        p_orig = np.array([phi1_orig[0], phi2_orig[0], phi3_orig[0], phi4_orig[0]])
        # solve ode for the offset line
        _, phi1_offset, phi2_offset, phi3_offset, phi4_offset = double_pendulum.solve_double_pendulum(init_cond = p0+p0d, t_min= t, t_max = t+t_step, eval_points = [t + t_step])
        p_offset = np.array([phi1_offset[0], phi2_offset[0], phi3_offset[0], phi4_offset[0]])
        # calculate the lyapunov exponent
        l = lyapunov_single(t_step, length, p_orig, p_offset)
        # calculate running mean of lyapunov exponent
        t += t_step
        l_mean = ((t - t_step)/t_step * l_mean + l) / (t / t_step)
        # set values ready for next loop
        p0d_not_scale = p_offset - p_orig
        p0d = length/np.sqrt(np.dot(p0d_not_scale, p0d_not_scale))*p0d_not_scale
        p0 = p_orig
        lyapunovs.append(l)
        lyapunov_means.append(l_mean)
        # if you change d0 or t/delta t slightly it won't change at all
    return (lyapunovs, lyapunov_means)

def lyapunov_mean(lyapunov_exps, N):
    '''calculates the mean lyapunov exponent'''
    return sum(lyapunov_exps)/N

def lyapunov_error(lyapunov_exps, N):
    '''calculates the error in the lyapunov exponent'''
    return 1/np.sqrt(N-1)*np.std(lyapunov_exps)
