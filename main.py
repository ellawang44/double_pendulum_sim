import scipy.integrate
import scipy
import numpy as np
import matplotlib.pyplot as plt
import argparse
import init
import double_pendulum
import lyapunov
import time
# typing LaTeX code in matplotlib graphs
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

# parsers to change inputs
parser = argparse.ArgumentParser(description='Double pendulum chaotic system.')

# double pendulum
# initial conditions and parameters
parser.add_argument('--init_cond', metavar = 'initial_condition', type = str, default = '0 0 0 0', help='Change the initial conditions of the system, default value is 0 for all angles and angular velocities. Inputs are strings, with spacebars between each variable. Variables are positioned as: "theta_1 \dot{theta}_1 theta_2 \dot{theta}_2". Accepts "n*pi" as inputs as well, where n is some float. E.g: "6 8.1 9.2 1.5*pi".')
parser.add_argument('--m1', metavar = 'm1', type = float, default = 1, help='Change the mass of the inner pendulum.')
parser.add_argument('--m2', metavar = 'm2', type = float, default = 1, help='Change the mass of the outer pendulum.')
parser.add_argument('--l1', metavar = 'l1', type = float, default = 1, help='Change the length of the inner pendulum.')
parser.add_argument('--l2', metavar = 'l2', type = float, default = 1, help='Change the length of the outer pendulum.')
parser.add_argument('--I1', metavar = 'I1', type = float, default = 1, help='Change the angular momentum of the inner pendulum.')
parser.add_argument('--I2', metavar = 'I2', type = float, default = 1, help='Change the angular momentum of the outer pendulum.')
# ODE solver details
parser.add_argument('--t_min', metavar = 't_min', type = float, default = 0, help='Change the start integration time.')
parser.add_argument('--t_max', metavar = 't_max', type = float, default = 100, help='Change the end integration time.')
parser.add_argument('-n', '--number', metavar = 'number', type = float, default = 1000, help='Change the number of points at which the integral is evaluated between the start and end integration time')

# lyapunov exponent
parser.add_argument('--lyapunov', action = 'store_true', default = False, help='Calculates and prints the lyapunov exponent to the screen.')
parser.add_argument('--lyapunov_max_t', metavar = 'lyapunov_max_t', type = float, default = 100, help='Change the end time used when calculating the lyapunov exponent.')
parser.add_argument('--lyapunov_t_step', metavar = 'lyapunov_t_step', type = float, default = 0.01, help='Change the time step used in the lyapunov exponent calculation.')
parser.add_argument('--lyapunov_length', metavar = 'lyapunov_length', type = float, default = 1e-8, help='Change the d0 length used in the lyapunov exponent calculation.')

# lyapunov exponent values for different variables
parser.add_argument('--stability', action = 'store_true', default = False, help='runs the code used to calculate the lyapunov exponent for varying variables to determine how variables affect the stability of the system. If run with the argument "--xaxis", the lyapunov exponent will be calaculated for different values of that variable with the results saved in a csv file. If run with the arguments "--twoD", "--xaxis", and "--yaxis", the lyapunov exponent will be calculated for two varying variables. You can still set initial conditions or change the parameters using the designated flags to change those arguments.')
parser.add_argument('--xaxis', metavar = 'xaxis', type = str, default = 'theta1', help='Changes the variable changed and plotted on the x-axis.')
parser.add_argument('--x_range', metavar = 'x_range', type = str, default = '-1*pi 1*pi', help='Changes the sampling range of the x variables. Accepts "n*pi" as inputs - read the help documentation for "--init_cond" for more information.')
parser.add_argument('--x_num', metavar = 'x_num', type = int, default = 10, help='Changes the number of x variables sampled.')
parser.add_argument('--twoD', action = 'store_true', default = False, help='Calculates the lyapunov exponent over 2 changing variables.')
parser.add_argument('--yaxis', metavar = 'yaxis', type = str, default = 'theta2', help='Changes the variable changed and plotted on the y-axis.')
parser.add_argument('--y_range', metavar = 'y_range', type = str, default = '-1*pi 1*pi', help='Changes the sampling range of the y variables. Accepts "n*pi" as inputs - read the help documentation for "--init_cond" for more information.')
parser.add_argument('--y_num', metavar = 'y_num', type = int, default = 10, help='Changes the number of y variables sampled.')

# plots
parser.add_argument('--time_angles', action = 'store_true', default = False, help='plots the angles against time.')
parser.add_argument('--cartesian', action = 'store_true', default = False, help='plots the cartesian frame of reference for the motion of the double pendulum.')
parser.add_argument('--phase_space', action = 'store_true', default = False, help='plots the phase space frame of reference for the motion of the double pendulum.')
parser.add_argument('--lyapunov_exponent', action = 'store_true', default = False, help='plots the running mean of the lyapunov exponent')

if __name__ == '__main__':
    # set parameter based upon initial args
    args = parser.parse_args()
    init.m1 = args.m1
    init.m2 = args.m2
    init.l1 = args.l1
    init.l2 = args.l2
    init.I1 = args.I1
    init.I2 = args.I2
    init_cond = []
    for i in args.init_cond.split():
        if i[-2:] == 'pi':
            factor = float(i[:-3])
            init_cond.append(factor*np.pi)
        else:
            init_cond.append(float(i))
    init.init_cond = np.array(init_cond)

    if not args.stability:
        # times and initial conditions
        t_min = args.t_min
        t_max = args.t_max
        eval_points = np.linspace(t_min, t_max, args.number)

        # evaluate the system of odes
        xvals, phi1, phi2, phi3, phi4 = double_pendulum.solve_double_pendulum(init_cond = init.init_cond, t_min= t_min, t_max = t_max, eval_points = eval_points)

        # plots
        if args.time_angles: # plots time against angles
            plt.plot(xvals, phi1, label = r'$\theta_1$')
            plt.plot(xvals, phi3, label = r'$\theta_2$')
            plt.xlabel('time (s)')
            plt.ylabel(r'$\theta$')
            plt.legend()
            plt.savefig('time_angles.png', bbox = 'tight')
            plt.close()

        if args.phase_space: # plots the phase space
            plt.plot(phi1, phi3)
            plt.xlabel(r'$\theta_1$')
            plt.ylabel(r'$\theta_3$')
            plt.savefig('phase_space.png', bbox = 'tight')
            plt.close()

        if args.cartesian: # plots the pendulum in Cartesian coords
            x1 = init.l1*np.sin(phi1)
            y1 = -init.l1*np.cos(phi1)
            x2 = init.l1*np.sin(phi1) + init.l2*np.sin(phi3)
            y2 = -init.l1*np.cos(phi1) - init.l2*np.cos(phi3)
            plt.plot(x1, y1, color = 'aqua', label = 'inner motion')
            plt.plot(x2, y2, color = 'b', label = 'total motion')
            plt.xlabel('x')
            plt.ylabel('y')
            plt.legend()
            plt.savefig('cartesian.png', bbox = 'tight')
            plt.close()

        if args.lyapunov: # calculates the lyapunov exponent and prints it to screen
            start = time.time()
            lyapunov_max_t = args.lyapunov_max_t
            t_step = args.lyapunov_t_step
            lyaps, lyapunov_means = lyapunov.lyapunov_exp(p0 = init.init_cond, max_t = lyapunov_max_t, t_step = t_step)

            # print the lyapunov exponent and error
            print('lyapunov exponent', lyapunov.lyapunov_mean(lyaps,N = (lyapunov_max_t/t_step)))
            print('error', lyapunov.lyapunov_error(lyaps, N = (lyapunov_max_t/t_step)))
            end = time.time()
            print(end - start)

        if args.lyapunov_exponent: # plot running mean of lyapunov exponent
            plt.plot(np.arange(0, lyapunov_max_t, t_step), lyapunov_means)
            plt.xlabel('time (s)')
            plt.ylabel(r'running mean $\lambda$')
            plt.savefig('lyapunov_running_mean.png', bbox = 'tight')
            plt.close()

    if args.stability:
        # get the different x values that we need to calculate the lyapunov exponent for
        x_range = []
        for i in args.x_range.split():
            if i[-2:] == 'pi':
                factor = float(i[:-3])
                x_range.append(factor*np.pi)
            else:
                x_range.append(float(i))
        x_range = np.array(x_range)
        xvals = np.linspace(x_range[0], x_range[1], args.x_num)

        # if the change is in the initial conditions
        # generate initial conditions
        if args.xaxis == 'theta1':
            x_init_conds = [np.array([i, init_cond[1], init_cond[2], init_cond[3]]) for i in xvals]
        elif args.xaxis == 'theta1d':
            x_init_conds = [np.array([init_cond[0], i, init_cond[2], init_cond[3]]) for i in xvals]
        elif args.xaxis == 'theta2':
            x_init_conds = [np.array([init_cond[0], init_cond[1], i, init_cond[3]]) for i in xvals]
        elif args.xaxis == 'theta2d':
            x_init_conds = [np.array([init_cond[0], init_cond[1], init_cond[2], i]) for i in xvals]
        # if the change is not in the initial conditions, then set the right values to change in the init file
        elif args.xaxis == 'm1':
            xval = init.m1
        elif args.xaxis == 'm2':
            xval = init.m2
        elif args.xaxis == 'l1':
            xval = init.l1
        elif args.xaxis == 'l2':
            xval = init.l2
        elif args.xaxis == 'I1':
            xval = init.I1
        elif args.xaxis == 'I2':
            xval = init.I2

        # if the twD flag is active, then get the different y values that we need to calculate the lyapunov exponent for
        if args.twoD:
            y_range = []
            for i in args.y_range.split():
                if i[-2:] == 'pi':
                    factor = float(i[:-3])
                    y_range.append(factor*np.pi)
                else:
                    y_range.append(float(i))
            y_range = np.array(y_range)
            yvals = np.linspace(y_range[0], y_range[1], args.y_num)

            # if the change is in the initial conditions
            # generate initial conditions
            if args.yaxis == 'theta1':
                y_init_conds = [np.array([i, init_cond[1], init_cond[2], init_cond[3]]) for i in yvals]
            elif args.yaxis == 'theta1d':
                y_init_conds = [np.array([init_cond[0], i, init_cond[2], init_cond[3]]) for i in yvals]
            elif args.yaxis == 'theta2':
                y_init_conds = [np.array([init_cond[0], init_cond[1], i, init_cond[3]]) for i in yvals]
            elif args.xaxis == 'theta2d':
                y_init_conds = [np.array([init_cond[0], init_cond[1], init_cond[2], i]) for i in yvals]
            # if the change is not in the initial conditions, then set the right values to change in the init file
            elif args.yaxis == 'm1':
                yval = init.m1
            elif args.yaxis == 'm2':
                yval = init.m2
            elif args.yaxis == 'l1':
                yval = init.l1
            elif args.yaxis == 'l2':
                yval = init.l2
            elif args.yaxis == 'I1':
                yval = init.I1
            elif args.yaxis == 'I2':
                yval = init.I2

        # set lyapunov calculation values
        lyapunov_max_t = args.lyapunov_max_t
        t_step = args.lyapunov_t_step
        N = (lyapunov_max_t/t_step)

        if not args.twoD: # whilst I know there is a lot of code I can reuse between 1D and 2D, I can't be bothered. This is most definately not the elegant way to do things
            # calculate all the lyapunov values
            if args.xaxis == 'theta1' or args.xaxis == 'theta1d' or args.xaxis == 'theta2' or args.xaxis == 'theta2d':
                lyap_vals, lyap_err = ([], [])
                for i in x_init_conds:
                    init.init_cond = i
                    lyaps, _ = lyapunov.lyapunov_exp(p0 = init.init_cond, max_t = lyapunov_max_t, t_step = t_step)
                    err = lyapunov.lyapunov_error(lyaps, N)
                    l = lyapunov.lyapunov_mean(lyaps, N)
                    lyap_vals.append(l)
                    lyap_err.append(err)
            else:
                lyap_vals, lyap_err = ([], [])
                for i in xvals:
                    xval = i
                    lyaps, _ = lyapunov.lyapunov_exp(p0 = init.init_cond, max_t = lyapunov_max_t, t_step = t_step)
                    err = lyapunov.lyapunov_error(lyaps, N)
                    l = lyapunov.lyapunov_mean(lyaps, N)
                    lyap_vals.append(l)
                    lyap_err.append(err)
            x_lyap = np.array([xvals, lyap_vals])
            np.savetxt(fname = 'lyapunov_'+args.xaxis+'.csv', X = x_lyap.T, delimiter = ',')
            x_lyap_err = np.array([xvals, lyap_err])
            np.savetxt(fname = 'lyapunov_' + args.xaxis + '_err.csv', X = x_lyap_err.T, delimiter = ',')

        if args.twoD:
            # accumulate all the lyapunov values
            lyap_vals, lyap_err = ([], [])
            # there's 4 permutations of paramters/variables
            if (args.xaxis == 'theta1' or args.xaxis == 'theta1d' or args.xaxis == 'theta2' or args.xaxis == 'theta2d') and (args.yaxis == 'theta1' or args.yaxis == 'theta1d' or args.yaxis == 'theta2' or args.yaxis == 'theta2d'):
                if args.yaxis == 'theta1':
                    var = 0
                elif args.yaxis == 'theta1d':
                    var = 1
                elif args.yaxis == 'theta2':
                    var = 2
                elif args.yaxis == 'theta2d':
                    var = 3
                for j in yvals:
                    lyap_x_vals, lyap_x_err = ([], [])
                    for i in x_init_conds:
                        np.put(i, var, j)
                        init.init_cond = i
                        lyaps, _ = lyapunov.lyapunov_exp(p0 = init.init_cond, max_t = lyapunov_max_t, t_step = t_step)
                        err = lyapunov.lyapunov_error(lyaps, N)
                        l = lyapunov.lyapunov_mean(lyaps, N)
                        lyap_x_vals.append(l)
                        lyap_x_err.append(err)
                    lyap_vals.append(lyap_x_vals)
                    lyap_err.append(lyap_x_err)
            elif (args.xaxis == 'theta1' or args.xaxis == 'theta1d' or args.xaxis == 'theta2' or args.xaxis == 'theta2d') and (args.yaxis == 'm1' or args.yaxis == 'm2' or args.yaxis == 'l1' or args.yaxis == 'l2' or args.yaxis == 'I1' or args.yaxis == 'I2'):
                for j in yvals:
                    yval = j
                    lyap_x_vals, lyap_x_err = ([], [])
                    for i in xvals:
                        init.init_cond = i
                        lyaps, _ = lyapunov.lyapunov_exp(p0 = init.init_cond, max_t = lyapunov_max_t, t_step = t_step)
                        err = lyapunov.lyapunov_error(lyaps, N)
                        l = lyapunov.lyapunov_mean(lyaps, N)
                        lyap_x_vals.append(l)
                        lyap_x_err.append(err)
                    lyap_vals.append(lyap_x_vals)
                    lyap_err.append(lyap_x_err)
            elif ((args.xaxis == 'm1' or args.xaxis == 'm2' or args.xaxis == 'l1' or args.xaxis == 'l2' or args.xaxis == 'I1' or args.xaxis == 'I2') and (args.yaxis == 'theta1' or args.yaxis == 'theta1d' or args.yaxis == 'theta2' or args.yaxis == 'theta2d')):
                for j in y_init_conds:
                    init.init_cond = j
                    lyap_x_vals, lyap_x_err = ([], [])
                    for i in xvals:
                        xval = i
                        lyaps, _ = lyapunov.lyapunov_exp(p0 = init.init_cond, max_t = lyapunov_max_t, t_step = t_step)
                        err = lyapunov.lyapunov_error(lyaps, N)
                        l = lyapunov.lyapunov_mean(lyaps, N)
                        lyap_x_vals.append(l)
                        lyap_x_err.append(err)
                    lyap_vals.append(lyap_x_vals)
                    lyap_err.append(lyap_x_err)
            else:
                for j in yvals:
                    yval = j
                    lyap_x_vals, lyap_x_err = ([], [])
                    for i in xvals:
                        xval = i
                        lyaps, _ = lyapunov.lyapunov_exp(p0 = init.init_cond, max_t = lyapunov_max_t, t_step = t_step)
                        err = lyapunov.lyapunov_error(lyaps, N)
                        l = lyapunov.lyapunov_mean(lyaps, N)
                        lyap_x_vals.append(l)
                        lyap_x_err.append(err)
                    lyap_vals.append(lyap_x_vals)
                    lyap_err.append(lyap_x_err)
            # print the lyapunov values into a file
            x_lyap_vals = np.append([xvals], lyap_vals, axis = 0)
            nan_yvals = np.append([np.nan], yvals)
            x_y_lyap_vals = np.append(np.array([nan_yvals]).T, x_lyap_vals, axis = 1)
            np.savetxt(fname = 'lyapunov_' + args.xaxis + '_'+args.yaxis + '.csv', X = x_y_lyap_vals, delimiter = ',')
            x_lyap_err = np.append([xvals], lyap_err, axis = 0)
            x_y_lyap_err = np.append(np.array([nan_yvals]).T, x_lyap_err, axis = 1)
            np.savetxt(fname = 'lyapunov_' + args.xaxis + '_' + args.yaxis + '_err.csv', X = x_y_lyap_err, delimiter = ',')
