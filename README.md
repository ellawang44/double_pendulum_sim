# Description
A numerical simulation of the motion of the double pendulum with arbitrary parameters and initial conditions. The Lyapunov exponent of varying initial conditions and parameters can be calculated to determine whether the system is chaotic or not.

# Running the Program
## Inputs
There are many inputs available that all change different parts of the simulation, they can be loosely grouped into: parameters and initial conditions, numerical integration, single Lyapunov exponent, multiple Lyapunov exponents, and plots. Whilst all of these flags have default values and are optional, I recommend changing the initial conditions using the `--init_cond` flag - you'll have a very boring system otherwise. I also recommend toggling on some output flags, these are the ones that fall under plot and lyapunov exponent (both single and double). If you don't use output flags the code will run but not produce any outputs, which would be rather dull. Some flags do produce either `png` or `csv` output files and will replace any files under the same names. Run `--help` to see detailed documentation on what each flag does.

### Parameters and Initial conditions
Changes the parameters and initial conditions of the double pendulum. `--init_cond`, `--m1`, `--m2`, `--l1`, `--l2`, `--I1`, and `--I2` are the available flags.

### Numerical Integration
Changes the numerical integration parameters in the code. `--t_min`, `--t_max`, and `--number` are the available flags.

### Single Lyapunov Exponent
`--lyapunov` toggles on the calculation of the Lyapunov exponent. `--lyapunov_max_t`, `--lyapunov_t_step`, and `--lyapunov_length` change the parameters used in calculating the lyapunov exponent.

### Multiple Lyapunov Exponents
`--stability` toggles on the calculation of multiple Lyapunov exponents. If `--stability` is toggled on, then `--xaxis`, `--x_range`, and `--x_num` are used to determine the changing parameter/initial condition and how it changes. `--twoD` allows two variables to change as the Lyapunov exponent is calculated (2D). When `--twoD` is toggled on, then `--xaxis`, `--x_range`, `--x_num`, `--yaxis`, `--y_range`, and `--y_num` are used to determine the two parameters/initial conditions that change and how they change as the Lyapunov exponent is calculated.

### Plots
Purely output flags, each flag generates a plot. `--time_angles`, `--cartesian`, `--phase_space`, and `--lyapunov_exponent` are the available flags.

## Sample Inputs
The following are some example runs of the simulation.  

Produce the cartesian plot of a system with initial conditions theta1 = 0.5*pi, d theta1 / dt = 0.001, and theta2 = 0, d theta2 / dt = -2*pi. All parameters are 1.
```
python main.py --init_cond "0.5*pi 0.01 0 -2*pi" --cartesian  
```

Produces the cartesian, phase space and angles against time plots of a system with initial conditions theta1 = 0.001, d theta1 / dt = 0, theta2 = -0.001, and d theta2 / dt = 0. All parameters are 1.
```
python main.py --init_cond "0.001 0 -0.001 0" --time_angles --cartesian --phase_space  
```

Calculates and prints to terminal the Lyapunov exponent for a system with unitary paramters and initial conditions given by: theta1 = 0.1*pi, d theta1 / dt = 0, theta2 = 0.2*pi, and d theta2 / dt = 1. The Lyapunov exponent is calcaulated when the time step reaches 500 (by default time step reaches 100 only).
```
python main.py --init_cond "0.1*pi 0 0.2*pi 1" --lyapunov --lyapunov_max_t 500
```

Produces a csv file output containing the Lyapunov exponents for a system with theta1 values varying from -2*pi to 0, 50 equally spaced theta1 values are sampled. The other initial conditions are set to 0, all parameters are unitary.
```
python main.py --stability --xaxis 'theta1' --x_range '-2*pi 0' --init_cond '0 0 -0.001 0' --x_num 50  
```

Produces a csv file output containing the Lyapunov exponents for a system with varying `m1` and `m2` values. `m1` is sampled 50 times from 5 to 1000. `m2` is sampled 20 times from 1 to 200. The initial conditions of this system are given by: theta1 = 0.001, d theta1 / dt = 0, theta2 = -0.001, and d theta2 / dt = 0. All parameters are 1.
```
python main.py --stability --twoD --xaxis 'm1' --yaxis 'm2' --x_range '5 1000' --y_range '1 200' --init_cond '0.001 0 -0.001 0' --x_num 50 --y_num 20
```

# Dependencies
Package dependencies include: `scipy`, `numpy`, and `matplotlib`.
