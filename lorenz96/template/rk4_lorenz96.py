import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import shutil
import copy
from mpl_toolkits.mplot3d import Axes3D
import pyemu


def lorenz96(x, t, N, F):
    # The Lorenz 1996 model (Lorenz E., 1996. ``Predictability: a problem partly solved'')

    dxdt = np.zeros(N)
    # edges first
    dxdt[0] = (x[1] - x[N-2]) * x[N-1] - x[0]
    dxdt[1] = (x[2] - x[N-1]) * x[0] - x[1]
    dxdt[N-1] = (x[0] - x[N-3]) * x[N-2] - x[N-1]
    for i in range(2, N-1):  # now the rest
        dxdt[i] = (x[i+1] - x[i-2]) * x[i-1] - x[i]
    dxdt = dxdt + F

    # Alternatively can use np.roll
    #dxdt = (-np.roll(x, 1) * (np.roll(x, 2) - np.roll(x, -1)) - x + F)

    return dxdt

def rK4(xt, h, N, F):
    xt = xt.copy()
    k1 = h * lorenz96(xt, h, N, F)
    k2 = h * lorenz96(xt + k1 / 2, h, N, F)
    k3 = h * lorenz96(xt + k2 / 2, h, N, F)
    k4 = h * lorenz96(xt + k3, h, N, F)
    xt += 1 / 6 * (k1 + 2 * (k2 + k3) + k4)

    return xt

def run_L96(N, h, tt, homegrown, random_initial=False):
    '''
    homegrown : True means to use homegrown fourth-order Runge-Kutta solution;
        False means to use scipy's odeint (which uses RK4)
    random_initial : True means to initialize Lorenz state vector randomly;
        False means to initialize as uniform with some perturbation.
    '''

    #wd = os.path.join("lorenz96", "template")

    # treat forcing as 'parameter'
    in_f = "l96_in.csv"
    #if not os.path.isfile(in_f):
    #    in_f = os.path.join(wd, in_f)  # for debugging locally only
    F = pd.read_csv(os.path.join(in_f))['F'][0]
    assert isinstance(F, float)
    print("atmospheric forcing par value: {}".format(F))

    t = np.arange(0.0, tt * h, h)
    x = np.zeros(shape=(t.shape[0], N))
    if random_initial:
        np.random.seed(123)
        x[0, :] = np.random.rand(N)
    else:
        x[0, :] = F * np.ones(N)
        x[0, 19] += 0.01  # add perturbation

    out_f = "l96_out.csv"
    #if not os.path.isfile(out_f):
    #    out_f = os.path.join(wd, out_f)  # for debugging locally only
    with open(os.path.join(out_f), 'w') as f:
        for col in range(N):
            f.write(",x_{}".format(col))
    if homegrown:
        for ti in range(1, x.shape[0]):  # t:
            x[ti, :] = rK4(xt=x[ti - 1], h=h, N=N, F=F)
            with open(os.path.join(out_f), 'a') as f:
                f.write("\n{}".format(ti - 1))
                [f.write(",{}".format(a)) for a in x[ti, :]]
    else:
        from scipy.integrate import odeint
        x = odeint(lorenz96, x[0], t, args=(N, F))

    print("simulation complete")

    return x

def L96(homegrown=True):
    N = 40  # dims
    #F = 2.0  # forcing  # handled within pst
    h = 0.05  # dt
    tt = 1000  # total time

    # run
    x = run_L96(N=N, h=h, tt=tt, homegrown=homegrown)

    # process results
    #plot_hovmoeller(x, h=h, tt=tt, homegrown=homegrown)
    #plot_time_series(x, tt=tt, homegrown=homegrown)
    #plot_profiles(x, tt=tt, homegrown=homegrown)

def plot_time_series(x, tt, homegrown):
    fig, ax = plt.subplots()
    start_time, end_time = tt - 100, tt
    start_loc, end_loc = 5, 8
    lines = ax.plot(np.arange(start_time, end_time),
                    x[start_time:end_time, start_loc:end_loc])
    ax.legend(lines, list(np.arange(start_loc, end_loc)), loc='upper right')
    ax.set_ylabel("$x_{0},..,x_{1}$".format(start_loc, end_loc - 1))
    ax.set_xlabel("time step index")
    #plt.savefig("lorenz_ts_homegrown{}.pdf".format(homegrown))

def plot_profiles(x, tt, homegrown):
    fig, ax = plt.subplots()
    time = tt - 1  #start_time, end_time = 4000, 4100
    lines = ax.plot(np.arange(x.shape[1]), x[time, :])
    #ax.legend(lines, list(np.arange(1)), loc='upper right')
    ax.set_ylabel("$x$ at time step {0}".format(time))
    ax.set_xlabel("site along lat circle")
    #plt.savefig("lorenz_profile_homegrown{}.pdf".format(homegrown))

def plot_hovmoeller(x, h, tt, homegrown, plot_all_times=True):
    fig1, ax2 = plt.subplots()
    if plot_all_times:
        assert tt == x.shape[0]
        t_range = np.arange(0, x.shape[0] * h, h)
        xxx, yyy = np.meshgrid(np.arange(x.shape[1]), t_range)
        CS = ax2.contourf(xxx, yyy, x[:])
    else:
        last_few = 200
        t_range = np.arange((x.shape[0] - last_few) * h, x.shape[0] * h, h)
        xxx, yyy = np.meshgrid(np.arange(x.shape[1]), t_range)
        CS = ax2.contourf(xxx, yyy, x[-t_range.shape[0]:])
    c = plt.colorbar(CS)
    ax2.set_xlabel("equidistant sites along latitude circle", fontsize=12)
    if plot_all_times:
        ax2.set_ylabel("time", fontsize=12)
    else:
        ax2.set_ylabel("time step index", fontsize=12)
    #ax2.set_title("$F$ = {}".format(F), fontsize=12)
    #plt.savefig("l96_hovmoeller_homegrown{}.pdf".format(homegrown))

if (__name__ == '__main__'):
    L96()  # with homegrown RK4
