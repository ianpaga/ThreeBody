import yt
import yt.units as u
import os
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from scipy.optimize import fsolve


def L1(q):
    def f(x):
        return -1 / (x - q / (1 + q))**2 + q / (x + 1 /
                                                (1 + q))**2 - (1 + q) * x

    return (fsolve(f, 0)[0])


def Theta(q):
    mu = q / (1 + q)
    x = L1(q)
    A = mu / abs(-x - 1 + mu)**3 + (1 - mu) / abs(-x + mu)**3
    sin = -(8 / 9 / A)**(1 / 2) * (1 - 2 / A + 3 *
                                   (1 - 8 / 9 / A)**(1 / 2))**(1 / 2)
    theta = np.arcsin(sin) / 2
    return (theta)


def circle(x, y, r, ax=None, label=None, color=None):
    if ax == None:
        f, ax = plt.subplots(figsize=(12, 8))
    theta = np.linspace(0, np.pi * 2, 1000)
    ax.plot(x + r * np.cos(theta),
            y + r * np.sin(theta),
            label=label,
            color=color)


def Radius(m):
    a, b = m / 1.44, m / 0.00057
    return (0.0114 * (a**(-2 / 3) - a**(2 / 3))**(1 / 2) *
            (1 + 3.5 * b**(-2 / 3) + b**(-1))**(-2 / 3))


class WDWD:
    def __init__(self, file_name):
        Rdap = np.loadtxt('./data/R_{}'.format(file_name))
        Jdap = np.loadtxt('./data/J_{}'.format(file_name))
        Para = np.loadtxt('./data/Para_{}'.format(file_name))
        with open(file='./data/Type_{}'.format(file_name)) as f:
            self.type = f.readline()[:2]

        self.md, self.ma, self.mp = Para[0], Para[1], Para[2]
        self.q = self.md / self.ma

        self.P = Para[3] # Orbital period in minute
        self.Omega = Para[4]

        self.x = L1(self.q)
        self.thetas = Theta(self.q)

        self.Radiusd = Radius(self.md)
        self.Radiusa = Radius(self.ma)

        self.a = self.Radiusd / (1 / (1 + self.q) + self.x)
        self.RL = self.a * (1 / (1 + self.q) + self.x)
        self.RL1 = self.a * (self.q / (1 + self.q) - self.x)
        self.Rc = (1 + self.q) * self.RL1**4 / self.a**3

        self.t = Rdap[:, 0] / np.pi / 2
        self.Rd = Rdap[:, 1:3] * self.a
        self.Ra = Rdap[:, 3:5] * self.a
        self.Rp = Rdap[:, 5:7] * self.a

        self.Jdorb, self.Jdsp = Jdap[:, 1], Jdap[:, 2]
        self.Jaorb, self.Jasp = Jdap[:, 3], Jdap[:, 4]
        self.Jp, self.Jtot = Jdap[:, 5], Jdap[:, 6]
        self.Jorb = self.Jdorb + self.Jaorb

    def plot_trajectory(self, ax=None):
        if ax == None:
            f, ax = plt.subplots(figsize=(12, 8))

        ax.plot(self.Rd[:, 0], self.Rd[:, 1], linestyle='-', label='Donor')
        ax.plot(self.Ra[:, 0], self.Ra[:, 1], linestyle='--', label='Accretor')
        ax.plot(self.Rp[:, 0], self.Rp[:, 1], linestyle=':')

        circle(self.Ra[-1, 0], self.Ra[-1, 1], self.Rc, ax=ax,
               label='Circularization Radius')
        circle(self.Ra[-1, 0], self.Ra[-1, 1], self.Radiusa, ax=ax, color='Black')
        circle(self.Rd[0, 0], self.Rd[0, 1], self.Radiusd, ax=ax, color='Black')

        ax.axis('equal')
        ax.set_xlabel('X ($R_{\odot}$)', fontsize=25)
        ax.set_ylabel('Y ($R_{\odot}$)', fontsize=25)
        ax.tick_params(labelsize=25)
        plt.legend(prop={'size': 20})
        plt.tight_layout()
        # plt.show()

    def plot_angular_m(self):
        f, ax = plt.subplots(2, 1, figsize=(12, 16), sharex=True)
        m = self.mp/(self.ma+self.md)

        ax[0].plot(self.t, (self.Jdorb - self.Jdorb[0])/self.Jtot[0]/m, label='Donor orbital')
        ax[0].plot(self.t, (self.Jaorb - self.Jaorb[0])/self.Jtot[0]/m, label='Accretor orbital')
        ax[0].plot(self.t, (self.Jtot - self.Jtot[0])/self.Jtot[0]/m, label='Total', color='k', linestyle='--')
        ax[0].plot(self.t, (self.Jp - self.Jp[0])/self.Jtot[0]/m, label='Particle')

        ax[1].plot(self.t, (self.Jdsp - self.Jdsp[0])/self.Jtot[0]/m, label='Donor spin')
        ax[1].plot(self.t, (self.Jasp - self.Jasp[0])/self.Jtot[0]/m, label='Accretor spin')

        ax[1].set_xlabel('T ($P_{orb}$)', fontsize=25)
        ax[0].set_ylabel(r'$\Delta J_{orb}/m_P$ ($J_{tot}$)', fontsize=25)
        ax[1].set_ylabel(r'$\Delta J_{spin}/m_P$ ($J_{tot}$)', fontsize=25)
        for a in ax:
            a.tick_params(labelsize=25)
            a.legend(prop={'size': 20})
        f.tight_layout()
        # plt.show()

    def plot_delta_J(self, ax=None):
        m = self.mp/(self.ma+self.md)
        if ((self.Jtot[-1]-self.Jtot[0])/self.Jtot[0] > 1e-8) or self.type != 'DI' or self.q>0.8:
            #print(self.type)
            return
        if ax == None:
            f, ax = plt.subplots(figsize=(12, 8))
        ax.scatter(self.md, (self.Jdorb[-1] - self.Jdorb[0])/self.Jtot[0]/m/self.t[-1], color='black', label='Donor orbital', s=2*np.exp(self.q*5))
        ax.scatter(self.md, (self.Jdsp[-1] - self.Jdsp[0])/self.Jtot[0]/m/self.t[-1], color='cyan', label='Donor spin', s=2*np.exp(self.q*5))
        ax.scatter(self.md, (self.Jaorb[-1] - self.Jaorb[0])/self.Jtot[0]/m/self.t[-1], color='blue', label='Accretor orbital', s=2*np.exp(self.q*5))
        ax.scatter(self.md, (self.Jasp[-1] - self.Jasp[0])/self.Jtot[0]/m/self.t[-1], color='pink', label='Accretor spin', s=2*np.exp(self.q*5))
        ax.scatter(self.md, (self.Jorb[-1] - self.Jorb[0])/self.Jtot[0]/m/self.t[-1], color='green', label='Total orbital', s=2*np.exp(self.q*5))
        ax.set_xlabel(r'$M_D$ ($M_{\odot}$)', fontsize=25)
        ax.set_ylabel(r'$\Delta J/m_P$ ($J_{tot}$)', fontsize=25)
