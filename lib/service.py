# coding: utf-8
# Copyright (C) 2021  Nguyen Ngoc Sang, <https://github.com/SangVn>

from numpy import sin, cos, deg2rad, loadtxt, sqrt
import matplotlib.pyplot as plt

Rgas = 287.052873836    # hằng số chất khí
gamma = 1.4             # số mũ đoạn nhiệt
gamma_m1 = 0.4
cp = Rgas*gamma/(gamma_m1)
Pr = 0.72               # số Prandtl
c_mu = 1.72e-5*(273+122.)/273.**1.5 # hằng số để tính đột nhớt

def inch2m(x_inch):
    return x_inch*0.0254

def P_US2SI(M=1.0, T_rankine=1.0, p_psi=1.0, alf=0.0):
    T_kelvin = T_rankine*5.0/9.0
    p_pascal = p_psi*6894.75728
    return Pm2P(M, T_kelvin, p_pascal, alf)

def PRe_US2SI(M=1.0, T_rankine=1.0, Re=1.0, L=1.0, alf=0.0):
    T_kelvin = T_rankine*5.0/9.0
    mu = muT(T_kelvin)
    p_pascal = Re*mu*sqrt(Rgas*T_kelvin/gamma)/(M*L)
    return Pm2P(M, T_kelvin, p_pascal, alf)

def rho(T, p):
    return p/(Rgas*T)

def Temperature(P):
    return P[3]/(Rgas*P[0])

# vận tốc âm thanh
def VSound_pr(p, rho):
    return (gamma*p/rho)**0.5

def VSound(P):
    return (gamma*P[3]/P[0])**0.5

# hàm tính số mach
def Mach(P):
    u = (P[1]*P[1] + P[2]*P[2])**0.5
    return u/VSound(P)

def muT(T):
    return c_mu*(T**1.5)/(T+122.0)

def muP(P):
    T = Temperature(P)
    return c_mu*(T**1.5)/(T+122.)

def reynolds_number(P, L=1.0):
    r = P[0]
    T = Temperature(P)
    mu = muT(T)
    V = (P[1]*P[1] + P[2]*P[2])**0.5
    Re_L = r*V*L/mu
    return Re_L

# t - total
def Pt2P(M=0.0, Tt=293.15, pt=101325.0, alf=0.0):
    m = 1.0/(1.0+M*M*gamma_m1/2.0)
    T = Tt*m
    p = pt*(m**(gamma/gamma_m1))
    r = rho(T, p)
    V = M*VSound_pr(p,r)
    alf = deg2rad(alf)
    return [r, V*cos(alf), V*sin(alf), p]

def Pm2P(M=0.0, T=293.15, p=101325.0, alf=0.0):
    r = rho(T, p)
    V = M*VSound_pr(p,r)
    alf = deg2rad(alf)
    return [r, V*cos(alf), V*sin(alf), p]

def P2Pt(P):
    r, u, v, p = P[0], P[1], P[2], P[3]
    M = Mach(P)
    m = (1.0+M*M*gamma_m1/2.0)
    T = Temperature(P)
    Tt = T*m
    pt = p*(m**(gamma/gamma_m1))
    rt = r*(m**(1.0/gamma_m1))
    return [rt, pt, Tt]

def plot_state(files, x="iter", y="misclosure", legend=None, save_path=None):
    fd = {"iter": 0, "time": 1, "misclosure": 2}
    for file in files:
        data = loadtxt(file, skiprows=2).T
        plt.plot(data[fd[x]], data[fd[y]])
    if legend is not None:
        plt.legend(legend)
    plt.xlabel(x)
    plt.ylabel(y)
    if save_path is not None:
        plt.savefig(save_path, bbox_inches="tight")
    plt.show()