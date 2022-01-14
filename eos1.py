# 存放状态方程函数的地方
#%%
# import numpy as np
import pandas as pd
# import matplotlib.pyplot as plt
from scipy import optimize as op
from scipy import integrate as inte
import math

#%%
# 常数
G = 6.67259e-11
#%%
# 用来定于积分函数的函数
def m(x,rho):
    return 4*math.pi*(x**2)*rho

def f(x,a_T,b_T,c_T):
    return a_T+b_T*x-c_T*x**(-2)

def p(x, rho):
    return G*m(x, rho)*rho/(x**2)
#%%
# 第一个状态方程，分别对应求密度和求压强两种情况
def eos_density1(x,pre,tem,mat):
    
    data = pd.read_csv(mat+'.txt')
    a_T = data.iloc[0,0]
    b_T = data.iloc[0,1]
    c_T = data.iloc[0,2]
    K_0 = data.iloc[0,3]
    a_p = data.iloc[0,4]
    T_0 = data.iloc[0,5]
    K_p_0 = data.iloc[0,6]
    rho_0 = data.iloc[0,7]

    xmid_rho = inte.quad(f,300,tem,args=(a_T,b_T,c_T))[0]

    K_0_T0 = K_0+a_p*(tem-T_0)
    K_p_T0 = K_p_0
    rho_T0 = rho_0*math.exp(xmid_rho)

    return 3/2*K_0_T0*((x/rho_T0)**(7/3)-(x/rho_T0)**(5/3))*(1-3/4*(4-K_p_T0)*((x/rho_T0)**(2/3)-1))-pre

def eos_pressure1(rho,tem,mat):
    
    data = pd.read_csv(mat+'.txt')
    a_T = data.iloc[0,0]
    b_T = data.iloc[0,1]
    c_T = data.iloc[0,2]
    K_0 = data.iloc[0,3]
    a_p = data.iloc[0,4]
    T_0 = data.iloc[0,5]
    K_p_0 = data.iloc[0,6]
    rho_0 = data.iloc[0,7]

    xmid_rho = inte.quad(f,300,tem,args=(a_T,b_T,c_T))[0]

    K_0_T0 = K_0+a_p*(tem-T_0)
    K_p_T0 = K_p_0
    rho_T0 = rho_0*math.exp(xmid_rho)

    pre = 3/2*K_0_T0*((rho/rho_T0)**(7/3)-(rho/rho_T0)**(5/3))*(1-3/4*(4-K_p_T0)*((rho/rho_T0)**(2/3)-1))
    return pre

#%%
# 第二个状态方程，分别对应求密度和求压强两种情况

#%%
# 计算过程中需要使用的包装函数

def getdensity(pres,temp,mate,rho_initial):
    density_i = op.fsolve(eos_density1, rho_initial, args=(pres,temp,mate))
    return density_i

def getpressure(dens,temp,mate):
    pressure_i = eos_pressure1(dens,temp, mate)
    return pressure_i

def gettemperature(pres_f,pres_b,temp_f,rho_f,rho_b,mate):
    data = pd.read_csv(mate+'.txt')
    # a_T = data.iloc[0,0]
    # b_T = data.iloc[0,1]
    # c_T = data.iloc[0,2]
    # K_0 = data.iloc[0,3]
    # a_p = data.iloc[0,4]
    # T_0 = data.iloc[0,5]
    # K_p_0 = data.iloc[0,6]
    rho_0 = data.iloc[0,7]
    gam_0 = data.iloc[0,1]
    q = data.iloc[0,2]

    dpre = pres_b-pres_f
    drho = rho_b-rho_f
    gam = gam_0*(rho_0/rho_f)**(q)
    phi = dpre/drho

    temp = (dpre)*(gam*temp_f)/(rho_f*phi)
    return temp

def compute_mass(rad_f,rad_b,rho):
    mass_i = inte.quad(m, rad_b, rad_f,args=(rho))[0]
    return mass_i

def compute_pressure(rad_f,rad_b,rho):
    pressure_i = inte.quad(p, rad_b, rad_f, args=(rho))[0]
    return pressure_i

# %%
