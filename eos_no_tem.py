# 存放状态方程函数的地方
#%%
import math

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pandas.io.parsers import read_table
from scipy import integrate as inte
from scipy import optimize as op
from scipy.integrate.odepack import odeint

import parameter

#%%
# 星球参数设置区域

# 常数
G = 6.67259e-11

# 星球相关参数
M_p = parameter.M_p # 行星质量，（0.008个地球质量）kg
# 输入模型参数
M_h2o = parameter.M_h2o # 壳的质量分数
# M_fe = parameter.M_fe # 核的质量分数
# 参数计算
# M_mgfe2sio4 = 1-M_h2o-M_fe # 硅酸盐幔的质量分数
step = parameter.step
#%%
# 用来定于积分函数的函数
def m(x,rho):
    return 4*math.pi*(x**2)*rho

def f(x,a_T,b_T,c_T):
    return a_T+b_T*x-c_T*x**(-2)

def g(x,rho):
    return (x**2)*rho

def p(x, rho,):
    return G*rho*4*math.pi*(x**2)*rho/(x**2)
#%%
# 第一个状态方程，分别对应求密度和求压强两种情况
def eos_density1(x,pre,mat):
    
    data = pd.read_csv(mat+'.txt')

    K_0 = data.iloc[0,0]
    K_p_0 = data.iloc[0,1]
    rho_0 = data.iloc[0,2]

    return 3/2*K_0*((x/rho_0)**(7/3)-(x/rho_0)**(5/3))*(1-3/4*(4-K_p_0)*((x/rho_0)**(2/3)-1))-pre

def eos_pressure1(rho,mat):
    
    data = pd.read_csv(mat+'.txt')

    K_0 = data.iloc[0,0]
    K_p_0 = data.iloc[0,1]
    rho_0 = data.iloc[0,2]

    pre = 3/2*K_0*((rho/rho_0)**(7/3)-(rho/rho_0)**(5/3))*(1-3/4*(4-K_p_0)*((rho/rho_0)**(2/3)-1))
    return pre

#%%
# 超级简化版本状态方程
def eos_density2(pre,mat):
    data = pd.read_csv(mat+'.txt')

    K_0 = data.iloc[0,0]
    K_p_0 = data.iloc[0,1]
    rho_0 = data.iloc[0,2]

    return rho_0*(1+(K_p_0*pre)/K_0)**(1/K_p_0)

#%%
# 计算过程中需要使用的包装函数

def getdensity(pres,mate):
    # if mate == 'curst':
    #     rho_initial = 1000
    # elif mate == 'mantle':
    #     rho_initial = 4260
    # else:
    #     rho_initial = 8340
    # density_i = op.fsolve(eos_density1, rho_initial, args=(pres,mate))
    # return density_i[0]
    density_i = eos_density2(pres,mate)
    return density_i

def getpressure(dens,mate):
    pressure_i = eos_pressure1(dens, mate)
    return pressure_i

def compute_mass(rad_f,rad_b,rho):
    mass_i = inte.quad(m, rad_b, rad_f,args=(rho))[0]
    return mass_i

def compute_gravity(rad_f,rad_b,rho):
    gravity_i = ((4*math.pi*G)/rad_f**2)*inte.quad(g, rad_b, rad_f, args=(rho))[0]
    return gravity_i

def compute_pressure(rad_f,rad_b,rho,gra_f):
    pressure_i = inte.quad(p, rad_b, rad_f, args=(rho))[0]
    return pressure_i

# %%
# 定义求解微分方程的方法

#%%
# 将控制微分方程封装成微分方程组,再每个步长的位置处以上一次的解为条件解一次微分方程组
def ELplanetmodel(var,dr,r,com_m_fe):
    m,p,rho = var
    if m <= 0:
        return 'stop'
    
    else:
        mass_fraction = m/M_p 
        if mass_fraction > 1-M_h2o:
            mate = 'curst'
        elif mass_fraction <= com_m_fe:
            mate = 'core'
        else:
            mate = 'mantle'
        dm = 4*math.pi*(r**2)*rho
        dp = G*(m)*rho/(r**2)
        next_m = m-dm*dr
        next_p = p+dp*dr
        next_rho = getdensity(p,mate)
        return [next_m,next_p,next_rho,mate]

def RKplanetmodel(var,dr,r,com_m_fe):
    m,p,rho = var
    if m <= 0:
        return 'stop'
    
    else:
        mass_fraction = m/M_p 
        if mass_fraction > 1-M_h2o:
            mate = 'curst'
        elif mass_fraction <= com_m_fe:
            mate = 'core'
        else:
            # if p <= 25e9:
            #     mate = 'up_mantle'
            # else:
            #     mate = 'low_mantle'

            mate = 'mantle'

        if r == 0:
            r = 1

        dm1 = 4*math.pi*(r**2)*rho
        dm2 = 4*math.pi*((r+dr/2)**2)*rho
        dm3 = 4*math.pi*((r+dr/2)**2)*rho
        dm4 = 4*math.pi*((r+dr)**2)*rho
        next_m = m-(dm1+2*dm2+2*dm3+dm4)*(dr/6)

        if mate == 'core':
            dp1 = G*(4/3)*math.pi*r*rho**2
            dp2 = G*(4/3)*math.pi*(r+dr/2)*rho**2
            dp3 = G*(4/3)*math.pi*(r+dr/2)*rho**2
            dp4 = G*(4/3)*math.pi*(r+dr)*rho**2
            next_p = p+(dp1+2*dp2+2*dp3+dp4)*(dr/6)
            # dp1 = G*(next_m)*rho/(r**2)
            # dp2 = G*(next_m)*rho/((r+dr/2)**2)
            # dp3 = G*(next_m)*rho/((r+dr/2)**2)
            # dp4 = G*(next_m)*rho/((r+dr)**2)
            # next_p = p+(dp1+2*dp2+2*dp3+dp4)*(dr/6)
        else:
            dp1 = G*(m)*rho/(r**2)
            dp2 = G*(m)*rho/((r+dr/2)**2)
            dp3 = G*(m)*rho/((r+dr/2)**2)
            dp4 = G*(m)*rho/((r+dr)**2)
            next_p = p+(dp1+2*dp2+2*dp3+dp4)*(dr/6)
        
        next_rho = getdensity(p,mate)
        return [next_m,next_p,next_rho,mate]

# %%
