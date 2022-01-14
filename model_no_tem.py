# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import optimize as op
from scipy import integrate as inte
import math
import eos_no_tem

# %%
# 星球参数设置区域
# P_s = 0.0 # 行星表面的压力
# P_center = [48.8,51.0] # 行星中心位置处的压强范围 kbar
Rho_s = 1000 # 行星表面的密度 kg/m^3
T_s = 102.0 # 行星表面的温度 K
d_T = 3.0 # 对流区域的地温梯度 K/km
R_p = 1565.0# 行星半径 km
M_p = 4.799844e22 # 行星质量，（0.008个地球质量）kg
g_s = 1.314 # 行星表面的重力, m/s^2

# 输入模型参数
M_h2o = 0.10 # 冰壳的质量分数
M_fe = 0.12 # 铁核的质量分数

step = 0.5 # 计算步长 km

# 参数计算
M_mgfe2sio4 = 1-M_h2o-M_fe # 硅酸盐幔的质量分数


# %%
# 开始计算

# initial

material = ['water']
mass = [M_p]
density = [Rho_s]
# temperature = [T_s]
radius = [R_p]
# gravity = [g_s]
pressrue_0 = eos_no_tem.getpressure(density[0], material[0])
pressrue = [pressrue_0]


i = 0 # 计数器

# compute

while radius[-1] != 0:
    i+=1
    radius.append(radius[i-1]-step)
    mass.append(mass[i-1]-eos_no_tem.compute_mass(radius[i-1],radius[i],density[i-1]))
    massfrac = mass[i]/M_p
    if massfrac > M_h2o:
        material.append('water')
    elif M_mgfe2sio4 <= massfrac < M_h2o:
        material.append('silicate')
    else:
        material.append('iron')
    pressrue.append(pressrue[i-1]+eos_no_tem.compute_pressure(radius[i-1],radius[i],density[i-1]))
    density.append(eos_no_tem.getdensity(pressrue[i],material[i-1],density[i-1]))
    

# %%
print(len(radius))

plt.plot(radius,density)
plt.title('density')
plt.show()
plt.plot(radius,mass)
plt.title('mass')
plt.show()
plt.plot(radius,pressrue)
plt.title('pressure')
plt.show()
# plt.plot(radius,temperature)
# plt.title('temperature')
# plt.show()


# %%
print(material)

# %%



