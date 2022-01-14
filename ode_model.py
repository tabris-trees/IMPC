# %%
from os import write
import numpy as np
# import pandas as pd
import matplotlib.pyplot as plt
# from scipy import optimize as op
# from scipy import integrate as inte
# import math
import eos_no_tem

# %%
# 星球参数设置区域
# P_s = 0.0 # 行星表面的压力
# P_center = [48.8,51.0] # 行星中心位置处的压强范围 kbar
Rho_s = 1000 # 行星表面的密度 kg/m^3
T_s = 102.0 # 行星表面的温度 K
# d_T = 3.0 # 对流区域的地温梯度 K/km
R_p = 1.565e6# 行星半径 m
M_p = 4.799844e22 # 行星质量，（0.008个地球质量）kg
g_s = 1.314 # 行星表面的重力, m/s^2

# 输入模型参数
M_h2o = 0.05 # 冰壳的质量分数
M_fe = 0.15 # 铁核的质量分数

step = 500 # 计算步长 km

# 参数计算
M_mgfe2sio4 = 1-M_h2o-M_fe # 硅酸盐幔的质量分数


# %%
# 开始计算

# initial

material = ['water']
# mass = [M_p]
# density = [Rho_s]
# # temperature = [T_s]
# radius = [R_p]
# gravity = [g_s]
# pressrue_0 = eos_no_tem.getpressure(density[0], material[0])
# pressrue = [pressrue_0]

# mass = []
# pressure = []
# density = []

# compute

# mate = 'water'
r_range = np.arange(1565000,1475000,-1.0)
# ans = odeint(planetmodel,(4.799844e22,0.0,1000.0),r,args=(mate,))
var_0 = [4.799844e22,0.0,1000.0]
ans = []

var = var_0
for i in r_range:
    vars = eos_no_tem.planetmodel(var,1.0,i)
    var = vars[0:3]
    material.append(vars[-1])
    ans.append(var)

    if i%1000 == 0:
        print(i)

ansarray = np.array(ans)

mass = ansarray[:,0]
pressure = ansarray[:,1]
density = ansarray[:,2]

plt.plot(r_range,mass)
plt.title('mass')
plt.show()
plt.plot(r_range,pressure)
plt.title('pressure')
plt.show()
plt.plot(r_range,density)
plt.title('density')
plt.show()

print('-'*60 )
print(ansarray[:,0][-1]/ansarray[:,0][0])
# %%

def save(file,ls):
    for i in ls:
        file.writelines(str(i)+'\n')

m = open('mass.txt','w')
save(m,mass)
p = open('pressure.txt','w')
save(p,pressure)
rh = open('density.txt','w')
save(rh,density)
# %%
