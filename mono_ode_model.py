# %%
from os import write
import numpy as np
# import pandas as pd
import matplotlib.pyplot as plt
# from scipy import optimize as op
# from scipy import integrate as inte
# import math
import mono_function
import parameter


# %%
# 星球参数设置区域
P_s = parameter.P_s # 行星表面的压力
# P_center = [48.8,51.0] # 行星中心位置处的压强范围 kbar
Rho_s = parameter.Rho_s # 行星表面的密度 kg/m^3
# T_s = 102.0 # 行星表面的温度 K
# d_T = 3.0 # 对流区域的地温梯度 K/km
R_p = parameter.R_p# 行星半径 m
M_p = parameter.M_p # 行星质量，（0.008个地球质量）kg
# g_s = 1.314 # 行星表面的重力, m/s^2

# 输入模型参数
# M_h2o = 0.10 # 壳的质量分数
M_fe = parameter.M_fe # 核的质量分数
# error_M_fe = parameter.error_M_fe # 核的质量分数的误差范围

step = parameter.step # 计算步长 km

# 参数计算
# M_mgfe2sio4 = 1-M_h2o-M_fe # 硅酸盐幔的质量分数



# %%
# 开始计算

# initial
material = ['curst']
r_range = np.arange(R_p,-step,-step)
var_0 = [M_p,P_s,Rho_s]
ans = []
ansls = []

var = var_0

# compute

for i in r_range:
    vars = mono_function.RKplanetmodel(var,step,i)
    if vars == 'stop':
        break
    var = vars[0:3]
    material.append(vars[-1])
    ans.append(var)

    if i%1000 == 0:
        print(i)
# if abs(var[-1]/M_p) <= 0.001:
#     ansls.append(ans)

# %%
# 作图
ansarray = np.array(ans)

mass = ansarray[:,0]
# print(mass)
pressure = ansarray[:,1]
density = ansarray[:,2]

r_stop = R_p-len(mass)*step
print('stop at {}m of radius'.format(r_stop+step))

r_compute = np.arange(R_p,r_stop,-step)
# print(r_compute)


plt.plot(r_compute,mass)
plt.title('mass')
plt.show()
plt.plot(r_compute,pressure)
plt.title('pressure')
plt.show()
plt.plot(r_compute,density)
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
m.close()
p = open('pressure.txt','w')
save(p,pressure)
p.close()
rh = open('density.txt','w')
save(rh,density)
rh.close()
ma = open('material.txt','w')
save(ma,material)
ma.close()

# %%



