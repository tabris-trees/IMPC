#%%
from scipy.integrate import odeint 
import numpy as np 

def lorenz(w, t, p, r, b): 
    # 给出位置矢量w，和三个参数p, r, b计算出
    # dx/dt, dy/dt, dz/dt的值
    x, y, z = w
    # 直接与lorenz的计算公式对应 
    return np.array([p*(y-x), x*(r-z)-y, x*y-b*z]) 

t = np.arange(0, 30, 0.01) # 创建时间点 
# 调用ode对lorenz进行求解, 用两个不同的初始值 
track1 = odeint(lorenz, (0.0, 1.00, 0.0), t, args=(10.0, 28.0, 3.0)) 
track2 = odeint(lorenz, (0.0, 30, 0.0), t, args=(10.0, 28.0, 3.0)) 

# 绘图
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt 

fig = plt.figure()
ax = Axes3D(fig)
ax.plot(track1[:,0], track1[:,1], track1[:,2])
ax.plot(track2[:,0], track2[:,1], track2[:,2])
plt.show()
# %%
# -*- coding: utf8 -*-
import numpy as np

"""
移动方程：
t时刻的位置P(x,y,z)
steps：dt的大小
sets：相关参数
"""


def move(P, steps, sets):
    x, y, z = P
    sgima, rho, beta = sets
    # 各方向的速度近似
    dx = sgima * (y - x)
    dy = x * (rho - z) - y
    dz = x * y - beta * z
    return [x + dx * steps, y + dy * steps, z + dz * steps]


# 设置sets参数
sets = [10., 28., 3.]
t = np.arange(0, 30, 0.01)

# 位置1：
P0 = [0., 1., 0.]

P = P0
d = []
for v in t:
    P = move(P, 0.01, sets)
    d.append(P)
dnp = np.array(d)

# 位置2：
P02 = [0., 1.01, 0.]

P = P02
d = []
for v in t:
    P = move(P, 0.01, sets)
    d.append(P)
dnp2 = np.array(d)
"""
画图
"""
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

fig = plt.figure()
ax = Axes3D(fig)
ax.plot(dnp[:, 0], dnp[:, 1], dnp[:, 2])
ax.plot(dnp2[:, 0], dnp2[:, 1], dnp2[:, 2])
plt.show()
# %%
# -*- coding: utf-8 -*-
import numpy as np
from scipy.integrate import odeint
"""
定义常微分方程，给出各方向导数,即速度
"""
def dmove(Point,t,sets):
    """
    p：位置矢量
    sets：其他参数
    """
    p,r,b = sets
    x,y,z = Point
    return np.array([p*(y-x),x*(r-z),x*y-b*z])
 
t = np.arange(30,60,0.01)
#调用odeint对dmove进行求解，用两个不同的初始值
P1 = odeint(dmove,(0.,1.,0.),t,args = ([10.,28.,3.],))  #(0.,1.,0.)是point的初值
#([10.,28.,3.],)以元祖的形式给出 point,t之后的参数
P2 = odeint(dmove,(0.,1.01,0.),t,args = ([10.,28.,3.],))
 
"""
画3维空间的曲线
"""
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
 
fig = plt.figure()
ax = Axes3D(fig)
ax.plot(P1[:,0],P1[:,1],P1[:,2])
ax.plot(P2[:,0],P2[:,1],P2[:,2])
plt.show()
# %%
dp1 = G*(4/3)*math.pi*r*rho**2
        dp2 = G*(4/3)*math.pi*(r+dr/2)*rho**2
        dp3 = G*(4/3)*math.pi*(r+dr/2)*rho**2
        dp4 = G*(4/3)*math.pi*(r+dr)*rho**2
        next_p = p+(dp1+2*dp2+2*dp3+dp4)*(dr/6)