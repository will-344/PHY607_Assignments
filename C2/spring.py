import numpy as np
import matplotlib.pyplot as plt

T = 50
dt = 0.1
x0 = 1
v0 = 0
k = 1
m = 1

x_EE = np.array([x0])
v_EE = np.array([v0])
x_ES = np.array([x0])
v_ES = np.array([v0])
x_RK2 = np.array([x0])
v_RK2 = np.array([v0])
nIter = int(T/dt)
print(nIter)

def position_velocity_update_EE(x, v, F, dt):
    a = F/m
    v_new = v + a*dt
    x_new = x + v*dt
    return x_new, v_new

def position_velocity_update_ES(x, v, F, dt):
    a = F/m
    v_new = v + a*dt
    x_new = x + v_new*dt
    return x_new, v_new

def position_velocity_update_RK2(x, v, F, dt, k):
    a = F/m
    k1x = v
    k1v = a
    k2x = v + k1v*dt
    k2v = force_update(k, x + k1x*dt)/m
    x_new = x + 1/2*(k1x + k2x)*dt
    v_new = v + 1/2*(k1v + k2v)*dt
    return x_new, v_new

def force_update(k, x):
    F_new = -k*x
    return F_new

for i in range(nIter):
    F_EE = force_update(k, x_EE[i])
    x_new_EE, v_new_EE = position_velocity_update_EE(x_EE[i], v_EE[i], F_EE, dt)
    x_EE = np.append(x_EE, x_new_EE)
    v_EE = np.append(v_EE, v_new_EE)
    
    F_ES = force_update(k, x_ES[i])
    x_new_ES, v_new_ES = position_velocity_update_ES(x_ES[i], v_ES[i], F_ES, dt)
    x_ES = np.append(x_ES, x_new_ES)
    v_ES = np.append(v_ES, v_new_ES)
    
    F_RK2 = force_update(k, x_RK2[i])
    x_new_RK2, v_new_RK2 = position_velocity_update_RK2(x_RK2[i], v_RK2[i], F_RK2, dt, k)
    x_RK2 = np.append(x_RK2, x_new_RK2)
    v_RK2 = np.append(v_RK2, v_new_RK2)

fig, axs = plt.subplots(nrows=3, ncols=1)

axs[0].plot(np.linspace(0,T,nIter+1),x_EE)
axs[0].set_xlabel('Time (s)')
axs[0].set_ylabel('Position (m)')
axs[0].set_title('Spring Oscillator w/ Euler Explicit Method, dt=0.1')

axs[1].plot(np.linspace(0,T,nIter+1),x_ES)
axs[1].set_xlabel('Time (s)')
axs[1].set_ylabel('Position (m)')
axs[1].set_title('Spring Oscillator w/ Euler Symplectic Method, dt=0.1')

axs[2].plot(np.linspace(0,T,nIter+1),x_RK2)
axs[2].set_xlabel('Time (s)')
axs[2].set_ylabel('Position (m)')
axs[2].set_title("Spring Oscillator w/ RK2 (Heun's) Method, dt=0.1")

fig.tight_layout()
plt.show()
