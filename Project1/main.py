import matplotlib.pyplot as plt
import analytical as an
import numerical as num

temp_an, rad_an = an.ode1(10, 20, 1, 8, 0.1)
temp_eul, rad_eul = num.euler(10, 20, 1, 8, 0.5)
temp_rk4, rad_rk4 = num.rk4(10, 20, 1, 8, 0.5)

fig1, ax1 = plt.subplots()
ax1.plot(rad_an, temp_an, rad_eul, temp_eul, '.-', rad_rk4, temp_rk4, '.-')
ax1.set_xlabel('Radius')
ax1.set_ylabel('Temperature')
ax1.set_title('Cylindrical Temperature Distribution')
ax1.grid()
ax1.legend(['Exact', 'Euler', 'RK4'])
plt.show()
