import matplotlib.pyplot as plt
import analytical as an
import numerical as nu

temp_an, rad_an = an.ode1(10, 20, 1, 8, 0.1)
temp_eul, rad_eul = nu.euler(10, 20, 1, 8, 0.1)

fig1, ax1 = plt.subplots()
ax1.plot(rad_an, temp_an, rad_eul, temp_eul, '.-')
ax1.set_xlabel('Radius')
ax1.set_ylabel('Temperature')
ax1.set_title('Cylindrical Temperature Distribution')
ax1.grid()
plt.show()
