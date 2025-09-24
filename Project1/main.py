import matplotlib.pyplot as plt
import numpy as np
import analytical as an
import numerical as num

temp_an, rad_an = an.ode1(10, 20, 1, 8, 0.1)
temp_eul, rad_eul = num.euler(10, 20, 1, 8, 0.5)
temp_rk4, rad_rk4 = num.rk4(10, 20, 1, 8, 0.5)
temp_scipy, rad_scipy = num.scipy_ode_solve(10, 20, 1, 8, 0.5)

lify_intervals = 4
lift_an = an.int1(1.225, 10, 1, 5)
unit_lift_an, y_an = an.int1_curve(1.225, 10, 1, 5, 50)
lift_rie, unit_lift_rie, y_rie = num.riemann(1.225, 10, 1, 5, lify_intervals)
unit_lift_plot_rie, y_plot_rie = num.riemann_plot(unit_lift_rie, y_rie)
lift_trap, unit_lift_trap, y_trap = num.trapezoidal(1.225, 10, 1, 5, lify_intervals)
unit_lift_plot_trap, y_plot_trap = num.trapezoidal_plot(unit_lift_trap, y_trap)
lift_simp = num.simpson(1.225, 10, 1, 5, lify_intervals)
lift_scipy = num.scipy_int_solve(1.225, 10, 1, 5)
lift_trap_scipy = num.scipy_trapezoidal(1.225, 10, 1, 5, lify_intervals)
lift_simp_scipy = num.scipy_simpson(1.225, 10, 1, 5, lify_intervals)

lift_labels = np.array(['Exact', 'Riemann Sum', 'Trapezoidal Rule', 'Simpsons Rule', 'SciPy Generic', 'SciPy Trapezoidal', 'SciPy Simpsons'])
lift_values = np.array([lift_an, lift_rie, lift_trap, lift_simp, lift_scipy, lift_trap_scipy, lift_simp_scipy])
lifts = np.c_[lift_values, lift_labels]
print('Lift results, number of intervals =', lify_intervals, ':')
print(lifts)

fig1, ax1 = plt.subplots()
ax1.plot(rad_an, temp_an, rad_eul, temp_eul, '.:', rad_rk4, temp_rk4, '.:',
         rad_scipy, temp_scipy, '.:')
ax1.set_xlabel('Radius')
ax1.set_ylabel('Temperature')
ax1.set_title('Temperature Distribution in a Hollow Cylinder')
ax1.grid()
ax1.legend(['Exact', 'Euler', 'RK4', 'SciPy'])
plt.show()

fig2, ax2 = plt.subplots()
ax2.plot(y_an, unit_lift_an, y_plot_rie, unit_lift_plot_rie)
ax2.set_xlabel('Spanwise Position')
ax2.set_ylabel('Lift per Unit Span')
ax2.set_title('Spanwise Elliptical Lift Distribution')
ax2.grid()
ax2.legend(['Exact', 'Riemann Sum'], loc='upper left')
plt.show()

fig3, ax3 = plt.subplots()
ax3.plot(y_an, unit_lift_an, y_plot_trap, unit_lift_plot_trap)
ax3.set_xlabel('Spanwise Position')
ax3.set_ylabel('Lift per Unit Span')
ax3.set_title('Spanwise Elliptical Lift Distribution')
ax3.grid()
ax3.legend(['Exact', 'Trapezoidal Rule'], loc='upper left')
plt.show()
