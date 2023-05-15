import numpy as np
import matplotlib.pyplot as plt

# astronomical unit  [km]
au = 1.496e8
# solar wind speed [km/s]
v = 500/3600
# angular velocity sun [1/s]
omega = 2.9e-4 #2*np.pi/27.28/86400 # 2.9e-4
print(omega)

# step size
stepsize = 8e6

'''Irgendwas funktioniert mit dem Omega nicht!!!'''


r = np.arange(0, 2*au, step=stepsize)


# Sonnenwind bewegt sich mit v radial nach außen. Die Zeit berechnet sich als r/v
# Ankunft an 1 au:
t_arrival = au/v

# theta = 2*np.pi*omega*np.linspace(0,t_arrival,steps)
theta = omega*r/v
# theta = np.degrees(theta)
print(r)
print(theta)


fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
ax.plot(theta, r)
ax.set_rmax(2*au)
ax.set_rticks(np.array([0.5, 1, 1.5, 2])*au)  # Less radial ticks
ax.set_rlabel_position(-22.5)  # Move radial labels away from plotted line
ax.grid(True)

ax.set_title("A line plot on a polar axis", va='bottom')
plt.show()