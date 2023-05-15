import numpy as np
import matplotlib.pyplot as plt

# astronomical unit  [m]
au = 1.5e11
# solar wind speed [m/s]
v = 500*3.6
# angular velocity earth [1/s]
omega = 2*np.pi/86400

# number of steps
steps = 1000



r = np.linspace(0.0, au, steps)


# Sonnenwind bewegt sich mit v radial nach außen. Die Zeit berechnet sich als r/v
# Ankunft an 1 au:
t_arrival = au/v

# theta = 2*np.pi*omega*np.linspace(0,t_arrival,steps)
theta = omega*r/v

print(theta)



fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
ax.plot(theta, r)
ax.set_rmax(2)
ax.set_rticks([0.5, 1, 1.5, 2])  # Less radial ticks
ax.set_rlabel_position(-22.5)  # Move radial labels away from plotted line
ax.grid(True)

ax.set_title("A line plot on a polar axis", va='bottom')
plt.show()