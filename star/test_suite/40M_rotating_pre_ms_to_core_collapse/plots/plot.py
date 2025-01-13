import numpy as np
import matplotlib.pyplot as plt
import mesa_reader
from matplotlib import rc


model = mesa_reader.MesaData('profile46.data',file_type = None)

radius = model.data('radius_cm')
velocity = model.data('vel_km_per_s')

logr = np.log10(radius)

plt.plot(logr,velocity)

plt.xlim([6,11])
plt.xlabel(r'log10(Radius) (cm)')
plt.ylabel(r'Velocity (km/s)')
plt.savefig('velocity_at_cc2.pdf',dpi =300)
