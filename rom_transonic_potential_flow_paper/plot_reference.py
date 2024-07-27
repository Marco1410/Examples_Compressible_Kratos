import numpy as np
import matplotlib.pyplot as plt

#### CP PLOT
######################################################################
fig = plt.figure()
fig.set_figwidth(12.0)
fig.set_figheight(8.0)
data_filename = 'reference_data/flo36/cp_flo36_aoa_1_mach_72.dat'

x_1  = np.loadtxt('reference_data/flo36/cp_flo36_aoa_1_mach_72.dat', usecols=(0,))
cp_1 = np.loadtxt('reference_data/flo36/cp_flo36_aoa_1_mach_72.dat', usecols=(1,))
fig = plt.plot(x_1, cp_1, "o", markersize = 2.0, label = '1')

x_2  = np.loadtxt('reference_data/flo36/cp_flo36_aoa_1_mach_73.dat', usecols=(0,))
cp_2 = np.loadtxt('reference_data/flo36/cp_flo36_aoa_1_mach_73.dat', usecols=(1,))
fig = plt.plot(x_2, cp_2, "o", markersize = 2.0, label = '2')

x_3  = np.loadtxt('reference_data/flo36/cp_flo36_aoa_1_mach_75.dat', usecols=(0,))
cp_3 = np.loadtxt('reference_data/flo36/cp_flo36_aoa_1_mach_75.dat', usecols=(1,))
fig = plt.plot(x_3, cp_3, "o", markersize = 2.0, label = '3')

x_4  = np.loadtxt('reference_data/flo36/cp_flo36_aoa_2_mach_75.dat', usecols=(0,))
cp_4 = np.loadtxt('reference_data/flo36/cp_flo36_aoa_2_mach_75.dat', usecols=(1,))
fig = plt.plot(x_4, cp_4, "o", markersize = 2.0, label = '4')

fig = plt.title('Cp vs x')
fig = plt.ylabel('Cp')
fig = plt.xlabel('x')
fig = plt.grid()
fig = plt.legend()
fig = plt.gca().invert_yaxis()
fig = plt.tight_layout()
# fig = plt.show()
fig = plt.savefig("reference_data_naca0012_flo36.png", dpi=400)
fig = plt.close('all')