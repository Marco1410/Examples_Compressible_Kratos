import os
import numpy as np
import matplotlib.pyplot as plt

directory        = '150/'
case_name        = '1.99780167280766, 0.7408431899218146'
capture_filename = 'test.png'


#### CP PLOT
######################################################################
cp_min = cp_max = cp_hrom = cp_fom = cp_rom = cp_rbf = 0
fig, (axsrc, axzoom) = plt.subplots(1, 2, figsize=(12, 8))
axsrc.set(xlim=(0, 1), ylim=(0, 1), autoscale_on=False,
          title='Click to zoom')
axzoom.set(xlim=(0.45, 0.55), ylim=(0.4, 0.6), autoscale_on=False,
           title='Zoom window')

#### FOM ######
fom_name = f'{directory}FOM_Snapshots/{case_name}.npy'
if os.path.exists(fom_name):
    fom = np.load(fom_name)
    fom_skin_data_filename = f"{directory}FOM_Skin_Data/{case_name}.dat"
    x_fom  = np.loadtxt(fom_skin_data_filename, usecols=(0,))
    cp_fom = np.loadtxt(fom_skin_data_filename, usecols=(3,))
    axsrc.plot(x_fom, cp_fom, 'ob', markersize = 1.0, label = 'FOM')
#### ROM ######
rom_name = f'{directory}ROM_Snapshots/{case_name}.npy'
if os.path.exists(rom_name):
    rom = np.load(rom_name)
    rom_skin_data_filename = f"{directory}ROM_Skin_Data/{case_name}.dat"
    x_rom  = np.loadtxt(rom_skin_data_filename, usecols=(0,))
    cp_rom = np.loadtxt(rom_skin_data_filename, usecols=(3,))
    axsrc.plot(x_rom, cp_rom, 'xr', markersize = 2.0, label = f'ROM-FOM e: {(np.linalg.norm(fom-rom)/np.linalg.norm(fom)):.2E}')
#### RBF ####
rbf_name = f"{directory}RBF_Snapshots/{case_name}.npy"
if os.path.exists(rbf_name):
    rbf = np.load(rbf_name).T
    rbf_skin_data_filename = f"{directory}RBF_Skin_Data/{case_name}.dat"
    x_rbf  = np.loadtxt(rbf_skin_data_filename, usecols=(0,))
    cp_rbf = np.loadtxt(rbf_skin_data_filename, usecols=(3,))
    axsrc.plot(x_rbf, cp_rbf, '+m', markersize = 2.0, label = f'RBF-FOM e: {(np.linalg.norm(fom-rbf)/np.linalg.norm(fom)):.2E}')
cp_min = np.min([np.min(cp_hrom), np.min(cp_fom), np.min(cp_rom), np.min(cp_rbf)])
cp_max = np.max([np.max(cp_hrom), np.max(cp_fom), np.max(cp_rom), np.max(cp_rbf)])

def on_press(event):
    if event.button != 1:
        return
    x, y = event.xdata, event.ydata
    axzoom.set_xlim(x - 0.1, x + 0.1)
    axzoom.set_ylim(y - 0.1, y + 0.1)
    fig.canvas.draw()

fig.canvas.mpl_connect('button_press_event', on_press)
fig = plt.title('Cp vs x')
fig = plt.axis([-0.05,1.35,cp_max+0.1,cp_min-0.1])
fig = plt.ylabel('Cp')
fig = plt.xlabel('x')
fig = plt.grid()
fig = plt.legend()
fig = plt.tight_layout()
fig = plt.show()

# fig = plt.savefig(capture_filename)
# fig = plt.close('all')