import numpy as np
import torch
from scipy import spatial
from matplotlib import pyplot as plt

from ezyrb import POD, RBF, Database, Snapshot, Parameter, Linear, ANN
from ezyrb import ReducedOrderModel as ROM
from ezyrb.plugin import AutomaticShiftSnapshots

def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

def wave(t, res=256):   
    x = np.linspace(0, 10, res)
    return x, gaussian(x, t, 0.1)   # parameterizing mean value

n_params = 15
params = np.linspace(0.5, 4.5, n_params).reshape(-1, 1)
print(params.shape)


pod = POD() #rank=1 
rbf = RBF()
db = Database()
for param in params:
    space, values = wave(param)
    snap = Snapshot(values=values, space=space)
    db.add(Parameter(param), snap)
    plt.rcParams["figure.figsize"] = (10,3.5)
    plt.plot(space, values, label = param)    #################### mod
    plt.ylabel('$f_{g}$(X)') 
    plt.xlabel('X')
    plt.legend(prop={'size': 6})

snap = db.snapshots_matrix
pa = db.parameters_matrix
print(snap.shape, pa.shape)

db_train, db_test = db.split([len(db)-3, 3])
print("lenght of training data set :", len(db_train))
print("lenght of test data set :", len(db_test))

interp = ANN([10, 10], torch.nn.Softplus(), 1000, frequency_print=200, lr=0.03)
shift  = ANN([], torch.nn.LeakyReLU(), [1e-3, 4000], optimizer=torch.optim.Adam, frequency_print=200, l2_regularization=0,  lr=0.0023)
nnspod = AutomaticShiftSnapshots(shift, interp, Linear(fill_value=0.0), parameter_index=2, reference_index=2, barycenter_loss=20.)
rom = ROM(db_train, pod, rbf, plugins=[nnspod])

for _ in range(10):
    rom.fit()    # Calculate reduced space
    if rom.plugins[0].shift_network.loss_trend[-1] < 1e-3:
        break


modes = pod.modes
m = modes.transpose()
print(m.shape)
plt.plot(space, pod.modes)
plt.title('All Modes')
plt.ylabel('$f_{g}$(X)') 
plt.xlabel('X')
plt.show()

"""
# To plot each mode
for i in range(0, len(db_train)):
    plt.plot(space, m[i, :]*-1, label = i+1 )
    plt.legend()
    plt.title('Mode = {0}'.format(i+1))
    plt.show()
"""

pred = rom.predict(db_test.parameters_matrix) # Calculate predicted solution for given mu

for i in range(len(pred)):
    plt.plot(space, pred.snapshots_matrix[i], label = pred.parameters_matrix[i])
    plt.legend()
    plt.ylabel('$f_{g}$(X)') 
    plt.xlabel('X')
    plt.title('Snapshots corresponding to db_test set shifted to reference position')

error = 0.0

for (_, snap), (_, truth_snap) in zip(pred._pairs, db_test._pairs):
    tree = spatial.KDTree(truth_snap.space.reshape(-1, 1))
    for coord, value in zip(snap.space, snap.values):
        a = tree.query(coord)
        error += np.abs(value - truth_snap.values[a[1]])

assert error < 25

print(error)

new_params = np.linspace(5, 9.5, n_params).reshape(-1, 1)
pred_new = rom.predict(new_params)

for param in new_params:
    x, y = wave(param)
    plt.rcParams["figure.figsize"] = (10,3.5)
    plt.plot(x,y, label = param)
    plt.legend(prop={'size': 7})
    plt.ylabel('$f_{g}$(X)') 
    plt.xlabel('X')
    plt.title('Snapshots corresponding to new parameters in original position')

# prior to shift 
x, y = wave(new_params)
U, s = np.linalg.svd(y, full_matrices=False)[:2]
N_modes = np.linspace(1, len(new_params),len(new_params))
plt.plot(N_modes, s, ".-")
plt.ylabel('Singular values')
plt.xlabel('Modes')
plt.title('Singular Values obtained for snapshots in original position')
plt.show()

pred_new = rom.predict(new_params)
p = pred_new.parameters_matrix
l = pred_new.snapshots_matrix
#print(p.shape, l.shape)
for i in range(len(pred_new)):
    plt.rcParams["figure.figsize"] = (10,3.5)
    plt.plot(space, pred_new.snapshots_matrix[i], label = pred_new.parameters_matrix[i])
    plt.legend(prop={'size': 7})
    plt.ylabel('$f_{g}$(X)') 
    plt.xlabel('X')
    plt.title('Snapshots corresponding to new parameters shifted to reference position')

# After the shift-based pre-processing
U_new, s_new = np.linalg.svd(pred_new.snapshots_matrix, full_matrices=False)[:2]

N_modes = np.linspace(1, len(new_params),len(new_params))

plt.plot(N_modes, s_new, ".-")
plt.ylabel('Singular values')
plt.xlabel('Modes')
plt.title('Singular Values obtained after the shift-based pre-processing')
plt.show()