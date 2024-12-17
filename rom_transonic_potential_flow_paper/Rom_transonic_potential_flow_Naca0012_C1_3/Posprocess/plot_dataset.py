import os
import numpy as np
import matplotlib.pyplot as plt


if os.path.exists("mu_train.npy"):
    mu_train = np.load('mu_train.npy')
    mu_train_not_scaled = np.load('mu_train_not_scaled.npy')
    mu_train =  [mu.tolist() for mu in mu_train]
    mu_train_not_scaled =  [mu.tolist() for mu in mu_train_not_scaled]
else:
    mu_train = []
    mu_train_not_scaled = []
if os.path.exists("mu_test.npy"):
    mu_test = np.load('mu_test.npy')
    mu_test_not_scaled = np.load('mu_test_not_scaled.npy')
    mu_test =  [mu.tolist() for mu in mu_test]
    mu_test_not_scaled =  [mu.tolist() for mu in mu_test_not_scaled]
else:
    mu_test = []
    mu_test_not_scaled = []
if os.path.exists("mu_validation.npy"):
    mu_validation = np.load('mu_validation.npy')
    mu_validation_not_scaled = np.load('mu_validation_not_scaled.npy')
    mu_validation =  [mu.tolist() for mu in mu_validation]
    mu_validation_not_scaled =  [mu.tolist() for mu in mu_validation_not_scaled]
else:
    mu_validation = []
    mu_validation_not_scaled = []

# name = "MuValues"
name = "MuValuesNotScaled"

plt.figure(figsize=(10, 6)) 
if len(mu_train) > 0: plt.plot(np.array(mu_train)[:,1], np.array(mu_train)[:,0], 'bs', label="Train Values", markersize=2)
if len(mu_test ) > 0: plt.plot(np.array(mu_test)[:,1], np.array(mu_test)[:,0], 'rx', label="Test Values")
if len(mu_validation ) > 0: plt.plot(np.array(mu_validation)[:,1], np.array(mu_validation)[:,0], 'g*', label="Validation Values")
plt.title('Mu Values')
plt.ylabel('Alpha')
plt.xlabel('Mach')
plt.grid(True)
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
plt.tight_layout()
plt.savefig(f"{name}.pdf", bbox_inches='tight', dpi=400)
plt.close('all')