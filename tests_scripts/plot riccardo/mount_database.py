import numpy as np
import json
import os, glob

def MountDatabase(path):
    all_params = []
    all_velocities = []
    all_temperatures = []

    for filename in glob.glob(os.path.join(path, '*.json')):
        print("filename = ",filename)
        with open(os.path.join(os.getcwd(), filename), 'r') as f: # open in readonly mode
            d = json.load(f)
            all_params.append( np.array([d["inlet_mass_flow"],d["inlet_temperature"]],dtype=float) )
            all_velocities.append( np.array(d["velocity"],dtype=float) )
            all_temperatures.append( np.array(d["temperature"],dtype=float) )

    all_params = np.column_stack(all_params)
    all_velocities = np.column_stack(all_velocities)
    all_temperatures = np.column_stack(all_temperatures)
    return all_params,all_velocities,all_temperatures


