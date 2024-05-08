import subprocess
import time as time

start_time = time.time()

salome_cmd = "salome -t python"
salome_script_name = "generate_naca0012_2d.py"

salome_exe = " ".join([salome_cmd, salome_script_name])

sp = subprocess.Popen(["/bin/bash", "-i", "-c", salome_exe])

sp.communicate()

exe_time = time.time() - start_time

print(' Executing SALOME with "' + salome_script_name +
      '" took ' + str(round(exe_time, 2)) + ' sec')
