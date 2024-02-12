# 1. About

This repository contains a tool to automatically i) generate simple meshes ii) prepare numerical simulations input files through a parametrization iii) run the simulation locally. 

The current version supports 2D external aerodynamics cases, but it could be extended to more complex configurations. 

While it can be configured for a wide range of solvers, examples of simulations (to be found in `templates`) are provided for the particular cases:

- The popular panel method included in Xfoil software.
- The open-source code Kratos Multiphysics, with full potential formulation.

# 2. Requirements and setup

It is assumed you are using a Linux OS distribution. You will need a working Python 3 environment in your local machine, with `matplotlib`, `numpy` and `scipy` installed. 

You will need a tailored version of meshio, as it is used in connection to Kratos meshes manipulation:
 - Clone the repository somewhere in your file system (e.g. in $HOME/git) with the command `git clone https://github.com/marandra/meshio.git`. Note that these files are to be kept in your syste..
 - Everytime these scripts will be used, you will need to type `export PYTHONPATH=$HOME/git/meshio/src:$PYTHONPATH` so that the right version of meshio is loaded. Therefore, we recommend that you copy this line in your `.bashrc` file

>>>
**NOTE**: This solution is temporary, as the developments will be merged at some point with the meshio version that can be directly pip installed.
>>>

You will also need to install the mesh generator gmsh [https://gmsh.info/]. The syntax here depends on your Linux distribution. The python interface should be accessible from your Python distribution (check by trying `import gmsh`).

You will need to install all the solvers that will be run locally. 
- For Kratos, please follow the instructions of https://github.com/KratosMultiphysics/Kratos, after cloning it, to install it in your Python distribution. Choose the latest master for that. Remember to include the potential flow module as an app, as it is needed (in your configuration script: `add_app ${KRATOS_APP_DIR}/CompressiblePotentialFlowApplication`). 
- For XFOIL, you will need to install the `xfoil` package of your Linux distro.