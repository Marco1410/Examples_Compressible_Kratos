This repo contains two mini-apps concerning the numerical simulation of NACA airfoils, and developed in the framework of NextSim. The project was funded by the European Union’s Horizon 2020 research and innovation program, under the grant agreement No. 956104, and the Spanish Ministry of science and Innovation (MICINN), under the programme "Acciones Programación Conjunta Internacional (APCIN)" grant number PCI2021-121945.

# 1. Get ready

Go to cfddispatcher, a common set of utilities shared by both mini-apps, and follow the instructions provided by the Readme.

# 2. Neural Networks mini-app

This mini-app allows the user to run a full potential transonic simulation using the code Kratos Multiphysics, and to subsequently apply an Euler-based neural network correction generated with CODA. More details are to be found in [TODO CITE]. To run this mini-app, follow the steps below:

- Go to neural_networks_miniapp subfolder 
- Edit run_cases.py
 - Modify "INPUT DATA" with the airfoil you want to simulate
 - Modify "USER-RELATED INPUT DATA" to adapt to your setup
- Execute run_cases.py

The mini-app will generate a figure (pressure_coeffs.png), comparing the results of the original Kratos Multiphysics simulation with their corrected counterpart.

# 3. Optmizitation mini-app

This mini-app illustrates the development of a framework for multi-fidelity optimization, and represents the basis of the strategy described in [TODO CITE]. This is done by performing two simple numerical optimizations, both aiming at maximing the lift coefficient (objective function) by varying the airfoil thickness (design variable). In particular, the panel method implemented in XFOIL is used in what he have called a lower fidelity (LF) optimization. For the higher fidelity (HF), the full potential model of Kratos Multiphysics is employed (assuming a compressible element formulation). In order to use this mini-app, you should first install the genetic algorithm it is based on, by e.g. typing 'pip install mini-app geneticalgorithm'. Then follow the steps below:

- Go to optimizitation_miniapp subfolder 
- Modify "USER-RELATED INPUT DATA" in full_potential.py to adapt Kratos execution to your setup
- Modify "USER-RELATED INPUT DATA" in run_optimization.py to point to your Python interpreter
- Modify "INPUT DATA" in run_optimization.py with the optimization you would like to perform
- Execute run_optimization.py

The mini-app will generate a figure (optima_found.png) with all the simulations involved for both LF and HF optimizations, highlighting the optima which were found. As the trends are somehow trivial (the higher the thickness, the higher the lift), the optima should be found at one of the edges of the design space.