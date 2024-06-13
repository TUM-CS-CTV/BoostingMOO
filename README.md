# BoostingMOO
The pyomo codes used to produce results for the boosting by intermediate paper are given here. 

# Setup
Install Anaconda in your PC. PYOMO, GLPK and IPOPT 3.11.1 are required to run the optimization codes. These can be installed by running the following commands in the Powershell Prompt of anaconda navigator:

conda install -c conda-forge pyomo

conda install -c conda-forge glpk

conda install -c conda-forge ipopt=3.11.1 

# Run the simulation code
The simulation code used to produce Figure 2 is available to download and run (Simulation for CRN a.py). 

# Run the optimization codes
The optimization codes used to produce Figure 3 are available to download and run (Optimization1_CascadeMOO_BTO_variedIC.py & Optimization2_CascadeMOO_BTO_variedIC.py) Both files must be saved in the same directory. To produce the results one has to run the Optimization1_CascadeMOO_BTO_variedIC.py file. You can vary the values of the parameters to produce different pareto frontiers. 

# Publications
When using this work, please cite our paper:

Leandros Paschalidis, Daniela Fröschl, Manuel Ibañez, Samuel Sutiono, Volker Sieber, Jakob Burger, Boosting of enyzmatic cascades by intermediates: theoretical analysis and model-based optimization, Submited to Biochemical Engineering Journal


# Further information
For more information about how multi-objective optimization is applied in the design of enzymatic cascade reaction processes refer to the CascadeMOO folder. 
