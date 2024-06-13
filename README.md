# BoostingMOO
The pyomo codes used to produce results for the boosting by intermediate paper are given here. 

# Setup
Install Anaconda in your PC. PYOMO, GLPK and IPOPT 3.11.1 are required to run the optimization codes. These can be installed by running the following commands in the Powershell Prompt of anaconda navigator:

conda install -c conda-forge pyomo

conda install -c conda-forge glpk

conda install -c conda-forge ipopt=3.11.1 

# Run the optimization codes
You can use the the optimization codes (Optimization1_CascadeMOO_BTO_variedIC.py & Optimization2_CascadeMOO_BTO_variedIC.py) to produce all optimization results in our paper and more. Save both files in the same directory. Open both files in Spyder and run the Optimization1_CascadeMOO_BTO_variedIC.py file. You can vary the values of the parameters to produce different pareto frontiers. 

# Publications
When using this work, please cite our paper:

Leandros Paschalidis, Daniela Fröschl, Manuel Ibañez, Samuel Sutiono, Volker Sieber, Jakob Burger, Boosting of enyzmatic cascades by intermediates: theoretical analysis and model-based optimization, Submited to Biochemical Engineering Journal
