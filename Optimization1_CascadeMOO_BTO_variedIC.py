import pyomo.environ as pe
from pyomo.dae import *
import Optimization2_CascadeMOO_BTO_variedIC as FfPF
import pylab
import matplotlib.pyplot as plt

SpaceTimeYield=[]
CofactorConsumption=[]
EnzymeConsumption=[]
IntermediateConsumption=[]
FinalTime=[]
InitialS1Concentration=[]
InitialIntermediateConcentration=[]
InitialS7Concentration=[]
EHsXylDH1=[]
EPuDHT=[]
ELlKdcA=[]
EEcAdhZ3=[]
TotalEnzymeConcentration=[]
info1=[]
info2=[]



i = 1
w1 = 1

# a for loop for the epsilon-IC constraint
for ICi in [0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06]:
    # initial values for the control variables selected manually for the first two runs
    if ICi == 0: 
        tfi = 1000
        S1i = 10
        S2i = 10
        S3i = 10
        S4i = 10
        S5i = 10
        S7i = 10
        EHsXylDH1i = 10
        EPuDHTi = 10
        ELlKdcAi = 10
        EEcAdhZ3i = 10
    elif ICi == 0.01: 
        tfi = 1000
        S1i = 50
        S2i = 10
        S3i = 10
        S4i = 10
        S5i = 10
        S7i = 10
        EHsXylDH1i = 10
        EPuDHTi = 10
        ELlKdcAi = 10
        EEcAdhZ3i = 10
    else: 
        # the control variables are initialized from the solutions of the previous run
        # in order to help with faster conversion
        tfi = model.tf
        S1i = model.S1[model.tau.first()] 
        S2i = model.S2[model.tau.first()]
        S3i = model.S3[model.tau.first()]
        S4i = model.S4[model.tau.first()]
        S5i = model.S5[model.tau.first()]
        S7i = model.S7[model.tau.first()]
        EHsXylDH1i = model.EHsXylDH1
        EPuDHTi = model.EPuDHT
        ELlKdcAi = model.ELlKdcA
        EEcAdhZ3i = model.EEcAdhZ3
        
    model = FfPF.FunctionforPF(ICi,tfi,S1i,S2i,S3i,S4i,S5i,S7i,EHsXylDH1i,
                               EPuDHTi,ELlKdcAi,EEcAdhZ3i)
    # solver selection
    solver=pe.SolverFactory('ipopt')
    # maximum iteration count limit 
    solver.options['max_iter'] = 100000
    # selection of the acceptible tolerance to be stricter than default 
    solver.options['acceptable_tol'] = 10**(-10)
    results = solver.solve(model, tee=True)
    
    # saving the results of each optimization run inside the for loop 
    SpaceTimeYield.append(pe.value(model.OBJ[1]))  
    EnzymeConsumption.append(pe.value(model.EC))
    CofactorConsumption.append(pe.value(model.CC))
    IntermediateConsumption.append(pe.value(model.IC))
    FinalTime.append(pe.value(model.tf))
    InitialS1Concentration.append(pe.value(model.S1[model.tau.first()]))
    InitialIntermediateConcentration.append(pe.value(model.SumIntermediates))
    InitialS7Concentration.append(pe.value(model.S7[model.tau.first()]))
    EHsXylDH1.append(pe.value(model.EHsXylDH1))
    EPuDHT.append(pe.value(model.EPuDHT))
    ELlKdcA.append(pe.value(model.ELlKdcA))
    EEcAdhZ3.append(pe.value(model.EEcAdhZ3))
    TotalEnzymeConcentration.append(pe.value(model.SumEnzymes))
    info1.append(results.solver.termination_condition)
    info2.append(results.solver.status)


    # printing the results of each optimization on an individual text file 
    print('BTO cascade', file = open("ParetoOptimalPoint{}.txt".format(w1), "w+"))
    print('tau', file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print(list(model.tau), file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print('S1', file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print([pe.value(model.S1[jo]) for jo in model.tau], file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print('S2', file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print([pe.value(model.S2[jo]) for jo in model.tau], file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print('S3', file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print([pe.value(model.S3[jo]) for jo in model.tau], file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print('S4', file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print([pe.value(model.S4[jo]) for jo in model.tau], file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print('S5', file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print([pe.value(model.S5[jo]) for jo in model.tau], file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print('S6', file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print([pe.value(model.S6[jo]) for jo in model.tau], file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print('S7', file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print([pe.value(model.S7[jo]) for jo in model.tau], file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print('S8', file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print([pe.value(model.S8[jo]) for jo in model.tau], file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print('vI', file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print([pe.value(model.vI[jo]) for jo in model.tau], file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print('vII', file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print([pe.value(model.vII[jo]) for jo in model.tau], file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print('vIII', file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print([pe.value(model.vIII[jo]) for jo in model.tau], file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print('vIV', file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print([pe.value(model.vIV[jo]) for jo in model.tau], file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print('vV', file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print([pe.value(model.vV[jo]) for jo in model.tau], file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print('EpHsXylDH1', file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print(pe.value(model.EHsXylDH1), file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print('EpPuDHT', file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print(pe.value(model.EPuDHT), file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print('EpLlKdcA', file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print(pe.value(model.ELlKdcA), file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print('EpEcAdhZ3', file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print(pe.value(model.EEcAdhZ3), file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print('tf', file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print(pe.value(model.tf), file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print('SumEnzymes', file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print(pe.value(model.SumEnzymes), file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print('OBJF', file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print(pe.value(model.OBJ[1]), file = open("ParetoOptimalPoint{}.txt".format(w1), "a")) 
    print('Yield', file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print(pe.value(model.Yield), file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print('S1[0]', file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print(pe.value(model.S1[model.tau.first()]), file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print('S7[0]', file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print(pe.value(model.S7[model.tau.first()]), file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    
    # definition of the realtime variable for plotting and printing
    tfp = pe.value(model.tf)
    realtime = [element * tfp for element in list(model.tau)]
        
    EHsXylDH1PLOT= {} 
    for t in list(realtime):
        EHsXylDH1PLOT[t]= pe.value(model.EHsXylDH1)
        
    EPuDHTPLOT= {} 
    for t in list(realtime):
        EPuDHTPLOT[t]= pe.value(model.EPuDHT)
        
    ELlKdcAPLOT= {} 
    for t in list(realtime):
        ELlKdcAPLOT[t]= pe.value(model.ELlKdcA)
        
    EEcAdhZ3PLOT= {} 
    for t in list(realtime):
        EEcAdhZ3PLOT[t]= pe.value(model.EEcAdhZ3)

    # plotting the process schedules     
    plt.figure(figsize=(9,8))
    plt.subplot(211)
    plt.plot(realtime, [pe.value(model.S1[x]) for x in model.tau], label = 'S1', c = 'b', ls = '-', lw = '2')
    plt.plot(realtime, [pe.value(model.S2[x]) for x in model.tau], label = 'S2', c = 'c', ls = '-', lw = '2')
    plt.plot(realtime, [pe.value(model.S3[x]) for x in model.tau], label = 'S3', c = 'm', ls = '-', lw = '2')
    plt.plot(realtime, [pe.value(model.S4[x]) for x in model.tau], label = 'S4', c = 'y', ls = '-', lw = '2')
    plt.plot(realtime, [pe.value(model.S5[x]) for x in model.tau], label = 'S5', c = 'r', ls = '-', lw = '2')
    plt.plot(realtime, [pe.value(model.S6[x]) for x in model.tau], label = 'S6', c = 'k', ls = '-', lw = '2')
    plt.plot(realtime, [pe.value(model.S7[x]) for x in model.tau], label = 'S7', c = 'g', ls = '-', lw = '2')
    plt.plot(realtime, [pe.value(model.S8[x]) for x in model.tau], label = 'S8', c = '0.75',  ls = '-', lw = '2')
    plt.xlabel('$\it{t}$ / (min)', fontsize=17)
    plt.ylabel('$\it{S}$$_i$ / (mM)', fontsize=17)
    plt.grid(False)
    plt.xticks(fontsize=17)
    plt.yticks(fontsize=17)
    plt.legend()

    plt.subplot(212)
    plt.plot(realtime, [pe.value(EHsXylDH1PLOT[x]) for x in realtime], label = 'GlucD', c = 'r', ls = '-', lw = '2')
    plt.plot(realtime, [pe.value(EPuDHTPLOT[x]) for x in realtime], label = 'KdgD', c = 'c', ls = '-', lw = '2')
    plt.plot(realtime, [pe.value(ELlKdcAPLOT[x]) for x in realtime], label = 'KgsalDH', c = 'g', ls = '-', lw = '2')
    plt.plot(realtime, [pe.value(EEcAdhZ3PLOT[x]) for x in realtime], label = 'UDH', c = 'b', ls = '-', lw = '2')
    plt.xlabel('$\it{t}$ / (min)', fontsize=17)  
    plt.ylabel('$\it{E}$$^{\mathrm{j}}$ / (μM)', fontsize=17)
    plt.grid(False)
    plt.yticks(fontsize=17)
    plt.xticks(fontsize=17)
    
    plt.tight_layout()

    pylab.savefig('ParetoOptimalPoint{}.PNG'.format(i))
    pylab.savefig('ParetoOptimalPoint{}.PDF'.format(i))
    plt.show() 
    
    i = i + 1
    w1 = w1 + 1

# printing an overview of the results of all runs 
print('BTO cascade', file = open("Overview.txt", "w+"))
print('SpaceTimeYield', file = open("Overview.txt", "a"))
print(SpaceTimeYield, file = open("Overview.txt", "a"))
print('EnzymeConsumption', file = open("Overview.txt", "a"))
print(EnzymeConsumption, file = open("Overview.txt", "a"))
print('CofactorConsumption', file = open("Overview.txt", "a"))
print(CofactorConsumption, file = open("Overview.txt", "a"))
print('IntermediateConsumption', file = open("Overview.txt", "a"))
print(IntermediateConsumption, file = open("Overview.txt", "a"))
print('FinalTime', file = open("Overview.txt", "a"))
print(FinalTime, file = open("Overview.txt", "a"))
print('InitialS1Concentration', file = open("Overview.txt", "a"))
print(InitialS1Concentration, file = open("Overview.txt", "a"))
print('InitialIntermediateConcentration', file = open("Overview.txt", "a"))
print(InitialIntermediateConcentration, file = open("Overview.txt", "a"))
print('InitialS7Concentration', file = open("Overview.txt", "a"))
print(InitialS7Concentration, file = open("Overview.txt", "a"))
print('EHsXylDH1', file = open("Overview.txt", "a"))
print(EHsXylDH1, file = open("Overview.txt", "a"))
print('EPuDHT', file = open("Overview.txt", "a"))
print(EPuDHT, file = open("Overview.txt", "a"))
print('ELlKdcA', file = open("Overview.txt", "a"))
print(ELlKdcA, file = open("Overview.txt", "a"))
print('EEcAdhZ3', file = open("Overview.txt", "a"))
print(EEcAdhZ3, file = open("Overview.txt", "a"))

print('TotalEnzymeConcentration', file = open("Overview.txt", "a"))
print(TotalEnzymeConcentration, file = open("Overview.txt", "a"))
print('info1', file = open("Overview.txt", "a"))
print(info1, file = open("Overview.txt", "a"))
print('info2', file = open("Overview.txt", "a"))
print(info2, file = open("Overview.txt", "a"))