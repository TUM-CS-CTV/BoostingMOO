from pyomo.environ import *
from pyomo.dae import *

def FunctionforPF(ICi,tfi,S1i,S2i,S3i,S4i,S5i,S7i,EHsXylDH1i,EPuDHTi,ELlKdcAi,
                  EEcAdhZ3i): 
    
    
    model = ConcreteModel()


    # total batch running time (min)
    model.tf = Var(bounds = (10,10000), initialize = tfi)     
    # scaled batch running time (min/min)                 
    model.tau = ContinuousSet(bounds=(0,1)) 
    # the reaction rates (mM/min)
    model.vI = Var(model.tau, within=NonNegativeReals)
    model.vII = Var(model.tau, within=NonNegativeReals)
    model.vIII = Var(model.tau, within=NonNegativeReals)
    model.vIV = Var(model.tau, within=NonNegativeReals)
    model.vV = Var(model.tau, within=NonNegativeReals)
    # the substrate concentrations (mM)
    model.S1 = Var(model.tau, bounds = (0, 1000), initialize = S1i) 
    model.S2 = Var(model.tau, bounds = (0, 1000), initialize = S2i) 
    model.S3 = Var(model.tau, bounds = (0, 1000), initialize = S3i) 
    model.S4 = Var(model.tau, bounds = (0, 1000), initialize = S4i) 
    model.S5 = Var(model.tau, bounds = (0, 1000), initialize = S5i) 
    model.S6 = Var(model.tau, bounds = (0, 1000))                   
    model.S7 = Var(model.tau, bounds = (0, 500), initialize = S7i)
    model.S8 = Var(model.tau, bounds = (0, 500))
    # the enzyme concentrations (μM)
    model.EHsXylDH1 = Var(bounds = (0, 1000), initialize = EHsXylDH1i)
    model.EPuDHT = Var(bounds = (0, 1000), initialize = EPuDHTi)
    model.ELlKdcA = Var(bounds = (0, 1000), initialize = ELlKdcAi)
    model.EEcAdhZ3 = Var(bounds = (0, 1000), initialize = EEcAdhZ3i)
    # the first order derivatives of the substrate concentrations (mM/min)
    model.dS1dt = DerivativeVar(model.S1, wrt=model.tau, within=Reals)
    model.dS2dt = DerivativeVar(model.S2, wrt=model.tau, within=Reals)
    model.dS3dt = DerivativeVar(model.S3, wrt=model.tau, within=Reals)
    model.dS4dt = DerivativeVar(model.S4, wrt=model.tau, within=Reals)
    model.dS5dt = DerivativeVar(model.S5, wrt=model.tau, within=Reals)
    model.dS6dt = DerivativeVar(model.S6, wrt=model.tau, within=Reals)
    model.dS7dt = DerivativeVar(model.S7, wrt=model.tau, within=Reals)
    model.dS8dt = DerivativeVar(model.S8, wrt=model.tau, within=Reals)
    # the total concentration of all enzymes used during the batch (μΜ)
    model.SumEnzymes = Var(within=NonNegativeReals)
    # the total concentration of all intermediates added at tf=0 (mM)
    model.SumIntermediates = Var(within=NonNegativeReals)
    # the yield (mM/mM)
    model.Yield = Var(within=Reals)
    # the cofactor consumption (mM/min)
    model.CC = Var(within=NonNegativeReals)
    # the enzyme consumption (μΜ/min)
    model.EC = Var(within=NonNegativeReals)
    # the intermediate consumption (mM/min)
    model.IC = Var(within=NonNegativeReals, initialize = ICi)
    # the objective (Space-time yield) (mM/min)
    model.OBJ = Var(model.tau, within=NonNegativeReals)
    
    
    model.L = Set(initialize = ['EcAdhZ3','HsXylDH1','PuDHT','LlKdcA'])
    model.M = Set(initialize = [1,2,3,4,5,6,7,8,9])

    # the molecular weights of all enzymes (mg/mM)
    mw = {}
    mw['HsXylDH1']      =  31000
    mw['PuDHT']         =  70000
    mw['LlKdcA']        =  70000
    mw['EcAdhZ3']       =  45000
    model.mw = Param(model.L, initialize = mw)

    # the maximum reaction rates of all enzyme catalyzed reactions (U/mg)
    Vmax = {}
    Vmax['HsXylDH1']    = 20.4
    Vmax['PuDHT']       = 11.9
    Vmax['LlKdcA']      = 0.9
    Vmax['EcAdhZ3']     = 12      
    model.Vmax = Param(model.L, initialize = Vmax)
    
    # the kinetic parameters of all enzyme catalyzed reactions (mM)
    K = {}
    K[1]      = 0.2
    K[2]      = 1.4 
    K[3]      = 524.5
    K[4]      = 0.5
    K[5]      = 7.5 
    K[6]      = 16.7 
    K[7]      = 1.5 
    model.K = Param(model.M, initialize = K)
    
    # the epsilon constraint for the intermediate consumption (mΜ/min)
    model.ICc = Param(initialize = ICi)
    
    # the reaction rate kinetics     
    def rr1(m, tau):
        return model.vI[tau] == (model.Vmax['HsXylDH1']*model.S7[tau])/(model.K[1]+model.S7[tau]*(1+model.K[2]/model.S1[tau])*(1+model.S1[tau]/model.K[3]))*model.EHsXylDH1*model.mw['HsXylDH1']*10**(-6)
    model.rr1con = Constraint(model.tau, rule=rr1)    
    
    def rr2(m, tau):
        return model.vII[tau] == 10000*model.S2[tau] 
    model.rr2con = Constraint(model.tau, rule=rr2)
    
    def rr3(m, tau):   
        return model.vIII[tau] == (model.Vmax['PuDHT']*model.S3[tau])/ \
            (model.K[4]+model.S3[tau])*model.EPuDHT*model.mw['PuDHT']*10**(-6)
    model.rr3con = Constraint(model.tau, rule=rr3)
    
    def rr4(m, tau): 
        return model.vIV[tau] == (model.Vmax['LlKdcA']*model.S4[tau])/ \
            (model.K[5]+model.S4[tau])*model.ELlKdcA*model.mw['LlKdcA']*10**(-6)
    model.rr4con = Constraint(model.tau, rule=rr4)
    
    def rr5(m, tau):  
        return model.vV[tau] == (model.Vmax['EcAdhZ3']*model.S5[tau]*model.S8[tau])/(model.K[6]*model.K[7]+model.K[7]*model.S5[tau]+model.S5[tau]*model.S8[tau])*model.EEcAdhZ3*model.mw['EcAdhZ3']*10**(-6)
    model.rr5con = Constraint(model.tau, rule=rr5)
    
    
    # the material balances for all substrates 
    def d1(m, tau):
        return model.dS1dt[tau] / model.tf == -model.vI[tau]
    model.d1con = Constraint(model.tau, rule=d1)
    
    def d2(m, tau):  
        return model.dS2dt[tau] / model.tf == +model.vI[tau]-model.vII[tau]
    model.d2con = Constraint(model.tau, rule=d2)
    
    def d3(m, tau):
        return model.dS3dt[tau] / model.tf == +model.vII[tau]-model.vIII[tau]
    model.d3con = Constraint(model.tau, rule=d3)
    
    def d4(m, tau):
        return model.dS4dt[tau] / model.tf == +model.vIII[tau]-model.vIV[tau]
    model.d4con = Constraint(model.tau, rule=d4)
    
    def d5(m, tau):
        return model.dS5dt[tau] / model.tf == +model.vIV[tau]-model.vV[tau]
    model.d5con = Constraint(model.tau, rule=d5)
    
    def d6(m, tau):
        return model.dS6dt[tau] / model.tf == +model.vV[tau]
    model.d6con = Constraint(model.tau, rule=d6)
    
    def d7(m, tau):
        return model.dS7dt[tau] / model.tf == -model.vI[tau]+model.vV[tau]
    model.d7con = Constraint(model.tau, rule=d7)
    
    def d8(m, tau):
        return model.dS8dt[tau] / model.tf == +model.vI[tau]-model.vV[tau]   
    model.d8con = Constraint(model.tau, rule=d8)
    
    # calculation of the total enzyme concentration used 
    def c1(m):
        return model.SumEnzymes == model.EHsXylDH1 + model.EPuDHT + model.ELlKdcA + model.EEcAdhZ3
    model.c1con = Constraint(rule=c1) 
        
    # calculation of the total intermediate concentration added at t=0
    def c2(m):
        return model.SumIntermediates == model.S2[0] + model.S3[0] + model.S4[0] + model.S5[0]
    model.c2con = Constraint(rule=c2)
    
    # definition of the yield 
    def c3(m):
        return model.Yield == (model.S6[1]-model.S6[0])/(model.S1[0])
    model.c3con = Constraint(rule=c3)
    
    # a constraint on the yield
    def c4(m):
        return model.Yield >= 0.95
    model.c4con = Constraint(rule=c4)
    
    # definition of the cofactor consumption 
    def c5(m):
        return model.CC == (model.S7[0] + model.S8[0]) / (model.tf+30)
    model.c5con = Constraint(rule=c5)
    
    # definition of the enzyme consumption
    def c6(m):
        return model.EC == model.SumEnzymes / (model.tf+30) 
    model.c6con = Constraint(rule=c6)
    
    # definition of the intermediate consumption
    def c7(m):
        return model.IC == model.SumIntermediates / (model.tf+30)
    model.c7con = Constraint(rule=c7)
    
    # constraint on the cofactor consumption (changed manually)
    def c8(m):
        return model.CC <= 0.0005
    model.c8con = Constraint(rule=c8)
    
    # constraint on the enzyme consumption (changed manually)
    def c9(m, tau):
        return model.EC <= 0.1
    model.c9con = Constraint(model.tau, rule=c9)
    
    # constraint on the intermediate consumption (changed automatically)
    def c10(m, tau):
        return model.IC <= model.ICc
    model.c10con = Constraint(model.tau, rule=c10)
    
    # definition of the objective (space-time yield)
    def OBJ(m, tau):
        return model.OBJ[tau] == model.S6[tau]/(model.tf+30)
    model.OBJcon = Constraint(model.tau, rule=OBJ)
    
    # initial values for the substrates not added to the reactor
    model.ic = ConstraintList()
    # model.ic.add(model.S2[model.tau.first()] == 0.001)
    # model.ic.add(model.S3[model.tau.first()] == 0.001)
    # model.ic.add(model.S4[model.tau.first()] == 0.001)
    # model.ic.add(model.S5[model.tau.first()] == 0.001)
    model.ic.add(model.S6[model.tau.first()] == 0.001)
    model.ic.add(model.S8[model.tau.first()] == 0.001)
    #model.ic.add(model.EHsXylDH1 == 15.70
    #model.ic.add(model.EPuDHT    == 3.41)
    #model.ic.add(model.ELlKdcA   == 46.75)
    #model.ic.add(model.EEcAdhZ3  == 9.14)

    # selection of the objective for pyomo
    model.obj = Objective(expr=model.OBJ[1], sense=maximize)
    # selection of a discretization method
    discretizer = TransformationFactory('dae.finite_difference')
    # selection of the number of finite elements and discretization options 
    discretizer.apply_to(model, wrt=model.tau, nfe=50, scheme='BACKWARD')
   
    return model