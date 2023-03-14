import pyomo.environ as pyomo
import matplotlib.pyplot as plt
import math
pipesections_primary = [1,2] #number of pipe sections with same diameter in the primary heating network
pipesections_secondary = [1,2] #number of pipe sections with same diameter in the secondary heating network
circulating_pumps_primary = [1,2]
circulating_pumps_secondary = [1,2]

##Insulation parameters
dw = 0.291
dz = 0.355
lambdab = 0.03 
Lambdat = 2.5 #https://www.researchgate.net/figure/Soil-characterisation-of-a-site-in-Ostend-Belgium-showing-a-thermal-conductivity_fig18_321485524
H_var = 2
h = H_var - dz/2
b = 0.150 
## Insulation calculations
def Rb_cal(lambdab, dz, dw):
    Rb = (1/(2*math.pi*lambdab))*math.log(dz/dw)
    return Rb

def Rt_cal(lambdat, H, dz):
    Rt = (1/(2*math.pi*lambdat))*math.log((4*H)/dz)
    return Rt

def H_cal(h, lambdat, alpha):
    H = h + (lambdat/alpha)
    return H

def Rc_cal(lambdat, H, b):
    Rc = (1/(2*math.pi*lambdat))*math.log(math.sqrt(1+((2*H)/b)**2))
    return Rc

m = pyomo.ConcreteModel()
m.cons = pyomo.ConstraintList()
T = range(len(pipesections_primary))
#sets
m.n = pyomo.Set(initialize=range(len(pipesections_primary)))
m.m = pyomo.Set(initialize=range(len(pipesections_secondary)))
m.npp =  pyomo.Set(initialize=range(len(circulating_pumps_primary))) #number of circulating pumps of the primary side
m.nsp= pyomo.Set(initialize=range(len(circulating_pumps_secondary))) #number of circulating pumps of the secondary side

#Data parameters
tg = 18
Rb1 = [Rb_cal(lambdab,dz,dw),Rb_cal(lambdab*0.8,dz,dw)]
Rb2 = [Rb_cal(lambdab,dz,dw),Rb_cal(lambdab*0.8,dz,dw)]
Rt=[Rt_cal(Lambdat,H_var,dz), Rb_cal(lambdab*0.8,dz,dw)]
Rc= [Rc_cal(Lambdat,H_var,b), Rb_cal(lambdab*0.8,dz,dw)]
l = [200,180]
beta = 0.5
Ue = 0.45
Nr = [77.6, 100]
ne = [0.95,0.94]
nf = [0.96,0.98]
Qs = [184,140]
Gr = [170,160]
c = 4.182
Uh = 0.08
a_s = 5
Kh = 2000 #https://www.engineeringtoolbox.com/heat-transfer-coefficients-exchangers-d_450.html
Ah = 4
dtm = 59.9 #https://www.thermal-engineering.org/what-is-logarithmic-mean-temperature-difference-lmtd-definition/
tps_min = 60
tps_max = 120

tpr_min = 30
tpr_max = 70

tss_min = 30
tss_max = 85

tsr_min = 30
tsr_max = 60
tn = 21
ns = 1
alphas = 1 
bs = 1

#Variables
m.tps = pyomo.Var(domain=pyomo.Reals)
m.tpr = pyomo.Var(domain=pyomo.Reals)
m.tsr = pyomo.Var(domain=pyomo.Reals)
m.tss = pyomo.Var(domain=pyomo.Reals)
m.Cp = pyomo.Var(domain=pyomo.Reals)
m.Cl = pyomo.Var(domain=pyomo.Reals)
m.Qd = pyomo.Var(domain=pyomo.Reals)
m.Qlp_supply = pyomo.Var(domain=pyomo.Reals)
m.Qlp_return = pyomo.Var(domain=pyomo.Reals)
m.Qls_supply = pyomo.Var(domain=pyomo.Reals)
m.Qls_return = pyomo.Var(domain=pyomo.Reals)
m.Qlp= pyomo.Var(domain=pyomo.Reals)
m.Qls= pyomo.Var(domain=pyomo.Reals)
m.Qs= pyomo.Var(domain=pyomo.Reals)
m.Cp_primary= pyomo.Var(domain=pyomo.Reals)
m.Cp_secondary= pyomo.Var(domain=pyomo.Reals)
m.Cl_primary = pyomo.Var(domain=pyomo.Reals)
m.Cl_secondary= pyomo.Var(domain=pyomo.Reals)
#Objective
cost = m.Cp+m.Cl
cost = sum(sum((g.variable_costs_el+g.fuelpricing[t])*(m.P_CHPCOALht[t,g]+m.P_CHPCOALel[t,g]) - m.P_CHPCOALel[t,g]*Pel[t] +m.Cramp[t,g]*g.rampcost for g in m.plants) for t in m.time)

#Model
m.cost = pyomo.Objective(expr = cost, sense=pyomo.minimize)
m.cons.add(m.Qd == 0.0001*ns*a_s*((((m.tss+m.tsr)/2)-tn)**bs))

m.cons.add(m.Qlp_supply == sum((((m.tps-tg)*(Rb2[i]+Rt[i])-(m.tpr-tg)*Rc[i])/((Rb1[i]+Rt[i])*(Rb2[i]+Rt[i])-(Rc[i]**2)))*l[i] for i in m.n))
m.cons.add(m.Qlp_return == sum((((m.tpr-tg)*(Rb1[i]+Rt[i])-(m.tps-tg)*Rc[i])/((Rb1[i]+Rt[i])*(Rb2[i]+Rt[i])-(Rc[i]**2)))*l[i] for i in m.n))
m.cons.add(m.Qls_supply == sum((((m.tss-tg)*(Rb2[j]+Rt[j])-(m.tsr-tg)*Rc[j])/((Rb1[j]+Rt[j])*(Rb2[j]+Rt[j])-(Rc[j]**2)))*l[j] for j in m.m))
m.cons.add(m.Qls_return == sum((((m.tsr-tg)*(Rb1[j]+Rt[j])-(m.tss-tg)*Rc[j])/((Rb1[j]+Rt[j])*(Rb2[j]+Rt[j])-(Rc[j]**2)))*l[j] for j in m.m))
m.cons.add(m.Qlp == 10**(-3)*(1+beta)*(m.Qlp_supply+m.Qlp_return))
m.cons.add(m.Qls == 10**(-3)*(1+beta)*(m.Qls_supply+m.Qls_return))
m.cons.add(m.Qs== m.Qd+m.Qlp+m.Qls)
m.cons.add(Kh*Ah*dtm == m.Qd + m.Qls)
m.cons.add(m.Cp_primary == sum(((Ue*Nr[k])/(3600.0*ne[k]*nf[k]))*((m.Qs/(Gr[k]*c*(m.tps-m.tpr)))**3.0) for k in m.npp))
m.cons.add(m.Cp_secondary == sum(((Ue*Nr[q])/(3600*ne[q]*nf[q]))*((m.Qs/(Gr[q]*c*(m.tss-m.tsr)))**3) for q in m.nsp))
m.cons.add(m.Cp == m.Cp_primary + m.Cp_secondary)
m.cons.add(m.Cl_primary == ((10**(-3)*Uh*(1+beta))/3600)*(m.Qlp_supply+m.Qlp_return))
m.cons.add(m.Cl_secondary == ((10**(-3)*Uh*(1+beta))/3600)*(m.Qls_supply+m.Qls_return))
m.cons.add(m.Cl == m.Cl_primary+m.Cl_secondary)

m.cons.add(tps_min <= m.tps)
m.cons.add(tps_max >= m.tps)

m.cons.add(tpr_min <= m.tpr)
m.cons.add(tpr_max >= m.tpr)

m.cons.add(tss_min <= m.tss)
m.cons.add(tss_max >= m.tss)

m.cons.add(tsr_min <= m.tsr)
m.cons.add(tsr_max >= m.tsr)

m.cons.add(m.tsr <= m.tpr)
m.cons.add(m.tss <= m.tps)
m.cons.add(m.tpr <= m.tps)
m.cons.add(m.tsr <= m.tss)
m.cons.add(tn <= m.tsr)

solver = pyomo.SolverFactory('ipopt').solve(m, tee=True, options="halt_on_ampl_error=yes")

for v in m.component_objects(pyomo.Var, active=True):
    print ("Variable component object",v)
    print ("Type of component object: ", str(type(v))[1:-1]) # Stripping <> for nbconvert
    varobject = getattr(m, str(v))
    print ("Type of object accessed via getattr: ", str(type(varobject))[1:-1])