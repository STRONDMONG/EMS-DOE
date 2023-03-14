import math

pipesections_primary = [1,2] #number of pipe sections with same diameter in the primary heating network
pipesections_secondary = [1,2] #number of pipe sections with same diameter in the secondary heating network
circulating_pumps_primary = [1,2]
circulating_pumps_secondary = [1,2]
T = range(len(pipesections_primary))

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
t = [1,1]
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
tps = 60
tpr = 40
tg = 18
tss = 70
tsr = 50
Qlp_supply = sum((((tps-tg)*(Rb2[i]+Rt[i])-(tpr-tg)*Rc[i])/((Rb1[i]+Rt[i])*(Rb2[i]+Rt[i])-(Rc[i]**2)))*l[i] for i in T)
Qd = 0.0001*ns*a_s*(((tss+tsr)/2)-tn)**bs
print(Qd)
print(Qlp_supply)