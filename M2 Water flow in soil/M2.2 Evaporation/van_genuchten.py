""" Mualem - van Genuchten model, equations from van Genuchten, MT (1980) """
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from scipy import optimize

# class containing the van Genuchten parameters
class Parameters:
    def __init__(self, R, S, alpha, n, Ksat, l = 0.5):
        self.theta_R = R
        self.theta_S = S        
        self.alpha = alpha # [1/cm]         
        self.n = n
        self.m = 1.-1./n
        self.Ksat = Ksat   
        self.l = l     

# returns the volumetric water content at a given pressure head  according to the van genuchten model (Eqn 21)
def water_content(h, sp):
    return sp.theta_R + (sp.theta_S-sp.theta_R)/pow(1. + pow(sp.alpha*abs(h),sp.n),sp.m)

# returns pressure head at a given volumetric water content according to the van genuchten model
def pressure_head(theta, sp): 
    theta = min(theta,sp.theta_S) # saturated water conent is the maximum 
    return - pow( pow( (sp.theta_S - sp.theta_R)/(theta - sp.theta_R), (1./sp.m))-1., 1./sp.n) / sp.alpha

# returns the effective saturation according to the van genuchten model (dimensionless water content, Eqn 2)
def effective_saturation(h,sp):
    h = min(h,0) # pressure head is negative, zero the maximum
    theta = water_content(h,sp)
    se = (theta-sp.theta_R)/(sp.theta_S-sp.theta_R)
    return se

# returns the hydraulic conductivity according to the van genuchten model (Eqn 8)
def hydraulic_conductivity(h,sp):
    se = effective_saturation(h,sp) 
    K = sp.Ksat*np.sqrt(se)*( (1. - pow(1. - pow(se, 1. / sp.m),sp.m)) ** 2 )
    return K 

# returns the specific moisture storage according to the van genuchten model
def specific_moisture_storage(h,sp):
    C = -sp.alpha*sp.n*np.sign(h)*(1. / sp.n - 1.) * pow(sp.alpha*abs(h), sp.n-1.) * (sp.theta_R-sp.theta_S) * pow(pow(sp.alpha*abs(h),sp.n) + 1., 1./sp.n-2.)
    return C

# returns the water diffusivity (Eqn 11)
def water_diffusivity(TH, theta_i, theta_sur, sp):
    theta = TH * (theta_i - theta_sur) + theta_sur
    Se = (theta - sp.theta_R) / (sp.theta_S - sp.theta_R)
    m = sp.m
    D = (1 - m) * sp.Ksat / (sp.alpha * m * (sp.theta_S - sp.theta_R)) * pow(Se, 0.5 - 1. / m) * (pow(1 - pow(Se, 1. / m), -m) + pow(1 - pow(Se, 1 / m), m) - 2)
    return D

# returns the matric flux potential
def MFP(h,sp):
    K = lambda h: hydraulic_conductivity(h,sp) # integrand 
    MFP, err = integrate.quad(K,-15000, h)
    return MFP

# returns the matric potential from matric flux potential
def h(MFP_given,sp):
    MFP_root = lambda psi: MFP(psi,sp)-MFP_given
    h = optimize.brentq(MFP_root, -15000, 0)
    return h


# normalized root mean squared error
def nRMSE(y_i, y_hat):
    assert len(y_i)==len(y_hat), "number of observations y_i must equalt number of predictions y_hat"
    n = len(y_i)
    rmse = np.sqrt(np.sum(np.square(y_i-y_hat))/n)
    return rmse/np.abs(np.mean(y_i))

# Nash–Sutcliffe model efficiency coefficient
def nNSE(y_i, y_hat):
    assert len(y_i)==len(y_hat), "number of observation times must equalt number of predictions times"
    a = np.sum(np.square(y_i-y_hat))
    b = np.sum(np.square(y_i-np.mean(y_i) ))
    #return 1-a/b # NSE
    return 1/(2-(1-a/b)) # NNSE
    
    
""" settings for all plots """
    
plt.rc('font', size=18)          # controls default text sizes
plt.rc('axes', titlesize=18)     # fontsize of the axes title
plt.rc('axes', labelsize=18)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=18)    # fontsize of the tick labels
plt.rc('ytick', labelsize=18)    # fontsize of the tick labels
plt.rc('legend', fontsize=16)    # legend fontsize
plt.rc('figure', titlesize=18)  # fontsize of the figure title

# Color Scheme

# cc_inc = 5 # increment
# cmap = plt.get_cmap('tab20c')
# col = 0.75*cmap(range(0,100))

#cc_inc = 1 # increment
#col=np.array([(1,0.,0.),(0.,0.,1),(0.,1,0.),(1,0.,1),(0.5,1,1)])
#col = np.vstack((col, 0.5*col))

# cmap = plt.get_cmap('Dark2')
# cmap = plt.get_cmap('Accent') # ganz schlecht
cmap = plt.get_cmap('tab10')
col = cmap([0,1,2,4,5,6,7,8,9,10,11,12,13,14,15]) # avoid 'red'



