import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simpson 
# 

def rho_h(points):
    h = np.linspace(0,11e3,points)
    R = 287
    g = 9.80665
    t0 = 288.16 #K
    rho0 = 1.225 #kg/m3
    p0 = 1.01325e5 #Pa
    a = -6.5e-3 #K/m
    t = t0+a*h
    return rho0*(t/t0)**-(g/(a*R)+1)*0.001941
def aerodynamic_coeficients(ac_pars,q):
    AE,AE_max = {},{}
    CD0 = ac_pars["CD0"]
    CL = np.asarray((ac_pars["W"]["value"]/(q*ac_pars["S"]["value"])), dtype=np.float64)
    if "AR" in ac_pars:
        AR = ac_pars["AR"]["value"]
    else:
        AR = ac_pars["b"]["value"]**2/ac_pars["S"]["value"]
    k = 1/(np.pi*AR*ac_pars["e"])
    CDi = k*CL**2
    CD = CD0 + CDi
    AE = {
        "1":CL/CD,
        "1/2":CL**(0.5)/CD,
        "3/2":CL**(3/2)/CD
    }
    AE_max = {
        "1":np.sqrt(CD0/k)/(2*CD0),
        "1/2":np.sqrt(3*CD0/k)/(4*CD0),
        "3/2":np.sqrt(CD0/(3*k))/(4/3*CD0)
    }
    return AE,AE_max
def q (V,rho):
    """
    Docstring for q

    :param V: Description
    :param rho: Description
    """
    return 0.5 * rho * V**2
def thrust_required(ac_pars,AE):
    """
    Docstring for thrust_required
    """
    return ac_pars["W"]["value"]/AE["1"]
def power(ac_pars,Tr,V,rho):
    """
    Docstring for power_required
    rho = array
    """ 
    Pr = Tr*V
    if "Power available" in ac_pars:
        if ac_pars["type"].casefold() == "propeller":
            Pa = np.ones_like(V)*ac_pars["eta"]*ac_pars["Power available"]["value"] # Pa in lb*ft/s
    else: ## jett
        rho_0 = ac_pars["rho"]["value"]
        Ta = ac_pars["Thrust max"]["value"]*rho/rho_0
        Pa = Ta*V
    return Pr,Pa
def operation_speeds(ac_pars,AE_max):
    CD0 = ac_pars["CD0"]
    rho = ac_pars["rho"]["value"]
    S = ac_pars["S"]["value"]
    W = ac_pars["W"]["value"]
    if ac_pars["type"].casefold() == "jett":
        # V_end occurs when L/D its max -> c
        CL = AE_max["1"]*2*CD0
        V_endurance = np.sqrt(2*W/(rho*S*CL))
        # V_ra occurs whue L**(0.5)/D its max
        CL = AE_max["1/2"]*4/3*CD0
        V_range = np.sqrt(2*W/(rho*S*CL))
    else: # propeller
        # V_end occurs when L**(3/2)/D its max -> c
        CL = AE_max["3/2"]*4*CD0
        V_endurance = np.sqrt(2*W/(rho*S*CL))
        # V_ra occurs whue L/D its max
        CL = AE_max["1"]*2*CD0
        V_range = np.sqrt(2*W/(rho*S*CL))
    return V_endurance,V_range
def rate_of_climb(ac_pars,Pa,Pr,V):
    """
    Doc
    """
    rho_array = rho_h(len(V))
    RoC_array = []
    RoC_array_max = []
    for iterable in rho_array:
        q_inf = q(V,iterable)
        AE,_ = aerodynamic_coeficients(ac_pars,q_inf)
        Tr = thrust_required(ac_pars,AE)
        Pr,Pa = power(ac_pars,Tr,V,iterable)
        RoC = (Pa-Pr)/ac_pars["W"]["value"]
        roc_max = max(RoC)
        RoC_array_max.append(roc_max)
        RoC_array.append(RoC)
    RoC_array_max = np.asarray(RoC_array_max,dtype=float)
    roc = np.asarray(RoC_array[0])
    roc_max = max(roc)
    return roc,roc_max,RoC_array_max
def climb_time(ac_pars,roc,h):
    if ac_pars["type"].casefold() == "jett":
        obj_roc = 500
    else:
        obj_roc = 100
    h = h*3.2
    roc = roc*60
    idx = np.argmin(np.abs(roc-obj_roc))
    h = h[:idx+1]
    roc_inv = 1/(roc)
    roc_inv = roc_inv[:idx+1]
    t = simpson(roc_inv, h)
    return roc_inv,h,t,idx

#%%