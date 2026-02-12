import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simpson
#   1 GAL AVIAION = 5.64 LB
#   1 GAL DE KEROSENO = 6.67 LB
def rho_h(points):
    h = np.linspace(0,18e3,points)
    rho_array = np.zeros(len(h))
    R = 287
    g = 9.80665
    t0 = 288.16 #K
    rho0 = 1.225 #kg/m3
    p0 = 1.01325e5 #Pa
    for idx,iterable in enumerate(h):
        if iterable <= 11e3:
            a = -6.5e-3 #K/m
            t = t0+a*iterable
            rho = rho0*(t/t0)**-(g/(a*R)+1)
        else:
            h1 = 11e3
            t1 = t0+a*h1
            rho1 = rho0*(t1/t0)**-(g/(a*R)+1)
            rho = rho1*np.exp(-(g/(R*t1))*(iterable-h1))
        rho_array[idx] = rho
    return rho_array*0.001941,h*3
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
def power(ac_pars,Tr,V,rho=None):
    """
    Docstring for power_required
    rho = array
    """
    Pr = Tr*V
    if "power_available" in ac_pars:
        if ac_pars["type"].casefold() == "propeller":
            if rho is not None:
                rho_ratio = rho/ac_pars["rho"]["value"]
                Pa = np.ones_like(V)*ac_pars["eta_p"]*ac_pars["power_available"]["value"]*rho_ratio # Pa in lb*ft/s
            else:
                Pa = np.ones_like(V)*ac_pars["eta_p"]*ac_pars["power_available"]["value"]
    else: ## jett
        if rho is not None:
            rho_ratio = rho/ac_pars["rho"]["value"]
            Ta = ac_pars["thrust_max"]["value"]*rho_ratio
        else:
             Ta = ac_pars["thrust_max"]["value"]
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
def rate_of_climb(ac_pars,Pa0,Pr0,V):
    """
    Doc
    """
    rho_array,h = rho_h(len(V))
    roc_array_max = []
    roc0 = (Pa0-Pr0)/ac_pars["W"]["value"]
    roc0_max = max(roc0)
    for iterable in rho_array:
        q_var = q(V,iterable)
        AE_var,_=aerodynamic_coeficients(ac_pars,q_var)
        Tr_var = thrust_required(ac_pars,AE_var)
        Pr_var,Pa_var = power(ac_pars,Tr_var,V,iterable)
        roc_var = (Pa_var-Pr_var)/ac_pars["W"]["value"]
        roc_var_max = np.max(roc_var)
        roc_array_max.append(roc_var_max)
    
    roc_array_max=np.asarray(roc_array_max,dtype=float)
    auxVar = abs(roc_array_max)
    idx = np.argmin(auxVar)
    roc_array_max = roc_array_max[:idx+1]
    roc_array_max[-1] = 1e-15
    h = h[:len(roc_array_max)]
    celling = 500 if ac_pars["type"].casefold() == "jett" else  100 
    roc_array_max = roc_array_max*60
    auxVar = abs(roc_array_max-celling)
    idx = np.argmin(auxVar)
    h_celling = h[idx]
    return roc0,roc0_max,roc_array_max,h,h_celling
def climb_time(roc_max,h):
    roc_inv = 1/(roc_max)
    h
    return roc_inv,h
def endurance_range(ac,AE_max):
    S = ac["S"]["value"]
    W = ac["W"]['value']
    Wf = ac["Wf"]["value"]
    rho = ac["rho"]["value"]
    if ac["type"].casefold() == "propeller":
       eta_p = ac['eta_p']
       c = ac['SFC']*(1/(550*3600))
       E = (eta_p/c)*AE_max['3/2']*np.sqrt(2*rho*S)*(W**(-1/2)-Wf**(-1/2))
       R = eta_p/c*AE_max['1']*np.log(W/Wf)
    else:
       c = ac['TSFC']/3600
       E = AE_max['1']/c*np.log(W/Wf)
       R = 2*np.sqrt(2/(rho*S))*AE_max['1/2']*(1/c)*(W**(1/2)-Wf**(1/2))
    print(f'Endurance: {E} seconds')
    print(f'Range: {R} ft')
    return E,R

