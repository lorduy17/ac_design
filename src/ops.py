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
    return rho_array*0.001941,h*3.28084
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
    calca escalarmente si hay ro
    """
    if rho is not None: # Return vector Power = f(h)
        rho_ratio = rho/ac_pars["rho"]['value']
        Pr = Tr*V*rho_ratio
        if ac_pars['type'] == 'jett':
            Ta = ac_pars['thrust_max']*V*rho_ratio
            Pa = Ta*V
            return Pr,Pa,Ta
        else:
            Pa = ac_pars['power_available']['value']*ac_pars['eta_p']*np.ones_like(V)*rho_ratio
        return Pr,Pa,None
    else: ## Return vector at sea level
        Pr = Tr*V
        if ac_pars['type'] == 'jett':
            Pa =  ac_pars['thrust_max']['value']*V
        else:
            Pa = ac_pars['power_available']["value"]*ac_pars['eta_p']*np.ones_like(V)
    return Pr,Pa,None
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
def rate_of_climb(ac_pars,Pa0,Pr0,V,points):
    """
    Output:
    Rate of climb at sea level = array, ft/s
    Maximun rate of climb at sea level = float, ft/s
    Maximun rate of climb in f(h) = array, ft/min
    """
    rho_array, h = rho_h(points) 
    roc_array_max = []
    W = ac_pars["W"]["value"]
    S = ac_pars["S"]["value"]
    CD0 = ac_pars["CD0"]
    # rate of climb at sea level
    roc_sl = (Pa0-Pr0)/W
    roc_sl_max = np.max(roc_sl)
    # roc_max vector calc
    if ac_pars["type"].casefold() == "jett": 
        for iterable in rho_array:
            rho_ratio = iterable/ac_pars["rho"]['value']
            q_var = q(V,iterable)
            _,AE_max=aerodynamic_coeficients(ac_pars,q_var)
            Ta = ac_pars['thrust_max']['value']*rho_ratio
            z = 1+np.sqrt(1+3/(AE_max["1"]**2*((Ta)/W)**2))
            roc_max = np.sqrt(W*z/(3*iterable*CD0*S))*(Ta/W)**(3/2)*(1-z/6-3/(2*np.square(Ta/W)*np.square(AE_max["1"])*z))
            roc_array_max.append(roc_max)
    else:
        for iterable in rho_array:
            rho_ratio = iterable/ac_pars["rho"]['value']
            q_var = q(V,iterable)
            AE,_=aerodynamic_coeficients(ac_pars,q_var)
            Tr = thrust_required(ac_pars,AE)
            Pr,Pa,_ = power(ac_pars,Tr,V,iterable)
            roc = (Pa-Pr)/W
            roc_array_max.append(np.max(roc))
    roc_array_max = np.asarray(roc_array_max,dtype=float)
    roc_array_max = roc_array_max*60
    idx = np.argmin(abs(roc_array_max))
    roc_array_max = roc_array_max[:idx+1]
    h = h[:idx+1]
    return roc_sl,roc_sl_max,roc_array_max,h
def climb_time(roc_max,h):
    roc_inv = 1/(roc_max)
    h = h
    time_trapz = np.trapezoid(roc_inv, h)
    time_simpson = simpson(roc_inv, h)
    print(f"Trapezoid: {time_trapz:.4f} min")
    print(f"Simpson: {time_simpson:.4f} min")
    print(f"Diferencia: {abs(time_trapz-time_simpson):.4f} min")
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

