import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import root_scalar
from scipy.integrate import simpson
#   1 GAL AVIAION = 5.64 LB
#   1 GAL DE KEROSENO = 6.67 LB
def rho_h(points=None):
    h = np.linspace(0,18e3,points) if points is not None  else 1
    R = 287
    g = 9.80665
    t0 = 288.16 #K
    rho0 = 1.225 #kg/m3
    p0 = 1.01325e5 #Pa
    a = -6.5e-3
    if points is None:
        return rho0* 0.001941
    else:
        rho_array = np.zeros(len(h))
        for idx,iterable in enumerate(h):
            if iterable <= 11e3:
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
    if "AR" in ac_pars:
        AR = ac_pars["AR"]["value"]
    else:
        AR = ac_pars["b"]["value"]**2/ac_pars["S"]["value"]
    k = 1/(np.pi*AR*ac_pars["e"])
    AE_max = {
        "1":np.sqrt(CD0/k)/(2*CD0),
        "1/2":float(3*(CD0/(3*k))**(1/4)/(4*CD0)),
        "3/2":np.sqrt(CD0/(3*k))/(4/3*CD0)
    }
    if q is not None:
        CL = np.asarray((ac_pars["W"]["value"]/(q*ac_pars["S"]["value"])), dtype=np.float64)     
        CDi = k*CL**2
        CD = CD0 + CDi
        AE = {
            "1":CL/CD,
            "1/2":CL**(0.5)/CD,
            "3/2":CL**(3/2)/CD
        }
        return AE,AE_max
    return None,AE_max
def q (V,rho=None):
    """
    Docstring for q
    :param V: Description
    :param rho: Description
    """
    rho = rho_h() if rho is None else rho
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
    rho0 = rho_h(None)
    rho_ratio = 1 if rho is None else rho/rho0
    Pr = Tr*V
    if ac_pars['type'] == 'jett':
        Ta = ac_pars['thrust_max']['value']*rho_ratio
        Pa = Ta*V
        return Pr,Pa,Ta
    else:
        Pa = ac_pars['power_available']['value']*ac_pars['eta_p']*rho_ratio*np.ones_like(V)
    return Pr,Pa,None
def operation_speeds(ac_pars,AE_max,rho=None):
    CD0 = ac_pars["CD0"]
    S = ac_pars["S"]["value"]
    W = ac_pars["W"]["value"]
    rho = ac_pars['rho']['value'] if 'rho' in ac_pars else rho_h()
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
def rate_of_climb(ac_pars,V,points,rho=None):
    """
    Output:
    Rate of climb at sea level = array, ft/s
    Maximun rate of climb at sea level = float, ft/s
    Maximun rate of climb in f(h) = array, ft/min
    """

    roc_array_max = []
    W = ac_pars["W"]["value"]
    S = ac_pars["S"]["value"]
    CD0 = ac_pars["CD0"]
    _,AE_max=aerodynamic_coeficients(ac_pars,q=None)
    # rate of climb at user rho if None sea level
    rho0 = rho_h(None)
    q0 = q(V,rho0)
    AE,_=aerodynamic_coeficients(ac_pars,q0)
    Tr = thrust_required(ac_pars,AE)
    Pr0,Pa0,_=power(ac_pars,Tr,V)
    roc_sl = (Pa0-Pr0)/W
    roc_sl_max = np.max(roc_sl)
    # roc_max vector calc
    rho_array, h = rho_h(points)
    if ac_pars["type"].casefold() == "jett":
        for iterable in rho_array:
            rho_ratio = iterable/rho0
            q_var = q(V,iterable)
            Ta = ac_pars['thrust_max']['value']*rho_ratio
            z = 1+np.sqrt(1+3/(AE_max["1"]**2*((Ta)/W)**2))
            roc_max = np.sqrt(W*z/(3*iterable*CD0*S))*(Ta/W)**(3/2)*(1-z/6-3/(2*np.square(Ta/W)*np.square(AE_max["1"])*z))
            roc_array_max.append(roc_max)
    else:
        for iterable in rho_array:
            rho_ratio = iterable/rho0
            q_var = q(V,iterable)
            AE,_=aerodynamic_coeficients(ac_pars,q_var)
            Tr = thrust_required(ac_pars,AE)
            Pr,Pa,_ = power(ac_pars,Tr,V,iterable)
            roc = (Pa-Pr)/W
            roc_array_max.append(np.max(roc))
    roc_array_max = np.asarray(roc_array_max, dtype=float)
    mask = roc_array_max >= 0
    roc_array_max = roc_array_max[mask]
    h = h[mask]
    roc_array_max = roc_array_max*60
    f_ = interp1d(h,roc_array_max)
    h_celling = root_scalar(lambda x: f_(x)-100,bracket=[0, max(h)]).root
    if rho is not None: # For specific altitude
        rho_ratio = rho/rho0
        Pr_,Pa_ = Pr0*rho_ratio,Pa0*rho_ratio
        roc_ = (Pa_-Pr_)/W
        roc_max_ = np.max(roc_)
        return roc_,roc_max_
    return roc_sl,roc_sl_max,roc_array_max,h,h_celling
def climb_time(roc_max,h):
    h2 = (h >= 0) & (h <= 20e3)
    h_filt = h[h2]
    roc_fil = roc_max[h2]
    roc_inv = 1/roc_fil
    t = simpson(roc_inv, h_filt)
    print(f"Climb time: {t:.2f} minutes")
    return roc_inv, h_filt, t
def endurance_range(ac,AE_max):
    S = ac["S"]["value"]
    W = ac["W"]['value']
    Wf = ac["Wf"]["value"]
    rho = rho_h() if 'rho' not in ac else ac['rho']['value']
    if ac["type"].casefold() == "propeller":
       eta_p = ac['eta_p']
       c = ac['SFC']*(1/(550*3600))
       E = (eta_p/c)*AE_max['3/2']*np.sqrt(2*rho*S)*(Wf**(-1/2)-W**(-1/2))
       R = eta_p/c*AE_max['1']*np.log(W/Wf)
    else:
       c = ac['TSFC']/3600
       E = AE_max['1']/c*np.log(W/Wf)
       R = 2*np.sqrt(2/(rho*S))*AE_max['1/2']*(1/c)*(W**(1/2)-Wf**(1/2))
    print(f'Endurance: {E/3600} hours')
    print(f'Range: {R/5280} ft')
    return E,R
    