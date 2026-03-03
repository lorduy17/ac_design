import numpy as np
from data_base import coefficients

def fuel_fraction_mission(mission:dict,
                          type_driven:str,
                          weight_takeoff_guess:float=None,
                          weight_changes:dict=None) -> float:
    """
    Calculate the fuel fraction for the mission.

    Inputs:
    -------
    type_driven: str |
        Propulsion type: {'jet','propeller'}

    number_phases: int |
        Number of phases in the missio.

    weight_takeoff_guess: float | None
        Initial guess takeoff weight [lb].

    weight_changes: dict | None
        Keys are 'number_phase' and values are weight[lb] add or dropped:
       - 'number_phase': float 

    mission: dict |
            Ordered mission profile where keys are the type of phase and
            values are:
            - float: fuel fraction for standart phases.
            - dict: parameters for cruise or loiter phas with keys:
                - speed: float [mph]
                - c: float [lb/hp*hr for propeller or lb/lbf*hr for jet]
                - L/D: float
                - cruise: float [nautical miles] (for cruise phase)
                - endurance: float [hours] (for loiter phase)
                - eta_p: float (for propeller case)

    Output:
    -------
    m_ff_0: float | Total fuel fraction for the mission.
    """
    m_ff_0 = 1 # Initial fuel fraction
    # Calc fuel fraction for each phase.
    for phase_num,type_phase in enumerate(mission.keys(),start=1):
        if type_phase.casefold() != 'cruise' and type_phase.casefold() != 'loiter':
            m_ff_0 *= mission[type_phase]
        else:
            if type_driven.casefold() == 'jet': # Check type
                # Cruise case.
                if type_phase.casefold() == 'cruise':
                    V = mission[type_phase]['speed']
                    c = mission[type_phase]['c']
                    L_D = mission[type_phase]['L/D']
                    R = mission[type_phase]['range']
                    m_ff_range = np.exp(-R*c/(V*L_D))
                    m_ff_0 *= m_ff_range
                else: # loiter case.
                    c = mission[type_phase]['c']
                    L_D = mission[type_phase]['L/D']
                    t = mission[type_phase]['endurance']
                    m_ff_endurance = np.exp(-t*c/(L_D))
                    m_ff_0 *= m_ff_endurance
            else: # propeller case.
                if type_phase.casefold() == 'cruise':
                    c = mission[type_phase]['c']
                    L_D = mission[type_phase]['L/D']
                    R = mission[type_phase]['range']
                    eta_p = mission[type_phase]['eta_p']
                    m_ff_range = np.exp(-R*c/(375*eta_p*L_D))
                    m_ff_0 *= m_ff_range
                else: # loiter case.
                    V = mission[type_phase]['speed']
                    c = mission[type_phase]['c']
                    L_D = mission[type_phase]['L/D']
                    t = mission[type_phase]['endurance']
                    eta_p = mission[type_phase]['eta_p']
                    m_ff_endurance = np.exp(-t*V*c/(375*eta_p*L_D))
                    m_ff_0 *= m_ff_endurance
        if weight_changes is not None:
            if phase_num in np.asarray(list(weight_changes.keys()),dtype=int):
                initial_phase_weight = weight_takeoff_guess*m_ff_0
                final_phase_weight = initial_phase_weight - weight_changes[str(phase_num)]
                m_ff_0 *= final_phase_weight/initial_phase_weight
    return m_ff_0
def iterative_weight_estimation(w_takeoff:float,
                                weight_payload:float,
                                weight_crew:float,
                                mission:dict,
                                type_driven:str,
                                m_ff_s:dict,
                                A:float,
                                B:float,
                                tolerance:float=0.5/100,
                                lambda_v:float=0.5) -> float:
    """
    Determine the empty weight of the aircraft sibce Roskman literature.

    Inputs:
    -------
    w_takeoff: float |
        Initial guess of takeoff [lb].

    weight_payload: float |
        Payload weight [lb].

    weight_crew: float |
        Crew weight [lb].

    mission: dict |
        Ordered mission profile where keys are the type of phase and
        values are:
        - float: fuel fraction for standart phases.
        - dict: parameters for cruise or loiter phase with keys:
            - speed: float [mph]
            - c: float [lb/hp*hr for propeller or lb/lbf*hr for jet]
            - L/D: float
            - range: float [nautical miles] (for cruise phase)
            - endurance: float [hours] (for loiter phase)
            - eta_p: float (for propeller case)
       
    type_driven: str |
        Propulsion type: {'jet','propeller'}

    m_ff_s: dict |
        Additional fuel fractions:
        - 'liquids': float
        - 'reserve': float
    
    A: float |
        Coefficient of the power law.

    B: float |
        Exponent of the power law.

    max_iter: int |
        Maximum number of iterations to perform.

    tolerance: float |
        Relative converge tolerence

    lambda_v: float |
        relaxation factor (0 < lambda_v ≤ 1)

    Output:
    -------
    w_e_allowed: float |
        Estimate empty weight of the aircraft [lb].
    """
    # calculate W_empty tentative
    m_ff = fuel_fraction_mission(mission,type_driven,weight_takeoff_guess=w_takeoff)
    C = m_ff*(1+m_ff_s['reserve'])-m_ff_s['liquids']-m_ff_s['reserve']
    D = weight_payload + weight_crew
    w_e_tent = C*w_takeoff - D
    # calculate w_empty allowed
    w_e_allowed = 10**((np.log10(w_takeoff)-A)/B)
    # Check error and iterate if necessary
    diff = (w_e_tent-w_e_allowed)
    error = abs(diff)/w_e_allowed
    iteration = 0
    while error > tolerance: # Iterate until error is within tolerance or max iterations reached
        if diff > 0:
            w_takeoff += -0.1
        else:
            w_takeoff += 0.1 # Update takeoff weight guess
        m_ff = fuel_fraction_mission(mission,type_driven,weight_takeoff_guess=w_takeoff)
        C = m_ff*(1+m_ff_s['reserve'])-m_ff_s['liquids']-m_ff_s['reserve']
        w_e_tent = C*w_takeoff - D
        w_e_allowed = 10**((np.log10(w_takeoff)-A)/B)
        diff =w_e_tent-w_e_allowed
        error = abs(diff)/w_e_allowed 
        iteration += 1
        w_tof = 10**(A+B*np.log10(C*w_takeoff-D))
        if error < tolerance:
            break
    return w_e_allowed,w_takeoff,w_tof
def sensitivity_weights(A:float,B:float,C:float,D:float,
              weight_takeoff:float):
    """
    Calculated

    Inputs:
    -------
    A: float |
        Coefficient from weight empty allowed eq. ['adim']

    B: float |
        Exponent from weight empty allowed eq. ['adim']

    C: float |
        Constat from fuel fractions. [adim]

    D: float |
        Constant from tentative eq. [lb]
    
    weight_takeoff: float |
        Weight takeoff [lb].

    Output:
    -------
    d_Wto: dict |
        Key as derivate:
        - W_pl: Weight payload
        - W_e: Weight empty
        Value float.
    """
    W_to = weight_takeoff # Abr to calcs
    d_Wto_Wpl = B*W_to/(D-C*(1-B)*W_to)
    d_Wto_We = B*W_to/(10**((np.log10(W_to)-A)/B))
    d_Wto = {
        'W_pl':d_Wto_Wpl,
        'W_e':d_Wto_We,
    }
    return d_Wto
def sensitivity_4phase(A:float,B:float,C:float,D:float,
              weight_takeoff:float,parameter_mission_p:dict):
    """
    Describe fun

    Inputs:
    -------
    A: float |
        Coefficient from weight empty allowed eq. ['adim']

    B: float |
        Exponent from weight empty allowed eq. ['adim']

    C: float |
        Constat from fuel fractions. [adim]

    D: float |
        Constant from tentative eq. [lb]
    
    weight_takeoff: float |
        Weight takeoff [lb].
    """
    ## para mff debe ser la fracción hasta esa fase en concreto? o mff total
    w_to = weight_takeoff
    mff = None
    mff_r = None
    F = -B*np.square(w_to)*(1+mff_r)*mff/(C*w_to*(1-B)-D)
    # Take mission parameters.
    try:
        c = parameter_mission_p['c']
        L_D = parameter_mission_p['L/D']
        V = parameter_mission_p['V']
        n = parameter_mission_p['n']
        E = parameter_mission_p['E']
        R = parameter_mission_p['R']
    except ValueError as e:
        raise
    propeller_partials = {
            'range':{
                'R':c/(375*n*L_D),
                'Cp':R/(375*n*L_D),
                'eta_p':-R*c/(375*L_D*n**2),
                'L/D':-R*c/(375*n*L_D**2)
            },
            'endurance':{
                'E':V*c/(375*n*L_D),
                'Cp':E*V*c/(375*n*L_D),
                'eta_p':-E*V*c/(375*L_D*n**2),
                'V':E*c/(374*n*L_D),
                'L/D':-E*V*c/(375*n*L_D**2)
            }
    }
    jet_partials = {
        'range':{
            'R':c/(V*L_D),
            'Cj':R/(V*L_D),
            'V':-R*c/(L_D*V**2),
            'L/D':-R*c/(V*L_D**2)
        },
        'endurance':{
            'E':c/(L_D),
            'Cj':E/(L_D),
            'L/D':-E*c/(L_D**2)
        }
    }

if __name__ == "__main__":
    # E.X Propeller
    w_crew = 175*6 # lbs
    w_pl = 200 # lbs
    fuel_fractions = {
        'reserve':0.25,
        'liquids':0.5/100
    }
    driven_type = 'propeller'
    mission ={
    'start': 0.992,
    'taxi': 0.996,
    'takeoff': 0.996,
    'climb': 0.99,
    'cruise': {'speed': 250, 
            'c': 0.5, 
            'L/D': 11, 
            'range': 1000, 
            'eta_p': 0.82},
    'descent': 0.99278,
    'landing': 0.992
    }
    w_to_guess = 7000
    mff = fuel_fraction_mission(
        mission,
        driven_type
    )
    w_e,w_to,w_tof = iterative_weight_estimation(w_to_guess,
                                                 w_pl,
                                                 w_crew,mission,
                                                 driven_type,
                                                 fuel_fractions,
                                                 A=0.0966,
                                                 B=1.0298)
    a = 1231231
