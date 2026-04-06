import numpy as np
from data_base import coefficients
class OperationWeight:
    def isa_conditions(h=None):
        R = 287 # J/kg*K
        g = 9.80665 # m/s2
        t0 = 288.16 #K
        rho0 = 1.225 #kg/m3
        a = -6.5e-3 # K/m
        if h is None or h == 0:
            return rho0 * 0.001941, t0
        else: # Gradient
            if h <=11e3:
                t = t0+a*h
                rho = rho0*(t/t0)**-(g/(a*R)+1)
            else: # Isothermal
                h1 = 11e3
                t1 = t0+a*h1
                rho1 = rho0*(t1/t0)**-(g/(a*R)+1)
                rho = rho1*np.exp(-(g/(R*t1))*(h-h1)) 
        return rho*0.001941,t1 # rho [slugs/ft3], t [K]
    @staticmethod
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
        # dict for store weights ratio if weight changes are given.
        w_rs = {}
        var_aux= []
        m_ff_0 = 1 # Initial fuel fraction
        # Calc fuel fraction for each phase.
        for phase_num,type_phase in enumerate(mission.keys(),start=1):
            if 'cruise' not in type_phase.casefold() and 'loiter' not in type_phase.casefold():
                m_ff_0 *= mission[type_phase]
            else:
                if type_driven.casefold() == 'jet': # Check type
                    # Cruise case.
                    if 'cruise' in type_phase.casefold():
                        V = mission[type_phase]['speed']
                        c = mission[type_phase]['c']
                        L_D = mission[type_phase]['L/D']
                        R = mission[type_phase]['range']
                        m_ff_range = np.exp(-R*c/(V*L_D))
                        m_ff_0 *= m_ff_range
                    # loiter case.
                    elif 'loiter' in type_phase.casefold(): 
                        c = mission[type_phase]['c']
                        L_D = mission[type_phase]['L/D']
                        t = mission[type_phase]['endurance']
                        m_ff_endurance = np.exp(-t*c/(L_D))
                        m_ff_0 *= m_ff_endurance
                else: # propeller case.
                    if 'cruise' in type_phase.casefold():
                        c = mission[type_phase]['c']
                        L_D = mission[type_phase]['L/D']
                        R = mission[type_phase]['range']
                        eta_p = mission[type_phase]['eta_p']
                        m_ff_range = np.exp(-R*c/(375*eta_p*L_D))
                        m_ff_0 *= m_ff_range
                    elif 'loiter' in type_phase.casefold(): # loiter case.
                        V = mission[type_phase]['speed']
                        c = mission[type_phase]['c']
                        L_D = mission[type_phase]['L/D']
                        t = mission[type_phase]['endurance']
                        eta_p = mission[type_phase]['eta_p']
                        m_ff_endurance = np.exp(-t*V*c/(375*eta_p*L_D))
                        m_ff_0 *= m_ff_endurance
            if weight_changes is not None:
                if str(phase_num) in weight_changes.keys(): # Check if the current phase have weight change.
                    initial_weight = weight_takeoff_guess*m_ff_0
                    final_weight = initial_weight + weight_changes.get(str(phase_num),0)
                    weight_ratio = final_weight/initial_weight # These value is used in the next iteration
                    w_rs[str(phase_num)] = weight_ratio
                    var_aux.append(phase_num)
            # Check if the current phase need correction by weight ratio
            # for the fuel fraction calculation.
            if len(var_aux) > 0:
                if phase_num-1 == var_aux[-1]:
                    # Check actual case
                    if 'cruise' in type_phase.casefold():
                        m_ff_save = m_ff_range
                        m_ff_0 /= m_ff_range
                    elif 'loiter' in type_phase.casefold():
                        m_ff_save = m_ff_endurance
                        m_ff_0 /= m_ff_endurance
                    m_ff_c = (1-(1-m_ff_save)*w_rs[str(phase_num-1)])
                    m_ff_0 *= m_ff_c
        return m_ff_0
    def iterative_weight_estimation(w_takeoff:float,
                                    weight_payload:float,
                                    weight_crew:float,
                                    m_ff_s:dict,
                                    A:float,
                                    B:float,
                                    tolerance:float=0.1/100) -> float:
        """
        Determine the aircraft empty weight.

        Inputs:
        -------
        w_takeoff: float |
            Initial guess of takeoff [lb].

        weight_payload: float |
            Payload weight [lb].

        weight_crew: float |
            Crew weight [lb]

        m_ff_s: dict |
            Fuel fractions:
            - 'total': float
            - 'liquids': float
            - 'reserve': float
        
        A: float |
            Coefficient of the power law.

        B: float |
            Exponent of the power law.

        tolerance: float |
            Relative converge tolerence

        Output:
        -------
        w_e_allowed: float |
            Estimate empty weight of the aircraft [lb].
        """
        # calculate W_empty tentative
        m_ff0,m_ff_r,m_ff_l = m_ff_s['total'],m_ff_s['reserve'],m_ff_s['liquids']
        C = m_ff0*(1+m_ff_r)-m_ff_l-m_ff_r
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
            m_ff0,m_ff_r,m_ff_l = m_ff_s['total'],m_ff_s['reserve'],m_ff_s['liquids']
            C = m_ff0*(1+m_ff_r)-m_ff_l-m_ff_r 
            w_e_tent = C*w_takeoff - D
            w_e_allowed = 10**((np.log10(w_takeoff)-A)/B)
            diff =w_e_tent-w_e_allowed
            error = abs(diff)/w_e_allowed 
            iteration += 1
            w_tof = 10**(A+B*np.log10(C*w_takeoff-D))
            if error < tolerance:
                break
        return w_e_allowed,w_takeoff,w_tof,C,D,iteration
   