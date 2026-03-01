import numpy as np
from data_base import coefficients

def fraction_fuel(mission,type_driven,weight_takeoff_guess,weight_changes=None):
    """
    aircraft_type: str,
        either 'jet' or 'propeller'
    number_phases: int, 
        number of phases in the mission, e.x. 3 for climb, cruise and descent.
    weight_takeoff_guess: float,
        in lb, the initial guess of the takeoff weight of the aircraft.
    weight_changes: dict, 
        with keys as 'number_phase' in order and values as weight(dropped or gained) 
        in that phase in lb.
        Format:
            weight_changes = {
            'number_phase':weight:float,
            ...
            }
    mission: dict,
        with keys as 'type_phase' in order and values as the fraction 
        of fuel used in that phase.
        Note:
            For range or loiter phase, the value is a dict with keys as 'speed', 
            'c', 'L/D', 'range' or 'endurance' and 'eta_p' for propeller case.
        Format:
            mission = {
            'type_phase':fracction_fuel_phase:float,
            'type_phase':fracction_fuel_phase:float,
            ...
            'type_phase':{  # for range or loiter phase
                'speed':float, # In mph
                'c':float, # consumption in lb/hp*hr for propeller or lb/lbf*hr for jet
                'L/D':float
                'range':float, # distance of the phase in nautical miles.
                'endurance':float, # time of the phase, if the phase is endurance
                'eta_p':float, # propeller efficiency, if the phase is propeller
                }
    """
    m_ff_0 = 1 # Initial fuel fraction
    # Calc fuel fraction for each phase.
    for phase_num,type_phase in enumerate(mission.keys(),start=1):
        if type_phase.casefold() != 'range' and type_phase.casefold() != 'endurance':
            m_ff_0 *= mission[type_phase]
        else:
            if type_driven.casefold() == 'jet': # Check type
                # Ranges case.
                if type_phase.casefold() == 'range':
                    V = mission[type_phase]['speed']
                    c = mission[type_phase]['c']
                    L_D = mission[type_phase]['L/D']
                    R = mission[type_phase]['range']
                    m_ff_range = np.exp(R*c/(V*L_D))
                    m_ff_0 *= m_ff_range
                else: # loiter case.
                    c = mission[type_phase]['c']
                    L_D = mission[type_phase]['L/D']
                    t = mission[type_phase]['endurance']
                    m_ff_endurance = np.exp(t*c/(L_D))
                    m_ff_0 *= m_ff_endurance
            else: # propeller case.
                if type_phase.casefold() == 'range':
                    c = mission[type_phase]['c']
                    L_D = mission[type_phase]['L/D']
                    R = mission[type_phase]['range']
                    eta_p = mission[type_phase]['eta_p']
                    m_ff_range = np.exp(R*c/(375*eta_p*L_D))
                    m_ff_0 *= m_ff_range
                else: # loiter case.
                    V = mission[type_phase]['speed']
                    c = mission[type_phase]['c']
                    L_D = mission[type_phase]['L/D']
                    t = mission[type_phase]['endurance']
                    eta_p = mission[type_phase]['eta_p']
                    m_ff_endurance = np.exp(t*V*c/(375*eta_p*L_D))
                    m_ff_0 *= m_ff_endurance
        if weight_changes is not None:
            if phase_num in np.asarray(list(weight_changes.keys()),dtype=int):
                initial_phase_weight = weight_takeoff_guess*m_ff_0
                final_phase_weight = initial_phase_weight - weight_changes[str(phase_num)]
                m_ff_0 *= final_phase_weight/initial_phase_weight
    return m_ff_0
def weight_empty_tentaive(weight_takeoff_guess,weight_payload,weight_crew,fuel_fractions):
    """
    Determine a tentative empty weight of the aircraft based on the fuel fractions 
    and the payload and crew weight.
    Input:
    fuel_fractions: dict,
        description
        Format:
            fuel_fractions = {
                'used':float,
                'liquids':float,
                'reserve':float
            }
    weight_takeoff_guess: float,
        in lb, the initial guess of the takeoff weight of the aircraft.
    weight_payload: float,
        in lb, the weight of the payload of the aircraft.
    weight_crew: float,
        in lb, the weight of the crew of the aircraft.
    Output:
    W_e_tentative: float,
        in lb, the tentative empty weight of the aircraft.
    """
    D = weight_payload + weight_crew
    C = fuel_fractions['used']*(1+fuel_fractions['reserve'])-fuel_fractions['liquids']-fuel_fractions['liquids']
    W_e_tentative = C*weight_takeoff_guess - D
    return W_e_tentative

if __name__ == "__main__":
    None