import numpy as np

def _(weight_payload,weight_crew,):
    D = weight_payload + weight_crew
def weight_fuel(mission,type_driven,weight_takeoff_guess,weight_drop=None):
    """
    aircraft_type: str,
        either 'jet' or 'propeller'
    number_phases: int, 
        number of phases in the mission, e.x. 3 for climb, cruise and descent.
    weight_takeoff_guess: float,
        in lb, the initial guess of the takeoff weight of the aircraft.
    weight_drop: dict, 
        with keys as 'number_phase' in order and values as t weight drop in that phase.
        Format:
            weight_drop = {
            'number_phase':weight_drop_phase:float,
            ...
            }
    mission: dict,
        with keys as 'type_phase' in order and values as the fraction of fuel used in that phase, e.g. 0.2 for 20% fuel used in climb phase. For range or loiter phase, the value is a dict with keys as 'speed', 'c', 'L/D', 'range' or 'endurance' and 'eta_p' for propeller case.
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
    if type_driven.casefold() == 'jet': # Check type
        # Ranges case.
        if mission['type_phase'].casefold() == 'range':
            V = mission['type_phase']['speed']
            c = mission['type_phase']['c']
            L_D = mission['type_phase']['L/D']
            R = mission['type_phase']['range']
            m_ff_range = np.exp(R*c/(V*L_D))
        else: # loiter case.
            V = mission['type_phase']['speed']
            c = mission['type_phase']['c']
            L_D = mission['type_phase']['L/D']
            t = mission['type_phase']['endurance']
            m_ff_endurance = np.exp((V/c)*L_D*t)
    else:
        if mission['type_phase'].casefold() == 'range':
            c = mission['type_phase']['c']
            L_D = mission['type_phase']['L/D']
            R = mission['type_phase']['range']
            eta_p = mission['type_phase']['eta_p']
            m_ff_range = np.exp(R*c/(eta_p*L_D))
    m_ff_0 = 1
    for i in mission:
        if mission[i].casefold() != 'range' and mission[i].casefold() != 'endurance':
            m_ff_0 *= mission[i]
    m_ff_0 = m_ff_0*m_ff_range
    mff = m_ff_0*m_ff_endurance if mission['type_phase'].casefold() == 'endurance' else m_ff_0
    weight_fuel_used = weight_takeoff_guess*(1-mff)
    weight_fuel_reserve = 0.25*weight_fuel_used
    weight_fuel_total = weight_fuel_used + weight_fuel_reserve
    return mff, weight_fuel_total