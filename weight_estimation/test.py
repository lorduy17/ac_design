#%%
from weights_n_fuel_fractions import OperationWeight as OW
#%%
# CLASS WORK 
w_crew = 200*9 # lbs
w_pl = 420000 # lb
driven_type = 'jet'
part_range = 6000*0.43*0.621371
mission ={
'start': 0.99,
'taxi': 0.99,
'takeoff': 0.995,
'climb_1': 0.98,
'cruise_1': {'speed': 520*1.15, 
        'c': 0.8, # tiene by pass alto TENER EN CUENTA
        'L/D': 13, # puede ser hasta más bajo porque lleva otra aeronave encima
        'range': part_range},
'drop': 1,
'cruise_2': {'range': part_range,
            'c': 0.6, # tiene by pass alto TENER EN CUENTA
            'L/D': 14,
            'speed': 520*1.15},
'loiter':{'endurance':20/60,
          'c':0.5,
          'L/D':18},
'descent_2': 0.99,
'landing': 0.99
}
w_to_guess = 9e5 # lbs
"""
valor de wto por la velocidad más que nada
"""
mff = OW.fuel_fraction_mission(
    mission,
    driven_type,
    weight_takeoff_guess=w_to_guess,
    weight_changes={
        '6':-w_pl, # Drop payload
    }
)
fuel_fractions = {
    'total':mff,
    'reserve':0,
    'liquids':0.5/100
}
A=-0.2009 # from roskman table 2.15
B=1.1037 # from roskman table 2.15
w_e,w_to,w_tof,C,D,iters = OW.iterative_weight_estimation(
    w_to_guess,w_pl,w_crew,fuel_fractions,A,B,tolerance=0.001/100)
_print = f"""
Takeoff weight by eq.: {w_tof:.3f} lbs
Takeoff weight by iterative respect w_empty.: {w_to:.3f} lbs
Empty weight: {w_e:.3f} lbs
Fuel weight: {w_to * (1 - mff):.3f} lbs
iters: {iters}

ASSUMPTIONS:
------------

Initial weight guess: {w_to_guess:.3f} lbs:
 Below check the similar aircrafts estimated takeoff weight cause:
    - Ve' more cruise speed above aircrafts.
    - Respect B747-100 SCA (Same aplication) have 280000 lbs more. So can
      be guess that the weight takeoff is higher.
        ** note: The B747-100 SCA have 2 engine less.
    - Respect Antonov An-225, the payload is 80000 lbs less. So can
      be guess that the weight takeoff is lower.
    

Cruise 1:
 Fuel consumption,  c:{mission['cruise_1']['c']} not is the max because the aircraft has a by pass ratio high,
                    so it is more efficient but in this stage
                    carry the other aircraft. So for keep the 
                    same speed need more power -> more fuel.

 L/D: {mission['cruise_1']['L/D']} Have to lower for the drag of the other aircraft.

Cruise 2:
 Fuel consumption,  c:{mission['cruise_2']['c']} is lower than cruise 1 cause' the SFA dropped the payload,
                    for keep the same speed need less power -> less fuel.

 L/D:   {mission['cruise_2']['L/D']} Increase respect cruise 1, the SFA dropped the payload, so the drag polar decrease.
        The lift is the same because the speed is the same, but the drag decrease, so the L/D increase.

Loiter:
 Fuel consumption,  c:{mission['loiter']['c']} is low cause' the A/C not need be fast, just wait for descent.
                    So the engine can be at low power, so less fuel consumption with out
                    the high by pass ratio.
 L/D:   {mission['loiter']['L/D']} For this point, in addition to the above, the airplane flight in a max
        aerodynamic efficient.
"""
print(_print)
#%%
## SENSITIVITY ANALYSIS
dWto = OW.sensitivity_weights(A,B,C,D,
                           w_to)
# new range = 7350 km
# cruise 1:
parameter_cruise_1 = {
    'c':0.8,
    'L/D':13,
    'R': part_range,
    'V': 520*1.15
}
dwto_cruise_1 = OW.sensitivity_4phase(driven_type,B,D,w_to,fuel_fractions,parameter_cruise_1)
dWto_dR_1= dwto_cruise_1['dRange']['R']
dWto_dLD_1= dwto_cruise_1['dRange']['L/D']
# cruise 2:
parameter_cruise_2 = {
    'c':0.6,
    'L/D':14,
    'R': part_range,
    'V': 520*1.15
}
dwto_cruise_2 = OW.sensitivity_4phase(driven_type,B,D,w_to,fuel_fractions,parameter_cruise_2)
dWto_dR_2= dwto_cruise_2['dRange']['R']
dWto_dLD_2= dwto_cruise_2['dRange']['L/D']
# loiter phase:
parameter_loiter = {
    'E':20/60,
    'c':0.5,
    'L/D':18,
}
dwto_loiter = OW.sensitivity_4phase(driven_type,B,D,w_to,fuel_fractions,parameter_loiter)
dWto_dE = dwto_loiter['dEndurance']['E']
dWto_dLD_loiter = dwto_loiter['dEndurance']['L/D']
delta_range = 0.43*(7350-6000)*0.621371# sm
_p = f"""
SENSITIVITY ANALYSIS
--------------------
results:

    Phase       Parameter(p)            dW_TO/dp (lbs/unit)
    --------------------------------------------------------
    Cruise 1    Range (km)              {dWto_dR_1:.3f}
                L/D                     {dWto_dLD_1:.3f}
    Cruise 2    Range (km)              {dWto_dR_2:.3f}
                L/D                     {dWto_dLD_2:.3f}
    Loiter      Endurance (hr)          {dWto_dE:.3f}
                L/D                     {dWto_dLD_loiter:.3f}
    -           Weight payload (lbs)    {dWto['W_pl']:.3f}

for:

    phase       change                   delta W_to (lbs)
    --------------------------------------------------------
    -           Payload -22%            {dWto['W_pl']*w_pl*(-0.22):.3f}
    cruise 1    range +10% (7350 km)    {dWto_dR_1*delta_range:.3f}
                L/D +18%                {dWto_dLD_1*(13*0.18):.3f}
    cruise 2    range +10% (7350 km)    {dWto_dR_2*delta_range:.3f}
                L/D +18%                {dWto_dLD_2*(14*0.18):.3f}
    loiter      reduce 15 mins          {dWto_dE*((20-15)/60):.3f}
                L/D -23%                {dWto_dLD_loiter*(-18*0.23):.3f}

"""
print(_p)
# %%