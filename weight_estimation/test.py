#%%
import numpy as np
from weights_n_fuel_fractions import OperationWeight as OW


# mac = v/a; a = sqrt(gamma*R*T);
#               gamma = 1.4; R = 287; T@altitude.

w_crew = 200 # lbs
# E.X military A/C
w_crew = 175*2 # lbs
w_pl = 200*19+2000 # lbs
driven_type = 'jet'
mission ={
'start': 0.99,
'taxi': 0.99,
'takeoff': 0.99,
'climb_1': 0.971,
'cruise_1': {'speed': 459*1.15, 
        'c': 0.6, 
        'L/D': 7, 
        'range': 300-47},
'loiter_1': {'endurance': 0.5,
        'c': 0.6,
        'L/D': 9},
'descent_1': 0.99,
'cruise dash out':{'speed':400*1.15,
          'range':100,
          'c':0.9,
          'L/D':4.5},
'drop': 1,
'loiter strafe': {'endurance': 5/60,
            'c': 0.9,
            'L/D': 4.5},
'cruise dash in':{'speed':450*1.15,
          'range':100,
          'c':0.9,
          'L/D':5.5},
'climb_2': 0.969,
'cruise_4':{'speed':488*1.15,
          'c':0.6,
          'L/D':7.7,
          'range':300-47},
'descent_2': 0.99,
'landing': 0.995
}
w_to_guess = 64000

mff = OW.fuel_fraction_mission(
    mission,
    driven_type,
    weight_takeoff_guess=w_to_guess,
    weight_changes={
        '9':-20*500, # Drop payload
        '10':-2000 # Drop ammo
    }
)
fuel_fractions = {
    'total':mff,
    'reserve':0,
    'liquids':0.5/100
}
A=0.5091
B=0.9505

w_e,w_to,w_tof,C,D = OW.iterative_weight_estimation(
    w_to_guess,w_pl,w_crew,fuel_fractions,A,B)
dWto = OW.sensitivity_weights(A,B,C,D,
                           w_to)
parameter_mission_p = {
    'n':0.82,
    'c':0.5,
    'L/D':11,
    'R':1000
}
dWto_roe = OW.sensitivity_4phase(driven_type,B,D,w_to,fuel_fractions,parameter_mission_p)
# %%


# E.X 2 : propeller
payload = 200*19 # lbs
crew = 175*2 # lbs
mission ={
'start': 0.99,
'taxi': 0.995,
'takeoff': 0.995,
'climb': 0.985,
'cruise_1': {'speed': 210, 
        'c': 0.4, 
        'eta_p': 0.85,
        'L/D': 13, 
        'range': 700*0.47},
'drop':1,
'cruise_2':{'speed': 210, 
        'c': 0.4,
        'eta_p': 0.85,
        'L/D': 14, 
        'range': 700*0.53},
'loiter': {'endurance': 15/60,
        'c': 0.5,
        'eta_p': 0.77,
        'L/D': 14,
        'speed':0.7*210},
'descent': 0.985,
'finish': 0.995,
}
w_to_guess = 15800
driven_type = 'propeller'
mff = OW.fuel_fraction_mission(
    mission,
    'propeller',
    weight_takeoff_guess=w_to_guess,
    weight_changes={
        '6':-payload, # Drop payload
    }
)
fuel_fractions = {
    'total':mff,
    'reserve':0,
    'liquids':0.5/100
}
A=0.3774
B=0.9647
w_e,w_to,w_tof,C,D = OW.iterative_weight_estimation(
    w_to_guess,payload,crew,fuel_fractions,A,B)

