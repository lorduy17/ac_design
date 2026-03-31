import numpy as np

class mc_ops:
    ### MATH FUNCTIONS WORKS IN ENGLISH UNITS
    @staticmethod
    def stall_speed(vs,cl_max):
        """
        Calc required value of wing loading by stall speed requirement.
        """
        rho = 0.00237 # slugs/ft^3
        wing_load = 0.5*rho*cl_max*vs**2
        return wing_load
    @staticmethod # PUEDE CAMBIAR PORQUE PUEDE QUE AÑADA UN SELF (esp. para mantener el far, tipo de motor, etc)
    def takeoff_distance(
        FAR:int,
        distances:dict,
        cl_max:int,
        dens_ratio:int,
        wing_load:np.asarray,
        driven_by_jett=True
    )->np.asarray:
        """
        Calcs array of power loading or thrust-weight by distances required.

        Parameters
        ----------
            FAR : int
                Federal Aviation Regulation (FAR) part number, either 23 or 25.
                
            distances : dict
                A dictionary containing the relevant distances for the takeoff calculation. 
                - For FAR 23, it should include 's_tog' (takeoff ground roll distance) and 's_to' (total takeoff distance).
                - For FAR 25, it should include 's_tofl' (takeoff field length).
            cl_max : int
                Maximum lift coefficient of the aircraft.
            dens_ratio : int
                Density ratio of the air at the takeoff conditions compared to sea level standard conditions.
            wing_load : np.asarray
                Wing loading of the aircraft, typically in pounds per square foot.
            driven_by_jett : bool, optional
                A flag indicating whether the aircraft is driven by jet engines (True) or turboprops (False). 

                This affects the thrust-to-weight ratio calculation for FAR 25. Default is True.
        Returns
        -------
            np.asarray
                The calculated power loading or thrust-weight ratio for the aircraft based on the provided parameters and distances.

        Reference
        ---------

        """
        if FAR == 23:
            s_tog = distances.get('s_tog',0)
            s_to = distances.get('s_to',0)
            aux_var = s_to < 1.66*s_tog  # check which distance is limiting
            if aux_var and s_to != 0: # can use s_tog
                p = [0.009, 4.9, -s_tog]
            else: # have to use s_to
                p = [0.0149, 8.134, -s_to]
            top23_roots = np.roots(p)
            for r in top23_roots:
                if r.real < 0:
                    top23_roots.remove(r)
            top23 = top23_roots[0].real
            power_load = top23*dens_ratio*cl_max/wing_load
            return power_load
        else: # FAR 25
            s_tofl = distances.get('s_tofl',0)
            thrust_to_weight_ratio = 37.5*wing_load/(
                dens_ratio*cl_max*s_tofl)
            output = (thrust_to_weight_ratio if driven_by_jett
                       else thrust_to_weight_ratio/2.8) # PREGUNTAR
            return output
    @staticmethod # PUEDE CAMBIAR PORQUE PUEDE QUE AÑADA UN SELF (esp. para mantener el far, tipo de motor, etc)
    def landing_distance(
        FAR:int,
        w_ratio:int,
        s_land:int,
        rho:float=0.00237
    )->float:
        """
        def ...

        Parameters
        ----------
        FAR : int
            Federal Aviation Regulation (FAR) part number, either 23 or 25.
        w_ratio : int
            Weight ratio, defined as the landing weight divided by the maximum takeoff weight of the aircraft. By tables ref roskman
        s_land : int
            Landing distance.
        rho : float, optional
            Air density at sea level in slugs per cubic foot. Default is 0.00237 slugs/ft^3.
        
        Returns
        -------
        wing_load_takeoff : float
            Wing loading at takeoff per CL_max, calculated based on the landing distance and weight ratio.
        """
        if FAR == 23:
            v = np.sqrt(s_land/0.265)*1.687810
            wing_load_landing = 0.5*v**2*rho
            wing_load_takeoff = wing_load_landing/w_ratio
            return wing_load_takeoff # per CL_max
        else: # FAR 25
            va = np.sqrt(s_land/0.3)*1.687810
            v_stall = va/1.3
            wing_load_landing = 0.5*v_stall**2*rho
            wing_load_takeoff = wing_load_landing/w_ratio
            return wing_load_takeoff # per CL_max
    @staticmethod
    def areas_f_Wset(
            a:float,
            d:float,
            c:float,
            w_to:float,
            b:float=1.0
    )->dict:
        """
        Calculate the f,S_wet(equivalent parasite, wetted areas) of the aircraft(A/C) by climb performance requirement.

        note:
            The equation is empirical and proposed by Roskman. It takes the data base of different A/C's and
            build the equation. Theofore, the coefficients a,b,c,d change depending on the A/C type.

            The values coeffincients can be found in the Roskman book, Chapter 3.
            - eq:
                log_10(f) = a + b*log_10(S_wet)
                log_10(S_wet) = c + d*log_10(w_to)

            *** Can build a own database and fit the coefficients. Suggest look the graphics of Roskman ***

        Parameters
        ----------
        a,b,c,d : float
            Coefficients of empirical eq. proposed by Roskam.
        w_to : float
            Takeoff weight of the aircraft.
        
        Returns
        -------
        result : dict
            Dict with f and S_wet values.
        """
        s_wet = 10**c*w_to**d
        f = np.exp(a)*np.power(s_wet,b)
        result = dict(f=f, s_wet=s_wet)
        return result
    def drag_polar_estimation(
            wing_load:float,
            w_to:float,
            f:float,
            AR:float,
            e_owswald:dict,
            delta_CD0:dict    
    ):
        """
        Calulcate values of CD0 and k in different configurations.

        Parameters
        ----------
        delta_CD0 : dict
            Dict with values as increment of CD0 in different configurations. 
            
            Structure

                delta_CD0 = {
                    'clean': '0',
                    'TO_flaps': '.01 - .02',
                    'Landing_flaps': '.055 - .075',
                    'landing_gear': '.015 - .025'
                }
                
        e_owswald : dict
            Dict with values as Oswald efficiency factor in different configurations.

            Structure

                e_owswald = {
                    'clean': '0.8 - 0.85',
                    'TO_flaps': '0.75 - 0.8',
                    'Landing_flaps': '0.7 - 0.75',
                }
        wing_load : float
            Wing loading.
        w_to : float
            Takeoff weight.
        f : float
            Equivalent parasite area.
        AR : float
            Aspect ratio of the wing.

        

        Returns
        -------
        result : dict
            Dict with values of CD0 and k in different configurations. Structure
    
                result = {
                    'clean': dict(cd0=cd0_clean, k=k_clean),
                    'take-off':{
                        'gear_up': dict(cd0=cd0_to_flaps, k=k_to_flaps),
                        'gear_down': dict(cd0=cd0_to_flaps + delta_CD0.get('landing_gear',0), k=k_to_flaps)
                        },
                    'landing':{
                        'gear_up': dict(cd0=cd0_landing_flaps, k=k_landing_flaps),
                        'gear_down': dict(cd0=cd0_landing_gear+delta_CD0.get('landing_gear',0), k=k_landing_flaps)
                        }
                }
        
        **Note**
            the values are suggest by Roskman. Can use another values suggest user creteria.  
        """
        S = w_to/wing_load
        cd0_clean = f/S
        cd0_to_flaps = cd0_clean + delta_CD0.get('TO_flaps',0)
        cd0_landing_flaps = cd0_clean + delta_CD0.get('Landing_flaps',0)
        cd0_landing_gear = cd0_clean + delta_CD0.get('landing_gear',0)
        k_clean = 1/(np.pi*e_owswald.get('clean',0)*AR)
        k_to_flaps = 1/(np.pi*e_owswald.get('TO_flaps',0)*AR)
        k_landing_flaps = 1/(np.pi*e_owswald.get('Landing_flaps',0)*AR)
        result = {
            'clean': dict(cd0=cd0_clean, k=k_clean),
            'take-off':{
                'gear_up': dict(cd0=cd0_to_flaps, k=k_to_flaps),
                'gear_down': dict(cd0=cd0_to_flaps + delta_CD0.get('landing_gear',0), k=k_to_flaps)
            },
            'landing':{
                'gear_up': dict(cd0=cd0_landing_flaps, k=k_landing_flaps),
                'gear_down': dict(cd0=cd0_landing_gear+delta_CD0.get('landing_gear',0), k=k_landing_flaps)
            }
        }
        return result
    def rate_climb(
            FAR:int,
            wing_load:np.asarray,
            sigma:float,
            cd0:float,
            k:float,
            configurations:dict,
            rate_climb:float=0.0,
            grandient_climb:float=0.0,
            w_to:float=0.0,
            eta_p:float=None,
            cl_max:float=None,
            e_num:int=1,
            driven_type:str='jet'  
    ):
        """
        Def ...

        note:
            Configuration considered for cumplies with the regulation
            FAR 23
             AOE (All Engines Operating):
              - Rate of climb: gear-up, take off flaps, max. cont. power
              - Gradient of climb: gear-down, landing flaps, take-off. power
             One Engine Inoperative (OEI):
              - Rate of climb: gear-up, flaps most favorable, take-off power in operative engine, max. cont. power in operative engine
            FAR 25:

        Parameters
        ----------
        Far : int
            Federal Aviation Regulation (FAR) part number, either 23 or 25.
        wing_load : np.asarray
            Wing loading of the aircraft, typically in pounds per square foot.
        sigma : float
            Density ratio of the air at the climb conditions compared to sea level standard conditions.
        cd0 : float
            Zero-lift drag coefficient of the aircraft.
        k : float
            Induced drag factor of the aircraft, which is related to the aspect ratio and Oswald
        rate_climb : float, (ft/min)
            Desired rate of climb.
        eta_p : float | default None if FAR 25
            Propulsive efficiency of the aircraft's engines.
        driven_type : str | default 'jet'

        Returns
        -------
        power_load : np.asarray
        """
        
        if (FAR==25 and driven_type.casefold() == 'propeller') or FAR == 23:
            # requierments by regulation
            # AEO rate of climb requirement (FAR 23.65)
             
            if rate_climb > 0:
                AEO_rate_climb = rate_climb
            else: 
                AEO_rate_climb= 300

            cd0 = configurations["take-off"]["gear_up"][cd0]
            k = configurations["take-off"]["gear_up"][k]
            rcp = 33000/AEO_rate_climb
            cl_32 = (3*cd0/k)**(3/4)
            cd = 4*cd0
            ae_ef = cl_32/cd
            w_p = ((rcp+np.fraction(
                np.sqrt(wing_load),
                19*ae_ef*np.sqrt(sigma)
            ))/eta_p)**(-1)
            if rate_climb > 0:
                AEO_power_load = {f'by R/C':w_p}
            else:
                AEO_power_load = {f'.65 by R/C':w_p}
            # AEO climb gradient rate requieriment (FAR 23.65)
            cd0 = configurations["landing"]["gear_down"][cd0]
            k = configurations['landing']["gear_down"][k]
            if grandient_climb > 0:
                gradient_rate = grandient_climb
            else:
                gradient_rate = 1/12 # By regulation
            cl_32 = cl_32 - 0.2
            ae_ef_inv = 1/(cl_32/cd)
            try:
                cgrp = np.fraction(gradient_rate+ae_ef_inv,np.sqrt(cl))
                w_p = np.fraction(18.97*eta_p*np.sqrt(sigma),cgrp*wing_load)
                AEO_power_load.append({'.65 by gradient rate': w_p})
            except ValueError as e:
                cgrp = np.nan
                print(f'Error {e}')

            # AEO climb gradient requirement (FAR 23.77)
            cd0 = configurations["landing"]["gear_down"][cd0]
            k = configurations['landing']["gear_down"][k]
            if gradient_rate > 1/30:
                gradient_rate = grandient_climb
            else:   
                gradient_rate = 1/30 # By regulation
            try:
                cgrp = np.fraction(gradient_rate+ae_ef_inv,np.sqrt(cl))
            except ValueError as e:
                cgrp = np.nan
                print(f'Error {e}')
            w_p = np.fraction(18.97*eta_p*np.sqrt(sigma),cgrp*wing_load)
            AEO_power_load.append({'.77 by gradient rate': w_p})

            # OEI rate climb requirement (FAR 23.67)
            try:
                cd0 = configurations["clean"][cd0]
                k = configurations['clean'][k]
                w_to = w_to.get('w_to',np.nan)
                rho = 0.00237
                cl_max = cl_32 + 0.2
                v_stall = np.sqrt(2*w_to/(rho*cl_max))*60
                rate_climb = 0.027*v_stall**2
                rcp = 33000/rate_climb
                w_p = ((rcp+np.fraction(
                    np.sqrt(wing_load),
                    19*ae_ef*np.sqrt(sigma)
                ))/eta_p)**(-1)
                OEI_power_load = {
                    '.67 by R/C': w_p
                }
            except ValueError as e:
                print(f'Error {e}, returns OEI_power_load={np.nan}')
                OEI_power_load = np.nan

            power_load = dict(
                AEO=AEO_power_load,
                OEI=OEI_power_load
            )
            return power_load
         
        else: # FAR 25
            # OEI .111
            # conf gear up take-off flaps, ground effect, 1.2*V_s_take-off
            # requierments by regulation
            cd0 = configurations["take-off"]["gear_up"][cd0]
            k = configurations["take-off"]["gear_up"][k]
            if gradient_rate > 0.012:
                gradient_rate = grandient_climb
            else:
                gradient_rate = 0.012
            
            cl = np.sqrt(cd0/k)-0.2
            cd = cd0 + cd0
            ae_ef = cl/cd
            thurst_weight_ratio_aeo=np.fraction(1,ae_ef)+grandient_climb
            thurst_weight_ratio_oei = np.fraction(e_num,e_num-1)*thurst_weight_ratio_aeo
            thurst_weight_ratio = dict(
                aeo=thurst_weight_ratio_aeo,
                oei=thurst_weight_ratio_oei
            )
            return thurst_weight_ratio
            




            ##### Importante para jett se traba con empuje
            """
            LANDIING GEAR SIZING

            PARA PATIN LOS LOS ANGULOS DEBEN SER MENOR A ALPHA STALL
            6
            """
