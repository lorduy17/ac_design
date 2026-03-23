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
            rate_climb:float=0.0,
            grandient_climb:float=0.0,
            v:float=0.0,
            eta_p:float=None,
            cl_max:float=None,
            e_num:int=1,
            driven_type:str='jet'  
    ):
        """
        Def ...

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
        if FAR == 23:
            # power loading by rate of climb
            if rate_climb == 0:
                power_load_by_rate_climb = np.nan
            else:
                rcp = rate_climb/33000 # rate of climb in ft/min
                cl_32 = np.sqrt(3*cd0/k)
                cd = 4*cd0
                ae_ef = cl_32/cd
                power_load_by_rate_climb = ((rcp+np.fraction(
                    np.sqrt(wing_load),
                    19*ae_ef*np.sqrt(sigma)
                ))/eta_p)**(-1)
            # power loading by climb gradient rate
            """
            Preguntar porque para el caso de gradiente de aseceso, el cl que se usa
            es 0.2 menos que el cl max. (muy posiblemente la respuesta sea si
            porque en el ejemplo de Roskman se hace eso, pero no se especifica en el libro)
            """
            if grandient_climb == 0:
                power_load_by_gradient_climb = np.nan
            else:
                cl = cl_32 - 0.2
                cd = cd0 + k*cl**2
                ae_ef = cl/cd
                cgr_plus_ae = grandient_climb + np.fraction(1,ae_ef)
                power_load_by_gradient_climb = np.fraction(
                    np.sqrt(sigma*cl)*18.97*eta_p,
                    cgr_plus_ae*np.sqrt(wing_load)
                )
            power_loading=dict(
                rate_climb=power_load_by_rate_climb,
                gradient_climb=power_load_by_gradient_climb
            )
            return power_loading
        else: # FAR 25
            if driven_type.casefold() == 'jet':
                cd = cd0 + k*cl_max**2
                ae_ef = cl_max/cd
                thurst_weight_ratio_aeo=np.fraction(1,ae_ef)+grandient_climb
                thurst_weight_ratio_oei = np.fraction(e_num,e_num-1)*thurst_weight_ratio_aeo
                thurst_weight_ratio = dict(
                    aeo=thurst_weight_ratio_aeo,
                    oei=thurst_weight_ratio_oei
                )
                return thurst_weight_ratio
            
