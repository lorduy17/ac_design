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
    def  climb(
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
        # cd0 = f/S
        s_wet = 10**c*w_to**d
        f = np.exp(a)*np.power(s_wet,b)
        result = dict(f=f, s_wet=s_wet)
        return result
    def drag_polar_estimation(
            wing_load: float,
            w_to: float,
            f: float,
            S_wet: float,
            AR: float,
            e: dict,
            delta_CD0: dict
    ) -> dict:
        """
        Estimate the drag polar coefficients (CD0 and k) for multiple configurations.

        This method uses a simple drag polar model:
          CD = CD0 + k * CL^2
        where:
          CD0 (parasite drag at CL=0) is computed from f/S (clean condition) plus
          configuration-specific increments from delta_CD0.
          k = 1/(pi * AR * e).

        Parameters
        ----------
        wing_load : float
            Wing loading [lb/ft^2].
        w_to : float
            Takeoff weight [lb].
        f : float
            Equivalent parasite drag coefficient times area (unitless).
        S_wet : float
            Wetted area [ft^2] (included for completeness, not used in current model). 
        AR : float
            Wing aspect ratio.
        e : dict
            Oswald efficiency factors for configurations.
            Required key: 'clean'. Optional: 'TO_flaps', 'Landing_flaps', 'landing_gear'.
        delta_CD0 : dict
            Configuration drag-rise increments with same keys as `e`.

            Example:
            >>> delta_CD0 = {
            >>>     'clean': 0.0,
            >>>     'TO_flaps': 0.015,
            >>>     'Landing_flaps': 0.04,
            >>>     'landing_gear': 0.01
            >>> }

        Returns
        -------
        dict
            Dictionary containing:
              - S: wing area [ft^2]
              - cd0_clean: clean zero-lift drag coefficient
              - polar: nested config dictionary with 'CD0' and 'k'

            Example:
            {
                'S': ..., 'cd0_clean': ..., 
                'polar': {
                    'clean': {'CD0': ..., 'k': ...},
                    'TO_flaps': {...},
                    ...
                }
            }

        Notes
        -----
        * S_wet is accepted for compatibility with Roskam-style formulas,
          but this implementation uses `f/S` for CD0 (clean) and ignores S_wet directly.
        """
        # core computation
        S = w_to / wing_load
        cd0_clean = f / S

        polar = {}
        config_list = ['clean', 'TO_flaps', 'Landing_flaps', 'landing_gear']

        for cfg in config_list:
            delta = float(delta_CD0.get(cfg, 0.0))
            e_cfg = float(e.get(cfg, e.get('clean', 1.0)))
            k_cfg = 1.0 / (np.pi * AR * e_cfg)

            polar[cfg] = {
                'CD0': cd0_clean + delta,
                'k': k_cfg,
            }

        return {
            'S': S,
            'cd0_clean': cd0_clean,
            'polar': polar,
        }
