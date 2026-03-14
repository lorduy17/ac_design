import numpy as np
import matplotlib.pyplot as plt

# Speed sizing, with CL suggest

class mc_operations:
    def __init__(self,wing_loading):
        self.wing_loading = wing_loading

    def speed_sizing(self, CL, V_stall, rho):
        """
        Calculate the weight-to-power ratio based on speed sizing.
        """
        w_to = 0.5 * rho * V_stall**2 * CL
        return w_to
    def take_off_distance_far23(self,s_to,sigma,CL_max_TO):
        """
        docstring
        """
        s_relation = 1.66
        def top_23_value():
            c = -s_to
            b = 8.134
            a = 0.0149
            ans = np.roots([a,b,c])
            for index, x in enumerate(ans):
                if np.isreal(x) and x > 0:
                    return x
            ans = 0 if not np.isreal(x) else x
            return ans
        top23 = top_23_value()
        power_loading = top23*sigma*CL_max_TO/self.wing_loading
        return power_loading
    def take_off_distance_far25(self,s_tofl,sigma,CL_max_TO):
        """
        docstring
        """
        ## TOP25 = S_TOFL / 37.5
        top25 = s_tofl / 37.5
        power_landing = self.wing_loading / (top25*sigma*CL_max_TO)
        return power_landing