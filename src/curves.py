import numpy as np
import matplotlib.pyplot as plt
import ops as f

def power_required_plot(ax,V,Pr,Pa,ac_pars,AEmax):
    W = ac_pars["W"]["value"]
    vr,ve = f.operation_speeds(ac_pars,AEmax)
    ax.plot(V,Pr/550,color="mediumorchid",label=r"$P_{required}$")
    ax.plot(V,Pa/550,color="red",label=r"$P_{available}$")
    if ac_pars["type"].casefold() == "propeller":
        ax.axline((0,0),(vr,W/AEmax["1"]*vr/550),linestyle = "--",color="dimgray")
    ax.axvline(x=vr, linestyle ="-.",color="black")
    ax.axvline(x=ve, linestyle ="--",color="black")
    ax.set_ylabel(r"P, $HP$")
    ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
    ax.grid(True)
def thrust_required_plot(ax,V,Tr,ac_pars,AEmax):
    W = ac_pars["W"]["value"]
    vr,ve = f.operation_speeds(ac_pars,AEmax)
    ax.plot(V,Tr,color="mediumorchid",label=r"$T_{required}$")
    if ac_pars["type"].casefold() == "jett":
        ax.axline((0,0),(ve,W/AEmax["1/2"]),linestyle = "--",color="dimgray")
    ax.axvline(x=vr, linestyle ="-.",color="black",label=rf"$V_{{range}}$ = {vr:.2f} ft/s")
    ax.axvline(x=ve, linestyle ="--",color="black",label=rf"$V_{{endurance}}$ = {ve:.2f} ft/s")
    ax.set_ylabel("Tr, lbf")
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    ax.grid(True)
def aerodynamic_ef_plot(ax,V,AE,AEmax,ac_pars):
    ax.plot(V,AE["1"],label=r"$\frac{C_L}{C_D}$")
    ax.plot(V,AE["1/2"],label = r"$\frac{\sqrt{C_L}}{C_D}$")
    ax.plot(V,AE["3/2"],label=r"$\frac{\sqrt{C_L^3}}{C_D}$")
    vr,ve = f.operation_speeds(ac_pars,AEmax)
    ax.axvline(x=vr, linestyle ="-.",color="black")
    ax.axvline(x=ve, linestyle ="--")
    ax.set_xlabel(r"$V_{\infty}, \frac{ft}{s}$")
    ax.set_ylabel(r"$\frac{C_L}{C_D}$")
    ax.grid()
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
def rate_of_climb_plot(ax,RoC,RoC_max_value,V):
    ax.plot(V,RoC,color="mediumorchid",label=r"$R_C$")
    ax.axhline(y = RoC_max_value,linestyle="--",color="black",label=rf"$\frac{{R}}{{C}} = {RoC_max_value:.2f} \frac{{ft}}{{s}}$")
    ax.set_xlabel(r"$V_{\infty}, ft/s$")
    ax.set_ylabel(r"$\frac{R}{C}, ft/s$")
    ax.legend(bbox_to_anchor=(1.05, 1),loc='upper left')
    ax.grid(True)
def roc_height_plot(ax,roc,h):
    roc = roc*1e-3
    h = h*1e-3
    ax.plot(roc,h)
    ax.set_xlabel(r"$\frac{R}{C}_{max}, ft/min \cdot 10^3$")
    ax.set_ylabel(r"$Height, ft \cdot 10^3$")
    ax.legend()
    ax.grid(True)
def climb_time_plot(ax,roc,h,time):
    ax.plot(h,roc,color="darkblue",label=f"Time to climb = {time:.2f}")
    ax.fill_between(h,roc,0,alpha=0.3)
    ax.set_xlabel(r"$Altitude, ft \cdot 10^3$")
    ax.set_ylabel(r"$\frac{R}{C}^{-1}, ft/min^{-1} \cdot 10^{-3}$")
    ax.legend(bbox_to_anchor=(1.05, 1),loc='upper left')
    ax.grid(True)
    ax.legend()
def plot_ac_performance(CJ1,points):
    V = np.linspace(CJ1["v_range"]["value"][0],CJ1["v_range"]["value"][1],points)
    q_inf = f.q(V,CJ1["rho"]["value"])
    AE,AE_max = f.aerodynamic_coeficients(CJ1,q_inf)
    Tr = f.thrust_required(CJ1,AE)
    Pr,Pa = f.power(CJ1,Tr,V)
    RoC,RoC_max_val,RoC_max_values,h,celling = f.rate_of_climb(CJ1,Pa,Pr,V,points)
    # Figure 1
    fig, ax = plt.subplots(3,1)
    thrust_required_plot(ax[0],V,Tr,CJ1,AE_max)
    power_required_plot(ax[1],V,Pr,Pa,CJ1,AE_max)
    aerodynamic_ef_plot(ax[2],V,AE,AE_max,CJ1)
    for axs in ax:
        axs.set_xlim(CJ1["v_range"]["value"][0]-10,CJ1["v_range"]["value"][1]+10)
    plt.tight_layout()
    plt.show(block=False)
    # Figure 2
    fig, ax = plt.subplots(1)
    rate_of_climb_plot(ax,RoC,RoC_max_val,V)
    plt.tight_layout()
    plt.show(block=False)
    # Figure 3
    fig, ax = plt.subplots(1)
    roc_height_plot(ax,RoC_max_values,h)
    plt.tight_layout()
    plt.show(block = False)
    # Figure 4
    fig, ax = plt.subplots(1)
    roc_inv,h,time = f.climb_time(RoC_max_values,h,celling)
    climb_time_plot(ax,roc_inv,h,time)
    # Endurance and range
    E,R = f.endurance_range(CJ1,AE_max)
    plt.show()
CJ1 = {"type":"jett",
    "v_range":{"value":[100,1000],"unit":"ft/s"},
    "rho": {"value":0.0023769,"unit":"slugs/ft^3"},
    "S": {"value":318,"unit":"ft^2"},  # 
    "W": {"value":19815,"unit":"lbf"},
    "Wf": {"value":19815-1119*6.67,'unit':'lbf'},
    "CD0": 0.02,
    "b": {"value":53.3,"unit":"ft"},
    "e": 0.81,
    "thrust_max": {"value":3650*2,"unit":"lbs"},
    "eta_p":0.6,
    "TSFC":0.6
}
CP1 = {"type":"propeller",
    "v_range":{"value":[100,350],"unit":"ft/s"},
    "rho": {"value":0.0023769,"unit":"slugs/ft^3"},
    "S": {"value":174,"unit":"ft^2"},  
    "W": {"value":2950,"unit":"lbs"},
    "Wf": {"value":2950-5.64*65,'unit':'lbf'},
    "CD0": 0.025,
    "b": {"value":35.8,"unit":"ft"},
    "e": 0.8,
    "power_available": {"value":230*550,"unit":"HP"},
    "eta_p": 0.8,
    "SFC":0.45
    }
plot_ac_performance(CJ1,1000)
plot_ac_performance(CP1,1000)
