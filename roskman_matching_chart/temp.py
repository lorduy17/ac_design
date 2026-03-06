import numpy as np
import matplotlib.pyplot as plt

## FAR 35
def general_eq():
    None

#
top25 = 133.3
CLmax=1.2
def T_W(W_S,CL,sigma,top25):
    return W_S/(top25*CL*sigma)

plt.figure()
W_S = np.linspace(40,100, 100)
for i in [1.2,1.6,2,2.4]:
    plt.plot(W_S, T_W(W_S,i,1,top25))
plt.xlabel('W/S')
plt.ylabel('T/W')
plt.title('T/W vs W/S for different CL')
plt.legend(['CL=1.2','CL=1.6','CL=2','CL=2.4'])
plt.grid()
plt.show(block=False)

# Landing

def _(rho,V,CL):
    return (0.5*rho*V**2*CL)/0.95

ax = plt.figure().add_subplot(111)
W_S = np.linspace(40,100, 100)
for i in [1.2,1.6,2,2.4]:
    ax.axvline((0.00204834,69.8*1.688,i),0,1)
plt.xlabel('W/S')
plt.ylabel('T/W')
plt.title('T/W vs W/S for different CL')
plt.legend(['CL=1.2','CL=1.6','CL=2','CL=2.4'])
plt.grid()
plt.show()




