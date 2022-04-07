from transmission_coefficient import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import MultipleLocator



def g(E,Ee,Eh,V,d):
  return 1-abs(transmissionParams(E,Ee,Eh,V,d,'ref', 1))**2+abs(transmissionParams(E,Ee,Eh,V,d,'ref', 2))**2
d = 0.5e-10*Meter
V_arr = np.linspace(5*eV, 0, 100)
g_arr = [g(0,0,0,V,d) for V in V_arr]
plt.plot(V_arr/eV, g_arr)
plt.show()