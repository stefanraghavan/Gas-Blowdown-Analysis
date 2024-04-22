from math import cos
import numpy as np
import numpy
import matplotlib.pyplot as plt
import math
import math
from turtle import end_fill
from typing import Sized
import matplotlib.pyplot as plt
import numpy as np
import numpy
from pandas.core.indexes.base import ensure_index 
import xlrd 
from xlrd import open_workbook
import pandas
from tempfile import TemporaryFile
import os
import subprocess
import numpy as np
import pandas as pd
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

# Conversion Factors
psi2pa = 0.00014504**-1
lbm2kg = 2.20462**-1
cuft2cum = 35.3147**-1
ftsq2msq = 10.7639**-1

def unitConversion(value, quantity, targetUnit):
    if quantity == "pressure":
        if targetUnit == "english":
            return value*0.00014504 # convert to psi
        else:
            return value*0.00014504**-1 # convert to Pa
    if quantity == "temperature":
        if targetUnit == "english":
           return (value-273.15)*(9/5)+32 # convert to degF
        else:
           return (value-32)*(5/9)+273.15 # convert to K
    if quantity == "mass":
        if targetUnit == "english":
            return value*2.20462 # convert to lbm
        else:
            return value*2.20462**-1 # convert to kg                    
    if quantity == "length":
        if targetUnit == "english":
            return value*3.28084 # convert to ft   
        else:
            return value*3.28084**-1 # convert to m     
    if quantity == "volume":
        if targetUnit == "english":
            return value*35.3147 # convert to ft^3
        else: 
            return value*35.3147**-1 # convert to m^3             
    if quantity == "area":
        if targetUnit == "english":
            return value*10.7639 # convert to ft^2
        else:
            return value*10.7639**-1 # convert to m^2           
    if quantity == "force":
        if targetUnit == "english":
            return value*0.224809 # convert to lbf
        else:
            return value*0.224809**-1 # convert to N               

def gasPropertyLookup(fluid, property):
    gasPropTable = pd.read_excel("gasProperties.xlsx", skiprows=0)
    gas_list     = gasPropTable['Gas'].tolist()
    k_list       = gasPropTable['k'].tolist()
    R_list       = gasPropTable['R'].tolist()
    Cv_list      = gasPropTable['Cv'].tolist()
    Cp_list      = gasPropTable['Cp'].tolist()
    
    for i in range(len(gas_list)):
        if gas_list[i] == fluid:
            if property == "k":
               return k_list[i] # [] - k (Cp/Cv)
            if property == "R":
                return R_list[i]*1000 # [J/(kg*K)] - Gas Constant
            if property == "Cv":
                return Cv_list[i]*1000 # [J/(kg*K)] - Cv
            if property == "Cp":
                return Cp_list[i]*1000 # [J/(kg*K)] - Cp
            if property == "Pcrit":
                return (2/(k_list[i]+1))**(k_list[i]/(k_list[i]-1)) # [] - Pcrit

def mdotCalc(flowType, pi, po, rho, CdA, k, Pcrit):
    if flowType == "incompressible":
        return CdA*math.sqrt(2*rho*(pi-po))
    if flowType == "compressible":
       if po/pi <= Pcrit: # flow is choked
           return CdA*math.sqrt(k*rho*pi*(2/(k+1))**((k+1)/(k-1)))
       if po/pi > Pcrit: # flow is unchoked
        #    return CdA*math.sqrt((2*rho*pi*(k/(k-1))*(po/pi)**(2/k))-(po/pi)**((k+1)/k))
           value_under_sqrt = 2 * rho * pi * (k / (k - 1)) * ((po/pi)**(2 / k) - (po/pi)**((k + 1) / k))
           if value_under_sqrt < 0:
                print("Negative value under sqrt encountered in mdotCalc")
                return 0
           return CdA * math.sqrt(value_under_sqrt)

# Fluid Inputs
fluid = "Helium"
pi = 344738 # [pa] - Starting pressure
po = 103421 # [Pa] - Final pressure
Ti  = 293 # [K] - Starting temperature

# Hardware Inputs
V = 0.021*cuft2cum # [m^3] - Vent volume
CdA = 7.93e-5*ftsq2msq # [m^2] - Effective flow area

 # Blowdown
blowdown = ["Isothermal","Adiabatic"]
# Calculate gas properties        
k     = gasPropertyLookup(fluid, "k")
R     = gasPropertyLookup(fluid, "R")
Pcrit = gasPropertyLookup(fluid, "Pcrit")

# Define blowdown duration
#t = np.arange(0, tvent, 0.01)  # [s]

# Initialize Lists
m    = []    # [kg]     - Mass
mdot = []    # [kg/s]   - Mass flow rate
rho  = []    # [kg/m^3] - Density
p    = []    # [Pa]     - Pressure
T    = []    # [K]      - Temperature
t    = []    # [s]      - Time

for j in range(len(blowdown)):
   
   m.append([])
   mdot.append([])
   rho.append([])
   p.append([])
   T.append([])
   t.append([])

   m[j].append(pi*V/(R*Ti))
   rho[j].append(m[j][0]/V)
   mdot[j].append(mdotCalc("compressible", pi, po, rho[j][0], CdA, k, Pcrit))
   p[j].append(pi)
   T[j].append(Ti)
   t[j].append(0)

   # Initialize time variable
   time_counter = 0
   t_step = 0.001  # [s] - You can adjust the time step as needed
   while p[j][-1] > po+1:  # Continue loop until pressure p is less or equal to po
       time_counter += t_step
       t[j].append(time_counter)
       m[j].append(m[j][-1] - mdot[j][-1] * t_step)
       rho[j].append(m[j][-1] / V)

       if blowdown[j] == "Adiabatic":
           p[j].append(p[j][-1] * (rho[j][-1] / rho[j][-2])**k)
           T[j].append(T[j][-1] * (p[j][-1] / p[j][-2])**((k - 1) / k))

       elif blowdown[j] == "Isothermal":
           p[j].append(rho[j][-1] * R * Ti)
           T[j].append(Ti)

       mdot[j].append(mdotCalc("compressible", p[j][-1], po, rho[j][-1], CdA, k, Pcrit))
       if mdot[j][-1] == 0:
           break

color = ['b-','r-']
fig1, axs = plt.subplots(3, figsize=(7,10))

for i in range(len(blowdown)):
    axs[0].plot(t[i], numpy.array(p[i])/1000, color[i], label=blowdown[i])
    axs[1].plot(t[i], T[i], color[i])
    axs[2].plot(t[i], mdot[i], color[i])

axs[0].set( ylabel='Pressure [kPa]')
axs[0].grid(True)
#axs[0].legend(loc='best')

axs[1].set( ylabel='Temperature [K]')
axs[1].grid(True)
#axs[1].legend(loc='best')

axs[2].set(xlabel='Time [s]', ylabel='Mass Flow Rate [kg/s]')
axs[2].grid(True)
#axs[2].legend(loc='best')

lines = [] 
labels = [] 
  
for ax in fig1.axes: 
    Line, Label = ax.get_legend_handles_labels() 
    # print(Label) 
    lines.extend(Line) 
    labels.extend(Label)

fig1.legend(loc='upper center') 

lines = []
labels = []
plt.show()


# fig1 = plt.figure(1)
# plt.title('Pressure')
# fig1.set_facecolor("w")
# plt.ylabel('Pressure [Pa]')
# plt.xlabel('Time [s]')
# for i in range(len(blowdown)):
#     plt.plot(t, p[i], color[i], label=blowdown[i])
# plt.tight_layout()
# plt.legend(loc='best')
# plt.show()

# fig2 = plt.figure(2)
# plt.title('Temperature')
# fig2.set_facecolor("w")
# plt.ylabel('Temperature [K]')
# plt.xlabel('Time [s]')
# for i in range(len(blowdown)):
#     plt.plot(t, T[i], color[i], label=blowdown[i])
# plt.tight_layout()
# plt.legend(loc='best')
# plt.show()