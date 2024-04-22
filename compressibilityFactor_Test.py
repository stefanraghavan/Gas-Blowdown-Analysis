from math import cos
import numpy as np
import numpy
import matplotlib.pyplot as plt
import math
from turtle import end_fill
from typing import Sized
from pandas.core.indexes.base import ensure_index 
from xlrd import open_workbook
from tempfile import TemporaryFile
import numpy as np
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import tkinter as tk
from tkinter import ttk
import pandas as pd
from PIL import Image, ImageTk
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from tkinter import filedialog
from xlsxwriter import Workbook
import tkinter.messagebox as messagebox
import webbrowser
import ttkbootstrap as ttk
from ttkbootstrap.constants import *
import TKinterModernThemes as TKMT
from ttkthemes import ThemedTk
from IPython.display import display, Math
import os
import sys

def resource_path(relative_path):
    """ Get absolute path to resource, works for dev and for PyInstaller """
    base_path = getattr(sys, '_MEIPASS', os.path.dirname(os.path.abspath(__file__)))
    return os.path.join(base_path, relative_path)

def gasPropertyLookup(fluid, property):
    gas_path = resource_path('gasProperties.xlsx')
    gasPropTable = pd.read_excel(gas_path, skiprows=0)
    gas_list     = gasPropTable['Gas'].tolist()
    k_list       = gasPropTable['k'].tolist()
    R_list       = gasPropTable['R'].tolist()
    Cv_list      = gasPropTable['Cv'].tolist()
    Cp_list      = gasPropTable['Cp'].tolist()
    Pc_list      = gasPropTable['Pc'].tolist()
    Tc_list      = gasPropTable['Tc'].tolist()

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
                return (2/(k_list[i]+1))**(k_list[i]/(k_list[i]-1)) # [] - Critical Pressure Ratio for Choked Flow
            if property == "Pc":
                return Pc_list[i] # [Pa] - Critical Pressure
            if property == "Tc":
                return Tc_list[i] # [Pa] - Critical Temperature

def compressibilityFactor(p, T, fluid):

    Pc = gasPropertyLookup(fluid, "Pc")
    Tc = gasPropertyLookup(fluid, "Tc")
    Pr = p/Pc # [] - Reduced Pressure
    Tr = T/Tc # [] - Reduced Temperature
    Ai = 0.42747*(Pr/Tr**(5/2))
    Bi = 0.08664*(Pr/Tr)
    
    q = Bi**2 + Bi -Ai
    r = Ai*Bi

    a3 = 1
    a2 = -1
    a1 = -1*q
    a0 = -1*r

    A = a2/a3
    B = a1/a3
    C = a0/a3

    pi = (-A**2 / 3 + B) / 3
    qi = (9 * A * B - 2 * A**3 - 27 * C) / 54
    Disc = qi**2 + pi**3

    if Disc > 0:
       h = qi + Disc**(1 / 2)
       
       if h < 0:
           y = -(abs(h))**(1 / 3)
       if h >= 0:
           y = (abs(h))**(1 / 3)
       
       return y - pi / y - A / 3
       
    if Disc <= 0:
        theta = math.atan((-Disc) ** (1 / 2) / qi)
        c1 = math.cos(theta / 3) 

        if qi < 0:
           s1 = math.sin(theta / 3)
           c1 = (c1 - s1 * 3 ** (1 / 2)) / 2

        z1 = 2 * (-pi) ** (1 / 2) * c1 - A / 3
        m = A + z1
        r = (m ** 2 - 4 * (B + m * z1)) ** (1 / 2)
        z2 = (-m + r) / 2
        z3 = (-m - r) / 2
        
        if z2 > z1:
            return z2
        if z3 > z1:
            return z3
        if z1 >= z3 or z1 >= z2:
            return z1

# def compressibilityFactor(p, T, fluid, tolerance=1e-6, max_iterations=1000):
#     Pc = gasPropertyLookup(fluid, "Pc")
#     Tc = gasPropertyLookup(fluid, "Tc")
#     Pr = p/Pc
#     Tr = T/Tc
#     A = 0.42747*(Pr/Tr**5/2)
#     B = 0.08664*(Pr/Tr)
#     q = B**2 + B -A
#     r = A*B

#     # Define the cubic function and its derivative
#     def f(z):
#         return z**3 - z**2 - q*z - r

#     def df(z):
#         return 3*z**2 - 2*z - q

#     # Start at a point (arbitrary choice, could be improved)
#     z = 1.0
#     for i in range(max_iterations):
#         # Calculate the value of the function and its derivative at z
#         fz = f(z)
#         dfz = df(z)

#         # Avoid division by zero
#         if dfz == 0:
#             return None

#         # Newton-Raphson formula
#         z_new = z - fz / dfz

#         # Check for convergence
#         if abs(z_new - z) < tolerance:
#             return z_new

#         z = z_new

#     return None  # Return None if no convergence within max_iterations

print(compressibilityFactor(700000, 293, "Helium"))

# press = np.arange(100000, 1000000, 100)  # [Pa]
# Temp = 293

# z = []

# for i in range(len(press)):
#     z.append(compressibilityFactor(press[i], Temp, "Helium"))
# END

# fig1 = plt.figure(1)
# plt.title('Pressure')
# fig1.set_facecolor("w")
# plt.ylabel('Compressibility Factor')
# plt.xlabel('Pressure [Pa]')
# plt.plot(press, z)
# plt.tight_layout()
# plt.legend(loc='best')
# plt.show()