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

def resource_path(relative_path):
    """ Get absolute path to resource, works for dev and for PyInstaller """
    base_path = getattr(sys, '_MEIPASS', os.path.dirname(os.path.abspath(__file__)))
    return os.path.join(base_path, relative_path)

def export_data():
    global t, p, T, m, mdot, rho  # Declare these variables as global

    # Ask user for file location and name
    file_path = filedialog.asksaveasfilename(defaultextension=".xlsx", 
                                             filetypes=[("Excel files", "*.xlsx"), ("All files", "*.*")])
    if not file_path:
        # User cancelled the dialog, exit the function
        return
    
    unit = selected_unit.get()
    if unit == "SI":
    # Gather data
       data1 = {
           'Time [s]': t[0], 
           'Pressure [Pa]': p[0], 
           'Temperature [K]': T[0], 
           'Mass [kg]': m[0],
           'Flow Rate [kg/s]': mdot[0]
       }
   
       data2 = {
           'Time [s]': t[1], 
           'Pressure [Pa]': p[1], 
           'Temperature [K]': T[1], 
           'Mass [kg]': m[1],
           'Flow Rate [kg/s]': mdot[1]
       }
    if unit == "English":
    # Gather data
       data1 = {
           'Time [s]': t[0], 
           'Pressure [psia]': unitConversion(numpy.array(p[0]), "pressure","english"), 
           'Temperature [degF]': unitConversion(numpy.array(T[0]), "temperature", "english"), 
           'Mass [lbm]': unitConversion(numpy.array(m[0]), "mass", "english"),
           'Flow Rate [kg/s]': unitConversion(numpy.array(mdot[0]),"mass", "english")
       }
   
       data2 = {
           'Time [s]': t[1], 
           'Pressure [psia]': unitConversion(numpy.array(p[1]), "pressure","english"), 
           'Temperature [degF]': unitConversion(numpy.array(T[1]), "temperature", "english"), 
           'Mass [lbm]': unitConversion(numpy.array(m[1]), "mass", "english"),
           'Flow Rate [kg/s]': unitConversion(numpy.array(mdot[1]),"mass", "english")
       }
    # Create a DataFrame
    df1 = pd.DataFrame(data1)
    df2 = pd.DataFrame(data2)
    
   # Save to Excel, each DataFrame in a separate sheet
    with pd.ExcelWriter(file_path, engine='xlsxwriter') as writer:
        df1.to_excel(writer, sheet_name='Isothermal', startrow=0, index=False)
        df2.to_excel(writer, sheet_name='Adiabatic', startrow=0, index=False)

        # Access the xlsxwriter workbook and worksheet objects
        workbook  = writer.book
        worksheet1 = writer.sheets['Isothermal']
        worksheet2 = writer.sheets['Adiabatic']

        # Define a format for removing the border
        format = workbook.add_format({'border': 0})

        # Apply the format to the first row (excluding the header row)
        worksheet1.set_row(1, None, format)
        worksheet2.set_row(1, None, format)

def open_pdf():
    #pdf_path = 'Theory_Assumptions.pdf'  # Replace with the actual path to your PDF file
    pdf_path = resource_path('Theory & Assumptions Real Gas.pdf')
    webbrowser.open(pdf_path)

def unitConversion(value, quantity, targetUnit):
    if quantity == "pressure":
        if targetUnit == "english":
            return value*0.00014504 # convert to psi
        if targetUnit == "si":
            return value*0.00014504**-1 # convert to Pa
    if quantity == "temperature":
        if targetUnit == "english":
           return (value-273.15)*(9/5)+32 # convert to degF
        if targetUnit == "si":
           return (value-32)*(5/9)+273.15 # convert to K
    if quantity == "mass":
        if targetUnit == "english":
            return value*2.20462 # convert to lbm
        if targetUnit == "si":
            return value*2.20462**-1 # convert to kg                    
    if quantity == "length":
        if targetUnit == "english":
            return value*3.28084 # convert to ft   
        if targetUnit == "si":
            return value*3.28084**-1 # convert to m     
    if quantity == "volume":
        if targetUnit == "english":
            return value*35.3147 # convert to ft^3
        if targetUnit == "si": 
            return value*35.3147**-1 # convert to m^3             
    if quantity == "area":
        if targetUnit == "english":
            return value*10.7639 # convert to ft^2
        if targetUnit == "si":
            return value*10.7639**-1 # convert to m^2           
    if quantity == "force":
        if targetUnit == "english":
            return value*0.224809 # convert to lbf
        if targetUnit == "si":
            return value*0.224809**-1 # convert to N     

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

def mdotCalc(flowType, pi, po, rho, CdA, k, Pcrit):
    if flowType == "incompressible":
        return CdA*math.sqrt(2*rho*(pi-po))
    if flowType == "compressible":
       if po/pi <= Pcrit: # flow is choked
           value_under_sqrt = k*rho*pi*(2/(k+1))**((k+1)/(k-1))
        #    if value_under_sqrt < 0:
        #         messagebox.showwarning("Warning", "Invalid CdA, V, and T entry")
        #         return 
           return CdA*math.sqrt(value_under_sqrt)
       if po/pi > Pcrit: # flow is unchoked
        #    return CdA*math.sqrt((2*rho*pi*(k/(k-1))*(po/pi)**(2/k))-(po/pi)**((k+1)/k))
           value_under_sqrt = 2 * rho * pi * (k / (k - 1)) * ((po/pi)**(2 / k) - (po/pi)**((k + 1) / k))
           if value_under_sqrt < 0:
                #print("Negative value under sqrt encountered in mdotCalc")
                return 0
           return CdA * math.sqrt(value_under_sqrt)
    
def update_units_label():
    unit = selected_unit.get()
    if unit == "English":
        # Updating Text widget content for English units
        pi_unit_label.configure(state="normal")
        pi_unit_label.delete(1.0, "end")
        pi_unit_label.insert("end", "psia")
        pi_unit_label.configure(state="disabled")
        pi_unit_label.configure(font=("Helvetica", 11, "bold"))  # Change "12" to your desired font size

        po_unit_label.configure(state="normal")
        po_unit_label.delete(1.0, "end")
        po_unit_label.insert("end", "psia")
        po_unit_label.configure(state="disabled")
        po_unit_label.configure(font=("Helvetica", 11, "bold"))  # Change "12" to your desired font size


        Ti_unit_label.configure(state="normal")
        Ti_unit_label.delete(1.0, "end")
        Ti_unit_label.insert("end", "°F")
        Ti_unit_label.configure(state="disabled")
        Ti_unit_label.configure(font=("Helvetica", 11, "bold"))  # Change "12" to your desired font size

        CdA_unit_label.configure(state="normal")
        CdA_unit_label.delete(1.0, "end")
        CdA_unit_label.insert("end", "in²")
        CdA_unit_label.configure(state="disabled")
        CdA_unit_label.configure(font=("Helvetica", 11, "bold"))  # Change "12" to your desired font size

        V_unit_label.configure(state="normal")
        V_unit_label.delete(1.0, "end")
        V_unit_label.insert("end", "in³")
        V_unit_label.configure(state="disabled")
        V_unit_label.configure(font=("Helvetica", 11, "bold"))  # Change "12" to your desired font size

        axs[0].set(ylabel='Pressure [psia]')
        axs[1].set(ylabel='Temperature [°F]')
        axs[2].set(ylabel='Flow Rate [lbm/s]')

    elif unit == "SI":
        # Updating Text widget content for SI units
        pi_unit_label.configure(state="normal")
        pi_unit_label.delete(1.0, "end")
        pi_unit_label.insert("end", "kPa")
        pi_unit_label.configure(state="disabled")
        pi_unit_label.configure(font=("Helvetica", 11, "bold"))  # Change "12" to your desired font size

        po_unit_label.configure(state="normal")
        po_unit_label.delete(1.0, "end")
        po_unit_label.insert("end", "kPa")
        po_unit_label.configure(state="disabled")
        po_unit_label.configure(font=("Helvetica", 11, "bold"))  # Change "12" to your desired font size

        Ti_unit_label.configure(state="normal")
        Ti_unit_label.delete(1.0, "end")
        Ti_unit_label.insert("end", "K")
        Ti_unit_label.configure(state="disabled")
        Ti_unit_label.configure(font=("Helvetica", 11, "bold"))  # Change "12" to your desired font size


        CdA_unit_label.configure(state="normal")
        CdA_unit_label.delete(1.0, "end")
        CdA_unit_label.insert("end", "cm²")
        CdA_unit_label.configure(state="disabled")
        CdA_unit_label.configure(font=("Helvetica", 11, "bold"))  # Change "12" to your desired font size

        V_unit_label.configure(state="normal")
        V_unit_label.delete(1.0, "end")
        V_unit_label.insert("end", "cm³")
        V_unit_label.configure(state="disabled")
        V_unit_label.configure(font=("Helvetica", 11, "bold"))  # Change "12" to your desired font size

        axs[0].set(ylabel='Pressure [kPa]')
        axs[1].set(ylabel='Temperature [K]')
        axs[2].set(ylabel='Flow Rate [kg/s]')

    canvas.draw()  # Redraw the canvas to update the plot

global canvas, toolbar  # Declare these as global if they are not already

def plot_results():
    global t, p, T, m, mdot, rho  # Declare these variables as global

    update_units_label()
    
    # Retrieve selections from dropdowns
    fluid = selected_fluid.get()
    unit = selected_unit.get()

    #tvent = float(tvent_entry.get()) # [s] - Vent duration

    if unit == "SI": # Keep SI Units
        pi = float(pi_entry.get())*1000 # convert from kPa to Pa
        po = float(po_entry.get())*1000 # convert from kPa to Pa
        Ti = float(Ti_entry.get())
        CdA = float(CdA_entry.get())*0.0001 # convert from cm^2 to m^2
        V = float(V_entry.get())*.000001 # convert from cm^3 to m^3

    if unit == "English": # Convert to SI Units for calcs
        pi = unitConversion(float(pi_entry.get()), "pressure", "si")
        po = unitConversion(float(po_entry.get()), "pressure", "si")
        Ti = unitConversion(float(Ti_entry.get()), "temperature", "si")
        CdA = unitConversion(float(CdA_entry.get())*0.00694444, "area", "si") # convert from in^2 to ft^2
        V = unitConversion(float(V_entry.get())*0.000578704, "volume", "si") # convert from in^3 to ft^3

    # Check if starting pressure is less than reservoir pressure
    if pi <= po:
        messagebox.showwarning("Warning", "Initial pressure must be greater than outlet pressure.")
        return  # Exit the function if the condition is met        
    # Blowdown
    blowdown = ["Isothermal", "Adiabatic"]

    # Calculate gas properties        
    k     = gasPropertyLookup(fluid, "k")
    R     = gasPropertyLookup(fluid, "R")
    Pcrit = gasPropertyLookup(fluid, "Pcrit")
    z     = compressibilityFactor(pi, Ti, fluid)
    
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

       m[j].append(pi*V/(z*R*Ti))
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
           
           if m[j][-1] < 0:
              messagebox.showwarning("Warning", "Invalid CdA: CdA is too large relative to vent volume")
              return 
           
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

    # Plot the results using matplotlib
           
    # Example: Clear existing plot
    for widget in frame_plot.winfo_children():
        widget.destroy()

    # Example: Create a new figure and plot
    color = ['b-','r-']
    fig1, axs = plt.subplots(3, figsize=(5,8))
    fig1.set_tight_layout(True)
    legendloc = "upper right"

    axs[0].clear()
    axs[1].clear()
    axs[2].clear()

    # Clear previous figure
    for widget in frame_plot.winfo_children():
        widget.destroy()


    if unit == "SI": # plot in SI units
       
       for i in range(len(blowdown)):
           axs[0].plot(t[i], numpy.array(p[i])/1000, color[i], label=blowdown[i])
           axs[1].plot(t[i], T[i], color[i], label=blowdown[i])
           axs[2].plot(t[i], mdot[i], color[i], label=blowdown[i])
       
       axs[0].set( ylabel='Pressure [kPa]')
       axs[0].grid(True)
       axs[0].legend(loc=legendloc)
       
       axs[1].set( ylabel='Temperature [K]')
       axs[1].grid(True)
       axs[1].legend(loc=legendloc)
       
       axs[2].set(xlabel='Time [s]', ylabel='Flow Rate [kg/s]')
       axs[2].grid(True)
       axs[2].legend(loc=legendloc)

    if unit == "English": # Convert and plot in English units
       
       for i in range(len(blowdown)):
          axs[0].plot(t[i], unitConversion(numpy.array(p[i]), "pressure", "english"), color[i], label=blowdown[i])
          axs[1].plot(t[i], unitConversion(numpy.array(T[i]), "temperature", "english"), color[i], label=blowdown[i])
          axs[2].plot(t[i], unitConversion(numpy.array(mdot[i]), "mass", "english") , color[i], label=blowdown[i])

       axs[0].set( ylabel='Pressure [psia]')
       axs[0].grid(True)
       axs[0].legend(loc=legendloc)
       
       axs[1].set( ylabel='Temperature [°F]')
       axs[1].grid(True)
       axs[1].legend(loc=legendloc)
       
       axs[2].set(xlabel='Time [s]', ylabel='Flow Rate [lbm/s]')
       axs[2].grid(True)
       axs[2].legend(loc=legendloc)       
    
    # Remove top and right ticks
    for ax in axs:
        # Only remove the tick marks, not the entire spine
        ax.tick_params(top=False, right=False, which='both')  # which='both' refers to both major and minor ticks

    canvas = FigureCanvasTkAgg(fig1, master=frame_plot)  
    canvas.draw()
    canvas.get_tk_widget().pack()
    
    # Update the canvas with the new figure
    canvas.figure = fig1  
    canvas.draw()

    # Recreate the toolbar linked to the updated canvas
    for widget in toolbar_frame.winfo_children():
        widget.destroy()
    toolbar = NavigationToolbar2Tk(canvas, toolbar_frame)
    toolbar.update()
    
# Main window
# root = ttk.Window(themename="darkly")
root = ThemedTk(theme="park")
root.title("Blowdown Analysis")

style = ttk.Style()
#style.configure("Custom.TLabel", background="darkly", foreground="white", borderwidth=0)
style.theme_use('litera')  # Replace 'darkly' with your desired theme name

bg_color = style.lookup('TFrame', 'background')  # Get background color of the theme
style.configure("Borderless.TEntry", bordercolor=bg_color)

# Create the outer frame
outer_frame = tk.Frame(root, borderwidth=2, bg="white")
outer_frame.grid(row=0, column=0, padx=5, pady=5)

# Create an inner frame inside the outer frame
inner_frame = tk.Frame(outer_frame, borderwidth=2, bg="white")
inner_frame.grid(row=0, column=0, padx=5, pady=5)

# Dropdown for fluid
fluidlabel = tk.Label(inner_frame, text="Gas", bg="white",font='Helvetica 11 bold')
fluidlabel.grid(row=0, column=0,  sticky='w')

fluid_options = ["Helium",
                 "Nitrogen",
                 "Oxygen",
                 "Air",
                 "Water Vapor",
                 "Methane",
                 "Hydrogen",
                 "Carbon dioxide",
                 "Carbon monoxide",
                 "Acetone", 
                 "Acetylene",
                 "Ethanol",
                 "Methanol",
                 "Ammonia",
                 "Argon",
                 "Benzene",
                 "Bromine",
                 "Butane",
                ]

selected_fluid = tk.StringVar()
selected_fluid.set(fluid_options[0])

fluidentry = tk.OptionMenu(inner_frame, selected_fluid, *fluid_options)
fluidentry.grid(row=0, column=2, padx=10,pady=10)

# Dropdown for Units
unitlabel = tk.Label(inner_frame, text="Unit System", bg="white",font='Helvetica 11 bold')
unitlabel.grid(row=1, column=0,  sticky='w')

unit_options = ["SI","English"]

selected_unit = tk.StringVar()
selected_unit.set(unit_options[0])
selected_unit.trace("w", lambda *args: update_units_label())

unitentry = tk.OptionMenu(inner_frame, selected_unit, *unit_options)
unitentry.grid(row=1, column=2, padx=10,pady=10)

# Entry widgets for pi, po, Ti, CdA, V, t_vent
pi_label = ttk.Label(inner_frame, text="Initial Pressure", background=bg_color, font='Helvetica 11 bold')
pi_label.grid(row=2, column=0,  sticky='w')
# pi_label.bind("<Enter>", on_hover)


pi_variable_label = ttk.Label(inner_frame, text="Pᵢ", background=bg_color, font='Helvetica 11 bold')
pi_variable_label.grid(row=2, column=1,  sticky='E')

pi_entry = ttk.Entry(inner_frame, width=10)
pi_entry.grid(row=2, column=2, padx=10,pady=10)
pi_entry.insert(0, "700")

po_label = ttk.Label(inner_frame, text="Outlet Pressure", background=bg_color, font='Helvetica 11 bold')
po_label.grid(row=3, column=0,  sticky='w')

po_variable_label = ttk.Label(inner_frame, text="Pₒ", background=bg_color, font='Helvetica 11 bold')
po_variable_label.grid(row=3, column=1,  sticky='E')

po_entry = ttk.Entry(inner_frame, width=10)
po_entry.grid(row=3, column=2, padx=10,pady=10)
po_entry.insert(0, "100")

Ti_label = ttk.Label(inner_frame, text="Initial Temperature", background=bg_color, font='Helvetica 11 bold')
Ti_label.grid(row=4, column=0,  sticky='w')

Ti_variable_label = ttk.Label(inner_frame, text="Tᵢ", background=bg_color, font='Helvetica 11 bold')
Ti_variable_label.grid(row=4, column=1,  sticky='E')

Ti_entry = ttk.Entry(inner_frame, width=10)
Ti_entry.grid(row=4, column=2, padx=10,pady=10)
Ti_entry.insert(0, "293")

CdA_label = ttk.Label(inner_frame, text="Effective Flow Area", background=bg_color, font='Helvetica 11 bold')
CdA_label.grid(row=5, column=0,  sticky='w')

CdA_variable_label = ttk.Label(inner_frame, text="CdA", background=bg_color, font='Helvetica 11 bold')
CdA_variable_label.grid(row=5, column=1,  sticky='E')

CdA_entry = ttk.Entry(inner_frame, width=10)
CdA_entry.grid(row=5, column=2, padx=10,pady=10)
CdA_entry.insert(0, "0.07")

V_label = ttk.Label(inner_frame, text="Vent Volume", background=bg_color, font='Helvetica 11 bold')
V_label.grid(row=6, column=0,  sticky='w')

V_variable_label = ttk.Label(inner_frame, text="V", background=bg_color, font='Helvetica 11 bold')
V_variable_label.grid(row=6, column=1,  sticky='E')

V_entry = ttk.Entry(inner_frame, width=10)
V_entry.grid(row=6, column=2, padx=10,pady=10)
V_entry.insert(0, "600")

# Unit labels

pi_unit_label = tk.Text(inner_frame, height=1, width=5, bg=bg_color, bd=0, highlightthickness=0, font='Helvetica 11 bold')
pi_unit_label.grid(row=2, column=3, sticky='w')
pi_unit_label.insert(tk.END, "kPa")
pi_unit_label.tag_configure("bold", font=("Helvetica", 11, "bold"))
pi_unit_label.configure(state="disabled")


# po_unit_label = tk.Label(inner_frame, text="Pa", bg="white",font='Helvetica 8 bold')
# po_unit_label.grid(row=3, column=2)

po_unit_label = tk.Text(inner_frame, height=1, width=5, bg="white", bd=0, font='Helvetica 11 bold')
po_unit_label.grid(row=3, column=3, sticky='w')
po_unit_label.insert(tk.END, "kPa")
po_unit_label.tag_configure("bold", font=("Helvetica", 11, "bold"))
po_unit_label.configure(state="disabled")

# Ti_unit_label = tk.Label(inner_frame, text="K", bg="white",font='Helvetica 8 bold')
# Ti_unit_label.grid(row=4, column=2)

Ti_unit_label = tk.Text(inner_frame, height=1, width=5, bg="white", bd=0, font='Helvetica 11 bold')
Ti_unit_label.grid(row=4, column=3, sticky='w')
Ti_unit_label.insert(tk.END, "K")
Ti_unit_label.tag_configure("bold", font=("Helvetica", 11, "bold"))
Ti_unit_label.configure(state="disabled")

# CdA_unit_label = tk.Label(inner_frame, text="m^2", bg="white",font='Helvetica 8 bold')
# CdA_unit_label.grid(row=5, column=2)

CdA_unit_label = tk.Text(inner_frame, height=1, width=5, bg="white", bd=0, font='Helvetica 11 bold')
CdA_unit_label.grid(row=5, column=3,  sticky='w')
CdA_unit_label.insert(tk.END, "cm", "bold")
CdA_unit_label.insert(tk.END, "2", "bold_superscript")  # Marking "2" as superscript
CdA_unit_label.tag_configure("bold", font=("Helvetica", 11, "bold"))
CdA_unit_label.tag_configure("bold_superscript", offset=4, font=("Helvetica", 8))  # Adjust font size and offset as needed
CdA_unit_label.configure(state="disabled")

# V_unit_label = tk.Label(inner_frame, text="m^3", bg="white",font='Helvetica 8 bold')
# V_unit_label.grid(row=6, column=2)

V_unit_label = tk.Text(inner_frame, height=1, width=5, bg="white", bd=0, font='Helvetica 11 bold')
V_unit_label.grid(row=6, column=3,  sticky='w')
V_unit_label.insert(tk.END, "cm", "bold")
V_unit_label.insert(tk.END, "3", "bold_superscript")  # Marking "2" as superscript
V_unit_label.tag_configure("bold", font=("Helvetica", 11, "bold"))
V_unit_label.tag_configure("bold_superscript", offset=4, font=("Helvetica", 8))  # Adjust font size and offset as needed
V_unit_label.configure(state="disabled")

pi_unit_label.config(bg=bg_color, bd=0, highlightthickness=0)
po_unit_label.config(bg=bg_color, bd=0, highlightthickness=0)
Ti_unit_label.config(bg=bg_color, bd=0, highlightthickness=0)
CdA_unit_label.config(bg=bg_color, bd=0, highlightthickness=0)
V_unit_label.config(bg=bg_color, bd=0, highlightthickness=0)

# Button to execute and plot
button = tk.Button(inner_frame, text="Plot", command=plot_results,font='Helvetica 11 bold', height= 2, width=8)
button.grid(row=8, column=2, padx=10,pady=10)

# Button to export data to excel
button2 = tk.Button(inner_frame, text="Export", bg='light green', command=plot_results,font='Helvetica 11 bold',height= 2, width=8)
button2.config(command=export_data)
button2.grid(row=9, column=2, padx=10,pady=10)

# Theory and Assumptions Button
pdf_button = tk.Button(inner_frame, text="Documentation", command=open_pdf, font='Helvetica 11 bold', height=5, width=18,bg="#A9A9A9")
pdf_button.grid(row=8, column=0, rowspan=2, padx=10, pady=10, sticky='w')
pdf_button.grid(row=8, column=0, rowspan=2,columnspan=2, padx=(15, 10), pady=10,sticky='w')

# Display Empty plots
fig1, axs = plt.subplots(3, figsize=(5,8))
fig1.set_tight_layout(True)

axs[0].plot([], [])
axs[1].plot([], [])
axs[2].plot([], [])
       
axs[0].set( ylabel='Pressure [kPa]')
axs[0].grid(True)
       
axs[1].set( ylabel='Temperature [K]')
axs[1].grid(True)
       
axs[2].set(xlabel='Time [s]', ylabel='Flow Rate [kg/s]')
axs[2].grid(True)

# Creating the Matplotlib canvas
canvas = FigureCanvasTkAgg(fig1, master=outer_frame)  
canvas_widget = canvas.get_tk_widget()
canvas_widget.grid(row=0, column=4, padx=10, pady=10)

# Create a separate frame for the toolbar
toolbar_frame = tk.Frame(outer_frame)
toolbar_frame.grid(row=1, column=4, padx=10, pady=10)

# Add the Matplotlib navigation toolbar
toolbar = NavigationToolbar2Tk(canvas, toolbar_frame)
toolbar.update()

# Frame for plot
frame_plot = tk.Frame(outer_frame)
frame_plot.grid(row=0, column=4, padx=10,pady=10)

# Load the image with PIL
diagram_path = resource_path('diagram.jpg')
original_image = Image.open(diagram_path)  # Replace with your image path

# Desired new width (change this as needed)
new_width = 450

# Calculate the new height to maintain aspect ratio
original_width, original_height = original_image.size
aspect_ratio = original_height / original_width
new_height = int(new_width * aspect_ratio)

# Resize the image
resized_image = original_image.resize((new_width, new_height))

# Convert to PhotoImage
image = ImageTk.PhotoImage(resized_image)

# Assuming you have already created 'inner_frame' in your Tkinter application
image_label = tk.Label(inner_frame, image=image)
image_label.grid(row=10, column=0, columnspan=5, padx=10, pady=10)

# Important: Keep a reference to the image to prevent garbage-collection
image_label.image = image

# Load the image
logo_path = resource_path('logo2.jpg')
icon_image = Image.open(logo_path)  # Replace with the path to your image file

# Convert the image to PhotoImage
icon_photo = ImageTk.PhotoImage(icon_image)

# Set the window icon
root.iconphoto(False, icon_photo)

plot_results()

root.mainloop()
