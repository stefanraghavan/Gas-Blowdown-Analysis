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
import xlwings as xw
from xlrd import open_workbook
import pandas
import xlwt
from tempfile import TemporaryFile
import os
import subprocess
import numpy as np
import pandas as pd
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import *


root = tk.Tk()

root.title("Simple Calculator")

fluidlabel = tk.Label(root)
e = tk.Entry(root, width=35, borderwidth=5)
e.grid(row=0, column=0, columnspan=3, padx=10,pady=10)

def button_click(number):
    #e.delete(0,END)
    current = e.get()
    e.delete(0, END)
    e.insert(0,str(current) + str(number))

def button_click_clear():
    e.delete(0, END)

def button_add():
    first_number =e.get()
    global f_num
    global math
    math = "addition"
    f_num = int(first_number)
    e.delete(0, END)

def button_equal():
    second_number = e.get()
    e.delete(0, END)

    if math == "addition":
        e.insert(0, f_num + int(second_number))
    if math == "subtraction":
        e.insert(0, f_num - int(second_number))
    if math == "multiplication":
        e.insert(0, f_num * int(second_number))
    if math == "division":
        e.insert(0, f_num / int(second_number))

def button_subtract():
    first_number =e.get()
    global f_num
    global math
    math = "subtraction"
    f_num = int(first_number)
    e.delete(0, END)

def button_multiply():
    first_number =e.get()
    global f_num
    global math
    math = "multiplication"
    f_num = int(first_number)
    e.delete(0, END)

def button_divide():
    first_number =e.get()
    global f_num
    global math
    math = "division"
    f_num = int(first_number)
    e.delete(0, END)

# Define buttons
button1 = tk.Button(root, text="1", padx=40, pady=20, command=lambda: button_click(1))
button2 = tk.Button(root, text="2", padx=40, pady=20, command=lambda: button_click(2))
button3 = tk.Button(root, text="3", padx=40, pady=20, command=lambda: button_click(3))
button4 = tk.Button(root, text="4", padx=40, pady=20, command=lambda: button_click(4))
button5 = tk.Button(root, text="5", padx=40, pady=20, command=lambda: button_click(5))
button6 = tk.Button(root, text="6", padx=40, pady=20, command=lambda: button_click(6))
button7 = tk.Button(root, text="7", padx=40, pady=20, command=lambda: button_click(7))
button8 = tk.Button(root, text="8", padx=40, pady=20, command=lambda: button_click(8))
button9 = tk.Button(root, text="9", padx=40, pady=20, command=lambda: button_click(9))
button0 = tk.Button(root, text="0", padx=40, pady=20, command=lambda: button_click(0))
button_add = tk.Button(root, text="+", padx=39, pady=20, command=button_add)
button_equal = tk.Button(root, text="=", padx=91, pady=20, command=button_equal)
button_clear = tk.Button(root, text="Clear", padx=79, pady=20, command=lambda: button_click_clear())

button_subtract = tk.Button(root, text="-", padx=41, pady=20, command=button_subtract)
button_multiply = tk.Button(root, text="*", padx=40, pady=20, command=button_multiply)
button_divide = tk.Button(root, text="/", padx=41, pady=20, command=button_divide)

# Put the buttons on screeen
button1.grid(row=3, column=0)
button2.grid(row=3, column=1)
button3.grid(row=3, column=2)

button4.grid(row=2, column=0)
button5.grid(row=2, column=1)
button6.grid(row=2, column=2)

button7.grid(row=1, column=0)
button8.grid(row=1, column=1)
button9.grid(row=1, column=2)

button0.grid(row=4, column=0)

button_add.grid(row=5, column=0)
button_equal.grid(row=5, column=1, columnspan=2)
button_clear.grid(row=4, column=1, columnspan=2)

button_subtract.grid(row=6, column=0)
button_multiply.grid(row=6, column=1)
button_divide.grid(row=6, column=2)




root.mainloop() 