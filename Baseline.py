from tkinter import *
from tkinter.ttk import *
import tkinter as tk
from tkinter import ttk, Menu, filedialog
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal, sparse
from scipy.signal import savgol_filter
from scipy.sparse.linalg import spsolve
import os


##import seaborn as sns

##sns.set()

##def update_oS(val):
##    global O
##    O = float(val)
####    txt_lam_var.set(round(O, 5))
##    Rplot()
##
##def update_tS(val):
##    global T
##    T = float(val) * 0.001
####    txt_lam1_var.set(round(T, 5))
##    Rplot()

def update_lamE(event):
    global LAM
    LAM = float(txt_lam.get())
    Rplot()


def update_lam1E(event):
    global LAM1
    LAM1 = float(txt_lam1.get())
    Rplot()


def update_pE(event):
    global P
    P = float(txt_p.get())
    Rplot()


def update_lamS(val):
    global LAM
    LAM = float(val)
    txt_lam_var.set(round(LAM, 4))
    Rplot()


def update_lam1S(val):
    global LAM1
    LAM1 = float(val) * 0.00001
    txt_lam1_var.set(round(LAM1, 4))
    Rplot()


def update_pS(val):
    global P
    P = float(val) * 0.0001
    txt_p_var.set(round(P, 4))
    Rplot()


def update_Wl1S(val):
    global Wl1
    Wl1 = int(float(val))
    line_C11.set_data((WL[Wl1], WL[Wl1]), (0, max(I) + 250))
    line_C21.set_data((WL[Wl1], WL[Wl1]), (0, max(I) + 250))
    canvas.draw()


def update_Wl2S(val):
    global Wl2
    Wl2 = int(float(val))
    line_C12.set_data((WL[Wl2], WL[Wl2]), (0, max(I) + 250))
    line_C22.set_data((WL[Wl2], WL[Wl2]), (0, max(I) + 250))
    canvas.draw()


def Rplot():
    global Wl1, Wl2
    if Wl1 > Wl2:
        tWl = Wl1
        Wl1 = Wl2
        Wl2 = tWl

    tI = I[Wl1:Wl2]
    z = B_als(tI, LAM, LAM1, P)
    for i in range(0, (Wl2 - Wl1)):
        I_F[Wl1 + i] = z[i]

    line.set_data(WL, I_F)
    line_Final.set_data(WL, (I - I_F))
    canvas.draw()

    Config[str(WL[Wl1]) + "," + str(WL[Wl2])] = LAM, LAM1, P


def B_als(y, lam, lam1, p, niter=10):
    L = len(y)
    D = sparse.diags([1, -2, 1], [0, -1, -2], shape=(L, L - 2))
    D1 = sparse.diags([1, -1], [0, -1], shape=(L, L - 1))
    w = np.ones(L)
    for i in range(niter):
        W = sparse.spdiags(w, 0, L, L)
        Z = W.dot(W.transpose()) + lam1 * D1.dot(D1.transpose()) + lam * D.dot(D.transpose())
        b = (W.dot(W.transpose()) + lam1 * D1.dot(D1.transpose())).dot(y)
        z = spsolve(Z, b)
        w = p * (y > z) + (1 - p) * (y < z)
    return z


def file_save():
    name = filedialog.asksaveasfile(mode='w', initialfile=os.path.basename(file_path))
    B = (I - I_F)
    for i in range(len(B)):
        name.write(f"{WL[i]}\t{B[i]}\n")
    name.close()


def Ex_plot():
    plt.figure(1, figsize=(13, 3))
    plt.subplots_adjust(left=0.071, bottom=0.214, right=0.976, top=0.974)
    plt.plot(WL, (I - I_F))
    plt.xlabel('Ramanshift $(cm^{-1})$')
    plt.ylabel('Intensity(counts)')
    plt.grid();
    plt.figure(1).show()


def Ex_Para():
    name = filedialog.asksaveasfile(mode='w', initialfile='Config_' + os.path.basename(file_path))
    matrix = []
    for key, value in Config.items():
        row = [float(x) if '.' in x or x[0] == '-' and '.' in x[1:] else int(x) for x in key.split(",")] + list(
            map(float, value))
        matrix.append(row)

    for i in range(len(matrix)):
        name.write(f"{matrix[i][0]}:{matrix[i][1]}\t{matrix[i][2]},{matrix[i][3]},{matrix[i][4]}\n")

    name.close()


def update_plot():
    global file_path, WL, I, Wl1, Wl2, I_F

    file_path = filedialog.askopenfilename()
    data = np.loadtxt(file_path)
    lb4.configure(text=os.path.basename(file_path))

    WL = (data[:, 0] * m) + c
    I = data[:, 1]

    lenWL = len(WL)

    ##    b, a = signal.butter(1, 0.08)
    ####    b, a = signal.butter(5, 0.25)
    ##    b, a = signal.butter(1, 0.9999)
    ##    I = signal.filtfilt(b, a, I)

    ##    I=savgol_filter(I,51, 3)
    ##    I=I
    ##    I=savgol_filter(I,41, 2)

    Wl1 = 0
    Wl2 = len(WL) - 1

    ##    I_F=savgol_filter(I,51, 3)
    I_F = savgol_filter(I, 1, 0)

    ##    I_F=savgol_filter(I,41, 2)

    ax1.set_xlim((min(WL) - 50), (max(WL) + 50))
    ax1.set_ylim(min(I), (max(I) + 100))
    ax2.set_xlim((min(WL) - 50), (max(WL) + 50))
    ax2.set_ylim(-2, (max(I) - min(I)))

    line_In.set_data(WL, I)
    line_C11.set_data((WL[Wl1], WL[Wl1]), (0, max(I) + 250))
    line_C12.set_data((WL[Wl2], WL[Wl2]), (0, max(I) + 250))
    line_C21.set_data((WL[Wl1], WL[Wl1]), (0, max(I) + 250))
    line_C22.set_data((WL[Wl2], WL[Wl2]), (0, max(I) + 250))
    canvas.draw()


window = tk.Tk()
window.wm_title("B_als")

tabControl = ttk.Notebook(window)

frame1 = ttk.Frame(window)

frame3 = ttk.Frame(tabControl)
tabControl.add(frame3, text='Calibration')

frame4 = ttk.Frame(tabControl)
tabControl.add(frame4, text='Filter')

frame2 = ttk.Frame(tabControl)
tabControl.add(frame2, text='Base line')

frame1_1 = ttk.Frame(frame1)
frame2_1 = ttk.Frame(frame2)
frame2_2 = ttk.Frame(frame2)
frame2_3 = ttk.Frame(frame2)

Config = {}

m = 1
c = 0

LAM = 1.0
LAM1 = 0.001
P = 0.00496
lenWL = 1600

txt_lam_var = tk.StringVar()
txt_lam_var.set(round(LAM, 5))
txt_lam1_var = tk.StringVar()
txt_lam1_var.set(round(LAM1, 5))
txt_p_var = tk.StringVar()
txt_p_var.set(round(P, 5))

plt.rc('xtick', labelsize=15)
plt.rc('ytick', labelsize=15)
plt.rc('axes', labelsize=15)
plt.rc('legend', fontsize=12)

fig = plt.figure(0, figsize=(10, 5))
ax1 = plt.subplot(2, 1, 1)
plt.subplots_adjust(left=0.1, bottom=0.126, right=0.976, top=0.974)
ax1.set_xlim(0, 3510)
ax1.set_ylim(1000, 10000)
line, = ax1.plot([], [])
line_C11, = ax1.plot([], '--')
line_C12, = ax1.plot([], '--')
line_In, = ax1.plot([], [])
ax1.set_ylabel('Intensity(counts)')
ax1.grid()

ax2 = plt.subplot(2, 1, 2)
ax2.set_xlim(0, 3510)
ax2.set_ylim(-2, 10000)
line_Final, = ax2.plot([], [])
line_C21, = ax2.plot([], '--r')
line_C22, = ax2.plot([], '--g')
ax2.set_xlabel('Ramanshift $(cm^{-1})$')
ax2.set_ylabel('Intensity(counts)')
ax2.grid()

canvas = FigureCanvasTkAgg(fig, frame1_1)
toolbar = NavigationToolbar2Tk(canvas, frame1_1)
toolbar.update()

lb1 = ttk.Label(frame2_1, text='lam value')
lb2 = ttk.Label(frame2_2, text='lam1 value')
lb3 = ttk.Label(frame2_3, text='p value')

lb4 = tk.Label(frame1, text='')

txt_lam = Entry(frame2_1, width=9, textvariable=txt_lam_var)
txt_lam.bind("<Return>", update_lamE)

txt_lam1 = Entry(frame2_2, width=9, textvariable=txt_lam1_var)
txt_lam1.bind("<Return>", update_lam1E)

txt_p = Entry(frame2_3, width=9, textvariable=txt_p_var)
txt_p.bind("<Return>", update_pE)

slider_lam = ttk.Scale(frame2_1, from_=5000, to=1, length=550,
                       orient=tk.VERTICAL, command=update_lamS, value=100)
slider_lam1 = ttk.Scale(frame2_2, from_=50000, to=100, length=550,
                        orient=tk.VERTICAL, command=update_lam1S, value=100)
slider_p = ttk.Scale(frame2_3, from_=1000, to=1, length=550,
                     orient=tk.VERTICAL, command=update_pS, value=10)

##slider_o = ttk.Scale(window, from_=1, to=15, length=1000,
##                     orient=tk.HORIZONTAL, command=update_oS, value=100)
####slider_o.pack(side=tk.TOP)
####
##slider_t = ttk.Scale(window, from_=1, to=10000, length=1000,
##                     orient=tk.HORIZONTAL, command=update_tS, value=10)
####slider_t.pack(side=tk.TOP)

slider_Wl1 = ttk.Scale(frame1, from_=0, to=(lenWL - 1), length=(lenWL - 600),
                       orient=tk.HORIZONTAL, command=update_Wl1S, value=0)
slider_Wl2 = ttk.Scale(frame1, from_=0, to=(lenWL - 1), length=(lenWL - 600),
                       orient=tk.HORIZONTAL, command=update_Wl2S, value=(lenWL - 1))

canvas.get_tk_widget().pack(padx=5, side=tk.TOP, fill=tk.BOTH, expand=1)
toolbar.pack(padx=5, side=tk.BOTTOM, fill=tk.X)

slider_lam.pack(padx=0, pady=10, side=tk.TOP, fill=tk.Y, expand=1)
txt_lam.pack(padx=0, side=tk.TOP)
lb1.pack(padx=0, pady=5, side=tk.TOP)

slider_lam1.pack(padx=0, pady=10, side=tk.TOP, fill=tk.Y, expand=1)
txt_lam1.pack(padx=0, side=tk.TOP)
lb2.pack(padx=0, pady=5, side=tk.TOP)

slider_p.pack(padx=0, pady=10, side=tk.TOP, fill=tk.Y, expand=1)
txt_p.pack(padx=0, side=tk.TOP)
lb3.pack(padx=0, pady=5, side=tk.TOP)

lb4.pack(ipadx=5, side=tk.TOP, fill=tk.BOTH)

frame1_1.pack(side=tk.TOP, fill=tk.BOTH, expand=1)

slider_Wl1.pack(padx=5, pady=1, side=tk.TOP, fill=tk.X)
slider_Wl2.pack(padx=5, pady=1, side=tk.TOP, fill=tk.X)

frame2_1.pack(side=tk.LEFT, fill=tk.Y, expand=1)
frame2_2.pack(side=tk.LEFT, fill=tk.Y, expand=1)
frame2_3.pack(side=tk.LEFT, fill=tk.Y, expand=1)

frame1.pack(ipadx=10, side=tk.LEFT, fill=tk.BOTH, expand=1)

tabControl.pack(ipadx=10, side=tk.LEFT, fill=tk.Y)

menu = Menu(window)
filemenu = Menu(menu, tearoff=0)
filemenu.add_command(label='Open', command=update_plot)
filemenu.add_command(label="Save as", command=file_save)
filemenu.add_command(label="External plot", command=Ex_plot)
filemenu.add_command(label="Export Parameters", command=Ex_Para)
menu.add_cascade(label='File', menu=filemenu)
window.config(menu=menu)

window.mainloop()