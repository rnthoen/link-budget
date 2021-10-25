import tkinter as tk
from tkinter import ttk
from tkinter.messagebox import showerror
import math
import numpy as np

# Setup of GUI
root = tk.Tk()
root.title('Link Budget Tool Group 36')
root.geometry('400x300')
root.resizable(True, True)
frame = ttk.Frame(root)
options = {'padx': 5, 'pady': 5}

# Conversion function from scalar to decibels
def decibels(value):    
	db = 10 * np.log10(value)
	return db

# Calculate results when button pressed
def button_clicked():
	# Gather input values
	P_t = decibels(float(P_t_var.get()))
	
	L_t = decibels(float(L_t_var.get()))
	L_r = decibels(float(L_r_var.get()))
	
	f = float(f_var.get())
	
	D_t = float(D_t_var.get())
	D_r = float(D_r_var.get())
	
	e_t = float(e_t_var.get())
	
	T_s = float(T_s_var.get())
	
	L_a = float(L_a_var.get())
	
	inter = int(inter_var.get())
	
	d_space_sun = int(d_space_sun_var.get())
	theta = int(theta_var.get())
	h = int(h_var.get())
	h_or = int(h_or_var.get())
	
	updown = int(updown_var.get())
	
	B_p = float(B_p_var.get())
	S_wa = float(S_wa_var.get())
	P_s = float(P_s_var.get())
	R = float(R_var.get())
	GM = float(GM_var.get())
	D_c = float(D_c_var.get()) / 100
	T_dl = float(T_dl_var.get()) / 24
	
	v = np.sqrt(GM / R)
	
	R_req = float(R_req_var.get())
	
	#Constants
	k = 1.38 * 10 ** -23    #Boltzmann constant [J/K]
	mu = 0.55               #Parabolic antenna efficiency [-]
	R_E = 6371000           #Earth's radius in [m]
	c = 3e+08               #Speed of light [m/s]
	d_moon = 385000000      #Distance Earth-moon
	d_sun = 149597900000    #Distance Earth-Sun
	
	#Calculations
	lmbda = c / f           #Wavelength in [m]

	G_peak_t = 20 * np.log10(D_t) + 20 * np.log10(f/(10**9)) + 17.8     #Peak Gain in [dB] for the transmitting antenna
	G_r = decibels((np.pi**2 * D_r**2 * mu) / (lmbda**2))               #Gain in [dB] for the receiving antenna

	a_t = 21 / ((f/(10**9)) * D_t)                                      #Half-power angle for the transmitting antenna
	a_r = 21 / ((f/(10**9)) * D_r)                                      #Half-power angle for the receiving antenna
	L_pt = -12 * ((e_t / a_t)**2)                                       #Pointing loss transmitting antenna [dB]
	L_pr = -12 * ((e_t / a_r)**2)                                       #Pointing loss receiving antenna [dB]
	L_p = L_pt + L_pr                                                   #Total pointing loss [dB]

	#Determine how the ground station-spacecraft distance is calculated for the Earth orbit, interplanetary and moon orbit cases
	if inter==0:
		S = np.sqrt((R_E + h)**2 - R_E**2)
		h_or = h
	elif inter==1:
		S = np.sqrt(d_sun**2 + d_space_sun**2 - (2 * d_sun * d_space_sun * np.cos(np.radians(theta))))
	elif inter==2:
		S = d_moon

	if updown==1:
		S_w = 2 * np.tan(np.radians(0.5 * S_wa)) * h_or     #Swath width, in case of the downlink situation [m]
		R_G = (B_p * S_w * v) / (P_s ** 2)      #Generated bit rate [bit/s]
		R_req = (R_G * D_c) / T_dl              #Required bit rate [bit/s]
	#In case of uplink, the required bit rate equals the previously requirement-input.
	elif updown==0:
		R_req = R_req

	L_s = decibels((lmbda/(4 * np.pi * S))**2)  #Space losses [dB]

	R_reg_db = decibels(1 / R_req)      #Inverse of required bit rate converted to [dB]
	k_db = decibels(1 / k)              #Inverse of the Boltzmann-constant converted to [dB]
	T_s_db = decibels(1 / T_s)          #Inverse of the system noise temperature converted to [dB]

	Gain_total = G_peak_t + G_r                 #Total Gain in [dB]
	Loss_total = L_t + L_r + L_a + L_s + L_p    #Total losses in [dB]
	Noise_total = R_reg_db + k_db + T_s_db      #Total noise in [dB]

	SNR = P_t + Gain_total + Loss_total + Noise_total #Actual SNR-ratio.
	
	# Output results	
	P_t_label.config(text="P_t = " + str(round(P_t, 2)) + " [dBW]")
	Gain_total_label.config(text="Gain_total = " + str(round(Gain_total, 2)) + " [dBW]")
	Loss_total_label.config(text="Loss_total = " + str(round(Loss_total, 2)) + " [dBW]")
	Noise_total_label.config(text="Noise_total = " + str(round(Noise_total, 2)) + " [dBW]")
	SNR_label.config(text="SNR = " + str(round(SNR, 2)) + " [dBW]")
	
		

# Input variables
P_t_var		= tk.StringVar()
P_t_label 	= ttk.Label(frame, text='P_t')			.grid(column=0, row=0, sticky='W', **options)
P_t_entry 	= ttk.Entry(frame, textvariable=P_t_var)	.grid(column=1, row=0, sticky='W', **options)
P_t_unit 	= ttk.Label(frame, text='[W]')			.grid(column=2, row=0, sticky='W', **options)

L_t_var 	= tk.StringVar()
L_t_label 	= ttk.Label(frame, text='L_t')			.grid(column=0, row=1, sticky='W', **options)
L_t_entry 	= ttk.Entry(frame, textvariable=L_t_var)	.grid(column=1, row=1, sticky='W', **options)
L_t_unit 	= ttk.Label(frame, text='[-]')			.grid(column=2, row=1, sticky='W', **options)

L_r_var 	= tk.StringVar()
L_r_label 	= ttk.Label(frame, text='L_r')			.grid(column=0, row=2, sticky='W', **options)
L_r_entry 	= ttk.Entry(frame, textvariable=L_r_var)	.grid(column=1, row=2, sticky='W', **options)
L_r_unit 	= ttk.Label(frame, text='[-]')			.grid(column=2, row=2, sticky='W', **options)

f_var 		= tk.StringVar()
f_label 	= ttk.Label(frame, text='f')			.grid(column=0, row=3, sticky='W', **options)
f_entry 	= ttk.Entry(frame, textvariable=f_var)		.grid(column=1, row=3, sticky='W', **options)
f_unit 		= ttk.Label(frame, text='[Hz]')			.grid(column=2, row=3, sticky='W', **options)

D_t_var 	= tk.StringVar()
D_t_label 	= ttk.Label(frame, text='D_t')			.grid(column=0, row=4, sticky='W', **options)
D_t_entry 	= ttk.Entry(frame, textvariable=D_t_var)	.grid(column=1, row=4, sticky='W', **options)
D_t_unit	= ttk.Label(frame, text='[m]')			.grid(column=2, row=4, sticky='W', **options)

D_r_var 	= tk.StringVar()
D_r_label 	= ttk.Label(frame, text='D_r')			.grid(column=0, row=5, sticky='W', **options)
D_r_entry 	= ttk.Entry(frame, textvariable=D_r_var)	.grid(column=1, row=5, sticky='W', **options)
D_r_unit	= ttk.Label(frame, text='[m]')			.grid(column=2, row=5, sticky='W', **options)

e_t_var 	= tk.StringVar()
e_t_label 	= ttk.Label(frame, text='e_t')			.grid(column=0, row=6, sticky='W', **options)
e_t_entry 	= ttk.Entry(frame, textvariable=e_t_var)	.grid(column=1, row=6, sticky='W', **options)
e_t_unit 	= ttk.Label(frame, text='[deg]')		.grid(column=2, row=6, sticky='W', **options)

T_s_var 	= tk.StringVar()
T_s_label 	= ttk.Label(frame, text='T_s')			.grid(column=0, row=7, sticky='W', **options)
T_s_entry 	= ttk.Entry(frame, textvariable=T_s_var)	.grid(column=1, row=7, sticky='W', **options)
T_s_unit 	= ttk.Label(frame, text='[K]')			.grid(column=2, row=7, sticky='W', **options)

L_a_var 	= tk.StringVar()
L_a_label 	= ttk.Label(frame, text='L_a')			.grid(column=0, row=8, sticky='W', **options)
L_a_entry 	= ttk.Entry(frame, textvariable=L_a_var)	.grid(column=1, row=8, **options)
L_a_unit 	= ttk.Label(frame, text='[dB]')			.grid(column=2, row=8, sticky='W', **options)

inter_var 	= tk.StringVar()
inter_label 	= ttk.Label(frame, text='INTER')		.grid(column=0, row=9, sticky='W', **options)
inter_entry 	= ttk.Entry(frame, textvariable=inter_var)	.grid(column=1, row=9, sticky='W', **options)
inter_unit 	= ttk.Label(frame, text='[-]')			.grid(column=2, row=9, sticky='W', **options)

d_space_sun_var 	= tk.StringVar()
d_space_sun_label 	= ttk.Label(frame, text='d_space_sun')			.grid(column=0, row=10, sticky='W', **options)
d_space_sun_entry 	= ttk.Entry(frame, textvariable=d_space_sun_var)	.grid(column=1, row=10, sticky='W', **options)
d_space_sun_unit 	= ttk.Label(frame, text='[m]')				.grid(column=2, row=10, sticky='W', **options)

theta_var 	= tk.StringVar()
theta_label 	= ttk.Label(frame, text='theta')		.grid(column=0, row=11, sticky='W', **options)
theta_entry 	= ttk.Entry(frame, textvariable=theta_var)	.grid(column=1, row=11, sticky='W', **options)
theta_unit 	= ttk.Label(frame, text='[deg]')		.grid(column=2, row=11, sticky='W', **options)

h_var 		= tk.StringVar()
h_label 	= ttk.Label(frame, text='h')			.grid(column=0, row=12, sticky='W', **options)
h_entry 	= ttk.Entry(frame, textvariable=h_var)		.grid(column=1, row=12, sticky='W', **options)
h_unit 		= ttk.Label(frame, text='[m]')			.grid(column=2, row=12, sticky='W', **options)

h_or_var 	= tk.StringVar()
h_or_label 	= ttk.Label(frame, text='h_or')			.grid(column=0, row=13, sticky='W', **options)
h_or_entry 	= ttk.Entry(frame, textvariable=h_or_var)	.grid(column=1, row=13, sticky='W', **options)
h_or_unit 	= ttk.Label(frame, text='[m]')			.grid(column=2, row=13, sticky='W', **options)

updown_var 	= tk.StringVar()
updown_label 	= ttk.Label(frame, text='UP/DOWN')		.grid(column=0, row=14, sticky='W', **options)
updown_entry 	= ttk.Entry(frame, textvariable=updown_var)	.grid(column=1, row=14, sticky='W', **options)
updown_unit 	= ttk.Label(frame, text='[-]')			.grid(column=2, row=14, sticky='W', **options)

B_p_var 	= tk.StringVar()
B_p_label 	= ttk.Label(frame, text='B_p')			.grid(column=0, row=15, sticky='W', **options)
B_p_entry 	= ttk.Entry(frame, textvariable=B_p_var)	.grid(column=1, row=15, sticky='W', **options)
B_p_unit 	= ttk.Label(frame, text='[bit/px]')		.grid(column=2, row=15, sticky='W', **options)

S_wa_var 	= tk.StringVar()
S_wa_label 	= ttk.Label(frame, text='S_wa')			.grid(column=0, row=16, sticky='W', **options)
S_wa_entry 	= ttk.Entry(frame, textvariable=S_wa_var)	.grid(column=1, row=16, sticky='W', **options)
S_wa_unit 	= ttk.Label(frame, text='[deg]')		.grid(column=2, row=16, sticky='W', **options)

P_s_var 	= tk.StringVar()
P_s_label 	= ttk.Label(frame, text='P_s')			.grid(column=0, row=17, sticky='W', **options)
P_s_entry 	= ttk.Entry(frame, textvariable=P_s_var)	.grid(column=1, row=17, sticky='W', **options)
P_s_unit 	= ttk.Label(frame, text='[m]')			.grid(column=2, row=17, sticky='W', **options)

R_var 		= tk.StringVar()
R_label 	= ttk.Label(frame, text='R')			.grid(column=0, row=18, sticky='W', **options)
R_entry 	= ttk.Entry(frame, textvariable=R_var)		.grid(column=1, row=18, sticky='W', **options)
R_unit 		= ttk.Label(frame, text='[m]')			.grid(column=2, row=18, sticky='W', **options)

GM_var 		= tk.StringVar()
GM_label 	= ttk.Label(frame, text='GM')			.grid(column=0, row=19, sticky='W', **options)
GM_entry 	= ttk.Entry(frame, textvariable=GM_var)		.grid(column=1, row=19, sticky='W', **options)
GM_unit 	= ttk.Label(frame, text='[m^3/s^2]')		.grid(column=2, row=19, sticky='W', **options)

D_c_var 	= tk.StringVar()
D_c_label 	= ttk.Label(frame, text='D_c')			.grid(column=0, row=20, sticky='W', **options)
D_c_entry 	= ttk.Entry(frame, textvariable=D_c_var)	.grid(column=1, row=20, sticky='W', **options)
D_c_unit 	= ttk.Label(frame, text='[%]')			.grid(column=2, row=20, sticky='W', **options)

T_dl_var 	= tk.StringVar()
T_dl_label 	= ttk.Label(frame, text='T_dl')			.grid(column=0, row=21, sticky='W', **options)
T_dl_entry 	= ttk.Entry(frame, textvariable=T_dl_var)	.grid(column=1, row=21, sticky='W', **options)
T_dl_unit 	= ttk.Label(frame, text='[h]')			.grid(column=2, row=21, sticky='W', **options)

R_req_var 	= tk.StringVar()
R_req_label 	= ttk.Label(frame, text='R_req')		.grid(column=0, row=22, sticky='W', **options)
R_req_entry 	= ttk.Entry(frame, textvariable=R_req_var)	.grid(column=1, row=22, sticky='W', **options)
R_req_unit 	= ttk.Label(frame, text='[bit/s]')		.grid(column=2, row=22, sticky='W', **options)

# Outputs
P_t_label = ttk.Label(frame)
P_t_label.grid(row=101, column=0, sticky='W', **options)

Gain_total_label = ttk.Label(frame)
Gain_total_label.grid(row=102, column=0, sticky='W', **options)

Loss_total_label = ttk.Label(frame)
Loss_total_label.grid(row=103, column=0, sticky='W', **options)

Noise_total_label = ttk.Label(frame)
Noise_total_label.grid(row=104, column=0, sticky='W', **options)

SNR_label = ttk.Label(frame)
SNR_label.grid(row=105, column=0, sticky='W', **options)

	
# Button
convert_button = ttk.Button(frame, text='Calculate')
convert_button.grid(column=1, row=100, sticky='W', **options)
convert_button.configure(command=button_clicked)

# Run GUI
frame.grid(padx=5, pady=5)
root.mainloop()

