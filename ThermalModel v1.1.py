# -*- coding: utf-8 -*-
"""
Created on Mon May 10 13:40:04 2021

@author: ybshi
"""

import numpy as np
import sympy as sp
import math as mt
import timeit
import pandas as pd

from scipy.integrate import odeint
from matplotlib import pyplot as plt
from CellVoltage import Cell_Voltage as CV
from HydrogenProduction import Hydrogen_Production as HP

class Thermal_Model():
    
    def __init__(self, T, P):
        
        self.F = 96485.34
        self.R = 82.06 # atm cm3 / mol K
        self.T0 = 298
        self.T = T
        self.P = P
        
    def Thermoneutral_Voltage(self):
        """ Function for calculating the thermoneutral cell voltage
            Ref. K. Onda et al., Journal of Power Sources, 132 (2004), pp.64-70.
        """
        
        # Coefficients for specfic heat of water, hydrogen and oxygen [a, b, c, e]
        Coefficients = np.array(([72.39, 9.38, 0, 0], [26.57, 3.77, 1.17, 0], [34.35, 1.92, -18.45, 4.06]))
        # Enthalpy calculation according to temperature with coefficients
        ΔH = np.array(([(self.T-self.T0)], [0.5*10**(-3)*((self.T**2)-(self.T0**2))], [-(10**5)*((1/self.T)-(1/self.T0))], [-0.5*(10**8)*((1/(self.T**2))-(1/(self.T0**2)))]))        
        
        # Enthalpy of water splitting (Temperature change)
        Enthalpy_T = -285830 - (Coefficients.dot(ΔH)[1,0] + 0.5*Coefficients.dot(ΔH)[2,0] - Coefficients.dot(ΔH)[0,0])
        
        # Define the symbol of temperature (T)
        T = sp.symbols('T', real=True)
        
        # Coefficients for virial constant of hydrogen and oxygen [b1, b2, c1, c2]
        Coefficient_B = np.array(([20.5, -1857], [42.6, -17400]))
        Coefficient_C = np.array(([-351, 12760], [-2604, 61457]))
        
        # Virial constants
        B = np.array(([1], [1/T]))
        C = np.array(([1], [1/(T**0.5)]))
        
        # ∂B/∂T and ∂C/∂T
        Derivative_B = sp.derive_by_array(Coefficient_B.dot(B), T)
        Derivative_C = sp.derive_by_array(Coefficient_C.dot(C), T)
        
        # Enthalpy of water splitting (Pressure change)
        Enthalpy_P = 0.101325 * (((Coefficient_B.dot(B)[0,0].subs(T, self.T) - self.T * Derivative_B[0,0].subs(T, self.T)) * self.P 
                        + (self.P**2) * (Coefficient_C.dot(C)[0,0].subs(T, self.T) - Coefficient_B.dot(B)[0,0].subs(T, self.T) * Coefficient_B.dot(B)[0,0].subs(T, self.T)
                        - 0.5 * self.T * (Derivative_C[0,0].subs(T, self.T) - 2 * Coefficient_B.dot(B)[0,0].subs(T, self.T) * Derivative_B[0,0].subs(T, self.T)))/(self.R*self.T))
                      + 0.5 * ((Coefficient_B.dot(B)[1,0].subs(T, self.T) - self.T * Derivative_B[1,0].subs(T, self.T)) * self.P 
                        + (self.P**2) * (Coefficient_C.dot(C)[1,0].subs(T, self.T) - Coefficient_B.dot(B)[1,0].subs(T, self.T) * Coefficient_B.dot(B)[1,0].subs(T, self.T)
                        - 0.5 * self.T * (Derivative_C[1,0].subs(T, self.T) - 2 * Coefficient_B.dot(B)[1,0].subs(T, self.T) * Derivative_B[1,0].subs(T, self.T)))/(self.R*self.T)))
        
        # Enthalpy of water splitting according to temperature and pressure
        Enthalpy_T_P = Enthalpy_T - Enthalpy_P
        
        # Thermoneutral cell voltage
        ThermoVoltage = (-Enthalpy_T_P)/(2*self.F)
        
        return ThermoVoltage
    
    def Temperature(self, current, Time_final):
        """ Function for calculating the temperature
            Ref. O. Ulleberg, International Journal of Hydrogen Energy, 28 (2003), pp.21-33.
            Results: Time and Temperature
        """
        
        def model(T, t, I):
            """ Define ordinary differential equation """

            # Parameters
            C = 625000 # Overall thermal capacity of electrolyzer
            
            # Heat generation (Q_gen)
            Nc = 21 # Number of cell in series per stack
            V = CV(2, T[0], self.P, 30, I).CV_Jang(0.25, 0.2, 0) # Cell voltage
            V_tn = Thermal_Model(T[0], self.P).Thermoneutral_Voltage() # Thermoneutral voltage

            # Heat loss to air (Q_loss)
            R = 0.167 # Resistance of the elctrolyzer
            Ta = 293 # Ambient temperature
            
            # Cooling (Q_cool)
            h_conduction = 7 # Conduction heat transfer
            h_convection = 0.02 # Convection heat transfer
            Cw = 697.9 # Heat capacity of cooling water (In case of 0.6 Nm3/h, heat capacity is 0.6975 J/C s)
            UA = h_conduction + (h_convection * I)
            T_w_i = 287.5
            T_w_o = T_w_i + (T - T_w_i) * (1 - mt.exp(-(UA/Cw)))
            LMTD = ((T - T_w_i) - (T - T_w_o)) / mt.log(((T - T_w_i) / (T - T_w_o)))
            
            # Heat balance
            Q_gen = int((Nc * V * I) * (1 - (V_tn/V)))
            Q_loss = (1/R) * (T - Ta)
            Q_cool = UA * LMTD
            
            dTdt = (Q_gen - Q_loss - Q_cool) / C

            return dTdt
        
        # Initial condition of temperature
        Initial_Tem = self.T
        SimulationTime = Time_final
        
        # Time, current and temperature list for data save
        TimeList = []
        Tem = []
        
        for i in range(SimulationTime):
            
            # Simulation time
            Time_in = i
            Time_fi = i + 1
            Time_step = 0.1
            NumDataPoint = int(((Time_fi-Time_in)/Time_step))
            Time = np.linspace(Time_in, Time_fi, NumDataPoint)
            
            # Current variation of electrolyzer
            I = current
            
            # Solution of ordinary differential equation
            y = odeint(model, Initial_Tem, Time, args=(I[i],))
            
            for j in range(len(Time) - 1):
                TimeList.append(Time[j])
                Tem.append(y[j, 0])
            
            # Initial temperature
            Initial_Tem = y[-1, -1]
            
        TimeList.append(Time[-1])
        Tem.append(y[-1, -1])

        return TimeList, Tem

""" Simulation Setup """

# Simulation Parameters
SimulationTime = 25200 # Unit [s]
InitialTemperature = 303 # Unit [K]
Pressure = 7 # Unit [bar]
NumOfCell = 21

X = np.full((1800), 100)
X1 = np.linspace(0,1*np.pi,3600)
# X2 = np.linspace(0,0.5*np.pi,600)
# X3 = np.linspace(0,0.5*np.pi,500)
# X4 = np.linspace(0,0.5*np.pi,500)
# X5 = np.linspace(0,4*np.pi,500)
# X6 = np.linspace(0,0.5*np.pi,500)
X7 = np.full((3600), 100)

Profile1 = X
Profile2 = 100 + 700 * np.sin(X1)
# Profile3 = 600 + 350 * np.sin(-X2)
# Profile4 = 400 + 300 * np.sin(X3)
# Profile5 = 700 + 300 * np.sin(-X4)
# Profile6 = 400 + 200 * np.sin(X5)
# Profile7 = 400 + 300 * np.sin(-X6)
Profile8 = X7
Profile = np.concatenate((Profile1, Profile2, Profile1, Profile2, Profile1, Profile2, Profile1, Profile2, Profile8))
# plt.plot(Profile)
CurrentProfile = Profile

""" Simulation """

start_time = timeit.default_timer() # 타이머 시작
Results = Thermal_Model(InitialTemperature, Pressure).Temperature(CurrentProfile, SimulationTime)
terminate_time = timeit.default_timer() # 타이머 종료

SimTime = np.array(Results[0])
Temperature = np.array(Results[1]) - 273 # K to C

HydrogenProduction = []
ThermoneutralVoltage = []
Cellvoltage = []
Current = []
Current.append(CurrentProfile[0])

# Current profile list
for c in range(len(CurrentProfile)):
        
    for t in range(len(SimTime)):

        Simultaion = SimTime[t]
        
        if Simultaion > c + 1:
            break
        
        elif c < Simultaion <= c + 1:
            Current.append(CurrentProfile[c])

# Cell voltage profile list
for i in range(len(Results[1])):
    
    Tem = np.array(Results[1])
    V = CV(2, Tem[i], Pressure, 30, Current[i]).CV_Jang(0.25, 0.2, 0) * NumOfCell # Cell voltage
    V_tn = Thermal_Model(Tem[i], Pressure).Thermoneutral_Voltage() * NumOfCell # Thermoneutral voltage
    
    # Append the simulation data to list
    Cellvoltage.append(V)
    ThermoneutralVoltage.append(V_tn)

for i in range(len(Results[1])):
    
    Tem = np.array(Results[1])
    Hydrogen = HP(Tem[i], Pressure, Current[i], 21).HP_Jang()
    HydrogenProduction.append(Hydrogen)

# Faraday efficiency
f = 2.5 * Tem - 632.5
f2 = -0.00075 * Tem + 1.22475
η = ((((np.array(Current)*1000)/2500)**2) / (f + (((np.array(Current)*1000))/2500)**2)) * f2 * 100 
        
time = "{0:.2f}".format(terminate_time - start_time)
print("Time :", time, "초")

# Time - Temperature curve
plt.plot(SimTime, Temperature, color='indigo', marker='o', markevery=1500)
plt.xlim(0, SimulationTime)
plt.xlabel("Time [s]")
plt.ylabel("Temperature [°C]")
plt.grid()
plt.show()

# # Time - Hydrogen Production curve
# plt.plot(SimTime, HydrogenProduction, color='b')
# plt.xlim(0, SimulationTime)
# plt.xlabel("Time [s]")
# plt.ylabel("Hydrogen Production [Nm$^3$/h]")
# plt.grid()
# plt.show()

blank = []
fig, ax1 = plt.subplots()
ax1.plot(SimTime, HydrogenProduction, label="Hydrogen Production rate", color='blueviolet', marker='s', markevery=1500)
ax1.plot(blank, label="Faraday efficiency", color='darkblue', marker='^', markevery=1500)
ax1.set_xlabel("Time [s]")
ax1.set_ylabel("Hydrogen Production [Nm$^3$/h]")
ax1.set_xlim([0, SimulationTime])
ax1.legend(loc='lower right')
ax1.set_ylim([0, 10])
ax1.xaxis.grid(True, color='gray', linestyle='dashed', linewidth=0.5)
ax1.yaxis.grid(True, color='gray', linestyle='dashed', linewidth=0.5)
ax2 = ax1.twinx()
ax2.plot(SimTime, η, label="Faraday efficiency", color='darkblue', marker='^', markevery=1500)
ax2.set_ylabel("Efficiency [%]")
ax2.set_ylim([0, 100])
plt.show()

# Current and cell voltage plotting
blank = []
fig, ax1 = plt.subplots()
ax1.plot(SimTime, Cellvoltage, label="Cell voltage", color='maroon', marker='s', markevery=1500)
ax1.plot(SimTime, ThermoneutralVoltage, label="Thermoneutral voltage", color='violet', marker='v', markevery=1500)
ax1.plot(blank, label="Current", color='darkgreen', marker='H', markevery=1500)
ax1.set_xlabel("Time [s]")
ax1.set_ylabel("Voltage [V]")
ax1.set_xlim([0, SimulationTime])
ax1.set_ylim([30, 55])
ax1.xaxis.grid(True, color='gray', linestyle='dashed', linewidth=0.5)
ax1.yaxis.grid(True, color='gray', linestyle='dashed', linewidth=0.5)
ax1.legend(loc='upper right')
ax2 = ax1.twinx()
ax2.plot(SimTime, Current, label="Current", color='darkgreen',  marker='H', markevery=1500)
ax2.set_ylabel("Current [A]")
ax2.set_ylim([0, 1000])
# ax2.legend(loc='upper left')
plt.show()

HydrogenProduction = np.array(HydrogenProduction)
Cellvoltage = np.array(Cellvoltage)
ThermoneutralVoltage = np.array(ThermoneutralVoltage)
Current = np.array(Current)

TimeArray = SimTime.reshape(-1,1)
TemArray = Temperature.reshape(-1,1)
HydrogeArray = HydrogenProduction.reshape(-1,1)
EffiArray = η.reshape(-1,1)
VoltageArray = Cellvoltage.reshape(-1,1)
ThermoneutralArray = ThermoneutralVoltage.reshape(-1,1)
CurrentArray = Current.reshape(-1,1)

Data = np.hstack([TimeArray, TemArray, HydrogeArray, EffiArray, VoltageArray, ThermoneutralArray, CurrentArray])
DataSave = pd.DataFrame(Data, columns = ["Time", "Temperature", "Hydrogen", "Efficiency", "Voltage", "Thermoneutral", "Current"])
DataSave.to_excel("G:/연구/Digital Twin for Alkaline Water Eletrolysis/Simulation Results/SimResult(Case1).xlsx")
