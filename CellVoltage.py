# -*- coding: utf-8 -*-
"""
Created on Fri Apr  9 14:16:33 2021

@author: ybshi
"""

import math
from ReversibleVoltage import Reversible_Voltage as RV
from ReversibleVoltage import Gibbs_Energy as GE
from matplotlib import pyplot as plt

""" 
알칼라인 수전해 시스템의 Electrochemical 모델링을 위한 온도, 압력에 따른 Cell voltage를 계산하는 코드
Cell Voltage = Reversible voltage (Open-circuit voltage) + Activation overpotential + Ohmic overpotential + (Concentration overpotential)
"""

class Cell_Voltage(RV):
    """ * Reversible voltage equation number 
    1. Basic Model (Constant P) 2. LeRoy Model 3. Onda model 4. Hammoudi model """
    
    def __init__(self, Equation, T, P):
        
        super(Cell_Voltage, self).__init__()
        
        if Equation == 1: # Basic Model (Constant P)
            
            self.T = T
            self.RCV = RV().RV_ConstantP(self.T)
            
        elif Equation == 2: # LeRoy Model
            
            w = float(input("KOH 전해질의 wt%를 입력하시오. \nValue (%): "))
            self.T = T
            self.P = P
            self.RCV = RV().RV_LeRoy(self.T, self.P, w)
            
        elif Equation == 3: # Onda model
            
            self.T = T
            self.P = P
            self.RCV = RV().RV_Onda(self.T, self.P)
            
        elif Equation == 4: # Hammoudi model
            
            w = float(input("KOH 전해질의 wt%를 입력하시오. \nValue (%): "))
            self.T = T
            self.P = P
            self.RCV = RV().RV_Hammoudi(self.T, self.P, w)
            
        else:
            print("Out of bound")
            
    def CV_Ulleberg(self, I, A, r1, r2, t1, t2, t3, s):
        """
        2003년 Ulleberg 의해 제안된 모델 
        Ref. : O. Ulleberg, "Modeling of Advanced Alkaline Electrolyzers: a system simulation approach, Int. J. Hydrog. Energy, 28, pp. 21-33 (2003)"
        * 모델의 특징: Stationary electrolyte를 가정하며, ohmic resistance parameter와 overvoltage coefficients로 비교적 간단하게 나타낸 모델.
        * 필요 parameters
        I: 전류 [A], A: 전극의 면적 [m2], r: ohmic resistance parameters, s and t: coefficient for overvoltage on electrode
        * Parameters의 단위
        r1: ohm m2, r2: ohm m2/c, s: V, t1: m2/A, t2: m2 C/A, t3: m2 c2/A
        """
        
        Tem = self.T - 273
        CV = self.RCV + ((r1 + r2*Tem)/A) * I + s * math.log10((((t1 + (t2/Tem) + (t3/(Tem**2)))/A) * I) + 1 )
        
        return CV
    
    def CV_Henao(self, I, A, w, L_a, L_c, d_a, d_c):
        """
        2014년 Henao 의해 제안된 모델 
        Ref. : C. Henao et al., "Simulation Tool Based on a Physics Model and an Electrical Analogy for an Alkaline Electrolyser", J. Power Sources, 250, pp. 58-67 (2014).
        * 필요 parameters
        I: 전류[A], A: 전극의 총 표면적 [m2], A_m: Membrane 표면적 [m2], w: 전해질 질량% [%], L_a and L_c: 전극 길이 (Anode, Cathode) [cm]
        d_a and d_c: Anode-membrane and Cathode-membrane 거리 [cm]
        """
        
        # Theta: bubble coverage
        θ = (-97.25 + 182 * (self.T/self.T0) - 84 * ((self.T/self.T0)**2)) * ((I/(0.5*A)) / 300000)**0.3 # Fraction of electrode coverage by the gas bubbles
        # Effective electrode surface [cm2]
        S_eff = (0.5 * A * 1000) * (1 - θ) # Electrode surface except coverage [cm2]
        # Current density [mA/cm2]
        J = (I*1000) / S_eff
            
        def AO_Cathode():
            """ Activation Overpotential on Cathode """
            
            # Transfer coefficient
            a = 0.1175 + 0.00095 * self.T # 음극의 전하전달계수
            b = (2.303 * self.T * self.R) / (self.z * self.F * a)
            # Exchange current densities [mA/cm2]
            J0 = 13.72491 - 0.09055 * self.T + 0.009055 * (self.T**2)
            
            # Activation overpotential of Cathode [V]
            if I == 0: 
                AO_Ca = 0
                
            else:
                AO_Ca = b * math.log10(J/J0) - b * math.log10(1-θ)
            
            return AO_Ca
        
        def AO_Anode():
            """ Activation Overpotential on Anode """
            
            # Transfer coefficient
            a = 0.0675 + 0.00095 * self.T # 양극의 전하전달계수
            b = (2.303 * self.T * self.R) / (self.z * self.F * a)
            # Exchange current densities [mA/cm2]
            J0 = 30.4 - 0.206 * self.T + 0.00035 * (self.T**2) # mA/cm2
            
            # Activation overpotential of Anode [V]
            if I == 0: 
                AO_An = 0
                
            else:
                AO_An = b * math.log10(J/J0) - b * math.log10(1-θ)
            
            return AO_An

        def Ohm_Resistance_Electrodes():
            
            # Electrical conductivity of the Ni [S/cm]
            σ = (60000000 - 279650 * self.T + 532 * self.T**2 - 0.38057 * self.T**3) * 0.01
            
            # Ohm resistance of Cathode and Anode [S]
            Ohm_Cathode = (1/σ) * (L_c/(0.5*A*10000))
            Ohm_Anode = (1/σ) * (L_a/(0.5*A*10000))
            
            # Ohm resistance of electrodes [S]
            Ohm_Res_Electrodes = Ohm_Cathode + Ohm_Anode
            
            return Ohm_Res_Electrodes
        
        def Ohm_Resistance_Electrolyte():
            
            # wt% KOH soultion molar concentration [mol/m3]
            molarity = ((w * 10) * 997.03 * math.exp(0.0086*w)) / 56.1
            # Ionic conductivity of the KOH [S/cm]
            σ = -2.04 * molarity - 0.0028 * (molarity**2) + 0.005332 * molarity * self.T + 207.2 * (molarity/self.T) + 0.001043 * (molarity**3) - 0.0000003 * (molarity**2) * (self.T**2)
            ε = (2/3) * θ
            
            # Ohm resistance of bubble-free electrolyte [S]
            Ohm_Ele_free = (1/σ) * ((d_a/(0.5 * A * 10000)) + (d_c/(0.5 * A * 10000)))
            
            # Ohm resistance of bubble coverage electrolyte [S]
            if I == 0:
                Ohm_Ele_bubble = 0
                
            else:
                Ohm_Ele_bubble = Ohm_Ele_free * ((1/((1-ε)**1.5)) - 1)
             
            # Ohm resistance of electrolyte [S]
            Ohm_Electrolyte = Ohm_Ele_free + Ohm_Ele_bubble
            
            return Ohm_Electrolyte
        
        def Ohm_Resistance_Membrane():
            
            # Ohm resistance of electrolyte [S]
            Ohm_Res_Membrane = (0.060 + 80 * (0.5**(self.T/50))) / (10000 * A)
            
            return Ohm_Res_Membrane
        
        # Ohm resistance [S]
        Resistance = Ohm_Resistance_Electrodes() + Ohm_Resistance_Electrolyte() + Ohm_Resistance_Membrane()
        
        CV = self.RCV + AO_Cathode() + AO_Anode() + Resistance * I
        
        return CV
    
    def CV_Abdin(self, I, A, w, α_Ca, α_An, r_Ca, r_An, i0_Ca_Ref, i0_An_Ref, G_Ca, G_An, ε_Electrode, ε_Sep, 𝛿_Ele, 𝛿_Sep, l, β, τ, ω):
        """
        2017년 Abdin 의해 제안된 모델 
        Ref. : Z. Abdin et al., "Modelling and Simulation of an Alkaline Electrolyser Cell", Energy, 138, pp. 316-331 (2017).
        * 필요 parameters
        I: 전류[A], A: 총 전극의 표면적 [m2], w: KOH 전해질 농도 [%], α_Ca and α_An: Electrode charge transfer coefficient,
        r_Ca and r_An: Roughness factor, i0_Ca_Ref and i0_An_Ref: reference exchange current density [A/cm2],
        G_Ca, G_An: Gibbs free energy [kJ/mol], ε_Electrode and ε_Electrolyte: Porosity,
        𝛿_Ele and 𝛿_Sep: thickness [cm], l: length between electrode and seperator [cm], β: width of bubble zone [cm], τ: Effective length, ω: Wettability
        """
        
        # i: Current density[A/cm2]
        i = I / (A * 1000)
        # molality: KOH molality [mol/kg]
        molality = ((w *10) / 56.1056) / 0.7
        # Vapour_Pressure: Vapour pressure of the KOH solution [bar]
        Vapour_Pressure = 10**(-0.01508 * molality - 0.0016788 * (molality**2) + 2.25887 * (10**-5) * (molality**3) 
                           + (1 - 0.0012062 * molality + 5.6024 * (10**-4) * (molality**2) - 7.8228 * (10**-6) * (molality**3)) 
                           * (35.4462 - (3343.93/self.T) - 10.9 * math.log10(self.T) + 0.0041645*self.T))
        # Theta: bubble coverage
        θ = (-97.25 + 182 * (self.T/self.T0) - 84 * ((self.T/self.T0)**2)) * (((I/(0.5*A)) / 300000)**0.3) * (self.P / (self.P-Vapour_Pressure))
        
        def Activation_Overpotential():
            # Exchange current densities [mA/cm2]
            i0_Cathode = r_Ca * math.exp((-G_Ca*1000/self.R)*((1/self.T) - (1/self.T0))) * i0_Ca_Ref
            i0_Anode = r_An * math.exp((-G_An*1000/self.R)*((1/self.T) - (1/self.T0))) * i0_An_Ref
            
            if I == 0:
                AO_An = 0
                AO_Ca = 0
                
            else:
                AO_Ca = ((self.R * self.T) / (α_Ca * self.F)) * math.log(i/(i0_Cathode*(1-θ)))
                AO_An = ((self.R * self.T) / (α_An * self.F)) * math.log(i/(i0_Anode*(1-θ)))
                
            return AO_An + AO_Ca
        
        def Ohm_Resistance():
            
            # Resistivity of the Nickel (Material) [S cm]
            ρ_Electrode = 6.99 * (10**-8)
            # Temperature coefficient of Nickel
            κ = 0.006
            
            # Resistance of Electrode [S]
            Resistance_Cathode = (ρ_Electrode / ((1-ε_Electrode)**1.5)) * (𝛿_Ele/(A * 500)) * (1+κ*(self.T-293))
            Resistance_Anode = (ρ_Electrode / ((1-ε_Electrode)**1.5)) * (𝛿_Ele/(A * 500)) * (1+κ*(self.T-293))
            Resistance_Electrode = Resistance_Cathode + Resistance_Anode
            
            # Resistance of Electrolyte [S]
            # Ohm resistance of bubble-free electrolyte [S]
            
            # wt% KOH soultion molar concentration [mol/m3]
            molarity = ((w * 10) * 997.03 * math.exp(0.0086*w)) / 56.1
            # Ionic conductivity of the KOH [S/cm]
            σ = -2.04 * molarity - 0.0028 * (molarity**2) + 0.005332 * molarity * self.T + 207.2 * (molarity/self.T) + 0.001043 * (molarity**3) - 0.0000003 * (molarity**2) * (self.T**2)
            ε = (2/3) * θ
            
            Ohm_Ele_free = (1/σ) * ((l/(0.5 * A * 10000)) + (l/(0.5 * A * 10000)))
            
            # Ohm resistance of bubble coverage electrolyte [S]
            if I == 0:
                Ohm_Ele_bubble = 0
                
            else:
                Ohm_Ele_bubble = Ohm_Ele_free * ((1/((1-ε)**1.5)) - 1)
             
            # Ohm resistance of electrolyte [S]
            Resistance_Electrolyte = Ohm_Ele_free + Ohm_Ele_bubble
            
            # Resistance of Separator [S]
            Resistance_Separator = (0.060 + 80 * (0.5**(self.T/50))) / (10000 * A)
            
            ResistanceOverpotential = (Resistance_Electrode + Resistance_Separator + Resistance_Electrolyte) * I
            
            return ResistanceOverpotential
        
        CV = self.RCV + Activation_Overpotential() + Ohm_Resistance()
        
        return CV
    
    # def CV_Jang(self):
    #     """
    #     2021년 Jang 의해 제안된 모델 
    #     Ref. : D. Jang et al., "Numerical Modeling and Analysis of the Effect of Pressure on the Performance of an Alkaline Water Electrolysis System", Applied Energy, 287, 116554 (2021).
    #     * 필요 parameters
        
    #     * Parameters의 단위
        
    #     """
        
    #     CV = self.RCV
        
    #     return CV