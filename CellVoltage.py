# -*- coding: utf-8 -*-
"""
Created on Fri Apr  9 14:16:33 2021

@author: ybshi
"""

import math
from ReversibleVoltage import Reversible_Voltage as RV
import scipy as sp

""" 
알칼라인 수전해 시스템의 Electrochemical 모델링을 위한 온도, 압력에 따른 Cell voltage를 계산하는 코드
Cell Voltage = Reversible voltage (Open-circuit voltage) + Activation overpotential + Ohmic overpotential + (Concentration overpotential)
"""

class Cell_Voltage(RV):
    """ * Reversible voltage equation number 
    1. Basic Model (Constant P) 2. LeRoy Model 3. Onda model 4. Hammoudi model """
    
    def __init__(self, Equation, T, P, W, I):
        
        super(Cell_Voltage, self).__init__()
        
        if Equation == 1: # Basic Model (Constant P)
            
            self.T = T
            self.W = W
            self.I = I
            self.RCV = RV().RV_ConstantP(self.T)
            
        elif Equation == 2: # LeRoy Model
            
            self.T = T
            self.P = P
            self.W = W
            self.I = I
            self.RCV = RV().RV_LeRoy(self.T, self.P, W)
            
        elif Equation == 3: # Onda model
            
            self.T = T
            self.P = P
            self.W = W
            self.I = I
            self.RCV = RV().RV_Onda(self.T, self.P)
            
        elif Equation == 4: # Hammoudi model
            
            self.T = T
            self.P = P
            self.W = W
            self.I = I
            self.RCV = RV().RV_Hammoudi(self.T, self.P, W)
            
        else:
            print("Out of bound")
            
    def CV_Ulleberg(self, A, r1, r2, t1, t2, t3, s):
        """
        2003년 Ulleberg 의해 제안된 모델 
        Ref. : O. Ulleberg, "Modeling of Advanced Alkaline Electrolyzers: a system simulation approach, Int. J. Hydrog. Energy, 28, pp. 21-33 (2003)"
        * 모델의 특징: Stationary electrolyte를 가정하며, ohmic resistance parameter와 overvoltage coefficients로 비교적 간단하게 나타낸 모델.
        * 필요 parameters
        I: 전류 [A], A: 전극의 면적 [m2], r: ohmic resistance parameters, s and t: coefficient for overvoltage on electrode
        * Parameters의 단위
        r1: ohm m2, r2: ohm m2/c, s: V, t1: m2/A, t2: m2 C/A, t3: m2 c2/A
        """
        
        CV = self.RCV + ((r1 + r2*self.T)/A) * self.I + s * math.log10((((t1 + (t2/self.T) + (t3/(self.T**2)))/A) * self.I) + 1 )
        
        return CV
    
    def CV_Henao(self, A, L_a, L_c, d_a, d_c):
        """
        2014년 Henao 의해 제안된 모델 
        Ref. : C. Henao et al., "Simulation Tool Based on a Physics Model and an Electrical Analogy for an Alkaline Electrolyser", J. Power Sources, 250, pp. 58-67 (2014).
        * 필요 parameters
        I: 전류[A], A: 전극의 총 표면적 [m2], A_m: Membrane 표면적 [m2], w: 전해질 질량% [%], L_a and L_c: 전극 길이 (Anode, Cathode) [cm]
        d_a and d_c: Anode-membrane and Cathode-membrane 거리 [cm]
        """
        
        # Theta: bubble coverage (Fraction of electrode coverage by the gas bubbles)
        θ = (-97.25 + 182 * (self.T/self.T0) - 84 * ((self.T/self.T0)**2)) * ((self.I/A) / 300000)**0.3
        
        # Effective electrode surface [cm2]
        S_eff = (A * 10000) * (1 - θ) # Electrode surface except coverage [cm2]
        
        # Current density [mA/cm2]
        j = (self.I * 1000) / S_eff
            
        def Activation_Overportential():
            """ Activation Overpotential of Electrodes """
            
            # Charge Transfer Coefficient
            α_An = 0.0675 + 0.00095 * self.T
            α_Ca = 0.1175 + 0.00095 * self.T
            
            # Coefficient
            b_An = (2.303 * self.T * self.R) / (self.z * self.F * α_An)
            b_Ca = (2.303 * self.T * self.R) / (self.z * self.F * α_Ca)
            
            # Exchange current densities [mA/cm2]
            j_An = 30.4 - 0.206 * self.T + 0.00035 * (self.T**2)
            j_Ca = 13.72491 - 0.09055 * self.T + 0.0001555 * (self.T**2)
            
            if self.I == 0:
                AO = 0
                
            else:
                # Activation Overpotential (Anode + Cathode) [V]
                AO_An = b_An * math.log10(j/j_An) - b_An * math.log10(1-θ)
                AO_Ca = b_Ca * math.log10(j/j_Ca) - b_An * math.log10(1-θ)
                AO = AO_An + AO_Ca
                 
            return AO

        def Ohm_Resistance():
            
            # Ohm Resistance of Electrodes
            
            # Electrical conductivity of the Ni [S/cm]
            σ_Ni = 60000000 - 279650 * self.T + 532 * self.T**2 - 0.38057 * self.T**3 
            # Ohm resistance of Cathode and Anode [S]
            Ohm_Cathode = (1/σ_Ni) * (L_c/(A*10000))
            Ohm_Anode = (1/σ_Ni) * (L_a/(A*10000)) 
            # Ohm resistance of electrodes [S]
            Ohm_Res_Electrodes = Ohm_Cathode + Ohm_Anode

            # Ohm Resistance of Electrolyte
            
            # Molar Concentration of KOH Solution
            m = (self.W * (183.1221 - 0.56845 * self.T + 984.5679 * math.exp(self.W/115.96277))) / 5610.5
            # Conductivity of the KOH Solution [S/cm]
            σ_KOH = -2.04 * m - 0.0028 * (m**2) + 0.005332 * m * self.T + 207.2 * (m/self.T) + 0.001043 * (m**3) - 0.0000003 * (m**2) * (self.T**2)
            ε = (2/3) * θ
            # Ohm resistance of bubble-free electrolyte [S]
            Ohm_Ele_free = (1/σ_KOH) * ((d_a/(A * 10000)) + (d_c/(A * 10000)))
            # Ohm resistance of bubble coverage electrolyte [S]
            Ohm_Electrolyte = Ohm_Ele_free * ((1/((1-ε)**(2/3))) - 1)
            
            # Ohm resistance of electrolyte [S]
            
            Ohm_Res_Membrane = (0.060 + 80 * (0.5**(self.T/50))) / (10000 * A)
            
            # Ohm Overpotential (V=IR)
            Ohm_Overpotential = (Ohm_Res_Electrodes + Ohm_Electrolyte + Ohm_Res_Membrane) * self.I
            
            return Ohm_Overpotential
        
        CV = self.RCV + Activation_Overportential() + Ohm_Resistance()

        return CV
    
    def CV_Abdin(self, A, α_Ca, α_An, r_Ca, r_An, i0_Ca_Ref, i0_An_Ref, G_Ca, G_An, ε_Electrode, ε_Sep, 𝛿_Ele, 𝛿_Sep, l, β, τ, ω):
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
        i = self.I / (A * 10000)
        # molality: KOH molality [mol/kg]
        molality = ((self.W *10) / 56.1056) / 0.7
        # Vapour_Pressure: Vapour pressure of the KOH solution [bar]
        Vapour_Pressure = 10**(-0.01508 * molality - 0.0016788 * (molality**2) + 2.25887 * (10**-5) * (molality**3) 
                           + (1 - 0.0012062 * molality + 5.6024 * (10**-4) * (molality**2) - 7.8228 * (10**-6) * (molality**3)) 
                           * (35.4462 - (3343.93/self.T) - 10.9 * math.log10(self.T) + 0.0041645*self.T))
        # Theta: bubble coverage
        θ = (-97.25 + 182 * (self.T/self.T0) - 84 * ((self.T/self.T0)**2)) * (((self.I/(0.5*A)) / 300000)**0.3) * (self.P / (self.P-Vapour_Pressure))
        
        def Activation_Overpotential():
            # Exchange current densities [mA/cm2]
            i0_Cathode = r_Ca * math.exp((-G_Ca*1000/self.R)*((1/self.T) - (1/self.T0))) * i0_Ca_Ref
            i0_Anode = r_An * math.exp((-G_An*1000/self.R)*((1/self.T) - (1/self.T0))) * i0_An_Ref
            
            if self.I == 0:
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
            
            # Molar Concentration of KOH Solution
            m = (self.W * (183.1221 - 0.56845 * self.T + 984.5679 * math.exp(self.W/115.96277))) / 5610.5
            # Conductivity of the KOH Solution [S/cm]
            σ = -2.04 * m - 0.0028 * (m**2) + 0.005332 * m * self.T + 207.2 * (m/self.T) + 0.001043 * (m**3) - 0.0000003 * (m**2) * (self.T**2)
            ε = (2/3) * θ
            
            Ohm_Ele_free = (1/σ) * ((l/(0.5 * A * 10000)) + (l/(0.5 * A * 10000)))
            
            # Ohm resistance of bubble coverage electrolyte [S]
            if self.I == 0:
                Ohm_Ele_bubble = 0
                
            else:
                Ohm_Ele_bubble = Ohm_Ele_free * ((1/((1-ε)**1.5)) - 1)
             
            # Ohm resistance of electrolyte [S]
            Resistance_Electrolyte = Ohm_Ele_free + Ohm_Ele_bubble
            
            # Resistance of Separator [S]
            Resistance_Separator = (0.06 + 80 *0.5**(self.T/50)) / (10000 * A)
            
            ResistanceOverpotential = (Resistance_Electrode + Resistance_Separator + Resistance_Electrolyte) * self.I
            
            return ResistanceOverpotential
        
        CV = self.RCV + Activation_Overpotential() + Ohm_Resistance()
        
        return CV
    
    def CV_Jang(self, A, L, l):
        """
        2021년 Jang 의해 제안된 모델 
        Ref. : D. Jang et al., "Numerical Modeling and Analysis of the Effect of Pressure on the Performance of an Alkaline Water Electrolysis System", Applied Energy, 287, 116554 (2021).
        * 필요 parameters
        I: 전류[A], A: 총 전극의 표면적 [m2], w: KOH 전해질 농도 [%], L: 전극의 두께 [cm], l: 전극과 분리막 사이의 거리 [cm]
        """
        
        # Current density[A/cm2]
        j = self.I / (A * 10000)
    
        # Reference Parameters
        P_ref = 1
        T_ref = 293
        
        # Covering Coefficient
        # θ = 0.023 * (j**0.3) * ((self.T/T_ref) * (P_ref/self.P))**(2/3)
        θ = (-97.25 + 182 * (self.T/self.T0) - 84 * ((self.T/self.T0)**2)) * (((self.I/(0.5*A)) / 300000)**0.3) * 2
            
        def Activation_Overpotential():
            
            # Charge Transfer Coefficient
            α_An = 0.07835 + 0.001 * self.T
            α_Ca = 0.1175 + 0.00095 * self.T
            
            # Exchange Current Density of Anode and Cathode [A/cm2]
            j_An = 0.9 * (10**(-3)) * (((self.P/P_ref))**0.1) * math.exp((-42000/(self.R*self.T)) * (1-(self.T/T_ref)))
            j_Ca = 1.5 * (10**(-3)) * (((self.P/P_ref))**0.1) * math.exp((-23000/(self.R*self.T)) * (1-(self.T/T_ref)))
            
            # Activation Overpotential (Anode + Cathode) [V]               
            AO_An = ((self.R * self.T) / (self.F * α_An)) * (sp.arcsinh(j/j_An) - sp.arcsinh(1-θ))
            AO_Ca = ((self.R * self.T) / (self.F * α_Ca)) * (sp.arcsinh(j/j_Ca) - sp.arcsinh(1-θ))
                
            AO = AO_An + AO_Ca
                
            return AO
        
        def Ohm_Resistance():
            
            # Conductivity of the Ni Electrode [S/cm]
            σ_Ni = 60000000 - 279650 * self.T + 532 * self.T**2 - 0.38057 * self.T**3
            
            # Molar Concentration of KOH Solution
            m = (self.W * (183.1221 - 0.56845 * self.T + 984.5679 * math.exp(self.W/115.96277))) / 5610.5
            # Conductivity of the KOH Solution [S/cm]
            σ_KOH = -2.04 * m - 0.0028 * (m**2) + 0.005332 * m * self.T + 207.2 * (m/self.T) + 0.001043 * (m**3) - 0.0000003 * (m**2) * (self.T**2)
            
            # Resistance of Electrode
            Resistance_Ele = (1/σ_Ni) * (L/(A*0.5*10000)) + (1/σ_Ni) * (L/(A*0.5*10000))
            
            # Resistance of Electrolyte
            ε = (2/3) * θ
            Resistance_KOH = (1/(σ_KOH*((1-ε)**1.5))) * (l/(A*10000))
            
            # Resistance of Membrane
            S_mem = A * (1-θ)
            Resistance_Mem = (0.06 + 80 * math.exp(-self.T/50)) / (10000 * S_mem)
            
            ResOver = (Resistance_Ele + Resistance_KOH + Resistance_Mem) * self.I
            
            return ResOver
            
        CV = self.RCV + Activation_Overpotential() + Ohm_Resistance()
        
        return CV
    