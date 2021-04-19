# -*- coding: utf-8 -*-
"""
Created on Thu Apr  8 11:37:12 2021

@author: ybshi
"""

import math

""" 
알칼라인 수전해 시스템의 Electrochemical 모델링을 위한 온도, 압력에 따른 Open-circuit voltage (Reversible cell voltage)를 계산하는 코드
"""

class Reversible_Voltage(object):
    """ Open Circuit Voltage를 계산하는 4가지 모듈 """
    
    def __init__(self):
        
        super(Reversible_Voltage, self).__init__()
        self.T0 = 298
        self.R = 8.314 # 기체상수 - 단위 (J/mol K)
        self.F = 96485.34 # Faraday 상수 - 단위 (C/mol = J/mol V)
        self.z = 2  # 반응 1몰당 전자 2개 이동
        
    def RV_ConstantP(self, T):
        """
        압력이 1 atm로 일정할 때, 온도 T의 변화에 따른 Open-circuit voltage by LeRoy (1980)
        * 필요 parameters
        T: 온도 (K)
        """
        
        RevCV_ConP = 1.5184 - 1.5421 * (10**-3) * T + 9.523 * (10**-5) * T * math.log(T) + 9.84 * (10**-8) * (T**2)
        
        return RevCV_ConP
          
    def RV_LeRoy(self, T, P, w):
        """ 
        1980년 LeRoy에 의해 제안된 모델 
        Ref. : R. L. LeRoy and C. T. Bowen "The Thermodynamics of Aqueous Water Electrolysis", J. Electrochem. Soc., 127, pp. 1954-1962 (1980).
        온도, 압력에 따른 Open-circuit voltage를 계산하며, 압력에 영향을 주는 기체는 이상기체의 거동을 따른다. (1-100 atm)
        * 필요 parameters
        T: 온도 (K), P: 압력 (bar), w: 전해질의 질량% (wt%)
        * 변수 설명
        m: 전해질의 몰농도, Pw: x% KOH 용액의 증기압, VPw: 순수물의 증기압
        """
        
        if P == 1:
            RV = Reversible_Voltage().RV_ConstantP(T)
        
        else:
            P = 0.986923 * P # bar to atm
            m = (w * (183.1221 - 0.56845 * T + 984.5679 * math.exp(w/115.96277))) / 5610.5
            
            lnVPw = 37.04 - 6276/T - 3.416 * math.log(T)
            lnPw = 0.01621 - 0.1380 * m + 0.1933 * (m**0.5) + 1.024 * lnVPw
            
            VPw = math.exp(lnVPw)
            Pw = math.exp(lnPw)
            RV = Reversible_Voltage().RV_ConstantP(T) + ((self.R * T)/(self.z*self.F)) * math.log((P-Pw)**(1.5) * (VPw/Pw))
        
        return RV

    def RV_Onda(self, T, P):
        """ 
        2004년 Onda에 의해 제안된 모델 
        Ref. : K. Onda et al., "Prediction of Production Power for High-pressure Hydrogen by High-pressure Water Electrolysis", J. Power Sources, 132, pp. 64-70 (2004).
        온도, 압력에 따른 Open-circuit voltage를 계산하며, 높은 압력에서 적합하다. (1-700 atm)
        * 필요 parameters
        T: 온도 (K), P: 압력 (bar), w: 전해질의 질량% (wt%)
        * 변수 설명
        m: 전해질의 몰농도, Pw: 물의 압력, VPw: 물의 증기압
        """
            
        GibbsE = Gibbs_Energy(T, P).GibbsE()
        
        RV = GibbsE / (self.z*self.F)
        
        return RV

    def RV_Hammoudi(self, T, P, w):
        """ 
        2012년 Hammoudi에 의해 제안된 모델
        Nernst Equation기반 온도와 압력에 따른 Open-circuit voltage를 계산하는 함수 
        Ref. : M. Hammoudi et al., "New Multi-physics Approch for Modelling and Design of Alkaline Electrolyzers", Int. J. Hydrog. Energy, 37, pp. 13895-13913 (2012).
        온도, 압력에 따른 Open-circuit voltage를 계산하며, 첫번째 항은 온도의 영향, 뒤 세개의 항은 압력의 영향으로, 기존 LeRoy가 제안한 모델에 실제 기체 거동에 의한 영향을 추가하였다.
        * 필요 parameters
        T: 온도 (K), P: 압력 (bar), w: 전해질의 질량% (wt%)
        * 변수 설명
        m: 전해질의 몰농도, Pw: 물의 압력, VPw: 물의 증기압
        """
    
        if P == 1:
            RV = Reversible_Voltage().RV_ConstantP(T)
        
        else:
            P = 0.986923 * P # bar to atm
            m = (w * (183.1221 - 0.56845 * T + 984.5679 * math.exp(w/115.96277))) / 5610.5
            Pw = (T**-3.498) * math.exp(37.93 - (6426.32/T)) * math.exp(0.016214 - 0.13802 * m + 0.1933 * (m**0.5)) # atm
            VPw = (T**-3.4159) * math.exp(37.043 - (6275.7/T)) # atm
            
            RV = (Reversible_Voltage().RV_ConstantP(T) + (self.R*T)/(self.z*self.F) * math.log((P-Pw)**(1.5) * (VPw/Pw)) 
                   + (P-Pw) * (21.661 * (10**-6) - (5.471 * (10**-3)) / T) + ((P-Pw)**2) 
                   * ((-6.289 * (10**-6)) / T + (0.135 * (10**-3)) / (T**1.5) + (2.547 * (10**-3)) / (T**2) - 0.4825 / (T**3)))
        
        return RV
    
class Gibbs_Energy(Reversible_Voltage):
    """ 
    Hydrogen, Oxygen, Water의 깁스에너지를 계산하는 모듈
    * 변수 설명
    Entalpy0: 표준상태의 enthalpy, Entropy0: 표준상태의 entropy, Gib_T_1: 압력이 1 atm일 때 특정 T에서의 깁스에너지
    """
    
    def __init__(self, T, P):
        
        super(Gibbs_Energy, self).__init__()
        self.T = T
        self.P = P

    def Hydrogen(self):
        """ 특정 T, P에서 Hydrogen의 깁스에너지 계산 """
        
        a, b, c, e = 26.57, 3.77, 1.17, 0
        Entalpy0 = 0
        Entropy0 = 130.59
        
        Entalpy = Entalpy0 + a * (self.T - self.T0) + (b/2) * (10**-3) * ((self.T**2) - (self.T0**2)) - c * (10**5) * ((1/self.T) - (1/self.T0)) - (e/2) * (10**8) * ((1/self.T**2) - (1/self.T0**2))
        Entropy = Entropy0 + a * (math.log(self.T) - math.log(self.T0)) + b * (10**-3) * (self.T - self.T0) - (c/2) * (10**5) * ((1/self.T**2) - (1/self.T0**2)) - (e/3) * (10**8) * ((1/self.T**3) - (1/self.T0**3))
        Gib_T_1 = Entalpy - self.T * Entropy
        
        R = 82.06 # atm cm3 / mol K
            
        b1, b2 = 20.5, -1857
        c1, c2 = -351, 12760
        B = b1 + b2/self.T
        C = c1 + c2/(self.T**0.5)
        
        GibbsE_Hy = Gib_T_1 + (R * self.T * math.log(self.P) + B * self.P + ((C - (B**2)) * self.P**2) / (2 * R * self.T)) * 0.101325 # Conversion factor (atm cm3/mol to J/mol)
        
        return GibbsE_Hy
    
    def Oxygen(self):
        """ 특정 T, P에서 Oxygen의 깁스에너지 계산 """
        
        a, b, c, e = 34.35, 1.92, -18.45, 4.06
        Entalpy0 = 0
        Entropy0 = 205.29
        
        Entalpy = Entalpy0 + a * (self.T - self.T0) + (b/2) * (10**-3) * ((self.T**2) - (self.T0**2)) - c * (10**5) * ((1/self.T) - (1/self.T0)) - (e/2) * (10**8) * ((1/self.T**2) - (1/self.T0**2))
        Entropy = Entropy0 + a * (math.log(self.T) - math.log(self.T0)) + b * (10**-3) * (self.T - self.T0) - (c/2) * (10**5) * ((1/self.T**2) - (1/self.T0**2)) - (e/3) * (10**8) * ((1/self.T**3) - (1/self.T0**3))
        Gib_T_1 = Entalpy - self.T * Entropy
        
        R = 82.06 # atm cm3 / mol K
        
        b1, b2 = 42.6, -17400
        c1, c2 = -2604, 61457
        B = b1 + b2/self.T
        C = c1 + c2/(self.T**0.5)
        
        GibbsE_Ox = Gib_T_1 + (R * self.T * math.log(self.P) + B * self.P + ((C - (B**2)) * self.P**2) / (2 * R * self.T)) * 0.101325 # Conversion factor (atm cm3/mol to J/mol)
        
        return GibbsE_Ox
    
    def Water(self):
        """ 특정 T, P에서 water의 깁스에너지 계산 """
        
        a, b, c, e = 72.39, 9.38, 0, 0
        Entalpy0 = -285830
        Entropy0 = 69.96
        
        Entalpy = Entalpy0 + a * (self.T - self.T0) + (b/2) * (10**-3) * ((self.T**2) - (self.T0**2)) - c * (10**5) * ((1/self.T) - (1/self.T0)) - (e/2) * (10**8) * ((1/self.T**2) - (1/self.T0**2))
        Entropy = Entropy0 + a * (math.log(self.T) - math.log(self.T0)) + b * (10**-3) * (self.T - self.T0) - (c/2) * (10**5) * ((1/self.T**2) - (1/self.T0**2)) - (e/3) * (10**8) * ((1/self.T**3) - (1/self.T0**3))
        Gib_T_1 = Entalpy - self.T * Entropy
        
        GibbsE_Wa = Gib_T_1
        
        return GibbsE_Wa
    
    def GibbsE(self):
        """ H2O -> H2 + 1/2 O2 반응의 깁스 에너지를 계산 """
        
        GibbsE = Gibbs_Energy.Hydrogen(self) + 0.5 * Gibbs_Energy.Oxygen(self) - Gibbs_Energy.Water(self)
        
        return GibbsE