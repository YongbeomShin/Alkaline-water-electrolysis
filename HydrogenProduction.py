# -*- coding: utf-8 -*-
"""
Created on Tue Apr 27 17:41:56 2021

@author: ybshi
"""

class Hydrogen_Production():
    
    def __init__(self, T, P, I, NumCell):
        """ I: 전류 [A], A: 전극면적 [m2] """
        
        self.I = I
        self.T = T
        self.P = P
        self.F = 96485.34 # Faraday 상수 - 단위 (C/mol = J/mol V)
        self.z = 2  # 반응 1몰당 전자 2개 이동
        self.R = 8.205 * 10**(-5) # atm m3 /mol K
        self.NumCell = NumCell
     
    def HP_Ulleberg(self, A, f, f2):
        
        I = self.I * 1000 # A to mA
        A = A * 10000 # m2 to cm2
        
        η = (((I/A)**2) / (f + (I/A)**2)) * f2
        n_Hydrogen = η * (self.I/(self.F*self.z))
        
        return n_Hydrogen
        
    def HP_Jang(self):
        
        I = self.I * 1000 # A to mA
        A = 0.25 * 10000 # m2 to cm2
        
        f = 2.5 * self.T - 632.5
        f2 = -0.00075 * self.T + 1.22475
        
        η = (((I/A)**2) / (f + (I/A)**2)) * f2
        
        Nc = self.NumCell # Number of cell
        n_Hydrogen = ((η*self.I)/(self.z*self.F)) * ((self.R*298)/1) * 3600 * Nc
        
        return n_Hydrogen
        