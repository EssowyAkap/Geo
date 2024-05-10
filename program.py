# -*- coding: utf-8 -*-
"""
Created on Fri May 10 22:04:21 2024

@author: szyme
"""

import numpy as np
import sys 

class Transformacje:
    
    def __init__(self, model: str = 'wgs84'):
        """
        ----------
        model : str, model elipsoidy, domyslny model to WGS84,
            dostepne modele:
                + WGS84 = 'wgs84'
                + GRS80 = 'grs80'
                + ELipsoida Krasowskiego = 'elipsoida_krasowskiego'
        ----------
        Parametry elipsoid:
            a - duza polos
            b - mala polos
            flat - splaszczenie
            e2 - mimosrod elipsoidy podniesiony do kwadratu
        """
        if model == "wgs84":
            self.a = 6378137.0
            self.b = 6356752.31424518
        elif model == "grs80":
            self.a = 6378137.0
            self.b = 6356752.31414036
        elif model == "elipsoida_krasowskiego":
            self.a = 6378245.0
            self.b = 6356863.019
        else:
            raise NotImplementedError(f"model {model} nie jest obslugiwany")
            
        self.flat = (self.a - self.b) / self.a
        self.e = np.sqrt(2 * self.flat - self.flat ** 2)
        self.e2 = (2 * self.flat - self.flat ** 2) 

    def obliczenie_N(self, f):
        """
        funkcja pomocnicza - obliczenie promienia krzywizny w I wertykale
        """
        N = self.a/((1 - self.e2 * np.sin(f)**2)**(1/2))
        return N
    