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

    def N(self, f):
        """
        funkcja pomocnicza - obliczenie promienia krzywizny w I wertykale
        """
        N = self.a/((1 - self.e2 * np.sin(f)**2)**(1/2))
        return(N)
    def d_poludnika(self, f):
        '''
        funkcja wykorzystujaca szerokosc geodezyjna punktu w celu obliczenia
        dlugosci luku danego poludnika
        '''
        A0 = 1 - (self.e2/4) - ((3*(self.e2**2))/64) - ((5*(self.e2**3))/256)
        A2 = (3/8) * (self.e2 + ((self.e2**2)/4) + ((15*(self.e2**3))/128))
        A4 = (15/256) * (self.e2**2 + ((3*(self.e2**3))/4))
        A6 = (35*(self.e2**3))/3072
        sigma = self.a * (A0*f - A2*np.sin(2*f) + A4*np.sin(4*f) - A6*np.sin(6*f))
        
        return(sigma)
    
    def radtodegrees(self, radiany):
        '''
        daje stopnie dziesietne z podanych radianow
        '''
        stopnie_dziesietne = radiany*(180/np.pi)
        return(stopnie_dziesietne)

    def radians_dms(self, radiany):
        '''
        zamienia radiany na stopnie, minuty, sekundy 
        '''
        stopnie_dziesietne = radiany*(180/np.pi)
        stp = np.trunc(stopnie_dziesietne)
        minuty = 60*(stopnie_dziesietne - stp)
        minuty2 = np.trunc(minuty)
        sek = 60*(minuty - minuty2)
        return(f'{stp:.0f}\xb0{minuty2:.0f}\'{sek:.5f}\"')

    def dms(self, stopnie_dziesietne):
        '''
        funkcja zamieniajaca stopnie dziesietne na stopnie minuty i sekundy
        (jako lancuch znakow)
        '''
        stp = np.trunc(stopnie_dziesietne)
        minuty = 60*(stopnie_dziesietne-stp)
        minuty2 = np.trunc(minuty)
        sek = 60*(minuty - minuty2)
        return(f'{stp:.0f}\xb0{minuty2:.0f}\'{sek:.5f}\"')
    
    def dms_rad(self, stopnie_minuty_sekundy):
        '''
        funkcja zamieniajaca wartosc stopnie minuty sekundy na radiany
        '''
        stp, minuty, sek = stopnie_minuty_sekundy
        radiany = (np.pi/180)*(stp + minuty/60 + sek/3600)
        return(radiany)
    
    def degrees_rad(self, stopnie_dziesietne):
        '''
        zamienia wielkosc stopni dziesietnych na wartosc w radianach
        '''
        radiany = (np.pi/180)*stopnie_dziesietne 
        return(radiany)


    def XYZ2flh(self, X, Y, Z, output = 'dec'):
        '''
        przelicza wspolrzedne kartezjanskie(x,y,z) na wspolrzedne geodezyjne(phi,lambda,h)
        przy wykorzystaniu algorytmu hirvonena. Zwraca je z dokladnoscia  okolo 1 cm 
        
        -----
        parametry:
            X, Y, Z [float] - wspolrzedne w ukladzie ortokartezjanskim
            
        -----
        returns:
            f [float] - szerokosc geodezyjna
            l [float] - dlugosc geodezyjna
            h [float] - wysokosc elipsoidalna [m]
            
        w zaleznosci od parametru output szerokosc i dlugosc geodezyjna
        moga byc podane w dwoch formatach:
                dec - stopnie dziesietne (wartosc domyslna)
                dms - stopnie, minuty, sekundy
                calc - zwrocone wartosci w radianach, np. jesli nie potrzeba
                    dokladnie tych wartosci, tylko do dalszych obliczen
                
        '''
        p = np.sqrt(X**2 + Y**2)
        f = np.atan(Z/(p*(1-self.e2)))
        while True:
            f_p = f
            # N = self.a/((1 - self.e2 * m.sin(f)**2)**(1/2))
            N = self.N(f_p)
            h = p/np.cos(f_p) - N
            f = np.atan(Z/p * ((N*(1-self.e2)+h)/(N+h))**(-1))
            if (0.000001/206265) > abs(f_p - f):
                break
        N = self.N(f)
        l = np.atan2(Y,X)
        h = p/np.cos(f) - N
                
        if output == "dziesietne":
            return(self.radtodegrees(f), self.radtodegrees(l), h)
             
        elif output == "dms":
            return(self.radians_dms(f), self.radians_dms(l), h)
             
        elif output == "calc":
             return(f, l, h) 
         
         
    def flh2XYZ(self, f, l, h, inp = 'dec'):
        '''
        przelicza wspolrzedne geodezyjne(phi, lambda, h) na wspolrzedne kartezjanskie(x,y,z)
        ---
        parametry:
            f - szerokosc geodezyjna
            lam - dlugosc geodezyjna
            ^ typ w zaleznosci od inp, mozliwosci:
                    inp = 'dec' -> [float] - stopnie dziesietne (wartosc domyslna)
                    inp = 'dms' -> [tuple] - krotka: (stopnie, minuty, sekundy)
            h [float] - wysokosc 
            
        ---
        returns: 
            X, Y, Z - wspolrzedne ortokartezjanskie 
            
        '''
        if inp == 'dec':
            f2 = self.degrees_rad(f)
            l2 = self.degrees_rad(l)
            
        elif inp == 'dms':
            f2 = self.dms_rad(f)
            l2 = self.dms_rad(l)
            
        N = self.N(f)
        X = (N + h) * np.cos(f2) *np.cos(l2)
        Y = (N + h) * np.cos(f2) * np.sin(l2)
        Z = (N*(1-self.e2) + h) * np.sin(f2)        
        return(X, Y, Z)
    