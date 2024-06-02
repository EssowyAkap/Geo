# -*- coding: utf-8 -*-
"""
Created on Fri May 10 22:04:21 2024

@author: szyme
"""

import numpy as np
from argparse import ArgumentParser

class Transformation:
    
    def __init__(self,elip):
        """
        definiuje parametry elipsoidy obrotowej.

       Parameters
       ----------
      
       elip [list] - lista parametrów elipsoidy obrotowej

       Returns
       -------
       *brak*

       """
        self.a = elip[0]
        self.e2 = elip[1]
        


    def N(self, f):
        """
        funkcja pomocnicza - obliczenie promienia krzywizny w I wertykale
        """
        N = self.a/((1 - self.e2 * np.sin(f)**2)**(1/2))
        return N

    def radtodegrees(self, radiany):
        '''
        daje stopnie dziesietne z podanych radianow
        '''
        stopnie_dziesietne = radiany*(180/np.pi)
        return stopnie_dziesietne

    def radians_dms(self, radiany):
        '''
        zamienia radiany na stopnie, minuty, sekundy
        '''
        stopnie_dziesietne = radiany*(180/np.pi)
        stp = np.trunc(stopnie_dziesietne)
        minuty = 60*(stopnie_dziesietne - stp)
        minuty2 = np.trunc(minuty)
        sek = 60*(minuty - minuty2)
        return f'{stp:.0f}\xb0{minuty2:.0f}\'{sek:.5f}\"'

    def dms(self, stopnie_dziesietne):
        '''
        funkcja zamieniajaca stopnie dziesietne na stopnie minuty i sekundy
        (jako lancuch znakow)
        '''
        stp = np.trunc(stopnie_dziesietne)
        minuty = 60*(stopnie_dziesietne-stp)
        minuty2 = np.trunc(minuty)
        sek = 60*(minuty - minuty2)
        return f'{stp:.0f}\xb0{minuty2:.0f}\'{sek:.5f}\"'
    
    def dms_rad(self, stopnie_minuty_sekundy):
        '''
        funkcja zamieniajaca wartosc stopnie minuty sekundy na radiany
        '''
        stp, minuty, sek = stopnie_minuty_sekundy
        radiany = (np.pi/180)*(stp + minuty/60 + sek/3600)
        return radiany
    
    def degrees_rad(self, stopnie_dziesietne):
        '''
        zamienia wielkosc stopni dziesietnych na wartosc w radianach
        '''
        
        radiany = (np.pi/180)*stopnie_dziesietne
        return radiany
    
    def XYZ2flh(self, X, Y, Z):
        '''
        przelicza wspolrzedne kartezjanskie(x,y,z) na wspolrzedne geodezyjne(f,lbda,h)
        przy wykorzystaniu algorytmu hirvonena
        -----
        Parameters:
            X, Y, Z [float] - wspolrzedne w ukladzie ortokartezjanskim
            
        -----
        Returns
        -------
        wyniki [list] - lista współrzędnych geodezyjnych po transformacji, podane w kolejnosci:
            f, l, h

        """
        '''
        wyniki = []
        for X, Y, Z in zip(X, Y, Z):
            p = np.sqrt(X**2+Y**2)
            f = np.arctan(Z/(p*(1-self.e2)))
            while True:
                N = self.N(f)
                h = (p/np.cos(f)) - N
                fp = f
                f = np.arctan((Z/p)/(1-((N*self.e2)/(N+h))))
                if abs(fp-f)<(0.000001/206265):
                    break
            N = self.N(f)
            h = p/np.cos(f)-N
            l = np.arctan(Y/X)
            wyniki.append([np.rad2deg(f), np.rad2deg(l), h])
            
        return wyniki

    def flh2XYZ(self, f, l, h):
        '''
        przelicza wspolrzedne geodezyjne(f, l, h) na wspolrzedne kartezjanskie(x,y,z)
        ---
        Parameters:
            f [float]- szerokosc geodezyjna
            l [float] - dlugosc geodezyjna
            h [float] - wysokosc elipsoidalna
            
            wielkosci f i l nalezy podac w radianach!!!!
        -----
        Returns
        -------
        wyniki [list] - lista wyników złożonych ze współrzędnych kartezjańskich ułożonych
            w kozlejnosci: X, Y, Z
            
        '''
        wyniki = []
        for f, l, h in zip(f, l, h):
            N = self.N(f)
            X = (N+h)*np.cos(f)*np.cos(l)
            Y = (N+h)*np.cos(f)*np.sin(l)
            Z = (N*(1-self.e2)+h)*np.sin(f)
            wyniki.append([X,Y,Z])
            
        return wyniki
    def Rneu(self, phi, lam):
        """
          Funkcja definiuje macierz obrotu układu NEU dla podanej elipsoidy, szerokosci i długosci geodezyjnej

         Parameters
         ----------
         phi : [float] szerokosć geodezyjna danego punktu
              jednostka: [rad]
         lam : [float] długosć geodezyjna danego punktu
              jednostka: [rad]

         Returns
         -------
         Rneu : [array] macierz obrotu układu NEU
                jednostka: brak

         """
        Rneu = np.array([[-np.sin(phi)*np.cos(lam), -np.sin(lam), np.cos(phi)*np.cos(lam)],
                         [-np.sin(phi)*np.sin(lam), np.cos(lam), np.cos(phi)*np.sin(lam)],
                         [np.cos(phi), 0, np.sin(phi)]])
        
        return Rneu
    
    def XYZ2neup(self, X, Y, Z, X0, Y0, Z0):
        
        """
        zamienia współrzędne geocentryczne XYZ dla danego punktu na współrzędne topocentryczne
        w układzie NEU dla punktu o srodku w X0, Y0, Z0
        Parameters
        ----------
        X [float] - współrzędna geocentryczna X
        Y [float] - współrzędna geocentryczna Y
        Z [float] - współrzędna geocentryczna Z
        X0 [float] - współrzędna geocentryczna X0 nadajnika
        Y0 [float] - współrzędna geocentryczna Y0 nadajnika
        Z0 [float] - współrzędna geocentryczna Y0 nadajnika

         wyniki
         -------
         wyniki [list] - lista współrzędnych punktow w układzie topocentrycznym NEU, podanych
         w kolejnsci: X, Y, Z

         """
        wyniki = []
        p = np.sqrt(X0**2+Y0**2)
        fi = np.arctan(Z0/(p*(1-self.e2)))
        while True:
            N = self.N(fi)
            h = (p/np.cos(fi)) - N
            fi_poprzednia = fi
            fi = np.arctan((Z0/p)/(1-((N*self.e2)/(N+h))))
            if abs(fi_poprzednia-fi)<(0.000001/206265):
                break 
        N = self.N(fi)
        h = p/np.cos(fi) - N
        lam = np.arctan(Y0/X0)
        
        R_neu = self.Rneu(fi, lam)
        for X, Y, Z in zip(X, Y, Z):
            X_sr = [X-X0, Y-Y0, Z-Z0] 
            X_rneu = R_neu.T@X_sr
            wyniki.append(X_rneu.T)
            
        return wyniki
        
    
    def PL2000(self,f,l,m=0.999923):
        """
        przelicza wspolrzedne geodezyjne f,l na współrzędne w układzie PL2000
         Parameters
         ----------
         f [float] - szerokosć geodezyjna danego punktu
            
         l [float] - długosć geodezyjna danego punktu
            
         m [float] - skala dla układu PL2000

         wyniki
         -------
         wyniki [list] - wspolrzedne w ukladzie PL2000 podane w kolejnosci: X, Y
         """
        wyniki = []
        for f, l in zip (f,l):
            l0 = 0
            strefa = 0
            if l >np.deg2rad(13.5) and l < np.deg2rad(16.5):
                strefa = 5
                l0 = np.deg2rad(15)
            elif l >np.deg2rad(16.5) and l < np.deg2rad(19.5):
                strefa = 6
                l0 = np.deg2rad(18)
            elif l >np.deg2rad(19.5) and l < np.deg2rad(22.5):
                strefa =7
                l0 = np.deg2rad(21)
            elif l >np.deg2rad(22.5) and l < np.deg2rad(25.5):
                strefa = 8
                l0 = np.deg2rad(24)
            
            b2 = self.a**2*(1-self.e2)
            ep2 = (self.a**2-b2)/b2
            tg = np.tan(f)
            dl = l - l0
            n2 = ep2*(np.cos(f)**2)
            N = self.N(f)
            A0 = 1- (self.e2/4)-(3*self.e2**2/64)-(5*self.e2**3/256)
            A2 = (3/8)*(self.e2+(self.e2**2/4)+(15*self.e2**3/128))
            A4 = (15/256)*(self.e2**2+((3*self.e2**3)/4))
            A6 = (35*self.e2**3)/3072
            sigma = self.a *(A0*f-A2*np.sin(2*f)+A4*np.sin(4*f)-A6*np.sin(6*f))
            
            xgk = sigma+(((dl**2/2)*N*np.sin(f)*np.cos(f))*(1+((dl**2/12)*(np.cos(f)**2)*(5-tg**2+9*n2+4*n2**2))+((dl**4/360)*(np.cos(f)**4)*(61-58*tg**2+tg**4+270*n2-330*n2*tg**2))))
            ygk = (dl*N*np.cos(f))*(1+((dl**2/6)*(np.cos(f)**2)*(1-tg**2+n2))+(((dl**4/120)*(np.cos(f)**4))*(5-(18*tg**2)+tg**4+(14*n2)-(58*n2*tg**2))))
            x2000 = xgk*m
            y2000 = ygk*m + (strefa *1000000) +500000
            wyniki.append([x2000, y2000])
        return wyniki
    
    
    def PL1992(self,f,l,m=0.9993):
        """
        przelicza wspolrzedne geodezyjne f,l na współrzędne w układzie PL1992

         Parameters
         ----------
         f [float] - szerokosć geodezyjna danego punktu
        
        l [float] - długosć geodezyjna danego punktu
        
        m [float] - skala dla układu PL1992


         wyniki
         -------
         wyniki [list] - wspolrzedne w ukladzie PL1992 podane w kolejnosci: X, Y

         """
        
        wyniki = []
        l0 = np.deg2rad(19)
        for f, l in zip (f,l):
                b2 = self.a**2*(1-self.e2)
                ep2 = (self.a**2-b2)/b2
                tg = np.tan(f)
                dl = l - l0
                n2 = ep2*(np.cos(f)**2)
                         
                A0 = 1- (self.e2/4)-(3*self.e2**2/64)-(5*self.e2**3/256)
                A2 = (3/8)*(self.e2+(self.e2**2/4)+(15*self.e2**3/128))
                A4 = (15/256)*(self.e2**2+((3*self.e2**3)/4))
                A6 = (35*self.e2**3)/3072
                N = self.N(f)
                sigma = self.a *(A0*f-A2*np.sin(2*f)+A4*np.sin(4*f)-A6*np.sin(6*f))
                xgk = sigma+( ((dl**2/2)*N*np.sin(f)*np.cos(f))*(1+((dl**2/12)*(np.cos(f)**2)*(5 - tg**2 + 9*n2 + 4*n2**2))+((dl**4/360)*(np.cos(f)**4)*(61 - 58*tg**2 + tg**4 + 270*n2 - 330*n2*tg**2))))
                ygk = (dl*N*np.cos(f))*(1+((dl**2/6)*(np.cos(f)**2)*(1 - tg**2 + n2))+(((dl**4/120)*(np.cos(f)**4)) * (5 - (18*tg**2) + tg**4 + (14 * n2) - (58*n2*tg**2))))
                x92 = (xgk*m)-5300000
                y92 = (ygk*m)+500000
                wyniki.append([x92, y92])
        return wyniki
            
    
        

    
    def zwroty(self,wejsciowy_plik, rodzaj_transformacji):
        """
         wykorzystuje dane z pliku wejsciowego do przeprowadzenia danego rodzaju transformacji,
        tworzy plik z wynikami transformacji: 'wyniki(rodzaj_transformacji)_(args.el).txt'

        Parameters
        ----------
        wejsciowy_plik [.txt file] - plik z danymi
        rodzaj_transformacji [str] - tekst; zaleznie od rodzaju transformacji jaki uzytkownik
        chce uzyskac stosuje sie: XYZ2flh, flh2XYZ, XYZ2neu, PL1992, PL2000

        Returns
        -------
        brak
         """
        dane = np.genfromtxt(wejsciowy_plik,delimiter = " ")
        if rodzaj_transformacji == 'XYZ2BLH':
            wyniki = self.XYZ2flh(dane[:,0], dane[:,1], dane[:,2])
            np.savetxt(f"obliczone_{rodzaj_transformacji}_{args.elp}.txt", wyniki, delimiter=' ', fmt='%0.10f %0.10f %0.3f')
        
        elif rodzaj_transformacji == 'BLH2XYZ':
            wyniki = self.flh2XYZ(np.deg2rad((dane[:,0])), np.deg2rad(dane[:,1]), dane[:,2])
            np.savetxt(f"obliczone_{rodzaj_transformacji}_{args.elp}.txt", wyniki, delimiter=' ', fmt='%0.3f %0.3f %0.3f')

        elif rodzaj_transformacji == 'XYZ2NEUP':
            wyniki = self.XYZ2neup(dane[1:,0], dane[1:,1], dane[1:,2], dane[0,0], dane[0,1], dane[0,2])
            np.savetxt(f"obliczone_{rodzaj_transformacji}._{args.elp}.txt", wyniki, delimiter=' ', fmt='%0.3f %0.3f %0.3f')

        elif rodzaj_transformacji == 'PL2000':
            wyniki = self.PL2000(np.deg2rad(dane[:,0]), np.deg2rad(dane[:,1]))
            np.savetxt(f"obliczone_{rodzaj_transformacji}_{args.elp}.txt", wyniki, delimiter=' ', fmt='%0.3f %0.3f')
        
        elif rodzaj_transformacji == 'PL1992':
            wyniki = self.PL1992(np.deg2rad(dane[:,0]), np.deg2rad(dane[:,1]))
            np.savetxt(f"obliczone_{rodzaj_transformacji}_{args.elp}.txt", wyniki, delimiter=' ', fmt='%0.3f %0.3f')
        
    
if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('-plik', type=str, help='przyjmuje dany plik')
    parser.add_argument('-rt', type=str, help='rodzaj transformacji')
    parser.add_argument('-elp', type=str, help='przyjmuje dana elipsoide')
    args = parser.parse_args()

    elipsoidy = {'GRS80':[6378137.000, 0.00669438002290],'WGS84':[6378137.000, 0.00669438002290],  'KRASOWSKI':[6378245.000, 0.00669342162296]}
    transformacje = {'BLH2XYZ': 'BLH2XYZ', 'XYZ2NEUP':'XYZ2NEUP', 'PL2000':'PL2000','PL1992':'PL1992', 'XYZ2BLH': 'XYZ2BLH'}
    
    wybor = "DALEJ"
    try:
        while wybor =="DALEJ":
            if args.elp==None:
                args.elp = input(str('Elipsoida:'))
            if args.plik==None:
                args.plik = input(str('Sciezka do pliku .txt z danymi:'))
            if args.rt==None:
                args.rt = input(str('Nazwa żądanej transformacji:'))
                    
            obiekt = Transformation(elipsoidy[args.elp.upper()])
            dane = obiekt.zwroty(args.plik, transformacje[args.rt.upper()])
                
            print('Utworzono plik z wynikami.')               
            wybor = input(str("Wcisnij ENTER, zeby zakonczyc; jesli chcesz dokonac kolejnej transformacji wpisz DALEJ")).upper()
            args.rt= None
            args.elp = None
            args.plik= None



    except IndexError:
        print('W pliku podano niewłasciwe dane!')
    except ValueError:
        print('W pliku podano niewłasciwe dane!')
    except FileNotFoundError:
        print('Nie znaleziono podanego pliku!')
    except KeyError:
        print('Elipsoida lub rodzaj transformacji zostały niewłasciwie podane!')
    finally:
        print('Koniec na dzis!')      
