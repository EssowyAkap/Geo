# -*- coding: utf-8 -*-
"""
Created on Fri May 10 22:04:21 2024

@author: szyme
"""

import numpy as np
from argparse import ArgumentParser as args

class Transformation:
    
    def __init__(self,elip):
        """
    definiuje parametry elipsoidy obrotowej.

       Parameters
       ----------
      
       elip [list] - lista parametrów elipsoidy obrotowej

       Returns
       -------
       None.

       """
        self.a = elip[0]
        self.e2 = elip[1]
        


    def N(self, f):
        """
        funkcja pomocnicza - obliczenie promienia krzywizny w I wertykale
        """
        N = self.a/((1 - self.e2 * np.sin(f)**2)**(1/2))
        return(N)

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
        przelicza wspolrzedne kartezjanskie(x,y,z) na wspolrzedne geodezyjne(f,lbda,h)
        przy wykorzystaniu algorytmu hirvonena. Zwraca je z dokladnoscia  okolo 1 cm 
        
        -----
        Parameters:
            X, Y, Z [float] - wspolrzedne w ukladzie ortokartezjanskim
            
        -----
        Returns
        -------
        returns [list] - lista współrzędnych geodezyjnych po transformacji, podane w kolejnosci:
            f, l, h

        """
        '''
        returns = []
        for X, Y, Z in zip(X, Y, Z):
            p = np.sqrt(X**2+Y**2)
            f = np.arctan(Z/(p*(1-self.e2)))
            while True:
                N = self.Npu(f)
                h = (p/np.cos(f)) - N
                fp = f
                f = np.arctan((Z/p)/(1-((N*self.e2)/(N+h))))
                if abs(fp-f)<(0.000001/206265):
                    break
            N = self.N(f)
            h = p/np.cos(f)-N
            l = np.arctan(Y/X)
            returns.append([np.rad2deg(f), np.rad2deg(l), h])
            
        return(returns)

    def flh2XYZ(self, f, l, h):
        '''
        przelicza wspolrzedne geodezyjne(f, l, h) na wspolrzedne kartezjanskie(x,y,z)
        ---
        Parameters:
            f - szerokosc geodezyjna
            l - dlugosc geodezyjna
            h [float] - wysokosc elipsoidalna
            
        -----
        Returns
        -------
        returns [list] - lista wyników złożonych ze współrzędnych kartezjańskich ułożonych
            w kozlejnosci: X, Y, Z
            
        '''
        returns = []
        for f, l, h in zip(f, l, h):
            N = self.Npu(f)
            X = (N+h)*np.cos(f)*np.cos(l)
            Y = (N+h)*np.cos(f)*np.sin(l)
            Z = (N*(1-self.e2)+h)*np.sin(f)   
            returns.append([X,Y,Z])
            
        return(returns)
    
        
    
    def PL2000(self,f,l,m=0.999923):
        """
        przelicza wspolrzedne geodezyjne f,l na współrzędne w układzie PL2000
         Parameters
         ----------
         f [float] - szerokosć geodezyjna danego punktu
            
         l [float] - długosć geodezyjna danego punktu
            
         m [float] - skala dla układu PL2000

         Returns
         -------
         returns [list] - wspolrzedne w ukladzie PL2000 podane w kolejnosci: X, Y
         """
        returns = []
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
            returns.append([x2000, y2000])
        return(returns)
    
    def XYZ2neu(self, X, Y, Z, X0, Y0, Z0):
        
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

         Returns
         -------
         returns [list] - lista współrzędnych punktow w układzie topocentrycznym NEU, podanych
         w kolejnsci: X, Y, Z

         """
        returns = []
        p = np.sqrt(X0**2+Y0**2)
        f = np.arctan(Z0/(p*(1-self.e2)))
        while True:
             N = self.Npu(f)
             h = (p/np.cos(f)) - N
             fp = f
             f = np.arctan((Z0/p)/(1-((N*self.e2)/(N+h))))
             if abs(fp-f)<(0.000001/206265):
                 break
             N = self.Npu(f)
             h = p/np.cos(f)-N
             l = np.arctan(Y0/X0)
             Rneu = np.array([[-np.sin(f)*np.cos(l), -np.sin(l), np.cos(f)*np.cos(l)],
                              [-np.sin(f)*np.sin(l), np.cos(l), np.cos(f)*np.sin(l)],
                              [np.cos(f), 0, np.sin(f)]])
             
             for X, Y, Z in zip(X, Y, Z):
                 Xs = [X-X0, Y-Y0, Z-Z0]
                 Xr = Rneu.T@Xs
                 returns.append([Xr.T])
        
        
        return(returns)
    
    def PL1992(self,f,l,m=0.9993):
        """
        przelicza wspolrzedne geodezyjne f,l na współrzędne w układzie PL1992

         Parameters
         ----------
         f [float] - szerokosć geodezyjna danego punktu
        
        l [float] - długosć geodezyjna danego punktu
        
        m [float] - skala dla układu PL1992


         Returns
         -------
         returns [list] - wspolrzedne w ukladzie PL1992 podane w kolejnosci: X, Y

         """
        
        returns = []
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
                returns.append([x92, y92])
        return(returns)
            
    
        

    
    def wynik(self, wejsciowy_plik, rodzaj_transformacji):
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
         *brak*

         """
        dane = np.genfromtxt(wejsciowy_plik,delimiter = " ")
        if rodzaj_transformacji == 'XYZ2flh':
            returns = self.XYZ2flh(dane[:,0], dane[:,1], dane[:,2])
            np.savetxt(f"wyniki{rodzaj_transformacji}_{args.el}.txt", returns, delimiter=' ', fmt='%0.10f %0.10f %0.3f')
        
        elif rodzaj_transformacji == 'flh2XYZ':
            returns = self.flh2XYZ(np.deg2rad((dane[:,0])), np.deg2rad(dane[:,1]), dane[:,2])
            np.savetxt(f"wyniki{rodzaj_transformacji}_{args.el}.txt", returns, delimiter=' ', fmt='%0.3f %0.3f %0.3f')

        elif rodzaj_transformacji == 'PL2000':
            returns = self.PL2000(np.deg2rad(dane[:,0]), np.deg2rad(dane[:,1]))
            np.savetxt(f"wyniki{rodzaj_transformacji}_{args.el}.txt", returns, delimiter=' ', fmt='%0.3f %0.3f')
        
        elif rodzaj_transformacji == 'PL1992':
            returns = self.PL1992(np.deg2rad(dane[:,0]), np.deg2rad(dane[:,1]))
            np.savetxt(f"wyniki{rodzaj_transformacji}_{args.el}.txt", returns, delimiter=' ', fmt='%0.3f %0.3f')
        
        elif rodzaj_transformacji == 'XYZ2neu':
            returns = self.XYZ2neu(dane[1:,0], dane[1:,1], dane[1:,2], dane[0,0], dane[0,1], dane[0,2])
            np.savetxt(f"wyniki{rodzaj_transformacji}._{args.el}.txt", returns, delimiter=' ', fmt='%0.3f %0.3f %0.3f')
    
if __name__ == '__main__':
    parser = args()
    parser.add_argument('-t', type=str, help='Przyjmuje nazwe wybranej transformacji. Dostepne: XYZ2flh, flh2XYZ, XYZ2neu, PL1992, PL2000 ')
    parser.add_argument('-el', type=str, help='Przyjmuje daną elipsoidę. Dostepne: WGS84, GRS80 lub KRASOWSKI')
    parser.add_argument('-p', type=str, help='przyjmuje dany plik')
    args = parser.parse_args()

    elipsoidy = {'WGS84':[6378137.000, 0.00669438002290], 'GRS80':[6378137.000, 0.00669438002290], 'KRASOWSKI':[6378245.000, 0.00669342162296]}
    transformacje = {'XYZ2flh': 'XYZ2flh','flh2XYZ': 'flh2XYZ','PL2000':'PL2000','PL1992':'PL1992', 'XYZ2neu':'XYZ2neu'}
 
    wybor = "TAK"
    try:
        while wybor =="TAK":
            if args.el==None:
                args.el = input(str('Nazwa elipsoidy: '))
            if args.p==None:
                args.p = input(str('Podaj sciezke dostępu pliku z danymi: '))
            if args.t==None:
                args.t = input(str('Podaj rodzaj żądanej transformacji: '))
                    
            obiekt = Transformation(elipsoidy[args.el.upper()])
            dane = obiekt.wynik(args.p, transformacje[args.t.upper()])
                
            print('utworzono plik z wynikami')
                
            wybor = input(str("Kliknij ENTER aby zakonczyc program. W celu wykonania kolejnej transformacji wpisz ,,TAK'' ")).upper()
            args.el = None
            args.p= None
            args.t= None

    except FileNotFoundError:
        print('Podany plik nie istnieje.')
    except KeyError:
        print('Niepoprawnie podana elipsoida lub transformacja.')
    except IndexError:
        print('Sprawdz format danych pliku.')
    except ValueError:
        print('Sprawdz format danych pliku.')
    finally:
        print('Koniec programu')
