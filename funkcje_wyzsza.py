# -*- coding: utf-8 -*-
"""
Created on Wed Oct 18 10:21:33 2023

@author: piotr
"""


from math import *
import numpy as np



a_grs = 6378137
e2_grs = 0.00669438002290


def Np(f, a, e2):
    N = a/np.sqrt(1 - e2*np.sin(f)**2)
    return(N)

def hirvonen(X, Y, Z, a, e2):
    p = np.sqrt(X**2 + Y**2)
    f = np.arctan(Z / (p * (1 - e2)))
    while True:
        N = Np(f, a, e2)
        h = (p / np.cos(f)) - N
        fs = f
        f = np.arctan(Z / (p * (1 - e2 * (N / (N + h)))))
        if np.abs(fs - f) < (0.000001 / 206265):
            break
    l = np.arctan2(Y, X)    
    return(f, l, h)


def dms(x):
    sig = ' '
    if x<0:
        sig = '-'
        x = abs(x)
    x = x * 180 / pi
    d = int(x)
    m = int(60 * (x - d))
    s = (x - d - m/60) * 3600
    print(f'{sig}{d:3d}{chr(176)}{abs(m):2d}\'{abs(s):7.5f}\"')
    
def XYZ(fi, lm, h, a, e2):
    N = Np(fi, a, e2)
    
    x = (N + h) * np.cos(fi) * np.cos(lm)
    y = (N + h) * np.cos(fi) * np.sin(lm)
    z = (N * (1 - e2) + h) * np.sin(fi)
    
    return(x, y, z)

    
def Mp(f, a, e2):
    M = (a * (1 - e2)) / np.sqrt((1 - (e2 * (np.sin(f))**2))**3)
    return(M)

def claiuraut(f, A, a, e2):
    N = Np(f, a, e2)
    C = N * np.cos(f) * np.sin(A)
    return(C)


def R(f, a, e2):
    R = np.sqrt((Np(f, a, e2) * (Mp(f, a, e2))))
    return(R)







def kivioj(f, l, A, s, a, e2):
    n = int(s/1000)
    ds = s/n
    for i in range(n):
        M = Mp(f, a, e2)
        N = Np(f, a, e2)
        df = ds * np.cos(A) / M 
        dA = (np.sin(A) * np.tan(f) * ds) / N
        fm = f + df / 2
        Am = A + dA / 2 
        Mm = Mp(fm, a, e2)
        Nm = Np(fm, a, e2)
        df = (ds * np.cos(Am)) / Mm
        dl = (ds * np.sin(Am)) / (Nm * np.cos(fm))
        dA = (np.sin(Am) * np.tan(fm) * ds) / Nm
        f = f + df
        l = l + dl
        A = A + dA
    A = A + pi
    if A > 2 * pi:
        A = A - 2 * pi
        
    return(f, l, A)


def degrees(stonie, minuty, sekundy):
    stopnie_dziesietne = stopnie + (minuty / 60) + (sekundy / 3600)
    return(stopnie_dziesietne)

def radians_dms(stopnie, minuty, sekundy):
    stopnie_dziesietne = stopnie + (minuty / 60) + (sekundy / 3600)
    radians = stopnie_dziesietne * pi / 180
    return(radians)


def vincenty(fa, la, fb, lb, a_grs, e2_grs):
    b = a_grs * np.sqrt(1 - e2_grs)
    fl = 1 - b/a_grs
    dL = lb - la
    Ua = np.arctan((1-fl) * np.tan(fa))
    Ub = np.arctan((1-fl) * np.tan(fb))
    L = dL
    while True:
            sin_sig = np.sqrt((np.cos(Ub) * np.sin(L))**2 + (np.cos(Ua) * np.sin(Ub) - np.sin(Ua) * np.cos(Ub) *
            np.cos(L))**2)
            cos_sig = (np.sin(Ua) * np.sin(Ub)) + (np.cos(Ua) * np.cos(Ub) * np.cos(L))
            sigma = np.arctan2(sin_sig, cos_sig)
            sin_alfa = (np.cos(Ua) * np.cos(Ub) * np.sin(L)) / sin_sig
            cos2_alfa = 1 - (sin_alfa)**2
            cos_2sig_m = cos_sig - ((2 * np.sin(Ua) * np.sin(Ub)) / cos2_alfa)
            C = fl/16 * cos2_alfa * (4 + fl * (4 - 3 * cos2_alfa))
            Lstare = L
            L = dL + (1 - C) * fl * sin_alfa * (sigma + C * sin_sig * (cos_2sig_m + C * cos_sig * (-1 + 2 * (cos_2sig_m)**2)))
            if np.abs(Lstare - L) < (0.000001/206265):
                break
    u2 = ((a_grs**2 - b**2) / b**2) * cos2_alfa
    A = 1 + (u2/16384) * (4096 + u2 * (-768 + u2 * (320 - 175 * u2)))
    B = (u2/1024) * (256 + u2 * (-128 + u2 * (74 - 47 * u2)))
    dsigma = B * sin_sig * (cos_2sig_m + (1/4 * B) * (cos_sig * (-1 + 2 * (cos_2sig_m)**2) - (1/6 * B * cos_2sig_m) * (-3
    + 4 * (sin_sig)**2) * (-3 + 4 * (cos_2sig_m)**2)))
    s = b * A *(sigma - dsigma)
    Aab = np.arctan2((np.cos(Ub) * np.sin(L)),((np.cos(Ua) * np.sin(Ub)) - (np.sin(Ua) * np.cos(Ub) *
    np.cos(L)))) 
    Aba = np.arctan2((np.cos(Ua) * np.sin(L)),(-np.sin(Ua) * np.cos(Ub)) + (np.cos(Ua) * np.sin(Ub) *
    np.cos(L))) + np.pi
    return(s, Aab, Aba)



def saz2neu(s, alfa, z):
    dx = np.array([s*np.sin(z)*np.cos(alfa),
                   s*np.sin(z)*np.sin(alfa),
                   s*np.cos(z)])
    return(dx)



def Rneu(f,l):
    R = np.array([[-np.sin(f)*np.cos(l),-np.sin(l),np.cos(f)*np.cos(l)],
                  [-np.sin(f)*np.sin(l),np.cos(l),np.cos(f)*np.sin(l)],
                  [np.cos(f),0,np.sin(f)]])
    return(R)


def neu2XYZ(dx, f, l):
    R = Rneu(f,l)
    dX = R @ dx
    return(dX)


def neu2saz(dx):
    s = np.sqrt(dx @ dx)
    alfa = np.arctan2(dx[1],dx[0])   
    z = np.arccos(dx[2]/s) 
    return(s,alfa, z)
 
    
def XYZ2neu(dX, f, l):
    R = Rneu(f, l)
    dx = R.T @ dX
    return(dx)


def b2(a, e2):
    b2 = a**2 * (1 - e2)
    return(b2)


def ep2(a, e2):
    ep2 = (a**2 - b2(a, e2)) / b2(a, e2)
    return(ep2)



def n_2(a, e2, f):
    n2 = ep2(a, e2) * (np.cos(f))**2
    return(n2)


def dla(l, l0):
    dl = l - l0
    return(dl)


def tp(f):
    t = np.tan(f)
    return(t)




def sigma1(f, a, e2):  # długosć południka
    A0 = 1 - (e2 / 4) - (3 * e2**2 / 64) - (5 * e2**3 / 256)
    A2 = (3 / 8) * (e2 + e2**2 / 4 + 15 * e2**3 / 128)
    A4 = (15 / 256) * (e2**2 + 3 * e2**3 / 4)
    A6 = 35 * e2**3 / 3072
    si = a * (A0 * f - A2 * np.sin(2 * f) + A4 * np.sin(4 * f) - A6 * np.sin(6 * f))
    return(si)


def fl2gk(fa,la,l0,a,e2):
    b2=a**2*(1-e2)
    e_2=(a**2-b2)/b2
    dl=la-l0
    t=tan(fa)
    n2=e_2*cos(fa)**2
    N=Np(fa,a,e2)
    sigma=sigma1(fa,a, e2)
    x=sigma+(dl**2/2)*N*sin(fa)*cos(fa)*(1+(dl**2/12)*(cos(fa))**2*(5-t**2+9*n2+4*n2**2)+(dl**4/360)*(cos(fa))**4*(61-58*t**2+t**4+270*n2-330*n2*t**2))
    y=dl*N*cos(fa)*(1+(dl**2/6)*(cos(fa))**2*(1-t**2+n2)+(dl**4/120)*(cos(fa))**4*(5-18*t**2+t**4+14*n2-58*n2*t**2))
    return(x,y)


# do przeliczenia odwrotnego


def f1(xgk, a, e2):
    A0 = 1 - (e2 / 4) - (3 * e2**2 / 64) - (5 * e2**3 / 256)
    f1 = xgk / (a * A0)
    while True:
        fs = f1
        s = sigma1(f1, a, e2)
        f1 = f1 + (xgk - s) / (a * A0)
        if abs(f1 - fs) < (0.000001 / 206265):
            break
    return(f1)


def Np1(f1, a, e2):
    N1 = Np(f1, a, e2)
    return(N1)


def Mp1(f1, a, e2):
    M1 = Mp(f1, a, e2)
    return(M1)


def t1(f1):
    t1 = np.tan(f1)
    return(t1)




def gk2fl(x,y,l0,a,e2):
    b2=a**2*(1-e2)
    e_2=(a**2-b2)/b2
    A0=1-e2/4-3*e2**2/64-5*e2**3/256
    fl = f1(x, a, e2)
    N=Np(fl, a, e2)
    M=Mp(fl, a, e2)
    t=tan(fl)
    n2=e_2*(cos(fl)**2)
    f= fl - (((y**2) * t)/(2 * M * N)) * (1 - ((y**2)/(12*N**2)) * (5 + 3 * t**2 + n2 - 9 * n2 * t**2 - 4 * n2**2) + ((y**4)/(360 * N**4)) * (61 + 90 * t**2 + 45 * t**4))
    l=l0+(y/(N*cos(fl)))*(1-((y**2)/(6*N**2))*(1+2*t**2+n2)+((y**4)/(120*N**4))*(5+28*t**2+24*t**4+6*n2+8*n2*t**2))
    return(f,l)



# przeliczanie do układów PL2000 i PL1992 i odwrotne

def gk2pl2000(xgk, ygk, nr):
    m0_2000 = 0.999923
    x_2000 = xgk * m0_2000
    y_2000 = ygk * m0_2000 + ((nr*1000000)+500000)
    return(x_2000, y_2000)


def gk2pl92(xgk, ygk):
    m0_92 = 0.9993
    x_92 = xgk * m0_92 - 5300000
    y_92 = ygk * m0_92 + 500000
    return(x_92, y_92)


def pl20002gk(x_2000, y_2000, nr):
    m0_2000 = 0.999923
    xgk = x_2000 / m0_2000
    ygk = (y_2000 - 500000 - (nr * 1000000)) / m0_2000
    return(xgk, ygk)


def pl922gk(x_92, y_92):
    m0_92 = 0.9993
    xgk = (x_92 + 5300000) / m0_92
    ygk = (y_92 - 500000) / m0_92
    return(xgk, ygk)

def fl2pl92(f, l, a, e2):
    l0 = np.radians(19)
    xgk, ygk = fl2gk(f, l, l0, a, e2)
    x92, y92 = gk2pl92(xgk, ygk)
    return(xgk, ygk, x92, y92)

def pl922fl(x92, y92, a, e2):
    l0 = np.radians(19)
    xgk, ygk = pl922gk(x92, y92)
    f, l = gk2fl(xgk, ygk, l0, a, e2)
    return(xgk, ygk, f, l)
##nr 1 cyfra wspolrzednej y
def fl2pl2000(f, l, nr, a, e2): 
    if nr == 5:
        l0 = np.radians(15)
    elif nr == 6:
        l0 = np.radians(18)
    elif nr == 7:
        l0 = np.radians(21)
    elif nr == 8:
        l0 = np.radians(24)
    xgk, ygk = fl2gk(f, l, l0, a, e2)
    x2000, y2000 = gk2pl2000(xgk, ygk, nr)
    return(xgk, ygk, x2000, y2000)

def pl20002fl(x2000, y2000, nr, a, e2):
    if nr == 5:
        l0 = np.radians(15)
    elif nr == 6:
        l0 = np.radians(18)
    elif nr == 7:
        l0 = np.radians(21)
    elif nr == 8:
        l0 = np.radians(24)
    xgk, ygk = pl20002gk(x2000, y2000, nr)
    f, l = gk2fl(xgk, ygk, l0, a, e2)
    return(xgk, ygk, f, l)


    
# redukcje odległosci    

def mgk(xgk, ygk, a, e2):
    fi1 = f1(xgk, a, e2)
    R1 = np.sqrt(Np(fi1, a, e2) * Mp(fi1, a, e2))
    mgk = 1 + (ygk**2 / (2 * R1**2)) + (ygk**4 / (24 * R1**4))
    return(mgk)

def mgk1(dl, f, a, e2):
    n2 = n_2(a, e2, f)
    t = np.tan(f)
    m_gk = 1 + (dl**2 / 2) * (np.cos(f))**2 * (1 + n2) + (dl**4 / 24) * (np.cos(f))**4 * (5 - 4 * t**2)
    return(m_gk)


def m_ukladu(mgk, m0):
    m_ukladu = mgk * m0
    return(m_ukladu)


def redukcja(s_ab, ya, yb, Rm):   # musi być srednie fi
    r = s_ab * ((ya**2 + ya * yb + yb**2) / (6 * Rm**2))
    return(r)

def sgk(s_el, r):
    s = s_el + r
    return(s)

def s0_el(ha, hb, s_pom, Rm):
    s0 = np.sqrt((s_pom**2 - (hb - ha)**2) / ((1 + (hb / Rm)) * (1 + (ha / Rm))))
    return(s0)


def s_elips(Rm, s0):
    s_el = 2 * Rm * np.arcsin(s0 / 2 * Rm)


# redukcja azymutu


def gamma(f, dl, a, e2):     # kiedy znamy l0
    t = np.tan(f)
    n2 = n_2(a, e2, f)
    gamma = dl * np.sin(f) + (dl**3 / 3) * np.sin(f) * (np.cos(f))**2 * (1 + 3 * n2 + 2 * n2**2) + (dl**5 / 15) * np.sin(f) * (np.cos(f))**4 * (2 - t**2)
    

def gamma1(xgk, ygk, a, e2):
    f = f1(xgk, a, e2)
    N1 = Np(f,a,e2)
    t1 = np.tan(f)
    b2 = a**2 * (1-e2)
    eprim2 = (a**2 - b2)/b2
    n2 = eprim2 * (np.cos(f))**2
    
    gamma = ((ygk/N1) * t1) *(1 - ((ygk**2)/(3 * N1**2)) * (1 + t1**2 - n2 - 2 * n2**2) \
            + ((ygk**4)/(15 * N1**4)) * (2 + 5*t1**2 + 3*t1**4))
    return gamma


def delta_AB(xa, ya, xb, yb, Rm):       # do redukcji kierunku, azymutu
    d_ab = ((xb - xa) * (2 * ya + yb)) / (6 * Rm**2)
    return(d_ab)


def delta_BA(xa, ya, xb, yb, Rm):     # do redukcji kierunku, azymutu
    d_ba = ((xa - xb) * (2 * yb + ya)) / (6 * Rm**2)
    return(d_ba)

def alfa(xa, ya, xb, yb):
    dx_ab = xb - xa
    dy_ab = yb - ya
    dx_ba = xa - xb
    dy_ba = ya - yb
    alfa_ab = np.arctan2(dy_ab,dx_ab)
    if alfa_ab < 0:
        alfa_ab = alfa_ab + 2 * np.pi
    alfa_ba = np.arctan2(dy_ba, dx_ba)
    if alfa_ba < 0:
        alfa_ba = alfa_ba + 2 * np.pi
    print('alfa_ab = ', dms(alfa_ab))
    print('alfa_ba = ', dms(alfa_ba))
    return(alfa_ab, alfa_ba)

def azymut_el(gamma, alfa, delta):
    az = gamma + alfa + delta
    if az < 0:
        az = az + 2 * np.pi
    return az
    
    


def transf(x, y, z, kx, ky, kz, alfa, beta, gamma, x0, y0, z0):
    r1 = np.array([x, y, z])
    M = np.array([[kx, gamma, -beta],
                  [-gamma, ky, alfa],
                  [beta, -alfa, kz]])
    r0 = np.array([x0, y0, z0])
    r2 = r1 + M @ r1 + r0
    return(r2)

    
# dla punktu A
#dX = np.array([(xb - xa), (yb - ya), (zb - za)])
    







        
     

















































        
        
        
        

        
    
    

        
        
    
    
    


    
    


