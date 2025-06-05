import numpy as np
from scipy.optimize import newton

# Funktion berechnet Höhe eines Kreissegments auf Basis des Kreisradius r und der Fläche A
def getHeight(A, r):
    eq = lambda h: A - r**2 * np.arccos(1 - h / r) + (r - h) * np.sqrt(2 * r * h - h**2)
    h0 = r / 2
    if A < 0:
        #print('Querschnitt kleiner Null: ' + str(A))
        return 0
    elif A > np.pi * r**2:
        #print('Querschnitt größer als zulässig: ' + str(A))
        return 2*r
    return newton(eq, h0)

def getVdis(Sim):
    l_in = Sim.l_in
    D = Sim.Sub.D
    h_p0 = Sim.h_p0
    h_p = Sim.Sub.h_p_sim
    x = Sim.Sub.x_sim
    dl = Sim.dl
    V_dis = np.zeros(len(x))
    V_dis[0] = h_p0*D*l_in
    for i in range(1, len(x)):
        if (h_p[i]<0): continue
        V_dis[i] = h_p[i]*D*dl
    V_dis = np.sum(V_dis)
    return V_dis 

def getHeightArray(A, r):
    h = np.zeros_like(A)
    for i in range(len(h)):
        h[i] = getHeight(A[i], r)
    return h

# Funktion berechnet die Fläche eines Kreissegments auf Basis des Kreisradiuses r und der Höhe h des Segments
def getArea(h, r):
    return r**2 * np.arccos(1 - h / r) - (r - h) * np.sqrt(2 * r * h - h**2)

# Funktion berechnet dyn. Viskosität in der dicht gepackten Schicht nach Modell von Yaron und Gal-Or
def yaron(eta_c, eta_d, eta_v, eps):
    al = eta_c / (eta_d + eta_v)
    ga = eps ** (1 / 3)
    omega = ((4 * ga ** 7 + 10 - (84 / 11) * ga ** 2 + 4 * al * (1 - ga ** 7)) /
             (10 * (1 - ga ** 10) - 25 * ga ** 3 * (1 - ga ** 4) + 10 * al * (1 - ga ** 3) * (1 - ga ** 7)))
    return eta_c * (1 + 5.5 * omega * eps)