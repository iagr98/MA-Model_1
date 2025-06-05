import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve
from scipy.optimize import brentq
from scipy.optimize import minimize_scalar
from scipy.optimize import minimize
import helper_functions as hf
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Diese Datei beinhaltet alle Berechnungen und Plotting-Funktionen
class Simulation():

    def __init__(self, Settings, Substance_System):   # inputs need to be Objects from classes Settings and Substance_System

        self.Set = Settings
        self.Sub = Substance_System
        
        self.x = np.linspace(0, self.Sub.L, self.Set.N_x)
        self.dl = self.x[1]-self.x[0]

        self.l_in = 0
        self.l_coal = 0
        self.dV_ges_flood = 0
        self.h_p0_flood = 0
        self.h_p0 = 0

    # Funktion für Koaleszenzzeit
    def tau(self, h, d_32, ID, sigma, r_s_star):

        La_mod = (self.Sub.g * self.Sub.delta_rho / sigma) ** 0.6 * d_32 * h ** 0.2

        R_F = d_32 * (1 - (4.7 / (4.7 + La_mod))) ** 0.5

        if ID == "d":
            R_F = 0.3025 * R_F
        else:
            R_F = 0.5240 * R_F

        R_a = 0.5 * d_32 * (1 - (1 - 4.7 / (4.7 + La_mod)) ** 0.5)

        tau = 7.65 * self.Sub.eta_c * (R_a ** (7 / 3)
            / (self.Sub.H_cd ** (1 / 6) * sigma ** (5 / 6) * R_F * r_s_star))

        return tau
    
    # Funktion für die Berechnung der Einlauflänge
    def calc_L_in(self, H_p0):
        dV_ges = self.Sub.dV_ges
        D_ein = self.Sub.D_ein
        eps_0 = self.Sub.eps_0
        phi_0 = self.Sub.phi_0
        rho_mid = eps_0 * self.Sub.rho_d + (1 - eps_0) * self.Sub.rho_c
        v_ein = dV_ges / (np.pi * D_ein ** 2 / 4)
        v_A = dV_ges / (np.pi * self.Sub.D ** 2 / 4)
        Re_ein = v_ein * D_ein * rho_mid / self.Sub.eta_dis
        Re_A = v_A * self.Sub.D * rho_mid / self.Sub.eta_dis
        Ar_dis = self.Sub.delta_rho * self.Sub.g * rho_mid * phi_0**3 / self.Sub.eta_dis**2
        L_in = phi_0 * 43.7 * (phi_0/(phi_0 + H_p0))**0.4 * (Re_ein*Re_A/Ar_dis)**0.5 * (1/(1-eps_0))**0.2 * (self.Sub.delta_rho/rho_mid)**0.2 * (phi_0/self.Sub.D)**0.1
        return L_in
    
    def calc_DPZ(self, report=False):
        # initial values
        h_p0 = self.Sub.D/2
        dh_p0 = h_p0/2
        phi_32 = np.zeros(self.Set.N_h)
        dl = self.dl
        
        for k in range(self.Set.N_dpz_loop):
            l = 0
            i_unten = 0
            disp_laenge = True
            dV_dis = self.Sub.dV_ges * self.Sub.eps_0/self.Sub.eps_p
            dh = h_p0/self.Set.N_h
            h_p = np.zeros(self.Set.N_x)
            h_p[0] = h_p0
            L_ein = self.calc_L_in(h_p0)
            L_eff = self.Sub.L - L_ein
            for i in range(self.Set.N_h):
                phi_32[i] = self.Sub.phi_0
            j = 0
            while disp_laenge:
                dt = h_p[j] * self.Sub.D * dl / dV_dis
                j += 1
                l = l + dl
                h_p[j] = (self.Set.N_h - i_unten + 0.5) * dh
                # Berechnung der Koaleszenz 
                for i in range(self.Set.N_h):
                    h_py = (self.Set.N_h - i + 0.5) * dh
                    if h_py < 0.01*phi_32[i]:
                        h_py = 0.01*phi_32[i]
                    tau_dd = self.tau(h_py, phi_32[i], "d", self.Sub.sigma, self.Sub.r_s_star)
                    phi_32[i] = phi_32[i] + dt * phi_32[i] / (6 * tau_dd)
                    if i==i_unten:
                        if h_p[j] < phi_32[i_unten]/2:
                            h_p[j] = phi_32[i_unten]/2
                        tau_di = self.tau(h_p[j], phi_32[i_unten], "i", self.Sub.sigma, self.Sub.r_s_star)
                        ddV_dis = dl * 2 * self.Sub.D * self.Sub.eps_di * phi_32[i_unten] / (3 * tau_di * self.Sub.eps_p)
                d_hp = ( 126 * (self.Sub.eta_c + self.Sub.eta_d) + 11.3 * self.Sub.s * self.Sub.eta_dis) * dV_dis * dl
                d_hp = d_hp / (h_p[j-1] * self.Sub.D**3 * self.Sub.g * self.Sub.eps_p * (1 - self.Sub.eps_p) * self.Sub.delta_rho)
                h_p[j] = h_p[j-1] - d_hp
                dV_dis = dV_dis - ddV_dis
                i_unten = round(self.Set.N_h - (self.Set.N_h - 1) * dV_dis * self.Sub.eps_p / (self.Sub.dV_ges * self.Sub.eps_0))
                dh = h_p[j] * self.Sub.dV_ges * self.Sub.eps_0 / (self.Set.N_h * dV_dis * self.Sub.eps_p)
                if dh > h_p[j]:
                    dh = h_p[j]
                if dV_dis <= 0 or h_p[j] <= 0 or l > L_eff:
                    disp_laenge = False
                    break
            # print iteration results
            if report:
                print("Iteration: ", k)
                print("dV_dis in L/h: ", dV_dis*1000*3600)
                print("h_p0: ", h_p0)
                print("l: ", l)
                print("L_ein: ", L_ein)
                print("h_p: ", h_p)
            if dV_dis >= 0:
                h_p0 = h_p0 + dh_p0
            else:
                h_p0 = h_p0 - dh_p0
            dh_p0 = dh_p0 / 2
            
        # safe results
        self.Sub.h_p_sim = h_p
        self.Sub.x_sim = self.x
        self.l_coal = l
        self.l_in = L_ein
        self.dV_dis_end = dV_dis
        self.h_p0 = h_p0
        
        return
    
    def calc_DPZ_brentq(self, report=False):
        # initial values
        h_p0 = self.Sub.D/2
        dh_p0 = h_p0/2
        phi_32 = np.zeros(self.Set.N_h)
        dl = self.dl
        
        def calc_dpz(h_p0):
            l = 0
            i_unten = 0
            disp_laenge = True
            dV_dis = self.Sub.dV_ges * self.Sub.eps_0/self.Sub.eps_p
            dh = h_p0/self.Set.N_h
            h_p = np.zeros(self.Set.N_x)
            h_p[0] = h_p0
            L_ein = self.calc_L_in(h_p0)
            L_eff = self.Sub.L - L_ein
            for i in range(self.Set.N_h):
                phi_32[i] = self.Sub.phi_0
            j = 0
            while disp_laenge:
                dt = h_p[j] * self.Sub.D * dl / dV_dis
                j += 1
                l = l + dl
                h_p[j] = (self.Set.N_h - i_unten + 0.5) * dh
                # Berechnung der Koaleszenz 
                for i in range(self.Set.N_h):
                    h_py = (self.Set.N_h - i + 0.5) * dh
                    if h_py < 0.01*phi_32[i]:
                        h_py = 0.01*phi_32[i]
                    tau_dd = self.tau(h_py, phi_32[i], "d", self.Sub.sigma, self.Sub.r_s_star)
                    phi_32[i] = phi_32[i] + dt * phi_32[i] / (6 * tau_dd)
                    if i==i_unten:
                        if h_p[j] < phi_32[i_unten]/2:
                            h_p[j] = phi_32[i_unten]/2
                        tau_di = self.tau(h_p[j], phi_32[i_unten], "i", self.Sub.sigma, self.Sub.r_s_star)
                        ddV_dis = dl * 2 * self.Sub.D * self.Sub.eps_di * phi_32[i_unten] / (3 * tau_di * self.Sub.eps_p)
                d_hp = ( 126 * (self.Sub.eta_c + self.Sub.eta_d) + 11.3 * self.Sub.s * self.Sub.eta_dis) * dV_dis * dl
                d_hp = d_hp / (h_p[j-1] * self.Sub.D**3 * self.Sub.g * self.Sub.eps_p * (1 - self.Sub.eps_p) * self.Sub.delta_rho)
                h_p[j] = h_p[j-1] - d_hp
                dV_dis = dV_dis - ddV_dis
                i_unten = round(self.Set.N_h - (self.Set.N_h - 1) * dV_dis * self.Sub.eps_p / (self.Sub.dV_ges * self.Sub.eps_0))
                dh = h_p[j] * self.Sub.dV_ges * self.Sub.eps_0 / (self.Set.N_h * dV_dis * self.Sub.eps_p)
                if dh > h_p[j]:
                    dh = h_p[j]
                if dV_dis <= 0 or h_p[j] <= 0 or l > L_eff:
                    disp_laenge = False
                    break
            # safe results
            self.Sub.h_p_sim = h_p
            self.Sub.x_sim = self.x
            self.l_coal = l
            self.l_in = L_ein
            self.dV_dis_end = dV_dis
            # print iteration results
            if report:
                print("dV_dis end in L/h: ", dV_dis*1000*3600)
                print("h_p0: ", h_p0)
                print("l: ", l)
                print("L_ein: ", L_ein)
                print("h_p: ", h_p)
            return dV_dis
        
        h_p0, r = brentq(calc_dpz, self.Sub.D/100, self.Sub.D, full_output=True, xtol=self.Sub.D/1000)
        # run final calculation
        if r.converged:
            self.h_p0 = h_p0
            calc_dpz(h_p0)
        else:
            print("Brentq did not converge")
        return
    
    def find_flooding_flowrate(self, h_p_max):
        """Find flooding flowrate for given h_p_max

        Args:
            h_p_max (float): critical height of DPZ height in m

        Returns:
            void: Finds flooding flowrate and stores it in self.dV_ges_flood
        """
        self.h_p0_flood = h_p_max
        def residual(dV_ges):
            self.Sub.set_operating_point(dV_ges[0], self.Sub.phi_0, self.Sub.eps_0)
            self.calc_DPZ_brentq()
            res = np.max(self.Sub.h_p_sim) - h_p_max
            return res**2 # square residual to make it positive
        res = minimize(residual, self.Sub.dV_ges, bounds=((self.Sub.dV_ges/10, self.Sub.dV_ges*10),), method='Nelder-Mead', tol=1e-8, options={'disp': True})
        # recalculation of DPZ with flooding flowrate
        self.Sub.set_operating_point(res.x[0], self.Sub.phi_0, self.Sub.eps_0)
        self.calc_DPZ_brentq()
        self.dV_ges_flood = res.x[0]
    
# Plot h_p vs x
def plot_h_p(Sim, label='1', henschkeData=True):
    fig, ax = plt.subplots(figsize=(8,6))
    if henschkeData:
        ax.plot(Sim.Sub.x_exp, Sim.Sub.h_p_exp, 'o', label='Henschke')
    ax.plot((np.concatenate((np.linspace(0,Sim.l_in,2),Sim.Sub.x_sim + Sim.l_in)))*1000, (np.concatenate((np.ones(2)*Sim.h_p0,Sim.Sub.h_p_sim)))*1000, label='Simulation ' + label)    
    ax.axvline(x=1000*Sim.l_in, color='r', linestyle='--', label='Einlauflänge')
    ax.axhline(y=0, color='k')
    ax.set_xlim(0, Sim.Sub.L * 1000)
    ax.set_xlabel('Abscheiderlänge / mm', size=12)
    ax.set_ylabel('Höhe der DGTS / mm', size=12)
    plt.legend()
    plt.tight_layout()
    plt.show()

# Kriegt mehrere Simulationsobjekte in Liste übergeben und ruft einzeln calc_comparison auf, um diese zu plotten
def plot_comparison(Sims, labels=['1','2','3','4', '5'], legend_title=None, title=None,
                    henschkeData=True, xlim=None, figsize=(8,6)):

    fig, ax = plt.subplots(figsize=figsize)
    xmax = 0
    ymax = 0
    single = False
    if len(Sims) == 1:
        single = True
    for i in range(len(Sims)):
        xmax_new, ymax_new = Sims[i].calc_comparison(ax, i, single, labels, henschkeData)
        if xmax_new > xmax:
            xmax = xmax_new
        if ymax_new > ymax:
            ymax = ymax_new
    if title != None:
        ax.set_title(title)
    # ax.set_xlim(0, xmax)
    if xlim == None:
        ax.set_xlim(0, Sims[0].Sub.L * 1000)
    else:
        ax.set_xlim(xlim)
    ax.set_ylim(0, ymax)
    ax.set_xlabel('Abscheiderlänge / mm', size=12)
    ax.set_ylabel('Höhe der DGTS / mm', size=12)
    plt.tick_params(axis='x', top=True, direction='in')
    plt.tick_params(axis='y', right=True, direction='in')
    if legend_title == None:
        plt.legend(frameon=False)
    else:
        plt.legend(title = legend_title, frameon=False)

    plt.tight_layout()
    plt.show()

# erstellt den Sensitivitätsplot
def plot_sensitivity(Sims, parameters, title=None, x_label='x-Achse', xlim=None):
    fig, ax = plt.subplots(figsize=(5, 4.4))
    V_tot = np.zeros(len(parameters))
    colors = ['b', 'r', 'g', 'm', 'k', 'y', 'c']
    for i in range(len(Sims)):
        Vdis_end = Sims[i].V_dis[:, -1]
        if np.where(Vdis_end < 1e-8)[0].size > 0:
            last_idx = np.where(Vdis_end < 1e-8)[0][0]
            Vdis_end = Vdis_end[:last_idx]
        Vdis_tot = np.sum(Vdis_end) * 1000
        V_tot[i] = Vdis_tot
        ax.scatter(parameters[i], Vdis_tot, color=colors[i])
    ax.set_xlabel(x_label, size=12)
    ax.set_ylabel('Gesamtvolumen der DGTS / L', size=12)
    if xlim == None:
        ax.set_xlim(left=0)
    else:
        ax.set_xlim(xlim)
    ax.set_ylim(0, 1)
    if title != None:
        ax.set_title(title)
    plt.tick_params(axis='x', direction='in', top=True)
    plt.tick_params(axis='y', right=True, direction='in')
    plt.tight_layout()
    plt.show()