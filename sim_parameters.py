import numpy as np
import pandas as pd
import os
import helper_functions as hf

# Settings-Objekt beinhaltet Simulationsparameter und Abscheidergeometrie

class Settings():

    def __init__(self):
        self.N_h = 50 # number of grid points in height direction
        self.N_x = 100 # number of grid points in x direction
        self.N_dpz_loop = 20 # number of iterations in dpz-loop

# Folgende Klasse beinhaltet weitere Parameter und liest Stoffwerte und Messdaten aus Excel-Dateien aus
class Substance_System():

    def __init__(self):
        # constants
        self.g = 9.81
        self.eps_p = 0.9 # hold up of dense-packed zone [-]
        self.eps_di = 1
        self.H_cd = 1e-20
        self.eta_v = 23e-3 # Grenzflächenverzögerungsviskosität [Pas] (Henschke1995)
        
        # settler geometry
        self.L = 1 # length of settler
        self.D = 0.2 # Durchmesser
        self.D_ein = 0.017 # Einlaufdurchmesser

        # Substance parameters
        self.rho_c = 0      # Density of conti phase [kg/m³]
        self.rho_d = 0      # Density of dispers phase [kg/m³]
        self.sigma = 0      # Surface tension [N/m]
        self.eta_c = 0      # Viscosity conti phase [Pas]
        self.eta_d = 0      # Viscosity disp phase [Pas]
        self.eps_0 = 0          # Hold-Up Feed [-]
        
        self.r_s_star = 0   # r_s_star parameter [-]
        self.h_p_star = 0   # h_p_star parameter [-]
        
        # operating parameters
        self.dV_ges = 0     # Total Volume Flow Rate [m^3/s]
        self.phi_0 = 0      # Sauter-Diameter at the beginning [m]

        # Calculated parameters
        self.delta_rho = 0      # Density difference [kg/m^3]
        self.s = 0              # slip parameter [-]
        self.eta_dis = 0        # Viscosity dpz [-]
        self.light_in_heavy = 0 # boolean whether light phase is dispersed in heavy phase

        # Data
        self.x_exp = []
        self.h_p_exp = []
        self.x_sim = []
        self.h_p_sim = []

    def update(self,excel_file):

        import_data = pd.read_excel(os.path.join('Input', excel_file), 'Parameters', index_col=0)

        # Update experimental parameter
        self.rho_c = import_data["ρ_c"]["Value"]
        self.rho_d = import_data["ρ_d"]["Value"]
        self.sigma = import_data["σ"]["Value"]
        self.eta_c = import_data["η_c"]["Value"]
        self.eta_d = import_data["η_d"]["Value"]
        self.eps_0 = import_data["eps_0"]["Value"]
        self.dV_ges = import_data["dV_ges"]["Value"] / 3.6 * 1e-6     # convert into SI unit
        self.phi_0 = import_data["phi_0"]["Value"]
        self.r_s_star = import_data["r_s_star"]["Value"]
        self.h_p_star = import_data["h_p_star"]["Value"]
        self.D = import_data["D"]["Value"]
        self.D_ein = import_data["D_ein"]["Value"]
        self.L = import_data["L"]["Value"]

        import_data = pd.read_excel(os.path.join('Input', excel_file), 'DataExp')

        # Update results experiment
        self.x_exp = np.array(import_data["x_exp"])
        self.h_p_exp = np.array(import_data["h_p_exp"])

        import_data = pd.read_excel(os.path.join('Input', excel_file), 'DataSim')

        # Update results experiment
        self.x_sim = np.array(import_data["x_sim"])
        self.h_p_sim = np.array(import_data["h_p_sim"])


        self.delta_rho = abs(self.rho_c - self.rho_d)
        self.s = 1 - np.exp(-12300 * self.rho_c / self.delta_rho * self.sigma ** 3)
        # self.s = 1
        self.eta_dis = hf.yaron(self.eta_c, self.eta_d, self.eta_v, self.eps_p)
        
        self.name = excel_file[:-5]

        print("Updated Parameters with Excel")
        
    def set_operating_point(self, dV_ges, phi_0, eps_0):
        self.dV_ges = dV_ges
        self.phi_0 = phi_0
        self.eps_0 = eps_0
    
    def set_coalescence_parameter(self, r_s_star):
        self.r_s_star = r_s_star
