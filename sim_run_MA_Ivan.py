import sim_model as sm
import sim_parameters as sp
import helper_functions as hf
import numpy as np

def init_sim(exp, phi_0, dV_ges, eps_0, N_h, N_x):
    if (exp == "ye"):
        filename = "Paraffin_flut_20C.xlsx"
        Set = sp.Settings(N_h=N_h, N_x=N_x)
        SubSys = sp.Substance_System(L=0.56, D=0.15, D_ein=0.04)
    elif(exp == "niba1" or exp == "niba2" or exp == "niba3" or exp == "niba4"):
        Set = sp.Settings(N_h=N_h, N_x=N_x)
        SubSys = sp.Substance_System(L=1.0, D=0.2, D_ein=0.03)
        filename = "niba_V1.xlsx" if exp == "niba1" else \
        "niba_V2.xlsx" if exp == "niba2" else \
        "niba_V3.xlsx" if exp == "niba3" else \
        "niba_V4.xlsx" if exp == "niba4" else None
    elif(exp == "2mmol_21C" or exp == "2mmol_30C" or exp == "5mmol_30C" or exp == "10mmol_21C" or exp == "10mmol_30C" or exp == "15mmol_20C" or exp == "15mmol_30C"):
        Set = sp.Settings(N_h=N_h, N_x=N_x)
        SubSys = sp.Substance_System(L=1.3, D=0.2, D_ein=0.05)
        filename = "2mmolNa2CO3_21C.xlsx" if exp == "2mmol_21C" else \
        "2mmolNa2CO3_30C.xlsx" if exp == "2mmol_30C" else \
        "5mmolNa2CO3_30C.xlsx" if exp == "5mmol_30C" else \
        "10mmolNa2CO3_21C.xlsx" if exp == "10mmol_21C" else \
        "10mmolNa2CO3_30C.xlsx" if exp == "10mmol_30C" else \
        "15mmolNa2CO3_20C.xlsx" if exp == "15mmol_20C" else \
        "15mmolNa2CO3_30C.xlsx" if exp == "15mmol_30C" else None
    else:
        print('Test does not belong to either Ye or Niba.')
    SubSys.update(filename)
    SubSys.phi_0 = phi_0
    SubSys.dV_ges = dV_ges / 3.6 * 1e-6
    SubSys.eps_0 = eps_0
    return sm.Simulation(Set, SubSys)

def run_sim(exp="ye", phi_0=610e-6, dV_ges=240, eps_0=0.2, N_h=100, N_x=700):
    Sim = init_sim(exp, phi_0, dV_ges, eps_0, N_h, N_x)
    # Sim.calc_DPZ_brentq(report=False)
    Sim.calc_DPZ(report=False) # Results in reproducicle results
    return Sim

if __name__ == '__main__':
    
    exp = "2mmol_21C"    
    phi_0 = 548e-6
    dV_ges = 991
    eps_0 = 0.5
    Sim = run_sim(exp, phi_0, dV_ges, eps_0)
    sm.plot_h_p(Sim, henschkeData=False)
    print('V_dis: ', Sim.V_dis, '[m^3]')