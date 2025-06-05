import sim_model as sm
import sim_parameters as sp
import helper_functions as hf
import pandas as pd
import numpy as np

def init_sim(filename):
    Set = sp.Settings()
    SubSys = sp.Substance_System()
    SubSys.update(filename)
    return sm.Simulation(Set, SubSys)

if __name__ == '__main__':

    # filename = 'Hexan_1_1_o_in_w.xlsx'
    filename = 'Paraffin_flut_20C.xlsx' 
    Sim = init_sim(filename)
    Sim.Sub.dV_ges = 280 / 3.6 * 1e-6  ############################################################################################# Ändern
    file_path = "c:\\Users\\Asus\\OneDrive\\Documentos\\RWTH\\Semester 4 (SS25)\\MA\\W9_Vergleichen_Vdis_sep.Effi\\Vdis_Ye.xlsx"
    sheet = "FilteredDataQ280"  #################################################################################################### Ändern
    df_input = pd.read_excel(file_path, sheet_name=sheet) 
    phi_in0 = df_input['d_32in'].to_numpy() 

    V_dis = np.zeros(len(phi_in0))
    H_DPZ = np.zeros(len(phi_in0))
    L_DPZ = np.zeros(len(phi_in0))
    L_ein = np.zeros(len(phi_in0))
    for i in range(len(phi_in0)):
        Sim.Sub.phi_0 = phi_in0[i]
        if (Sim.Sub.phi_0 == 0) :
            V_dis[i] = 0
            H_DPZ[i] = 0
            L_DPZ[i] = 0
            L_ein[i] = 0
        else:
            Sim.calc_DPZ_brentq(report=False)
            V_dis[i] = hf.getVdis(Sim)
            H_DPZ[i] = Sim.h_p0
            L_DPZ[i] = Sim.l_coal+Sim.l_in
            L_ein[i] = Sim.l_in


    df = pd.DataFrame({
        'H_DPZ': H_DPZ,
        'L_DPZ': L_DPZ,
        'L_ein': L_ein,
        'V_dis': V_dis
    })
    file_path_results = "c:\\Users\\Asus\\OneDrive\\Documentos\\RWTH\\Semester 4 (SS25)\\MA\W9_Vergleichen_Vdis_sep.Effi\\Vdis_Model_1_results.xlsx"
    sheet = "Results_Q280" ############################################################################################################################ Ändern
    with pd.ExcelWriter(file_path_results, engine='openpyxl', mode='a') as writer:
        df.to_excel(writer, sheet_name=sheet, index=False)