import numpy as np
import pandas as pd
import joblib
from sim_run_MA_Ivan import run_sim

N_CPU = 4

experiments = "detail_V_dis"   # "main" for niba and ye or "sozh" for tests from AVT.FVT lab.



df = pd.read_excel("Input/data_main.xlsx", sheet_name=experiments)
exp = df['exp'].tolist()
phi_0 = df['phi_0'].tolist()
dV_ges = df['dV_ges'].tolist()
eps_0 = df['eps_0'].tolist()


def parallel_simulation(params):
    exp, phi_0, dV_ges, eps_0 = params
    print(f"Start simulation with exp={exp}, phi_0={phi_0}, dV_ges={dV_ges}, eps_0={eps_0}")
    try:
        Sim = run_sim(exp, phi_0, dV_ges, eps_0)
        return {'exp': exp, 'phi_0': phi_0, 'dV_ges': dV_ges, 'eps_0': eps_0,
                'V_dis_total': Sim.V_dis,'h_dpz':Sim.h_dpz , 'l_dpz': Sim.l_dpz, 'status': 'success'}
    except Exception as e:
        print(f"Simulation failed for exp={exp}, phi_0={phi_0}, dV_ges={dV_ges}, eps_0={eps_0}: {str(e)}")
        return {'exp': exp, 'phi_0': phi_0, 'dV_ges': dV_ges, 'eps_0': eps_0, 'error': str(e), 'status': 'failed'}

if __name__ == "__main__":
    parameters = [(exp[i], phi_0[i], dV_ges[i], eps_0[i]) for i in range(len(exp))]
    
    results = joblib.Parallel(n_jobs=N_CPU, backend='multiprocessing')(joblib.delayed(parallel_simulation)(param) for param in parameters)
    
    # Save results
    df_results = pd.DataFrame(results)
    h_dpz_columns = pd.DataFrame(df_results['h_dpz'].tolist())   # Convert h_dpz (list of arrays) into separate columns
    h_dpz_columns.columns = [f'h_dpz_{i}' for i in range(h_dpz_columns.shape[1])]
    df_results = df_results.drop(columns=['h_dpz'])
    l_dpz_columns = pd.DataFrame(df_results['l_dpz'].tolist())   # Convert l_dpz (list of arrays) into separate columns
    l_dpz_columns.columns = [f'l_dpz_{i}' for i in range(l_dpz_columns.shape[1])]
    df_results = df_results.drop(columns=['l_dpz'])
    df_results = pd.concat([df_results, h_dpz_columns, l_dpz_columns], axis=1)  # Concatenate V_dis columns with the main result dataframe

    df_results.to_csv('simulation_results_parallel_evaluation_detail_new_fit.csv', index=False)
    print("Alle Simulationen abgeschlossen. Ergebnisse gespeichert.")