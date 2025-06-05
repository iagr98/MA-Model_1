import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import joblib
from sim_run_MA_Ivan import run_sim

N_CPU = 5
N_h = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]                   # Update
var = 'N_h'                                                                 # Update


def parallel_simulation(params): 
    N_h = params                                                            # Update parameters
    print(f"Start simulation with {var}={N_h}")                             # Update parameter in second {}
    try:
        Sim = run_sim(N_h=N_h)                                              # Update inputs (VERY IMPORTANT)
        return {f"{var}": N_h,                                              # Update parameter in second place
                'V_dis': Sim.V_dis,
                'status': 'success'}    
    except Exception as e:
        print(f"Simulation failed by {var}={N_h}: {str(e)}")                # Update parameter in second {}
        return {f"{var}": N_h, 'error': str(e), 'status': 'failed'}         # Update parameter in second place           

if __name__ == "__main__":
    parameters = [N_h_value for N_h_value in N_h]                           # Update parameter var_value, var_value & var 
    
    results = joblib.Parallel(n_jobs=N_CPU)(joblib.delayed(parallel_simulation)(param) for param in parameters)
    
    # Save results
    df_results = pd.DataFrame(results)
    df_results.to_csv('simulation_results_parallel_N_h.csv', index=False)   # update name of file
    print("Alle Simulationen abgeschlossen. Ergebnisse gespeichert.")

    # Plot results
    df = pd.read_csv("simulation_results_parallel_N_h.csv")                 # update name of file
    df.columns = df.columns.str.strip()
    plt.figure(figsize=(8, 5))
    plt.plot(df['N_h'], df['V_dis'], marker='o')                      # Update parameter in first place
    # plt.xscale('log')  # da atol logarithmisch skaliert ist
    # plt.yscale('log')  # da atol logarithmisch skaliert ist
    plt.xlabel('N_h')                                                       # Change x-label
    plt.ylabel('V_dis')                                         
    plt.title(f'Gitterunabhängigkeitsanalyse ({var})')
    plt.grid(True)
    plt.tight_layout()
    plt.show()