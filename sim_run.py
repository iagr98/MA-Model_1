import sim_model as sm
import sim_parameters as sp
import pandas as pd
import numpy as np

def init_sim(filename):
    Set = sp.Settings()
    SubSys = sp.Substance_System()
    SubSys.update(filename)
    return sm.Simulation(Set, SubSys)

def run_floodingsearch(filename, d32, eps_0, name_addition = ''):
    Sim = init_sim(filename)
    Sim.calc_DPZ_brentq() # calculates dpz for the simulation
    
    dV_ges_flood = []
    h_p = []
    l_ein = []
    # iteration
    for i in range(len(d32)):
        Sim.Sub.set_operating_point(Sim.Sub.dV_ges, d32[i]/1e6, eps_0[i])
        Sim.find_flooding_flowrate(0.030) # calculates flooding flowrate
        print(Sim.dV_ges_flood*3600)
        dV_ges_flood.append(Sim.dV_ges_flood*3600)
        h_p.append(np.max(Sim.Sub.h_p_sim))
        l_ein.append(Sim.l_in)

    # create data frame for the results
    df = pd.DataFrame({'d32': d32, 'eps_0': eps_0, 'dV_ges_flood': dV_ges_flood, 'h_p in m': h_p, 'l_ein in m': l_ein})
    
    # save results to excel with filename in output folder
    df.to_excel("Output/" + filename[:-5] + name_addition +'_results.xlsx')
    print('Results saved to Output/' + filename[:-5] +  name_addition + '_results.xlsx')

# run
if __name__ == '__main__':
    
    d32_sensitivity = 1.2
    # some plotting functions may only work for Simulation with a Sim Object calculated with simulate_upwind() and not ivp!
    filename = 'Octanol_flut_20C_endoscope.xlsx'
    d32 = [700, 880, 905, 950]
    d32 = [x*d32_sensitivity for x in d32]
    eps_0 = [0.3, 0.3, 0.5, 0.5]
    run_floodingsearch(filename, d32, eps_0, name_addition='_sensitivity_120d32')
    
    filename = 'Octanol_flut_30C_endoscope.xlsx'
    d32 = [608, 632, 738, 865]
    d32 = [x*d32_sensitivity for x in d32]
    eps_0 = [0.3, 0.3, 0.5, 0.5]
    run_floodingsearch(filename, d32, eps_0, name_addition='_sensitivity_120d32')
    
    filename = 'Octanol_flut_40C_endoscope.xlsx'
    d32 = [585, 628, 635, 707]
    d32 = [x*d32_sensitivity for x in d32]
    eps_0 = [0.3, 0.3, 0.5, 0.5]
    run_floodingsearch(filename, d32, eps_0, name_addition='_sensitivity_120d32')
    
    filename = 'Octanol_flut_50C_endoscope.xlsx'
    d32 = [630, 668, 673, 723]
    d32 = [x*d32_sensitivity for x in d32]
    eps_0 = [0.3, 0.3, 0.5, 0.5]
    run_floodingsearch(filename, d32, eps_0, name_addition='_sensitivity_120d32')
    
    # Sim = init_sim(filename)
    # Sim.calc_DPZ_brentq() # calculates dpz for the simulation
    

    # dV_ges_flood = []
    # h_p = []
    # l_ein = []
    # # iteration
    # for i in range(len(d32)):
    #     Sim.Sub.set_operating_point(Sim.Sub.dV_ges, d32[i]/1e6, eps_0[i])
    #     Sim.find_flooding_flowrate(0.030) # calculates flooding flowrate
    #     print(Sim.dV_ges_flood*3600)
    #     dV_ges_flood.append(Sim.dV_ges_flood*3600)
    #     h_p.append(np.max(Sim.Sub.h_p_sim))
    #     l_ein.append(Sim.l_in)

    # # create data frame for the results
    # df = pd.DataFrame({'d32': d32, 'eps_0': eps_0, 'dV_ges_flood': dV_ges_flood, 'h_p in m': h_p, 'l_ein in m': l_ein})
    
    # # save results to excel with filename in output folder
    # df.to_excel("Output/" + filename[:-5] + '_results.xlsx')
    # print('Results saved to Output/' + filename[:-5] + '_results.xlsx')
    
