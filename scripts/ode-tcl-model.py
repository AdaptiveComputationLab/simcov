import os
import sys
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FixedFormatter
import numpy as np
from scipy.integrate import odeint
from scipy.interpolate import make_interp_spline, BSpline

def ode_solver(t, initial_conditions, params, ode_model):
    initT1, initT2, initE1, initE2, initI1, initI2, initVT, initVS = initial_conditions
    beta_t, beta_s, sigma, gamma_t, gamma_s, pi_t, pi_s, c, transport = params
    res = odeint(ode_model, [initT1, initT2, initE1, initE2, initI1, initI2, initVT, initVS], t, args=(beta_t, beta_s, sigma, gamma_t, gamma_s, pi_t, pi_s, c, transport))
    return res

def main(initT1, initT2, initE1, initE2, initI1, initI2, initVT, initVS, beta_t, beta_s, sigma, gamma_t, gamma_s, pi_t, pi_s, c, transport, time, ode_model, name):
    initial_conditions = [initT1, initT2, initE1, initE2, initI1, initI2, initVT, initVS]
    params = [beta_t, beta_s, sigma, gamma_t, gamma_s, pi_t, pi_s, c, transport]
    tspan = np.arange(0, time, 1)
    sol = ode_solver(tspan, initial_conditions, params, ode_model)
    T1, T2, E1, E2, I1, I2, VT, VS = sol[:, 0], sol[:, 1], sol[:, 2], sol[:, 3], sol[:, 4], sol[:, 5], sol[:, 6], sol[:, 7]
    return T2, E2, I2, VS

def ode_model(z, t, infection_t, infection_s, leave, death_t, death_s, pi_virus_production_t, pi_virus_production_s, virus_death, transport):
    T1, T2, E1, E2, I1, I2, VT, VS = z
    dT1dt = -infection_t*T1*VT
    dE1dt = infection_t*T1*VT - leave*E1
    dI1dt = leave*E1 - death_t*I1
    dVTdt = pi_virus_production_t*I1 - virus_death*VT
    
    dT2dt = -infection_s*T2*VS
    dE2dt = infection_s*T2*VS - leave*E2
    dI2dt = leave*E2 - death_s*I2
    dVSdt = pi_virus_production_s*I2 - virus_death*VS + transport*VT
    return [dT1dt, dT2dt, dE1dt, dE2dt, dI1dt, dI2dt, dVTdt, dVSdt]

def plot_graph(ax, VS, time, letter):
    VS[1:] = np.log10(VS[1:])
    VS[:100]= 0
    
    tspan = np.arange(0, time, 1)
    
    xnew = np.linspace(tspan.min(), tspan.max(), 300)
    spl_s = make_interp_spline(tspan, VS, k=3)
    power_smooth_s = spl_s(xnew)
     
    ax.plot(xnew, power_smooth_s, c='blue', label='sputum virus', linewidth=3)
    
    ax.set_title(letter, loc='left')
    ax.set_ylim(bottom=0)
    ax.set_xlim(left=0)
  
    x_formatter = FixedFormatter(['0', '2', '4', '6', '8', '10', '12', '14'])
    x_locator = FixedLocator([0, 2880, 5760, 8640, 11520, 14400, 17280, 20160])
    ax.xaxis.set_major_formatter(x_formatter)
    ax.xaxis.set_major_locator(x_locator)

    y_formatter = FixedFormatter([0, 2, 4, 6, 8])
    y_locator = FixedLocator([0, 2, 4, 6, 8])
    ax.yaxis.set_major_formatter(y_formatter)
    ax.yaxis.set_major_locator(y_locator)
    
def create_ode_txt_file(T, E, I, V, name):
    T = 480000000 - T
    time = np.array(range(20160))
    A = np.zeros(T.shape[0])
    title = name + '.txt'
    np.savetxt(title, np.transpose([time,E,I,A,T,V]), delimiter="\t", header="time\tincubating\texpressing\tapoptotic\tdead\tode-virus")

time_scale = 1440 ## rates are per day, but we want per minute rates

initT1 = 4600000
initT2 = 480000000
initE1 = initE2 = 0
initI1 = initI2 = 1
initVT = initVS = 0
virus_death_rate = 10 / time_scale
transport_rate = 0.01 / time_scale
minutes = 20160 # 14 days
leave_incubation_rate = 4 / time_scale

params_dict = {
    'patient_A': {
        'epi_death_rate_throat': 2.48 / time_scale,
        'epi_death_rate_sputum': 0.55 / time_scale,
        'infection_rate_throat': 0.0000008515 / time_scale, 
        'infection_rate_sputum': 0.0000006978 / time_scale,
        'virus_production_rate_throat': 50.93 / time_scale,
        'virus_production_rate_sputum': 0.41 / time_scale
    },
    'patient_B': {
        'epi_death_rate_throat': 2.31 / time_scale,
        'epi_death_rate_sputum': 0.54 / time_scale,
        'infection_rate_throat': 0.0000005195 / time_scale,
        'infection_rate_sputum': 0.0000006886 / time_scale,
        'virus_production_rate_throat': 50.99 / time_scale,
        'virus_production_rate_sputum': 0.35 / time_scale
    },
    'patient_C': {
        'epi_death_rate_throat': 2.38 / time_scale,
        'epi_death_rate_sputum': 0.81 / time_scale,
        'infection_rate_throat': 0.0000013861 / time_scale,
        'infection_rate_sputum': 0.0000017973 / time_scale,
        'virus_production_rate_throat': 51.07 / time_scale,
        'virus_production_rate_sputum': 0.37 / time_scale
    },
    'patient_D': {
        'epi_death_rate_throat': 2.44 / time_scale,
        'epi_death_rate_sputum': 0.71 / time_scale,
        'infection_rate_throat': 0.0000005308 / time_scale,
        'infection_rate_sputum': 0.0000009590 / time_scale,
        'virus_production_rate_throat': 53.34 / time_scale,
        'virus_production_rate_sputum': 0.37 / time_scale
    },
    'patient_E': {
        'epi_death_rate_throat': 1.98 / time_scale,
        'epi_death_rate_sputum': 0.53 / time_scale,
        'infection_rate_throat': 0.0000005135 / time_scale,
        'infection_rate_sputum': 0.0000007067 / time_scale,
        'virus_production_rate_throat': 49.80 / time_scale,
        'virus_production_rate_sputum': 0.34 / time_scale
    },
    'patient_F': {
        'epi_death_rate_throat': 2.06 / time_scale,
        'epi_death_rate_sputum': 0.49 / time_scale,
        'infection_rate_throat': 0.0000005177 / time_scale,
        'infection_rate_sputum': 0.0000007530 / time_scale,
        'virus_production_rate_throat': 50.15 / time_scale,
        'virus_production_rate_sputum': 0.37 / time_scale
    },
    'patient_G': {
        'epi_death_rate_throat': 0.82 / time_scale,
        'epi_death_rate_sputum': 0.42 / time_scale,
        'infection_rate_throat': 0.0000010001 / time_scale,
        'infection_rate_sputum': 0.0000009611 / time_scale,
        'virus_production_rate_throat': 50.84 / time_scale,
        'virus_production_rate_sputum': 0.38 / time_scale,
    },
    'patient_H': {
        'epi_death_rate_throat': 2.15 / time_scale,
        'epi_death_rate_sputum': 1.62 / time_scale,
        'infection_rate_throat': 0.0000004446 / time_scale,
        'infection_rate_sputum': 0.0000016857 / time_scale,
        'virus_production_rate_throat': 51.83 / time_scale,
        'virus_production_rate_sputum': 0.4 / time_scale
    }
}

fig, axes = plt.subplots(2, 4, figsize=(20,7.5))
patient_IDs = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']

plt.rcParams.update({'font.size': 20})
plt.rcParams.update({'xtick.labelsize': 75})
plt.rcParams.update({'ytick.labelsize': 75})
plt.rcParams.update({'axes.xmargin': 0})

for i in range(8):
    letter = patient_IDs[i]
    
    T, E, I, VS = main(initT1,
                                initT2,
                                initE1,
                                initE2,
                                initI1,
                                initI2,
                                initVT,
                                initVS,
                                params_dict['patient_' + letter]['infection_rate_throat'],
                                params_dict['patient_' + letter]['infection_rate_sputum'],
                                leave_incubation_rate,
                                params_dict['patient_' + letter]['epi_death_rate_throat'],
                                params_dict['patient_' + letter]['epi_death_rate_sputum'],
                                params_dict['patient_' + letter]['virus_production_rate_throat'],
                                params_dict['patient_' + letter]['virus_production_rate_sputum'],
                                virus_death_rate,
                                transport_rate,
                                minutes,
                                ode_model,
                                letter)

    plot_graph(axes[i // 4][i % 4], VS, minutes, letter)
    create_ode_txt_file(T, E, I, VS, 'patient_' + letter + '_ode')

fig.suptitle("Viral Log Over 14 Days")
fig.text(0.5, 0.04, "Days post infection", ha='center', va='center')
fig.text(0.02, 0.5, "VL in Log10 (/swab & /mL)", ha='center', va='center', rotation='vertical')

plt.subplots_adjust(left=0.05, bottom=0.1, right=0.97, top=0.88, wspace=0.2, hspace=0.3)
plt.savefig('ode-LRT-plots.png')
