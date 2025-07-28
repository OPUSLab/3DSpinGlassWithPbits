import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes
import plot_config as pltc
from pathlib import Path

#matplotlib.rcParams['pdf.fonttype']=42

ebparams = {
    'capsize': 1.0,
    'elinewidth': 0.75,
    'capthick': 0.75,
}


SQA_result_dir = f'DT_SQA_results/extracted_data/'
SQA_dir_path = Path(__file__).resolve().parent.parent.joinpath(SQA_result_dir)

APT_result_dir = f'APT_results/extracted_data/'
APT_dir_path = Path(__file__).resolve().parent.parent.joinpath(APT_result_dir)

data1 = np.loadtxt(SQA_dir_path.joinpath('L15_DTSQA_orig_132Replicas_bs_ci_logical.txt'))  # DTSQA: L = 15, R = 132
data2 = np.loadtxt(APT_dir_path.joinpath('L15_ICM_w_ICM_4xReplicas_bs_ci_old.txt'))

data3 = np.loadtxt(APT_dir_path.joinpath('L08_ICM_w_ICM_4xReplicas_bs_ci_logical_SSR1.txt'))
data4 = np.loadtxt(APT_dir_path.joinpath('L10_ICM_w_ICM_4xReplicas_bs_ci_logical_SSR1.txt'))
data5 = np.loadtxt(APT_dir_path.joinpath('L12_ICM_w_ICM_4xReplicas_bs_ci_logical_SSR1.txt'))
data6 = np.loadtxt(APT_dir_path.joinpath('L16_ICM_w_ICM_4xReplicas_bs_ci_logical_SSR1.txt'))


fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8.0, 2.8), dpi=150, facecolor='w', edgecolor='k')

ax1.errorbar(data1[:, 0], data1[:, 1], yerr=(data1[:, 2], data1[:, 3]), 
                label='DT-SQA (R = 132)', marker='d', color=pltc.color[2], 
                markeredgecolor=pltc.color[2], markersize=6.0,
                markerfacecolor=pltc.color[2]+f'55', linestyle='', **ebparams )
ax1.errorbar(data2[:, 0], data2[:, 1], yerr=(data2[:, 2], data2[:, 3]), 
                label='APT (R = 132)', marker='s', color=pltc.color[1], 
                markeredgecolor=pltc.color[1], markersize=5.0,
                markerfacecolor=pltc.color[1]+f'55', linestyle='', **ebparams )

ax1.set_xlabel(f'Annealing time, $t_a$ (MCS)')
ax1.set_ylabel(r'Residual energy, $\rho_E^f$')
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.tick_params(axis='both', which='both', direction='in')
ax1.minorticks_on()
ax1.legend(frameon=False, shadow=False)
ax1.set_ylim([2e-4, 1.5])

mu_ext = 6.34
b_ext = 3.0

ax2.errorbar(data3[:-14, 0] / (np.power(float(8), mu_ext)), data3[:-14, 1] * (np.power(float(8), b_ext)),
                yerr=(data3[:-14, 2] * (np.power(float(8), b_ext)), data3[:-14, 3] * (np.power(float(8), b_ext))), 
                label='L = 8', marker='o', color=pltc.color[0], 
                markeredgecolor=pltc.color[0], markersize=5.0,
                markerfacecolor=pltc.color[0]+f'55', linestyle='', **ebparams )
ax2.errorbar(data4[:-20, 0] / (np.power(float(10), mu_ext)), data4[:-20, 1] * (np.power(float(10), b_ext)),
                yerr=(data4[:-20, 2] * (np.power(float(10), b_ext)), data4[:-20, 3] * (np.power(float(10), b_ext))), 
                label='L = 10', marker='s', color=pltc.color[1], 
                markeredgecolor=pltc.color[1], markersize=5.0,
                markerfacecolor=pltc.color[1]+f'55', linestyle='', **ebparams )
ax2.errorbar(data5[:-12, 0] / (np.power(float(12), mu_ext)), data5[:-12, 1] * (np.power(float(12), b_ext)),
                yerr=(data5[:-12, 2] * (np.power(float(12), b_ext)), data5[:-12, 3] * (np.power(float(12), b_ext))), 
                label='L = 12', marker='d', color=pltc.color[2], 
                markeredgecolor=pltc.color[2], markersize=5.0,
                markerfacecolor=pltc.color[2]+f'55', linestyle='', **ebparams )
ax2.errorbar(data6[:-9, 0] / (np.power(float(16), mu_ext)), data6[:-9, 1] * (np.power(float(16), b_ext)),
                yerr=(data6[:-9, 2] * (np.power(float(16), b_ext)), data6[:-9, 3] * (np.power(float(16), b_ext))), 
                label='L = 16', marker='^', color=pltc.color[3], 
                markeredgecolor=pltc.color[3], markersize=5.0,
                markerfacecolor=pltc.color[3]+f'55', linestyle='', **ebparams )

ax2.set_xlabel('$t_aL^{-\mu}$')
ax2.set_ylabel(r'$\rho_E^fL^{b}$')
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.tick_params(axis='both', which='both', direction='in')
ax2.minorticks_on()
ax2.legend(shadow=False, frameon=False, loc='lower left')

fig.tight_layout()
#plt.savefig('Fig_3.pdf', dpi=300, bbox_inches='tight', transparent=False)
plt.show(block=True)
