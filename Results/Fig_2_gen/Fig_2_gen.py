import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes
import plot_config as pltc
from pathlib import Path
import custom_fit as cfit
from scipy.special import erfcinv
import matplotlib.ticker as ticker

ebparams = {
    'capsize': 1.0,
    'elinewidth': 0.75,
    'capthick': 0.75,
}

slopes = []
ciAll = []
L = 15

SQA_result_dir = f'DT_SQA_results/extracted_data/'
SQA_dir_path = Path(__file__).resolve().parent.parent.joinpath(SQA_result_dir)

fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(8.0, 2.8), dpi=150, facecolor='w', edgecolor='k')

ttest = np.arange(10, 5001)  # for slope illustration


data50 = np.loadtxt(SQA_dir_path.joinpath('L15_DTSQA_orig_50Replicas_bs_ci_logical.txt'))  # DTSQA: L = 15, R = 132
fit_indices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
unused_indices = []

ax1.errorbar(data50[fit_indices, 0], data50[fit_indices, 1], yerr=(data50[fit_indices, 2], data50[fit_indices, 3]), 
                label='R = 50', marker='o', color=pltc.color[1], 
                markeredgecolor=pltc.color[1], markersize=6.0,
                markerfacecolor=pltc.color[1]+f'55', linestyle='', **ebparams )

pp, ci, zx, zy  = cfit.custom_fit(np.log10(data50[fit_indices,0]), np.log10(data50[fit_indices,1]))
slopes.append(pp[0])
ciAll.append(ci)

ax1.plot(ttest, 10**(pp[0]*np.log10(ttest)+pp[1]), linestyle='-',color=pltc.color[1])




data132 = np.loadtxt(SQA_dir_path.joinpath('L15_DTSQA_orig_132Replicas_bs_ci_logical.txt'))
fit_indices = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
unused_indices = [0, 1]

pp, ci, zx, zy  = cfit.custom_fit(np.log10(data132[fit_indices,0]), np.log10(data132[fit_indices,1]))
slopes.append(pp[0])
ciAll.append(ci)

data250 = np.loadtxt(SQA_dir_path.joinpath('L15_DTSQA_orig_250Replicas_bs_ci_logical.txt'))
fit_indices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
unused_indices = []

pp, ci, zx, zy  = cfit.custom_fit(np.log10(data250[fit_indices,0]), np.log10(data250[fit_indices,1]))
slopes.append(pp[0])
ciAll.append(ci)

data500 = np.loadtxt(SQA_dir_path.joinpath('L15_DTSQA_orig_500Replicas_bs_ci_logical.txt'))
fit_indices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
unused_indices = []

ax1.errorbar(data500[fit_indices, 0], data500[fit_indices, 1], yerr=(data500[fit_indices, 2], data500[fit_indices, 3]), 
                label='R = 500', marker='d', color=pltc.color[2], 
                markeredgecolor=pltc.color[2], markersize=6.0,
                markerfacecolor=pltc.color[2]+f'55', linestyle='', **ebparams )

pp, ci, zx, zy  = cfit.custom_fit(np.log10(data500[fit_indices,0]), np.log10(data500[fit_indices,1]))
slopes.append(pp[0])
ciAll.append(ci)

ax1.plot(ttest, 10**(pp[0]*np.log10(ttest)+pp[1]), linestyle='-',color=pltc.color[2])


data1000 = np.loadtxt(SQA_dir_path.joinpath('L15_DTSQA_orig_1000Replicas_bs_ci_logical.txt'))
fit_indices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
unused_indices = []

pp, ci, zx, zy  = cfit.custom_fit(np.log10(data1000[fit_indices,0]), np.log10(data1000[fit_indices,1]))
slopes.append(pp[0])
ciAll.append(ci)

data2850 = np.loadtxt(SQA_dir_path.joinpath('L15_DTSQA_orig_2850Replicas_bs_ci_logical.txt'))
fit_indices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
unused_indices = []

ax1.errorbar(data2850[fit_indices, 0], data2850[fit_indices, 1], yerr=(data2850[fit_indices, 2], data2850[fit_indices, 3]), 
                label='R = 2850', marker='s', color=pltc.color[0], 
                markeredgecolor=pltc.color[0], markersize=6.0,
                markerfacecolor=pltc.color[0]+f'55', linestyle='', **ebparams )

ax2.errorbar(data2850[fit_indices, 0], data2850[fit_indices, 1], yerr=(data2850[fit_indices, 2], data2850[fit_indices, 3]), 
                label='DT-SQA (R = 2850)', marker='s', color=pltc.color[0], 
                markeredgecolor=pltc.color[0], markersize=6.0,
                markerfacecolor=pltc.color[0]+f'55', linestyle='', **ebparams )

pp, ci, zx, zy  = cfit.custom_fit(np.log10(data2850[fit_indices,0]), np.log10(data2850[fit_indices,1]))
slopes.append(pp[0])
ciAll.append(ci)

ax1.plot(ttest, 10**(pp[0]*np.log10(ttest)+pp[1]), linestyle='-',color=pltc.color[0])
ax2.plot(ttest, 10**(pp[0]*np.log10(ttest)+pp[1]), linestyle='-',color=pltc.color[0])

data6754 = np.loadtxt(SQA_dir_path.joinpath('L15_DTSQA_orig_6754Replicas_bs_ci_logical.txt'))
fit_indices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
unused_indices = []

pp, ci, zx, zy  = cfit.custom_fit(np.log10(data6754[fit_indices,0]), np.log10(data6754[fit_indices,1]))
slopes.append(pp[0])
ciAll.append(ci)

data16000 = np.loadtxt(SQA_dir_path.joinpath('L15_DTSQA_orig_16000Replicas_bs_ci_logical.txt'))
fit_indices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
unused_indices = []

pp, ci, zx, zy  = cfit.custom_fit(np.log10(data16000[fit_indices,0]), np.log10(data16000[fit_indices,1]))
slopes.append(pp[0])
ciAll.append(ci)

ax1.set_xlabel(f'Annealing time, $t_a$ (MCS)')
ax1.set_ylabel(r'Residual energy, $\rho_E^f$')
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.tick_params(axis='both', which='both', direction='in')
ax1.minorticks_on()
ax1.legend(frameon=False, shadow=False)
ax1.set_ylim([2e-3, 1.0])
ax1.set_xlim([5, 5e3])


dataQPU = np.loadtxt('L15_QA_bs_ci_Embedded.txt')
inds = np.array([ True,  True,  True,  True,  True,  True,  True,  True,  True,  True,  True,  True,
  True,  True,  True,  True, False, False, False, False, False, False, False, False,
 False, False, False, False, False, False, False, False, False, False, False, False,
 False, False, False, False, False, False, False, False, False, False, False, False,
 False, False])

ninds = np.array([False, False, False, False, False, False, False, False, False, False, False, False,
 False, False, False, False,  True,  True,  True,  True,  True,  True,  True,  True,
  True,  True,  True,  True,  True,  True,  True,  True,  True,  True,  True,  True,
  True,  True, True,  True,  True,  True,  True,  True,  True,  True,  True,  True,
  True,  True])

ax2.errorbar(dataQPU[inds,0], 2*dataQPU[inds,1],
              yerr=(dataQPU[inds,2],dataQPU[inds,3]),
              marker='o',
              markersize=7,
              markerfacecolor='#A0A0A0'+f'55',
              markeredgecolor='#A0A0A0',
              color='#A0A0A0',
              label='QA',
              linestyle='',
              **ebparams, zorder=0)

ax2.errorbar(dataQPU[ninds,0], 2*dataQPU[ninds,1],
              yerr=(dataQPU[ninds,2], dataQPU[ninds,3]),
              marker='o',
              markersize=7,
              color='#A0A0A0',
              fillstyle='none',
              label='_',
              linestyle='',
              **ebparams, zorder=0)

pp, ci, zx, zy = cfit.custom_fit(np.log10(dataQPU[inds,0]), np.log10(2*dataQPU[inds,1]))
ttest2 = np.arange(30, 2001)
ttest3 = np.arange(1, 30)
ax2.plot(ttest3, 10**(pp[0]*np.log10(ttest3)+pp[1]), linestyle='-', color='#A0A0A0')
ax2.plot(ttest2, 10**(pp[0]*np.log10(ttest2)+pp[1]), linestyle='--', color='#A0A0A0')


ax2.set_xlabel(f'Annealing time, $t_a$ (QA: ns, DT-SQA: MCS)')
ax2.set_ylabel(r'Residual energy, $\rho_E^f$')
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.tick_params(axis='both', which='both', direction='in')
ax2.minorticks_on()
ax2.legend(frameon=False, shadow=False)
ax2.set_ylim([2e-3, 1.0])
ax2.set_xlim([5, 5e3])


### the rightmost plot with modified extreme value prediction


MCS = np.array([50, 100, 200])
params = []

a = np.array([0.1182,  0.073777,  0.0477])
b = np.array([0.0141,  0.010156,  0.0086])
c = np.array([5, 16, 25])

reps = [50, 132, 250, 500,  1000, 2850, 6754,  16000]

ciHigh = np.array([-x[0][0] for x in ciAll]) - (-np.array(slopes))
ciLow = - np.array([-x[0][1] for x in ciAll]) + (-np.array(slopes))

ax3.errorbar(np.array(reps), -np.array(slopes), yerr=(ciLow, ciHigh), label='Measured slopes', linestyle='', linewidth=2.0, color=pltc.color[0], marker='o', markersize=6, markerfacecolor=pltc.color[0]+f'55',
         markeredgecolor=pltc.color[0], **ebparams, zorder=10)

ax3.plot(np.array(reps), 0.785*np.ones(len(reps)), label=f'Slope from QA', linestyle='-.', linewidth=1.5, color='#A0A0A0', zorder=0)

for ii in range(50, 20001):
    xxp = a + np.sqrt(2)* b * erfcinv(2 ** (1 - c / ii))
    ppx = np.polyfit(np.log10(MCS), np.log10(xxp), 1)
    params.append(ppx)

paramsList = np.array(params)
ax3.plot(np.arange(50, 20001), -paramsList[:, 0], label='Extreme value theory', linestyle='--', linewidth=2.0, color=pltc.color[2])

ax3.set_xlabel(r'Number of Trotter replicas, $R$')
ax3.set_ylabel(r'Slope, $\kappa_f$')
ax3.set_xscale('log')
ax3.set_yscale('log')

ax3.yaxis.set_major_formatter(ticker.ScalarFormatter())
ax3.yaxis.set_minor_formatter(ticker.ScalarFormatter())

formatter = ticker.ScalarFormatter(useOffset=False)
formatter.set_scientific(False)
formatter.set_powerlimits((-2,2))
ax3.yaxis.set_major_formatter(formatter)
ax3.yaxis.set_minor_formatter(formatter)

def format_func(value, tick_number):
    return f'{value:.1f}'  

ax3.yaxis.set_major_formatter(ticker.FuncFormatter(format_func))
ax3.yaxis.set_minor_formatter(ticker.FuncFormatter(format_func))

ax3.tick_params(axis='both', which='both', direction='in')
ax3.minorticks_on()
ax3.legend(frameon=False, shadow=False, loc='best')
ax3.set_xlim([40, 1e5])


fig.tight_layout()
#plt.savefig('Fig_2.pdf', dpi=300, bbox_inches='tight', transparent=False)
plt.show(block=True)