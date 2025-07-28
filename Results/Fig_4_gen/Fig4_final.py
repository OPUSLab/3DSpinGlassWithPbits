import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat
import matplotlib.font_manager as fm
import plot_config as pltc
import matplotlib.gridspec as gridspec

# Clear any existing plots
plt.close('all')
ebparams = {
    'capsize': 1.5,
    'elinewidth': 1.0,
    'capthick': 1.0,
}

# Load variables from .mat files
subplot1_vars = loadmat('subplot1_vars.mat')
subplot2_vars = loadmat('subplot2_vars.mat')
subplot3_vars = loadmat('subplot3_vars.mat')

# Extract variables from the loaded data
# Note: You'll need to verify the exact variable names in your .mat files
spin_sizes = subplot1_vars['spin_sizes'].flatten()
mean_min_seq_time_CPU = subplot1_vars['mean_min_seq_time_CPU'].flatten()
ci_min_seq_time_CPU = subplot1_vars['ci_min_seq_time_CPU']
cpu_time_per_sweep = subplot2_vars['cpu_time_per_sweep'].flatten()
ci_time_per_sweep = subplot2_vars['ci_time_per_sweep']
fpga_time_per_sweep = subplot2_vars['fpga_time_per_sweep'].item()
mtj_time_per_MCS_us = subplot2_vars['mtj_time_per_MCS_us'].item()
mean_total_time_CPU = subplot3_vars['mean_total_time_CPU'].flatten()
ci_total_time_CPU = subplot3_vars['ci_total_time_CPU']
mean_total_time_FPGA = subplot3_vars['mean_total_time_FPGA'].flatten()
ci_total_time_FPGA = subplot3_vars['ci_total_time_FPGA']
total_time_MTJ_upper = subplot3_vars['total_time_MTJ_upper'].flatten()
total_time_MTJ_lower = subplot3_vars['total_time_MTJ_lower'].flatten()

# Define plot parameters
color_cpu = pltc.color[1]
color_fpga = pltc.color[2]
color_mtj = pltc.color[4]
error_bar_thickness = 1.0
marker_size_cpu = 5
marker_size_fpga = 5

# Create a single figure with three subplots

fig = plt.figure(figsize=(8, 2.8), dpi=150)
gs = gridspec.GridSpec(1, 3, width_ratios=[1.2, 1.2, 2.0])  # Adjust these ratios as needed
ax1 = fig.add_subplot(gs[0])
ax2 = fig.add_subplot(gs[1])
ax3 = fig.add_subplot(gs[2])

# Subplot 1: Minimum MCS vs. spin size
ax1.errorbar(spin_sizes, mean_min_seq_time_CPU,
            yerr=[mean_min_seq_time_CPU - ci_min_seq_time_CPU[:, 0],
                  ci_min_seq_time_CPU[:, 1] - mean_min_seq_time_CPU],
            marker='o', color=color_cpu, markerfacecolor=color_cpu+f'55',
            markeredgecolor=color_cpu, linewidth=error_bar_thickness,
            markersize=marker_size_cpu, linestyle='', **ebparams)

ax1.set_ylabel(r'Minimum MCS to $ρ_\mathrm{E0}^\mathrm{f}$')
ax1.set_xlabel(f'Number of p-bits, $n$')
ax1.tick_params(which='both', direction='in', axis='both')
#ax1.legend(loc='best', frameon=False, shadow=False)
ax1.set_xscale('linear')
ax1.set_yscale('log')
ax1.set_xticks(spin_sizes)
ax1.set_xlim(300, 4300)
ax1.set_ylim(100, 3100)

# Subplot 2: Time per sweep vs. spin size
ax2.errorbar(spin_sizes, cpu_time_per_sweep,
            yerr=[cpu_time_per_sweep - ci_time_per_sweep[:, 0],
                  ci_time_per_sweep[:, 1] - cpu_time_per_sweep], label=f'CPU',
            marker='o', color=color_cpu, markerfacecolor=color_cpu+f'55',
            markeredgecolor=color_cpu, linewidth=error_bar_thickness,
            markersize=marker_size_cpu, linestyle='', **ebparams)

ax2.errorbar(spin_sizes, np.ones_like(spin_sizes) * fpga_time_per_sweep,
            yerr=0, label=f'FPGA', marker='s', color=color_fpga, markerfacecolor=color_fpga+f'55',
            markeredgecolor=color_fpga, linewidth=error_bar_thickness,
            markersize=marker_size_fpga, linestyle='', **ebparams)

ax2.plot(spin_sizes, np.ones_like(spin_sizes) * mtj_time_per_MCS_us,
         '--', label=f'sMTJ', color=color_mtj, linewidth=1.5)

ax2.set_ylabel(f'Measured MCS time per replica (μs)')
ax2.set_xlabel(f'Number of p-bits, $n$')
ax2.legend(frameon=False, shadow=False)
ax2.tick_params(which='both', direction='in', axis='both')
ax2.set_xscale('linear')
ax2.set_yscale('log')
ax2.set_xticks(spin_sizes)
ax2.set_xlim(300, 4300)
ax2.set_ylim(0.0001, 1.4e6)

# Subplot 3: MTJ bounds with fitted lines
ax3.errorbar(spin_sizes, mean_total_time_CPU,
            yerr=[mean_total_time_CPU - ci_total_time_CPU[:, 0],
                  ci_total_time_CPU[:, 1] - mean_total_time_CPU], label=f'CPU',
            marker='o', color=color_cpu, markerfacecolor=color_cpu+f'55',
            markeredgecolor=color_cpu, linewidth=error_bar_thickness,
            markersize=marker_size_cpu, linestyle='', **ebparams)


print(mean_total_time_FPGA)

ax3.errorbar(spin_sizes, mean_total_time_FPGA,
            yerr=[mean_total_time_FPGA - ci_total_time_FPGA[:, 0],
                  ci_total_time_FPGA[:, 1] - mean_total_time_FPGA], label=f'FPGA',
             marker='s', color=color_fpga, markerfacecolor=color_fpga+f'55',
            markeredgecolor=color_fpga, linewidth=error_bar_thickness,
            markersize=marker_size_fpga, linestyle='', **ebparams)

# Fit lines to MTJ bounds in log scale
log_spin_sizes = np.log10(spin_sizes)
log_total_time_MTJ_upper = np.log10(total_time_MTJ_upper)
log_total_time_MTJ_lower = np.log10(total_time_MTJ_lower)

# Linear fit in log-log space
p_upper = np.polyfit(log_spin_sizes, log_total_time_MTJ_upper, 1)
p_lower = np.polyfit(log_spin_sizes, log_total_time_MTJ_lower, 1)

# Evaluate fitted lines
fitted_MTJ_upper = 10 ** np.polyval(p_upper, log_spin_sizes)
fitted_MTJ_lower = 10 ** np.polyval(p_lower, log_spin_sizes)

ax3.set_ylabel(f'Total estimated time (μs)')
ax3.set_xlabel(f'Number of p-bits, $n$')
ax3.legend( loc='best', frameon=False, shadow=False)
ax3.tick_params(which='both', direction='in', axis='both')
ax3.set_xscale('linear')
ax3.set_yscale('log')
ax3.set_xticks(spin_sizes)
ax3.set_xlim(300, 4300)
ax3.set_ylim(0.0001, 5e9)

plt.tight_layout()
# plt.savefig('Fig4_new.pdf', dpi=300, bbox_inches='tight', transparent=False)
plt.show()