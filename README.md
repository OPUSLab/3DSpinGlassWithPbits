# 3DSpinGlassWithPbits

This repository contains code and associated data used to generate the results of the paper titled:  
**“Pushing the Boundary of Quantum Advantage in Hard Combinatorial Optimization with Probabilistic Computers.”**

## ⚙️ Requirements

- MATLAB (developed and tested on recent versions)

## ⚙️ Usage

1. Clone this repository:
   ```
   https://github.com/OPUSLab/3DSpinGlassWithPbits.git
   ```
2. Open MATLAB and navigate to the cloned directory.

3. Run:
   ```matlab
   APT.m
   ``` 
   to use the adaptive parallel tempering (APT) or,

   ```matlab
   SQA.m
   ```
   to use simulated quantum annealing.

4. The codes will save the lowest energy among all the replicas at the end of annealing in a .mat file.


## ⚙️ Customization

You can modify the following parameters in the scripts (APT.m or SQA.m):

- `logical`: 1 for logical instances and 0 for embedded instances
- `instanceSize`: Size of the instances (8, 10, 12, 15 and 16 can be used for logical instances, 15 can be used for embedded instances)
- `instanceID`: Which instance to run (0 to 299)
- `run_per_instance`: Number of independent experiments to run per instance (50 is used in this work)
- `total_real_num_sweeps`: Total number of sweeps (MCS) to be used
- `num_sweeps_per_swap`: Number of sweeps to use before making a swap attempt (only used in APT.m)
- `num_sweeps_read_per_swap`: Number of sweeps to use in a swap attempt (1 is used in this work, only used in APT.m)
- `num_internal_replicas` : Number of ICM replicas to be used per temperature (only used in APT.m)
- `base_seed`: Seed to the random number generator
- `num_replicas`: Number of Trotter replicas to be used (only used in SQA.m)
- `GammaX`: Magnitude of the transverse field (only used in SQA.m)
- `betaAll`: Inverse temperature to number of replicas ratio (only used in SQA.m)


## ⚙️ Plotting

The repository also contains processed data and python codes to generate Figures 2, 3 and 4 of the main text.

## Contributing

Contributions to improve the code or extend its functionality are welcome. Please feel free to submit issues or pull requests.


## Citations

To cite this work, please cite the following paper: 
Chowdhury, S., Aadit, N.A., Grimaldi, A. et al. Pushing the boundary of quantum advantage in hard combinatorial optimization with probabilistic computers. Nat Commun 16, 9193 (2025). https://doi.org/10.1038/s41467-025-64235-y

## Contact

If you have any questions or suggestions, please open an issue in this repository or contact Shuvro Chowdhury (schowdhury@ucsb.edu).
