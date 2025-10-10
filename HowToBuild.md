# Build Instructions

## Follow these steps to build n2p2:

1. Navigate to the interface directory:
   ```
   cd <n2p2_directory>/src/interface/
   ```

2. Modify the `makefile` at line 58 to enable the necessary LAMMPS packages:
   ```
   cd lammps-hdnnp/src/ && $(MAKE) yes-ml-hdnnp && $(MAKE) yes-kspace && $(MAKE) yes-molecule && $(MAKE) yes-manybody && $(MAKE) yes-misc && $(MAKE) yes-rigid && $(MAKE) yes-qeq && $(MAKE) yes-spin && $(MAKE) yes-misc && $(MAKE) yes-tally && $(MAKE) yes-phonon && $(MAKE) mpi
   ```

3. Navigate to the source directory and build the LAMMPS interface:
   ```
   cd <n2p2_directory>/src/
   export LD_LIBRARY_PATH=~/anaconda3/lib:$LD_LIBRARY_PATH
   make clean && make -j8 lammps-hdnnp
   ```
   > **Note:** If fails, fix the issue and run `make clean-lammps-hdnnp && make -j8 lammps-hdnnp`

4. Install required libraries via conda:
   ```
   conda install eigen gsl
   ```

5. Update the `makefile.gnu` configuration:
   - Set line 8: `PROJECT_GSL=<conda_environment_path>/include/`
   - Set line 9: `PROJECT_EIGEN=<conda_environment_path>/include/eigen3/`
   - Set line 24: `PROJECT_LDFLAGS_BLAS=-L<conda_environment_path>/lib -lopenblas -lgsl -lgslcblas`

6. Build all targets:
   ```
   make -j8 all
   ```

## Optional: Control charge feature for Q-HDNNP and 4G-HDNNP

You can control whether the predicted charge is included as a feature in the energy neural network for Q-HDNNP and 4G-HDNNP types using the following flags:

- `INCLUDE_CHARGE_IN_ENERGY_Q=1` (default: enabled)
- `INCLUDE_CHARGE_IN_ENERGY_4G=1` (default: enabled)

To disable charge as a feature for either type, set the flag to 0 when building, e.g.:

```
make clean && make -j8 lammps-hdnnp INCLUDE_CHARGE_IN_ENERGY_Q=0 INCLUDE_CHARGE_IN_ENERGY_4G=1
```

This allows you to independently control charge inclusion for Q-HDNNP and 4G-HDNNP builds.

## Alternative method with system packages

If you have sudo permissions, you can install system packages instead of using conda:

```
sudo apt install -y libopenmpi-dev libgsl-dev libeigen3-dev libblas-dev liblapack-dev
```

Replace `<n2p2_directory>` with the actual path to your n2p2 installation and `<conda_environment_path>` with the path to your conda environment (e.g., `~/anaconda3`).