# MACE_ASE-NEB_tests
An example for testing and analysing different implementations of NEB with the MACE ML force-field.

Running script found in either neb_tester.py (with HAWK GPU submission script) or on your PC with Testrunner.ipynb. Paths to the working directory will need to be changed manually unfortunately (lazy coding on my part).

If you only care about the results, check out Testanalysis.ipynb!

# Summary

Comparisons are performed for the default ASE-NEB implementation and the ODE-based string/spline NEBOptimizer. NEB calculations for the latter have been performed with and without the Exp geometry optimization pre-conditioner. Either 5 or 6 interpolating images were used.

Tests were performed for reactions: occurring on the Cu100 and Cu110 surface (HCO2 formation (Reaction 1) and CO2 dissociation respectively (Reaction 2)); and the migration of a CH3 group from one oxygen site to another on a H-SSZ-13 zeolite (Reaction 3).

From the tests performed, the spline method combined with the Exp preconditioner gave the most reliable convergence of the tested methods, successfully converging below the (rather tight) Fmax threshold of 0.02 Ha/Bohr for all test examples. The string method was far less robust, failing to approach a looser convergence threshold of 0.04 Ha/Bohr. Convergence was achieved for the spline NEB Optimizer case within 161, 84, and 118 total steps for reactions 1-3 respectively. In the case where the default NEB implementation managed to converge, the preconditioned spline optimizer gave a speed-up of 8 steps overall.

In summary - the spline method with preconditioning converged more robustly than the other NEB methods within a reasonable number of steps. If you want to have a closer look at the data, some run information is stored in JSON files, which can be read and plotted in the Testanalysis Jupyter notebook in the git repository above. The running scripts in there should explain how to set-up these calculations as well.

![alt text](https://github.com/GabrielBram/MACE_ASE-NEB_tests/blob/main/figures/Cu100_CO2_dissoc_NEB-comparison.png?raw=true)
![alt text](https://github.com/GabrielBram/MACE_ASE-NEB_tests/blob/main/figures/Cu110_HCO2_formation_NEB-comparison.png?raw=true)
![alt text](https://github.com/GabrielBram/MACE_ASE-NEB_tests/blob/main/figures/H-SSZ-13_SMS_migration_NEB-comparison.png?raw=true)
