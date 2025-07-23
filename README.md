# AMBER Pre-Simulation Protocol
Preparation protocol and production scripts for simulating Protein-Ligand Complex in AMBER - a tutorial with a brief introduction to MD Simulations concepts and typical workflow in AMBER.

## Introduction to AMBER & MD Simulations
[Para on Purpose of MDSim]
### Force Fields: What exactly are they?
A force field in an MD Simulation software defines the characteristic potential energy of conventional/unconventional atoms and their bonds via equations and paramters which it stores in text files(`.dat` for parameters, `.frcmod` for modified parameters and `.lib` for topology template). In doing so, the force on each atom can be calculated by the gradient(in the direction of maximum change) of the potential energy. 
A Force Field therefore refers to the functional form and the parameter sets used to calculate the potential energy of a system on an atomistic level.
Some usual components of the force fields are as follows:
#### Bond Stretching
#### Angle 
#### Torsional Rotation
#### Non-bonded Interactions
[The parameter sets for different atoms are determined by various experimental techniques to ensure accurate predictions which are closer to reality - Raman spectroscopy...]
### Approximations in MD Simulations:
- Classical Mechanics is considered where atoms are approximated to point particles and newtonian mechanics is only considered.
- Only pairwise interactions are considered.
- Simple harmonic motion for bonds and angles considered.
- Fixed partial charges instead of assigning them dynamically.

### Structure of Force Field Files
The force field is stored in `.dat`, `.frcmod` and `.lib` files:
#### `.dat`
Main parameter file which stores all the numerical constants, pparameters categorised into MASS, BOND, ANGLE, DIHEDTRAL(torsional) and NONBOND sections.
#### `.frcmod`
Same format as `.dat` but for modifying parameters, especially used for non-standard residues, ligands, metal ions etc.
> Parameters written in `.frcmod` overwrite the ones written in `.dat` files.
#### `.lib`
Stores residue topology template in OFF(Object File Format - [What is OFF?]) and carries the following information:
1. Connection information for linking residues.
2. Partial charge assignment for each atom.
3. Complete residue definition with atom names, types and connectivity within residue.
