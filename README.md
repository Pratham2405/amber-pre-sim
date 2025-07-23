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
##### `.dat`
Main parameter file which stores all the numerical constants, pparameters categorised into MASS, BOND, ANGLE, DIHEDRAL(torsional) and NONBOND sections.
[Image of .dat file screenshot]
##### `.frcmod`
Same format as `.dat` but for modifying parameters, especially used for non-standard residues, ligands, metal ions etc.
> Parameters written in `.frcmod` overwrite the ones written in `.dat` files.
##### `.lib`
Stores residue topology template in OFF(Object File Format - [What is OFF?]) and carries the following information:
1. Connection information for linking residues.
2. Partial charge assignment for each atom.
3. Complete residue definition with atom names, types and connectivity within residue.
[Image of .lib file sreenshot]

#### Different types of Force Fields in AMBER
[Information from the AMBER manual about force fields(brief-one or two lines)]
1. ff19sb
2. ff14sb
3. gaff
4. olaf
### Important functions in AMBER
#### LEaP
Our standard structure files like `.pdb` or `.mol2` cannot be used directly by AMBER since they do not contain crucial data required for implementing an MD simulation. We have already discussed that an MD simulation requires force field files because that is what leads to force calculations - the 'dynamic' end of any MD simulation. However, even with this force field data at our disposal, we need a function to map those parameters to each atom of the protein/ligand of interest, enabling the force-computing functions(sander/pmemd) downstream with the right set of force-field parameters, connectivity(both bonded and non-bonded) and partial charges, such that accurate simulations become a natural result: LEaP is that function of AMBER. A very intuitive tutorial by Pengfei Li and David Cerutti[https://ambermd.org/tutorials/pengfei/index.php] might be of help to understand the fundamentals of LEaP in more detail.
[Image from the tutorial]
In summary, LEaP is the central preparatory program in AMBER which collates input structure files(`.pdb` and `.mol2`) and force field files(`.dat`, `.lib` and `.frcmod`) to make a complete description of molecules ready for dynamic calculations.
It has the following output files:
1. Topology files(`.prmtop`)
  - Contain complete molecular connectivity information(bonded and non-bonded)
  - All force-field parameters(bonds, angles, dihedrals, non-bonded terms)
  - Atom Types and Charges
2. Coordinate Files(`.inpcrd`,`.rst7`)
  - Atomic Coordinates
  - Box dimensions(for periodic boundary conditions, `ntb=1`)

#### sander/pmemd
These are the functions which carry out the simulations and the dynamic calculations. There are many kinds of simulations you can carry out with sander:
1. Minimisation
When water molecules and ions are added to the system, they are simply arranged randomly without accounting for any interactions and minimum distance between molecules; this needs to be corrected such that the solvent and ions acquire a relaxed and appropriately-spced configuration. Moreover, you are bound to receive many warnings from the `leap.log` file about gaps and improper contacts within the protein; minimisation also allows the protein to correct these improper contacts.
2. Heating
[Para about Heating]
3. Equilibriation
[Para about Equilibriation]
4. Production 
[Para about Production]

#### `reduce`
[Para about reduce]

#### `pdb4amber`
[Para about pdb4amber]

#### `ambpdb`
[Para about ambpdb]

### That which we won't cover in this repository
It is important to know the bounds of our knowledge when starting out, and to that effect, the following will not be discussed here and the reader is referred to primary sources about them:
1. Non-standard residues encountered in protein
2. Non-standard ligands covered in the [Corrections required]

### Overall Workflow
[Image of the flowchart]

### Steps of the Workflow
#### Protein Preparation
#### Ligand Preparation

#### LEaP
[Input Files Required -> Explanation of Arguments -> Output File and what they carry]
#### Minimisation
#### Heating(Thermalisation)
#### Equilibriation
#### Production
