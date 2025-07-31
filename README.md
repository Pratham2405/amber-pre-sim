# AMBER Pre-Simulation Protocol
Preparation protocol and production scripts for simulating Protein-Ligand Complex in AMBER - a tutorial with an in-depth introduction to MD Simulations concepts and typical workflow in AMBER.

## Introduction to AMBER & MD Simulations
### Objectives of an MD Simulation
A Molecular Dynamics Simulation is an attempt to solve nature on our computers by approximating how microscopic interactions work and how they affect the dynamics of atoms in a molecular system. Although MD Simulations, like any other solution, is also wrought with problems and challenges of its own, it is still the most predictive tool we have short of _in vitro_ experimentation.

### Force Fields: What exactly are they?
A force field in an MD Simulation software defines the characteristic potential energy of conventional/unconventional atoms and their bonds via equations and paramters which it stores in text files(`.dat` for parameters, `.frcmod` for modified parameters and `.lib` for topology template). In doing so, the force on each atom can be calculated by the gradient(in the direction of maximum change) of the potential energy. 
A Force Field therefore refers to the functional form and the parameter sets used to calculate the potential energy of a system on an atomistic level.
Some usual components of the force fields are as follows:

#### Bond Stretching  
$E_{\text{bond}} = \sum_{\text{bonds}} K_{b}\(r - r_{0})^{2}$

#### Angle  
$E_{\text{angle}} = \sum_{\text{angles}} K_{\theta}\(\theta - \theta_{0})^{2}$

#### Torsional Rotation  
$E_{\text{torsion}} = \sum_{\text{torsions}} \frac{V_{n}}{2}\\left[1 + \cos(n\phi - \gamma)\right]$

#### Non-bonded Interactions  
$E_{\text{nonbond}} = \sum_{i<j} \left[ \epsilon_{ij} \left( \left(\frac{R_{\min,ij}}{r_{ij}}\right)^{12} - 2 \left(\frac{R_{\min,ij}}{r_{ij}}\right)^6 \right) + \frac{q_i q_j}{4 \pi \epsilon_0 r_{ij}} \right]$


The parameter sets for different atoms are determined by various experimental techniques to ensure accurate predictions which are closer to reality - Raman spectroscopy, infrared (IR) spectroscopy, X-ray crystallography, nuclear magnetic resonance (NMR) spectroscopy and other methods.

### Approximations in MD Simulations:
- Classical Mechanics is considered where atoms are approximated to point particles and newtonian mechanics is only considered.
- Only pairwise interactions are considered.
- Simple harmonic motion for bonds and angles considered.
- Fixed partial charges instead of assigning them dynamically.

### Minimisation
Incremental movement of the atoms under the influence of the gradient of Potential Energy(PE) _in vacuo_ in search of a global minima of PE.
#### Rationale behind Minimisation
There are several reasons why it makes sense to minimise your system:
1. PDBs have structural artifacts like atomic clashes, steric overlap and unrealistic bond lengths which may lead to unrealistic potential energies and might crash the simulation.
2. Steric clashes of non-bonded atoms can lead to unusually high Van der Waals repulsion. 
3. Provides a reasonable starting value of PE for later MD runs.
#### Implementation in AMBER
Generally, minimisation is carried out in 2 steps:
1. **Restrained Minimisation(`ntr=1`)**: It is advisable to restrain the protein-ligand complex for the first minimisation run. This is useful because we don't need the protein to undergo drastic changes in its backbone structure due to the effects of the solvent first because the solvent itself has not been minimised and just packed into the solvent box, facing steric clashes among itself. All of this can lead to drastic minimisation steps which might move the protein atoms more than is required for minimisation.
AMBER applies a harmonic restraint term to the restrained atoms in the form of a quadratic penalty term which penalises deviation of atoms from their initial coordinates($r_0$).
$E_{\text{restraint}} = \sum_{\text{i}} (K/2)\(r_{i} - r_{i,0})^{2}$
2 sander parameters are crucial when carrying out restrained minimisation(`ntr=1`):
- `restraint_wt`: It is equivalent to `K` in the above equation and determines the strength of penalty. Described in kcal/molÅ^2 and ranges from 5-10(weak) to 20-50(moderate) and strong(100-500).
- `restraintmask`: Chooses the atoms you want to restrain.
2. **Unrestrained Minimisation(`ntr=2`)**: Remove all restraints and minimise the entire system. Alternatively, you could stage the restrained minimisation in steps, slowly decreasing the restraint weight before going for the unrestrained minimisation.

#### Differences from Docking
Although Minimisation and Docking are both minima search/optimisation problems which happen when the atoms have 0 kinetic energy, they differ in their approach and their final objectives:
- Minimisation aims to decrease the potential energy of the entire system(protein + ligand + solvent), whereas docking doesn't account for the solvent and assumes the receptor as a mostly rigid entity.
- Minimisation optimises the PE as determined by the force field and docking focusses on other scoring functions which prioritises speed over realistic behavior. It predicts a binding affinity for finding the best ligand pose.
- Minimisation makes small increments/displacements to original coordinated to find a minima in the topographical vicinity, thereby assuming that the structure is already close to one. Docking is more useful when trying to find the best binding pose.
> It is for this very reason that we would be using docked poses as our starting coordinates for our MD simulation.

### Thermalisation(Heating)
#### Rationale behind Thermalisation
Minimisation achieves a low PE state but the kinetic energy of atoms is 0. This contradicts with the kinetic energy(KE) at room temperature(300K). Hence, in order to carry out simulation we need to bring the molecules at 300K(essentially giving them a velicty) slowly from 0K.
This increase of temperature is achieved by the use of thermostats like Berendsen(`ntt=1`), Andersen(`ntt=2`), Langevin(`ntt=3`) and Nose-Hoover(`ntt=4`). Langevin thermostat is the most widely used and that is the one we would be using too.
When a system of molecules is subjected to a heat bath, it turns energy in the form of KE to molecules of the system. To impart KE in computational terms is like trying to emulate this phenomena. Langevin thermostat does this by choosing a subset of atoms randomly and gives them a 'push'(Random force simulating transfer of KE from heat bath) and a dissipation force(similar to friction/damping) which balances the the kicks, thus leading to temperature stabilisation eventually.
Following equation defines the Langevin Dynamics:

$m \frac{d^2 \mathbf{r}}{dt^2} = - \gamma \frac{d\mathbf{r}}{dt} + \mathbf{F}(\mathbf{r}) + \mathbf{R}(t)$

#### Implementation in AMBER(Fluctuation-Dissipation Theorem)
The thermal kick has been defined as follows:

$\langle R_{\alpha}(t)\, R_{\beta}(t') \rangle = 2\gamma k_B T\ \delta_{\alpha\beta} \delta(t-t')$

The Kronecker-Delta term tells and the Dirac Delta function imply that the expected/average value of R(t) has no relation with [Equation here]. Which, in essence, means that the thermal impulses/pushes/kicks are completely random, emulating molecular collisions.
Moreover, the Random force increases with T, $m_i$(heavier atoms need stronger pushes) and $\gamma_{i}$, the collision frequency, which decides the frequency of introduction of the kicks. It is also called coupling constant, since it also measures how affected/coupled is the system with the thermodynamic bath. 
It is advisable to heat the system with a restraint on the protein-ligand complex to prevent from the protein changing its structure too much.

#### Relevant sander parameters
- **`ntt`**: Thermostat selection.
- **gamma_ln**: Collision Frequency(default = 2). Lower values take more time to reach T0 but give more realistic behaviour.

### Equilibration

#### Rationale behind Equilibration
The step of equilibriation is to make the system converge to a stable temperature nad pressure suitable for production runs. Equilibriation stabilises the temperature of the system to a target temperature via a thermostat with the same initial and final temperature. The solvent molecules still be cramped up in the box unevenly. A pressure equilibriation allows for the system to expand the box under constant pressure to execute a density equilibration.
#### Implementation in AMBER
It is generally carried out by an initial NVT(using a thermostat just like in thermalisation but `tempi` and `temp0` are equal) followed by a series of NPT runs with decreasing `restraint_wt`. NPT is carried out using a thermostat and a barostat(Pressure coupling) together like Monte Carlo, Parinello-Rahman and Berendsen. These work by adjusting the volume of the system under control pressure. The system's instantaneous pressure is calculated by the virial equation:

$pV = Nk_BT + \frac{1}{3} \left\langle \sum_{i=1}^N \mathbf{r}_i \cdot \mathbf{F}_i \right\rangle$

All barostats use this fundamnental equation to convert microscopic forces and system variables to compute macroscopic pressure of the system. The barostats differ by how much complexity(tensor pressure accounts for anisotropic stresses like shear stress on the box) they allow in pressure consideration, how often they compute the virial pressure, how they respond to pressure fluctuations from the target pressure and how they introduce changes to box volume and atom coordinates. If you have had a long, stable NVT in your heating step, it's admissable to skip the NVT equilibriation and go for NPT equilibriation with decreasing `restraint_wt`.
#### Relevant sander parameters
- **`ntb`**: 
  - `=1`: constant volume condition(NVT).
  - `=2`: constant pressure condition(NPT); allows the box to change the dimensions
- **`ntt`**: Thermostat selection
  - `=1`: Berendsen
  - `=2`: Andersen
  - `=3`: Langevin
  - `=4`: Nose-Hoover
- **`ntp`**: Barostat selection
  - `=0`: No pressure coupling(NVT)
  - `=1`: Berendsen 
  - `=2`: Monte-Carlo 
  - `=3`: Parinello-Rahman

### Production
The production run is when the actual 'simulation' takes place and Newton's equations of motion are solved for the force-field along with verlet increments in coordinates and velocity.
### File Types in AMBER
One approach to know about the various types of files with overlapping information is to know them by what they do:
#### Defining Molecular Structure and Composition  
**Extensions:** `.pdb`, `.mol2`, `.prepi`, `.lib` / `.off`  
These files define atoms, residues, connectivity, coordinates, and charges for proteins, ligands, and non-standard residues. They are loaded by LEaP during system setup to build the initial molecular model.

#### Providing Force Field Parameters  
**Extensions:** `.frcmod`, `.prmtop` (also `.top`, `.parm7`), `parm10.dat` (internal)   
They contain bonded and nonbonded interaction parameters (bonds, angles, dihedrals, charges, atom types). FRCMOD supplements missing parameters; PRMTOP combines the full topology and force field parameters for MD engines.

#### Supplying Initial Coordinates and Velocities  
**Extensions:** `.inpcrd`, `.rst`, `.crd` (legacy), `.rst7`   
These files store starting atomic positions, optional velocities, and box dimensions. They serve as input for launching or restarting molecular dynamics simulations.

#### Recording Simulation Trajectories  
**Extensions:** `.mdcrd`, `.nc` 
Trajectory files record atomic coordinates frame-by-frame. MDCRD is ASCII; NC is compressed NetCDF binary developed for efficiency. They capture system evolution for visualization and analysis.

#### Controlling Simulation Parameters  
**Extension:** `.mdin`  
This control file specifies MD protocol parameters (time step, temperature coupling, restraints, output frequency) and is read by `sander`/`pmemd` to run simulations.

#### Storing Analysis and Miscellaneous Data  
**Extensions:** `.dat`, `.out`, `.log` 
General-purpose data, log, and analysis output files. Most `.dat` files store analysis results, but files like `parm10.dat` carry core force field tables loaded by LEaP internally.


### Important Tools in AMBER
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
As we have already covered most of the pivotal sander arguments, here are some miscellaneous arguments which are also used in sander/pmemd:
- **ntx:**  
  Specifies the format of the input coordinate file: `1` means only coordinates, `5` means coordinates plus velocities.  
  Used by `sander` to interpret whether velocities are read from the restart file for continued simulations.

- **ntc:**  
  Determines how bond constraints are applied: `1` constrains bonds involving hydrogen, `2` constrains all bonds.  
  Frequently set to `2` to allow a larger integration time step by fixing fast bond vibrations.

- **ntf:**  
  Specifies which forces are excluded for constrained bonds: `1` excludes none, `2` excludes forces on bonds constrained by `ntc=2`.  
  Ensures consistency between constraint application (`ntc`) and force calculation to maintain energy conservation.

- **nstlim:**  
  Number of MD steps to perform in the simulation.  
  Combined with `dt` to define the total simulated time (`nstlim × dt`).

- **dt:**  
  Time step in picoseconds for each MD integration step (e.g., `0.002` for 2 fs).  
  Chosen based on system stability and constraint settings to ensure accurate energy conservation.

- **irest:**  
  If `0`, start a new simulation; if `1`, restart from existing velocities and positions in the input coordinate file.  
  Controls whether `sander` reads velocities (`irest=1`) or initializes them (`irest=0`).

- **ntpr** & **ntwx:**  
  `ntpr` sets how frequently (in steps) energy and temperature data are printed to the output file.  
  `ntwx` sets how frequently (in steps) coordinates are written to the trajectory file.

- **cut:**  
  Defines the nonbonded interaction cutoff distance in angstroms.  
  Only interactions within `cut` are computed explicitly, improving computational efficiency.  

#### `pdb4amber`
Prepares the protein pdb file for LEaP and can also add missing atoms and fixing residue numbers(31-150 becomes 1-120).

#### `ambpdb`
Converts AMBER-compatible topology files(`.prmtop`) and coordinate files(`.crd`, `.inpcrd`, `.rst7`) into PDB format. Helpful when viewing final views of trajectory runs on PyMol.

#### `antechamber`
Automatically assigns GAFF/GAFF2 atom types, bond types, and AM1-BCC charges to small molecules, generating MOL2 topology files.

#### `parmchk2`
Analyzes a MOL2 file’s atom types and bond connectivity to identify missing GAFF/GAFF2 parameters, then generates an `.frcmod` file with estimated bond, angle, dihedral, and nonbonded parameters. It ensures non-standard residues receive reasonable force field parameter values before LEaP integration 


### Overall Workflow
[Image of the flowchart]

### That which we won't cover in this repository
It is important to know the bounds of our knowledge when starting out, and to that effect, the following will not be discussed here and the reader is referred to more primary sources:
1. Non-standard residues encountered in protein. These would require custom preparation using parmchk2 and antechamber.
2. Systems with more than one ligand, multiple bound small molecules, or complex environments (e.g., cofactors, substrates, and ions together) are not addressed.
3. We have used explicit solvent in this protocol. However, implicit solvent(`igb=1`) also gives accurate results if you want to reduce your computation load.
4. Post simulation analysis is left for some later repository.

### Steps of the Workflow
#### Protein Preparation
The only data needed in a PDB file to set up AMBER simulations are atom names, coordinates of heavy atoms and chain identifier(in case there is more than one chain in the protein). If there is any extra feature except the protein like important water molecules near the binding site or any other moiety would require you to load extra library files.
The challenge in protein preparation arises because by default, PDB files carry information irrelevant to AMBER like ANISOU, REMARK, HEADER lines and CONECT records(tleap implies connection from atom distances). Moreover, PDBs suffer with issues like missing hydrogen atoms(due to X-Ray crystallography limitations), steric clashes among atoms of protein backbone, lack of protonation states of residues and insonsistent naming convention for residues. These and other irrelevant information needs to be removed from the PDB file. Due to limitations in experimental structure determination there can be many structural artifacts in the PDB file like gaps and the protein of interest might not even have an experimental structure to begin with. To circumvent this issue, we can take the assistance of deep-learning structure prediction tools like [I-TASSER], [Swissmodel], [MODELER] etc. We will be using I-TASSER for our tutorial to fill in the gaps in our target protein. I-TASSER may take upto 2 days to generate the predicted structure and you are advised to choose the structure with the highest C-value. 
- Once you have obtained the predicted structure, inspect the PDB for any unnecessary lines. Generally, the PDB file would have no irrelevant lines such as ANISOU, CONECT and no HOH atoms or H atoms either.
- Take the PDB structure and give it as input to [H++]. Choose the pH as per your use case(physiological pH = 7.4). The resulting output would contain many files but we would be needing the .crd and the .top files for our purposes.
- Use the following command to use `ambpdb` to convert the topology file and the coordinate files into PDB format to be used in LEaP:
```
ambpdb -p hpp_output.top -c hpp_output.crd > hpp_output_2WZX.pdb
```
The output pdb file would have different protonation states for each residue, especially histidine, H would be added appropriately.
- After we know the correct protonation states for our residues at working pH, we can use `pdb4amber` for our final protein preparation by the following command:
```
pdb4amber -i hpp_output_2WZX.pdb -o 2WZX_prepared.pdb
```
- Inspect both the pdb file and its contents and also view the output using visualisation tools like PyMol or VMD and look for any abnormalities.

#### Ligand Preparation
We have discussed in ([Minimisation](#minimisation)) how docking output poses can help us give us a very good starting point for initial coordinates for our ligand w.r.t the protein. This implies that we would be using the highest docked pose from the `output.pdbqt` file of our docking run we carried out in [Vina Docking Protocol](https://github.com/Pratham2405/Vina_Docking_Protocol). 
- With the `output.pdbqt`, simply delete everything below the first `ENDMDL` line.
- Before proceeding further, we need to calculate `-nc` of the ligand for the coming antechamber step:
```
grep -E '^(ATOM|HETATM)' 0.160_out.pdbqt | awk '{sum += $(NF-1)} END {print sum}'
```
Note down this charge to the nearest integer for future use.
- Convert to `.sdf`  and then to `.pdb`using basic Openbabel commands:
```
obabel -ipdbqt 0.160_out.pdbqt -osdf -O 0.160_MD.sdf
obabel -isdf 0.160_MD.sdf -opdb -O 0.160_MD.pdb
```
- The file now has irrelevant lines that need to be filtered out to a new line, moreover, the `UNL` tag has which is often found in ligand pdb files has to be changed to `LIG` for compatibility with AMBER conventions:
```
sed 's/UNL/LIG/g' 0.160_MD.pdb | grep -E '^HETATM' > 0.160_MD_clean.pdb
```
This command renames the `UNL` lines to `LIG`, pipes the output to next command, which then filters all the `HETATM` lines and copy them to a new pdb file. You might want to preserve the `CONECT` records for inferred connectivity. In that case, you might want to replace `HETATM` with `(HETATM|CONECT)`. However, from experience, antechamber is quite good at inferring connectivity provided all the H are provided.
- Obabel, antechamber and reduce are all tools capable of adding hydrogens to a ligand. However, from experience, these tools can't add H even if there is one other H in the ligand already. To overcome this challenge, we can remove all H atoms from the ligand and then add them again(I prefer Obabel for addition at required pH):
```
obabel -ipdb 0.160_MD_clean.pdb -opdb -O 0.160_MD_noH.pdb -d
obabel 0.160_MD_clean.pdb -O 0.160_MD_H.pdb -p 7.4
```
- Now we use antechamber to prepare the parameter file for the ligand, which uses the GAFF force field to assign atom types:
```
antechamber -i 0.160_MD_H.pdb -fi pdb -o 0.160_MD.mol2 -fo mol2 -c bcc -nc -2 -m 1 -at gaff2 -rn LIG -s 2 -j 4
```
- `parmchk2` identifies any atom types and parameters(bonds, angles, dihedral and non-bonded) not identified by GAFF2 and generates the `.frcmod` file with estimated values for those atom types and parameters:
```
parmchk2 -i 0.160_MD.mol2 -f mol2 -o 0.160_MD.frcmod -at gaff2
```
- In conjunction with a `.frcmod` file, we need a `.prepi` file for the coordinates, connectivities and atom information like atom charges to incorporate the non-standard residue in the simulation:
```
antechamber -i 0.160_MD.mol2 -fi mol2 -o 0.160_MD.prepi -fo prepi
```

#### LEaP
We have discussed in length about the purpose and function of `tleap`(alias for LEaP) in ([LEaP](#leap)). 
Input files required by tleap: 
- Force Field Parameter Files(`.dat` and `.frcmod`)
- Residue Topology files(`.lib`)
- Prepared Structure Files(`.pdb` and `.mol2`)
You can view `leap.in` to view the arguments for the tleap program which are pretty self-explanatory. 
The file is run by the following command:
```
tleap -f leap.in > leap.log
```
Output files provided:
- `.prmtop`(Topology File): Contains molecular connectivity information, all force field parameters, atom types and charges
- `.inpcrd`(Coordinate File): Contains atomic coordinates and box dimension(for periodic boundary condition). 
- `leap.log`(Log File): Carries a detailed log of the tleap run. `leap.log` gives out many warnings of steric clashes, while these are normal, any error encountered could fatally kill tleap execution and `leap.log` helps us pinpoint the exact cause for the fatal error in the input files.

> It is important to view `system.prmtop` on PyMol to check for any structural inconsistencies.

#### Minimisation
Minimisation is implemented in two stages - restrained and unrestrained - as discussed in ([Minimisation](#minimisation)). The input `.in` files for both can be found in this repository.
The command for both the minimisations using `pmemd.MPI` are as follows(the syntax for sander remains exactly the same save for starting with `sander`):
```
pmemd.MPI -O -i min_restrained.in -o min1.out -p system.prmtop -c system.inpcrd -r min1.rst -x min1.nc -inf min1.mdinfo -ref system.inpcrd
memd.MPI -O -i min_unrestrained.in -o min2.out -p system.prmtop -c min1.rst -r min2.rst -x min2.nc -inf min2.mdinfo -ref min1.rst
```
> `-c`, `-p` and `-i` are the input arguments and the rest are output. This stays true for all the subsequent operations of sander and pmemd.

#### Post-Minimisation
The commands for Thermalisation, Equilibriation and Production runs are very similar to minimisation:
```
pmemd.MPI -O -i heat.in -o heat.out -p system.prmtop -c min2.rst -r heat.rst -x heat.nc -inf mdinfo -ref min2.rst
pmemd.MPI -O -i eq1.in -o eq1.out -p system.prmtop -c heat.rst -r eq1.rst -x eq1.nc -inf mdinfo -ref heat.rst
pmemd.MPI -O -i eq2.in -o eq2.out -p system.prmtop -c eq1.rst -r eq2.rst -x eq2.nc -inf mdinfo -ref eq1.rst
pmemd.MPI -O -i eq3.in -o eq3.out -p system.prmtop -c eq2.rst -r eq3.rst -x eq3.nc -inf mdinfo -ref eq2.rst
pmemd.MPI -O -i eq4.in -o eq4.out -p system.prmtop -c eq3.rst -r eq4.rst -x eq4.nc -inf mdinfo -ref eq3.rst
pmemd.MPI -O -i eq5.in -o eq5.out -p system.prmtop -c eq4.rst -r eq5.rst -x eq5.nc -inf mdinfo -ref eq4.rst
pmemd.MPI -O -i eq.in -o eq.out -p system.prmtop -c eq5.rst -r eq.rst -x eq.nc -inf mdinfo -ref eq5.rst
pmemd.MPI -O -i prod.in -o prod.out -p system.prmtop -c eq.rst -r prod.rst -x prod.nc -inf mdinfo
```
> Notice how the `-c` and `-i` arguments change sequentially with each run. Also note that the tleap output `system.prmtop` is being used in all `-p` flags.

#### Post-Production Analysis
For the sake of brevity of this repository, there will be a separate repository to document the post-simulation analysis very soon.
