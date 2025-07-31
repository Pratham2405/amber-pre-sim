# This script collates the important commands require for AMBER simulations of a protein-ligand complex, including 
# parameterization and complex formation with a protein.

#!/bin/bash

#OBabel Command for converting PDBQT to SDF format.
obabel -ipdbqt 0.160_out.pdbqt -osdf -O 0.160_MD.sdf

#Obabel Command for converting SDF to PDB format.
obabel -isdf 0.160_MD.sdf -opdb -O 0.160_MD.pdb

#Obabel Command for removing hydrogens from the PDB file.
obabel -ipdb 0.160_MD_clean.pdb -opdb -O 0.160_MD_noH.pdb -d

#Calculates the net charge of a ligand from its pdbqt file.
grep -E '^(ATOM|HETATM)' 0.160_out.pdbqt | awk '{sum += $(NF-1)} END {print sum}'

# Remove CONECT records and let antechamber infer connectivity
sed 's/UNL/LIG/g' 0.160_MD.pdb | grep -E '^HETATM' > 0.160_MD_clean.pdb

#Adding hydrogens to the ligand.
reduce -build 0.160_MD_noH.pdb > 0.160_MD_H.pdb
antechamber -i 0.160_MD_H.pdb -fi pdb -o 0.160_MD_withH.pdb -fo pdb -c bcc -nc -2 -m 1 -at gaff2

#Antechamber parameter file preparation for ligand.
antechamber -i 0.160_MD_H.pdb -fi pdb -o 0.160_MD.mol2 -fo mol2 -c bcc -nc -2 -m 1 -at gaff2 -rn LIG -s 2 -j 4

#parmchk2: parmchk2 generates missing force field parameters that aren't available in the standard GAFF parameter set.
parmchk2 -i 0.160_MD.mol2 -f mol2 -o 0.160_MD.frcmod -at gaff2

#prepi file preparation: PREPI files contain residue connectivity and internal coordinate information for non-standard residues.
antechamber -i 0.160_MD.mol2 -fi mol2 -o 0.160_MD.prepi -fo prepi

#Combine protein and ligand in tleap.in as below:
tleap
source leaprc.protein.ff19SB
source leaprc.gaff2
source leaprc.water.tip3p
loadamberparams frcmod.ionsjc_tip3p
loadamberprep ligand.prepi    # Correct command for prepi files
loadamberparams ligand.frcmod # Correct command for frcmod files
ligand = loadmol2 ligand.mol2
protein = loadpdb protein.pdb
complex = combine {protein ligand}
check complex
solvateoct complex TIP3PBOX 12.0
addIons complex Cl- 6
saveamberparm complex system.prmtop system.inpcrd
savepdb complex system.pdb
quit

tleap -f leap.in > leap.log


 
# #!/bin/bash

# # Script for AMBER ligand parametrization and complex preparation
# set -e  # Exit on any error

# # Variables
# LIGAND_NAME="LIG"
# NET_CHARGE=$(grep -E '^(ATOM|HETATM)' ligand.pdbqt | awk '{sum += $(NF-1)} END {printf "%.0f", sum}')
# PROTEIN_PDB="protein.pdb"
# LIGAND_PDB="ligand.pdb"

# echo "Calculated net charge: $NET_CHARGE"

# # Step 1: Add hydrogens to ligand
# echo "Adding hydrogens to ligand..."
# antechamber -i $LIGAND_PDB -fi pdb -o ligand_withH.pdb -fo pdb -c bcc -nc $NET_CHARGE -m 1 -at gaff2

# # Step 2: Generate MOL2 file
# echo "Generating MOL2 file..."
# antechamber -i ligand_withH.pdb -fi pdb -o ligand.mol2 -fo mol2 -c bcc -nc $NET_CHARGE -m 1 -at gaff2 -rn $LIGAND_NAME -s 2

# # Step 3: Generate force field parameters
# echo "Generating force field parameters..."
# parmchk2 -i ligand.mol2 -f mol2 -o ligand.frcmod -at gaff2 -a Y

# # Step 4: Generate PREPI file
# echo "Generating PREPI file..."
# antechamber -i ligand.mol2 -fi mol2 -o ligand.prepi -fo prepi

# # Step 5: Create LEaP input file and run
# echo "Creating complex with LEaP..."
# # [Insert the enhanced LEaP input file creation from above]

# # Step 6: Run LEaP
# tleap -f leap.in > leap.log 2>&1

# # Step 7: Validate output
# echo "Validating output files..."
# # [Insert validation steps from above]

# echo "Complex preparation completed successfully!"
