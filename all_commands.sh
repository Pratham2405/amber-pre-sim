# This file contains all the commands used in this tutorial.

#OBabel Command for converting PDBQT to SDF format.
obabel -ipdbqt 0.160_out.pdbqt -osdf -O 0.160_MD.sdf

#Obabel Command for converting SDF to PDB format.
obabel -isdf 0.160_MD.sdf -opdb -O 0.160_MD.pdb

#Obabel Command for removing hydrogens from the PDB file.
obabel -ipdb 0.160_MD_clean.pdb -opdb -O 0.160_MD_noH.pdb -d

#Calculates the net charge of a ligand from its pdbqt file.
grep -E '^(ATOM|HETATM)' 0.160_out.pdbqt | awk '{sum += $(NF-1)} END {print sum}'

# Remove CONECT records and let antechamber infer connectivity. Only take the hetatom lines.
sed 's/UNL/LIG/g' 0.160_MD.pdb | grep -E '^HETATM' > 0.160_MD_clean.pdb

#Adding hydrogens to the ligand(optional).
reduce -build 0.160_MD_noH.pdb > 0.160_MD_H.pdb
antechamber -i 0.160_MD_H.pdb -fi pdb -o 0.160_MD_withH.pdb -fo pdb -c bcc -nc -2 -m 1 -at gaff2

#Antechamber parameter file preparation for ligand.
antechamber -i 0.160_MD_H.pdb -fi pdb -o 0.160_MD.mol2 -fo mol2 -c bcc -nc -2 -m 1 -at gaff2 -rn LIG -s 2 -j 4

#parmchk2: parmchk2 generates missing force field parameters that aren't available in the standard GAFF parameter set.
parmchk2 -i 0.160_MD.mol2 -f mol2 -o 0.160_MD.frcmod -at gaff2

#prepi file preparation: PREPI files contain residue connectivity and internal coordinate information for non-standard residues.
antechamber -i 0.160_MD.mol2 -fi mol2 -o 0.160_MD.prepi -fo prepi

#Run leap.in
tleap -f leap.in > leap.log

#Commands for Minimisation

#Commands for Post-minimisation runs
