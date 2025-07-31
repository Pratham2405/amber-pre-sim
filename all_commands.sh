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
pmemd.MPI -O -i min_restrained.in -o min1.out -p system.prmtop -c system.inpcrd -r min1.rst -x min1.nc -inf min1.mdinfo -ref system.inpcrd
memd.MPI -O -i min_unrestrained.in -o min2.out -p system.prmtop -c min1.rst -r min2.rst -x min2.nc -inf min2.mdinfo -ref min1.rst

#Commands for Post-minimisation runs
pmemd.MPI -O -i heat.in -o heat.out -p system.prmtop -c min2.rst -r heat.rst -x heat.nc -inf mdinfo -ref min2.rst
pmemd.MPI -O -i eq1.in -o eq1.out -p system.prmtop -c heat.rst -r eq1.rst -x eq1.nc -inf mdinfo -ref heat.rst
pmemd.MPI -O -i eq2.in -o eq2.out -p system.prmtop -c eq1.rst -r eq2.rst -x eq2.nc -inf mdinfo -ref eq1.rst
pmemd.MPI -O -i eq3.in -o eq3.out -p system.prmtop -c eq2.rst -r eq3.rst -x eq3.nc -inf mdinfo -ref eq2.rst
pmemd.MPI -O -i eq4.in -o eq4.out -p system.prmtop -c eq3.rst -r eq4.rst -x eq4.nc -inf mdinfo -ref eq3.rst
pmemd.MPI -O -i eq5.in -o eq5.out -p system.prmtop -c eq4.rst -r eq5.rst -x eq5.nc -inf mdinfo -ref eq4.rst
pmemd.MPI -O -i eq.in -o eq.out -p system.prmtop -c eq5.rst -r eq.rst -x eq.nc -inf mdinfo -ref eq5.rst
pmemd.MPI -O -i prod.in -o prod.out -p system.prmtop -c eq.rst -r prod.rst -x prod.nc -inf mdinfo
