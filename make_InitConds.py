# The purpose of this file is to produce the "initconds"
#  for input into the QD_PLDM method after both an
#  optimization and normal modes calculation.
###      ~ Braden Weight 9/19/2020 ~       ###

# Syntax: python3 make_InitConds.py NAME.xyz
# NAME.xyz should look like regular xyz file. I will optimize for you.
#   NAtoms
#   (SPACE)
#   Type x y z
#   Type x y z
#   .
#   .
#   .

import numpy as np
import sys
import subprocess as sp

# Where is MOLPRO? I have mine set in bashrc.
MOLPRO = "$MOLPRO/molpro.exe"

# Where is the file "wigner.py"?
WIGNER = "$SHARC/wigner.py"

# How many trajectories will you run in QD_PLDM?
NTRAJ = 500

### BEGIN MAIN PROGRAM ###

# Read in geometry
name, fileType = sys.argv[1].split(".")
file01 = open(name + "." + fileType,"r")
atoms = [] # Store atom types
coords = [] # Store atom coordinates
NAtoms = 0 # Stores the number of atoms
for count, line in enumerate(file01):
    s = line.split()
    if (count == 0):
        NAtoms = int(s[0])
    if (count >= 2 and s != []):
        typ = str(s[0])
        x = s[1]
        y = s[2]
        z = s[3]
        atoms.append(typ)
        coords.append([x,y,z])
file01.close()

if (NAtoms != len(atoms)):
    print ("\nThe number of atoms does not match with geometry. Please check.\n")
    exit()

# Clean up directory
sp.call("rm *.log*", shell=True)
sp.call("rm *.out*", shell=True)
sp.call("rm *.xml*", shell=True)
sp.call("rm *.molden*", shell=True)

# Make MOLPRO input with new coordinates read from file
sp.call("cat templates/input_template_TOP > " + name + ".com", shell=True)
file01 = open(name + ".com","a")
file01.write(str(NAtoms) + "\n")
file01.write("\t" + name + " CARTESIAN COORDINATES\n")
for n in range(len(atoms)):
    file01.write(atoms[n] + "\t" + coords[n][0] + "\t" + coords[n][1] + "\t" + coords[n][2] + "\n")
file01.close()
sp.call("cat templates/input_template_BOTTOM >> " + name + ".com", shell=True)
sp.call("sed -i 's/$NAME/" + name + "/g' " + name + ".com", shell=True)

# Now we run the optimization and normal modes job.
sp.call(MOLPRO + " " + name + ".com", shell=True)

# Now we getthe wigner distribution of the normal mode geometries
sp.call("python2.7 " + WIGNER + " -n " + str(NTRAJ) + " " + name.lower() + ".molden", shell=True)




