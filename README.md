# PART - Plumed Automatic Restraining Tool

## Software Requirements

You will need a molecular dynamics (MD) engine patched with PLUMED. Pick an MD engine of your choice (for example GROMACS) and patch it with PLUMED ([for example PLUMED v2.7 install instructions](https://www.plumed.org/doc-v2.7/user-doc/html/_installation.html)). Make sure that the version of the MD engine is compatible with the PLUMED version. Some MD engines, such as AMBER, do not need to be patched with PLUMED as the default install already includes PLUMED.

You will need a Python distribution with the `numpy` module. The code was tested using the following versions, but we expect other versions to work as well:

* Python: 3.8.10
* numpy: 1.19.2

## Installation

Get the code using:

```
git clone 
```

Now either add PART to your PATH and give `PART.py` or move the `PART.py` script to the directory where you wish to use it.


## Usage

Basic usage uitilizes the following syntax:

```
PART.py -f <.pdb/.gro of input system> -rp <restraint pairs string> -dp <distance point string> > -o <.dat, output file name> 
```

The restraints pair string lists all residue type pairs you wish to restrain. Seperate the two pair members' residue name by a colon and seperate all pairs by a comma, for example `'P1:P1'` to prevent aggregation of P1 molecules. The distance point string should contain all residue names included in the restraints pair list, together with the atom on which to base the intermolecular distance calculation. Seperate the residue name and seperate eqch residue type by commas. Use for example `'P1:C1'` to indicate that P1 intermolecular distances distances should be calculated by calculating the distances between the C1 atoms. It is also possible to use center of masses of molecules, for example `'P1:COM:C1,C2,C3'` to use the center of mass of the C1, C2 and C3 carbons. See the examples under Example Usage for further info.


To list the full description of all the keyword arguments use:

```
PART.py -h
```

#### Note that PART program units are different from the units used in the PART paper.


## Example Usage

Example systems and output can be found under the `Examples` folder.


#### Example 1

To generate the PLUMED file of example 1 (.pdb file of an Amber system, pdb file generated using ambpdb):

```
PART.py -f complex.solvated.pdb -rp 'P1:P1,P1:P2,P2:P2' -dp 'P1:C1,P2:C1' -o test.dat
```

This command will create a PLUMED file `test.dat` with restraints applied between all P1 and all other P1 molecules, between P1 and P2 molecules and between all P2 and all other P2 molecules. The distance calculation will be based on the C1 atoms of P1 and P2. The restraint potential will have default PART parameters.

#### Example 2

To generate the PLUMED file of example 1 (.gro file of a GROMACS system):
```
PART.py -f solvated_complex.gro -rp 'BNZ:BNZ,PRP:PRP,BNZ:PRP,MAM:ACE' -dp 'BNZ:COM:C01-C02-C03-C04-C05-C06-H07-H08-H09-H10-H11-H12,PRP:C01,MAM:N01,ACE:C02' -o plumed.dat
```

This command will create a `plumed.dat` file that makes sure restraints are applied to prevent aggregation between benzenes (BNZs) themselves, between PRPs themselves, between PRP and BNZ, and between MAM and ACE. The distances are calculated using the COM of all benzene (BNZ) atoms, the C01 atom of PRP, the N01 of MAM and C02 atom of ACE.

