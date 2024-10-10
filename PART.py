#!/usr/bin/env python3

from dataclasses import dataclass
import argparse
import numpy as np
import sys

@dataclass
class AtomGro:
    """
    This dataclass contains all atom information for one atom  provided  by a .gro file. Attributes are descriptive :). Note that every atom has a unique ID and resid (gro file format start with 0 again for atomid/resid 100,000)

    This "overflow" is controlled by atomid_overflow and resid_overflow. 

    Parameters
    -----------
    :grofile_line : the raw line of the grofile
    :atomid_overflow : counter on how many times 100,000 should be added to the atomid
    :resid_overflow : counter on how many times 100,000 should be added to the resid

    Attributes 
    -----------
    :resid 
    :resname
    :atomname
    :atomid
    :position
    :velocity_info: Boolean that stores wether velocity info is present
    :velocity

    """
    
    
    def __init__(self, grofileline,atomid_overflow,resid_overflow):
        self.resid = int(int(grofileline[0:5]) + 100000*resid_overflow)
        self.resname = grofileline[5:10].strip()
        self.atomname = grofileline[10:15].strip()
        self.atomid = int(int(grofileline[15:20])+atomid_overflow*100000)

        posx = float(grofileline[20:28].strip())
        posy = float(grofileline[28:36].strip())
        posz = float(grofileline[36:44].strip())
        
        self.velocity_info = True
        if grofileline[44:52].strip() == '':
            self.velocity_info = False
        else:
            velx = float(grofileline[44:52].strip())
            vely = float(grofileline[52:60].strip())
            velz = float(grofileline[60:68].strip())
            self.velocity = [velx, vely, velz]

        self.position = [posx, posy, posz]


class ParserGro:
    '''
    This class reads the .gro file and can perform parsing operations such as getting all atomids for specidic resids.


    Attributes
    -----------
    :natoms : number of atoms
    :boxvectors: vectors of the box
    :header : title line
    :atoms: list of atom instances

    Methods
    ----------
    print_atom_info(atomid):
        Prints all known informatio of an atom (aFtomid provided as input) in a fancy format

    resids_with_resname(resname)
        Returns all resids with a given resname

    atomids_with_resid(resid)
        Returns all atomids with a given resid
    
    atomids_with_atomname(self,atomname)
        Returns all atomids with a given atomname
        
    atomids_with_resname_and_atomname
        Returns all atomids with a given resname and atomname combo

    atomid_to_atomname
        Returns the atomname that is found for the given atomid

    write_file
        Writes out the gro file of the system stored in the class

    recount_atom_ids
        Function that recalculates all the atomids in the class, starting from 1.
    '''

    def __init__(self, grofile):
        '''
        Initialize the atoms list with instances of the atom class.
        '''

        self.atoms = []

        self.original_filename = grofile

        # 1. Calculate number of lines
        nlines = sum(1 for line in open(grofile,'r'))


        # 2. Build the atoms list
        with open(grofile,'r') as myfile: 
            linecounter = 0
            atomid_overflow = 0
            resid_overflow = 0
            
            all_lines = myfile.readlines()
            for line in all_lines:
                # All lines should contain atoms, except the first two and the last line
                if linecounter != 0 and linecounter != 1 and linecounter != (nlines-1):
                    self.atoms.append(AtomGro(line,atomid_overflow,resid_overflow))

                    if self.atoms[-1].atomid ==int(99999+100000*atomid_overflow):
                        atomid_overflow +=1
                    
                    # When to up resid: next line is an atom (not the last line!) with resid 0
                    if self.atoms[-1].resid == int(99999+100000*resid_overflow) and linecounter+1 != nlines -1 and AtomGro(all_lines[linecounter+1], atomid_overflow,resid_overflow).resid == 0:
                        resid_overflow += 1
                elif linecounter == 0:
                    self.fileheader = line
                
                elif linecounter == 1:
                    natoms = int(line)
                
                elif linecounter == (nlines-1) and natoms == nlines-3:
                    self.boxvectors_unprocessed = line
                    self.boxvectors = [float(i) for i in line.split()]

                elif linecounter == (nlines-1) and natoms == nlines-2:
                    print('Found no boxvector, please check if this is correct')
                    self.atoms.append(AtomGro(line,atomid_overflow,resid_overflow))

                    if self.atoms[-1].atomid ==int(99999+100000*atomid_overflow):
                        atomid_overflow +=1
                    
                    # When to up resid: next line is an atom (not the last line!) with resid 0
                    if self.atoms[-1].resid == int(99999+100000*resid_overflow) and linecounter+1 != nlines -1 and AtomGro(all_lines[linecounter+1], atomid_overflow,resid_overflow).resid == 0:
                        resid_overflow += 1

                    self.boxvectors_unprocessed = ''
                    self.boxvectors = None

                else:
                    raise ValueError('ERROR1: something is wrong with ParserGro linecounter')
                
                linecounter += 1
        
        # 3. Consistency check
        if natoms != len(self.atoms):
            raise ValueError('ERROR2: Number of atoms in atoms list does not correspond with natoms value')


    def resids_with_resname(self,resname):
        '''
        Returns a list of resids of atoms that match a given residue name
        '''
        
        result_list = []
        
        for atom in self.atoms:
            if atom.resname == resname and atom.resid not in result_list:
                result_list.append(atom.resid)

        return result_list

    def atomids_with_resid(self,resid):
        '''
        Returns a list of atomids that match a given resid
        '''
        result_list = []

        for atom in self.atoms:
            if atom.resid == resid:
                result_list.append(atom.atomid)

        

        return result_list

    def atomids_with_atomname(self,atomname):
        '''
        Returns a list of atomids that match a given atom name
        '''
        
        result_list = []

        for atom in self.atoms:
            if atom.atomname == atomname:
                result_list.append(atom.atomid)
        
        return result_list

    def atomids_with_resname_and_atomname(self, resname, atomname):
        '''
        Return a list of atomids that match a given residue name and atomname
        '''
        
        result_list = []

        for atom in self.atoms:
            if atom.resname == resname and atom.atomname == atomname:
                result_list.append(atom.atomid)
        
        return result_list

    def atomid_to_atomname(self,atomid):
        '''
        Returns the atomname that is found for the given atomid

        '''

        for atom in self.atoms:
            if atom.atomid == atomid:
                myatomname = atom.atomname
        return myatomname

    def write_file(self, filename):
        '''
        Writes out the gro file of the system stored in the class
        '''
        def _build_resid_id_str(atom):
            id_for_file = atom.resid
            while id_for_file > 99999:
                id_for_file = id_for_file - 100000
            
            id_for_file = str(id_for_file)

            while len(id_for_file) != 5:
                id_for_file = f' {id_for_file}'
            
            return id_for_file

        def _build_atom_resname_str(atom):
            resn_for_file = atom.resname
            
            while len(resn_for_file) != 5:
                resn_for_file = f'{resn_for_file} '            
            
            return resn_for_file

        def _build_atom_atomname_str(atom):
            atom_atomname_str = atom.atomname
            
            while len(atom_atomname_str) != 5:
                atom_atomname_str = f' {atom_atomname_str}'
            return atom_atomname_str

        def _build_atom_id_str(atom):
            id_for_file = atom.atomid
            while id_for_file > 99999:
                id_for_file = id_for_file - 100000
            
            id_for_file = str(id_for_file)

            while len(id_for_file) != 5:
                id_for_file = f' {id_for_file}'
            
            return id_for_file
        
        def _build_position_str(atom):
            x = str(atom.position[0])
            y = str(atom.position[1])
            z = str(atom.position[2])



            while len(x) != 8:
                x = f' {x}'
            while len(y) != 8:
                y = f' {y}'
            while len(z) != 8:
                z = f' {z}'

            while x[-4] != '.':
                x = f'{x}0'
                x = x[1:]
            while y[-4] != '.':
                y = f'{y}0'
                y = y[1:]
            while z[-4] != '.':
                z = f'{z}0'
                z = z[1:]

            return f'{x}{y}{z}'

        def _build_velocity_str(atom):
            if atom.velocity_info == False:
                x = ''
                y = ''
                z = ''
            
            else:
                x = str(atom.velocity[0])
                y = str(atom.velocity[1])
                z = str(atom.velocity[2])



                while len(x) != 8:
                    x = f' {x}'
                while len(y) != 8:
                    y = f' {y}'
                while len(z) != 8:
                    z = f' {z}'

                while x[-5] != '.':
                    x = f'{x}0'
                    x = x[1:]
                while y[-5] != '.':
                    y = f'{y}0'
                    y = y[1:]
                while z[-5] != '.':
                    z = f'{z}0'
                    z = z[1:]

            return f'{x}{y}{z}'


        # 1. Calculat natoms
        n_atoms = len(self.atoms)
        
        # 2. Build the file 
        myfile_str = f'{self.fileheader}{n_atoms}\n'

        for atom in self.atoms:
            # 2.1 Build atom id str
            atom_id_str = _build_atom_id_str(atom)
            atom_resname_str = _build_atom_resname_str(atom)
            atom_atomname_str = _build_atom_atomname_str(atom)
            atom_resid_str = _build_resid_id_str(atom)
            position_str =  _build_position_str(atom)
            velocity_str = _build_velocity_str(atom)
            
            myfile_str += f'{atom_resid_str}{atom_resname_str}{atom_atomname_str}{atom_id_str}{position_str}{velocity_str}\n'

        myfile_str += self.boxvectors_unprocessed

        # 3. Write out the file
        if filename == self.original_filename:
            cont = input('Filename for new file is the same as the original filename. File will be overwritten, proceed? (y/n)')
            if cont == 'y':
                pass
            
            else:
                sys.exit('Exit was requestred by the user to not overwrite the file')

        with open(filename,'w') as myfile:
            myfile.write(myfile_str)

    def recount_atom_ids(self):
        '''
        Function that recalculates all the atomids in the class, starting from 1. If followed by write_file, a .gro file with wrong atomids can be corrected.
        '''
        mycounter = 1
        for atom in self.atoms:
            atom.atomid = mycounter
            mycounter += 1


    def print_atom_info(self,atomid):
        '''
        Print some atoms information of a given atomid
        '''
        thisatom = self.atoms[atomid-1]
        
        if thisatom.atomid != atomid:
            raise ValueError('Something is wrong with atomids')
        
        print('-----------')
        print('Fancy print of the requested atom information:')
        print('Atomid:', thisatom.atomid)
        print('Atomname:', thisatom.atomname)
        print('Resid:', thisatom.resid)
        print('Resname:', thisatom.resname)
        print('Position:', thisatom.position)
        print('Velocities:', thisatom.velocity)
        print('-----------')
    

@dataclass
class AtomPdb:
    """
    This dataclass contains all atom information for one atom  provided  by a .pdb file. Attributes are descriptive :). Note that every atom has a unique ID and resid (pdb file format start with 0 again for atomid/resid 100,000)

    This "overflow" is controlled by atomid_overflow and resid_overflow. 

    Parameters
    -----------
    :pdbfile_line : the raw line of the pdbfile
    :atomid_overflow : counter on how many times 100,000 should be added to the atomid
    :resid_overflow : counter on how many times 100,000 should be added to the resid

    Attributes 
    -----------
    :resid 
    :resname
    :atomname
    :atomid
    :chain
    :position

    """
    
    
    def __init__(self, pdbfileline,atomid_overflow,resid_overflow):
        self.resid = int(int(pdbfileline[22:26]) + 100000*resid_overflow) #c
        self.resname = pdbfileline[17:20].strip() #c
        self.atomname = pdbfileline[12:16].strip() #c
        self.atomid = int(int(pdbfileline[6:11])+atomid_overflow*100000) #c
        
        self.chain = pdbfileline[21] #c
        
        x = float(pdbfileline[30:38].strip())
        y = float(pdbfileline[38:46].strip())
        z = float(pdbfileline[46:54].strip())

        self.position = [x,y,z] #c


class ParserPdb:
    '''
    This class reads the .pdb file and can perform parsing operations such as getiing all atomids for specidic resids.


    Attributes
    -----------

    :atoms: list of atom instances

    Methods
    ----------
    print_atom_info(atomid):
        Prints all known informatio of an atom (atomid provided as input) in a fancy format

    resids_with_resname(resname)
        Returns all resids with a given resname

    atomids_with_resid(resid)
        Returns all atomids with a given resid
    
    atomids_with_atomname(self,atomname)
        Returns all atomids with a given atomname
        
    atomids_with_resname_and_atomname
        Returns all atomids with a given resname and atomname combo
    
    atomid_to_atomname
        Returns the atomname that is found for the given atomid
    '''

    def __init__(self, pdbfile):
        '''
        Initialize the atoms list with instances of the atom class.
        '''

        self.atoms = []

        # 1. Calculate number of lines
        nlines = sum(1 for line in open(pdbfile,'r'))


        # 2. Build the atoms list
        with open(pdbfile,'r') as myfile: 
            linecounter = 0
            atomid_overflow = 0
            resid_overflow = 0
            
            all_lines = myfile.readlines()
            for line in all_lines:
                # All lines that start with ATOM or HETATM
                if 'ATOM' in line or 'HETATM' in line:
                    self.atoms.append(AtomPdb(line,atomid_overflow,resid_overflow))

                    #Corrction 1
                    if self.atoms[-1].atomid ==int(99999+100000*atomid_overflow):
                        atomid_overflow +=1
                    
                    #Correction 2
                    # When to up resid: next line is an atom (not the last line!) with resid 0
                    if self.atoms[-1].resid == int(99999+100000*resid_overflow) and linecounter+1 != nlines -1 and AtomPdb(all_lines[linecounter+1], atomid_overflow,resid_overflow).resid == 0:
                        resid_overflow += 1
                
                linecounter += 1
        

    def resids_with_resname(self,resname):
        '''
        Returns a list of all resids with a given resname
        '''
        result_list = []
        
        for atom in self.atoms:
            if atom.resname == resname and atom.resid not in result_list:
                result_list.append(atom.resid)

        return result_list

    def atomids_with_resid(self,resid):
        '''
        Returns a list of all atomids with given resids
        '''
        result_list = []

        for atom in self.atoms:
            if atom.resid == resid:
                result_list.append(atom.atomid)

        

        return result_list

    def atomids_with_atomname(self,atomname):
        '''
        Returns a list of atomids with given atomnames
        '''
        result_list = []

        for atom in self.atoms:
            if atom.atomname == atomname:
                result_list.append(atom.atomid)
        
        return result_list

    def atomids_with_resname_and_atomname(self, resname, atomname):
        '''
        Returns a list of atomnames of the atoms that have the given resname and atomname
        '''
        
        result_list = []

        for atom in self.atoms:
            if atom.resname == resname and atom.atomname == atomname:
                result_list.append(atom.atomid)
        
        return result_list

    def atomnames_with_resname(self,resname):
        '''
        Returns a list of atomnames that are found for the given resname
        Contains duplicates for each occurence of the chosen residue

        '''
        result_list = []

        for atom in self.atoms:
            if atom.resname == resname:
                result_list.append(atom.atomname)

        return result_list

    def atomid_to_atomname(self,atomid):
        '''
        Returns the atomname that is found for the given atomid

        '''

        for atom in self.atoms:
            if atom.atomid == atomid:
                myatomname = atom.atomname
        return myatomname

    def print_atom_info(self,atomid):
        '''
        Print some atom information of an atom with a given atomid
        '''
        thisatom = self.atoms[atomid-1]
        
        if thisatom.atomid != atomid:
            raise ValueError('Something is wrong with atomids')
        
        print('-----------')
        print('Fancy print of the requested atom information:')
        print('Atomid:', thisatom.atomid)
        print('Atomname:', thisatom.atomname)
        print('Resid:', thisatom.resid)
        print('Resname:', thisatom.resname)
        print('Position:', thisatom.position)
        print('-----------')

def find_atom_group_gro(grofile,residue_type,central_atom_for_probe):
    '''
    Returns all atoms ids of the grofile with the wanted residue type and central_atom_for_probe in a comma separated string

    Inputs
    ---------
    grofile: string of the filename
    residue_type: string of residue name
    central_atom_for_probe: string of the atom name

    Outputs
    ---------
    atomids_str: string of all atomids matching the residue type and taomtype, comma separated
    '''


    grofile_parsed = ParserGro(grofile)
    atom_id_list = grofile_parsed.atomids_with_resname_and_atomname(residue_type,central_atom_for_probe)
    atomids_str = ','.join([str(x) for x in atom_id_list])

    return atomids_str

def find_atom_group_pdb(pdbfile,residue_type,central_atom_for_probe):
    '''
    Returns all atoms ids of the pdbfile with the wanted residue type and central_atom_for_probe in a comma separated string

    Inputs
    ---------
    pdbfile: string of the filename
    residue_type: string of residue name
    central_atom_for_probe: string of the atom name    

    Outputs
    ---------
    atomids_str: string of all atomids matching the residue type and taomtype, comma separated
    '''


    pdbfile_parsed = ParserPdb(pdbfile)
    atom_id_list = pdbfile_parsed.atomids_with_resname_and_atomname(residue_type,central_atom_for_probe)
    atomids_str = ','.join([str(x) for x in atom_id_list])

    return atomids_str


def find_atom_group_COM(strucutrefile, mode,residue_type, atomnames_for_com):
    '''
    Produces two necessary strings to be included in the plumed restraint file.

    Inputs
    -------
    structurefile: string, filename of boxed protein .gro/.pdb file
    mode: string, pdb mode or gro mode
    atomnames_for_com: list of strings with atomnames to be selected

    Outputs
    --------
    atomgroup_str: str with the atomgroup of all benzene COMS
    com_calculation_str: str with all the COM calculation statements. 
    '''

    if mode  == 'pdb':
        file_parsed = ParserPdb(strucutrefile)
        com_resid_list = file_parsed.resids_with_resname(residue_type)

    elif mode  == 'gro':
        file_parsed = ParserGro(strucutrefile)
        com_resid_list = file_parsed.resids_with_resname(residue_type)


    else:
        raise ValueError(" Invalid file mode")

    counter =1

    atomgroup_str =''
    com_calculation_str = '\n'

    for resid in com_resid_list:
    
        atomids_raw = file_parsed.atomids_with_resid(resid)

        atomids = []
        for atomid in atomids_raw:
            if file_parsed.atomid_to_atomname(atomid) in atomnames_for_com:
                atomids.append(atomid)

        if resid != com_resid_list[-1]:
            atomgroup_str += f'{residue_type.lower()}{counter},'
        else:
            atomgroup_str += f'{residue_type.lower()}{counter}'

        str_thiscom = f'{residue_type.lower()}{counter}: COM ATOMS='
        
        
        for atomid in atomids:
            if atomid != atomids[-1]:
                str_thiscom+=str(atomid)+','
            else:
                str_thiscom+=str(atomid)+'\n'
        
        com_calculation_str += str_thiscom
        counter +=1

    return atomgroup_str, com_calculation_str



# Used in main
def command_line_parser():
    """
    Performs argument parsing
       
    """

    my_parser = argparse.ArgumentParser(description='''
This script builds a plumed file that will prevent aggregation between specified residues based on a strucure file provided by the user.
The format of the structure file must be either a .pdb or .gro file of the full system (all probes present, atom identifiers equal to those in the MD simulation). 
The user can choose to which residue pairs the artificial repulsion potential should be applied.
The user should specify, for every residue type that is subjugate to an artificial repulsion, which atom of the residue should be used to calculate the interresidue distances. It is also possible to use a virtual atom for the distance calculation, where the virtual atom is located at the center of mass of some chosen atoms of a residue.
The plumed file can be customized after running the script, for example to include additional RMSD restraints.
    ''')
   
    my_parser.add_argument("-f", "--StructureFile", help = 'Structure file of the full system. Supported formats: .gro or .pdb', type = str) 

    my_parser.add_argument("-rp", "--RestraintPairs", help = """Residue pairs to which an artificial repulsion potential should be applied. Residue couples should be joined by a colon, while several pairs should be separated by ','.
Examples: 'BNZ:BNZ' will apply the arificial repulsion potential only between all benzene residue pairs with residue names BNZ.
'BNZ:BNZ,PRP:BNZ,PRP:PRP' will apply the repulsion potential to all benzene pairs, to all benzene-propane pairs, where the residue name of propane is PRP and to all propane pairs. """ \
            , type = str)
    
    my_parser.add_argument("-dp", "--DistancePoints", help = """Input argumant that will determine which atoms will be used for the interresidue distance calculation. Every residue present in the RestraintPairs should be present in this argument. Central atoms of each residue present in a residue restraint pair can by passed to the script by typing the residue name, a colon, and then the atomname of the chosen atom.
It is is possible to use virtual atoms that are at the center of mass of a chosen set of atoms using the residue name, a colon, COM and colon followed by the atom names of which the COM should be calculated, separated by a hyphen. Examples:
'PRP:C01' will base the distance calculation on the C01 atom of propane zith residue name PRP. 'BNZ:COM:C01-C02-C03-C04-C05-C06,PRP:C01' will base the distance calculation of benzene (residue name BNZ) on the center of mass of the atoms C01 to C06 and the distance calculation of PRP will be based on the C01 atom of PRP.
 """ \
                , type=str) 
    
    my_parser.add_argument("-o", "--Output", help = "Name of the output plumed file. Default: plumed.dat", type=str, default='plumed.dat') 

    my_parser.add_argument("-ul", "--UpperLimit", help = "Optional input argument. Changes the upper limit of when the distance restraints become active (LWALLS AT keyword). Units of nm. Default is 0.8.", type=float, default =0.8) 

    my_parser.add_argument("-k", "--Kappa", help = "Optional input argument. Force constant for the restraint potential (LWALLS KAPPA keyword). Units of kJ/(mol*nm**x), with x being the value of the -exp argument. Default is 20000.0 ",default=20000.0,type=float)

    my_parser.add_argument("-exp", "--EXP", help = "Optional input argument. Exponent for the restraint potential (LWALLS EXP keyword). Default is 4.0",default=4.0,type=float)

    my_parser.add_argument("-eps", "--EPS", help = "Optional input argument. Rescaling facor for the restraint potential (LWALLS EPS keyword). Default is 1.0",default=1.0,type=float)

    # When no arghuments given
    if len(sys.argv)==1:
        my_parser.print_help()
        my_parser.exit()


    # Parse the arguments
    args = my_parser.parse_args()

    myfile =  args.StructureFile
    rp =  args.RestraintPairs
    dp =  args.DistancePoints
    output =  args.Output
    upper_limit =  args.UpperLimit
    kappa =  args.Kappa
    exp =  args.EXP
    eps =  args.EPS

    # Process the arguments
    ## Process RestraintPairs
    self_restraint_list = []
    nonself_restraint_list = []
    
    mypairs = rp.split(',')
    for entry in mypairs:
        thispair = entry.split(':')
        if thispair[0] == thispair [1]:
            self_restraint_list.append(thispair[0])
        
        else:
            nonself_restraint_list.append([thispair[0],thispair[1]])

    ## Process DistancePoints
    distance_info = {}
    com_info = {}      
    individualinfos = dp.split(',')

    for entry in individualinfos:
        thisinfo = entry.split(':') 

        if thisinfo [1] == 'COM' and len(thisinfo) == 3:
            distance_info[thisinfo[0]] = thisinfo[1]
            atomlist = thisinfo[2].split('-')
            com_info[thisinfo[0]] = atomlist

        elif len(thisinfo) == 2:
            distance_info[thisinfo[0]] = thisinfo[1]

        else:
            raise ValueError("ERROR: Could not parse the -dp input properly. Please check the input format.")


    # Sanity check
    ## All in selfrestraintlist and nonselfrestraintlist in distinfo?
    for entry1 in self_restraint_list:
            if entry1 in distance_info:
                pass
            else:
                raise ValueError('ERROR: Found entries in -rp that were not present in -dp. Please provide all necessary information')

    for entry1 in nonself_restraint_list:
        for entry in entry1:    
            if entry in distance_info:
                pass
            else:
                raise ValueError('ERROR: Found entries in -rp that were not present in -dp. Please provide all necessary information')

    ## Also check the reverse
    for entry in distance_info:
        if entry in self_restraint_list:
            pass
        else:
            if entry in [item for sublist in nonself_restraint_list for item in sublist]:
                pass
            else:
                raise ValueError('ERROR: Found entries in -dp that were not present in -rp. Please provide all necessary information')

    return myfile, self_restraint_list, nonself_restraint_list, distance_info, com_info, output, upper_limit, kappa, exp, eps


def write_plumed_file( mode,myfile, self_restraint_list, nonself_restraint_list, distance_info, com_info, output, upper_limit, kappa, exp, eps):
    """
    Writes the plumed multibias file
    """


    # 1. Build the atom groups and the com calculation statement
    atomgroups_str = {}
    plumed_file_com_calcs = ''
    for residue_type in distance_info:
        if distance_info[residue_type] == 'COM':
            atomgroups_str[residue_type], this_com_calc_string = find_atom_group_COM(myfile, mode,residue_type, com_info[residue_type])

            if this_com_calc_string not in plumed_file_com_calcs:
                plumed_file_com_calcs += this_com_calc_string

        else:
            if mode == 'pdb':
                atomgroups_str[residue_type] = find_atom_group_pdb(myfile,residue_type, distance_info[residue_type])
            elif mode == 'gro':
                atomgroups_str[residue_type] = find_atom_group_gro(myfile,residue_type, distance_info[residue_type])

            else:
                raise ValueError("Invalid file mode")

    # 2. Compose the plumed file
    plumed_file = plumed_file_com_calcs
    plumed_print_statement = '\n'
    bias_counter = 1
    
    for restype in self_restraint_list:
        str_distances = f'\nd{bias_counter}: DISTANCES GROUP={atomgroups_str[restype]}'
        str_restraint = f'\nb{bias_counter}: LWALLS DATA=d{bias_counter} AT={upper_limit} KAPPA={kappa} EXP={exp} EPS={eps} \n'
        plumed_file +=  str_distances + str_restraint
        plumed_print_statement += f'PRINT ARG=b{bias_counter}.* FILE=COLVAR-BIAS{bias_counter} STRIDE=500 \n'
        bias_counter += 1

    for pair in nonself_restraint_list:
        restype1 = pair[0]
        restype2 = pair[1]
        str_distances = f'\nd{bias_counter}: DISTANCES GROUPA={atomgroups_str[restype1]} GROUPB={atomgroups_str[restype2]}'
        str_restraint = f'\nb{bias_counter}: LWALLS DATA=d{bias_counter} AT={upper_limit} KAPPA={kappa} EXP={exp} EPS={eps} \n'
        plumed_file +=  str_distances + str_restraint
        plumed_print_statement += f'PRINT ARG=b{bias_counter}.* FILE=COLVAR-BIAS{bias_counter} STRIDE=500 \n'
        bias_counter += 1

    plumed_file += plumed_print_statement

    with open(output, 'w') as myfile:
        myfile.write(plumed_file)

    # Print success message
    print(f'PART created the following output file: {output}')

def main ():
    """
    Controls the execution of the script
    """

    # 1. Parsing the command line
    myfile, self_restraint_list, nonself_restraint_list, distance_info, com_info, output, upper_limit, kappa, exp, eps = command_line_parser()

    # 2. Determine the correct mode (pdb/gro)
    if myfile[-3:] == 'pdb' or myfile[-3:] == 'gro':
        mode = myfile[-3:]
    else:
        raise ValueError('ERROR: Please use a supported structure file format')

    # 3. Writing the plumed file
    write_plumed_file(mode, myfile, self_restraint_list, nonself_restraint_list, distance_info, com_info, output, upper_limit, kappa, exp, eps)

main()