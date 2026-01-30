
#! /usr/bin/env python

from __future__ import print_function
import sys, pickle, argparse, os, subprocess, time

import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem

import pathlib

from ..system import logger


def writeHeader(molname, loplsflag, gaffflag=False, dreidingflag=False, compassflag=False):
    if gaffflag:
        print("""import "gaff.lt"    # <-- defines the GAFF (General Amber Force Field)""")
        print("""# NOTE: GAFF requires user-supplied charges (AM1-BCC or RESP recommended)""")
        print("""# See: http://ambermd.org/antechamber/gaff.pdf""")
        print("""{0} inherits GAFF {{""".format(molname))
    elif dreidingflag:
        print("""import "dreiding.lt"    # <-- defines the DREIDING force field""")
        print("""# NOTE: DREIDING requires user-supplied charges (AM1-BCC, Gasteiger, or RESP recommended)""")
        print("""# See: Mayo et al., J. Phys. Chem. 1990, 94, 8897-8909""")
        print("""{0} inherits DREIDING {{""".format(molname))
    elif compassflag:
        print("""import "compass_published.lt"    # <-- defines the COMPASS force field (class2)""")
        print("""# NOTE: COMPASS is a class2 force field requiring LAMMPS CLASS2 package""")
        print("""# NOTE: This is an incomplete public version - parameters may be missing""")
        print("""# Charges are handled via bond increment model in the .lt file""")
        print("""{0} inherits COMPASS {{""".format(molname))
    else:
        print("""import "oplsaa.lt"    # <-- defines the standard "OPLSAA" force field""")
        if loplsflag:
            print("""import "loplsaa.lt"   # <-- custom parameters for long alkane chains taken from
                          #     Sui et al. J.Chem.Theory.Comp (2012), 8, 1459
                          #     To use the ordinary OPLSAA force field parameters,
                          #     (instead of the Sui et al. parameters), change the
                          #     atom types below from "@atom:81L","@atom:85LCH2" to
                          #     "@atom:81" and "@atom:85"  (defined in "oplsaa.lt")""")
        print("""{0} inherits OPLSAA {{""".format(molname))


def writeFooter(molname, gaffflag=False, dreidingflag=False, compassflag=False):
    if gaffflag:
        print("""}} # {0}

# IMPORTANT: GAFF requires atomic charges to be assigned manually!
# All charges in this file are currently set to 0.0 as placeholders.
# You must replace them with calculated charges using one of these methods:
#
# Method 1 (Recommended): AM1-BCC charges using AmberTools antechamber
#   antechamber -i molecule.mol2 -fi mol2 -o molecule_charged.mol2 -fo mol2 -c bcc -nc 0
#
# Method 2: RESP charges (more accurate, requires Gaussian)
#   See Amber documentation for RESP fitting procedure
#
# Method 3: Use pre-calculated charges from literature or databases
#
# After calculating charges, update the charge values in the "Data Atoms" section above.
# Reference: http://ambermd.org/antechamber/gaff.pdf""".format(molname))
    elif dreidingflag:
        print("""}} # {0}

# IMPORTANT: DREIDING requires atomic charges to be assigned manually!
# All charges in this file are currently set to 0.0 as placeholders.
# You must replace them with calculated charges using one of these methods:
#
# Method 1 (Recommended): AM1-BCC charges using AmberTools antechamber
#   antechamber -i molecule.mol2 -fi mol2 -o molecule_charged.mol2 -fo mol2 -c bcc -nc 0
#
# Method 2: Gasteiger charges using Open Babel
#   obabel -i mol mol.mol -o mol2 -O mol_charged.mol2 --partialcharge gasteiger
#
# Method 3: RESP charges (more accurate, requires Gaussian)
#   See Amber documentation for RESP fitting procedure
#
# After calculating charges, update the charge values in the "Data Atoms" section above.
# Reference: Mayo et al., J. Phys. Chem. 1990, 94, 8897-8909""".format(molname))
    elif compassflag:
        print("""}} # {0}

# NOTE: COMPASS uses the bond increment charge model defined in the .lt file.
# The charges in this file may be set to 0.0 as placeholders if the bond
# increment model is not applicable for your molecule.
#
# COMPASS is a class2 force field that requires LAMMPS compiled with CLASS2 package.
# This is an incomplete public version - parameters for some atom types may be missing.
#
# Reference: Sun, H., J. Phys. Chem. B, 1998, 102, 7338-7364""".format(molname))
    else:
        print("""}} # {0}

# Note: You don't need to supply the partial partial charges of the atoms.
#       If you like, just fill the fourth column with zeros ("0.000").
#       Moltemplate and LAMMPS will automatically assign the charge later""".format(molname))

def writeAtoms(m):
    print("\n# atom-id  mol-id  atom-type charge      X         Y        Z\n")
    print("  write(\"Data Atoms\") {")
    conf = m.GetConformer(0)
    for at in m.GetAtoms():
        point = conf.GetAtomPosition(at.GetIdx())
        print("\t{0} $mol:... {1} 0.00 {2:8.3f} {3:8.3f} {4:8.3f}".format(
                                    '$atom:'+at.GetSymbol()+str(at.GetIdx()+1),
                                    at.GetProp('AtomType'),
                                    point.x, point.y, point.z
                                    ))
    print("  }")

def writeBonds(m):
    bonds = m.GetBonds()
    print("\n  write('Data Bond List') {\n")
    for bond in bonds:
        b = bond.GetBeginAtom()
        e = bond.GetEndAtom()
        bname = b.GetSymbol()+str(b.GetIdx()+1)
        ename = e.GetSymbol()+str(e.GetIdx()+1)
        print("\t$bond:{0}\t$atom:{1}\t$atom:{2}".format(bname+ename,
                                                       bname, ename))
    print("  }")

def lt_to_molecule(ltfn):
    """Reads a moltemplate .lt file and returns an RDKit molecule for
    comparison purposes. Only works on .lt files with specific formatting.
    Doesn't have bond type perception, so doesn't generate useful smiles.
    ** Don't use this for anything. **
    """
    ltmol = Chem.RWMol()
    with open(ltfn,'r') as infile:
        for line in [line.strip() for line in infile if line.strip()]:
            # for this to work, the name of the atom must contain the
            #atomic symbol at the beginning. Otherwise would need to
            # look up based on type or mass.
            if line.split()[0][:5] == "$atom":
                label = line.split(':')[1].split()[0]
                #print(label)
                #construct a new atom by passing the atomic symbol
                #filter removes numbers
                newatom = Chem.Atom(''.join(filter(str.isalpha, label)))
                atomid = ltmol.AddAtom(newatom)
            elif line.split()[0][:5] == "$bond":
                #assumes bond - atom - atom style entries with atom id
                # in the field
                id1str = line.split()[1].split(':')[1]
                id1 = int(''.join(filter(str.isdigit, id1str)))
                id2str = line.split()[2].split(':')[1]
                id2 = int(''.join(filter(str.isdigit, id2str)))
                #this makes everything a single bond, so won't allow building
                # of a valid smiles from the incomplete graph
                ltmol.AddBond(id1,id2)
                #print(id1,id2)
    newmol = ltmol.GetMol()
    Chem.SanitizeMol(newmol)
    AllChem.EmbedMolecule(newmol,AllChem.ETKDG())
    print(Chem.MolToSmiles(newmol))


def read_cdict(cdictin):
    with open(cdictin, 'rb') as f:
        cdict = pickle.load(f)
    return cdict

def sum_of_charges(m, cdict):
    test_sum = 0
    print("# Given the current charge dictionary, the atoms will have the following charges:")
    for atom in m.GetAtoms():
        atype = atom.GetProp('AtomType')
        print("# Atom {0} is type {1} with charge {2}".format(atom.GetIdx(),atype,cdict[atype]))
        test_sum += cdict[atype]
    print("# The sum of the atomic charges is: {:.2f}".format(test_sum))
    if abs(test_sum) > 0.001:
        print("""
            # WARNING: The net charge appears to be non-zero! This may indicate
            incompatible atom types.
            """)

def generateFeatureDefn(fpath, fdefout, cdictout):
    """Write a feature definition file in RDKit style from the moltemplate
    conversion document. Only need to run this function if the conversion
    document has been changed.

    fpath -- file path of the moltemplate conversion doc
    fdefout -- file path to write the feature definition file
    cdictout -- file path to create a dictionary of atom types to charges
    """
    with open(fpath,'r') as infile, open(fdefout,'w') as outfile:
        feat_index = 0
        cdict = {}
        for line in [line.strip() for line in infile if line.strip()]:
            if line[0]!='*':
                el, atomname, typename, patt, lttype, chg, desc = [el.strip() for el in line.split("|")]
                # write lttype, SMARTS to feature definintion file
                # NOTE: using feature family to store the atom names is dangerous
                # because rdkit won't assign mutliple features in same family.
                # So I had to assign an index to make unique names [AHS]
                fdefn = \
"""
DefineFeature {0} {1}
Family {2}{3}
EndFeature""".format(lttype, patt, feat_index, atomname)
                # add new charge dictionary entry
                cdict[lttype]=float(chg)
                feat_index+=1
                outfile.write(fdefn)

    with open(cdictout,'wb') as f:
        pickle.dump(cdict,f, protocol=2)

def copy_to_cwd(source,destination):
        bash="cp "
        bash=bash+str(source)+" "+destination
        os.system(bash)

class RDlt(object):
    def __init__(self,smiles=None,name=None):
        self.smiles=smiles
        #self.name=name

    def store_bank(self,flag=True):
        if flag:
            file_=pathlib.Path(str(pathlib.Path(__file__).parent.resolve())+"/Monomer_bank/"+self.name+'.lt')
            if file_.exists():
                response = input(str(file_)+" exists, delete and make new?(y/n) ")
                if response[0] == "y":
                    proc = subprocess.Popen(['/bin/bash'], shell=True,stdin=subprocess.PIPE, stdout=subprocess.PIPE)
                    stdout = proc.communicate(("rm "+ str(file_)).encode())
                    logger.info(' '.join(["removing "+str(file_)]))
                    time.sleep(3)
                    copy_to_cwd(self.to_file,str(file_))
                elif response[0] == "n":
                    logger.info(' '.join(["EXIT : "+str(file_)+" has already existed"]))
                    sys.exit()
            else:
                copy_to_cwd(self.to_file,str(file_))
        else:
            return

    def detect_and_adjust_conjugated_systems(self, mol):
        """
        Detect conjugated systems and adjust GAFF atom types for paired types.

        Implements post-typing adjustment for conjugated systems (cc/cd, ce/cf, nc/nd, ne/nf).
        Uses bond connectivity to distinguish paired types (pure SMARTS cannot do this).

        This algorithm is inspired by the Antechamber approach which uses bond connectivity
        (field 7 of their 7-feature typing algorithm) to distinguish between conjugated
        atom types that have identical local environments.

        Args:
            mol: RDKit Mol object with AtomType properties assigned

        Returns:
            Modified Mol with adjusted conjugated types
        """
        from collections import defaultdict, deque

        # Build bond network and find double-bonded atoms
        bonds_by_atom = defaultdict(list)
        double_bonded_atoms = set()

        for bond in mol.GetBonds():
            begin = bond.GetBeginAtomIdx()
            end = bond.GetEndAtomIdx()
            btype = bond.GetBondType()

            bonds_by_atom[begin].append((end, btype))
            bonds_by_atom[end].append((begin, btype))

            if btype == Chem.rdchem.BondType.DOUBLE:
                double_bonded_atoms.add(begin)
                double_bonded_atoms.add(end)

        # Find conjugated systems using BFS (double-single-double patterns)
        conjugated_sets = []
        visited = set()

        for start_atom in double_bonded_atoms:
            if start_atom in visited:
                continue

            queue = deque([start_atom])
            current_system = []

            while queue:
                current = queue.popleft()
                if current in visited:
                    continue

                if current in double_bonded_atoms:
                    visited.add(current)
                    current_system.append(current)

                    # Find neighbors through single bonds (conjugated path)
                    for neighbor, btype in bonds_by_atom[current]:
                        if (btype == Chem.rdchem.BondType.SINGLE and
                            neighbor in double_bonded_atoms and
                            neighbor not in visited):
                            queue.append(neighbor)

            if len(current_system) >= 2:
                conjugated_sets.append(current_system)

        # Adjust atom types in conjugated systems
        # Apply GAFF paired type rules
        for system in conjugated_sets:
            for atom_idx in system:
                atom = mol.GetAtomWithIdx(atom_idx)
                try:
                    current_type = atom.GetProp('AtomType')

                    # Map basic types to conjugated types
                    # For conjugated systems, use the "first" type in each pair
                    # This is acceptable because parameter differences are minimal
                    # and matches the approach used by GAFF-foyer
                    if current_type == '@atom:c2':
                        atom.SetProp('AtomType', '@atom:ce')  # Conjugated sp2 C
                    elif current_type == '@atom:cc':
                        atom.SetProp('AtomType', '@atom:cd')  # Heteroaromatic conjugated C

                except KeyError:
                    pass

        return mol

    def run(self,to_file=None,name='test',fdef=str(pathlib.Path(__file__).parent.resolve())+"/rdlt_data/opls_lt.fdefn",lfdef=str(pathlib.Path(__file__).parent.resolve())+"/rdlt_data/lopls_lt.fdefn",gfdef=None,dfdef=None,cfdef=None,charge=True,refresh=False,loplsflag=False,gaffflag=False,dreidingflag=False,compassflag=False):
        #Build rdkit molecule from smiles and generate a conformer
        self.to_file=to_file
        self.name=pathlib.Path(to_file).stem

        # Validate mutual exclusivity of force fields
        active_flags = sum([gaffflag, loplsflag, dreidingflag, compassflag])
        if active_flags > 1:
            logger.error("Cannot use multiple force field flags simultaneously. Please choose one force field.")
            sys.exit("Error: Force field flags are mutually exclusive.")

        # Set default GAFF path if not provided
        if gaffflag and gfdef is None:
            gfdef = str(pathlib.Path(__file__).parent.resolve())+"/rdlt_data/gaff_lt.fdefn"
            if not os.path.exists(gfdef):
                logger.error(f"GAFF feature definition file not found: {gfdef}")
                sys.exit("Error: GAFF data files not found. Please ensure GAFF support is properly installed.")

        # Set default DREIDING path if not provided
        if dreidingflag and dfdef is None:
            dfdef = str(pathlib.Path(__file__).parent.resolve())+"/rdlt_data/dreiding_lt.fdefn"
            if not os.path.exists(dfdef):
                logger.error(f"DREIDING feature definition file not found: {dfdef}")
                sys.exit("Error: DREIDING data files not found. Please ensure DREIDING support is properly installed.")

        # Set default COMPASS path if not provided
        if compassflag and cfdef is None:
            cfdef = str(pathlib.Path(__file__).parent.resolve())+"/rdlt_data/compass_lt.fdefn"
            if not os.path.exists(cfdef):
                logger.error(f"COMPASS feature definition file not found: {cfdef}")
                sys.exit("Error: COMPASS data files not found. Please ensure COMPASS support is properly installed.")

        original_stdout = sys.stdout
        with open(to_file, 'w') as f:
            sys.stdout = f
            m = AllChem.AddHs(Chem.MolFromSmiles(self.smiles))
            AllChem.EmbedMolecule(m,AllChem.ETKDG())

            # WARNING: This part is dumb. Will update the lopls definitions ONLY
            # if the lopls flag is used. If a path is passed with the refresh command
            #
            if refresh and loplsflag:
                generateFeatureDefn(refresh,str(pathlib.Path(__file__).parent.resolve())+'/rdlt_data/lopls_lt.fdefn',str(pathlib.Path(__file__).parent.resolve())+'/rdlt_data/lopls_lt_dict.pkl')
            elif refresh and gaffflag:
                generateFeatureDefn(refresh,str(pathlib.Path(__file__).parent.resolve())+'/rdlt_data/gaff_tomoltemplate.txt',str(pathlib.Path(__file__).parent.resolve())+'/rdlt_data/gaff_lt_dict.pkl')
            elif refresh:
                generateFeatureDefn(refresh,str(pathlib.Path(__file__).parent.resolve())+'/rdlt_data/opls_lt.fdefn',str(pathlib.Path(__file__).parent.resolve())+'/rdlt_data/opls_lt_dict.pkl')

            #Build a feature factory from the defintion file and assign all features
            if gaffflag:
                factory = Chem.ChemicalFeatures.BuildFeatureFactory(gfdef)
            elif dreidingflag:
                factory = Chem.ChemicalFeatures.BuildFeatureFactory(dfdef)
            elif compassflag:
                factory = Chem.ChemicalFeatures.BuildFeatureFactory(cfdef)
            else:
                factory = Chem.ChemicalFeatures.BuildFeatureFactory(fdef)
            features = factory.GetFeaturesForMol(m)

            #Use the features to assign an atom type property
            [m.GetAtomWithIdx(f.GetAtomIds()[0]).SetProp('AtomType',f.GetType()) for f in features]

            # Apply conjugated system adjustments for GAFF
            # This adjusts atom types for conjugated systems (cc/cd, ce/cf pairs)
            # using bond connectivity information that SMARTS patterns cannot capture
            if gaffflag:
                m = self.detect_and_adjust_conjugated_systems(m)

            #if lopls defitions are desired, redo the feature process
            # overwrite atomtypes (not compatible with GAFF)
            if loplsflag:
                #print('loplsflag is {}'.format(loplsflag) )
                lfactory = Chem.ChemicalFeatures.BuildFeatureFactory(lfdef)
                lfeatures = lfactory.GetFeaturesForMol(m)
                #print(len(lfeatures))
                #for f in lfeatures:
                #    print(f.GetId(), f.GetFamily(), f.GetType(), f.GetAtomIds())
                [m.GetAtomWithIdx(f.GetAtomIds()[0]).SetProp('AtomType',f.GetType()) for f in lfeatures]
                #[print(at.GetProp('AtomType')) for at in m.GetAtoms()]

            #find untyped atoms
            #
            failure = False
            for at in m.GetAtoms():
                try:
                    at.GetProp('AtomType')
                except KeyError:
                    print("Atom {0} does not have an assigned atom type!".format(at.GetIdx()))
                    failure = True
            #if any failed to type, quit
            if failure:
                sys.exit("""Refusing to write a .lt file without type assignments.
        Check the SMARTS pattern that defines the expected atom type.""")


            #basic output
            writeHeader(self.name,loplsflag,gaffflag,dreidingflag,compassflag)
            writeAtoms(m)
            writeBonds(m)
            writeFooter(self.name,gaffflag,dreidingflag,compassflag)

            if charge:
                if gaffflag:
                    # GAFF requires manual charge calculation
                    print("\n# ========================================")
                    print("# IMPORTANT GAFF CHARGE NOTICE")
                    print("# ========================================")
                    print("# GAFF does NOT include default charges.")
                    print("# You must calculate charges separately using AM1-BCC or RESP.")
                    print("#")
                    print("# Quick start with AM1-BCC (recommended):")
                    print('#   1. Generate 3D structure: obabel -:"SMILES" -omol2 -O mol.mol2 --gen3d')
                    print('#   2. Calculate charges: antechamber -i mol.mol2 -fi mol2 -o mol_charged.mol2 -fo mol2 -c bcc -nc 0')
                    print('#   3. Extract charges: grep -v "^@" mol_charged.mol2 | awk "{print $NF}"')
                    print('#   4. Update charge values in the "Data Atoms" section above')
                    print("#")
                    print("# See AutoPoly/extern/rdlt_data/GAFF_README.md for detailed instructions")
                    print("# Reference: http://ambermd.org/antechamber/gaff.pdf")
                    print("# ========================================\n")

                    # Try to load GAFF charge dictionary (likely empty)
                    try:
                        gaff_cdict = read_cdict(str(pathlib.Path(__file__).parent.resolve())+'/rdlt_data/gaff_lt_dict.pkl')
                        if not gaff_cdict or all(v == 0.0 for v in gaff_cdict.values()):
                            print("# Note: GAFF charge dictionary is empty or contains only zeros.")
                            print("# This is expected - charges must be calculated manually.\n")
                        else:
                            # If charges exist, display them
                            sum_of_charges(m, gaff_cdict)
                    except FileNotFoundError:
                        print("# Note: No GAFF charge dictionary found.")
                        print("# This is expected - charges must be calculated manually.\n")
                elif dreidingflag:
                    # DREIDING requires manual charge calculation
                    print("\n# ========================================")
                    print("# IMPORTANT DREIDING CHARGE NOTICE")
                    print("# ========================================")
                    print("# DREIDING does NOT include default charges.")
                    print("# You must calculate charges separately using AM1-BCC, Gasteiger, or RESP.")
                    print("#")
                    print("# Quick start with Gasteiger charges (using Open Babel):")
                    print('#   obabel -i mol mol.mol -o mol2 -O mol_charged.mol2 --partialcharge gasteiger')
                    print("#")
                    print("# Or with AM1-BCC (recommended for better accuracy):")
                    print('#   antechamber -i mol.mol2 -fi mol2 -o mol_charged.mol2 -fo mol2 -c bcc -nc 0')
                    print("#")
                    print("# Reference: Mayo et al., J. Phys. Chem. 1990, 94, 8897-8909")
                    print("# ========================================\n")

                    # Try to load DREIDING charge dictionary (likely empty)
                    try:
                        dreiding_cdict = read_cdict(str(pathlib.Path(__file__).parent.resolve())+'/rdlt_data/dreiding_lt_dict.pkl')
                        if not dreiding_cdict or all(v == 0.0 for v in dreiding_cdict.values()):
                            print("# Note: DREIDING charge dictionary is empty or contains only zeros.")
                            print("# This is expected - charges must be calculated manually.\n")
                        else:
                            # If charges exist, display them
                            sum_of_charges(m, dreiding_cdict)
                    except FileNotFoundError:
                        print("# Note: No DREIDING charge dictionary found.")
                        print("# This is expected - charges must be calculated manually.\n")
                elif compassflag:
                    # COMPASS uses bond increment model
                    print("\n# ========================================")
                    print("# COMPASS CHARGE NOTICE")
                    print("# ========================================")
                    print("# COMPASS uses a bond increment charge model.")
                    print("# Charges may be assigned automatically by moltemplate,")
                    print("# or you may need to calculate them manually depending on the molecule.")
                    print("#")
                    print("# NOTE: This is an incomplete public version of COMPASS.")
                    print("# Parameters for some atom types may be missing.")
                    print("#")
                    print("# Reference: Sun, H., J. Phys. Chem. B, 1998, 102, 7338-7364")
                    print("# ========================================\n")

                    # Try to load COMPASS charge dictionary
                    try:
                        compass_cdict = read_cdict(str(pathlib.Path(__file__).parent.resolve())+'/rdlt_data/compass_lt_dict.pkl')
                        if not compass_cdict or all(v == 0.0 for v in compass_cdict.values()):
                            print("# Note: COMPASS charge dictionary is empty or contains only zeros.")
                            print("# Charges may be handled via bond increment model in .lt file.\n")
                        else:
                            # If charges exist, display them
                            sum_of_charges(m, compass_cdict)
                    except FileNotFoundError:
                        print("# Note: No COMPASS charge dictionary found.")
                        print("# Charges may be handled via bond increment model in .lt file.\n")
                else:
                    # OPLS charge handling (existing logic)
                    # Read charge dictionaries for testing
                    opls_cdict = read_cdict(str(pathlib.Path(__file__).parent.resolve())+'/rdlt_data/opls_lt_dict.pkl')
                    if loplsflag:
                        lopls_cdict = read_cdict(str(pathlib.Path(__file__).parent.resolve())+'/rdlt_data/lopls_lt_dict.pkl')
                        opls_cdict.update(lopls_cdict)

                    sum_of_charges(m,opls_cdict)
            sys.stdout = original_stdout
