
#! /usr/bin/env python

from __future__ import print_function
import sys, pickle, argparse, os, subprocess, time

import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem

import pathlib

from ..system import logger


def writeHeader(molname, loplsflag):
    print("""import "oplsaa.lt"    # <-- defines the standard "OPLSAA" force field""")
    if loplsflag:
        print("""import "loplsaa.lt"   # <-- custom parameters for long alkane chains taken from
                      #     Sui et al. J.Chem.Theory.Comp (2012), 8, 1459
                      #     To use the ordinary OPLSAA force field parameters,
                      #     (instead of the Sui et al. parameters), change the
                      #     atom types below from "@atom:81L","@atom:85LCH2" to
                      #     "@atom:81" and "@atom:85"  (defined in "oplsaa.lt")""")
    print("""{0} inherits OPLSAA {{""".format(molname))


def writeFooter(molname):
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


    def run(self,to_file=None,name='test',fdef=str(pathlib.Path(__file__).parent.resolve())+"/rdlt_data/opls_lt.fdefn",lfdef=str(pathlib.Path(__file__).parent.resolve())+"/rdlt_data/lopls_lt.fdefn",charge=True,refresh=False,loplsflag=False):
        #Build rdkit molecule from smiles and generate a conformer
        self.to_file=to_file
        self.name=pathlib.Path(to_file).stem


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
            elif refresh:
                generateFeatureDefn(refresh,str(pathlib.Path(__file__).parent.resolve())+'/rdlt_data/opls_lt.fdefn',str(pathlib.Path(__file__).parent.resolve())+'/rdlt_data/opls_lt_dict.pkl')

            #Build a feature factory from the defintion file and assign all features
            factory = Chem.ChemicalFeatures.BuildFeatureFactory(fdef)
            features = factory.GetFeaturesForMol(m)

            #Use the features to assign an atom type property
            [m.GetAtomWithIdx(f.GetAtomIds()[0]).SetProp('AtomType',f.GetType()) for f in features]

            #if lopls defitions are desired, redo the feature process
            # overwrite atomtypes
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
            writeHeader(self.name,loplsflag)
            writeAtoms(m)
            writeBonds(m)
            writeFooter(self.name)

            if charge:
                # Read charge dictionaries for testing
                opls_cdict = read_cdict(str(pathlib.Path(__file__).parent.resolve())+'/rdlt_data/opls_lt_dict.pkl')
                if loplsflag:
                    lopls_cdict = read_cdict(str(pathlib.Path(__file__).parent.resolve())+'/rdlt_data/lopls_lt_dict.pkl')
                    opls_cdict.update(lopls_cdict)

                sum_of_charges(m,opls_cdict)
            sys.stdout = original_stdout
