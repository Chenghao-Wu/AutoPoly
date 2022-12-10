import sys
sys.path.append('/home/zwq2834/development/AutoPoly')
import AutoPoly

rdlt=AutoPoly.RDlt(smiles='c1c2ccccc2ccc1')
rdlt.run(to_file='Naphthalene.lt')
rdlt.store_bank()

# Define the system
# out is the folder name for the output
system=AutoPoly.System(out="test_molecule_from_RDkit")

# create polymers
# Just an example of polypropylene and polyethylene with all-atom and united-atom resolutions

naphthalene=AutoPoly.Polymer(ChainNum=50,Sequence=["Naphthalene"])

# polymerization
# Name is the output folder for this polymer
poly=AutoPoly.Polymerization(Name="Naphthalene",System=system,Model=[naphthalene],run=True)