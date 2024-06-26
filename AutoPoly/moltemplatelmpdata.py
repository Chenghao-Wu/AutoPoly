#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 21 12:19:08 2018

@author: zwu
"""
import sys
import os
from pathlib import Path

import numpy as np
from .system import logger

class Moltemplatelmpdata(object):
    def __init__(self, System=None, Model=None, run=True, path_MonomerBank=None):
        self.System=System
        self.path_cwd=self.System.get_FolderPath+"/"+self.Name+"/moltemplate/"
        self.path_master=str(Path(__file__).parent.resolve())+"/extern/"
        logger.info(' '.join(["\n'you are now using extern path of "+self.path_master+"\n"]))
        self.path_moltemplatesrc=self.path_master+"moltemplate/src/"
        self.path_oplsaaprm=self.path_master+"moltemplate/oplsaa.prm"
        if path_MonomerBank==None:
            self.path_MonomerBank=self.path_master+"Monomer_bank/"
        else:
            self.path_MonomerBank=path_MonomerBank

        self.Model=Model

        self.rotate=90.0
        self.offset_spacing=2.0
        self.offset=4.0
        self.packingL_spacing=5.0
        self.moltemplateBoxSize=400.0 # OPTIMIZED: will change accord. to actual packing 

        self.seuqenceSet = []
        self.chem_mono = []

        self.create_Folder()

        if run==True:
            self.make_LmpDataFilebyMoltemplate()

    def create_Folder(self):
        if self.System.MadeFolder==True:
            path=Path(self.path_cwd)
            if path.exists():
                logger.error(' '.join([self.path_cwd," already exist! Please have a check"]))
                sys.exit()
            else:
                path.mkdir(parents=True, exist_ok=True)

    def set_tacticity(self,tacticity):
        self.tacticity=tacticity

    def make_atacticSequenceSet(self):
        SequenceSet = []
        return SequenceSet


    def n_monomerAtoms(self,merltfile):
        n_monomerAtoms=0
        MonomerBank=Path(self.path_MonomerBank)
        merltfile_Path=MonomerBank/merltfile
        if merltfile_Path.is_file():
            is_inside_block=False
            with open(merltfile_Path) as f:
                while True: 
                    line = f.readline() 
                    if line.strip()=="write(\"Data Atoms\") {":
                        is_inside_block=True
                        line = f.readline() 
                    elif line.strip()=="}":
                        is_inside_block=False
                    
                    if is_inside_block:
                        n_monomerAtoms=n_monomerAtoms+1

                    if not line: 
                        break
        else:
            logger.error(' '.join(["in MoltemplateLmpData::n_monomerAtoms():\n",merltfile_Path," file cannot open.\n"]))
            sys.exit()
        return n_monomerAtoms
    
    def make_SequenceSet_adhoc(self):
        sequenceSet = []
        merSet = []
        counter=0
        n_chainAtoms=0
        n_chains=0
        chain=0
        mer=0

        for indexi in range(len(self.chem_mono)):
            counter += self.mono_beads[indexi]
            if self.sequenceLen>1:
                n_chainAtoms += self.n_monomerAtoms(self.chem_mono[indexi]+"le.lt")*self.mono_beads[indexi]
            else:
                n_chainAtoms += self.n_monomerAtoms(self.chem_mono[indexi]+".lt")*self.mono_beads[indexi]
        # /** sanity check **/
        if counter!=self.sequenceLen:
            logger.error(' '.join(["in MoltemplateLmpData::make_SequenceSet_adhoc():\nNumber of mers per chain",str(counter),"!= sequenceLen (",str(self.sequenceLen),"). Please check.\n"]))
            sys.exit()
        # /** Determine number of chains **/
        if self.n_total_atoms==0: #{/** controlled by sequenceNum **/
            n_chains=self.sequenceNum
        else: #{/** controlled by n_total_atoms **/
            if n_chainAtoms>self.n_total_atoms:
                logger.error(' '.join(["in MoltemplateLmpData::make_SequenceSet_adhoc():\nn_chainAtoms",str(n_chainAtoms),") > n_total_atoms(",str(self.n_total_atoms),")\nThe upper threshold of total number of atoms in system is reached.\nPlease either increase n_total_atoms or choose a shorter chain.\n"]))
                sys.exit()
            n_chains=np.floor(self.n_total_atoms/n_chainAtoms)
        # /** sanity check **/
        if n_chains<=0:
            logger.error(' '.join(["in MoltemplateLmpData::make_SequenceSet_adhoc():\nn_chain=",str(n_chains),", which is invalid. Please check.\n"]))
            sys.exit()

        #/** For generating small molecule **/
        # //--------------------------------------------------------------------------
        if self.sequenceLen==1:
        
            if self.chem_mono.size()>1:
                logger.error(' '.join(["sequenceLen = ",str(self.sequenceLen),"merSet should only have one mer type! Please check.\n"]))
                sys.exit()
            for chain in range(n_chains):
                merSet = []
                for mer in range(self.sequenceLen):
                    merSet.append(self.chem_mono[0]+".lt")
                self.sequenceSet.append(merSet)
            return sequenceSet
        #/** Homopolymer **/
        import random

        if len(self.chem_mono) == 1:
            if bool(random.getrandbits(1)):
                chosenTac="_T1.lt"
            else:
                chosenTac=".lt"

            for chain in range(n_chains):
                merSet = []
                
                if self.tacticity=="atactic": #/** random chirality **/
                    for mer in range(self.sequenceLen):
                        if mer==0:# {//1st monomer: le ltfile
                            if bool(random.getrandbits(1)):# {//true: T1
                                merSet.append(self.chem_mono[0]+"le_T1.lt")
                            else:# {//false: normal
                                merSet.append(self.chem_mono[0]+"le.lt")
                        elif mer==(self.sequenceLen-1): #{//last monomer: re ltfile
                            if bool(random.getrandbits(1)):# {//true: T1
                                merSet.append(self.chem_mono[0]+"re_T1.lt")
                            else:# {//false: normal
                                merSet.append(self.chem_mono[0]+"re.lt")
                        else: #{//intermediate monomers: i ltfiles
                            if bool(random.getrandbits(1)):# {//true: T1
                                merSet.append(self.chem_mono[0]+"i_T1.lt")
                            else:# {//false: normal
                                merSet.append(self.chem_mono[0]+"i.lt")

                elif self.tacticity=="isotactic":# /** monotonic chirality **/
                    for mer in range(self.sequenceLen):
                        if mer==0: #{//1st monomer: le ltfile
                            merSet.append(self.chem_mono[0]+"le"+chosenTac)
                        elif mer==(self.sequenceLen-1):# {//last monomer: re ltfile
                            merSet.append(self.chem_mono[0]+"re"+chosenTac)
                        else:# {//intermediate monomers: i ltfiles
                            merSet.append(self.chem_mono[0]+"i"+chosenTac)

                elif self.tacticity=="syndiotactic":# /** alternating chirality **/
                    #string startTac,nextTac,currentTac;
                    #/** coin flip to decide start tacticity **/
                    if bool(random.getrandbits(1)):
                        startTac="_T1.lt"
                        nextTac =".lt"
                    else: 
                        startTac=".lt"
                        nextTac ="_T1.lt"
                    for mer in range(self.sequenceLen):
                        if mer%2==0:
                            currentTac=startTac
                        else:
                            currentTac=nextTac
                        if mer==0: #{//1st monomer: le ltfile
                            merSet.append(self.chem_mono[0]+"le"+currentTac)
                        elif mer==(self.sequenceLen-1):# {//last monomer: re ltfile
                            merSet.append(self.chem_mono[0]+"re"+currentTac)
                        else:# {//intermediate monomers: i ltfiles
                            merSet.append(self.chem_mono[0]+"i"+currentTac)
                else:
                    logger.error(' '.join(["in MoltemplateLmpData::make_SequenceSet_adhoc()\nPlease choose tacticity from {atactic,syndiotactic,isotactic}\n"]))
                    sys.exit()
                #/** load in tacticity of current chain **/
                sequenceSet.append(merSet)
        else:
            if self.copolymerType=="random":
                #//TODO
                return
            
            if self.copolymerType=="block":# { /** Block Copolymer **/

                if bool(random.getrandbits(1)):
                    chosenTac="_T1.lt"
                else:
                    chosenTac=".lt"

                for chain in range(n_chains):
                    merSet = []
                    counter_mer=0
                    if self.tacticity=="atactic": #/** random chirality **/
                        for mertype in range(len(self.chem_mono)):
                            for mer in range(self.sequenceLen):
                                if counter_mer==0:# {//1st monomer: le ltfile
                                    if bool(random.getrandbits(1)):# {//true: T1
                                        merSet.append(self.chem_mono[mertype]+"le_T1.lt")
                                    else:# {//false: normal
                                        merSet.append(self.chem_mono[mertype]+"le.lt")
                                elif counter_mer==(self.sequenceLen-1): #{//last monomer: re ltfile
                                    if bool(random.getrandbits(1)):# {//true: T1
                                        merSet.append(self.chem_mono[mertype]+"re_T1.lt")
                                    else:# {//false: normal
                                        merSet.append(self.chem_mono[mertype]+"re.lt")
                                else: #{//intermediate monomers: i ltfiles
                                    if bool(random.getrandbits(1)):# {//true: T1
                                        merSet.append(self.chem_mono[mertype]+"i_T1.lt")
                                    else:# {//false: normal
                                        merSet.append(self.chem_mono[mertype]+"i.lt")
                                counter_mer=counter_mer+1

                    elif self.tacticity=="isotactic":# /** monotonic chirality **/
                        for mertype in range(len(self.chem_mono)):
                            for mer in range(self.sequenceLen):
                                if counter_mer==0: #{//1st monomer: le ltfile
                                    merSet.append(self.chem_mono[mertype]+"le"+chosenTac)
                                elif counter_mer==(self.sequenceLen-1):# {//last monomer: re ltfile
                                    merSet.append(self.chem_mono[mertype]+"re"+chosenTac)
                                else:# {//intermediate monomers: i ltfiles
                                    merSet.append(self.chem_mono[mertype]+"i"+chosenTac)
                                counter_mer=counter_mer+1

                    elif self.tacticity=="syndiotactic":# /** alternating chirality **/
                        #string startTac,nextTac,currentTac;
                        #/** coin flip to decide start tacticity **/
                        if bool(random.getrandbits(1)):
                            startTac="_T1.lt"
                            nextTac =".lt"
                        else: 
                            startTac=".lt"
                            nextTac ="_T1.lt"
                        for mertype in range(len(self.chem_mono)):
                            for mer in range(self.sequenceLen):
                                if counter_mer%2==0:
                                    currentTac=startTac
                                else:
                                    currentTac=nextTac
                                if counter_mer==0: #{//1st monomer: le ltfile
                                    merSet.append(self.chem_mono[mertype]+"le"+currentTac)
                                elif counter_mer==(self.sequenceLen-1):# {//last monomer: re ltfile
                                    merSet.append(self.chem_mono[mertype]+"re"+currentTac)
                                else:# {//intermediate monomers: i ltfiles
                                    merSet.append(self.chem_mono[mertype]+"i"+currentTac)
                                counter_mer=counter_mer+1
                    else:
                        logger.error(' '.join(["in MoltemplateLmpData::make_SequenceSet_adhoc()\nPlease choose tacticity from {atactic,syndiotactic,isotactic}\n"]))
                        sys.exit()
                    #/** load in tacticity of current chain **/
                    sequenceSet.append(merSet)

        return sequenceSet
    
    def getridof_ljcutcoullong(self):
        in_=self.path_cwd+"system.in.settings"
        out=self.path_cwd+"tmp.data"
        write_f=open(out, "w")
        with open(in_,'r') as read_f:
            
            in_path=Path(in_)
            if in_path.is_file():
                
                while True:
                    line = read_f.readline() 
                    if line.strip()=="":
                        write_f.write("\n")
                    elif line.strip().split()[0]=="pair_coeff":
                        #write_f.write("    pair_coeff ")
                        space_i=0
                        for ii in line.split():
                            if ii=="lj/cut/coul/long":
                                continue
                            else:
                                if space_i==0:
                                    write_f.write("    ")
                                    space_i=space_i+1
                                write_f.write(ii+" ")
                        write_f.write("\n")
                    else:
                        write_f.write(line)
                        
                    if not line: 
                        break
            else:
                logger.error(' '.join(["system.in.setting does not exist plase check ",in_]))
                sys.exit()
        write_f.close()
        mv="rm "+in_+";mv "+out+" "+in_
        os.system(mv)

    def make_LmpDataFilebyMoltemplate(self):
        logger.info(' '.join(["\nlmpdata prepared by Moltemplate.\n"]))
        
        poly_index=0
        for modelii in self.Model:
            logger.info(' '.join(["\nNumber of Molecules = ",str(len(modelii.sequenceSet))]))
            # loop through all polymers and make corresponding polymer lt files
            for moleii in range(len(modelii.sequenceSet)):
                # check degrees of polymerization (DOP) of current polymer
                if modelii.DOP>0:
                    if len(modelii.sequenceSet[moleii])!=modelii.DOP:
                        logger.warning(' '.join(["\nWarning: At molecule# ",str(moleii)," DOP=",str(len(modelii.sequenceSet[moleii])),
                                                " != ",str(modelii.DOP),"\n"]))
                else:
                    logger.error(' '.join(["\nWarning: At molecule#",str(moleii+1),",",
                                            "DOP=",str(len(modelii.sequenceSet[moleii])),str(modelii.DOP)]))

                # check monomer.lt's in the monomer bank
                for indexii in range(len(modelii.sequenceSet[moleii])):
                    if self.check_monomerbank(modelii.sequenceSet[moleii][indexii]):
                        source=Path(self.path_MonomerBank)/modelii.sequenceSet[moleii][indexii]
                        self.copy_to_cwd(source)
                    else:
                        logger.error(' '.join(["\nError: ","At monomer#" +str(indexii+1)+" ("+modelii.sequenceSet[moleii][indexii]+") ","of molecule#"+str(moleii+1) +": ",
                                            "\nCan't find corresponding lt file in the monomer bank.\n\n",
                                            "Automation terminated.\n\n"]))
                        sys.exit()

                if modelii.DOP>1:
                    # make poly.lt file
                    self.make_polylt(poly_index,modelii.sequenceName[moleii])
                    poly_index=poly_index+1

        #make oplsaa.lt (requiring oplsaa_subset.prm)
        #NOTE: the oplsaa.lt is shared by the current polymer and all its
        #constituent monomers; to make it more general, this function should
        #generate an unique oplsaa.lt file for each polymer and its own
        #monomers. This feature is NOT supported in the current version
    
        self.make_oplsaalt()

        # make system.lt file
        self.make_systemlt()

        # invoke moltemplate to generate LAMMPS datafile
        self.invoke_moltemplate()

        self.getridof_ljcutcoullong()

        # move files to working directory
        self.mv_files()

        logger.info("\nmoltemplate-enabled lmpdata made!\n")


    def mv_files(self):
        datafile="system.data"
        incharge="system.in.charges"
        insetting="system.in.settings"
        data="cd "+self.path_cwd+";"+"cp "+datafile+" ../; cd ../;"
        os.system(data)
        init="cd "+self.path_cwd+";"+"cp "+incharge+" "+insetting+" system.in system.in.init ../;"+"cd ..;"
        os.system(init)
        output="cd "+self.path_cwd+";"
        output+="mkdir output; mv system.in* system*data output_ttree output/"
        os.system(output)
        input_="cd "+self.path_cwd+"; mkdir input; mv *.lt *.prm input/"
        os.system(input_)

    def evaluate_boxLen(self):
        in_=path_cwd+"system.data"
        dubVar=0
        lmin=0
        lmax=0

    def invoke_moltemplate(self):
        # NOTE: system.lt is in cwd
        bash="cd "+self.path_cwd+"; "+self.path_moltemplatesrc+"moltemplate.sh ./system.lt"
        os.system(bash)

    def make_systemlt(self):
        output=self.path_cwd+"/system.lt"
        with open(output, "w") as write_f:
            polyindex=0
            for modelii in self.Model:
                n_poly=len(modelii.sequenceSet)
                if modelii.DOP>1:
                    for indexi in range(n_poly):
                        write_f.write("import \"poly_"+str(polyindex+1)+".lt\"\n")
                        polyindex=polyindex+1
                    write_f.write("\n")
                else:
                    if len(modelii.merSet)>1:
                        logger.error(' '.join(["sequenceLen = "+str(modelii.DOP)+", "
                                                    , " merSet should only have one mer type! Please check.\n"]))
                        sys.exit()
                    #import constituent monomer.lt's
                    unique_Sequence=[i[0] for i in modelii.sequenceSet]
                    print(unique_Sequence)
                    for sequenceii in range(len(unique_Sequence)):
                        write_f.write("import \""+unique_Sequence[sequenceii]+"\"\n")
                    write_f.write("\n")

            polyindex=0

            index=0

            packingL=self.offset+self.packingL_spacing
            
            offset_x=-50

            offset_x_increment=10.0
            
            first=True
            coordinates=[]
            
            #max_DOP=max([i.DOP for i in self.Model])

            for modelii in self.Model:
                counter=0
                n=0
                bndl=0
                bndh=0
                n_now=0
                n_pre=0
                signy=1
                signz=-1
                timey=0
                timez=0
                valy=0
                valz=0
                n_poly=len(modelii.sequenceSet)

                # Pack molecules in square spiral shape
                if modelii.DOP>1:
                    offset_x+=offset_x_increment
                    write_f.write("polymer_"+str(index+1)+" = new poly_"+str(polyindex+1)+".move("+"{:.4f}".format(offset_x)+","+"{:.4f}".format(valy)+","+"{:.4f}".format(valz)+")"+ "\n")
                    offset_x_increment=self.offset*(modelii.DOP+2)
                    polyindex+=1
                else:
                    #packingL=10
                    offset_x+=offset_x_increment
                    write_f.write("molecule_"+str(index+1)+" = new "+modelii.merSet[0]+".move("+"{:.4f}".format(offset_x)+","+"{:.4f}".format(valy)+","+"{:.4f}".format(valz)+")"+ "\n")
                    offset_x_increment=self.offset*(modelii.DOP+2)

                index=index+1

                for indexi in range(1,n_poly):
                    n=0
                    while True:
                        n=n+1
                        bndl=(n-1)*n
                        bndh=n*(n+1)
                        if bndl<index and indexi<=bndh:
                            break
                    n_now=n
                    if n_now!=n_pre:
                        counter=0
                        signy*=-1
                        signz*=-1
                    if counter<n_now:
                        timey=1
                        valy+=packingL*signy*timey
                        timez=0
                    else:
                        timey=0
                        timez=1 
                        valz+=packingL*signz*timez
                    if modelii.DOP>1:
                        #offset_x=-50.0
                        write_f.write("polymer_"+str(index+1)+" = new "+"poly_"+str(polyindex+1)+".move("+"{:.4f}".format(offset_x)+","+"{:.4f}".format(valy)+","+"{:.4f}".format(valz)+")"+ "\n")
                        polyindex+=1
                    else:
                        #offset_x=-60.0
                        write_f.write("molecule_"+str(index+1)+" = new "+modelii.merSet[0]+".move("+"{:.4f}".format(offset_x)+","+"{:.4f}".format(valy)+","+"{:.4f}".format(valz)+")"+ "\n")
                    n_pre=n_now
                    counter+=1
                    index=index+1
                    

                write_f.write("\n")
                #first=True

            hbox=self.moltemplateBoxSize*0.5
            fbox=self.moltemplateBoxSize

            if True:
                write_f.write("write_once(\"Data Boundary\") {\n")
                write_f.write("   -"+str(hbox)+"  "+str(hbox)+"  xlo xhi\n")
                write_f.write("   -"+str(hbox)+"  "+str(hbox)+"  ylo yhi\n")
                write_f.write("   -"+str(hbox)+"  "+str(hbox)+"  zlo zhi\n")
                write_f.write("}\n")
                write_f.write("\n")
            else:
                write_f.write("write_once(\"Data Boundary\") {\n")
                write_f.write("   0.0  "+fbox+"  xlo xhi\n")
                write_f.write("   0.0  "+fbox+"  ylo yhi\n")
                write_f.write("   0.0  "+fbox+"  zlo zhi\n")
                write_f.write("}\n")
                write_f.write("\n")
            write_f.close()



    def make_polylt(self,polyindex,monomerSet):
        output=self.path_cwd+"/poly_"+str(polyindex+1)+".lt"

        with open(output, "w") as write_f:
            write_f.write("import \"oplsaa.lt\"\n")


            unique_monomers=list(dict.fromkeys(monomerSet))
            for monomerii in range(len(unique_monomers)):
                write_f.write("import \""+unique_monomers[monomerii]+".lt"+"\"\n")

            write_f.write("\n")

            #Define combined molecule (ex.polymer)
            write_f.write("poly_"+str(polyindex+1)+" inherits OPLSAA {\n\n")
            write_f.write("    "+ "create_var {$mol}\n\n")

            monomerSet_copy=monomerSet
            offset_cum=0
            for indexii in range(len(monomerSet)):
                # erase .lt from name string
                #del monomerSet_copy[-3:]
                # pack monomers along x-axis and rotate accordingly (1,0,0)
                write_f.write("    "+"monomer["+str(indexii)+"] = new "+monomerSet[indexii])
                if indexii>0:
                    write_f.write(".rot(" +str(self.rotate*(indexii%2))+",1,0,0)"+".move("+"{:.4f}".format(offset_cum)+",0,0)")
                write_f.write("\n")
                # evaluate offset distance based on C1-C2 of the pre-mer

                self.evaluate_offset(monomerSet[indexii]+".lt")
                offset_cum+=self.offset
                #print(offset_cum)
            # add a list of bonds connecting propagating carbons
            write_f.write("\n    write('Data Bond List') {\n")
            for indexii in range(len(monomerSet)-1):
                write_f.write("      "+"$bond:b"+str(indexii+1)+"  "+"$atom:monomer["+str(indexii)+"]/C2"+"  "+"$atom:monomer["+str(indexii+1)+"]/C1"+"  "+"\n")
            write_f.write("    }\n")
            # end cap of poly.lt scope
            write_f.write("\n} # poly_"+str(polyindex+1)+ "\n") 
            write_f.close()
    

    def evaluate_offset(self,merltfile):
        MonomerBank=Path(self.path_MonomerBank)
        merltfile_Path=MonomerBank/merltfile
        if merltfile_Path.is_file():
            C1=[]
            C2=[]
            dubVar=0
            with open(merltfile_Path) as f:
                while True: 
                    line = f.readline() 
                    if line.strip()=="write(\"Data Atoms\") {":
                        # C1 coordinates
                        line = f.readline() 
                        for i in range(3):
                            C1.append(float(line.split()[i+4]))
                        # C2 coordinates
                        line = f.readline() 
                        for i in range(3):
                            C2.append(float(line.split()[i+4]))
                        
                        # calculate C1-C2 distance
                        self.offset=np.linalg.norm(np.array(C1)-np.array(C2))+self.offset_spacing
                        
                        return
            
    def make_oplsaalt(self):
        
        self.make_oplsaa_subset()

        # invoke oplsaa_moltemplate.py to make oplsaa.lt 
        oplsaa_subset=self.path_cwd+"oplsaa_subset.prm"
        oplsaa_py=self.path_moltemplatesrc+"oplsaa_moltemplate.py "+oplsaa_subset
        bash="cd "+self.path_cwd+"; "+oplsaa_py
        os.system(bash)
    
    def make_oplsaa_subset(self):
        # path to oplsaa_subset.prm file
        opls_subset_file = self.path_cwd+"oplsaa_subset.prm"

        atom_keys=[]
        for modelii in self.Model:
            for monomerii in range(len(modelii.sequenceSet)):
                monomerSet=modelii.sequenceSet[monomerii]
                # vector to store all atom types including the repeats
                for vecii in range(len(monomerSet)):
                    # path to monomer.lt in monomer bank
                    MonomerBank=Path(self.path_MonomerBank)
                    merltfile_Path=MonomerBank/monomerSet[vecii]
                    if merltfile_Path.is_file():
                        mono = self.path_MonomerBank+monomerSet[vecii]
                        read_switch= False 
                        with open(mono) as f:
                            while True: 
                                line = f.readline() 
                                
                                if line.strip()=="write(\"Data Atoms\") {":
                                    read_switch=True
                                    continue
                                elif line.strip()=="}":
                                    read_switch=False
                                    break
                                
                                # Determine atom types, element names and the raw_charges as
                                # given in the opls table
                                
                                if read_switch:
                                    load_line=""
                                    stringvector=line.split()
                                    
                                    load_switch=False
                                    for readii in range(len(stringvector[2])):
                                        if stringvector[2][readii]==":":
                                            load_switch = True
                                            continue
                                        if load_switch:
                                            load_line += stringvector[2][readii]
                                    atom_keys.append(load_line)
                                
                                if not line: 
                                    break
                    else:
                        logger.error(' '.join(["Monomer ("+ monomerSet[vecii] + ") does NOT exist. \n",
                                                "Please check the following path to the file\n" + merltfile_Path + "\n"]))
                        sys.exit()
                        
        # Cleaning up the stored data. Remove duplicate atoms types
        atom_types=list(dict.fromkeys(atom_keys))
        #print(atom_types)
        # Convert the vectors string to vector int in order to sort the atom_types in ascending order
        atom_types=sorted([int(i) for i in atom_types])
        # Read the master opls file and store the ones that match the atom_types into new subset file
        write_f=open(opls_subset_file, "w")
        
        with open(self.path_oplsaaprm,'r') as read_f:
            path_oplsaaprm=Path(self.path_oplsaaprm)
            if path_oplsaaprm.is_file():
                check_switch=False
                while True:
                    prm_line = read_f.readline()
                    if len(prm_line.strip())!=0:
                       
                        if prm_line.strip() == "##  Atom Type Definitions  ##":
                            check_switch = True
                            write_f.write(prm_line+"\n")
                            prm_line = read_f.readline()
                            write_f.write(prm_line+"\n")
                            prm_line = read_f.readline()
                            write_f.write(prm_line+"\n")
                            continue
                        elif prm_line.strip()=="################################":
                            check_switch = False
                            write_f.write(prm_line+"\n")
                            continue
                        elif check_switch:
                            
                            stringvector=prm_line.split()
                            
                            
                            for checkii in range(len(atom_types)):
                                
                                if atom_types[checkii]==int(stringvector[1]):
                                    write_f.write(prm_line+"\n")
                                    break
                        else:
                            write_f.write(prm_line+"\n")
                    else:
                        write_f.write(prm_line+"\n")

                    if not prm_line: 
                        break
        write_f.close()

    def check_monomerbank(self,monomer):
        monomer_path=Path(self.path_MonomerBank)/monomer
        return monomer_path.is_file()

    def copy_to_cwd(self,source):
        bash="cp "
        bash=bash+str(source)+" "+self.path_cwd
        os.system(bash)

