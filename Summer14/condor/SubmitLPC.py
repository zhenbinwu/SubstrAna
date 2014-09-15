#!/usr/bin/env python
# coding: utf-8

import os
import re
import shutil
import glob

DelDir    = '/uscms_data/d3/benwu/JetMET_TP/Puppi/CMSSW_6_2_0_SLHC16_patch1'
DelExe    = 'runPuppi'
Directory = '/eos/uscms/store/user/benwu/TP/Puppi/Sep12'
UserEMAIL = 'benwu@fnal.gov'
Fraction  = 1
Detectors = [
    'PhaseI',
    #'PhaseII3',
    #'PhaseII4'
]
PileUps   = [
    '50PU_noaged',
    '140PU_aged',
    #'50PileUp',
    #'140PileUp',
]
Projects  = {
    'QCD':[1, 0, 0],
    'DYMM':[0, 1, 1],
}


def Condor_Sub():
    for dec in Detectors:
        for pu in PileUps:
            for pro in Projects.keys():
                for splitpro in SplitPro(dec, pu, pro):
                    cond_file = splitpro.split('/')[-1].split('.')[0]
                    cond_file += "_condor"
                    print cond_file
                    filelist = ','.join("root://cmseos.fnal.gov:1094/"+x.strip() for x in open(splitpro, "r"))
                    print filelist
                    with open(cond_file, "wt") as out:
                        for line in open("Puppi_condor", "r"):
                            line = line.replace("USER@FNAL.GOV", UserEMAIL)
                            if dec == "PhaseI" and pu == "50PU_noaged":
                                line = line.replace("TAG", "DES19_V1_MC")
                            elif dec == "PhaseI" and pu == "140PU_aged":
                                line = line.replace("TAG", "AGE1K_V1_MC")
                            line = line.replace("GEN", str(Projects[pro][0]))
                            line = line.replace("SelZ", str(Projects[pro][1]))
                            line = line.replace("MET", str(Projects[pro][2]))
                            line = line.replace("List", filelist)
                            out.write(line)


                    os.system("condor_submit " + cond_file)
                    #return


def my_CheckFile():
    ## Check the Delphes Dir
    if os.path.isdir(DelDir):
        pass
    else:
        print "Please input the path to Delphes"
        quit()


    ## Check RunHT.csh file
    if os.path.isfile("RunHT.csh") and os.access("RunHT.csh", os.X_OK):
        #print "Found RunHT.csh"
        pass
    else:
        print "Please locate RunHT.csh"
        quit()

    ## Check DelFill to be execute
    DelFill = DelDir + "/" + DelExe
    if os.path.isfile(DelFill) and os.access(DelFill, os.X_OK):
        #print "Found DelFill"
        pass
    else:
        print "Please locate %s" % DelFill
        quit()

    ## Check HTadd
    if os.path.isfile("HTadd") and os.access("HTadd", os.X_OK):
        #print "Found HTadd"
        pass
    else:
        print "Please compile HTadd"
        return None

    ## Check Delphes_condor
    if os.path.isfile("Delphes_condor"):
        #print "Found HTadd"
        pass
    else:
        print "Please compile HTadd"
        quit()

def SplitPro(detector, pileup, pro):
    globout=glob.glob('./FileList/%s_%s_%s_*.list'  % (pro, detector,pileup))
    return sorted(globout[:len(globout)/Fraction])
    testout=[]
    for out in globout:
        file = out.split('/')[-1]
        testout.append(file.split(".")[0])
    return testout

def my_process():
    ## Create the output directory
    outdir = Directory + "/condor"
    try:
        os.makedirs(outdir)
    except OSError:
        pass

    ## Update RunHT.csh with DelDir and pileups
    RunHTFile = outdir + "/" + "RunPuppi.csh"
    with open(RunHTFile, "wt") as outfile:
        for line in open("RunPuppi.csh", "r"):
            line = line.replace("DELDIR", DelDir)
            line = line.replace("DELEXE", DelExe)
            outfile.write(line)

    ## Update condor files
    shutil.copy2("Puppi_condor", outdir)
    os.system("tar -czf %s/FileList.tgz FileList" % outdir )
    os.chdir(outdir)
    os.system("tar -xzvf FileList.tgz")
    Condor_Sub()


if __name__ == "__main__":
    my_process()
