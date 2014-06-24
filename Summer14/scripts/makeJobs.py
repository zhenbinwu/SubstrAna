import os
import sys
import numpy
import glob

from optparse import OptionParser
parser = OptionParser()
parser.add_option("-i","--inputdir"   , dest="inputdir"   , type="string", help="Path where ntuples are located. Example: /store/group/phys_jetmet/ntran/PUPPI/miniSamples/62x/qcd300-470_62x_PU40BX50")
parser.add_option("-w","--workdir"    , dest="workdir"    , type="string", default="mydir",help="Name of the directory for jobs")
parser.add_option("-o","--outputname" , dest="outputname" , type="string", default="outtree",help="Name of the output file. Default is: outtree")
parser.add_option(""  ,"--eosdir"     , dest="eosdir"     , type="string", default="",help="Name of the eos output directory for jobs")
parser.add_option("-c","--config"     , dest="config"     , type="string", default="minintuplizer_cfg.py",help="Ntuplizer config file. Default is: minintuplizer_cfg.py")
parser.add_option("",  "--puppiConfig", dest="puppiConfig", type="string", default="Puppi_cff.py",help="Puppi config file. Default is: Puppi_cff.py")
parser.add_option("-n","--njobs"      , dest="njobs"      , type="int"   , help="Number of jobs")
parser.add_option("-e","--executable" , dest="executable" , type="string", default="MiniNtuplizer",help="Name of the executable. Default is: MiniNtuplizer")
parser.add_option("-q","--queue"      , dest="queue"      , type="string", default="1nh",help="Name of the queue on lxbatch")
parser.add_option(""  ,"--checkJobs"  , dest="checkJobs"  , action="store_true", default=False,help="Checks job status")
parser.add_option(""  ,"--resubmit"   , dest="resubmit"   , action="store_true", default=False,help="Resubmit job ")
parser.add_option("-j","--job"        , dest="job"        , type="int", help="Job number (for resubmission). You can resubmit one job at time for now.")
parser.add_option(""  ,"--dryRun"     , dest="dryRun"     , default=False, action="store_true",help="Do not submit jobs")


(options,args)=parser.parse_args()

eos = '/afs/cern.ch/project/eos/installation/pro/bin/eos.select'


def makeFilesList(indir,wdir):
    list = []
    command = ('%s find -f %s | grep root > %s/list.txt' % (eos,indir,wdir))
    #print command
    os.system(command)
    file = open('%s/list.txt'%wdir, 'r')
    list =[line.replace('/eos/cms/','root://eoscms//eos/cms/').replace('\n','') for line in file ]
    print 'Found %d files' %len(list)
    #print list             
    return list


def writeJobs(wdir, analysis, config, puppiconfig, indir, output, eosoutdir, njobs):
    #---------------------------------------------
    # --- prepare the list of files to be analyzed
    #---------------------------------------------
    listoffiles = []
    listoffiles = makeFilesList(indir,wdir)

    #---------------------------------------------
    # --- now split the jobs
    #---------------------------------------------
    for job in range(njobs):
        print 'job %d' %job
        jobdir = '%s/JOB_%d'%(wdir,job)
        os.mkdir(jobdir)
        
        #--- prepare the list of files for each job
        f = open('%s/input_%d.txt'%(jobdir,job), 'w')
        sublist = [file for i,file in enumerate(listoffiles) if (i%njobs==job)]
        for fname in sublist:
            f.write('%s \n'%fname)

        #--- prepare the jobs scripts
        jobscript = open('%s/sub_%d.sh'%(jobdir,job),'w')
        jobscript.write('cd %s \n'%jobdir)
        #jobscript.write('export SCRAM_ARCH=slc5_amd64_gcc472 \n')
        jobscript.write('export SCRAM_ARCH=slc6_amd64_gcc472 \n')
        jobscript.write('eval ` scramv1 runtime -sh ` \n')
        jobscript.write('cd - \n')
        jobscript.write('cp %s ./ \n'%(puppiconfig)) 
        jobscript.write('if ( \n')
        jobscript.write('\t touch %s/sub_%d.run \n'%(jobdir,job))
        jobscript.write('\t %s %s %s/input_%d.txt %s_%d.root'%(analysis, config, jobdir, job, output, job))
        jobscript.write(') then \n')
        if (eosoutdir == ''):
            jobscript.write('\t cp ./%s_%d.root %s \n'%(output,job,jobdir))
        else:
            jobscript.write('\t cmsStage ./%s_%d.root %s/ \n'%(output,job,eosoutdir))
        jobscript.write('\t touch %s/sub_%d.done \n'%(jobdir,job))
        jobscript.write('else \n')
        jobscript.write('\t touch %s/sub_%d.fail \n'%(jobdir,job))
        jobscript.write('fi \n')
        os.system('chmod a+x %s/sub_%d.sh'%(jobdir,job))
   

def submitJobs(wdir, njobs, queue):
    for job in range(njobs):
        print 'job %d' %job
        jobdir = '%s/JOB_%d'%(wdir,job)
        jobname = '%s/sub_%d.sh'%(jobdir,job)
        print 'bsub -q %s -o %s/sub_%d.log %s'%(queue,jobdir,job,jobname )
        os.system('bsub -q %s -o %s/sub_%d.log %s'%(queue,jobdir,job,jobname ))
        


def checkJobs(wdir, output, queue, eosoutdir):
    jobs = glob.glob( '%s/JOB_*/sub*.sh'% (wdir) )
    print 'Total number of jobs: %d' %len(jobs)
    
    listdone = [j for j in range(len(jobs)) if os.path.isfile('%s/JOB_%d/sub_%d.done' % (wdir,j,j))]
    print 'Total number of DONE jobs: %s ' % len(listdone)
    print '  %s' %listdone
    for j in listdone:
        f = '%s/JOB_%d/sub_%d.run'%(wdir,j,j)
        #print 'rm %s/JOB_%d/sub_%d.run'%(wdir,j,j)
        if (os.path.isfile(f)):
            os.system('rm %s'%f)
            
    listrun = [j for j in range(len(jobs)) if os.path.isfile('%s/JOB_%d/sub_%d.run' % (wdir,j,j))]
    print 'Total number of RUNNING jobs: %d ' %len(listrun)
    print '   %s' %listrun
            
    listfailed = [j for j in range(len(jobs)) if os.path.isfile('%s/JOB_%d/sub_%d.fail' % (wdir,j,j))]
    print 'Failed jobs: %s ' % listfailed
    print '   %s' %listfailed
                        
    if (len(listdone) == len(jobs)):
        print "All jobs successful! Merging output files..."
        os.chdir(wdir)
        #        os.system('ls JOB*/*.root')
        # add automathically only if saved locally
        if (eosoutdir == ''):
            os.system('hadd -f %s.root JOB_*/*.root'%output)
    else:
        for j in listfailed:
            print 'bsub -q %s -o %s/JOB_%d/sub_%d.log %s/JOB_%d/sub_%d.sh'%(queue,wdir,j,j,wdir,j,j )
                    


#-----------------------------------------
#--- MAIN
#-----------------------------------------

path = os.getcwd()
conf = path+'/'+options.config
puppiconf = path+'/'+options.puppiConfig
workingdir = path+'/'+options.workdir

eosoutdir = ''

if not options.checkJobs and not options.resubmit:

    # -- write jobs scripts
    os.mkdir(workingdir)
    if (options.eosdir !=''):
        eosoutdir = options.eosdir+'/'+options.workdir
        mkdir = 'cmsMkdir'
        command = '%s %s'%(mkdir,eosoutdir)
        print command
        os.system(command)
    writeJobs(workingdir, options.executable, conf, puppiconf, options.inputdir, options.outputname, eosoutdir, options.njobs)
    
    # -- submit jobs 
    if not options.dryRun:
        submitJobs(workingdir, options.njobs, options.queue)
        
elif options.resubmit and options.job >-1 :
    print 'Resubmitting job %d ' %options.job
    resubcmd = 'bsub -q %s -o %s/JOB_%d/sub_%d.log %s/JOB_%d/sub_%d.sh'%(options.queue,workingdir,options.job,options.job,workingdir,options.job,options.job )
    #print resubcmd
    os.system(resubcmd)

elif options.checkJobs:
    checkJobs(workingdir,options.outputname, options.queue, eosoutdir)
    
