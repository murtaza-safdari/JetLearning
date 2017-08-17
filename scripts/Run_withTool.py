import ROOT
import logging
import shutil
import os
import pdb

logging.basicConfig(level=logging.INFO)
from optparse import OptionParser

parser = OptionParser()

# job configuration
parser.add_option("--submitDir", help="dir to store the output", default="submit_dir_z1and2")
parser.add_option("--inputDQ2", dest='use_scanDQ2', action='store_true', help="input files from DQ2", default=False )
parser.add_option("--inputFileList", dest='use_inputFileList', action='store_true', help="input files from list", default=False )
parser.add_option("--inputFiles", type=str, dest='search_directories', help="list of input files", default = "/atlas/local/acukierm/dijetz1and2/")
parser.add_option("--driver", help="select where to run", choices=("direct", "prooflite", "grid", "lsf"), default="direct")
parser.add_option("--files_per_job", help="how many files to submit per job", type=int, default=1)
parser.add_option("--nevents", type=int, help="number of events to process for all the datasets")
parser.add_option("--skip-events", type=int, help="skip the first n events")
parser.add_option("-w", "--overwrite", action='store_true', default=False, help="overwrite previous submitDir")
# input configuration
parser.add_option("--dataSource", help="undefined=-1, data=0, FullSim=1, AF-II=2 ", type="int", default=1)
parser.add_option("--doLC", dest='doLC', action='store_true', help="Use LC Weights", default=False )
parser.add_option("--doClusterInfo", dest='doClusterInfo', action='store_true', help="Store all cluster info for each jet", default=False )
parser.add_option("--doJetReclustering", dest='doJetReclustering', action='store_true', help="Make custom jets", default=False )
parser.add_option("--doTracks", dest='doTracks', action='store_true', help="Store track info", default=False )
parser.add_option("--doAllClusters", dest='doAllClusters', action='store_true', help="Store all cluster info from event", default=False )

(options, args) = parser.parse_args()

# check submission directory
if os.path.exists(options.submitDir):
  if options.overwrite:
    logging.info("removing {0:s}".format(options.submitDir))
    import shutil
    shutil.rmtree(options.submitDir, True)
  else:
    raise OSError('Output directory {0:s} already exists. Either re-run with -w (overwrite), choose a different --submitDir, or rm -rf it yourself.'.format(options.submitDir))

import atexit
@atexit.register
def quite_exit():
    ROOT.gSystem.Exit(0)

logging.info("loading packages")
ROOT.gROOT.Macro("$ROOTCOREDIR/scripts/load_packages.C")

if options.overwrite:
    import shutil
    shutil.rmtree(options.submitDir, True)

#Set up the job for xAOD access:
ROOT.xAOD.Init().ignore();

use_scanDQ2 = options.use_scanDQ2 or options.driver in ['grid']
if use_scanDQ2:
  if os.getenv('XRDSYS') is None:
    raise EnvironmentError('xrootd client is not setup. Run localSetupFAX or equivalent.')

if options.driver == 'lsf':
  if getattr(ROOT.EL, 'LSFDriver') is None:
    raise KeyError('Cannot load the LSF driver from EventLoop. Did you not compile it?')

# create a new sample handler to describe the data files we use
logging.info("creating new sample handler")
sh_all = ROOT.SH.SampleHandler()

# scan for datasets in the given directories
directory = options.search_directories
if options.use_inputFileList:
  ROOT.SH.readFileList(sh_all, "sample", directory)
else:
  if use_scanDQ2:
    #ROOT.SH.scanDQ2(sh_all, directory)
    ROOT.SH.scanRucio(sh_all, directory)
  else:
    ROOT.SH.scanDir(sh_all, directory)

# print out the samples we found
logging.info("%d different datasets found scanning all directories", len(sh_all))

# set the name of the tree in our files
sh_all.setMetaString("nc_tree", "CollectionTree")

# this is the basic description of our job
logging.info("creating new job")
job = ROOT.EL.Job()
job.sampleHandler(sh_all)

if options.driver == 'lsf':
  job.options().setBool(ROOT.EL.Job.optResetShell, False);
  job.options().setDouble(ROOT.EL.Job.optFilesPerWorker, options.files_per_job)

#Set the xAOD access mode of the job:
#job.options().setString( ROOT.EL.Job.optXaodAccessMode,ROOT.EL.Job.optXaodAccessMode_branch );
job.options().setString( ROOT.EL.Job.optXaodAccessMode,ROOT.EL.Job.optXaodAccessMode_athena );
        
if options.nevents:
    logging.info("processing only %d events", options.nevents)
    job.options().setDouble(ROOT.EL.Job.optMaxEvents, options.nevents)

if options.skip_events:
    logging.info("skipping first %d events", options.skip_events)
    job.options().setDouble(ROOT.EL.Job.optSkipEvents, options.skip_events)


# add our algorithm to the job
logging.info("creating algorithms")
analysis = ROOT.JetLearningAnalysis()

setattr(analysis, 'm_doLC', options.doLC)
setattr(analysis, 'm_doClusterInfo', options.doClusterInfo)
setattr(analysis, 'm_doAllClusters', options.doAllClusters)
setattr(analysis, 'm_doJetReclustering', options.doJetReclustering)
if not options.doJetReclustering: setattr(analysis,'m_jets','AntiKt4EMTopoJets')
setattr(analysis,'m_doTracks',options.doTracks)
if options.doJetReclustering: setattr(analysis,'m_doTracks',False) #no ghost tracks

logging.info("adding algorithms")
job.algsAdd(analysis)

# make the driver we want to use:
# this one works by running the algorithm directly
logging.info("creating driver")
driver = None
if (options.driver == "direct"):
  logging.info("running on direct")
  driver = ROOT.EL.DirectDriver()
  logging.info("submit job")
  driver.submit(job, options.submitDir)
elif (options.driver == "prooflite"):
  logging.info("running on prooflite")
  driver = ROOT.EL.ProofDriver()
  logging.info("submit job")
  driver.submit(job, options.submitDir)
elif (options.driver == "grid"):
  logging.info("running on Grid")
  driver = ROOT.EL.PrunDriver()   
  driver.options().setString("nc_outputSampleName", "user.%acukierm%.%in:name[2]%.%in:name[3]%.%in:name[6]%_11") 
  driver.options().setDouble("nc_nGBPerJob", 2)
  #driver.options().setDouble("nc_disableAutoRetry", 1)
  #driver.options().setString("nc_cmtConfig", "x86_64-slc6-gcc48-opt")
  #driver.options().setDouble("nc_mergeOutput", 1);  
  #driver.options().setString("nc_excludedSite", "TRIUMF-LCG2_DATADISK")
  #driver.options().setString("nc_EventLoop_SubmitFlags", "--allowTaskDuplication");
  driver.options().setString("EL::Job::optGridDestSE","CERN-PROD_LOCALGROUPDISK")
  driver.options().setString("nc_optGridDestSE","CERN-PROD_LOCALGROUPDISK")
  driver.options().setDouble(ROOT.EL.Job.optGridMergeOutput, 1)
  logging.info("submit job")
  driver.submitOnly(job, options.submitDir)
elif (options.driver == "lsf"):
  logging.info("running on LSF")
  driver = ROOT.EL.LSFDriver()
  #driver.options().setString(ROOT.EL.Job.optSubmitFlags, '-q medium -W200')
  driver.options().setString(ROOT.EL.Job.optSubmitFlags, '-q atlas-t3 -W300')
  driver.submitOnly(job, options.submitDir)
