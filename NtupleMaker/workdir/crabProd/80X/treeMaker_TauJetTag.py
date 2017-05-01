import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras

import os
import sys
import re
import tarfile


## --------------------------
## -- Command line options --
## --------------------------

import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing ('standard')

options.register('era', "Run2_25ns", VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "Run2_25ns or Run2_50ns")
options.register('ntpVersion', "Ntp_80X_12Jul2016_v8.0", VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "ntpVersion: to be same as the tag of the release. But can be used to produce 72X ntuple as well!")
options.register('GlobalTag', "80X_mcRun2_asymptotic_2016_miniAODv2", VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "74X PromptReco: 74X_dataRun2_Prompt_v0")
options.register('cmsswVersion', '80X', VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "'36X' for example. Used for specific MC fix")
options.register('specialFix', 'JEC', VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "special fixes ==>   JEC : use external JEC; IVF : fix IVF; BADMUON : bad muon filters")
options.register('jecDBname', "Spring16_25nsV6_MC", VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "Summer15_25nsV6_DATA for data")
options.register('hltName', 'HLT', VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "HLT menu to use for trigger matching")

options.register('mcInfo', True, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.bool, "process MonteCarlo data, default is data")

options.register('doPDFs', False, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.bool, "switch to enable the production of PDF weights for NNPDF3.0, CT10, MMHT2014, n.b. you need to do `scram setup lhapdf` before running (default=False)")

options.register('addJetsForZinv', True, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.bool, "switch to add top projected jets for Zinv")

options.register('externalFilterList', '', VarParsing.VarParsing.multiplicity.list, VarParsing.VarParsing.varType.string, "event list for filters")

options.register('jetCorrections', 'L2Relative', VarParsing.VarParsing.multiplicity.list, VarParsing.VarParsing.varType.string, "Level of jet corrections to use: Note the factors are read from DB via GlobalTag")
options.jetCorrections.append('L3Absolute')

options.register('mcVersion', '', VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "'36X' for example. Used for specific MC fix")
options.register('dataVersion', '', VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "'36X' for example. Used for specific DATA fix")

options.register('jetTypes', 'AK4PF', VarParsing.VarParsing.multiplicity.list, VarParsing.VarParsing.varType.string, "Additional jet types that will be produced (AK4Calo and AK4PF, cross cleaned in PF2PAT, are included anyway)")
options.register('hltSelection', '*', VarParsing.VarParsing.multiplicity.list, VarParsing.VarParsing.varType.string, "hlTriggers (OR) used to filter events. for data: ''HLT_Mu9', 'HLT_IsoMu9', 'HLT_IsoMu13_v*''; for MC, HLT_Mu9")
options.register('addKeep', '', VarParsing.VarParsing.multiplicity.list, VarParsing.VarParsing.varType.string, "Additional keep and drop statements to trim the event content")

options.register('debug', False, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.bool, "switch on/off debug mode")
options.register('verbose', 0, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "verbose of debug")

options.register('test', False, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.bool, "switch on/off debug mode")

options.register('doPtHatWeighting', False, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.bool, "PtHat weighting for QCD flat samples, default is False")

options.register('fileslist', '', VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "name of a file with input source list")

options.register('fastsim', False, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.bool, "fastsim sample or not, default is False")

options.register('doTopTagger', False, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.bool, "do top tagger or not, default is True")

options.register('usePhiCorrMET', False, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.bool, "use phi corrected MET or not, default is False")

options.register('reducedfilterTags', True, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.bool, "use phi corrected MET or not, default is True")

options.register('smsModel', 'T1tttt', VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "SMS model name")
options.register('smsMotherMass',  -1, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "SMS mother mass")
options.register('smsDaughterMass',  -1, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "SMS daughter mass")
options.register('selSMSpts', False, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.bool, "select model pobools")

options.register('pythia8', True, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.bool, "pythi8 or not, default True")

options.parseArguments()
options._tagOrder =[]

print options

## -------------------------
## -- Check CMSSW version --
## -------------------------

procCMSSWver = os.environ['CMSSW_RELEASE_BASE'].split("/")[-1]
print "procCMSSWver : ", procCMSSWver, "\n"

if not "CMSSW_8_0" in procCMSSWver:
   print "You should be using CMSSW 80X!! Please change your release area"
   sys.exit("ERROR: Not using 80X release")

if not options.cmsswVersion == "80X":
   print "You should be using CMSSW 80X as option!! Please change to be consistent with the release area"
   sys.exit("ERROR: Not using 80X option")


## ------------------------
## -- Define the process --
## ------------------------

process = cms.Process("SUSY")

if options.era == "Run2_25ns":
   process = cms.Process("SUSY", eras.Run2_25ns)
elif options.era == "Run2_50ns":
   process = cms.Process("SUSY", eras.Run2_50ns)

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

## -- MessageLogger --
process.MessageLogger.cerr.FwkReport.reportEvery = 100
if options.debug:
   process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.MessageLogger.suppressWarning = cms.untracked.vstring('ecalLaserCorrFilter','manystripclus53X','toomanystripclus53X')

## -- Options and Output Report --
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.options.allowUnscheduled = cms.untracked.bool(True)

## -- Maximal Number of Events --
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

if options.debug and options.verbose ==1:
   process.SimpleMemoryCheck = cms.Service('SimpleMemoryCheck',
                                        ignoreTotal=cms.untracked.int32(0),
                                        oncePerEventMode = cms.untracked.bool(False)
                                        )
   process.Timing = cms.Service("Timing")

## -- Input Source --
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
      '/store/relval/CMSSW_3_8_0_pre8/RelValTTbar/GEN-SIM-RECO/START38_V6-v1/0004/847D00B0-608E-DF11-A37D-003048678FA0.root'
    )
)

if options.files:
   process.source.fileNames = options.files
elif options.fileslist:
   inputfiles = cms.untracked.vstring()
   if os.path.exists(options.fileslist) == False or os.path.isfile(options.fileslist) == False:
      sys.exit(5)
   else:
      ifile = open(options.fileslist, 'r')
      for line in ifile.readlines():
         inputfiles.append(line)
   print "inputfiles : \n", inputfiles, "\n"
   process.source.fileNames = inputfiles
else:
   process.source.fileNames = [
#        '/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_M-50_HT-1200to2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/00000/00B2B39D-5D4D-E611-8BD4-002590D9D8B6.root',
        '/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-1200to2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/60000/00D97021-CFBE-E611-AD3F-0025901D08B8.root',
#        '/store/mc/RunIISpring16MiniAODv2/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/10000/7CE5EA6A-F132-E611-9E20-008CFA1660A8.root',
#       '/store/mc/RunIISpring16MiniAODv2/TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/50000/041F3A63-431E-E611-9E1E-008CFA1112CC.root',
#       '/store/mc/RunIISpring16MiniAODv2/SMS-T1tttt_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16Fast_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/00000/004A27F0-5132-E611-A936-02163E016171.root',
#       '/store/mc/RunIISpring16MiniAODv2/SMS-T1tttt_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16Fast_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/00000/00775AA9-5132-E611-A4FE-001E675049F5.root',
#       '/store/mc/RunIISpring16MiniAODv2/SMS-T1tttt_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16Fast_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/00000/020ABCF8-5B32-E611-A85E-02163E017932.root',
#       '/store/mc/RunIISpring16MiniAODv2/SMS-T1tttt_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16Fast_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/00000/0246C91A-A131-E611-95E8-02163E00E646.root',

#       '/store/data/Run2016C/HTMHT/MINIAOD/PromptReco-v2/000/275/420/00000/4AD126B0-F539-E611-AD77-02163E013390.root',
#       '/store/data/Run2016C/HTMHT/MINIAOD/PromptReco-v2/000/275/476/00000/5C70C94E-0D3A-E611-AE2E-02163E01346C.root',
#       '/store/data/Run2016C/HTMHT/MINIAOD/PromptReco-v2/000/275/588/00000/BE1D644B-393A-E611-905F-02163E0124F5.root',
   ]

## ---------------------
## -- Calibration tag --
## ---------------------

from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag

if options.GlobalTag:
   process.GlobalTag = GlobalTag(process.GlobalTag, options.GlobalTag, '')

if options.mcInfo == False: 
   options.jetCorrections.append('L2L3Residual')
options.jetCorrections.insert(0, 'L1FastJet')

## ---------------------------------------
## -- Create all needed jet collections --
## ---------------------------------------


## -- Other Analysis related configuration --
import HLTrigger.HLTfilters.hltHighLevel_cfi as hlt
if options.hltSelection:
   process.hltFilter = hlt.hltHighLevel.clone(
      TriggerResultsTag = cms.InputTag("TriggerResults","",options.hltName),
      HLTPaths = cms.vstring(options.hltSelection),
      throw = True, # Don't throw?!
      andOr = True
   )

############################# EDN SUSYPAT specifics ####################################

process.dummyCounter = cms.EDProducer("EventCountProducer")

# Other sequence
process.comb_seq = cms.Sequence(
  # All cleaning && all basic variables, e.g., mht, ht...     
  process.cleanpatseq *
  # hlt requirement
  process.hltFilter *
  process.prodMuons * process.prodElectrons * 
  process.QGTagger * process.QGTaggerOther * process.QGTaggerNoLep * process.QGAK4PFCHS *
  process.weightProducer *
  process.trackIsolation * process.loosetrackIsolation * process.refalltrackIsolation * 
  process.stopPFJets * process.stopBJets *
  process.ra2Objects *
  process.prepareCutvars_seq
)

process.load("SusyAnaTools.SkimsAUX.prodIsoTrks_cfi")

process.stopTreeMaker.vectorTLorentzVector.append(cms.InputTag("prodIsoTrks:trksForIsoVetoLVec"))

process.ak4Stop_Path = cms.Path(
                                   process.comb_seq * 
                                   process.printDecayPythia8 * process.prodGenInfo * 
                                   process.prodMuonsNoIso * process.prodElectronsNoIso * process.prodIsoTrks *  
                                   process.prodJets * process.prodMET * process.prodEventInfo * process.trig_filter_seq * 
                                   process.type3topTagger *
                                   process.stopTreeMaker
)

###-- Dump config ------------------------------------------------------------
if options.debug:
   file = open('allDump_cfg.py','w')
   file.write(str(process.dumpPython()))
   file.close()
