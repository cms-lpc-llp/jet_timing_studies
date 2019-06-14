import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
#------ Setup ------#

#Options
#options = VarParsing ('analysis')
#options.register('isdata',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,"")
#options.register('isfourjet',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,"")
#options.register('isqcd',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,"")
#options.outputFile = 'jetNtuple.root'
#options.inputFiles = 'file:/mnt/hadoop/store/group/phys_llp/RunIISummer17_QCD/RunIISummer17DRPremix_QCD_HT300-500_AODSIM_100.txt'
#options.inputFiles = '/store/mc/RunIISummer17DRPremix/QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/AODSIM/92X_upgrade2017_realistic_v10-v2/00000/BCA5EDA1-50AC-E711-BAB8-0CC47A4C8E8A.root'
#options.maxEvents = -1

#options.parseArguments()


#initialize the process
process = cms.Process("JetTimingStudies")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff")
process.load("Configuration.EventContent.EventContent_cff")

#load input files
process.source = cms.Source("PoolSource",
   # fileNames = cms.untracked.vstring(options.inputFiles),
    fileNames = cms.untracked.vstring(
#'file:/mnt/hadoop/store/mc/RunIIFall17DRPremix/WJetsToLNu_1J_TuneCP5_13TeV-amcatnloFXFX-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11_ext1-v1/110000/34256060-4155-E911-9B57-002590FD5694.root',
#'file:/mnt/hadoop/store/mc/RunIIFall17DRPremix/WJetsToLNu_1J_TuneCP5_13TeV-amcatnloFXFX-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11_ext1-v1/110000/501F0F4A-7455-E911-AD16-002590D9D9DA.root',
#'file:/mnt/hadoop/store/mc/RunIIFall17DRPremix/WJetsToLNu_1J_TuneCP5_13TeV-amcatnloFXFX-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11_ext1-v1/110000/8058CCE1-6055-E911-BC5B-AC1F6BB17832.root',
#'file:/mnt/hadoop/store/mc/RunIIFall17DRPremix/WJetsToLNu_1J_TuneCP5_13TeV-amcatnloFXFX-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11_ext1-v1/110000/80B6D4DC-F257-E911-8DC5-0CC47A57D1F8.root',
#'file:/mnt/hadoop/store/mc/RunIIFall17DRPremix/WJetsToLNu_1J_TuneCP5_13TeV-amcatnloFXFX-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11_ext1-v1/110000/8CEB8CD7-AA56-E911-948F-0CC47AB0B704.root',
#'file:/mnt/hadoop/store/mc/RunIIFall17DRPremix/WJetsToLNu_1J_TuneCP5_13TeV-amcatnloFXFX-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11_ext1-v1/110000/906D20A0-F155-E911-A247-00259048A862.root',
#'file:/mnt/hadoop/store/mc/RunIIFall17DRPremix/WJetsToLNu_1J_TuneCP5_13TeV-amcatnloFXFX-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11_ext1-v1/110000/C47D3025-0857-E911-B792-00259048A8F4.root',
#'file:/mnt/hadoop/store/mc/RunIIFall17DRPremix/WJetsToLNu_1J_TuneCP5_13TeV-amcatnloFXFX-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11_ext1-v1/110000/CA371157-8B56-E911-A413-0CC47A57CCF4.root',
#'file:/mnt/hadoop/store/mc/RunIIFall17DRPremix/WJetsToLNu_1J_TuneCP5_13TeV-amcatnloFXFX-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11_ext1-v1/110000/E841590F-2756-E911-8DD3-0CC47A0AD742.root',
#'file:/mnt/hadoop/store/mc/RunIIFall17DRPremix/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v3/260000/36BB19B7-E373-E911-BF36-0CC47A0AD630.root',
#'file:/mnt/hadoop/store/mc/RunIIFall17DRPremix/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v3/260000/62643DEC-3C74-E911-B041-002590D9D8D4.root',
#'file:/mnt/hadoop/store/mc/RunIIFall17DRPremix/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v3/260000/C8299DC5-3C74-E911-9C58-002590D9D990.root',

#'file:/mnt/hadoop/store/mc/RunIIFall17DRPremix/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v3/110000/0C2D4B43-2B5E-E911-9E54-002590791DA2.root',
#'file:/mnt/hadoop/store/mc/RunIIFall17DRPremix/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v3/110000/28EDCBB2-F55D-E911-8A64-0CC47A0AD74E.root',
#'file:/mnt/hadoop/store/mc/RunIIFall17DRPremix/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v3/110000/C036D6D1-C15E-E911-8C48-AC1F6BB17832.root',

#'file:/mnt/hadoop/store/mc/RunIIFall17DRPremix/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v3/70000/00886039-1259-E911-86DC-002590FD5A48.root',
#'file:/mnt/hadoop/store/mc/RunIIFall17DRPremix/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v3/70000/1813DDA4-2C59-E911-993E-0CC47A0AD3BC.root',
#'file:/mnt/hadoop/store/mc/RunIIFall17DRPremix/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v3/70000/206F691F-1259-E911-955B-0CC47AD24D28.root',
#'file:/mnt/hadoop/store/mc/RunIIFall17DRPremix/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v3/70000/2095DEFF-1C5A-E911-B3B7-0025907D2502.root',
#'file:/mnt/hadoop/store/mc/RunIIFall17DRPremix/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v3/70000/2E5D0ACE-2A59-E911-82B3-0CC47A0AD4A4.root',
#'file:/mnt/hadoop/store/mc/RunIIFall17DRPremix/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v3/70000/3E918620-F85A-E911-B4E6-0CC47A57D036.root',
#'file:/mnt/hadoop/store/mc/RunIIFall17DRPremix/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v3/70000/466A9267-2E5A-E911-BF85-002590D9D8BE.root',
#'file:/mnt/hadoop/store/mc/RunIIFall17DRPremix/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v3/70000/508DDB42-E358-E911-A3EF-0CC47A57CBCC.root',
#'file:/mnt/hadoop/store/mc/RunIIFall17DRPremix/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v3/70000/54111326-945A-E911-8A79-00259029ED0E.root',
#'file:/mnt/hadoop/store/mc/RunIIFall17DRPremix/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v3/70000/5C492124-1259-E911-9124-AC1F6BB17836.root',
#'file:/mnt/hadoop/store/mc/RunIIFall17DRPremix/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v3/70000/AA35FD00-1E5A-E911-B581-0025907D2502.root',
#'file:/mnt/hadoop/store/mc/RunIIFall17DRPremix/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v3/70000/B49B67CD-1C5A-E911-8E50-0CC47A57D086.root',
#'file:/mnt/hadoop/store/mc/RunIIFall17DRPremix/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v3/70000/D0311E5B-2C59-E911-8910-0CC47AA478BE.root',
#'file:/mnt/hadoop/store/mc/RunIIFall17DRPremix/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v3/70000/E2319D5C-F458-E911-9BC7-002590D9D8AE.root',
'file:/mnt/hadoop/store/user/christiw/RunIISummer16_withISR/ppTohToSS1SS2_SS1Tobb_SS2Tobb_vh_withISR_MC_prod/ppTohToSS1SS2_SS1Tobb_SS2Tobb_vh_withISR_mh125_mx50_pl1000_ev100000/crab_CMSSW_8_0_21_ppTohToSS1SS2_SS1Tobb_SS2Tobb_vh_ISR_mh125_mx50_pl1000_ev100000_AOD_CaltechT2_v1/190417_225744/0000/ppTohToSS1SS2_SS1Tobb_SS2Tobb_ggh_withISR_step2_2.root'
#'file:/mnt/hadoop/store/mc/RunIIFall17DRPremix/QCD_Pt_170to300_TuneCP5_13TeV_pythia8/AODSIM/94X_mc2017_realistic_v11_ext1-v1/100000/18BB89F2-8E32-E811-80D1-0025902BD8CE.root',

    ),
)

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100

#TFileService for output
process.TFileService = cms.Service("TFileService",
	fileName = cms.string('jet_timing_studies_ntuple_RunIIFall17DRPremix_WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_1.root'),
    closeFileFast = cms.untracked.bool(True)
)

#load run conditions
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.Geometry.GeometryIdeal_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')

#------ Declare the correct global tag ------#


process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_v3'
#process.GlobalTag.globaltag = '92X_upgrade2017_realistic_v10'

#------ If we add any inputs beyond standard miniAOD event content, import them here ------#

#process.load('CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi')
#process.HBHENoiseFilterResultProducer.minZeros = cms.int32(99999)
#process.HBHENoiseFilterResultProducer.IgnoreTS4TS5ifJetInLowBVRegion=cms.bool(False)
#process.HBHENoiseFilterResultProducer.defaultDecision = cms.string("HBHENoiseFilterResultRun2Loose")

#process.load('RecoMET.METFilters.BadChargedCandidateFilter_cfi')
#process.BadChargedCandidateFilter.muons = cms.InputTag("slimmedMuons")
#process.BadChargedCandidateFilter.PFCandidates = cms.InputTag("packedPFCandidates")
#process.BadChargedCandidateFilter.taggingMode = cms.bool(True)

process.load('RecoMET.METFilters.BadPFMuonFilter_cfi')
#process.BadPFMuonFilter.muons = cms.InputTag("slimmedMuons")
#process.BadPFMuonFilter.PFCandidates = cms.InputTag("packedPFCandidates")
#process.BadPFMuonFilter.taggingMode = cms.bool(True)

#------ Analyzer ------#

#list input collections
process.ntuples = cms.EDAnalyzer('jet_timing_studies',
    isData = cms.bool(False),
    useGen = cms.bool(True),
    isFastsim = cms.bool(False),
    enableTriggerInfo = cms.bool(True),
    enableRecHitInfo = cms.bool(True),
    readGenVertexTime = cms.bool(True),#needs to be false for glueball samples
    isQCD = cms.bool(True),
    model = cms.int32(1),
    genParticles_t0 = cms.InputTag("genParticles", "t0", ""),
    triggerPathNamesFile = cms.string("cms_lpc_llp/jet_timing_studies/data/trigger_names_llp_v1.dat"),
    eleHLTFilterNamesFile = cms.string("SUSYBSMAnalysis/RazorTuplizer/data/RazorElectronHLTFilterNames.dat"),
    muonHLTFilterNamesFile = cms.string("SUSYBSMAnalysis/RazorTuplizer/data/RazorMuonHLTFilterNames.dat"),
    photonHLTFilterNamesFile = cms.string("SUSYBSMAnalysis/RazorTuplizer/data/RazorPhotonHLTFilterNames.dat"),

    vertices = cms.InputTag("offlinePrimaryVerticesWithBS"),  # for non-timing case
    muons = cms.InputTag("muons"),
    electrons = cms.InputTag("gedGsfElectrons"),
    taus = cms.InputTag("hpsPFTauProducer"),
    photons = cms.InputTag("gedPhotons"),
    jetsCalo = cms.InputTag("ak4CaloJets","","RECO"),
    jets = cms.InputTag("ak4PFJetsCHS"),
    jetsPuppi = cms.InputTag("ak4PFJets"),
    jetsAK8 = cms.InputTag("ak8PFJetsCHS"),
    mets = cms.InputTag("pfMet"),
    #metsNoHF = cms.InputTag("pfMet30"),
    metsPuppi = cms.InputTag("pfMet"),
    pfCands = cms.InputTag("particleFlow","","RECO"),

    #packedPfCands = cms.InputTag("packedPFCandidates"),

    genParticles = cms.InputTag("genParticles"),

    #packedGenParticles = cms.InputTag("packedGenParticles"),
    #prunedGenParticles = cms.InputTag("prunedGenParticles"),
    genMetsCalo = cms.InputTag("genMetCalo"),
    genMetsTrue = cms.InputTag("genMetTrue"),
    genJets = cms.InputTag("ak4GenJets"),

    triggerBits = cms.InputTag("TriggerResults","","HLT"),
    #triggerBits = cms.InputTag("TriggerResults","","RECO"),
    hepMC = cms.InputTag("generatorSmeared", "", "SIM"),

    #triggerPrescales = cms.InputTag("patTrigger"),
    #triggerObjects = cms.InputTag("selectedPatTrigger"),

    metFilterBits = cms.InputTag("TriggerResults", "", "RECO"),

    #hbheNoiseFilter = cms.InputTag("HBHENoiseFilterResultProducer","HBHENoiseFilterResult"),
    #hbheTightNoiseFilter = cms.InputTag("HBHENoiseFilterResultProducer","HBHENoiseFilterResultRun2Tight"),
    #hbheIsoNoiseFilter = cms.InputTag("HBHENoiseFilterResultProducer","HBHEIsoNoiseFilterResult"),

    #BadChargedCandidateFilter = cms.InputTag("BadChargedCandidateFilter",""),
    #BadMuonFilter = cms.InputTag("BadPFMuonFilter",""),

    #lheInfo = cms.InputTag("externalLHEProducer", "", ""),
    genInfo = cms.InputTag("generator", "", "SIM"),

    tracks = cms.InputTag("generalTracks", "", "RECO"),
    #trackTime = cms.InputTag("trackTimeValueMapProducer","generalTracksConfigurableFlatResolutionModel"),
    #trackTimeReso = cms.InputTag("trackTimeValueMapProducer","generalTracksConfigurableFlatResolutionModelResolution"),

    puInfo = cms.InputTag("addPileupInfo", "", "HLT"), #uncomment if no pre-mixing
    #puInfo = cms.InputTag("mixData", "", "HLT"), #uncomment for samples with pre-mixed pileup
    #hcalNoiseInfo = cms.InputTag("hcalnoise", "", "RECO"),

    secondaryVertices = cms.InputTag("inclusiveSecondaryVertices", "", "RECO"),

    rhoAll = cms.InputTag("fixedGridRhoAll", "", "RECO"),

    rhoFastjetAll = cms.InputTag("fixedGridRhoFastjetAll", "", "RECO"),
    rhoFastjetAllCalo = cms.InputTag("fixedGridRhoFastjetAllCalo", "", "RECO"),
    rhoFastjetCentralCalo = cms.InputTag("fixedGridRhoFastjetCentralCalo", "", "RECO"),
    rhoFastjetCentralChargedPileUp = cms.InputTag("fixedGridRhoFastjetCentralChargedPileUp", "", "RECO"),
    rhoFastjetCentralNeutral = cms.InputTag("fixedGridRhoFastjetCentralNeutral", "", "RECO"),

    beamSpot = cms.InputTag("offlineBeamSpot", "", "RECO"),
    pfClusters = cms.InputTag("particleFlowClusterECAL","","RECO"),
    ebRecHits = cms.InputTag("reducedEcalRecHitsEB", "","RECO"),
    #ebRecHits = cms.InputTag("EcalRecHit", "reducedEcalRecHitsEB", "RECO"),
    eeRecHits  = cms.InputTag("reducedEcalRecHitsEE", "","RECO"),
    esRecHits = cms.InputTag("reducedEcalRecHitsES", "","RECO"),
    #ebeeClusters = cms.InputTag("reducedEgamma", "reducedEBEEClusters", "RECO"),
    ebeeClusters = cms.InputTag("particleFlowEGamma", "EBEEClusters", "RECO"),
    esClusters = cms.InputTag("particleFlowEGamma", "ESClusters", "RECO"),
    #conversions = cms.InputTag("reducedEgamma", "reducedConversions", "RECO"),
    conversions = cms.InputTag("allConversions", "", "RECO"),

    #singleLegConversions = cms.InputTag("reducedEgamma", "reducedSingleLegConversions", "RECO"),
    singleLegConversions = cms.InputTag("particleFlowEGamma", "", "RECO"),

    gedGsfElectronCores = cms.InputTag("gedGsfElectronCores", "", "RECO"),
    gedPhotonCores = cms.InputTag("gedPhotonCore", "", "RECO"),
    #superClusters = cms.InputTag("reducedEgamma", "reducedSuperClusters", "RECO"),

    #lostTracks = cms.InputTag("lostTracks", "", "RECO")
)

#run
process.p = cms.Path( #process.HBHENoiseFilterResultProducer*
                      #process.BadChargedCandidateFilter*
                      #process.BadPFMuonFilter*
                      process.ntuples)
