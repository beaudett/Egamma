import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")
process.load('Configuration.StandardSequences.Services_cff')
process.load("Geometry.CaloEventSetup.CaloTopology_cfi")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Geometry.CaloEventSetup.CaloTowerConstituents_cfi")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')


process.TFileService = cms.Service("TFileService", fileName = cms.string("histo.root") )

# Source : general definition
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(),
                            secondaryFileNames = cms.untracked.vstring(),
                            noEventSort = cms.untracked.bool(True),
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
                            )

# Update input files
#print dbs_discovery.search()
#process.source.fileNames=['/store/relval/CMSSW_3_3_0/RelValSingleElectronPt35/GEN-SIM-RECO/MC_31X_V9-v1/0009/98C1C609-75B7-DE11-941B-001D09F295FB.root','/store/relval/CMSSW_3_3_0/RelValSingleElectronPt35/GEN-SIM-RECO/MC_31X_V9-v1/0008/100559BC-EDB6-DE11-BA68-000423D6006E.root']
#process.source.fileNames=['root://cms-xrd-global.cern.ch///store/relval/CMSSW_9_4_4/RelValTTbar_13/GEN-SIM-RECO/PU25ns_94X_mc2017_realistic_v10For2017H_v2-v1/10000/F6E3F27C-352D-E811-BF29-0242AC130002.root']
process.source.fileNames=['/store/relval/CMSSW_9_4_4/RelValTTbar_13/GEN-SIM-RECO/PU25ns_94X_mc2017_realistic_v10For2017H_v2-v1/10000/F6E3F27C-352D-E811-BF29-0242AC130002.root',
'/store/relval/CMSSW_9_4_4/RelValTTbar_13/GEN-SIM-RECO/PU25ns_94X_mc2017_realistic_v10For2017H_v2-v1/10000/CA997FEE-392D-E811-96DC-0242AC130002.root',
'/store/relval/CMSSW_9_4_4/RelValTTbar_13/GEN-SIM-RECO/PU25ns_94X_mc2017_realistic_v10For2017H_v2-v1/10000/CA6E1F48-3E2D-E811-8258-0242AC130002.root',
'/store/relval/CMSSW_9_4_4/RelValTTbar_13/GEN-SIM-RECO/PU25ns_94X_mc2017_realistic_v10For2017H_v2-v1/10000/82639551-442D-E811-8B22-0242AC130002.root',
'/store/relval/CMSSW_9_4_4/RelValTTbar_13/GEN-SIM-RECO/PU25ns_94X_mc2017_realistic_v10For2017H_v2-v1/10000/8261D7BB-3F2D-E811-9334-0242AC130002.root',
'/store/relval/CMSSW_9_4_4/RelValTTbar_13/GEN-SIM-RECO/PU25ns_94X_mc2017_realistic_v10For2017H_v2-v1/10000/783EA421-352D-E811-B0A6-0242AC130002.root',
'/store/relval/CMSSW_9_4_4/RelValTTbar_13/GEN-SIM-RECO/PU25ns_94X_mc2017_realistic_v10For2017H_v2-v1/10000/380F4C0E-472D-E811-B740-0242AC130002.root',
'/store/relval/CMSSW_9_4_4/RelValTTbar_13/GEN-SIM-RECO/PU25ns_94X_mc2017_realistic_v10For2017H_v2-v1/10000/1AA59F88-422D-E811-88A0-0242AC130002.root',
'/store/relval/CMSSW_9_4_4/RelValTTbar_13/GEN-SIM-RECO/PU25ns_94X_mc2017_realistic_v10For2017H_v2-v1/10000/06E2103A-382D-E811-9A12-0242AC130002.root']

process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(-1)
        )



# Calo geometry service model
#process.load("Configuration.StandardSequences.GeometryRecoDB_cff")

#process.load("Geometry.CaloEventSetup.CaloTopology_cfi")
from FastSimulation.Calorimetry.Calorimetry_cff import FamosCalorimetryBlock
#process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.pfAllElectrons = cms.EDFilter("PdgIdPFCandidateSelector",
    pdgId = cms.vint32(11, -11),
    src = cms.InputTag("pfNoPileUp")
)


process.gensource = cms.EDProducer("GenParticlePruner",
    src = cms.InputTag("genParticles"),
    select = cms.vstring('drop *',
        'keep pdgId = 11',
        'keep pdgId = -11')
)


process.pfPileUp = cms.EDProducer("PFPileUp",
    PFCandidates = cms.InputTag("particleFlow"),
    verbose = cms.untracked.bool(False),
	Enable = cms.bool(True),
    Vertices = cms.InputTag("offlinePrimaryVerticesWithBS")
)


process.pfNoPileUp = cms.EDProducer("TPPileUpPFCandidatesOnPFCandidates",
    bottomCollection = cms.InputTag("particleFlow"),
    topCollection = cms.InputTag("pfPileUp"),
    name = cms.untracked.string('pileUpOnPFCandidates'),
    verbose = cms.untracked.bool(False),
	enable = cms.bool(True)
)
process.pfNoPileUpSequence = cms.Sequence(process.pfPileUp+process.pfNoPileUp)


process.trackExtra = cms.EDAnalyzer("EcalTrackExtra",FamosCalorimetryBlock,
                           GSFTracks=cms.InputTag("electronGsfTracks"),
                           BarrelRecHits=cms.InputTag("reducedEcalRecHitsEB"),
                           EndcapRecHits=cms.InputTag("reducedEcalRecHitsEE"),
                           PFCandidate = cms.InputTag("particleFlow"),
                           InputTruthLabel = cms.InputTag("gensource"),
                           DeltaRMatch = cms.double(0.2),
                           MatchWithMC = cms.bool(True),
						   maxPt = cms.double(10.),
                           TestParticleFilter = cms.PSet(
    # Particles with |eta| > etaMax (momentum direction at primary vertex)
    # are not simulated
    etaMax = cms.double(5.0),
    # Charged particles with pT < pTMin (GeV/c) are not simulated
    pTMin = cms.double(0.0),
    # Particles with energy smaller than EMin (GeV) are not simulated
    EMin = cms.double(0.0),
    # Protons with energy in excess of this value (GeV) will kept no matter what
    EProton = cms.double(99999.0),
    chargedPtMin = cms.double(0.1),
    invisibleParticles = cms.vint32(),
	protonEMin = cms.double(5000.0),
    ))

process.p =cms.Path(
#process.pfNoPileUpSequence+process.pfAllElectrons+
                    process.gensource+
                    process.trackExtra)


#process.out = cms.OutputModule("PoolOutputModule",
#                               outputCommands = cms.untracked.vstring('keep *'),
#                               outputFile = cms.string(os.environ['TEST_OUTPUT_FILE'])
#                               )
#process.outpath = cms.EndPath(process.out)

process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.MessageLogger.cerr.FwkReport.reportEvery = 100
