// -*- C++ -*-
// Class:      jet_timing_studies
/*
  Description: Base class for AOD analysis with CRAB
*/
//         Author:  Cristian Pena
//         Created:  Tue, 26 March 2019 15:00:06 GMT

#include "jet_timing_studies.h"
//------ Constructors and destructor ------//
jet_timing_studies::jet_timing_studies(const edm::ParameterSet& iConfig):
  //get inputs from config file
  isData_(iConfig.getParameter<bool> ("isData")),
  model_(iConfig.getParameter<int> ("model")),
  useGen_(iConfig.getParameter<bool> ("useGen")),
  isFastsim_(iConfig.getParameter<bool> ("isFastsim")),
  isQCD_(iConfig.getParameter<bool> ("isQCD")),
  enableTriggerInfo_(iConfig.getParameter<bool> ("enableTriggerInfo")),
  enableRecHitInfo_(iConfig.getParameter<bool> ("enableRecHitInfo")),
  readGenVertexTime_(iConfig.getParameter<bool> ("readGenVertexTime")),
  triggerPathNamesFile_(iConfig.getParameter<string> ("triggerPathNamesFile")),
  eleHLTFilterNamesFile_(iConfig.getParameter<string> ("eleHLTFilterNamesFile")),
  muonHLTFilterNamesFile_(iConfig.getParameter<string> ("muonHLTFilterNamesFile")),
  photonHLTFilterNamesFile_(iConfig.getParameter<string> ("photonHLTFilterNamesFile")),
  verticesToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
  tracksTag_(consumes<edm::View<reco::Track> >(iConfig.getParameter<edm::InputTag>("tracks"))),
  cscSegmentInputToken_(consumes<CSCSegmentCollection>(edm::InputTag("cscSegments"))),
  dtSegmentInputToken_(consumes<DTRecSegment4DCollection>(edm::InputTag("dt4DCosmicSegments"))),
  rpcRecHitInputToken_(consumes<RPCRecHitCollection>(edm::InputTag("rpcRecHits"))),
  muonsToken_(consumes<reco::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
  electronsToken_(consumes<reco::GsfElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
  tausToken_(consumes<reco::PFTauCollection>(iConfig.getParameter<edm::InputTag>("taus"))),
  photonsToken_(consumes<reco::PhotonCollection>(iConfig.getParameter<edm::InputTag>("photons"))),
  jetsCaloToken_(consumes<reco::CaloJetCollection>(iConfig.getParameter<edm::InputTag>("jetsCalo"))),
  jetsToken_(consumes<reco::PFJetCollection>(iConfig.getParameter<edm::InputTag>("jets"))),
  jetsPuppiToken_(consumes<reco::PFJetCollection>(iConfig.getParameter<edm::InputTag>("jetsPuppi"))),
  jetsAK8Token_(consumes<reco::PFJetCollection>(iConfig.getParameter<edm::InputTag>("jetsAK8"))),
  PFCandsToken_(consumes<reco::PFCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCands"))),
  PFClustersToken_(consumes<reco::PFClusterCollection>(iConfig.getParameter<edm::InputTag>("pfClusters"))),
  //genParticlesToken_(consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genParticles"))),
  //genParticlesToken_(consumes<edm::View<pat::PackedGenParticle> >(iConfig.getParameter<edm::InputTag>("genParticles"))),
  genParticlesToken_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticles"))),
  genJetsToken_(consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("genJets"))),
  triggerBitsToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerBits"))),
  hepMCToken_(consumes<edm::HepMCProduct>(iConfig.getParameter<edm::InputTag>("hepMC"))),
  //triggerObjectsToken_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("triggerObjects"))),
  //triggerPrescalesToken_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("triggerPrescales"))),
  genMetCaloToken_(consumes<reco::GenMETCollection>(iConfig.getParameter<edm::InputTag>("genMetsCalo"))),
  genMetTrueToken_(consumes<reco::GenMETCollection>(iConfig.getParameter<edm::InputTag>("genMetsTrue"))),
  metToken_(consumes<reco::PFMETCollection>(iConfig.getParameter<edm::InputTag>("mets"))),
//  metNoHFToken_(consumes<reco::PFMETCollection>(iConfig.getParameter<edm::InputTag>("metsNoHF"))),
  metPuppiToken_(consumes<reco::PFMETCollection>(iConfig.getParameter<edm::InputTag>("metsPuppi"))),
  metFilterBitsToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("metFilterBits"))),
  //hbheNoiseFilterToken_(consumes<bool>(iConfig.getParameter<edm::InputTag>("hbheNoiseFilter"))),
  //hbheTightNoiseFilterToken_(consumes<bool>(iConfig.getParameter<edm::InputTag>("hbheTightNoiseFilter"))),
  //hbheIsoNoiseFilterToken_(consumes<bool>(iConfig.getParameter<edm::InputTag>("hbheIsoNoiseFilter"))),
  //badChargedCandidateFilterToken_(consumes<bool>(iConfig.getParameter<edm::InputTag>("BadChargedCandidateFilter"))),
  //badMuonFilterToken_(consumes<bool>(iConfig.getParameter<edm::InputTag>("BadMuonFilter"))),
//  lheRunInfoTag_(iConfig.getParameter<edm::InputTag>("lheInfo")),
//  lheRunInfoToken_(consumes<LHERunInfoProduct,edm::InRun>(lheRunInfoTag_)),
//  lheInfoToken_(consumes<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("lheInfo"))),
  genInfoToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genInfo"))),
  genLumiHeaderToken_(consumes<GenLumiInfoHeader,edm::InLumi>(edm::InputTag("generator",""))),
  puInfoToken_(consumes<std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("puInfo"))),
  //hcalNoiseInfoToken_(consumes<HcalNoiseSummary>(iConfig.getParameter<edm::InputTag>("hcalNoiseInfo"))),
  secondaryVerticesToken_(consumes<vector<reco::VertexCompositePtrCandidate> >(iConfig.getParameter<edm::InputTag>("secondaryVertices"))),
  rhoAllToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoAll"))),
  rhoFastjetAllToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoFastjetAll"))),
  rhoFastjetAllCaloToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoFastjetAllCalo"))),
  rhoFastjetCentralCaloToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoFastjetCentralCalo"))),
  rhoFastjetCentralChargedPileUpToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoFastjetCentralChargedPileUp"))),
  rhoFastjetCentralNeutralToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoFastjetCentralNeutral"))),
  beamSpotToken_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"))),
  ebRecHitsToken_(consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > >(iConfig.getParameter<edm::InputTag>("ebRecHits"))),
  eeRecHitsToken_(consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > >(iConfig.getParameter<edm::InputTag>("eeRecHits"))),
  esRecHitsToken_(consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > >(iConfig.getParameter<edm::InputTag>("esRecHits"))),
  ebeeClustersToken_(consumes<vector<reco::CaloCluster> >(iConfig.getParameter<edm::InputTag>("ebeeClusters"))),
  esClustersToken_(consumes<vector<reco::CaloCluster> >(iConfig.getParameter<edm::InputTag>("esClusters"))),
  conversionsToken_(consumes<vector<reco::Conversion> >(iConfig.getParameter<edm::InputTag>("conversions"))),
  singleLegConversionsToken_(consumes<vector<reco::Conversion> >(iConfig.getParameter<edm::InputTag>("singleLegConversions"))),
  gedGsfElectronCoresToken_(consumes<vector<reco::GsfElectronCore> >(iConfig.getParameter<edm::InputTag>("gedGsfElectronCores"))),
  gedPhotonCoresToken_(consumes<vector<reco::PhotonCore> >(iConfig.getParameter<edm::InputTag>("gedPhotonCores"))),
  generalTrackToken_(consumes<std::vector<reco::Track>>(edm::InputTag("generalTracks")))
  //superClustersToken_(consumes<vector<reco::SuperCluster> >(iConfig.getParameter<edm::InputTag>("superClusters"))),
  //  lostTracksToken_(consumes<vector<reco::PFCandidate> >(iConfig.getParameter<edm::InputTag>("lostTracks")))
{
  //declare the TFileService for output
  edm::Service<TFileService> fs;

  //set up output tree
  llpTree = fs->make<TTree>("llp", "selected AOD information for llp analyses");
  //llpTree = new TTree("Jets", "selected AOD information");
  NEvents = fs->make<TH1F>("NEvents",";;NEvents;",1,-0.5,0.5);
  //*****************************************************************************************
  //Read in HLT Trigger Path List from config file
  //*****************************************************************************************
  for (int i = 0; i<NTriggersMAX; ++i) triggerPathNames[i] = "";
  ifstream myfile (edm::FileInPath(triggerPathNamesFile_.c_str()).fullPath().c_str()) ;
  if (myfile.is_open())
  {
    std::string line;
    int index;
    std::string hltpathname;

    while(myfile>>index>>hltpathname)
    {
      if (index < NTriggersMAX)
      {
        triggerPathNames[index] = hltpathname;
      }
    }
    myfile.close();
  }
  else
  {
    std::cout << "ERROR!!! Could not open trigger path name file : " << edm::FileInPath(triggerPathNamesFile_.c_str()).fullPath().c_str() << "\n";
  }

  if(enableTriggerInfo_)
  {
    std::cout << "\n";
    std::cout << "****************** Trigger Paths Defined For Razor Ntuple ******************\n";
    for (int i = 0; i<NTriggersMAX; ++i)
    {
      if (triggerPathNames[i] != "") std::cout << "Trigger " << i << " " << triggerPathNames[i] << "\n";
    }
    std::cout << "****************************************************************************\n";
    std::cout << "\n";
  }
  if(readGenVertexTime_) genParticles_t0_Token_ = consumes<float>(iConfig.getParameter<edm::InputTag>("genParticles_t0"));
  /*
  fJetPhotonRecHitEta = new std::vector<float>; fJetPhotonRecHitEta->clear();
  fJetPhotonRecHitPhi = new std::vector<float>; fJetPhotonRecHitPhi->clear();
  fJetPhotonRecHitE = new std::vector<float>; fJetPhotonRecHitE->clear();
  fJetPhotonRecHitTime = new std::vector<float>; fJetPhotonRecHitTime->clear();
*/
}

jet_timing_studies::~jet_timing_studies()
{
};

//------ Enable the desired set of branches ------//
void jet_timing_studies::setBranches(){

  llpTree->Branch("isData", &isData, "isData/O");
  llpTree->Branch("model", &model, "model/I");
  llpTree->Branch("isQCD", &isQCD, "isQCD/O");
  llpTree->Branch("runNum", &runNum, "runNum/i");
  llpTree->Branch("lumiNum", &lumiNum, "lumiNum/i");
  llpTree->Branch("eventNum", &eventNum, "eventNum/i");
  llpTree->Branch("pvX", &pvX, "pvX/F");
  llpTree->Branch("pvY", &pvY, "pvY/F");
  llpTree->Branch("pvZ", &pvZ, "pvZ/F");
  llpTree->Branch("nPV", &nPV, "nPV/I");
  llpTree->Branch("Rho", &Rho, "Rho/F");
  llpTree->Branch("nPU", &nPU, "nPU/I");
  llpTree->Branch("nPUmean", &nPUmean, "nPUmean/F");
  llpTree->Branch("PV_x", PV_x, "PV_x[nPV]/F");
  llpTree->Branch("PV_y", PV_y, "PV_y[nPV]/F");
  llpTree->Branch("PV_z", PV_z, "PV_z[nPV]/F");

  llpTree->Branch("nJets", &nJets,"nJets/I");
  llpTree->Branch("jetE", jetE,"jetE[nJets]/F");
  llpTree->Branch("jetEt", jetEt,"jetEt[nJets]/F");
  llpTree->Branch("jetPt", jetPt,"jetPt[nJets]/F");
  llpTree->Branch("jetEta", jetEta,"jetEta[nJets]/F");
  llpTree->Branch("jetPhi", jetPhi,"jetPhi[nJets]/F");
  llpTree->Branch("jetCISV", jetCISV,"jetCISV[nJets]/F");
  llpTree->Branch("jetMass", jetMass, "jetMass[nJets]/F");
  llpTree->Branch("jetAlphaMax",jetAlphaMax,"jetAlphaMax[nJets]/F");
  llpTree->Branch("jetBetaMax",jetBetaMax,"jetBetaMax[nJets]/F");

  llpTree->Branch("jetGammaMax_EM",jetGammaMax_EM,"jetGammaMax_EM[nJets]/F");
  llpTree->Branch("jetGammaMax_Hadronic",jetGammaMax_Hadronic,"jetGammaMax_Hadronic[nJets]/F");
  llpTree->Branch("jetGammaMax",jetGammaMax,"jetGammaMax[nJets]/F");
  llpTree->Branch("jetGammaMax_ET",jetGammaMax_ET,"jetGammaMax_ET[nJets]/F");
  llpTree->Branch("jetGammaMax_P",jetGammaMax_P,"jetGammaMax_P[nJets]/F");

  llpTree->Branch("jetPtAllTracks",jetPtAllTracks,"jetPtAllTracks[nJets]/F");
  llpTree->Branch("jetPtAllPVTracks",jetPtAllPVTracks,"jetPtAllPVTracks[nJets]/F");
  llpTree->Branch("jetMedianTheta2D",jetMedianTheta2D,"jetMedianTheta2D[nJets]/F");
  llpTree->Branch("jetMedianIP",jetMedianIP,"jetMedianIP[nJets]/F");
  llpTree->Branch("jetMinDeltaRAllTracks",jetMinDeltaRAllTracks,"jetMinDeltaRAllTracks[nJets]/F");
  llpTree->Branch("jetMinDeltaRPVTracks",jetMinDeltaRPVTracks,"jetMinDeltaRPVTracks[nJets]/F");
  llpTree->Branch("jetptPVTracks","std::vector<vector<double>>(nJets)",&jetptPVTracks);
  llpTree->Branch("jetnMatchedPVTracks","std::vector<vector<int> >(nJets)",&jetnMatchedPVTracks);
  llpTree->Branch("jetJetArea", jetJetArea, "jetJetArea[nJets]/F");
  llpTree->Branch("jetPileupE", jetPileupE, "jetPileupE[nJets]/F");
  llpTree->Branch("jetPileupId", jetPileupId, "jetPileupId[nJets]/F");
  llpTree->Branch("jetPileupIdFlag", jetPileupIdFlag, "jetPileupIdFlag[nJets]/I");
  llpTree->Branch("jetPassIDLoose", jetPassIDLoose, "jetPassIDLoose[nJets]/O");
  llpTree->Branch("jetPassIDTight", jetPassIDTight, "jetPassIDTight[nJets]/O");
  llpTree->Branch("jetPassMuFrac", jetPassMuFrac, "jetPassMuFrac[nJets]/O");
  llpTree->Branch("jetPassEleFrac", jetPassEleFrac, "jetPassEleFrac[nJets]/O");
  llpTree->Branch("jetPartonFlavor", jetPartonFlavor, "jetPartonFlavor[nJets]/I");
  llpTree->Branch("jetHadronFlavor", jetHadronFlavor, "jetHadronFlavor[nJets]/I");
  llpTree->Branch("jetChargedEMEnergyFraction", jetChargedEMEnergyFraction, "jetChargedEMEnergyFraction[nJets]/F");
  llpTree->Branch("jetNeutralEMEnergyFraction", jetNeutralEMEnergyFraction, "jetNeutralEMEnergyFraction[nJets]/F");
  llpTree->Branch("jetChargedHadronEnergyFraction", jetChargedHadronEnergyFraction, "jetChargedHadronEnergyFraction[nJets]/F");
  llpTree->Branch("jetNeutralHadronEnergyFraction", jetNeutralHadronEnergyFraction, "jetNeutralHadronEnergyFraction[nJets]/F");
  llpTree->Branch("jet_charged_hadron_multiplicity", jet_charged_hadron_multiplicity, "jet_charged_hadron_multiplicity[nJets]/F");
  llpTree->Branch("jet_neutral_hadron_multiplicity", jet_neutral_hadron_multiplicity, "jet_neutral_hadron_multiplicity[nJets]/F");
  llpTree->Branch("jet_photon_multiplicity", jet_photon_multiplicity, "jet_photon_multiplicity[nJets]/F");
  llpTree->Branch("jet_electron_multiplicity", jet_electron_multiplicity, "jet_electron_multiplicity[nJets]/F");
  llpTree->Branch("jet_muon_multiplicity", jet_muon_multiplicity, "jet_muon_multiplicity[nJets]/F");
  llpTree->Branch("jet_HF_hadron_multiplicity", jet_HF_hadron_multiplicity, "jet_HF_hadron_multiplicity[nJets]/F");
  llpTree->Branch("jet_HF_em_multiplicity", jet_HF_em_multiplicity, "jet_HF_em_multiplicity[nJets]/F");
  llpTree->Branch("jet_charged_multiplicity", jet_charged_multiplicity, "jet_charged_multiplicity[nJets]/F");
  llpTree->Branch("jet_neutral_multiplicity", jet_neutral_multiplicity, "jet_neutral_multiplicity[nJets]/F");
  llpTree->Branch("jetMatchedGenPt", jetMatchedGenPt,"jetMatchedGenPt[nJets]/F");
  llpTree->Branch("jetMatchedGenEta", jetMatchedGenEta,"jetMatchedGenEta[nJets]/F");
  llpTree->Branch("jetMatchedGenPhi", jetMatchedGenPhi,"jetMatchedGenPhi[nJets]/F");
  llpTree->Branch("jetMatchedGenMass", jetMatchedGenMass, "jetMatchedGenMass[nJets]/F");
  llpTree->Branch("pfMetEta",&pfMetEta,"pfMetEta/F");
  llpTree->Branch("pfMetE",&pfMetE,"pfMetE/F");
  llpTree->Branch("pfMetPt",&pfMetPt,"pfMetPt/F");
  llpTree->Branch("pfMetPhi",&pfMetPhi,"pfMetPhi/F");



  if( enableRecHitInfo_ )
  {
    llpTree->Branch("event_n_rechits", &event_n_rechits, "event_n_rechits/I");

    llpTree->Branch("jet_n_rechits", jet_n_rechits, "jet_n_rechits[nJets]/I");
    llpTree->Branch("good_jet", good_jet, "good_jet[nJets]/O");
    llpTree->Branch("good_jet0p5", good_jet0p5, "good_jet0p5[nJets]/O");

    llpTree->Branch("jet_n_rechits_Ecut1", jet_n_rechits_Ecut1, "jet_n_rechits_Ecut1[nJets]/I");
    llpTree->Branch("jet_n_rechits_Ecut0p5", jet_n_rechits_Ecut0p5, "jet_n_rechits_Ecut0p5[nJets]/I");

    llpTree->Branch("jet_rechit_T_rms_Ecut1", jet_rechit_T_rms_Ecut1, "jet_rechit_T_rms_Ecut1[nJets]/F");
    llpTree->Branch("jet_rechit_T_rms_Ecut0p5", jet_rechit_T_rms_Ecut0p5, "jet_rechit_T_rms_Ecut0p5[nJets]/F");

    llpTree->Branch("jet_rechit_E_Ecut3", jet_rechit_E_Ecut3, "jet_rechit_E_Ecut3[nJets]/F");
    llpTree->Branch("jet_rechit_T_Ecut3", jet_rechit_T_Ecut3, "jet_rechit_T_Ecut3[nJets]/F");
    llpTree->Branch("jet_rechit_E_Ecut4", jet_rechit_E_Ecut4, "jet_rechit_E_Ecut4[nJets]/F");
    llpTree->Branch("jet_rechit_T_Ecut4", jet_rechit_T_Ecut4, "jet_rechit_T_Ecut4[nJets]/F");
    llpTree->Branch("jet_rechit_E_Ecut2", jet_rechit_E_Ecut2, "jet_rechit_E_Ecut2[nJets]/F");
    llpTree->Branch("jet_rechit_T_Ecut2", jet_rechit_T_Ecut2, "jet_rechit_T_Ecut2[nJets]/F");
    llpTree->Branch("jet_rechit_E_Ecut1p5", jet_rechit_E_Ecut1p5, "jet_rechit_E_Ecut1p5[nJets]/F");
    llpTree->Branch("jet_rechit_T_Ecut1p5", jet_rechit_T_Ecut1p5, "jet_rechit_T_Ecut1p5[nJets]/F");
    llpTree->Branch("jet_rechit_E_Ecut1", jet_rechit_E_Ecut1, "jet_rechit_E_Ecut1[nJets]/F");
    llpTree->Branch("jet_rechit_T_Ecut1", jet_rechit_T_Ecut1, "jet_rechit_T_Ecut1[nJets]/F");
    llpTree->Branch("jet_rechit_E_Ecut0p5", jet_rechit_E_Ecut0p5, "jet_rechit_E_Ecut0p5[nJets]/F");
    llpTree->Branch("jet_rechit_T_Ecut0p5", jet_rechit_T_Ecut0p5, "jet_rechit_T_Ecut0p5[nJets]/F");
    llpTree->Branch("jet_rechit_E", jet_rechit_E, "jet_rechit_E[nJets]/F");
    llpTree->Branch("jet_rechit_T", jet_rechit_T, "jet_rechit_T[nJets]/F");


    llpTree->Branch("jet_pv_rechit_T_Ecut3", jet_pv_rechit_T_Ecut3, "jet_rechit_T_Ecut3[nJets]/F");
    llpTree->Branch("jet_pv_rechit_T_Ecut4", jet_pv_rechit_T_Ecut4, "jet_rechit_T_Ecut4[nJets]/F");
    llpTree->Branch("jet_pv_rechit_T_Ecut2", jet_pv_rechit_T_Ecut2, "jet_rechit_T_Ecut2[nJets]/F");
    llpTree->Branch("jet_pv_rechit_T_Ecut1p5", jet_pv_rechit_T_Ecut1p5, "jet_rechit_T_Ecut1p5[nJets]/F");
    llpTree->Branch("jet_pv_rechit_T_Ecut1", jet_pv_rechit_T_Ecut1, "jet_rechit_T_Ecut1[nJets]/F");
    llpTree->Branch("jet_pv_rechit_T_Ecut0p5", jet_pv_rechit_T_Ecut0p5, "jet_pv_rechit_T_Ecut0p5[nJets]/F");
    llpTree->Branch("jet_pv_rechit_T", jet_pv_rechit_T, "jet_rechit_T[nJets]/F");
    // llpTree->Branch("jet_rechits_phi", jet_rechits_phi, "jet_rechits_phi[event_n_rechits]/F");
    // llpTree->Branch("jet_rechits_eta", jet_rechits_eta, "jet_rechits_eta[event_n_rechits]/F");
    // llpTree->Branch("jet_rechits_E", jet_rechits_E, "jet_rechits_E[event_n_rechits]/F");
    // llpTree->Branch("jet_rechits_T", jet_rechits_T, "jet_rechits_T[event_n_rechits]/F");
    llpTree->Branch("jet_rechits_phi","std::vector<vector<double>>(nJets)",&jet_rechits_phi);
    llpTree->Branch("jet_rechits_eta","std::vector<vector<double>>(nJets)",&jet_rechits_eta);
    llpTree->Branch("jet_rechits_E","std::vector<vector<double>>(nJets)",&jet_rechits_E);
    llpTree->Branch("jet_rechits_T","std::vector<vector<double>>(nJets)",&jet_rechits_T);

    // llpTree->Branch("jet_rechits_E", jet_rechits_E, "jet_rechits_E[nJets][500]/F");
    // llpTree->Branch("jet_rechits_T", jet_rechits_T, "jet_rechits_T[nJets][500]/F");
    llpTree->Branch("jet_pv_rechits_T", jet_pv_rechits_T, "jet_pv_rechits_T[event_n_rechits]/F");
  }

  llpTree->Branch("nPhotons", &fJetNPhotons,"nPhotons/I");
  llpTree->Branch("phoPt", fJetPhotonPt,"phoPt[nPhotons]/F");
  llpTree->Branch("phoEta", fJetPhotonEta,"phoEta[nPhotons]/F");
  llpTree->Branch("phoPhi", fJetPhotonPhi,"phoPhi[nPhotons]/F");
  llpTree->Branch("phoSeedRecHitEta", fJetPhotonSeedRecHitEta, "phoSeedRecHitEta[nPhotons]/F");
  llpTree->Branch("phoSeedRecHitPhi", fJetPhotonSeedRecHitPhi, "phoSeedRecHitPhi[nPhotons]/F");
  llpTree->Branch("phoSeedRecHitE", fJetPhotonSeedRecHitE, "phoSeedRecHitE[nPhotons]/F");
  llpTree->Branch("phoSeedRecHitT", fJetPhotonSeedRecHitTime, "phoSeedRecHitT[nPhotons]/F");

  // llpTree->Branch("fJetPhotonRecHitEta", "std::vector<float>",&fJetPhotonRecHitEta);
  // llpTree->Branch("fJetPhotonRecHitPhi", "std::vector<float>",&fJetPhotonRecHitPhi);
  // llpTree->Branch("fJetPhotonRecHitE", "std::vector<float>",&fJetPhotonRecHitE);
  // llpTree->Branch("fJetPhotonRecHitTime", "std::vector<float>",&fJetPhotonRecHitTime);

  cout << "BRANCHES\n";
  enablePVTracksBranches();
  enableFatJetBranches();
  enableMCBranches();
  enableGenParticleBranches();
  enableCaloJetBranches();
  enableMuonSystemBranches();
  if (enableTriggerInfo_) enableTriggerBranches();
  if (isQCD_)enableQCDBranches();
};
void jet_timing_studies::enablePVTracksBranches()
{
  llpTree->Branch("nPVTracks", &nPVTracks,"nPVTracks/I");
  llpTree->Branch("pvTrackPt", pvTrackPt,"pvTrackPt[nPVTracks]/F");
  llpTree->Branch("pvTrackEta", pvTrackEta,"pvTrackEta[nPVTracks]/F");
  llpTree->Branch("pvTrackPhi", pvTrackPhi,"pvTrackPhi[nPVTracks]/F");
};
void jet_timing_studies::enableMuonSystemBranches()
{

    // csc_Phi = new std::vector<float>;
    // csc_Eta = new std::vector<float>;
    // csc_X = new std::vector<float>;
    // csc_Y = new std::vector<float>;
    // csc_Z = new std::vector<float>;
    // csc_NRecHits = new std::vector<float>;
    // csc_T = new std::vector<float>;
    // csc_Chi2 = new std::vector<float>;

    llpTree->Branch("nCsc",&nCsc,"nCsc/I");
    llpTree->Branch("cscPhi",cscPhi,"cscPhi[nCsc]");
    llpTree->Branch("cscEta",cscEta,"cscEta[nCsc]");
    llpTree->Branch("cscX",cscX,"cscX[nCsc]");
    llpTree->Branch("cscY",cscY,"cscY[nCsc]");
    llpTree->Branch("cscZ",cscZ,"cscZ[nCsc]");
    llpTree->Branch("cscNRecHits",cscNRecHits,"cscNRecHits[nCsc]");
    llpTree->Branch("cscT",cscT,"cscT[nCsc]");
    llpTree->Branch("cscChi2",cscChi2,"cscChi2[nCsc]");

    llpTree->Branch("nRpc",&nRpc,"nRpc/I");
    llpTree->Branch("rpcPhi",rpcPhi,"rpcPhi[nRpc]");
    llpTree->Branch("rpcEta",rpcEta,"rpcEta[nRpc]");
    llpTree->Branch("rpcX",rpcX,"rpcX[nRpc]");
    llpTree->Branch("rpcY",rpcY,"rpcY[nRpc]");
    llpTree->Branch("rpcZ",rpcZ,"rpcZ[nRpc]");
    llpTree->Branch("rpcT",rpcT,"rpcT[nRpc]");
    llpTree->Branch("rpcTError",rpcTError,"rpcTError[nRpc]");

    llpTree->Branch("nDt",&nDt,"nDt/I");
    llpTree->Branch("dtPhi",dtPhi,"dtPhi[nDt]");
    llpTree->Branch("dtEta",dtEta,"dtEta[nDt]");
    llpTree->Branch("dtX",dtX,"dtX[nDt]");
    llpTree->Branch("dtY",dtY,"dtY[nDt]");
    llpTree->Branch("dtZ",dtZ,"dtZ[nDt]");
    llpTree->Branch("dtDirX",dtDirX,"dtDirX[nDt]");
    llpTree->Branch("dtDirY",dtDirY,"dtDirY[nDt]");
    llpTree->Branch("dtDirZ",dtDirZ,"dtDirZ[nDt]");
    llpTree->Branch("dtT",dtT,"dtT[nDt]");
    llpTree->Branch("dtTError",dtTError,"dtTError[nDt]");
};
void jet_timing_studies::enableFatJetBranches()
{
  llpTree->Branch("n_fat_Jets", &n_fat_Jets,"n_fat_Jets/I");
  llpTree->Branch("fat_jetE", fat_jetE,"fat_jetE[n_fat_Jets]/F");
  llpTree->Branch("fat_jetPt", fat_jetPt,"fat_jetPt[n_fat_Jets]/F");
  llpTree->Branch("fat_jetEta", fat_jetEta,"fat_jetEta[n_fat_Jets]/F");
  llpTree->Branch("fat_jetPhi", fat_jetPhi,"fat_jetPhi[n_fat_Jets]/F");
  llpTree->Branch("fat_jetCISV", fat_jetCISV,"fat_jetCISV[n_fat_Jets]/F");
  llpTree->Branch("fat_jetMass", fat_jetMass, "fat_jetMass[n_fat_Jets]/F");
  llpTree->Branch("fat_jetJetArea", fat_jetJetArea, "fat_jetJetArea[n_fat_Jets]/F");
  llpTree->Branch("fat_jetPileupE", fat_jetPileupE, "fat_jetPileupE[n_fat_Jets]/F");
  llpTree->Branch("fat_jetPileupId", fat_jetPileupId, "fat_jetPileupId[n_fat_Jets]/F");
  llpTree->Branch("fat_jetPileupIdFlag", fat_jetPileupIdFlag, "fat_jetPileupIdFlag[n_fat_Jets]/I");
  llpTree->Branch("fat_jetPassIDLoose", fat_jetPassIDLoose, "fat_jetPassIDLoose[n_fat_Jets]/O");
  llpTree->Branch("fat_jetPassIDTight", fat_jetPassIDTight, "fat_jetPassIDTight[n_fat_Jets]/O");
  llpTree->Branch("fat_jetPassMuFrac", fat_jetPassMuFrac, "fat_jetPassMuFrac[n_fat_Jets]/O");
  llpTree->Branch("fat_jetPassEleFrac", fat_jetPassEleFrac, "fat_jetPassEleFrac[n_fat_Jets]/O");
  llpTree->Branch("fat_jetPartonFlavor", fat_jetPartonFlavor, "fat_jetPartonFlavor[n_fat_Jets]/I");
  llpTree->Branch("fat_jetHadronFlavor", fat_jetHadronFlavor, "fat_jetHadronFlavor[n_fat_Jets]/I");
  llpTree->Branch("fat_jetChargedEMEnergyFraction", fat_jetChargedEMEnergyFraction, "fat_jetChargedEMEnergyFraction[n_fat_Jets]/F");
  llpTree->Branch("fat_jetNeutralEMEnergyFraction", fat_jetNeutralEMEnergyFraction, "fat_jetNeutralEMEnergyFraction[n_fat_Jets]/F");
  llpTree->Branch("fat_jetChargedHadronEnergyFraction", fat_jetChargedHadronEnergyFraction, "fat_jetChargedHadronEnergyFraction[n_fat_Jets]/F");
  llpTree->Branch("fat_jetNeutralHadronEnergyFraction", fat_jetNeutralHadronEnergyFraction, "fat_jetNeutralHadronEnergyFraction[n_fat_Jets]/F");
  llpTree->Branch("fat_jet_charged_hadron_multiplicity", fat_jet_charged_hadron_multiplicity, "fat_jet_charged_hadron_multiplicity[n_fat_Jets]/F");
  llpTree->Branch("fat_jet_neutral_hadron_multiplicity", fat_jet_neutral_hadron_multiplicity, "fat_jet_neutral_hadron_multiplicity[n_fat_Jets]/F");
  llpTree->Branch("fat_jet_photon_multiplicity", fat_jet_photon_multiplicity, "fat_jet_photon_multiplicity[n_fat_Jets]/F");
  llpTree->Branch("fat_jet_electron_multiplicity", fat_jet_electron_multiplicity, "fat_jet_electron_multiplicity[n_fat_Jets]/F");
  llpTree->Branch("fat_jet_muon_multiplicity", fat_jet_muon_multiplicity, "fat_jet_muon_multiplicity[n_fat_Jets]/F");
  llpTree->Branch("fat_jet_HF_hadron_multiplicity", fat_jet_HF_hadron_multiplicity, "fat_jet_HF_hadron_multiplicity[n_fat_Jets]/F");
  llpTree->Branch("fat_jet_HF_em_multiplicity", fat_jet_HF_em_multiplicity, "fat_jet_HF_em_multiplicity[n_fat_Jets]/F");
  llpTree->Branch("fat_jet_charged_multiplicity", fat_jet_charged_multiplicity, "fat_jet_charged_multiplicity[n_fat_Jets]/F");
  llpTree->Branch("fat_jet_neutral_multiplicity", fat_jet_neutral_multiplicity, "fat_jet_neutral_multiplicity[n_fat_Jets]/F");
  llpTree->Branch("fat_jetMatchedGenPt", fat_jetMatchedGenPt,"fat_jetMatchedGenPt[n_fat_Jets]/F");
  llpTree->Branch("fat_jetMatchedGenEta", fat_jetMatchedGenEta,"fat_jetMatchedGenEta[n_fat_Jets]/F");
  llpTree->Branch("fat_jetMatchedGenPhi", fat_jetMatchedGenPhi,"fat_jetMatchedGenPhi[n_fat_Jets]/F");
  llpTree->Branch("fat_jetMatchedGenMass", fat_jetMatchedGenMass, "fat_jetMatchedGenMass[n_fat_Jets]/F");

  llpTree->Branch("fat_jet_rechit_E_Ecut3", fat_jet_rechit_E_Ecut3, "fat_jet_rechit_E_Ecut3[n_fat_Jets]/F");
  llpTree->Branch("fat_jet_rechit_T_Ecut3", fat_jet_rechit_T_Ecut3, "fat_jet_rechit_T_Ecut3[n_fat_Jets]/F");
  llpTree->Branch("fat_jet_rechit_E_Ecut4", fat_jet_rechit_E_Ecut4, "fat_jet_rechit_E_Ecut4[n_fat_Jets]/F");
  llpTree->Branch("fat_jet_rechit_T_Ecut4", fat_jet_rechit_T_Ecut4, "fat_jet_rechit_T_Ecut4[n_fat_Jets]/F");
  llpTree->Branch("fat_jet_rechit_E_Ecut2", fat_jet_rechit_E_Ecut2, "fat_jet_rechit_E_Ecut2[n_fat_Jets]/F");
  llpTree->Branch("fat_jet_rechit_T_Ecut2", fat_jet_rechit_T_Ecut2, "fat_jet_rechit_T_Ecut2[n_fat_Jets]/F");
  llpTree->Branch("fat_jet_rechit_E_Ecut1p5", fat_jet_rechit_E_Ecut1p5, "fat_jet_rechit_E_Ecut1p5[n_fat_Jets]/F");
  llpTree->Branch("fat_jet_rechit_T_Ecut1p5", fat_jet_rechit_T_Ecut1p5, "fat_jet_rechit_T_Ecut1p5[n_fat_Jets]/F");
  llpTree->Branch("fat_jet_rechit_E_Ecut1", fat_jet_rechit_E_Ecut1, "fat_jet_rechit_E_Ecut1[n_fat_Jets]/F");
  llpTree->Branch("fat_jet_rechit_T_Ecut1", fat_jet_rechit_T_Ecut1, "fat_jet_rechit_T_Ecut1[n_fat_Jets]/F");
  llpTree->Branch("fat_jet_rechit_E_Ecut0p5", fat_jet_rechit_E_Ecut0p5, "fat_jet_rechit_E_Ecut0p5[n_fat_Jets]/F");
  llpTree->Branch("fat_jet_rechit_T_Ecut0p5", fat_jet_rechit_T_Ecut0p5, "fat_jet_rechit_T_Ecut0p5[n_fat_Jets]/F");
  llpTree->Branch("fat_jet_rechit_E", fat_jet_rechit_E, "fat_jet_rechit_E[n_fat_Jets]/F");
  llpTree->Branch("fat_jet_rechit_T", fat_jet_rechit_T, "fat_jet_rechit_T[n_fat_Jets]/F");

  llpTree->Branch("fat_jet_pv_rechit_T_Ecut3", fat_jet_pv_rechit_T_Ecut3, "fat_jet_rechit_T_Ecut3[n_fat_Jets]/F");
  llpTree->Branch("fat_jet_pv_rechit_T_Ecut4", fat_jet_pv_rechit_T_Ecut4, "fat_jet_rechit_T_Ecut4[n_fat_Jets]/F");
  llpTree->Branch("fat_jet_pv_rechit_T_Ecut2", fat_jet_pv_rechit_T_Ecut2, "fat_jet_rechit_T_Ecut2[n_fat_Jets]/F");
  llpTree->Branch("fat_jet_pv_rechit_T_Ecut1p5", fat_jet_pv_rechit_T_Ecut1p5, "fat_jet_rechit_T_Ecut1p5[n_fat_Jets]/F");
  llpTree->Branch("fat_jet_pv_rechit_T_Ecut1", fat_jet_pv_rechit_T_Ecut1, "fat_jet_rechit_T_Ecut1[n_fat_Jets]/F");
  llpTree->Branch("fat_jet_pv_rechit_T_Ecut0p5", fat_jet_pv_rechit_T_Ecut0p5, "fat_jet_pv_rechit_T_Ecut0p5[n_fat_Jets]/F");
  llpTree->Branch("fat_jet_pv_rechit_T", fat_jet_pv_rechit_T, "fat_jet_rechit_T[n_fat_Jets]/F");

  if( enableRecHitInfo_ )
  {
    llpTree->Branch("fat_jet_n_rechits", fat_jet_n_rechits, "fat_jet_n_rechits[n_fat_Jets]/I");
    llpTree->Branch("fat_jet_rechits_E", fat_jet_rechits_E, "fat_jet_rechits_E[n_fat_Jets][1000]/F");
    llpTree->Branch("fat_jet_rechits_T", fat_jet_rechits_T, "fat_jet_rechits_T[n_fat_Jets][1000]/F");
    llpTree->Branch("fat_jet_pv_rechits_T", fat_jet_pv_rechits_T, "fat_jet_rechits_T[n_fat_Jets][1000]/F");
  }


  return;
};

void jet_timing_studies::enableMCBranches(){
  llpTree->Branch("nGenJets", &nGenJets, "nGenJets/I");
  llpTree->Branch("genJetE", genJetE, "genJetE[nGenJets]/F");
  llpTree->Branch("genJetPt", genJetPt, "genJetPt[nGenJets]/F");
  llpTree->Branch("genJetEta", genJetEta, "genJetEta[nGenJets]/F");
  llpTree->Branch("genJetPhi", genJetPhi, "genJetPhi[nGenJets]/F");
  llpTree->Branch("genJetME", genJetME, "genJetME[nGenJets]/F");
  llpTree->Branch("genMetPtCalo", &genMetPtCalo, "genMetPtCalo/F");
  llpTree->Branch("genMetPhiCalo", &genMetPhiCalo, "genMetPhiCalo/F");
  llpTree->Branch("genMetECalo", &genMetECalo, "genMetECalo/F");
  llpTree->Branch("genMetEtaCalo", &genMetEtaCalo, "genMetEtaCalo/F");
  llpTree->Branch("genMetPtTrue", &genMetPtTrue, "genMetPtTrue/F");
  llpTree->Branch("genMetPhiTrue", &genMetPhiTrue, "genMetPhiTrue/F");
  llpTree->Branch("genMetETrue", &genMetETrue, "genMetETrue/F");
  llpTree->Branch("genMetEtaTrue", &genMetEtaTrue, "genMeTEtaTrue/F");
  llpTree->Branch("genVertexX", &genVertexX, "genVertexX/F");
  llpTree->Branch("genVertexY", &genVertexY, "genVertexY/F");
  llpTree->Branch("genVertexZ", &genVertexZ, "genVertexZ/F");
  llpTree->Branch("genVertexT", &genVertexT, "genVertexT/F");
  llpTree->Branch("genWeight", &genWeight, "genWeight/F");
  llpTree->Branch("genSignalProcessID", &genSignalProcessID, "genSignalProcessID/i");
  llpTree->Branch("genQScale", &genQScale, "genQScale/F");
  llpTree->Branch("genAlphaQCD", &genAlphaQCD, "genAlphaQCD/F");
  llpTree->Branch("genAlphaQED", &genAlphaQED, "genAlphaQED/F");
  llpTree->Branch("genJet_match_jet_index", genJet_match_jet_index, "genJet_match_jet_index[nGenJets]/i");
  llpTree->Branch("genJet_min_delta_r_match_jet", genJet_min_delta_r_match_jet, "genJet_min_delta_r_match_jet[nGenJets]/F");

  /*scaleWeights = new std::vector<float>; scaleWeights->clear();
  pdfWeights = new std::vector<float>; pdfWeights->clear();
  alphasWeights = new std::vector<float>; alphasWeights->clear();
  if (isFastsim_) {
    llpTree->Branch("lheComments", "std::string",&lheComments);
  }
  llpTree->Branch("scaleWeights", "std::vector<float>",&scaleWeights);
  llpTree->Branch("pdfWeights", "std::vector<float>",&pdfWeights);
  llpTree->Branch("alphasWeights", "std::vector<float>",&alphasWeights);
  */
};
void jet_timing_studies::enableQCDBranches()
{
  //QCD BRANCHES
  llpTree->Branch("nGenQCDParticles", &nGenQCDParticles, "nGenQCDParticles/I");
  llpTree->Branch("genQCD_e", genQCD_e, "genQCD_e[nGenQCDParticles]/F");
  llpTree->Branch("genQCD_EB", genQCD_EB, "genQCD_EB[nGenQCDParticles]/O");
  llpTree->Branch("genQCD_ETL", genQCD_ETL, "genQCD_ETL[nGenQCDParticles]/O");
  llpTree->Branch("genQCD_pt", genQCD_pt, "genQCD_pt[nGenQCDParticles]/F");
  llpTree->Branch("genQCD_eta", genQCD_eta, "genQCD_eta[nGenQCDParticles]/F");
  llpTree->Branch("genQCD_phi", genQCD_phi, "genQCD_phi[nGenQCDParticles]/F");
  llpTree->Branch("genQCD_prod_vertex_x", genQCD_prod_vertex_x, "genQCD_prod_vertex_x[nGenQCDParticles]/F");
  llpTree->Branch("genQCD_prod_vertex_y", genQCD_prod_vertex_y, "genQCD_prod_vertex_y[nGenQCDParticles]/F");
  llpTree->Branch("genQCD_prod_vertex_z", genQCD_prod_vertex_y, "genQCD_prod_vertex_z[nGenQCDParticles]/F");
  llpTree->Branch("genQCD_decay_vertex_x", genQCD_decay_vertex_x, "genQCD_decay_vertex_x[nGenQCDParticles]/F");
  llpTree->Branch("genQCD_decay_vertex_y", genQCD_decay_vertex_y, "genQCD_decay_vertex_y[nGenQCDParticles]/F");
  llpTree->Branch("genQCD_decay_vertex_z", genQCD_decay_vertex_y, "genQCD_decay_vertex_z[nGenQCDParticles]/F");
  llpTree->Branch("genQCD_time", genQCD_time, "genQCD_time[nGenQCDParticles]/F");
  llpTree->Branch("genQCD_time_ETL", genQCD_time_ETL, "genQCD_time_ETL[nGenQCDParticles]/F");
  llpTree->Branch("genQCD_travel_time", genQCD_travel_time, "genQCD_travel_time[nGenQCDParticles]/F");
  llpTree->Branch("genQCD_travel_time_ETL", genQCD_travel_time_ETL, "genQCD_travel_time_ETL[nGenQCDParticles]/F");
  llpTree->Branch("genQCD_photon_travel_time", genQCD_photon_travel_time, "genQCD_photon_travel_time[nGenQCDParticles]/F");
  llpTree->Branch("genQCD_photon_travel_time_ETL", genQCD_photon_travel_time_ETL, "genQCD_photon_travel_time_ETL[nGenQCDParticles]/F");
  llpTree->Branch("genParticleQCD_match_jet_index", genParticleQCD_match_jet_index, "genParticleQCD_match_jet_index[nGenQCDParticles]/i");
  llpTree->Branch("genParticleQCD_min_delta_r_match_jet", genParticleQCD_min_delta_r_match_jet, "genParticleQCD_min_delta_r_match_jet[nGenQCDParticles]/F");
};
void jet_timing_studies::enableTriggerBranches()
{
  nameHLT = new std::vector<std::string>; nameHLT->clear();
  llpTree->Branch("HLTDecision", &triggerDecision, ("HLTDecision[" + std::to_string(NTriggersMAX) +  "]/O").c_str());
  //llpTree->Branch("HLTPrescale", &triggerHLTPrescale, ("HLTPrescale[" + std::to_string(NTriggersMAX) +  "]/I").c_str());
  //llpTree->Branch("HLTMR", &HLTMR, "HLTMR/F");
  //llpTree->Branch("HLTRSQ", &HLTRSQ, "HLTRSQ/F");
};
void jet_timing_studies::enableCaloJetBranches()
{
  llpTree->Branch("nCaloJets", &nCaloJets,"nCaloJets/I");
  llpTree->Branch("calojetE", calojetE,"calojetE[nCaloJets]/F");
   llpTree->Branch("calojetEt", calojetEt,"calojetEt[nCaloJets]/F");
  llpTree->Branch("calojetPt", calojetPt,"calojetPt[nCaloJets]/F");
  llpTree->Branch("calojetEta", calojetEta,"calojetEta[nCaloJets]/F");
  llpTree->Branch("calojetPhi", calojetPhi,"calojetPhi[nCaloJets]/F");
  llpTree->Branch("calojet_HadronicEnergyFraction", calojet_HadronicEnergyFraction,"calojet_HadronicEnergyFraction[nCaloJets]/F");
  llpTree->Branch("calojet_EMEnergyFraction", calojet_EMEnergyFraction,"calojet_EMEnergyFraction[nCaloJets]/F");
  llpTree->Branch("calojetGammaMax_EM",calojetGammaMax_EM,"calojetGammaMax_EM[nCaloJets]/F");
  llpTree->Branch("calojetGammaMax_Hadronic",calojetGammaMax_Hadronic,"calojetGammaMax_Hadronic[nCaloJets]/F");
  llpTree->Branch("calojetGammaMax",calojetGammaMax,"calojetGammaMax[nCaloJets]/F");
  llpTree->Branch("calojetGammaMax_ET",calojetGammaMax_ET,"calojetGammaMax_ET[nCaloJets]/F");
  llpTree->Branch("calojetGammaMax_P",calojetGammaMax_P,"calojetGammaMax_P[nCaloJets]/F");

  llpTree->Branch("calojetMass", calojetMass, "calojetMass[nCaloJets]/F");
  llpTree->Branch("calojetAlphaMax",calojetAlphaMax,"calojetAlphaMax[nCaloJets]/F");
  llpTree->Branch("calojetBetaMax",calojetBetaMax,"calojetBetaMax[nCaloJets]/F");
  llpTree->Branch("calojetPtAllTracks",calojetPtAllTracks,"calojetPtAllTracks[nCaloJets]/F");
  llpTree->Branch("calojetPtAllPVTracks",calojetPtAllPVTracks,"calojetPtAllPVTracks[nCaloJets]/F");
  llpTree->Branch("calojetMedianTheta2D",calojetMedianTheta2D,"calojetMedianTheta2D[nCaloJets]/F");
  llpTree->Branch("calojetMedianIP",calojetMedianIP,"calojetMedianIP[nCaloJets]/F");
  llpTree->Branch("calojetMinDeltaRAllTracks",calojetMinDeltaRAllTracks,"calojetMinDeltaRAllTracks[nCaloJets]/F");
  llpTree->Branch("calojetMinDeltaRPVTracks",calojetMinDeltaRPVTracks,"calojetMinDeltaRPVTracks[nCaloJets]/F");
  llpTree->Branch("calojetptPVTracks","std::vector<vector<double>>(nCaloJets)",&calojetptPVTracks);
  llpTree->Branch("calojetnMatchedPVTracks","std::vector<vector<int>>(nCaloJets)",&calojetnMatchedPVTracks);
  llpTree->Branch("calojetJetArea", calojetJetArea, "calojetJetArea[nCaloJets]/F");
  llpTree->Branch("calojetPileupE", calojetPileupE, "calojetPileupE[nCaloJets]/F");
  llpTree->Branch("calojetPileupId", calojetPileupId, "calojetPileupId[nCaloJets]/F");
  llpTree->Branch("calojetPileupIdFlag", calojetPileupIdFlag, "calojetPileupIdFlag[nCaloJets]/I");
  llpTree->Branch("calojetPassIDLoose", calojetPassIDLoose, "calojetPassIDLoose[nCaloJets]/O");
  llpTree->Branch("calojetPassIDTight", calojetPassIDTight, "calojetPassIDTight[nCaloJets]/O");
  // llpTree->Branch("calojetPassMuFrac", calojetPassMuFrac, "calojetPassMuFrac[nCaloJets]/O");
  // llpTree->Branch("calojetPassEleFrac", calojetPassEleFrac, "calojetPassEleFrac[nCaloJets]/O");
  // llpTree->Branch("calojetPartonFlavor", calojetPartonFlavor, "calojetPartonFlavor[nCaloJets]/I");
  // llpTree->Branch("calojetHadronFlavor", calojetHadronFlavor, "calojetHadronFlavor[nCaloJets]/I");
  // llpTree->Branch("calojetChargedEMEnergyFraction", calojetChargedEMEnergyFraction, "calojetChargedEMEnergyFraction[nCaloJets]/F");
  // llpTree->Branch("calojetNeutralEMEnergyFraction", calojetNeutralEMEnergyFraction, "calojetNeutralEMEnergyFraction[nCaloJets]/F");
  // llpTree->Branch("calojetChargedHadronEnergyFraction", calojetChargedHadronEnergyFraction, "calojetChargedHadronEnergyFraction[nCaloJets]/F");
  // llpTree->Branch("calojetNeutralHadronEnergyFraction", calojetNeutralHadronEnergyFraction, "calojetNeutralHadronEnergyFraction[nCaloJets]/F");
  // llpTree->Branch("calojetMuonEnergyFraction", calojetMuonEnergyFraction, "calojetMuonEnergyFraction[nCaloJets]/F");
  // llpTree->Branch("calojetHOEnergyFraction", calojetHOEnergyFraction, "calojetHOEnergyFraction[nCaloJets]/F");
  // llpTree->Branch("calojetHFHadronEnergyFraction", calojetHFHadronEnergyFraction, "calojetHFHadronEnergyFraction[nCaloJets]/F");
  // llpTree->Branch("calojetHFEMEnergyFraction",calojetHFEMEnergyFraction, "calojetHFEMEnergyFraction[nCaloJets]/F");
  // llpTree->Branch("calojetAllMuonPt", calojetAllMuonPt,"calojetAllMuonPt[nCaloJets]/F");
  // llpTree->Branch("calojetAllMuonEta", calojetAllMuonEta,"calojetAllMuonEta[nCaloJets]/F");
  // llpTree->Branch("calojetAllMuonPhi", calojetAllMuonPhi,"calojetAllMuonPhi[nCaloJets]/F");
  // llpTree->Branch("calojetAllMuonM", calojetAllMuonM,"calojetAllMuonM[nCaloJets]/F");
  // llpTree->Branch("calojetPtWeightedDZ", calojetPtWeightedDZ,"calojetPtWeightedDZ[nCaloJets]/F");
  llpTree->Branch("calojetNRechits", calojetNRechits,"calojetNRechits[nCaloJets]/I");
  llpTree->Branch("calojetRechitE", calojetRechitE,"calojetRechitE[nCaloJets]/F");
  llpTree->Branch("calojetRechitT", calojetRechitT,"calojetRechitT[nCaloJets]/F");
  llpTree->Branch("calojetRechitT_rms", calojetRechitT_rms,"calojetRechitT_rms[nCaloJets]/F");

  llpTree->Branch("calojet_match_track_index",calojet_match_track_index,"calojet_match_track_index[nCaloJets]/i");
  llpTree->Branch("calojet_min_delta_r_match_track",calojet_min_delta_r_match_track,"calojet_min_delta_r_match_track[nCaloJets]/F");


};

void jet_timing_studies::enableGenParticleBranches(){
  llpTree->Branch("gLLP_prod_vertex_x", gLLP_prod_vertex_x, "gLLP_prod_vertex_x[2]/F");
  llpTree->Branch("gLLP_prod_vertex_y", gLLP_prod_vertex_y, "gLLP_prod_vertex_y[2]/F");
  llpTree->Branch("gLLP_prod_vertex_z", gLLP_prod_vertex_z, "gLLP_prod_vertex_z[2]/F");
  llpTree->Branch("gLLP_decay_vertex_x", gLLP_decay_vertex_x, "gLLP_decay_vertex_x[2]/F");
  llpTree->Branch("gLLP_decay_vertex_y", gLLP_decay_vertex_y, "gLLP_decay_vertex_y[2]/F");
  llpTree->Branch("gLLP_decay_vertex_z", gLLP_decay_vertex_z, "gLLP_decay_vertex_z[2]/F");
  llpTree->Branch("gLLP_beta", gLLP_beta, "gLLP_beta[2]/F");
  llpTree->Branch("gLLP_e", gLLP_e, "gLLP_e[2]/F");
  llpTree->Branch("gLLP_pt", gLLP_pt, "gLLP_pt[2]/F");
  llpTree->Branch("gLLP_eta", gLLP_eta, "gLLP_eta[2]/F");
  llpTree->Branch("gLLP_phi", gLLP_phi, "gLLP_phi[2]/F");

  llpTree->Branch("gLLP_travel_time", gLLP_travel_time, "gLLP_travel_time[2]/F");

  llpTree->Branch("gLLP_daughter_travel_time", gLLP_daughter_travel_time, "gLLP_daughter_travel_time[4]/F");
  llpTree->Branch("gLLP_daughter_travel_time_ETL", gLLP_daughter_travel_time_ETL, "gLLP_daughter_travel_time_ETL[4]/F");
  llpTree->Branch("gLLP_daughter_pt", gLLP_daughter_pt, "gLLP_daughter_pt[4]/F");
  llpTree->Branch("gLLP_daughter_pz", gLLP_daughter_pz, "gLLP_daughter_pz[4]/F");
  llpTree->Branch("gLLP_daughter_eta", gLLP_daughter_eta, "gLLP_daughter_eta[4]/F");
  llpTree->Branch("gLLP_daughter_phi", gLLP_daughter_phi, "gLLP_daughter_phi[4]/F");
  llpTree->Branch("gLLP_daughter_eta_ecalcorr", gLLP_daughter_eta_ecalcorr, "gLLP_daughter_eta_ecalcorr[4]/F");
  llpTree->Branch("gLLP_daughter_phi_ecalcorr", gLLP_daughter_phi_ecalcorr, "gLLP_daughter_phi_ecalcorr[4]/F");
  llpTree->Branch("gLLP_daughter_eta_hcalcorr", gLLP_daughter_eta_hcalcorr, "gLLP_daughter_eta_hcalcorr[4]/F");
  llpTree->Branch("gLLP_daughter_phi_hcalcorr", gLLP_daughter_phi_hcalcorr, "gLLP_daughter_phi_hcalcorr[4]/F");
  llpTree->Branch("gLLP_daughter_e", gLLP_daughter_e, "gLLP_daughter_e[4]/F");
  llpTree->Branch("gLLP_daughter_EB", gLLP_daughter_EB, "gLLP_daughter_EB[4]/O");
  llpTree->Branch("gLLP_daughter_ETL", gLLP_daughter_ETL, "gLLP_daughter_ETL[4]/O");
  llpTree->Branch("photon_travel_time", photon_travel_time, "photon_travel_time[4]/F");
  llpTree->Branch("photon_travel_time_ETL", photon_travel_time_ETL, "photon_travel_time_ETL[4]/F");
  llpTree->Branch("gen_time", gen_time, "gen_time[4]/F");
  llpTree->Branch("gen_time_ETL", gen_time_ETL, "gen_time_ETL[4]/F");
  llpTree->Branch("gen_time_pv", gen_time_pv, "gen_time_pv[4]/F");


  llpTree->Branch("gLLP_daughter_match_genJet_index", gLLP_daughter_match_genJet_index, "gLLP_daughter_match_genJet_index[4]/i");
  llpTree->Branch("gLLP_min_delta_r_match_genJet", gLLP_min_delta_r_match_genJet, "gLLP_min_delta_r_match_genJet[4]/F");
  llpTree->Branch("gLLP_daughter_match_jet_index_hcal", gLLP_daughter_match_jet_index_hcal, "gLLP_daughter_match_jet_index_hcal[4]/i");
  llpTree->Branch("gLLP_min_delta_r_match_jet_hcal", gLLP_min_delta_r_match_jet_hcal, "gLLP_min_delta_r_match_jet_hcal[4]/F");
  llpTree->Branch("gLLP_daughter_match_jet_index_hcal_loose", gLLP_daughter_match_jet_index_hcal_loose, "gLLP_daughter_match_jet_index_hcal_loose[4]/i");
  llpTree->Branch("gLLP_min_delta_r_match_jet_hcal_loose", gLLP_min_delta_r_match_jet_hcal_loose, "gLLP_min_delta_r_match_jet_hcal_loose[4]/F");
  llpTree->Branch("gLLP_daughter_match_jet_index_loose", gLLP_daughter_match_jet_index_loose, "gLLP_daughter_match_jet_index_loose[4]/i");
  llpTree->Branch("gLLP_min_delta_r_match_jet_loose", gLLP_min_delta_r_match_jet_loose, "gLLP_min_delta_r_match_jet_loose[4]/F");
  llpTree->Branch("gLLP_daughter_match_jet_index", gLLP_daughter_match_jet_index, "gLLP_daughter_match_jet_index[4]/i");
  llpTree->Branch("gLLP_min_delta_r_match_jet", gLLP_min_delta_r_match_jet, "gLLP_min_delta_r_match_jet[4]/F");
  llpTree->Branch("gLLP_daughter_match_calojet_index", gLLP_daughter_match_calojet_index, "gLLP_daughter_match_calojet_index[4]/i");
  llpTree->Branch("gLLP_min_delta_r_match_calojet", gLLP_min_delta_r_match_calojet, "gLLP_min_delta_r_match_calojet[4]/F");
  llpTree->Branch("gLLP_min_delta_r_nocorr_match_jet", gLLP_min_delta_r_nocorr_match_jet, "gLLP_min_delta_r_nocorr_match_jet[4]/F");


  llpTree->Branch("nGenParticle", &nGenParticle, "nGenParticle/I");
  llpTree->Branch("gParticleMotherId", gParticleMotherId, "gParticleMotherId[nGenParticle]/I");
  llpTree->Branch("gParticleMotherIndex", gParticleMotherIndex, "gParticleMotherIndex[nGenParticle]/I");
  llpTree->Branch("gParticleId", gParticleId, "gParticleId[nGenParticle]/I");
  llpTree->Branch("gParticleStatus", gParticleStatus, "gParticleStatus[nGenParticle]/I");
  llpTree->Branch("gParticleE", gParticleE, "gParticleE[nGenParticle]/F");
  llpTree->Branch("gParticlePt", gParticlePt, "gParticlePt[nGenParticle]/F");
  llpTree->Branch("gParticlePx", gParticlePx, "gParticlePx[nGenParticle]/F");
  llpTree->Branch("gParticlePy", gParticlePy, "gParticlePy[nGenParticle]/F");
  llpTree->Branch("gParticlePz", gParticlePz, "gParticlePz[nGenParticle]/F");
  llpTree->Branch("gParticleEta", gParticleEta, "gParticleEta[nGenParticle]/F");
  llpTree->Branch("gParticlePhi", gParticlePhi, "gParticlePhi[nGenParticle]/F");
  llpTree->Branch("gParticleDecayVertexX", gParticleDecayVertexX, "gParticleDecayVertexX[nGenParticle]/F");
  llpTree->Branch("gParticleDecayVertexY", gParticleDecayVertexY, "gParticleDecayVertexY[nGenParticle]/F");
  llpTree->Branch("gParticleDecayVertexZ", gParticleDecayVertexZ, "gParticleDecayVertexZ[nGenParticle]/F");
}



//------ Load the miniAOD objects and reset tree variables for each event ------//
void jet_timing_studies::loadEvent(const edm::Event& iEvent){//load all miniAOD objects for the current event
  iEvent.getByToken(triggerBitsToken_, triggerBits);
  iEvent.getByToken(hepMCToken_, hepMC);
  iEvent.getByToken(triggerBitsToken_, triggerBits);
  iEvent.getByToken(metFilterBitsToken_, metFilterBits);
  iEvent.getByToken(verticesToken_, vertices);
  iEvent.getByToken(cscSegmentInputToken_,cscSegments);
  iEvent.getByToken(dtSegmentInputToken_,dtSegments);
  iEvent.getByToken(rpcRecHitInputToken_,rpcRecHits);
  iEvent.getByToken(tracksTag_,tracks);
  iEvent.getByToken(PFCandsToken_, pfCands);
  iEvent.getByToken(PFClustersToken_, pfClusters);
  iEvent.getByToken(muonsToken_, muons);
  iEvent.getByToken(electronsToken_, electrons);
  iEvent.getByToken(photonsToken_, photons);
  iEvent.getByToken(tausToken_, taus);
  iEvent.getByToken(jetsCaloToken_, jetsCalo);
  iEvent.getByToken(jetsToken_, jets);
  iEvent.getByToken(jetsPuppiToken_, jetsPuppi);
  iEvent.getByToken(jetsAK8Token_, jetsAK8);
  iEvent.getByToken(genMetCaloToken_, genMetsCalo);
  iEvent.getByToken(genMetTrueToken_, genMetsTrue);
  iEvent.getByToken(metToken_, mets);
  //iEvent.getByToken(metNoHFToken_, metsNoHF);
  iEvent.getByToken(metPuppiToken_, metsPuppi);
//  iEvent.getByToken(hcalNoiseInfoToken_,hcalNoiseInfo);
  iEvent.getByToken(secondaryVerticesToken_,secondaryVertices);
  iEvent.getByToken(rhoAllToken_,rhoAll);
  iEvent.getByToken(rhoFastjetAllToken_,rhoFastjetAll);
  iEvent.getByToken(rhoFastjetAllCaloToken_,rhoFastjetAllCalo);
  iEvent.getByToken(rhoFastjetCentralCaloToken_,rhoFastjetCentralCalo);
  iEvent.getByToken(rhoFastjetCentralChargedPileUpToken_,rhoFastjetCentralChargedPileUp);
  iEvent.getByToken(rhoFastjetCentralNeutralToken_,rhoFastjetCentralNeutral);
  iEvent.getByToken(beamSpotToken_,beamSpot);
  iEvent.getByToken(ebRecHitsToken_,ebRecHits);
  iEvent.getByToken(eeRecHitsToken_,eeRecHits);
  iEvent.getByToken(esRecHitsToken_,esRecHits);
  iEvent.getByToken(ebeeClustersToken_,ebeeClusters);
  iEvent.getByToken(esClustersToken_,esClusters);
  iEvent.getByToken(conversionsToken_,conversions);
  iEvent.getByToken(singleLegConversionsToken_,singleLegConversions);
  iEvent.getByToken(gedGsfElectronCoresToken_,gedGsfElectronCores);
  iEvent.getByToken(gedPhotonCoresToken_, gedPhotonCores);
  iEvent.getByToken(generalTrackToken_,generalTracks);
//  iEvent.getByToken(superClustersToken_,superClusters);
//  iEvent.getByToken(lostTracksToken_,lostTracks);
//  iEvent.getByToken(hbheNoiseFilterToken_, hbheNoiseFilter);
//  iEvent.getByToken(hbheTightNoiseFilterToken_, hbheTightNoiseFilter);
//  iEvent.getByToken(hbheIsoNoiseFilterToken_, hbheIsoNoiseFilter);
  //iEvent.getByToken(badChargedCandidateFilterToken_, badChargedCandidateFilter);
  //iEvent.getByToken(badMuonFilterToken_, badMuonFilter);
  if(readGenVertexTime_) iEvent.getByToken(genParticles_t0_Token_,genParticles_t0);
  if (useGen_) {
//    iEvent.getByToken(genParticlesToken_,genParticles);
    iEvent.getByToken(genParticlesToken_,genParticles);
    iEvent.getByToken(genJetsToken_,genJets);

    //for Spring16 fastsim, this has been changed and removed
//    if (!isFastsim_) iEvent.getByToken(lheInfoToken_, lheInfo);

    iEvent.getByToken(genInfoToken_,genInfo);
    iEvent.getByToken(puInfoToken_,puInfo);
  }


}

//called by the loadEvent() method
void jet_timing_studies::resetBranches(){
    //reset tree variables
    reset_event_variables();
    resetPVTracksBranches();

    reset_photon_variable();
    reset_jet_variables();
    resetCaloJetBranches();
    reset_gen_llp_variable();
    reset_gen_jet_variable();
    reset_qcd_variables();
    resetMuonSystemBranches();
}

void jet_timing_studies::reset_event_variables()
{
  eventNum = 0;
  lumiNum = 0;
  runNum = 0;
  pvX = -99.0;
  pvY = -99.0;
  pvZ = -99.0;
  nPV = -1;
  Rho = -99.0;
  nPUmean = -1;
  nPU = -1;
  for(int i = 0; i < OBJECTARRAYSIZE; i++)
  {
    PV_x[i]  = -999.;
    PV_y[i] = -999.;
    PV_z[i] = -999.;
  }
  return;
};
void jet_timing_studies::resetPVTracksBranches()
{
  nPVTracks = 0;
  for(int i = 0; i < OBJECTARRAYSIZE; i++)
  {
    pvTrackPt[i]  = -999.;
    pvTrackEta[i] = -999.;
    pvTrackPhi[i] = -999.;
  }
};
void jet_timing_studies::resetMuonSystemBranches()
{
    nCsc = 0;
    for ( int i = 0; i < OBJECTARRAYSIZE; i++)
    {
      cscPhi[i] = 0.0;
      cscEta[i] = 0.0;
      cscX[i] = 0.0;
      cscY[i] = 0.0;
      cscZ[i] = 0.0;
      cscNRecHits[i] = 0.0;
      cscT[i] = 0.0;
      cscChi2[i] = 0.0;
    }
    nRpc = 0;
    for ( int i = 0; i < OBJECTARRAYSIZE; i++)
    {
      rpcPhi[i] = 0.0;
      rpcEta[i] = 0.0;
      rpcX[i] = 0.0;
      rpcY[i] = 0.0;
      rpcZ[i] = 0.0;
      rpcT[i] = 0.0;
      rpcTError[i] = 0.0;
    }
    nDt = 0;
    for ( int i = 0; i < OBJECTARRAYSIZE; i++)
    {
      dtPhi[i] = 0.0;
      dtEta[i] = 0.0;
      dtX[i] = 0.0;
      dtY[i] = 0.0;
      dtZ[i] = 0.0;
      dtDirX[i] = 0.0;
      dtDirY[i] = 0.0;
      dtDirZ[i] = 0.0;
      dtT[i] = 0.0;
      dtTError[i] = 0.0;
    }
    return;
};
void jet_timing_studies::findTrackingVariables(const TLorentzVector &jetVec,const edm::EventSetup& iSetup,vector<double>& ptPVTrack, vector<int>& nMatchedPVTracks, float &pPVTracksMax, float &alphaMax,float &medianTheta2D,float &medianIP, int &nTracksPV,float &ptAllPVTracks,float &ptAllTracks,float &minDeltaRAllTracks, float &minDeltaRPVTracks)
{
  // int ptPVTracksMax_pvindex = 99;
  int nTracksAll = 0;
  //Displaced jet stuff
  double ptPVTracksMax = 0.;
  minDeltaRAllTracks = 15;
  minDeltaRPVTracks = 15;
  reco::Vertex primaryVertex = vertices->at(0);
  std::vector<double> theta2Ds;
  std::vector<double> IP2Ds;



  for (unsigned int iTrack = 0; iTrack < generalTracks->size(); iTrack ++){
    reco::Track generalTrack = generalTracks->at(iTrack);
    TLorentzVector generalTrackVecTemp;
    generalTrackVecTemp.SetPtEtaPhiM(generalTrack.pt(),generalTrack.eta(),generalTrack.phi(),0);

    if (generalTrack.pt() > 1) {
      if (minDeltaRAllTracks > generalTrackVecTemp.DeltaR(jetVec))
      {
    	    minDeltaRAllTracks =  generalTrackVecTemp.DeltaR(jetVec);
      }
      if (generalTrackVecTemp.DeltaR(jetVec) < 0.4){
      	nTracksAll ++;
      	//tot pt for alpha
      	ptAllTracks += generalTrack.pt();

    		// theta 2d
    		// ROOT::Math::XYZPoint innerPos = generalTrack.innerPosition();
    		// ROOT::Math::XYZPoint vertexPos = primaryVertex.position();
    		// ROOT::Math::XYZVector deltaPos = innerPos - vertexPos;
    		// ROOT::Math::XYZVector momentum = generalTrack.innerMomentum();
    		// double mag2DeltaPos = TMath::Sqrt((deltaPos.x()*deltaPos.x()) + (deltaPos.y()*deltaPos.y()));
    		// double mag2Mom = TMath::Sqrt((momentum.x()*momentum.x()) + (momentum.y()*momentum.y()));
    		// double theta2D = TMath::ACos((deltaPos.x()*momentum.x()+deltaPos.y()*momentum.y())/(mag2Mom*mag2DeltaPos));
    		// theta2Ds.push_back(theta2D);

    		//IP sig
    		edm::ESHandle<TransientTrackBuilder> theB;
    		iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);
    		reco::TransientTrack transTrack = theB->build(generalTrack);
    		TrajectoryStateClosestToBeamLine traj = transTrack.stateAtBeamLine();
    		Measurement1D meas = traj.transverseImpactParameter();
    		std::pair<bool, Measurement1D> ip2d = IPTools::absoluteTransverseImpactParameter(transTrack,primaryVertex);
    		IP2Ds.push_back(ip2d.second.value()/ip2d.second.error());
    		// IP2Ds.push_back(fabs(generalTrack.dxy()/generalTrack.dxyError()));
      }
     }
    }
    // int i = 0;
    if (ptAllTracks > 0.9){
	//No matched jets

	   for (auto vertex = vertices->begin(); vertex != vertices->end(); vertex++){
      double ptPVTracks = 0.;
      double pPVTracks = 0.;
      int nTracksPVTemp = 0;
      if(!vertex->isValid())continue;
      if (vertex->isFake())continue;
	    for(auto pvTrack=vertex->tracks_begin(); pvTrack!=vertex->tracks_end(); pvTrack++){
    		TLorentzVector pvTrackVecTemp;
    		pvTrackVecTemp.SetPtEtaPhiM((*pvTrack)->pt(),(*pvTrack)->eta(),(*pvTrack)->phi(),0);
    		//If pv track associated with jet add pt to ptPVTracks
	      if ((*pvTrack)->pt() > 1) {
  		    if (minDeltaRPVTracks > pvTrackVecTemp.DeltaR(jetVec))
  		    {
			       minDeltaRPVTracks =  pvTrackVecTemp.DeltaR(jetVec);
  		    }
	        if (pvTrackVecTemp.DeltaR(jetVec) < 0.4){
            pPVTracks += (*pvTrack)->p();
      			ptPVTracks += (*pvTrack)->pt();
      			ptAllPVTracks += (*pvTrack)->pt();
      			nTracksPVTemp++;
	        }
	      }



	    }
      ptPVTrack.push_back(ptPVTracks);//pt sum for each vertex
      nMatchedPVTracks.push_back(nTracksPVTemp);//number of matched tracks for each vertex


      if(pPVTracks > pPVTracksMax){
        pPVTracksMax = pPVTracks;
      }
	    if (ptPVTracks > ptPVTracksMax) {
    		ptPVTracksMax = ptPVTracks;
    		nTracksPV = nTracksPVTemp;
        // ptPVTracksMax_pvindex = i;
	    }
	    alphaMax = ptPVTracksMax/ptAllTracks;
      // i++;

	   }
     // std::cout<<nPV<<", "<<vertex_i<<std::endl;
    }
    // std::cout<<"alphamax index: "<<ptPVTracksMax_pvindex<<std::endl;
    std::sort(IP2Ds.begin(),IP2Ds.end());
    if (IP2Ds.size() > 0){
     medianIP = IP2Ds[IP2Ds.size()/2];

    }
    std::sort(theta2Ds.begin(),theta2Ds.end());
    if (theta2Ds.size() > 0){
     medianTheta2D = theta2Ds[theta2Ds.size()/2];
    }
};

void jet_timing_studies::reset_photon_variable()
{
  fJetNPhotons = 0;
  for (int i=0; i< OBJECTARRAYSIZE; i++) {
    fJetPhotonPt[i] = 0.0;
    fJetPhotonEta[i] = 0.0;
    fJetPhotonPhi[i] = 0.0;
    fJetPhotonSeedRecHitE[i]      = -99.0;
    fJetPhotonSeedRecHitEta[i]      = -99.0;
    fJetPhotonSeedRecHitPhi[i]      = -99.0;
    fJetPhotonSeedRecHitTime[i]      = -99.0;
  }
  return;
};
void jet_timing_studies::resetCaloJetBranches()
{
  nCaloJets = 0;
  calojetptPVTracks.clear();
  calojetnMatchedPVTracks.clear();
  for ( int i = 0; i < OBJECTARRAYSIZE; i++)
  {
    calojetE[i] = 0.0;
    calojetPt[i] = 0.0;
    calojetEta[i] = 0.0;
    calojetEt[i] = 0.0;
    calojetPhi[i] = 0.0;
    // calojetCSV[i] = 0.0;
    // calojetCISV[i] = 0.0;
    calojetMass[i] =  -99.0;
    calojetAlphaMax[i] = -99.0;
    calojetBetaMax[i] = -99.0;
    calojetGammaMax[i] = -99.0;
    calojetGammaMax_ET[i] = -99.0;
    calojetGammaMax_P[i] = -99.0;
    calojetGammaMax_EM[i] = -99.0;
    calojetGammaMax_Hadronic[i] = -99.0;

    calojetPtAllTracks[i] = -99.0;
    calojetPtAllPVTracks[i] = -99.0;
    calojetMedianTheta2D[i] = -99.0;
    calojetMedianIP[i] = -99.0;
    calojetMinDeltaRAllTracks[i] =-99.0;
    calojetMinDeltaRPVTracks[i] = -99.0;
    calojetJetArea[i] = -99.0;
    calojetPileupE[i] = -99.0;
    calojetPileupId[i] = -99.0;
    calojetPileupIdFlag[i] = -1;
    calojetPassIDLoose[i] = false;
    calojetPassIDTight[i] = false;
    calojet_match_track_index[i] = 666;
    calojet_min_delta_r_match_track[i] = -666.;
    calojet_HadronicEnergyFraction[i] = -666.;
    calojet_EMEnergyFraction[i] = -666.;

    // calojetPassMuFrac[i] = false;
    // calojetPassEleFrac[i] = false;
    // calojetPartonFlavor[i] = 0;
    // calojetHadronFlavor[i] = 0;
    // calojetChargedEMEnergyFraction[i] = -99.0;
    // calojetNeutralEMEnergyFraction[i] = -99.0;
    // calojetChargedHadronEnergyFraction[i] = -99.0;
    // calojetNeutralHadronEnergyFraction[i] = -99.0;
    // calojetMuonEnergyFraction[i] = -99.0;
    // calojetHOEnergyFraction[i] = -99.0;
    // calojetHFHadronEnergyFraction[i] = -99.0;
    // calojetHFEMEnergyFraction[i] = -99.0;
    // calojetAllMuonPt[i] = 0.0;
    // calojetAllMuonEta[i] = 0.0;
    // calojetAllMuonPhi[i] = 0.0;
    // calojetAllMuonM[i] = 0.0;
    // calojetPtWeightedDZ[i] = 0.0;
    calojetNRechits[i] = 0;
    calojetRechitE[i] = 0.0;
    calojetRechitT[i] = 0.0;
    calojetRechitT_rms[i] = 0.0;

  }
  return;
};


void jet_timing_studies::reset_jet_variables()
{
  nJets = 0;
  jetptPVTracks.clear();
  jetnMatchedPVTracks.clear();
  jet_rechits_E.clear();
  jet_rechits_T.clear();
  jet_rechits_phi.clear();
  jet_rechits_eta.clear();

  for ( int i = 0; i < OBJECTARRAYSIZE; i++)
  {
    jetE[i] = 0.0;
    jetEt[i] = 0.0;
    jetPt[i] = 0.0;
    jetEta[i] = 0.0;
    jetPhi[i] = 0.0;
    jetCISV[i] = 0.0;
    jetMass[i] =  -99.0;
    jetGammaMax[i] = -99.0;
    jetGammaMax_ET[i] = -99.0;
    jetGammaMax_P[i] = -99.0;
    jetGammaMax_EM[i] = -99.0;
    jetGammaMax_Hadronic[i] = -99.0;
    jetAlphaMax[i] = -99.0;
    jetBetaMax[i] = -99.0;
    jetPtAllTracks[i] = -99.0;
    jetPtAllPVTracks[i] = -99.0;
    jetMedianTheta2D[i] = -99.0;
    jetMedianIP[i] = -99.0;
    jetMinDeltaRAllTracks[i] =-99.0;
    jetMinDeltaRPVTracks[i] = -99.0;
    jetJetArea[i] = -99.0;
    jetPileupE[i] = -99.0;
    jetPileupId[i] = -99.0;
    jetPileupIdFlag[i] = -1;
    jetPassIDLoose[i] = false;
    jetPassIDTight[i] = false;
    jetPassMuFrac[i] = false;
    jetPassEleFrac[i] = false;
    jetPartonFlavor[i] = 0;
    jetHadronFlavor[i] = 0;
    jetChargedEMEnergyFraction[i] = -99.0;
    jetNeutralEMEnergyFraction[i] = -99.0;
    jetChargedHadronEnergyFraction[i] = -99.0;
    jetNeutralHadronEnergyFraction[i] = -99.0;
    jet_charged_hadron_multiplicity[i] = -99;
    jet_neutral_hadron_multiplicity[i] = -99;
    jet_photon_multiplicity[i] = -99;
    jet_electron_multiplicity[i] = -99;
    jet_muon_multiplicity[i] = -99;
    jet_HF_hadron_multiplicity[i] = -99;
    jet_HF_em_multiplicity[i] = -99;
    jet_charged_multiplicity[i] = -99;
    jet_neutral_multiplicity[i] = -99;
    jetMatchedGenPt[i] = 0.0;
    jetMatchedGenEta[i] = 0.0;
    jetMatchedGenPhi[i] = 0.0;
    jetMatchedGenMass[i] = 0.0;
    jetMatchedGenTime[i] = 0.0;
    jet_n_rechits[i] = 0;
    event_n_rechits = 0;
    jet_n_rechits_Ecut1[i] = 0;
    jet_n_rechits_Ecut0p5[i] = 0;
    good_jet[i] = true;
    good_jet0p5[i] = true;
    jet_rechit_E[i] = 0.0;
    jet_rechit_T[i] = 0.0;
    jet_rechit_T_rms_Ecut1[i] = 0.0;
    jet_rechit_T_rms_Ecut0p5[i] = 0.0;
    jet_rechit_E_Ecut3[i] = 0.0; //energy with a 2 GeV cut
    jet_rechit_T_Ecut3[i] = 0.0;

    jet_rechit_E_Ecut4[i] = 0.0; //energy with a 2 GeV cut
    jet_rechit_T_Ecut4[i] = 0.0;
    jet_rechit_E_Ecut2[i] = 0.0; //energy with a 2 GeV cut
    jet_rechit_T_Ecut2[i] = 0.0;
    jet_rechit_E_Ecut1p5[i] = 0.0; //energy with a 2 GeV cut
    jet_rechit_T_Ecut1p5[i] = 0.0;
    jet_rechit_E_Ecut1[i] = 0.0; //energy with a 2 GeV cut
    jet_rechit_T_Ecut1[i] = 0.0;
    jet_rechit_E_Ecut0p5[i] = 0.0; //energy with a 2 GeV cut
    jet_rechit_T_Ecut0p5[i] = 0.0;

    jet_pv_rechit_T[i] = 0.0;
    jet_pv_rechit_T_Ecut4[i] = 0.0;
    jet_pv_rechit_T_Ecut2[i] = 0.0;
    jet_pv_rechit_T_Ecut1p5[i] = 0.0;
    jet_pv_rechit_T_Ecut1[i] = 0.0;
    jet_pv_rechit_T_Ecut0p5[i] = 0.0;
  }
  // for(int i =0; i < RECHITARRAYSIZE;i++)
  // {
  //   // jet_rechits_E[i][j] = -666.;
  //   // jet_rechits_T[i][j] = -666.;
  //   jet_rechits_E[i] = -666.;
  //   jet_rechits_T[i] = -666.;
  //   jet_rechits_phi[i] = -666.;
  //   jet_rechits_eta[i] = -666.;
  //   jet_pv_rechits_T[i] = -666.;
  //
  // }
  pfMetPt = 0.0;
  pfMetPhi = 0.0;
  pfMetEta = 0.0;
  pfMetE = 0.0;
  return;
};

void jet_timing_studies::reset_gen_llp_variable()
{
  for ( int i = 0; i < LLP_ARRAY_SIZE; i++ )
  {
    gLLP_prod_vertex_x[i] = -666.;
    gLLP_prod_vertex_y[i] = -666.;
    gLLP_prod_vertex_z[i] = -666.;
    gLLP_decay_vertex_x[i] = -666.;
    gLLP_decay_vertex_y[i] = -666.;
    gLLP_decay_vertex_z[i] = -666.;
    gLLP_beta[i] = -666.;
    gLLP_pt[i] = -666.;
    gLLP_e[i] = -666.;
    gLLP_eta[i] = -666.;
    gLLP_phi[i] = -666.;
    gLLP_travel_time[i] = -666.;
  }

  for ( int i = 0; i < LLP_DAUGHTER_ARRAY_SIZE; i++ )
  {
    gLLP_daughter_pt[i] = -666.;
    gLLP_daughter_pz[i] = -666.;
    gLLP_daughter_eta[i] = -666.;
    gLLP_daughter_phi[i] = -666.;
    gLLP_daughter_eta_ecalcorr[i] = -666.;
    gLLP_daughter_phi_ecalcorr[i] = -666.;
    gLLP_daughter_eta_hcalcorr[i] = -666.;
    gLLP_daughter_phi_hcalcorr[i] = -666.;
    gLLP_daughter_e[i] = -666.;
    gLLP_daughter_EB[i] = false;
    gLLP_daughter_ETL[i] = false;
    gLLP_daughter_travel_time[i] = -666.;
    gLLP_daughter_travel_time_ETL[i] = -666.;
    gen_time[i] = -666.;
    gen_time_ETL[i] = -666.;
    gen_time_pv[i] = -666.;
    photon_travel_time[i] = -666.;
    photon_travel_time_ETL[i] = -666.;
    photon_travel_time_pv[i] = -666.;
    gLLP_daughter_match_calojet_index[i] = 666;
    gLLP_daughter_match_jet_index[i] = 666;
    gLLP_daughter_match_jet_index_hcal[i] = 666;
    gLLP_daughter_match_jet_index_hcal_loose[i] = 666;
    gLLP_daughter_match_jet_index_loose[i] = 666;
    gLLP_min_delta_r_match_calojet[i] = -666.;
    gLLP_min_delta_r_match_jet[i] = -666.;
    gLLP_min_delta_r_match_jet_hcal[i] = -666.;
    gLLP_min_delta_r_match_jet_loose[i] = -666.;
    gLLP_min_delta_r_match_jet_hcal_loose[i] = -666.;
    gLLP_min_delta_r_nocorr_match_jet[i] = -666.;
    gLLP_daughter_match_genJet_index[i] = 666;
    gLLP_min_delta_r_match_genJet[i] = -666.;

  }
  return;
};

void jet_timing_studies::reset_gen_jet_variable()
{
  nGenJets = 0;
  for ( int i = 0; i < OBJECTARRAYSIZE; i++ )
  {
    genJetE[i] = -666.;
    genJetPt[i] = -666.;
    genJetEta[i] = -666.;
    genJetPhi[i] = -666.;
    genJetME[i] = -666.;
    genJet_match_jet_index[i] = 666;
    genJet_min_delta_r_match_jet[i] = -666.;
  }
  genMetPtCalo  = -666.;
  genMetPhiCalo = -666.;
  genMetECalo  = -666.;
  genMetEtaCalo = -666.;
  genMetPtTrue  = -666.;
  genMetPhiTrue = -666.;
  genMetETrue  = -666.;
  genMetEtaTrue = -666.;
  return;
};
void jet_timing_studies::reset_qcd_variables()
{
  nGenQCDParticles = 0;
  for (int i = 0; i < GENPARTICLEARRAYSIZE; i++ )
  {
    genQCD_pt[i] = -666;
    genQCD_e[i] = -666;
    genQCD_eta[i] = -666;
    genQCD_phi[i] = -666;
    genQCD_EB[i] = false;
    genQCD_ETL[i] = false;
    genQCD_prod_vertex_x[i] = -666;
    genQCD_prod_vertex_y[i] = -666;
    genQCD_prod_vertex_z[i] = -666;
    genQCD_decay_vertex_x[i] = -666;
    genQCD_decay_vertex_y[i] = -666;
    genQCD_decay_vertex_z[i] = -666;
    genQCD_time[i] = -666;
    genQCD_time_ETL[i] = -666;
    genQCD_travel_time[i] = -666;
    genQCD_travel_time_ETL[i] = -666;
    genQCD_photon_travel_time[i] = -666;
    genQCD_photon_travel_time_ETL[i] = -666;
    genParticleQCD_match_jet_index[i] = 666;
    genParticleQCD_min_delta_r_match_jet[i] = -666.;
  }
  return;
};
//------ Methods to fill tree variables ------//




//------ Method called for each run ------//

void jet_timing_studies::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup) {


}


//------ Method called for each lumi block ------//
void jet_timing_studies::beginLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const&) {

}


//------ Method called for each event ------//

void jet_timing_studies::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  using namespace edm;

  //initialize
  loadEvent(iEvent); //loads objects and resets tree branches
  NEvents->Fill(0); //increment event count
  //resetting output tree branches
  resetBranches();
  //*************************************
  //Fill Event-Level Info
  //*************************************

  //store basic event info
  isData = isData_;
  model = model_;
  isQCD = isQCD_;
  runNum = iEvent.id().run();
  lumiNum = iEvent.luminosityBlock();
  eventNum = iEvent.id().event();

 //select the primary vertex, if any
  nPV = 0;
  myPV = &(vertices->front());

  bool foundPV = false;
  for(unsigned int i = 0; i < vertices->size(); i++)
  {
    if(vertices->at(i).isValid() && !vertices->at(i).isFake())
    {
      if (!foundPV)
      {
        myPV = &(vertices->at(i));
        // std::cout<<"pv index: "<<i<<std::endl;
        foundPV = true;
      }
      PV_x[nPV] = vertices->at(i).x();
      PV_y[nPV] = vertices->at(i).y();
      PV_z[nPV] = vertices->at(i).z();
      nPV++;
    }
  }

  pvX = myPV->x();
  pvY = myPV->y();
  pvZ = myPV->z();

  //get rho
  Rho = *rhoFastjetAll;

  //Fill Pileup info
  if (!isData)
  {
    for(const PileupSummaryInfo &pu : *puInfo)
    {
      if ( pu.getBunchCrossing() == 0)
      {
        nPU = pu.getPU_NumInteractions();
        nPUmean = pu.getTrueNumInteractions();
      }
    }
  }

  int i_jet = 0;
  int n_rechits = 0;
  for (const reco::PFJet &j : *jets)
  {
    //resetBranches();
    if (j.pt() < 20) continue;
    if (fabs(j.eta()) > 2.4) continue;



    //*************************************
    //Fill Jet-Level Info
    //*************************************
    jetE[i_jet] = j.energy();
    jetEt[i_jet] = j.et();
    jetPt[i_jet] = j.pt();
    jetEta[i_jet] = j.eta();
    jetPhi[i_jet] = j.phi();
    jetMass[i_jet] = j.mass();

    TLorentzVector thisJet;
    thisJet.SetPtEtaPhiE(jetPt[i_jet], jetEta[i_jet], jetPhi[i_jet], jetE[i_jet]);
    float alphaMax(0.0),medianTheta2D(0.0),medianIP(0.0),minDeltaRAllTracks(0.0),minDeltaRPVTracks(0.0),ptAllTracks(0.0), ptAllPVTracks(0.0);
    int nTracksPV(0);
    float pPVTracksMax(0.0);
    std::vector<double> ptPVTrack;
    std::vector<int> nMatchedPVTracks;


    findTrackingVariables(thisJet,iSetup,ptPVTrack,nMatchedPVTracks, pPVTracksMax,alphaMax,medianTheta2D,medianIP,nTracksPV,ptAllPVTracks,ptAllTracks, minDeltaRAllTracks, minDeltaRPVTracks);
    //jetCISV = j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
    jetAlphaMax[i_jet] = alphaMax;
    jetBetaMax[i_jet] = alphaMax * ptAllTracks/(j.pt());
    jetGammaMax[i_jet] = alphaMax * ptAllTracks/(j.energy());
    jetGammaMax_EM[i_jet] = alphaMax * ptAllTracks/(j.energy()*(j.chargedEmEnergyFraction()+j.neutralEmEnergyFraction()));
    jetGammaMax_Hadronic[i_jet] = alphaMax * ptAllTracks/(j.energy()*(j.chargedHadronEnergyFraction()+j.neutralHadronEnergyFraction()));
    jetGammaMax_ET[i_jet] = alphaMax * ptAllTracks/(j.et());
    jetGammaMax_P[i_jet] = pPVTracksMax/j.energy();

    jetMedianTheta2D[i_jet] = medianTheta2D;
    jetMedianIP[i_jet] = medianIP;
    jetPtAllPVTracks[i_jet] = ptAllPVTracks;
    jetPtAllTracks[i_jet] = ptAllTracks;
    jetMinDeltaRAllTracks[i_jet] = minDeltaRAllTracks;
    jetMinDeltaRPVTracks[i_jet] = minDeltaRPVTracks;
    jetptPVTracks.push_back(ptPVTrack);
    jetnMatchedPVTracks.push_back(nMatchedPVTracks);

    jetJetArea[i_jet] = j.jetArea();
    jetPileupE[i_jet] = j.pileup();

    jetPileupIdFlag[i_jet] = 0;
    jetPassIDLoose[i_jet] = passJetID(&j, 0);
    jetPassIDTight[i_jet] = passJetID(&j, 1);
    jetPassMuFrac[i_jet]  = ( j.muonEnergyFraction() < 0.80 );
    jetPassEleFrac[i_jet]  = ( j.electronEnergyFraction() < 0.90 );


    // if (useGen_) {
    //   jetPartonFlavor = j.partonFlavour();
    //   jetHadronFlavor = j.hadronFlavour();
    // }

    jetChargedEMEnergyFraction[i_jet] = j.chargedEmEnergyFraction();
    jetNeutralEMEnergyFraction[i_jet] = j.neutralEmEnergyFraction();
    jetChargedHadronEnergyFraction[i_jet] = j.chargedHadronEnergyFraction();
    jetNeutralHadronEnergyFraction[i_jet] = j.neutralHadronEnergyFraction();
    jet_charged_hadron_multiplicity[i_jet] = j.chargedHadronMultiplicity();
    jet_neutral_hadron_multiplicity[i_jet] = j.neutralHadronMultiplicity();
    jet_photon_multiplicity[i_jet] = j.photonMultiplicity();
    jet_electron_multiplicity[i_jet] = j.electronMultiplicity();
    jet_muon_multiplicity[i_jet] = j.muonMultiplicity();
    jet_HF_hadron_multiplicity[i_jet] = j.HFHadronMultiplicity();
    jet_HF_em_multiplicity[i_jet] = j.HFEMMultiplicity();
    jet_charged_multiplicity[i_jet] = j.chargedMultiplicity();
    jet_neutral_multiplicity[i_jet] = j.neutralMultiplicity();


    //*************************************
    //find photons inside the jet
    //*************************************
    for (const reco::Photon &pho : *photons) {
      //cout << "Nphoton: " << fJetNPhotons << "\n";

      if (!(deltaR(pho.eta(), pho.phi() , j.eta(), j.phi()) < 0.5)) continue;


      fJetPhotonPt[fJetNPhotons]  = pho.pt();
      fJetPhotonEta[fJetNPhotons] = pho.eta(); //correct this for the vertex
      fJetPhotonPhi[fJetNPhotons] = pho.phi(); //correct this for the vertex

      fJetPhotonSeedRecHitE[fJetNPhotons]      = pho.superCluster()->seed()->x();
      fJetPhotonSeedRecHitEta[fJetNPhotons]      = pho.superCluster()->seed()->y();
      fJetPhotonSeedRecHitPhi[fJetNPhotons]      = pho.superCluster()->seed()->z();
      fJetPhotonSeedRecHitTime[fJetNPhotons]      = pho.superCluster()->seed()->energy();

      // //get time coordinate for the seed
      // for (const reco::PFCluster &pfcluster : *pfClusters) {
      // 	if(pfcluster.seed() == pho.superCluster()->seed()->seed())
      // 	  {
      // 	    pho_superClusterSeedT[fJetNPhotons] = pfcluster.time();
      // 	    pho_pfClusterSeedE[fJetNPhotons]      = pfcluster.energy();
      // 	  }
      // }

      //*************************************
      //fill all rechits inside photons
      //*************************************

      fJetNPhotons++;

    }
    //***************************
    //Find RecHits Inside the Jet
    //***************************
    // geometry (from ECAL ELF)

    edm::ESHandle<CaloGeometry> geoHandle;
    iSetup.get<CaloGeometryRecord>().get(geoHandle);
    const CaloSubdetectorGeometry *barrelGeometry = geoHandle->getSubdetectorGeometry(DetId::Ecal, EcalBarrel);
    //const CaloSubdetectorGeometry *endcapGeometry = geoHandle->getSubdetectorGeometry(DetId::Ecal, EcalEndcap);
    double ecal_radius = 129.0;
    double hcal_radius = 179.0;
    int n_matched_rechits = 0;
    int n_matched_rechits_Ecut0p5 = 0;
    int n_matched_rechits_Ecut1 = 0;
    std::vector<double> rechitphi;
    std::vector<double> rechiteta;
    std::vector<double> rechite;
    std::vector<double> rechitt;

    for (EcalRecHitCollection::const_iterator recHit = ebRecHits->begin(); recHit != ebRecHits->end(); ++recHit)
    {
      // if (recHit->checkFlag(EcalRecHit::kSaturated) || recHit->checkFlag(EcalRecHit::kLeadingEdgeRecovered) || recHit->checkFlag(EcalRecHit::kPoorReco) || recHit->checkFlag(EcalRecHit::kWeird) || recHit->checkFlag(EcalRecHit::kDiWeird)) continue;
      // if (recHit->timeError() < 0 || recHit->timeError() > 100) continue;

      const DetId recHitId = recHit->detid();
      const auto recHitPos = barrelGeometry->getGeometry(recHitId)->getPosition();

      if ( deltaR(jetEta[i_jet], jetPhi[i_jet], recHitPos.eta(), recHitPos.phi())  < 0.4)
      {
        // if (good_jet[i_jet])
        // {
        //   std::cout << "good jet: " << good_jet[i_jet] << ","<<i_jet<<","<<nJets<<std::endl;
        //   std::cout<< "rechit energy" <<recHit->energy() <<std::endl;
        //   std::cout<< "time error flag" << (recHit->timeError() < 0 || recHit->timeError() > 100 )<<std::endl;
        //   std::cout<< "flag" << (recHit->checkFlag(EcalRecHit::kSaturated) || recHit->checkFlag(EcalRecHit::kLeadingEdgeRecovered) || recHit->checkFlag(EcalRecHit::kPoorReco) || recHit->checkFlag(EcalRecHit::kWeird) || recHit->checkFlag(EcalRecHit::kDiWeird))<<std::endl;
        //   std::cout<< "flag"<< (!recHit->checkFlag(0)) <<std::endl;
        // }

        if ( !recHit->checkFlag(0)){
          good_jet[i_jet] = false;
          if (recHit->energy()>0.5){
            good_jet0p5[i_jet] = false;
          }
          continue;
        }
        if (recHit->timeError() < 0 || recHit->timeError() > 100)
        {
          good_jet[i_jet] = false;
          if (recHit->energy()>0.5){
            good_jet0p5[i_jet] = false;
          }
          continue;
        }
        if (recHit->checkFlag(EcalRecHit::kSaturated) || recHit->checkFlag(EcalRecHit::kLeadingEdgeRecovered) || recHit->checkFlag(EcalRecHit::kPoorReco) || recHit->checkFlag(EcalRecHit::kWeird) || recHit->checkFlag(EcalRecHit::kDiWeird))
        {
          good_jet[i_jet] = false;
          if (recHit->energy()>0.5){
            good_jet0p5[i_jet] = false;
          }
          continue;
        }

        jet_rechit_E[i_jet] += recHit->energy();
        jet_rechit_T[i_jet] += recHit->time()*recHit->energy();
        // if (i_jet != 0 && (recHit->energy() >= 1.0 ||recHit->energy() == 0.0  )){
        //   std::cout << "before: "<< i_jet <<", "<<n_matched_rechits<<", "<< n_matched_rechits_Ecut1<<", "<<jet_rechits_E[i_jet][n_matched_rechits] << ", " <<recHit->energy()<< std::endl;
        // }
        // jet_rechits_E[i_jet][n_matched_rechits] = recHit->energy();
  	    // jet_rechits_T[i_jet][n_matched_rechits] = recHit->time();
        // jet_rechits_phi[n_rechits] = recHitPos.phi();
        // jet_rechits_eta[n_rechits] = recHitPos.eta();
        // jet_rechits_E[n_rechits]= recHit->energy();
        // jet_rechits_T[n_rechits] = recHit->time();
        rechitphi.push_back(recHitPos.phi());
        rechiteta.push_back(recHitPos.eta());
        rechite.push_back(recHit->energy());
        rechitt.push_back(recHit->time());

        // if (i_jet != 0 && (recHit->energy() >= 1.0||recHit->energy() == 0.0)){
        //   std::cout << "after: "<< i_jet <<", "<< n_matched_rechits<<", "<<n_matched_rechits_Ecut1<<", "<<jet_rechits_E[i_jet][n_matched_rechits] << ", " <<recHit->energy()<< std::endl;
        // }


  	    double rechit_x = ecal_radius * cos(recHitPos.phi());
  	    double rechit_y = ecal_radius * sin(recHitPos.phi());
  	    double rechit_z = ecal_radius * sinh(recHitPos.eta());
  	    double photon_pv_travel_time = (1./30) * sqrt(pow(pvX-rechit_x,2)+pow(pvY-rechit_y,2)+pow(pvZ-rechit_z,2));
        jet_pv_rechits_T[n_rechits] = recHit->time()+(1./30)*ecal_radius*cosh(recHitPos.eta()) - photon_pv_travel_time;
  	    jet_pv_rechit_T[i_jet] += recHit->energy()*jet_pv_rechits_T[n_rechits];
        if (recHit->energy() > 0.5)
  	    {
      		jet_rechit_E_Ecut0p5[i_jet] += recHit->energy();
      		jet_rechit_T_Ecut0p5[i_jet] += recHit->time()*recHit->energy();
          jet_pv_rechit_T_Ecut0p5[i_jet] += jet_pv_rechits_T[n_rechits] *recHit->energy();
          jet_rechit_T_rms_Ecut0p5[i_jet] += recHit->time()*recHit->time();
          n_matched_rechits_Ecut0p5++;
  	    }
        if (recHit->energy() > 1.0)
  	    {
  		    jet_rechit_E_Ecut1[i_jet] += recHit->energy();
  		    jet_rechit_T_Ecut1[i_jet] += recHit->time()*recHit->energy();
          jet_pv_rechit_T_Ecut1[i_jet] += jet_pv_rechits_T[n_rechits] *recHit->energy();
          jet_rechit_T_rms_Ecut1[i_jet] += recHit->time()*recHit->time();
          // std::cout << "rechit time, with pv"<<jet_rechit_T_Ecut1[i_jet]<< jet_pv_rechit_T_Ecut1[i_jet]<< std::endl;
          // std::cout << "rechits with pv, without" <<jet_pv_rechits_T[i_jet][n_matched_rechits] << jet_rechits_T[i_jet][n_matched_rechits] << std::endl;
          // std::cout << "rechit energy and time"<<recHit->energy()<< ", "<<recHit->time()<< std::endl;
          n_matched_rechits_Ecut1++;
  	    }
        if (recHit->energy() > 1.5)
  	    {
  		    jet_rechit_E_Ecut1p5[i_jet] += recHit->energy();
  		    jet_rechit_T_Ecut1p5[i_jet] += recHit->time()*recHit->energy();
          jet_pv_rechit_T_Ecut1p5[i_jet] += jet_pv_rechits_T[n_rechits] *recHit->energy();

  	    }
        if (recHit->energy() > 2.0)
  	    {
  		    jet_rechit_E_Ecut2[i_jet] += recHit->energy();
  		    jet_rechit_T_Ecut2[i_jet] += recHit->time()*recHit->energy();
          jet_pv_rechit_T_Ecut2[i_jet] += jet_pv_rechits_T[n_rechits] *recHit->energy();

  	    }
  	    if (recHit->energy() > 3.0)
        {
          jet_rechit_E_Ecut3[i_jet] += recHit->energy();
          jet_rechit_T_Ecut3[i_jet] += recHit->time()*recHit->energy();
          jet_pv_rechit_T_Ecut3[i_jet] += jet_pv_rechits_T[n_rechits] *recHit->energy();

        }

  	    if (recHit->energy() > 4.0)
        {
          jet_rechit_E_Ecut4[i_jet] += recHit->energy();
          jet_rechit_T_Ecut4[i_jet] += recHit->time()*recHit->energy();
          jet_pv_rechit_T_Ecut4[i_jet] += jet_pv_rechits_T[n_rechits] *recHit->energy();

        }
        n_matched_rechits++;
        n_rechits++;
      }

    }
    //cout << "Last Nphoton: " << fJetNPhotons << "\n";
    //std::cout << "n: " << n_matched_rechits << std::endl;
    jet_n_rechits[i_jet] = n_matched_rechits;
    jet_n_rechits_Ecut1[i_jet] = n_matched_rechits_Ecut1;
    jet_n_rechits_Ecut0p5[i_jet] = n_matched_rechits_Ecut0p5;
    jet_rechits_phi.push_back(rechitphi);
    jet_rechits_eta.push_back(rechiteta);
    jet_rechits_E.push_back(rechite);
    jet_rechits_T.push_back(rechitt);

    if (n_matched_rechits_Ecut1 > 0){
      jet_rechit_T_rms_Ecut1[i_jet] = sqrt(jet_rechit_T_rms_Ecut1[i_jet])/n_matched_rechits_Ecut1;

    }
    else{
      jet_rechit_T_rms_Ecut1[i_jet] = -666.0;

    }
    if (n_matched_rechits_Ecut0p5 > 0){
      jet_rechit_T_rms_Ecut0p5[i_jet] = sqrt(jet_rechit_T_rms_Ecut0p5[i_jet])/n_matched_rechits_Ecut0p5;

    }
    else{
      jet_rechit_T_rms_Ecut0p5[i_jet] = -666.0;

    }
    jet_rechit_T[i_jet] = jet_rechit_T[i_jet]/jet_rechit_E[i_jet];
    jet_rechit_T_Ecut4[i_jet] = jet_rechit_T_Ecut4[i_jet]/jet_rechit_E_Ecut4[i_jet];
    jet_rechit_T_Ecut3[i_jet] = jet_rechit_T_Ecut3[i_jet]/jet_rechit_E_Ecut3[i_jet];
    jet_rechit_T_Ecut2[i_jet] = jet_rechit_T_Ecut2[i_jet]/jet_rechit_E_Ecut2[i_jet];
    jet_rechit_T_Ecut1p5[i_jet] = jet_rechit_T_Ecut1p5[i_jet]/jet_rechit_E_Ecut1p5[i_jet];
    jet_rechit_T_Ecut1[i_jet] =  jet_rechit_T_Ecut1[i_jet]/jet_rechit_E_Ecut1[i_jet];
    jet_rechit_T_Ecut0p5[i_jet] = jet_rechit_T_Ecut0p5[i_jet]/jet_rechit_E_Ecut0p5[i_jet]; //incrementing jet counter
    jet_pv_rechit_T[i_jet] = jet_pv_rechit_T[i_jet]/jet_rechit_E[i_jet];
    jet_pv_rechit_T_Ecut4[i_jet] = jet_pv_rechit_T_Ecut4[i_jet]/jet_rechit_E_Ecut4[i_jet];
    jet_pv_rechit_T_Ecut3[i_jet] = jet_pv_rechit_T_Ecut3[i_jet]/jet_rechit_E_Ecut3[i_jet];
    jet_pv_rechit_T_Ecut2[i_jet] =  jet_pv_rechit_T_Ecut2[i_jet]/jet_rechit_E_Ecut2[i_jet];
    jet_pv_rechit_T_Ecut1p5[i_jet] = jet_pv_rechit_T_Ecut1p5[i_jet]/jet_rechit_E_Ecut1p5[i_jet];
    jet_pv_rechit_T_Ecut1[i_jet] = jet_pv_rechit_T_Ecut1[i_jet]/jet_rechit_E_Ecut1[i_jet];
    jet_pv_rechit_T_Ecut0p5[i_jet] = jet_pv_rechit_T_Ecut0p5[i_jet]/jet_rechit_E_Ecut0p5[i_jet]; //incrementing jet counter
    nJets++;
    i_jet++;


  } //loop over jets
  event_n_rechits = n_rechits;

  const reco::PFMET &Met = mets->front();
  pfMetPt = Met.pt();
  pfMetPhi = Met.phi();
  pfMetE = Met.energy();
  pfMetEta = Met.eta();
  //MC AND GEN LEVEL INFO
  if(!isData){
    fillMC();
    fillGenParticles();

  }
  // fillPVTracks();
  // fillCaloJets( iSetup );
  fillMuonSystem(iEvent, iSetup);
  //fill_fat_jet( iSetup );
  /*if(readGenVertexTime_)
  {
    genVertexT = *genParticles_t0; //std::cout << genVertexT << std::endl;
  }
  */

  if ( enableTriggerInfo_ ) fillTrigger( iEvent );
  llpTree->Fill();
}

//------ Method called once each job just before starting event loop ------//
void jet_timing_studies::beginJob(){
  setBranches();
}

//------ Method called once each job just after ending the event loop ------//
void jet_timing_studies::endJob(){
}
bool jet_timing_studies::fillMuonSystem(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    edm::ESHandle<CSCGeometry> cscG;
    edm::ESHandle<DTGeometry> dtG;
    edm::ESHandle<RPCGeometry> rpcG;

    iSetup.get<MuonGeometryRecord>().get(cscG);
    iSetup.get<MuonGeometryRecord>().get(dtG);
    iSetup.get<MuonGeometryRecord>().get(rpcG);

    for (const CSCSegment cscSegment : *cscSegments) {
	float globPhi   = 0.;
	float globX = 0.;
	float globY = 0.;
	float globZ = 0.;
	float globEta = 0.;
	CSCDetId id  = (CSCDetId)(cscSegment).cscDetId();
	LocalPoint segPos = (cscSegment).localPosition();
	const CSCChamber* cscchamber = cscG->chamber(id);
	if (cscchamber) {
	    GlobalPoint globalPosition = cscchamber->toGlobal(segPos);
	    globPhi   = globalPosition.phi();
	    globEta   = globalPosition.eta();
	    globX = globalPosition.x();
	    globY = globalPosition.y();
	    globZ = globalPosition.z();
	    // globR = pow(globX*globX+globY*globY,0.5);
	    cscNRecHits[nCsc] = cscSegment.nRecHits();
	    cscX[nCsc] = globX;
	    cscY[nCsc] = globY;
	    cscZ[nCsc] = globZ;
	    cscPhi[nCsc] = globPhi;
	    cscEta[nCsc] = globEta;
	    cscT[nCsc] = cscSegment.time();
	    cscChi2[nCsc] = cscSegment.chi2();
	    nCsc++;
	}
    }
    for (const RPCRecHit rpcRecHit : *rpcRecHits){
	LocalPoint  rpcRecHitLocalPosition       = rpcRecHit.localPosition();
	// LocalError  segmentLocalDirectionError = iDT->localDirectionError();
	DetId geoid = rpcRecHit.geographicalId();
	RPCDetId rpcdetid = RPCDetId(geoid);
	const RPCChamber * rpcchamber = rpcG->chamber(rpcdetid);
	if (rpcchamber) {
	    GlobalPoint globalPosition = rpcchamber->toGlobal(rpcRecHitLocalPosition);
	    rpcX[nRpc] = globalPosition.x();
	    rpcY[nRpc] = globalPosition.y();
	    rpcZ[nRpc] = globalPosition.z();
	    rpcPhi[nRpc] = globalPosition.phi();
	    rpcEta[nRpc] = globalPosition.eta();
	    rpcT[nRpc] = rpcRecHit.time();
	    rpcTError[nRpc] = rpcRecHit.timeError();
	    nRpc++;
	}
    }
    for(DTRecSegment4D dtSegment : *dtSegments){
	LocalPoint  segmentLocalPosition       = dtSegment.localPosition();
	LocalVector segmentLocalDirection      = dtSegment.localDirection();
	// LocalError  segmentLocalPositionError  = iDT->localPositionError();
	// LocalError  segmentLocalDirectionError = iDT->localDirectionError();
	DetId geoid = dtSegment.geographicalId();
	DTChamberId dtdetid = DTChamberId(geoid);
	const DTChamber * dtchamber = dtG->chamber(dtdetid);
	if (dtchamber) {
	    GlobalPoint globalPosition = dtchamber->toGlobal(segmentLocalPosition);
	    GlobalVector globalDirection = dtchamber->toGlobal(segmentLocalDirection);

	    dtPhi[nDt] = globalPosition.phi();
	    dtEta[nDt] = globalPosition.eta();
	    dtX[nDt] = globalPosition.x();
	    dtY[nDt] = globalPosition.y();
	    dtZ[nDt] = globalPosition.z();
	    dtDirX[nDt] = globalDirection.x();
	    dtDirY[nDt] = globalDirection.y();
	    dtDirZ[nDt] = globalDirection.z();
	    dtT[nDt] = 0;//dtSegment.time();
	    dtTError[nDt] = -1;//dtSegment.timeError();
	    nDt++;
	}

    }


    return true;
}
bool jet_timing_studies::fillCaloJets(const edm::EventSetup& iSetup)
{
  for (const reco::CaloJet &j : *jetsCalo)
  {
    if (j.pt() < 20) continue;
    if (fabs(j.eta()) > 2.4) continue;
    //-------------------
    //Fill Jet-Level Info
    //-------------------
    calojetE[nCaloJets] = j.energy();
    calojetEt[nCaloJets] = j.et();
    calojetPt[nCaloJets] = j.pt();
    calojetEta[nCaloJets] = j.eta();
    calojetPhi[nCaloJets] = j.phi();
    calojetMass[nCaloJets] = j.mass();
    calojet_HadronicEnergyFraction[nCaloJets] = j.energyFractionHadronic();
    calojet_EMEnergyFraction[nCaloJets] = j.emEnergyFraction();

    TLorentzVector thisJet;
    thisJet.SetPtEtaPhiE(calojetPt[nCaloJets], calojetEta[nCaloJets], calojetPhi[nCaloJets], calojetE[nCaloJets]);
    //calojetCISV = j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
    float alphaMax(0.0),medianTheta2D(0.0),medianIP(0.0),minDeltaRAllTracks(0.0),minDeltaRPVTracks(0.0),ptAllTracks(0.0), ptAllPVTracks(0.0);
    int nTracksPV(0);
    float pPVTracksMax(0.0);
    std::vector<double> ptPVTrack;
    std::vector<int> nMatchedPVTracks;
    // int ptPVTracksMax_pvindex(99);
    // findTrackingVariables(thisJet,iSetup,alphaMax,medianTheta2D,medianIP,nTracksPV,ptPVTracksMax_pvindex,ptAllPVTracks,ptAllTracks, minDeltaRAllTracks, minDeltaRPVTracks);
    // std::cout<<"pv index for alpha max: "<<ptPVTracksMax_pvindex<<" pt: "<<calojetPt[nCaloJets]<<std::endl;
    findTrackingVariables(thisJet,iSetup,ptPVTrack, nMatchedPVTracks, pPVTracksMax,alphaMax,medianTheta2D,medianIP,nTracksPV,ptAllPVTracks,ptAllTracks, minDeltaRAllTracks, minDeltaRPVTracks);
    //jetCISV = j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
    calojetAlphaMax[nCaloJets] = alphaMax;
    calojetBetaMax[nCaloJets] = alphaMax * ptAllTracks/j.pt();
    calojetGammaMax[nCaloJets] = alphaMax * ptAllTracks/(j.energy());
    calojetGammaMax_EM[nCaloJets] = alphaMax * ptAllTracks/(j.energy()*j.emEnergyFraction());
    calojetGammaMax_Hadronic[nCaloJets] =  alphaMax * ptAllTracks/(j.energy()*j.energyFractionHadronic());
    calojetGammaMax_ET[nCaloJets] = alphaMax * ptAllTracks/(j.et());
    calojetGammaMax_P[nCaloJets] = pPVTracksMax/j.energy();
    calojetMedianTheta2D[nCaloJets] = medianTheta2D;
    calojetMedianIP[nCaloJets] = medianIP;
    calojetPtAllPVTracks[nCaloJets] = ptAllPVTracks;
    calojetPtAllTracks[nCaloJets] = ptAllTracks;
    calojetMinDeltaRAllTracks[nCaloJets] = minDeltaRAllTracks;
    calojetMinDeltaRPVTracks[nCaloJets] = minDeltaRPVTracks;
    calojetptPVTracks.push_back(ptPVTrack);
    calojetnMatchedPVTracks.push_back(nMatchedPVTracks);

    // calojetptPVTracks[nCaloJets] = ptPVTrack;
    // calojetnMatchedPVTracks[nCaloJets] = nMatchedPVTracks;

    calojetJetArea[nCaloJets] = j.jetArea();
    calojetPileupE[nCaloJets] = j.pileup();

    calojetPileupIdFlag[nCaloJets] = 0;
    calojetPassIDLoose[nCaloJets] = passCaloJetID(&j, 0);
    calojetPassIDTight[nCaloJets] = passCaloJetID(&j, 1);
    //---------------------------
    //Find PV tracks close to calojet
    //---------------------------


    unsigned int match_track_index = 666;
    double min_delta_r = 666.;

    for (int i_track = 0; i_track < nPVTracks; i_track++)
    {

      double current_delta_r = deltaR(calojetEta[nCaloJets],calojetPhi[nCaloJets] , pvTrackEta[i_track], pvTrackPhi[i_track]);

      if ( current_delta_r < min_delta_r )
      {
        min_delta_r = current_delta_r;
        match_track_index = i_track;
      }
     }//end matching to jets
     if ( min_delta_r < 0.3 )
     {
       calojet_match_track_index[nCaloJets] = match_track_index;
       calojet_min_delta_r_match_track[nCaloJets] = min_delta_r;
     }




    //---------------------------
    //Find RecHits Inside the Jet
    //---------------------------
    // geometry (from ECAL ELF)

    edm::ESHandle<CaloGeometry> geoHandle;
    iSetup.get<CaloGeometryRecord>().get(geoHandle);
    const CaloSubdetectorGeometry *barrelGeometry = geoHandle->getSubdetectorGeometry(DetId::Ecal, EcalBarrel);
    //const CaloSubdetectorGeometry *endcapGeometry = geoHandle->getSubdetectorGeometry(DetId::Ecal, EcalEndcap);
    //double ecal_radius = 129.0;
    int n_matched_rechits = 0;
    for (EcalRecHitCollection::const_iterator recHit = ebRecHits->begin(); recHit != ebRecHits->end(); ++recHit)
    {
      if (recHit->checkFlag(EcalRecHit::kSaturated) || recHit->checkFlag(EcalRecHit::kLeadingEdgeRecovered) || recHit->checkFlag(EcalRecHit::kPoorReco) || recHit->checkFlag(EcalRecHit::kWeird) || recHit->checkFlag(EcalRecHit::kDiWeird)) continue;
      if (recHit->timeError() < 0 || recHit->timeError() > 100) continue;
      if ( recHit->checkFlag(0) )
      {
        const DetId recHitId = recHit->detid();
        const auto recHitPos = barrelGeometry->getGeometry(recHitId)->getPosition();
        if ( deltaR(calojetEta[nCaloJets], calojetPhi[nCaloJets], recHitPos.eta(), recHitPos.phi())  < 0.4)
        {
          //double rechit_x = ecal_radius * cos(recHitPos.phi());
          //double rechit_y = ecal_radius * sin(recHitPos.phi());
          //double rechit_z = ecal_radius * sinh(recHitPos.eta());
          //double photon_pv_travel_time = (1./30) * sqrt(pow(pvX-rechit_x,2)+pow(pvY-rechit_y,2)+pow(pvZ-rechit_z,2));

          if (recHit->energy() > 0.5)
          {
            calojetRechitE[nCaloJets] += recHit->energy();
            calojetRechitT[nCaloJets] += recHit->time()*recHit->energy();
            calojetRechitT_rms[nCaloJets] += recHit->time()*recHit->time();

          }
          n_matched_rechits++;
        }
      }
    }
    //cout << "Last Nphoton: " << fJetNPhotons << "\n";
    //std::cout << "n: " << n_matched_rechits << std::endl;
    calojetNRechits[nCaloJets] = n_matched_rechits;
    calojetRechitT[nCaloJets] = calojetRechitT[nCaloJets]/calojetRechitE[nCaloJets];
    calojetRechitT_rms[nCaloJets] = sqrt(calojetRechitT_rms[nCaloJets]);

    nCaloJets++;
  } //loop over calojets

  return true;
};


bool jet_timing_studies::fill_fat_jet(const edm::EventSetup& iSetup)
{
  int i_fat_jet = 0;
  for (const reco::PFJet &j : *jetsAK8)
  {
    //resetBranches();
    if (j.pt() < 20) continue;
    if (fabs(j.eta()) > 2.4) continue;
    //*************************************
    //Fill Jet-Level Info
    //*************************************
    fat_jetE[i_fat_jet] = j.energy();
    fat_jetPt[i_fat_jet] = j.pt();
    fat_jetEta[i_fat_jet] = j.eta();
    fat_jetPhi[i_fat_jet] = j.phi();
    fat_jetMass[i_fat_jet] = j.mass();

    TLorentzVector thisJet;
    thisJet.SetPtEtaPhiE(fat_jetPt[i_fat_jet], fat_jetEta[i_fat_jet], fat_jetPhi[i_fat_jet], fat_jetE[i_fat_jet]);
    //jetCISV = j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");

    fat_jetJetArea[i_fat_jet] = j.jetArea();
    fat_jetPileupE[i_fat_jet] = j.pileup();

    fat_jetPileupIdFlag[i_fat_jet] = 0;
    fat_jetPassIDLoose[i_fat_jet] = passJetID(&j, 0);
    fat_jetPassIDTight[i_fat_jet] = passJetID(&j, 1);
    fat_jetPassMuFrac[i_fat_jet]  = ( j.muonEnergyFraction() < 0.80 );
    fat_jetPassEleFrac[i_fat_jet]  = ( j.electronEnergyFraction() < 0.90 );


    // if (useGen_) {
    //   fat_jetPartonFlavor = j.partonFlavour();
    //   fat_jetHadronFlavor = j.hadronFlavour();
    // }

    fat_jetChargedEMEnergyFraction[i_fat_jet] = j.chargedEmEnergyFraction();
    fat_jetNeutralEMEnergyFraction[i_fat_jet] = j.neutralEmEnergyFraction();
    fat_jetChargedHadronEnergyFraction[i_fat_jet] = j.chargedHadronEnergyFraction();
    fat_jetNeutralHadronEnergyFraction[i_fat_jet] = j.neutralHadronEnergyFraction();
    fat_jet_charged_hadron_multiplicity[i_fat_jet] = j.chargedHadronMultiplicity();
    fat_jet_neutral_hadron_multiplicity[i_fat_jet] = j.neutralHadronMultiplicity();
    fat_jet_photon_multiplicity[i_fat_jet] = j.photonMultiplicity();
    fat_jet_electron_multiplicity[i_fat_jet] = j.electronMultiplicity();
    fat_jet_muon_multiplicity[i_fat_jet] = j.muonMultiplicity();
    fat_jet_HF_hadron_multiplicity[i_fat_jet] = j.HFHadronMultiplicity();
    fat_jet_HF_em_multiplicity[i_fat_jet] = j.HFEMMultiplicity();
    fat_jet_charged_multiplicity[i_fat_jet] = j.chargedMultiplicity();
    fat_jet_neutral_multiplicity[i_fat_jet] = j.neutralMultiplicity();

    //***************************
    //Find RecHits Inside the Jet
    //***************************
    // geometry (from ECAL ELF)

    edm::ESHandle<CaloGeometry> geoHandle;
    iSetup.get<CaloGeometryRecord>().get(geoHandle);
    const CaloSubdetectorGeometry *barrelGeometry = geoHandle->getSubdetectorGeometry(DetId::Ecal, EcalBarrel);
    //const CaloSubdetectorGeometry *endcapGeometry = geoHandle->getSubdetectorGeometry(DetId::Ecal, EcalEndcap);
    double ecal_radius = 129.0;
    int n_matched_rechits = 0;
    for (EcalRecHitCollection::const_iterator recHit = ebRecHits->begin(); recHit != ebRecHits->end(); ++recHit)
    {
      if (recHit->checkFlag(EcalRecHit::kSaturated) || recHit->checkFlag(EcalRecHit::kLeadingEdgeRecovered) || recHit->checkFlag(EcalRecHit::kPoorReco) || recHit->checkFlag(EcalRecHit::kWeird) || recHit->checkFlag(EcalRecHit::kDiWeird)) continue;
      if (recHit->timeError() < 0 || recHit->timeError() > 100) continue;
      if ( recHit->checkFlag(0) )
      {
        const DetId recHitId = recHit->detid();
        const auto recHitPos = barrelGeometry->getGeometry(recHitId)->getPosition();
        if ( deltaR(fat_jetEta[i_fat_jet], fat_jetPhi[i_fat_jet], recHitPos.eta(), recHitPos.phi())  < 0.4)
        {
          fat_jet_rechit_E[i_fat_jet] += recHit->energy();
          fat_jet_rechit_T[i_fat_jet] += recHit->time()*recHit->energy();
          fat_jet_rechits_E[i_fat_jet][n_matched_rechits] = recHit->energy();
          fat_jet_rechits_T[i_fat_jet][n_matched_rechits] = recHit->time();
          double rechit_x = ecal_radius * cos(recHitPos.phi());
          double rechit_y = ecal_radius * sin(recHitPos.phi());
          double rechit_z = ecal_radius * sinh(recHitPos.eta());
          double photon_pv_travel_time = (1./30) * sqrt(pow(pvX-rechit_x,2)+pow(pvY-rechit_y,2)+pow(pvZ-rechit_z,2));
          fat_jet_pv_rechits_T[i_fat_jet][n_matched_rechits] = recHit->time()+(1./30)*ecal_radius*cosh(recHitPos.eta()) - photon_pv_travel_time;
          fat_jet_pv_rechit_T[i_fat_jet] += recHit->energy()*fat_jet_pv_rechits_T[i_fat_jet][n_matched_rechits];
          // std::cout << fat_jet_pv_rechits_T[i_fat_jet][n_matched_rechits] << fat_jet_rechits_T[i_fat_jet][n_matched_rechits] << std::endl;
          if (recHit->energy() > 0.5)
          {
            fat_jet_rechit_E_Ecut0p5[i_fat_jet] += recHit->energy();
            fat_jet_rechit_T_Ecut0p5[i_fat_jet] += recHit->time()*recHit->energy();
            fat_jet_pv_rechit_T_Ecut0p5[i_fat_jet] += fat_jet_pv_rechits_T[i_fat_jet][n_matched_rechits] *recHit->energy();
          }
          if (recHit->energy() > 1.0)
          {
            fat_jet_rechit_E_Ecut1[i_fat_jet] += recHit->energy();
            fat_jet_rechit_T_Ecut1[i_fat_jet] += recHit->time()*recHit->energy();
            fat_jet_pv_rechit_T_Ecut1[i_fat_jet] += fat_jet_pv_rechits_T[i_fat_jet][n_matched_rechits] *recHit->energy();
            // std::cout << "rechit time, with pv"<<fat_jet_rechit_T_Ecut1[i_fat_jet]<< fat_jet_pv_rechit_T_Ecut1[i_fat_jet]<< std::endl;
            // std::cout << "rechits with pv, without" <<fat_jet_pv_rechits_T[i_fat_jet][n_matched_rechits] << fat_jet_rechits_T[i_fat_jet][n_matched_rechits] << std::endl;
            // std::cout << "rechit energy and time"<<recHit->energy()<< recHit->time()<< std::endl;


          }
          if (recHit->energy() > 1.5)
          {
            fat_jet_rechit_E_Ecut1p5[i_fat_jet] += recHit->energy();
            fat_jet_rechit_T_Ecut1p5[i_fat_jet] += recHit->time()*recHit->energy();
            fat_jet_pv_rechit_T_Ecut1p5[i_fat_jet] += fat_jet_pv_rechits_T[i_fat_jet][n_matched_rechits] *recHit->energy();

          }
          if (recHit->energy() > 2.0)
          {
            fat_jet_rechit_E_Ecut2[i_fat_jet] += recHit->energy();
            fat_jet_rechit_T_Ecut2[i_fat_jet] += recHit->time()*recHit->energy();
            fat_jet_pv_rechit_T_Ecut2[i_fat_jet] += fat_jet_pv_rechits_T[i_fat_jet][n_matched_rechits] *recHit->energy();

          }
          if (recHit->energy() > 3.0)
          {
            fat_jet_rechit_E_Ecut3[i_fat_jet] += recHit->energy();
            fat_jet_rechit_T_Ecut3[i_fat_jet] += recHit->time()*recHit->energy();
            fat_jet_pv_rechit_T_Ecut3[i_fat_jet] += fat_jet_pv_rechits_T[i_fat_jet][n_matched_rechits] *recHit->energy();

          }

          if (recHit->energy() > 4.0)
          {
            fat_jet_rechit_E_Ecut4[i_fat_jet] += recHit->energy();
            fat_jet_rechit_T_Ecut4[i_fat_jet] += recHit->time()*recHit->energy();
            fat_jet_pv_rechit_T_Ecut4[i_fat_jet] += fat_jet_pv_rechits_T[i_fat_jet][n_matched_rechits] *recHit->energy();

          }
          n_matched_rechits++;
        }
      }
    }
    //cout << "Last Nphoton: " << fJetNPhotons << "\n";
    //std::cout << "n: " << n_matched_rechits << std::endl;
    fat_jet_n_rechits[i_fat_jet] = n_matched_rechits;
    fat_jet_rechit_T[i_fat_jet] = fat_jet_rechit_T[i_fat_jet]/fat_jet_rechit_E[i_fat_jet];
    fat_jet_rechit_T_Ecut4[i_fat_jet] = fat_jet_rechit_T_Ecut4[i_fat_jet]/fat_jet_rechit_E_Ecut4[i_fat_jet];
    fat_jet_rechit_T_Ecut3[i_fat_jet] = fat_jet_rechit_T_Ecut3[i_fat_jet]/fat_jet_rechit_E_Ecut3[i_fat_jet];
    fat_jet_rechit_T_Ecut2[i_fat_jet] = fat_jet_rechit_T_Ecut2[i_fat_jet]/fat_jet_rechit_E_Ecut2[i_fat_jet];
    fat_jet_rechit_T_Ecut1p5[i_fat_jet] = fat_jet_rechit_T_Ecut1p5[i_fat_jet]/fat_jet_rechit_E_Ecut1p5[i_fat_jet];
    fat_jet_rechit_T_Ecut1[i_fat_jet] =  fat_jet_rechit_T_Ecut1[i_fat_jet]/fat_jet_rechit_E_Ecut1[i_fat_jet];
    fat_jet_rechit_T_Ecut0p5[i_fat_jet] = fat_jet_rechit_T_Ecut0p5[i_fat_jet]/fat_jet_rechit_E_Ecut0p5[i_fat_jet]; //incrementing fat_jet counter
    fat_jet_pv_rechit_T[i_fat_jet] = fat_jet_pv_rechit_T[i_fat_jet]/fat_jet_rechit_E[i_fat_jet];
    fat_jet_pv_rechit_T_Ecut4[i_fat_jet] = fat_jet_pv_rechit_T_Ecut4[i_fat_jet]/fat_jet_rechit_E_Ecut4[i_fat_jet];
    fat_jet_pv_rechit_T_Ecut3[i_fat_jet] = fat_jet_pv_rechit_T_Ecut3[i_fat_jet]/fat_jet_rechit_E_Ecut3[i_fat_jet];
    fat_jet_pv_rechit_T_Ecut2[i_fat_jet] =  fat_jet_pv_rechit_T_Ecut2[i_fat_jet]/fat_jet_rechit_E_Ecut2[i_fat_jet];
    fat_jet_pv_rechit_T_Ecut1p5[i_fat_jet] = fat_jet_pv_rechit_T_Ecut1p5[i_fat_jet]/fat_jet_rechit_E_Ecut1p5[i_fat_jet];
    fat_jet_pv_rechit_T_Ecut1[i_fat_jet] = fat_jet_pv_rechit_T_Ecut1[i_fat_jet]/fat_jet_rechit_E_Ecut1[i_fat_jet];
    fat_jet_pv_rechit_T_Ecut0p5[i_fat_jet] = fat_jet_pv_rechit_T_Ecut0p5[i_fat_jet]/fat_jet_rechit_E_Ecut0p5[i_fat_jet]; //incrementing fat_jet counter
    n_fat_Jets++;
    i_fat_jet++;

  } //loop over jets
  return true;
};
bool jet_timing_studies::passCaloJetID( const reco::CaloJet *jetCalo, int cutLevel) {
  bool result = false;

  return result;
}//passJetID CaloJet


bool jet_timing_studies::passJetID( const reco::PFJet *jet, int cutLevel) {
  bool result = false;

  double NHF = jet->neutralHadronEnergyFraction();
  double NEMF = jet->neutralEmEnergyFraction();
  int NumConst = jet->chargedMultiplicity() + jet->neutralMultiplicity() ;
  double CHF = jet->chargedHadronEnergyFraction();
  double MUF = jet->muonEnergyFraction();
  double CEMF = jet->chargedEmEnergyFraction();
  int NumNeutralParticles =jet->neutralMultiplicity();
  int CHM = jet->chargedMultiplicity();

  //Loose
  if (cutLevel == 0) {
    if ( fabs(jet->eta()) <= 2.4) {
      if ( NHF  < 0.99 && NEMF < 0.99 && NumConst > 1
	   && CHF > 0 && CHM > 0 && CEMF < 0.99 ) result = true;
    } else if( fabs(jet->eta()) <= 3.0)  {
      if ( NHF  < 0.99 && NEMF < 0.99 && NumConst > 1 ) result = true;
    } else {
      if ( NEMF < 0.90 && NumNeutralParticles > 10 ) result = true;
    }
  }

  //Tight
  else if (cutLevel == 1) {
    if ( fabs(jet->eta()) <= 2.4) {
      if ( NHF  < 0.90 && NEMF < 0.90 && NumConst > 1
	   && CHF > 0 && CHM > 0 && CEMF < 0.99 ) result = true;
    } else if( fabs(jet->eta()) <= 3.0)  {
      if ( NHF  < 0.90 && NEMF < 0.90 && NumConst > 1 ) result = true;
    } else {
      if ( NEMF < 0.90 && NumNeutralParticles > 10 ) result = true;
    }
  }

  //Tight Lep Veto
  else if (cutLevel == 2) {
    if ( fabs(jet->eta()) <= 2.4) {
      if ( NHF  < 0.90 && NEMF < 0.90 && NumConst > 1
	   && CHF > 0 && CHM > 0 && CEMF < 0.99 && MUF < 0.8 ) result = true;
    } else if( fabs(jet->eta()) <= 3.0)  {
      if ( NHF  < 0.90 && NEMF < 0.90 && NumConst > 1 ) result = true;
    } else {
      if ( NEMF < 0.90 && NumNeutralParticles > 10 ) result = true;
    }
  }

  return result;
}

double jet_timing_studies::deltaPhi(double phi1, double phi2)
{
  double dphi = phi1-phi2;
  while (dphi > TMath::Pi())
  {
    dphi -= TMath::TwoPi();
  }
  while (dphi <= -TMath::Pi())
  {
    dphi += TMath::TwoPi();
  }
  return dphi;
};

double jet_timing_studies::deltaR(double eta1, double phi1, double eta2, double phi2)
{
double dphi = deltaPhi(phi1,phi2);
double deta = eta1 - eta2;
return sqrt( dphi*dphi + deta*deta);
};
bool jet_timing_studies::fillPVTracks()
{
  //select the primary vertex, if any
  //myPV = &(vertices->front());
  //bool foundPV = false;
  for(unsigned int i = 0; i < vertices->size(); i++)
  {
    if(vertices->at(i).isValid() && !vertices->at(i).isFake())
    {
      myPV = &(vertices->at(i));
      for(auto pvTrack = myPV->tracks_begin(); pvTrack != myPV->tracks_end(); pvTrack++)
      {
        if( (*pvTrack)->pt() > pvTrack_pt_cut )
        {
          pvTrackPt[nPVTracks]  = (*pvTrack)->pt();
          pvTrackEta[nPVTracks] = (*pvTrack)->eta();
          pvTrackPhi[nPVTracks] = (*pvTrack)->phi();
          nPVTracks++;
        }
      }
    }
  }

  return true;
};

bool jet_timing_studies::fillMC()
{
  for(const reco::GenJet &j : *genJets)
  {
    //std::cout << nGenJets << std::endl;
    genJetE[nGenJets] = j.energy();
    genJetPt[nGenJets] = j.pt();
    genJetEta[nGenJets] = j.eta();
    genJetPhi[nGenJets] = j.phi();
    genJetME[nGenJets] = j.invisibleEnergy();
    nGenJets++;
  }

  const reco::GenMET &GenMetCalo = genMetsCalo->front();
  genMetPtCalo  = GenMetCalo.pt();
  genMetPhiCalo = GenMetCalo.phi();
  genMetEtaCalo = GenMetCalo.eta();
  genMetECalo = GenMetCalo.energy();


  const reco::GenMET &GenMetTrue = genMetsTrue->front();
  genMetPtTrue  = GenMetTrue.pt();
  genMetPhiTrue = GenMetTrue.phi();
  genMetEtaTrue = GenMetTrue.eta();
  genMetETrue = GenMetTrue.energy();

  bool foundGenVertex = false;
  for(size_t i=0; i<genParticles->size();i++)
  {
    if (!foundGenVertex)
    {
      for (unsigned int j=0; j<(*genParticles)[i].numberOfDaughters(); ++j)
      {
        const reco::Candidate *dau = (*genParticles)[i].daughter(j);
        if (dau)
        {
          genVertexX = dau->vx();
          genVertexY = dau->vy();
          genVertexZ = dau->vz();
          if(readGenVertexTime_) genVertexT = *genParticles_t0;
          foundGenVertex = true;
          break;
        }
      }
    }
  }

  genWeight = genInfo->weight();
  genSignalProcessID = genInfo->signalProcessID();
  genQScale = genInfo->qScale();
  genAlphaQCD = genInfo->alphaQCD();
  genAlphaQED = genInfo->alphaQED();
  for ( int i_genJet = 0; i_genJet < nGenJets; i_genJet++ )
  {

    unsigned int match_jet_index = 666;
    double min_delta_r = 666.;

    for (int i_jet = 0; i_jet < nJets; i_jet++)
    {

      double current_delta_r = deltaR(genJetEta[i_genJet],genJetPhi[i_genJet] , jetEta[i_jet], jetPhi[i_jet]);

      if ( current_delta_r < min_delta_r )
      {
        min_delta_r = current_delta_r;
        match_jet_index = i_jet;
      }
     }//end matching to jets
     if ( min_delta_r < 0.3 )
     {
       genJet_match_jet_index[i_genJet] = match_jet_index;
       genJet_min_delta_r_match_jet[i_genJet] = min_delta_r;
     }
   }





    /*
    if (isFastsim_) {

      //get lhe weights for systematic uncertainties:
      double nomlheweight = genInfo->weights()[0];

      //fill scale variation weights
      if (genInfo->weights().size()>=10) {
	for (unsigned int iwgt=1; iwgt<10; ++iwgt) {
	  //normalize to
	  double wgtval = genInfo->weights()[iwgt]*genWeight/genInfo->weights()[1];
	  scaleWeights->push_back(wgtval);
	}
      }

      //fill pdf variation weights
      if (firstPdfWeight>=0 && lastPdfWeight>=0 && lastPdfWeight<int(genInfo->weights().size()) && (lastPdfWeight-firstPdfWeight+1)==100) {

	//fill pdf variation weights after converting with mc2hessian transformation
	std::array<double, 100> inpdfweights;
	for (int iwgt=firstPdfWeight, ipdf=0; iwgt<=lastPdfWeight; ++iwgt, ++ipdf) {
	  inpdfweights[ipdf] = genInfo->weights()[iwgt]/genInfo->weights()[firstPdfWeight-1];
	}

	std::array<double, 60> outpdfweights;
	pdfweightshelper.DoMC2Hessian(inpdfweights.data(),outpdfweights.data());

	for (unsigned int iwgt=0; iwgt<60; ++iwgt) {
	  double wgtval = outpdfweights[iwgt]*genWeight;
	  pdfWeights->push_back(wgtval);
	}

	//fill alpha_s variation weights
	if (firstAlphasWeight>=0 && lastAlphasWeight>=0 && lastAlphasWeight<int(genInfo->weights().size())) {
	  for (int iwgt = firstAlphasWeight; iwgt<=lastAlphasWeight; ++iwgt) {
	    double wgtval = genInfo->weights()[iwgt]*genWeight/nomlheweight;
	    alphasWeights->push_back(wgtval);
	  }
	}

      }
    } else {

      if (lheInfo.isValid() && lheInfo->weights().size()>0) {

	double nomlheweight = lheInfo->weights()[0].wgt;

	//fill scale variation weights
	if (lheInfo->weights().size()>=9) {
	  for (unsigned int iwgt=0; iwgt<9; ++iwgt) {
	    double wgtval = lheInfo->weights()[iwgt].wgt*genWeight/nomlheweight;
	    scaleWeights->push_back(wgtval);
	  }
	}

	//fill pdf variation weights
	if (firstPdfWeight>=0 && lastPdfWeight>=0 && lastPdfWeight<int(lheInfo->weights().size()) && (lastPdfWeight-firstPdfWeight+1)==100) {

	  //fill pdf variation weights after converting with mc2hessian transformation
	  std::array<double, 100> inpdfweights;
	  for (int iwgt=firstPdfWeight, ipdf=0; iwgt<=lastPdfWeight; ++iwgt, ++ipdf) {
	    inpdfweights[ipdf] = lheInfo->weights()[iwgt].wgt/nomlheweight;
	  }

	  std::array<double, 60> outpdfweights;
	  pdfweightshelper.DoMC2Hessian(inpdfweights.data(),outpdfweights.data());

	  for (unsigned int iwgt=0; iwgt<60; ++iwgt) {
	    double wgtval = outpdfweights[iwgt]*genWeight;
	    pdfWeights->push_back(wgtval);
	  }

	  //fill alpha_s variation weights
	  if (firstAlphasWeight>=0 && lastAlphasWeight>=0 && lastAlphasWeight<int(lheInfo->weights().size())) {
	    for (int iwgt = firstAlphasWeight; iwgt<=lastAlphasWeight; ++iwgt) {
	      double wgtval = lheInfo->weights()[iwgt].wgt*genWeight/nomlheweight;
	      alphasWeights->push_back(wgtval);
	    }
	  }
	}
      }
    }

    //fill sum of weights histograms
    sumWeights->Fill(0.,genWeight);

    for (unsigned int iwgt=0; iwgt<scaleWeights->size(); ++iwgt) {
      sumScaleWeights->Fill(double(iwgt),(*scaleWeights)[iwgt]);
    }
    for (unsigned int iwgt=0; iwgt<pdfWeights->size(); ++iwgt) {
      sumPdfWeights->Fill(double(iwgt),(*pdfWeights)[iwgt]);
    }
    for (unsigned int iwgt=0; iwgt<alphasWeights->size(); ++iwgt) {
      sumAlphasWeights->Fill(double(iwgt),(*alphasWeights)[iwgt]);
    }
*/
    return true;
};

bool jet_timing_studies::fillGenParticles(){
  std::vector<const reco::Candidate*> prunedV;//Allows easier comparison for mother finding
  //Fills selected gen particles
  const double pt_cut = 0.0;
  int llp_id;
  if (model == 0)
  {
    llp_id = 35; //fourjet model
  }
  else if (model == 1)
  {
    llp_id = 9000006; // glueball model

  }
  else{
    llp_id = 1000021; // gluino, where both particles have the same id

  }

  for(size_t i=0; i<genParticles->size();i++)
  {
    if( (abs((*genParticles)[i].pdgId()) >= 1 && abs((*genParticles)[i].pdgId()) <= 6 && ( (*genParticles)[i].status() < 30 ))
       || (abs((*genParticles)[i].pdgId()) >= 11 && abs((*genParticles)[i].pdgId()) <= 16)
       || (abs((*genParticles)[i].pdgId()) == 21 && (*genParticles)[i].status() < 30)
       || (abs((*genParticles)[i].pdgId()) >= 22 && abs((*genParticles)[i].pdgId()) <= 25 && ( (*genParticles)[i].status() < 30))
       || (abs((*genParticles)[i].pdgId()) >= 32 && abs((*genParticles)[i].pdgId()) <= 42)
       || (abs((*genParticles)[i].pdgId()) >= 1000001 && abs((*genParticles)[i].pdgId()) <= 1000039)
       || (abs((*genParticles)[i].pdgId()) == 9000006 || abs((*genParticles)[i].pdgId()) == 9000007))
       {
         if ((*genParticles)[i].pt()>pt_cut){
           prunedV.push_back(&(*genParticles)[i]);
         }
       }

  }

  //Total number of gen particles
  nGenParticle = prunedV.size();
  int llp_index = 0;
  //Look for mother particle and Fill gen variables
  for(unsigned int i = 0; i < prunedV.size(); i++)
  {
    gParticleId[i] = prunedV[i]->pdgId();
    gParticleStatus[i] = prunedV[i]->status();
    gParticleE[i] = prunedV[i]->energy();
    gParticlePt[i] = prunedV[i]->pt();
    gParticlePx[i] = prunedV[i]->px();
    gParticlePy[i] = prunedV[i]->py();
    gParticlePz[i] = prunedV[i]->pz();
    gParticleEta[i] = prunedV[i]->eta();
    gParticlePhi[i] = prunedV[i]->phi();
    gParticleMotherId[i] = 0;
    gParticleMotherIndex[i] = -1;

    //For Neutralinos we try to find the decay vertex locaton.
    //Algorithm: Find the first daughter particle that is not a neutralino,
    //and call that the daughter. get the creation vertex of that daughter.
    // if ( (gParticleId[i] == 1000022 && gParticleStatus[i] == 22) )
    // {
    //   const reco::Candidate *dau = 0;
    //   bool foundDaughter = false;
    //   bool noDaughter = false;
    //   const reco::Candidate *tmpParticle = prunedV[i];
    //
    //   while (!foundDaughter && !noDaughter)
    //   {
    //     if (tmpParticle->numberOfDaughters() > 0)
    //     {
    //       dau = tmpParticle->daughter(0);
    //       if (dau && dau->pdgId() != 1000022){
    //         foundDaughter = true;
    //       }
    //       else{
    //         tmpParticle = dau;
    //       }
    //     }
    //     else
    //     {
    //       noDaughter = true;
    //     }
    //   }
    //
    //   if (foundDaughter)
    //   {
    //     gParticleDecayVertexX[i] = dau->vx();
    //     gParticleDecayVertexY[i] = dau->vy();
    //     gParticleDecayVertexZ[i] = dau->vz();
    //   }
    // }


    if(prunedV[i]->numberOfMothers() > 0)
    {
      //find the ID of the first mother that has a different ID than the particle itself
      const reco::Candidate* firstMotherWithDifferentID = findFirstMotherWithDifferentID(prunedV[i]);
      if (firstMotherWithDifferentID)
      {
        gParticleMotherId[i] = firstMotherWithDifferentID->pdgId();
        gParticleDecayVertexX[i] = firstMotherWithDifferentID->vx();
        gParticleDecayVertexY[i] = firstMotherWithDifferentID->vy();
        gParticleDecayVertexZ[i] = firstMotherWithDifferentID->vz();
      }

      //find the mother and keep going up the mother chain if the ID's are the same
      const reco::Candidate* originalMotherWithSameID = findOriginalMotherWithSameID(prunedV[i]);
      for(unsigned int j = 0; j < prunedV.size(); j++)
      {
        if(prunedV[j] == originalMotherWithSameID)
        {
          gParticleMotherIndex[i] = j;
          break;
        }
      }
    }
    else
    {
      gParticleMotherIndex[i] = -1;
    }

    //***************************************
    //Find LLPs production and decay vertices
    //***************************************
    bool llp_condition = (gParticleId[i] == llp_id) && gParticleStatus[i] == 106;
    if (!(model == 2)){
      llp_condition = (abs(gParticleId[i]) == llp_id || abs(gParticleId[i]) == llp_id+1)&& gParticleStatus[i] == 22;
    }

    if ( llp_condition)
    {
      bool firstllp = (gParticleId[i] == llp_id);
      bool secondllp = (gParticleId[i] == llp_id+1) || (gParticleId[i] == -1*llp_id);
      if(model == 2){
        firstllp = (gParticleId[i] == llp_id && llp_index == 0);
        secondllp = (gParticleId[i] == llp_id && llp_index == 1);
      }

      if (firstllp)
      {
        gLLP_prod_vertex_x[0] = prunedV[i]->vx();
        gLLP_prod_vertex_y[0] = prunedV[i]->vy();
        gLLP_prod_vertex_z[0] = prunedV[i]->vz();
      }
      else if (secondllp)
      {
        gLLP_prod_vertex_x[1] = prunedV[i]->vx();
        gLLP_prod_vertex_y[1] = prunedV[i]->vy();
        gLLP_prod_vertex_z[1] = prunedV[i]->vz();
      }
      // std::cout << "llp is "<<i<<","<<firstllp<<","<<secondllp<< std::endl;
      const reco::Candidate *dau = 0;
      bool foundDaughter = false;
      bool noDaughter = false;
      const reco::Candidate *tmpParticle = prunedV[i];

      while (!foundDaughter && !noDaughter)
      {
        // std::cout << "llp has daughters, and pdgID "<<tmpParticle->numberOfDaughters()<<","<<tmpParticle->pdgId()<< ","<<tmpParticle->status()<< ","<<tmpParticle->daughter(0)->pdgId()<<std::endl;
        if (tmpParticle->numberOfDaughters() > 1){
          // std::cout << "llp has daughters, and pdgID(1) "<<tmpParticle->numberOfDaughters()<<","<<tmpParticle->pdgId()<< ","<<tmpParticle->status()<< ","<<tmpParticle->daughter(1)->pdgId()<<std::endl;

        }
        if (tmpParticle->numberOfDaughters() > 0)
        {
          dau = tmpParticle->daughter(0);
          if (dau && (dau->pdgId() != llp_id && dau->pdgId() != llp_id+1))
          {
            foundDaughter = true;
          } else
          {
            tmpParticle = dau;
          }
        }
        else
        {
          noDaughter = true;
        }
      }

      if (foundDaughter)
      {

        if (firstllp)
        {
          gLLP_decay_vertex_x[0] = dau->vx();
          gLLP_decay_vertex_y[0] = dau->vy();
          gLLP_decay_vertex_z[0] = dau->vz();
          gLLP_pt[0] = sqrt(gParticlePx[i]*gParticlePx[i]+gParticlePy[i]*gParticlePy[i]);
          gLLP_e[0] = gParticleE[i];
          gLLP_eta[0] = gParticleEta[i];
          gLLP_phi[0] = gParticlePhi[i];
          gLLP_beta[0] = sqrt(gParticlePx[i]*gParticlePx[i]+gParticlePy[i]*gParticlePy[i]+gParticlePz[i]*gParticlePz[i])/gParticleE[i];
          gLLP_travel_time[0] = sqrt(pow(gLLP_decay_vertex_x[0]-gLLP_prod_vertex_x[0],2)
                                  +pow(gLLP_decay_vertex_y[0]-gLLP_prod_vertex_y[0],2)
                                  +pow(gLLP_decay_vertex_z[0]-gLLP_prod_vertex_z[0],2))/(30. * gLLP_beta[0]);//1/30 is to convert cm to ns
          if (model==2)
          {
            gLLP_travel_time[0] = sqrt(pow(gLLP_decay_vertex_x[0]-genVertexX,2)
                    +pow(gLLP_decay_vertex_y[0]-genVertexY,2)
                    +pow(gLLP_decay_vertex_z[0]-genVertexZ,2))/(30. * gLLP_beta[0]);//1/30 is to convert cm to ns
          }
          double radius = sqrt( pow(gLLP_decay_vertex_x[0],2) + pow(gLLP_decay_vertex_y[0],2) );
          double z = gLLP_decay_vertex_z[0] ;
          double ecal_radius = 129.0;
          double hcal_radius = 179.0;
          double EB_z = 268.36447217; // 129*sinh(1.479)
          double EE_z = 298.5; //where Ecal Endcap starts in z direction
          double ETL_rmin = 30.54540032; //Eta = 3.0, Z = 306cm
          double ETL_rmax = 128.81130156; //Eta = 1.6, Z = 306cm 
          double ETL_z = 306.0;
          // std::cout << "first llp has "<<tmpParticle->numberOfDaughters()<< " daughters" << std::endl;
          // std::cout << "first llp is "<<i<<","<<gParticleId[i] <<","<< gParticleStatus[i] <<","<<tmpParticle->pdgId()<<","<<tmpParticle->status()<<std::endl;

          for (unsigned int id = 0; id < tmpParticle->numberOfDaughters(); id++ )
          {
          //std::cout << "====================" << std::endl;
          //std::cout << " -> "<< tmpParticle->daughter(id)->pdgId() << std::endl;
            // std::cout << "first llp daughter: "<<tmpParticle->daughter(id)->pdgId() << "," << tmpParticle->daughter(id)->status()<<std::endl;
            if( id > 1 ) break;
            TLorentzVector tmp;
            tmp.SetPxPyPzE(tmpParticle->daughter(id)->px(), tmpParticle->daughter(id)->py(), tmpParticle->daughter(id)->pz(), tmpParticle->daughter(id)->energy());
            if(tmp.Pt()<pt_cut) continue;
            gLLP_daughter_pt[id] = tmp.Pt();
            gLLP_daughter_pz[id] = tmp.Pz();
            gLLP_daughter_eta[id] = tmp.Eta();
            gLLP_daughter_phi[id] = tmp.Phi();
            gLLP_daughter_e[id]  = tmp.E();

            double gLLP_daughter_travel_time_hcal= (1./30.)*(hcal_radius-radius)/(tmp.Pt()/tmp.E());// - (1./30.) * ecal_radius * cosh(tmp.Eta());//1/30 is to convert cm to ns

            gLLP_daughter_travel_time[id] = (1./30.)*(ecal_radius-radius)/(tmp.Pt()/tmp.E());// - (1./30.) * ecal_radius * cosh(tmp.Eta());//1/30 is to convert cm to ns
          //Calculate dt from generation point to ECAL face
            double x_ecal = gLLP_decay_vertex_x[0] + 30. * (tmp.Px()/tmp.E())*gLLP_daughter_travel_time[id];
            double y_ecal = gLLP_decay_vertex_y[0] + 30. * (tmp.Py()/tmp.E())*gLLP_daughter_travel_time[id];
            double z_ecal = gLLP_decay_vertex_z[0] + 30. * (tmp.Pz()/tmp.E())*gLLP_daughter_travel_time[id];
            double x_hcal = gLLP_decay_vertex_x[0] + 30. * (tmp.Px()/tmp.E())*gLLP_daughter_travel_time_hcal;
            double y_hcal = gLLP_decay_vertex_y[0] + 30. * (tmp.Py()/tmp.E())*gLLP_daughter_travel_time_hcal;
            double z_hcal = gLLP_decay_vertex_z[0] + 30. * (tmp.Pz()/tmp.E())*gLLP_daughter_travel_time_hcal;
    	    //std::cout << " gLLP_daughter_travel_time = " << gLLP_daughter_travel_time[id] << "x, y, z" << x_ecal << y_ecal << z_ecal <<std::endl;

	    if(tmp.Eta()>=0)
            {
            gLLP_daughter_travel_time_ETL[id] = (1./30.)*fabs(ETL_z-z)/fabs(tmp.Pz()/tmp.E());// - (1./30.) * ecal_radius * cosh(tmp.Eta());//1/30 is to convert cm to ns
            }
	    else
            {
            gLLP_daughter_travel_time_ETL[id] = (1./30.)*fabs(ETL_z+z)/fabs(tmp.Pz()/tmp.E());// - (1./30.) * ecal_radius * cosh(tmp.Eta());//1/30 is to convert cm to ns
            }
          //Calculate dt from generation point to ETL face
            double x_etl = gLLP_decay_vertex_x[0] + 30. * (tmp.Px()/tmp.E())*gLLP_daughter_travel_time_ETL[id];
            double y_etl = gLLP_decay_vertex_y[0] + 30. * (tmp.Py()/tmp.E())*gLLP_daughter_travel_time_ETL[id];
            double z_etl = gLLP_decay_vertex_z[0] + 30. * (tmp.Pz()/tmp.E())*gLLP_daughter_travel_time_ETL[id];
            double r_etl = sqrt( pow(x_etl,2) + pow(y_etl,2) );

      //if( fabs(z_ecal) < 10 && radius <= 1)
            if( fabs(z_ecal) < EB_z && radius <= ecal_radius && fabs(z) < EE_z)
    	    {
	      gLLP_daughter_EB[id] = true;
    	      photon_travel_time[id] = (1./30) * sqrt(pow(ecal_radius,2)+pow(z_ecal,2));
              photon_travel_time_pv[id] = (1./30) * sqrt(pow(x_ecal-genVertexX,2) + pow(y_ecal-genVertexY,2) + pow(z_ecal-genVertexZ,2));
              gen_time_pv[id] =  gLLP_travel_time[0] + gLLP_daughter_travel_time[id] - photon_travel_time_pv[id] + genVertexT;
              gen_time[id] = gLLP_travel_time[0] + gLLP_daughter_travel_time[id] - photon_travel_time[id] + genVertexT;

            }
            else if( r_etl > ETL_rmin && r_etl < ETL_rmax && radius <= ecal_radius && fabs(z) < ETL_z )
    	    {
    	      //std::cout <<" z_etl = "<< z_etl  <<std::endl;
    	      //std::cout << " gLLP_daughter_travel time ETL = " << gLLP_daughter_travel_time_ETL[id] << " , z = " << z << " , Pz = " << tmp.Pz() <<std::endl;
	      gLLP_daughter_ETL[id] = true;
    	      //std::cout << " gLLP_daughter_ETL = " << gLLP_daughter_ETL[id] << "travel time" << gLLP_daughter_travel_time_ETL[id] << "z" << z << "Pz" << tmp.Pz() <<std::endl;
    	      photon_travel_time_ETL[id] = (1./30) * sqrt(pow(r_etl,2)+pow(ETL_z,2));
    	      //std::cout << " gLLP_daughter_ETL = " << gLLP_daughter_ETL[id] << "photon travel time" << photon_travel_time_ETL[id] << " ,  ETL_z =  " << ETL_z << ", r_etl = " << r_etl <<std::endl;
              gen_time_ETL[id] = gLLP_travel_time[0] + gLLP_daughter_travel_time_ETL[id] - photon_travel_time_ETL[id] + genVertexT;

            }
            else
            {
              gLLP_daughter_travel_time[id] = -666;
              gen_time_pv[id] = -666.;
              gen_time[id] = -666.;
              photon_travel_time[id] = -666.;
              photon_travel_time_pv[id] = -666.;
            }
            double min_delta_r_calo = 666.;
    	      double min_delta_r = 666.;
    	      double min_delta_r_nocorr = 666.;
            double min_delta_r_hcal = 666.;
    	      unsigned int match_jet_index = 666;
            unsigned int match_jet_index_hcal = 666;
            unsigned int match_calojet_index = 666;
    	      double genJet_min_delta_r = 666.;
    	      unsigned int match_genJet_index = 666;

        // Correction of eta and phi based on ecal points
    	      double phi = atan((y_ecal-genVertexY)/(x_ecal-genVertexX));
            if  (x_ecal < 0.0){
              phi = TMath::Pi() + phi;
    	      }
    	      phi = deltaPhi(phi,0.0);
    	      double theta = atan(sqrt(pow(x_ecal-genVertexX,2)+pow(y_ecal-genVertexY,2))/abs(z_ecal-genVertexZ));
            double eta = -1.0*TMath::Sign(1.0, z_ecal-genVertexZ)*log(tan(theta/2));
    	      gLLP_daughter_eta_ecalcorr[id] = eta;
            gLLP_daughter_phi_ecalcorr[id] = phi;

        // Correction of eta and phi based on hcal points
            phi = atan((y_hcal-genVertexY)/(x_hcal-genVertexX));
            if  (x_hcal < 0.0){
              phi = TMath::Pi()  + phi;
            }
            phi = deltaPhi(phi,0.0);
            theta = atan(sqrt(pow(x_hcal-genVertexX,2)+pow(y_hcal-genVertexY,2))/abs(z_hcal-genVertexZ));
            eta = -1.0*TMath::Sign(1.0, z_hcal-genVertexZ)*log(tan(theta/2));
            gLLP_daughter_eta_hcalcorr[id] = eta;
            gLLP_daughter_phi_hcalcorr[id] = phi;


    	      for ( int i_jet = 0; i_jet < nGenJets; i_jet++)
    	      {
    		       double genJet_current_delta_r = deltaR(gLLP_daughter_eta[id], gLLP_daughter_phi[id],  genJetEta[i_jet], genJetPhi[i_jet]);
    	        //std::cout << i_jet << " current dR = " << genJet_current_delta_r << eta<<phi<<theta<<tan(theta/2.0)<<log(tan(theta/2.0))<<std::endl;
            	if ( genJet_current_delta_r < genJet_min_delta_r )
            	{
                genJet_min_delta_r = genJet_current_delta_r;
            		match_genJet_index = i_jet;
            		  //std::cout << i_jet << " min dR = " << genJet_min_delta_r << std::endl;
            	}
    	      }//end matching to genJets
            for ( int i_jet = 0; i_jet < nCaloJets; i_jet++ )
    	      {
  		        double current_delta_r = deltaR(gLLP_daughter_eta_ecalcorr[id], gLLP_daughter_phi_ecalcorr[id], calojetEta[i_jet], calojetPhi[i_jet]);
  	          if ( current_delta_r < min_delta_r_calo )
    	        {
    	  	      // min_delta_r_nocorr = deltaR(gLLP_daughter_eta[id], gLLP_daughter_phi[id], jetEta[i_jet], jetPhi[i_jet]);
    		        min_delta_r_calo = current_delta_r;
    		        match_calojet_index = i_jet;
    		        // std::cout << i_jet << " min dR = " << min_delta_r_calo << std::endl;
    	        }
              // std::cout << i_jet << " min dR = " << min_delta_r_calo << std::endl;
    	      }//end matching to calojets using ECAL radius
    	      for ( int i_jet = 0; i_jet < nJets; i_jet++ )
    	      {
  		        double current_delta_r = deltaR(gLLP_daughter_eta_ecalcorr[id], gLLP_daughter_phi_ecalcorr[id], jetEta[i_jet], jetPhi[i_jet]);
  	          if ( current_delta_r < min_delta_r )
    	        {
    	  	      min_delta_r_nocorr = deltaR(gLLP_daughter_eta[id], gLLP_daughter_phi[id], jetEta[i_jet], jetPhi[i_jet]);
    		        min_delta_r = current_delta_r;
    		        match_jet_index = i_jet;
    		  //std::cout << i_jet << " min dR = " << min_delta_r << std::endl;
    	        }
    	      }//end matching to jets using ECAL radius
            for ( int i_jet = 0; i_jet < nJets; i_jet++ )
    	      {
  		        double current_delta_r = deltaR(gLLP_daughter_eta_hcalcorr[id], gLLP_daughter_phi_hcalcorr[id], jetEta[i_jet], jetPhi[i_jet]);
  	          if ( current_delta_r < min_delta_r_hcal )
    	        {
    		        min_delta_r_hcal = current_delta_r;
    		        match_jet_index_hcal = i_jet;
    	        }
    	      }//end matching to jets using HCAL radius
    	      if( fabs(z_ecal) < EB_z && radius <= ecal_radius && z < EE_z)
            {
              if (min_delta_r_calo < 0.3)
              {
                gLLP_daughter_match_calojet_index[id] = match_calojet_index;
                gLLP_min_delta_r_match_calojet[id] = min_delta_r_calo;
                // gLLP_min_delta_r_nocorr_match_jet[id] = min_delta_r_nocorr;

              }
              if ( min_delta_r < 0.3 )
      	      {
      	        gLLP_daughter_match_jet_index[id] = match_jet_index;
      	        gLLP_min_delta_r_match_jet[id] = min_delta_r;
      	        gLLP_min_delta_r_nocorr_match_jet[id] = min_delta_r_nocorr;
      	      }
              if ( min_delta_r < 0.45 )
              {
                gLLP_daughter_match_jet_index_loose[id] = match_jet_index;
                gLLP_min_delta_r_match_jet_loose[id] = min_delta_r;
              }
            }
            if( fabs(z_hcal) < 388.0 && radius <= hcal_radius)
            {
              if ( min_delta_r_hcal < 0.3 )
      	      {
      	        gLLP_daughter_match_jet_index_hcal[id] = match_jet_index_hcal;
      	        gLLP_min_delta_r_match_jet_hcal[id] = min_delta_r_hcal;
      	      }
              if ( min_delta_r_hcal < 0.45 )
      	      {
      	        gLLP_daughter_match_jet_index_hcal_loose[id] = match_jet_index_hcal;
      	        gLLP_min_delta_r_match_jet_hcal_loose[id] = min_delta_r_hcal;
      	      }
            }
    	      if ( genJet_min_delta_r < 0.3 )
    	      {
    	        gLLP_daughter_match_genJet_index[id] = match_genJet_index;
    	        gLLP_min_delta_r_match_genJet[id] = genJet_min_delta_r;
    	        //std::cout << "min dR = " << min_delta_r << " matched to jet index " << match_jet_index << std::endl;
    	      }

          }
          llp_index ++;
        }
    	  else if (secondllp)
    	  {
    	    gLLP_decay_vertex_x[1] = dau->vx();
    	    gLLP_decay_vertex_y[1] = dau->vy();
    	    gLLP_decay_vertex_z[1] = dau->vz();
          gLLP_pt[1] = sqrt(gParticlePx[i]*gParticlePx[i]+gParticlePy[i]*gParticlePy[i]);
          gLLP_e[1] = gParticleE[i];
          gLLP_eta[1] = gParticleEta[i];
          gLLP_phi[1] = gParticlePhi[i];
    	    gLLP_beta[1] = sqrt(gParticlePx[i]*gParticlePx[i]+gParticlePy[i]*gParticlePy[i]+gParticlePz[i]*gParticlePz[i])/gParticleE[i];
    	    gLLP_travel_time[1] = sqrt(pow(gLLP_decay_vertex_x[1]-gLLP_prod_vertex_x[1],2)
    				      +pow(gLLP_decay_vertex_y[1]-gLLP_prod_vertex_y[1],2)
    				      +pow(gLLP_decay_vertex_z[1]-gLLP_prod_vertex_z[1],2))/(30. * gLLP_beta[1]);//1/30 is to convert cm to ns
          if (model==2)
          {
            gLLP_travel_time[1] = sqrt(pow(gLLP_decay_vertex_x[1]-genVertexX,2)
      				      +pow(gLLP_decay_vertex_y[1]-genVertexY,2)
      				      +pow(gLLP_decay_vertex_z[1]-genVertexZ,2))/(30. * gLLP_beta[1]);//1/30 is to convert cm to ns
          }
    	  double radius = sqrt( pow(gLLP_decay_vertex_x[1],2) + pow(gLLP_decay_vertex_y[1],2) );
          double z = gLLP_decay_vertex_z[1] ;
    	  double ecal_radius = 129.0;
          double hcal_radius = 179.0;
          double EB_z = 268.36447217; // 129*sinh(1.479)
          double EE_z = 298.5;
          double ETL_rmin = 30.54540032; //Eta = 3.0, Z = 306cm
          double ETL_rmax = 128.81130156; //Eta = 1.6, Z = 306cm 
          double ETL_z = 306.0;
    	    /*
    	    Second two LLP daughters belong to LLP->pdgID()=36
          */
          // std::cout << "second llp has "<<tmpParticle->numberOfDaughters()<< " daughters" << std::endl;
          // std::cout << "second llp is "<<i<<","<<gParticleId[i] <<","<< gParticleStatus[i] <<","<<tmpParticle->pdgId()<<","<<tmpParticle->status()<<std::endl;
    	    for (unsigned int id = 0; id < tmpParticle->numberOfDaughters(); id++ )
    	    {
    	      //std::cout << " -> "<< tmpParticle->daughter(id)->pdgId() << std::endl;
            // std::cout << "second llp daughter: "<<tmpParticle->daughter(id)->pdgId() << "," << tmpParticle->daughter(id)->status()<<std::endl;

    	      if( id > 1 ) break;
    	      TLorentzVector tmp;
    	      tmp.SetPxPyPzE(tmpParticle->daughter(id)->px(), tmpParticle->daughter(id)->py(), tmpParticle->daughter(id)->pz(), tmpParticle->daughter(id)->energy());
              if(tmp.Pt()<pt_cut) continue;
              gLLP_daughter_pt[id+2] = tmp.Pt();
              gLLP_daughter_pz[id+2] = tmp.Pz();
    	      gLLP_daughter_eta[id+2] = tmp.Eta();
    	      gLLP_daughter_phi[id+2] = tmp.Phi();
    	      gLLP_daughter_e[id+2]  = tmp.E();
    	      //gLLP_daughter_travel_time[id+2] = (1./30.)*(ecal_radius-radius)/(tmp.Pt()/tmp.E()) - (1./30.) * ecal_radius * cosh(tmp.Eta());//1/30 is to convert cm to ns
              double gLLP_daughter_travel_time_hcal = (1./30.)*(hcal_radius-radius)/(tmp.Pt()/tmp.E());// - (1./30.) * ecal_radius * cosh(tmp.Eta());//1/30 is to convert cm to ns
              gLLP_daughter_travel_time[id+2] = (1./30.)*(ecal_radius-radius)/(tmp.Pt()/tmp.E());// - (1./30.) * ecal_radius * cosh(tmp.Eta());//1/30 is to convert cm to ns

    	      //Calculate dt from generation point to ECAL face
    	      double x_ecal = gLLP_decay_vertex_x[1] + 30. * (tmp.Px()/tmp.E())*gLLP_daughter_travel_time[id+2];
    	      double y_ecal = gLLP_decay_vertex_y[1] + 30. * (tmp.Py()/tmp.E())*gLLP_daughter_travel_time[id+2];
    	      double z_ecal = gLLP_decay_vertex_z[1] + 30. * (tmp.Pz()/tmp.E())*gLLP_daughter_travel_time[id+2];

              double x_hcal = gLLP_decay_vertex_x[1] + 30. * (tmp.Px()/tmp.E())*gLLP_daughter_travel_time_hcal;
    	      double y_hcal = gLLP_decay_vertex_y[1] + 30. * (tmp.Py()/tmp.E())*gLLP_daughter_travel_time_hcal;
    	      double z_hcal = gLLP_decay_vertex_z[1] + 30. * (tmp.Pz()/tmp.E())*gLLP_daughter_travel_time_hcal;

	      if(tmp.Eta()>=0)
              {
              gLLP_daughter_travel_time_ETL[id+2] = (1./30.)*fabs(ETL_z-z)/fabs(tmp.Pz()/tmp.E());// - (1./30.) * ecal_radius * cosh(tmp.Eta());//1/30 is to convert cm to ns
              }
	      else
              {
              gLLP_daughter_travel_time_ETL[id+2] = (1./30.)*fabs(ETL_z+z)/fabs(tmp.Pz()/tmp.E());// - (1./30.) * ecal_radius * cosh(tmp.Eta());//1/30 is to convert cm to ns
              }
              //Calculate dt from generation point to ETL face
              double x_etl = gLLP_decay_vertex_x[1] + 30. * (tmp.Px()/tmp.E())*gLLP_daughter_travel_time_ETL[id+2];
              double y_etl = gLLP_decay_vertex_y[1] + 30. * (tmp.Py()/tmp.E())*gLLP_daughter_travel_time_ETL[id+2];
              double z_etl = gLLP_decay_vertex_z[1] + 30. * (tmp.Pz()/tmp.E())*gLLP_daughter_travel_time_ETL[id+2];
              double r_etl = sqrt( pow(x_etl,2) + pow(y_etl,2) );

    	      if( fabs(z_ecal) < EB_z && radius <= ecal_radius && fabs(z) < EE_z)
    	      // if( fabs(z_ecal) < 10 && radius <= 0.1)
    	      {
    	      gLLP_daughter_EB[id+2] = true;
              photon_travel_time[id+2] = (1./30) * sqrt(pow(ecal_radius,2)+pow(z_ecal,2));
              photon_travel_time_pv[id+2] = (1./30) * sqrt(pow(x_ecal-genVertexX,2) + pow(y_ecal-genVertexY,2) + pow(z_ecal-genVertexZ,2));
              gen_time_pv[id+2] =  gLLP_travel_time[1] + gLLP_daughter_travel_time[id+2] - photon_travel_time_pv[id+2] + genVertexT;
              gen_time[id+2] = gLLP_travel_time[1] + gLLP_daughter_travel_time[id+2] - photon_travel_time[id+2] + genVertexT;
    	      }
              else if( r_etl > ETL_rmin && r_etl < ETL_rmax && radius <= ecal_radius && fabs(z) < ETL_z )
    	      {
	        gLLP_daughter_ETL[id+2] = true;
    	        photon_travel_time_ETL[id+2] = (1./30) * sqrt(pow(r_etl,2)+pow(ETL_z,2));
                gen_time_ETL[id+2] = gLLP_travel_time[1] + gLLP_daughter_travel_time_ETL[id+2] - photon_travel_time_ETL[id+2] + genVertexT;

              }
    	      else
    	      {
    	      gLLP_daughter_travel_time[id+2] = -666;
              gen_time_pv[id+2] = -666.;
              gen_time[id+2] = -666.;
              photon_travel_time[id+2] = -666.;
              photon_travel_time_pv[id+2] = -666.;
    	      }
    	      double genJet_min_delta_r = 666.;
            unsigned int match_genJet_index = 666;
            double min_delta_r_calo = 666.;
    	      double min_delta_r = 666.;
            double min_delta_r_hcal = 666.;
    	      double min_delta_r_nocorr = 666.;
            unsigned int match_calojet_index = 666;
    	      unsigned int match_jet_index = 666;
            unsigned int match_jet_index_hcal = 666;

            //Corrections for angles based on ECAL radius
    	      double phi = atan((y_ecal-genVertexY)/(x_ecal-genVertexX));
            if  (x_ecal < 0.0){
              phi = TMath::Pi() + phi;
            }
    	      phi = deltaPhi(phi,0.0);
    	      double theta = atan(sqrt(pow(x_ecal-genVertexX,2)+pow(y_ecal-genVertexY,2))/abs(z_ecal-genVertexZ));
    	      double eta = -1.0*TMath::Sign(1.0,z_ecal-genVertexZ)*log(tan(theta/2));
    	      gLLP_daughter_eta_ecalcorr[id+2] = eta;
    	      gLLP_daughter_phi_ecalcorr[id+2] = phi;

            //Corrections for angles based on HCAL radius
            phi = atan((y_hcal-genVertexY)/(x_hcal-genVertexX));
            if  (x_hcal < 0.0){
              phi = TMath::Pi() + phi;
            }
            phi = deltaPhi(phi,0.0);
            theta = atan(sqrt(pow(x_hcal-genVertexX,2)+pow(y_hcal-genVertexY,2))/abs(z_hcal-genVertexZ));
            eta = -1.0*TMath::Sign(1.0,z_ecal-genVertexZ)*log(tan(theta/2));
            gLLP_daughter_eta_hcalcorr[id+2] = eta;
            gLLP_daughter_phi_hcalcorr[id+2] = phi;
    	      for ( int i_jet = 0; i_jet < nGenJets; i_jet++ )
    	      {
      		    double genJet_current_delta_r = deltaR(gLLP_daughter_eta[id+2], gLLP_daughter_phi[id+2],genJetEta[i_jet], genJetPhi[i_jet]);
      	      if ( genJet_current_delta_r < genJet_min_delta_r )
      	      {
    	  	      genJet_min_delta_r = genJet_current_delta_r;
    		        match_genJet_index = i_jet;
      		  //std::cout << i_jet << " min dR = " << min_delta_r << std::endl;
      	      }
    	      }//end matching to genJets
            for ( int i_jet = 0; i_jet < nCaloJets; i_jet++ )
    	      {
  		        double current_delta_r = deltaR(gLLP_daughter_eta_ecalcorr[id+2], gLLP_daughter_phi_ecalcorr[id+2], calojetEta[i_jet], calojetPhi[i_jet]);
  	          if ( current_delta_r < min_delta_r_calo )
    	        {
    	  	      // min_delta_r_nocorr = deltaR(gLLP_daughter_eta[id], gLLP_daughter_phi[id], jetEta[i_jet], jetPhi[i_jet]);
    		        min_delta_r_calo = current_delta_r;
    		        match_calojet_index = i_jet;
              }
    		  //std::cout << i_jet << " min dR = " << min_delta_r << std::endl;
            }// end matching to calojets
    	      for ( int i_jet = 0; i_jet < nJets; i_jet++ )
    	      {
              double current_delta_r = deltaR(gLLP_daughter_eta_ecalcorr[id+2], gLLP_daughter_phi_ecalcorr[id+2] , jetEta[i_jet], jetPhi[i_jet]);
          		if ( current_delta_r < min_delta_r )
          		{
          		  min_delta_r_nocorr = deltaR(gLLP_daughter_eta[id+2], gLLP_daughter_phi[id+2], jetEta[i_jet], jetPhi[i_jet]);
          		  min_delta_r = current_delta_r;
          		  match_jet_index = i_jet;
          		}
    	      }//end matching to jets ecal
            for ( int i_jet = 0; i_jet < nJets; i_jet++ )
    	      {
              double current_delta_r = deltaR(gLLP_daughter_eta_hcalcorr[id+2], gLLP_daughter_phi_hcalcorr[id+2], jetEta[i_jet], jetPhi[i_jet]);
          		if ( current_delta_r < min_delta_r_hcal )
          		{
          		  min_delta_r_hcal = current_delta_r;
          		  match_jet_index_hcal = i_jet;
          		}
    	      }//end matching to jets hcal
            if( fabs(z_ecal) < EB_z && radius <= ecal_radius && z < EE_z)
            {
              if ( min_delta_r_calo < 0.3 )
              {
                gLLP_daughter_match_calojet_index[id+2] = match_calojet_index;
                gLLP_min_delta_r_match_calojet[id+2] = min_delta_r_calo;
              }
              if ( min_delta_r < 0.3 )
              {
                gLLP_daughter_match_jet_index[id+2] = match_jet_index;
                gLLP_min_delta_r_match_jet[id+2] = min_delta_r;
                gLLP_min_delta_r_nocorr_match_jet[id] = min_delta_r_nocorr;
              }
              if ( min_delta_r < 0.45 )
              {
                gLLP_daughter_match_jet_index_loose[id+2] = match_jet_index;
                gLLP_min_delta_r_match_jet_loose[id+2] = min_delta_r;
              }


            }
            if( fabs(z_hcal) < 388.0 && radius <= hcal_radius)
            {
              if ( min_delta_r_hcal < 0.3 )
              {
                gLLP_daughter_match_jet_index_hcal[id+2] = match_jet_index_hcal;
                gLLP_min_delta_r_match_jet_hcal[id+2] = min_delta_r_hcal;
              }
              if ( min_delta_r_hcal < 0.3 )
              {
                gLLP_daughter_match_jet_index_hcal_loose[id+2] = match_jet_index_hcal;
                gLLP_min_delta_r_match_jet_hcal_loose[id+2] = min_delta_r_hcal;
              }


            }


    	      if ( genJet_min_delta_r < 0.3 )
    	      {
    	        gLLP_daughter_match_genJet_index[id+2] = match_genJet_index;
    	        gLLP_min_delta_r_match_genJet[id+2] = genJet_min_delta_r;
    	        //std::cout << "min dR = " << min_delta_r << " matched to jet index " << match_jet_index << std::endl;
    	      }
    	    }//for daughters loop
    	  }//if particle ID = 36
    	}//if found daughters
    }
    //******************************
    //QCD Matching
    //******************************
    if (isQCD_) {
    if ( (abs(gParticleId[i])  <= 6 || abs(gParticleId[i]) == 21) && gParticleStatus[i] == 23)
	    {
		    const reco::Candidate *tmpParticle = prunedV[i];
		    TLorentzVector tmp;
		    tmp.SetPxPyPzE(tmpParticle->px(), tmpParticle->py(), tmpParticle->pz(), tmpParticle->energy());
		    genQCD_pt[nGenQCDParticles] = tmp.Pt();
		    genQCD_eta[nGenQCDParticles] = tmp.Eta();
		    genQCD_phi[nGenQCDParticles] = tmp.Phi();
		    genQCD_e[nGenQCDParticles]  = tmp.E();

		    genQCD_prod_vertex_x[nGenQCDParticles]  = gParticleProdVertexX[i];
		    genQCD_prod_vertex_y[nGenQCDParticles]  = gParticleProdVertexY[i];
		    genQCD_prod_vertex_z[nGenQCDParticles]  = gParticleProdVertexZ[i];

		    genQCD_decay_vertex_x[nGenQCDParticles]  = gParticleDecayVertexX[i];
		    genQCD_decay_vertex_y[nGenQCDParticles]  = gParticleDecayVertexY[i];
		    genQCD_decay_vertex_z[nGenQCDParticles]  = gParticleDecayVertexZ[i];

		    double r0 = sqrt( pow(gParticleProdVertexX[i],2) + pow(gParticleProdVertexY[i],2) );
		    double z0 = gParticleProdVertexZ[i];

		    double ecal_radius = 129.0;
          	    double hcal_radius = 179.0;
          	    double EB_z = 268.36447217; // 129*sinh(1.479)
         	    double EE_z = 298.5; //where Ecal Endcap starts in z direction
          	    double ETL_rmin = 30.54540032; //Eta = 3.0, Z = 306cm
          	    double ETL_rmax = 128.81130156; //Eta = 1.6, Z = 306cm 
          	    double ETL_z = 306.0;

		    double t0_ecal = (1./30.)*(ecal_radius-r0)/(tmp.Pt()/tmp.E()); 

		    genQCD_travel_time[nGenQCDParticles]  = t0_ecal;

		    double x0_ecal = gParticleProdVertexX[i] + 30. * (tmp.Px()/tmp.E())*t0_ecal;
		    double y0_ecal = gParticleProdVertexY[i] + 30. * (tmp.Py()/tmp.E())*t0_ecal;
		    double z0_ecal = gParticleProdVertexZ[i] + 30. * (tmp.Pz()/tmp.E())*t0_ecal;
	
		    double t0_etl = 0.;

		    if(tmp.Eta()>=0)
    		    {
		    t0_etl = (1./30.)*fabs(ETL_z-z0)/fabs(tmp.Pz()/tmp.E()); 
		    }
		    else
    		    {
		    t0_etl = (1./30.)*fabs(ETL_z+z0)/fabs(tmp.Pz()/tmp.E()); 
		    }

		    genQCD_travel_time_ETL[nGenQCDParticles]  = t0_etl;

		    double x0_etl = gParticleProdVertexX[i] + 30. * (tmp.Px()/tmp.E())*t0_etl;
		    double y0_etl = gParticleProdVertexY[i] + 30. * (tmp.Py()/tmp.E())*t0_etl;
		    double z0_etl = gParticleProdVertexZ[i] + 30. * (tmp.Pz()/tmp.E())*t0_etl;
		    double r0_etl = sqrt( pow(x0_etl,2) + pow(y0_etl,2) );

/*
 * 		    double hcal_radius = 179.0;
 * 		    		    double t0_hcal = (1./30.)*(hcal_radius-r0)/(tmp.Pt()/tmp.E()); 
 		    double x0_hcal = gParticleProdVertexX[i] + 30. * (tmp.Px()/tmp.E())*t0_hcal;
		    double y0_hcal = gParticleProdVertexY[i] + 30. * (tmp.Py()/tmp.E())*t0_hcal;
		    double z0_hcal = gParticleProdVertexZ[i] + 30. * (tmp.Pz()/tmp.E())*t0_hcal;
*/
		    if( fabs(z0_ecal) < EB_z && r0 <= ecal_radius && fabs(z0) < EE_z)
		    {
		        genQCD_EB[nGenQCDParticles]  = true;
			double tg_ecal = (1./30) * sqrt(pow(ecal_radius,2)+pow(z0_ecal,2));
		    	genQCD_photon_travel_time[nGenQCDParticles]  = tg_ecal;
		    	genQCD_time[nGenQCDParticles]  = genQCD_travel_time[nGenQCDParticles] - genQCD_photon_travel_time[nGenQCDParticles] + genVertexT;
		    }
		    else if( r0_etl > ETL_rmin && r0_etl < ETL_rmax && r0 <= ecal_radius &&  fabs(z0) < ETL_z )
		    {
		        genQCD_ETL[nGenQCDParticles]  = true;
			double tg_etl = (1./30) * sqrt(pow(r0_etl,2)+pow(ETL_z,2));
		    	genQCD_photon_travel_time_ETL[nGenQCDParticles]  = tg_etl;
		    	genQCD_time_ETL[nGenQCDParticles]  = genQCD_travel_time_ETL[nGenQCDParticles] - genQCD_photon_travel_time_ETL[nGenQCDParticles] + genVertexT;
		    }


		    double min_delta_r = 666.;
		    unsigned int match_jet_index = 666;

		    for ( int i_jet = 0; i_jet < nJets; i_jet++ )
		    {
			    double current_delta_r = deltaR(genQCD_eta[nGenQCDParticles], genQCD_phi[nGenQCDParticles], jetEta[i_jet], jetPhi[i_jet]);
			    if ( current_delta_r < min_delta_r )
			    {
				    min_delta_r = current_delta_r;
				    match_jet_index = i_jet;
			    }
		    }//end matching to jets
		    if ( min_delta_r < 0.3 )
		    {
			    genParticleQCD_match_jet_index[nGenQCDParticles] = match_jet_index;
			    genParticleQCD_min_delta_r_match_jet[nGenQCDParticles] = min_delta_r;
		    }
		    nGenQCDParticles ++;
	    }// quarks or gluons with status 23
    }// end QCD matching part

  }// for loop of genParticles
  return true;
};



bool jet_timing_studies::fillTrigger(const edm::Event& iEvent)
{

  //fill trigger information
  const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
  // std::cout << "\n === TRIGGER PATHS === " << std::endl;
  //********************************************************************
  //Option to save all HLT path names in the ntuple per event
  //Expensive option in terms of ntuple size
  //********************************************************************
  nameHLT->clear();
  for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i)
  {
    string hltPathNameReq = "HLT_";
    //if (triggerBits->accept(i))
    if ((names.triggerName(i)).find(hltPathNameReq) != string::npos) nameHLT->push_back(names.triggerName(i));
    /*
    std::cout << "Trigger " << names.triggerName(i) <<
    ", prescale " << triggerPrescales->getPrescaleForIndex(i) <<
    ": " << (triggerBits->accept(i) ? "PASS" : "fail (or not run)")
    << std::endl;
    if ((names.triggerName(i)).find(hltPathNameReq) != string::npos && triggerBits->accept(i)) std::cout << "Trigger " << names.triggerName(i) <<
    ": " << (triggerBits->accept(i) ? "PASS" : "fail (or not run)")
    << std::endl;
    */
  }
  //std::cout << "n triggers: " <<  nameHLT->size() << std::endl;
  //std::cout << "====================" << std::endl;
  //for ( unsigned int i = 0; i < nameHLT->size(); i++ )
  //{
  //  std::cout << i << " -> " << nameHLT->at(i) << std::endl;
  //}
  //********************************************************************
  // Save trigger decisions in array of booleans
  //********************************************************************

  for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i)
  {
    string hltPathNameReq = "HLT_";
    if ((names.triggerName(i)).find(hltPathNameReq) == string::npos) continue;
    if ((names.triggerName(i)).find_last_of("_") == string::npos) continue;
    int lastUnderscorePos = (names.triggerName(i)).find_last_of("_");
    string hltPathNameWithoutVersionNumber = (names.triggerName(i)).substr(0,lastUnderscorePos);

    for (unsigned int j = 0; j < NTriggersMAX; ++j)
    {
      if (triggerPathNames[j] == "") continue;
      if (hltPathNameWithoutVersionNumber == triggerPathNames[j])
      {
        triggerDecision[j] = triggerBits->accept(i);
        //triggerHLTPrescale[j] = triggerPrescales->getPrescaleForIndex(i);
      }
    }
  }

  //********************************************************************
  // Print Trigger Objects
  //********************************************************************
/*
  for (pat::TriggerObjectStandAlone trigObject : *triggerObjects)
  {
    //cout << "triggerObj: " << trigObject.pt() << " " << trigObject.eta() << " " << trigObject.phi() << "\n";
    //bool foundRazor = false;
    //Need to unpack the filter labels before checking
    trigObject.unpackFilterLabels(iEvent, *triggerBits);
    for(int j=0; j<int(trigObject.filterLabels().size());j++)
    {
      //if ((trigObject.filterLabels())[j] == "hltRsqMR200Rsq0p0196MR100Calo") foundRazor = true;
      // trigObject.unpackPathNames(names);
      // cout << "filter: " << (trigObject.pathNames())[j] << " " << (trigObject.filterLabels())[j] << "\n";
      //cout << "filter: " << (trigObject.filterLabels())[j] << "\n";
    }
  }
*/
//define this as a plug-in
  return true;
};
DEFINE_FWK_MODULE(jet_timing_studies);
