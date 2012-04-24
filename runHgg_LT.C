// $Id: runHgg.C,v 1.6 2012/04/03 08:23:23 mingyang Exp $
#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TSystem.h>
#include <TProfile.h>
#include "MitAna/DataUtil/interface/Debug.h"
#include "MitAna/Catalog/interface/Catalog.h"
#include "MitAna/TreeMod/interface/Analysis.h"
#include "MitAna/TreeMod/interface/HLTMod.h"
#include "MitAna/PhysicsMod/interface/RunLumiSelectionMod.h"
#include "MitAna/PhysicsMod/interface/MCProcessSelectionMod.h"
#include "MitAna/PhysicsMod/interface/PublisherMod.h"
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitAna/DataTree/interface/PFJetCol.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitPhysics/Mods/interface/GoodPVFilterMod.h"
#include "MitPhysics/Mods/interface/MuonIDMod.h"
#include "MitPhysics/Mods/interface/ElectronIDMod.h"
#include "MitPhysics/Mods/interface/ElectronCleaningMod.h"
#include "MitPhysics/Mods/interface/PhotonIDMod.h"
#include "MitPhysics/Mods/interface/PhotonPairSelector.h"
#include "MitPhysics/Mods/interface/PhotonTreeWriter.h"
#include "MitPhysics/Mods/interface/PhotonCleaningMod.h"
#include "MitPhysics/Mods/interface/MergeLeptonsMod.h"
#include "MitPhysics/Mods/interface/JetCorrectionMod.h"
#include "MitPhysics/Mods/interface/PhotonMvaMod.h"
#include "MitPhysics/Mods/interface/MVASystematicsMod.h"
#endif

//--------------------------------------------------------------------------------------------------
void runHgg_LT(const char *fileset    = "0000",
            const char *skim       = "noskim",
	    //const char *dataset = "f11--h123gg-gf-v14b-pu",
            const char *dataset = "r11a-pho-j16-v1",   
            const char *book       = "local/filefi/025",
            const char *catalogDir = "/home/cmsprod/catalog",
            const char *outputName = "hgg",
            int         nEvents    = 2000)
{
  //------------------------------------------------------------------------------------------------
  // some parameters get passed through the environment
  //------------------------------------------------------------------------------------------------
  char json[1024], overlap[1024];
  float overlapCut = -1;

  if (gSystem->Getenv("MIT_PROD_JSON"))
    sprintf(json,   "%s",gSystem->Getenv("MIT_PROD_JSON"));
  else {
    sprintf(json, "%s", "~");
    //printf(" JSON file was not properly defined. EXIT!\n");
    //return;
  } 

  TString jsonFile = TString("/home/fabstoec/cms/json/") + TString(json);
  //TString jsonFile = TString("/home/fabstoec/cms/json/") + TString("Cert_136033-149442_7TeV_Dec22ReReco_Collisions10_JSON_v4.txt");
  Bool_t  isData   = ( (jsonFile.CompareTo("/home/fabstoec/cms/json/~") != 0) );
  
  if (gSystem->Getenv("MIT_PROD_OVERLAP")) {
    sprintf(overlap,"%s",gSystem->Getenv("MIT_PROD_OVERLAP"));
    if (EOF == sscanf(overlap,"%f",&overlapCut)) {
      printf(" Overlap was not properly defined. EXIT!\n");
      return;
    }
  }
  else {
     sprintf(overlap,"%s", "-1.0");
    //printf(" OVERLAP file was not properly defined. EXIT!\n");
    //return;
  } 

  printf("\n Initialization worked. \n\n");

  //------------------------------------------------------------------------------------------------
  // some global setups
  //------------------------------------------------------------------------------------------------
  using namespace mithep;
  gDebugMask  = Debug::kGeneral;
  gDebugLevel = 3;

  //------------------------------------------------------------------------------------------------
  // set up information
  //------------------------------------------------------------------------------------------------
  RunLumiSelectionMod *runLumiSel = new RunLumiSelectionMod;
  runLumiSel->SetAcceptMC(kTRUE);                          // Monte Carlo events are always accepted

  MCProcessSelectionMod *mcselmod = new MCProcessSelectionMod;
  
  MVASystematicsMod *sysMod = new MVASystematicsMod;
  sysMod->SetMCR9Scale(1.0035, 1.0035);  
  sysMod->SetIsData(isData);
  
  // only select on run- and lumisection numbers when valid json file present
  if ((jsonFile.CompareTo("/home/fabstoec/cms/json/~") != 0) &&
      (jsonFile.CompareTo("/home/fabstoec/cms/json/-") != 0)   ) {
    runLumiSel->AddJSONFile(jsonFile.Data());
  }
  if ((jsonFile.CompareTo("/home/fabstoec/cms/json/-") == 0)   ) {
    printf("\n WARNING -- Looking at data without JSON file: always accept.\n\n");
    runLumiSel->SetAbortIfNotAccepted(kFALSE);   // accept all events if there is no valid JSON file
  }

  printf("\n Run lumi worked. \n\n");

  //------------------------------------------------------------------------------------------------
  // HLT information
  //------------------------------------------------------------------------------------------------
  HLTMod *hltModP = new HLTMod("HLTModP");  
  hltModP->AddTrigger("HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v*",160000,161176);
  
  hltModP->AddTrigger("HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v*",161216,165633);
  hltModP->AddTrigger("HLT_Photon26_CaloIdL_IsoVL_Photon18_R9Id_v*",161216,165633);
  hltModP->AddTrigger("HLT_Photon26_R9Id_Photon18_CaloIdL_IsoVL_v*",161216,165633);
  hltModP->AddTrigger("HLT_Photon20_R9Id_Photon18_R9Id_v*",161216,165633);

  hltModP->AddTrigger("HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v*",165970,173198);
  hltModP->AddTrigger("HLT_Photon26_CaloIdL_IsoVL_Photon18_R9Id_v*",165970,173198);
  hltModP->AddTrigger("HLT_Photon26_R9Id_Photon18_CaloIdL_IsoVL_v*",165970,173198);
  hltModP->AddTrigger("HLT_Photon26_R9Id_Photon18_R9Id_v*",165970,173198);
  hltModP->AddTrigger("HLT_Photon36_CaloIdL_IsoVL_Photon22_CaloIdL_IsoVL_v*",165970,173198);
  hltModP->AddTrigger("HLT_Photon36_R9Id_Photon22_CaloIdL_IsoVL_v*",165970,173198);
  hltModP->AddTrigger("HLT_Photon36_R9Id_Photon22_R9Id_v*",165970,173198);
  
  hltModP->AddTrigger("HLT_Photon36_CaloId_IsoVL_Photon22_R9Id_v*",165970,166967);
  hltModP->AddTrigger("HLT_Photon36_CaloIdL_IsoVL_Photon22_R9Id_v*",167039,173198);
  
  hltModP->AddTrigger("HLT_Photon26_CaloIdXL_IsoXL_Photon18_CaloIdXL_IsoXL_v*",173236,178380);
  hltModP->AddTrigger("HLT_Photon26_CaloIdXL_IsoXL_Photon18_R9Id_v*",173236,178380);
  hltModP->AddTrigger("HLT_Photon26_R9Id_Photon18_CaloIdXL_IsoXL_v*",173236,178380);
  hltModP->AddTrigger("HLT_Photon26_R9Id_Photon18_R9Id_v*",173236,178380);
  hltModP->AddTrigger("HLT_Photon36_CaloIdL_IsoVL_Photon22_CaloIdL_IsoVL_v*",173236,178380);
  hltModP->AddTrigger("HLT_Photon36_CaloIdL_IsoVL_Photon22_R9Id_v*",173236,178380);
  hltModP->AddTrigger("HLT_Photon36_R9Id_Photon22_CaloIdL_IsoVL_v*",173236,178380);
  hltModP->AddTrigger("HLT_Photon36_R9Id_Photon22_R9Id_v*",173236,178380);
  
  hltModP->AddTrigger("HLT_Photon26_CaloIdXL_IsoXL_Photon18_CaloIdXL_IsoXL_Mass60_v*",178420,999999);
  hltModP->AddTrigger("HLT_Photon26_CaloIdXL_IsoXL_Photon18_R9IdT_Mass60_v*",178420,999999);
  hltModP->AddTrigger("HLT_Photon26_R9IdT_Photon18_CaloIdXL_IsoXL_Mass60_v*",178420,999999);
  hltModP->AddTrigger("HLT_Photon26_R9IdT_Photon18_R9IdT_Mass60_v*",178420,999999);
  hltModP->AddTrigger("HLT_Photon36_CaloIdL_IsoVL_Photon22_CaloIdL_IsoVL_v*",178420,999999);
  hltModP->AddTrigger("HLT_Photon36_CaloIdL_IsoVL_Photon22_R9Id_v*",178420,999999);
  hltModP->AddTrigger("HLT_Photon36_R9Id_Photon22_CaloIdL_IsoVL_v*",178420,999999);
  hltModP->AddTrigger("HLT_Photon36_R9Id_Photon22_R9Id_v*",178420,999999);
    
  hltModP->SetTrigObjsName("MyHltPhotObjs");
  hltModP->SetAbortIfNotAccepted(isData);
 
  //-----------------------------------
  // Lepton Selection
  //-----------------------------------
  ElectronIDMod* eleIdMod = new ElectronIDMod;
  eleIdMod -> SetEtMin(20.);
  eleIdMod -> SetEtaMax(2.5);
  eleIdMod -> SetApplyEcalFiducial(true);

  eleIdMod -> SetIDType("Hgg_LeptonTag_WP85Id");

  eleIdMod -> SetWhichVertex(0);
  eleIdMod -> SetD0Cut(0.02);
  eleIdMod -> SetDZCut(0.1);

  eleIdMod -> SetIsoType("CombinedRelativeConeAreaCorrected");
  //                                Barrel   Endcap
  eleIdMod -> SetCombRelativeIsoCut(0.053,   0.042  );
  eleIdMod -> SetOutputName("HggLeptonTagElectrons");


  MuonIDMod* muonIdMod = new MuonIDMod;
  // base kinematics
  muonIdMod -> SetPtMin(20.);
  muonIdMod -> SetEtaCut(2.4);
  // base ID
  muonIdMod -> SetIDType("WMuId");
  // distance to Vtx restriction
  muonIdMod -> SetWhichVertex(0); // this is a 'hack'.. but hopefully good enough...
  muonIdMod -> SetD0Cut(0.02);
  muonIdMod -> SetDZCut(0.1);
  // isolation
  muonIdMod -> SetIsoType("CombinedRelativeConeAreaCorrected");
  muonIdMod -> SetCombRelativeIsoCut(0.1);
  muonIdMod -> SetOutputName("HggLeptonTagMuons");

  //------------------------------------------------------------------------------------------------
  // select events with a good primary vertex
  //------------------------------------------------------------------------------------------------
  GoodPVFilterMod *goodPVFilterMod = new GoodPVFilterMod;
  goodPVFilterMod->SetMinVertexNTracks(0);
  goodPVFilterMod->SetMinNDof         (4.0);
  goodPVFilterMod->SetMaxAbsZ         (24.0);
  goodPVFilterMod->SetMaxRho          (2.0);
  goodPVFilterMod->SetAbortIfNotAccepted(kFALSE);
  goodPVFilterMod->SetIsMC(!isData);

  PublisherMod<PFJet,Jet> *pubJet = new PublisherMod<PFJet,Jet>("JetPub");
  pubJet->SetInputName("AKt5PFJets");
  pubJet->SetOutputName("PubAKt5PFJets");

  JetCorrectionMod *jetCorr = new JetCorrectionMod;
  jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/START42_V12_AK5PF_L1FastJet.txt")).Data())); 
  jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/START42_V12_AK5PF_L2Relative.txt")).Data())); 
  jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/START42_V12_AK5PF_L3Absolute.txt")).Data())); 
  if(isData){ 
    jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/START42_V12_AK5PF_L2L3Residual.txt")).Data())); 
  }
  jetCorr->SetInputName(pubJet->GetOutputName());
  jetCorr->SetCorrectedName("CorrectedJets");

  Bool_t excludedoubleprompt = kFALSE;
  if (TString(dataset).Contains("-pj")) {
    mcselmod->ExcludeProcess(18);
    mcselmod->ExcludeProcess(114);
    excludedoubleprompt = kTRUE;
  }

  if (TString(dataset).Contains("-qcd2em") || TString(dataset).Contains("-qcd-2em")) {
    excludedoubleprompt = kTRUE;
  }
  
  PhotonMvaMod *photreg = new PhotonMvaMod;
  photreg->SetOutputName("GoodPhotonsRegr");
  photreg->SetIsData(isData);

    
  PhotonPairSelector         *photcic = new PhotonPairSelector("PhotonPairSelectorCiC");
  photcic->SetOutputName("GoodPhotonsCIC");
  photcic->SetPhotonSelType("CiCSelection");
  photcic->SetVertexSelType("CiCMVASelection");
  photcic->DoMCSmear(kTRUE);
  photcic->DoDataEneCorr(kTRUE);
  photcic->SetPhotonsFromBranch(kFALSE);
  photcic->SetInputPhotonsName(photreg->GetOutputName());
  photcic->SetMCSmearFactors(0.0067,0.0077,0.0096,0.0141,0.0196,0.0268,0.0279,0.0293,0.0301);//ming:smear(sigE/E)
  photcic->AddEnCorrPerRun(160431,167913,0.9941,0.9941,1.0004,0.9916,1.0045,1.0033,1.0082,0.9958,1.0064);//ming:Emc/Edata
  photcic->AddEnCorrPerRun(170000,172619,0.9954,0.9954,1.0016,0.9937,1.0066,0.9976,1.0025,0.9940,1.0046);
  photcic->AddEnCorrPerRun(172620,173692,0.9955,0.9955,1.0017,0.9929,1.0058,0.9986,1.0035,0.9923,1.0029);
  photcic->AddEnCorrPerRun(175830,177139,0.9958,0.9958,1.0021,0.9944,1.0073,0.9968,1.0017,0.9933,1.004);
  photcic->AddEnCorrPerRun(177140,178421,0.9962,0.9962,1.0025,0.9946,1.0075,0.9960,1.0010,0.9944,1.005);
  photcic->AddEnCorrPerRun(178424,180252,0.9961,0.9961,1.0024,0.9942,1.0071,0.9921,0.9970,0.9953,1.0059); 
  photcic->SetDoMCR9Scaling(kTRUE);
  photcic->SetMCR9Scale(1.0048, 1.00492);
  photcic->SetDoMCErrScaling(kTRUE);
  photcic->SetMCErrScale(1.07, 1.045);    
  photcic->SetIsData(isData);

  PhotonPairSelector         *photpresel = new PhotonPairSelector("PhotonPairSelectorPresel");
  photpresel->SetOutputName("GoodPhotonsPresel");
  photpresel->SetPhotonSelType("MITSelection");
  photpresel->SetVertexSelType("CiCMVASelection");
  photpresel->DoMCSmear(kTRUE);
  photpresel->DoDataEneCorr(kTRUE);
  photpresel->SetPhotonsFromBranch(kFALSE);
  photpresel->SetInputPhotonsName(photreg->GetOutputName());
  //photpresel->SetMCSmearFactors(0.0045, 0.0084, 0.0109, 0.0156, 0.0203,0.0303,0.0326,0.0318,0.0331);
  photpresel->SetMCSmearFactors(0.0067,0.0077,0.0096,0.0141,0.0196,0.0268,0.0279,0.0293,0.0301);//ming:smear(sigE/E)
  //photpresel->AddEnCorrPerRun(160000,167913,0.9905,0.9905,0.9971,0.9976,1.0094,0.9994,1.0044,0.9968,1.0079);
  //photpresel->AddEnCorrPerRun(170249,172619,0.9909,0.9909,0.9975,0.9994,1.0112,0.9962,1.0012,0.9962,1.0072);
  //photpresel->AddEnCorrPerRun(172620,173692,0.9909,0.9909,0.9975,0.9977,1.0096,0.9963,1.0013,0.9947,1.0057);
  //photpresel->AddEnCorrPerRun(175860,177139,0.9911,0.9911,0.9977,0.9990,1.0109,0.9922,0.9973,0.9967,1.0077);
  //photpresel->AddEnCorrPerRun(177140,178421,0.9910,0.9910,0.9975,0.9987,1.0105,0.9921,0.9972,0.9975,1.0085);
  //photpresel->AddEnCorrPerRun(178424,999999,0.9903,0.9903,0.9969,0.9976,1.0095,0.9889,0.9940,0.9976,1.0086); 
  photpresel->AddEnCorrPerRun(160431,167913,0.9941,0.9941,1.0004,0.9916,1.0045,1.0033,1.0082,0.9958,1.0064);//ming:Emc/Edata
  photpresel->AddEnCorrPerRun(170000,172619,0.9954,0.9954,1.0016,0.9937,1.0066,0.9976,1.0025,0.9940,1.0046);
  photpresel->AddEnCorrPerRun(172620,173692,0.9955,0.9955,1.0017,0.9929,1.0058,0.9986,1.0035,0.9923,1.0029);
  photpresel->AddEnCorrPerRun(175830,177139,0.9958,0.9958,1.0021,0.9944,1.0073,0.9968,1.0017,0.9933,1.004);
  photpresel->AddEnCorrPerRun(177140,178421,0.9962,0.9962,1.0025,0.9946,1.0075,0.9960,1.0010,0.9944,1.005);
  photpresel->AddEnCorrPerRun(178424,180252,0.9961,0.9961,1.0024,0.9942,1.0071,0.9921,0.9970,0.9953,1.0059); 
  photpresel->SetDoMCR9Scaling(kTRUE);
  photpresel->SetMCR9Scale(1.0035, 1.0035);  
  photpresel->SetDoMCSigIEtaIEtaScaling(kTRUE);
  photpresel->SetDoMCWidthScaling(kTRUE);  
  photpresel->SetDoMCErrScaling(kTRUE);
  photpresel->SetMCErrScale(1.07, 1.045);    
  photpresel->SetIsData(isData);

  
  PhotonTreeWriter *phottreecic = new PhotonTreeWriter("PhotonTreeWriterCiC");
  phottreecic->SetPhotonsFromBranch(kFALSE);
  phottreecic->SetInputPhotonsName(photcic->GetOutputName());
  phottreecic->SetEnableJets(kTRUE);
  phottreecic->SetPFJetsFromBranch(kFALSE);
  phottreecic->SetPFJetName(jetCorr->GetOutputName());
  phottreecic->SetExcludeDoublePrompt(excludedoubleprompt);
  phottreecic->SetIsData(isData);
  
  phottreecic->SetApplyLeptonTag(true);
  phottreecic->SetLeptonTagElectronsName("HggLeptonTagElectrons");
  phottreecic->SetLeptonTagMuonsName    ("HggLeptonTagMuons");

  PhotonTreeWriter *phottreepresel = new PhotonTreeWriter("PhotonTreeWriterPresel");
  phottreepresel->SetPhotonsFromBranch(kFALSE);
  phottreepresel->SetInputPhotonsName(photpresel->GetOutputName());
  phottreepresel->SetEnableJets(kTRUE);
  phottreepresel->SetPFJetsFromBranch(kFALSE);
  phottreepresel->SetPFJetName(jetCorr->GetOutputName());  
  phottreepresel->SetExcludeDoublePrompt(excludedoubleprompt);  
  phottreepresel->SetIsData(isData);  
  
  phottreepresel->SetApplyLeptonTag(true);
  phottreepresel->SetLeptonTagElectronsName("HggLeptonTagElectrons");
  phottreepresel->SetLeptonTagMuonsName    ("HggLeptonTagMuons");

  //------------------------------------------------------------------------------------------------
  // making analysis chain
  //------------------------------------------------------------------------------------------------
  // this is how it always starts
  runLumiSel      ->Add(mcselmod);
    
  if (TString(dataset).Contains("-h")) {    
    mcselmod        ->Add(sysMod);
  }
  
  mcselmod         ->Add(hltModP);

  hltModP         ->Add(goodPVFilterMod);
  goodPVFilterMod->Add(photreg);
  photreg->Add(pubJet);
  pubJet->Add(jetCorr);
  jetCorr->Add(eleIdMod);
  eleIdMod->Add(muonIdMod);
  
  muonIdMod          ->Add(photcic);
  muonIdMod          ->Add(photpresel);  

  photcic         ->Add(phottreecic);
  photpresel      ->Add(phottreepresel);

  TFile::SetOpenTimeout(0);
  TFile::SetCacheFileDir("./rootfilecache",kTRUE,kTRUE);
  TFile::SetReadaheadSize(128*1024*1024);
  
  //------------------------------------------------------------------------------------------------
  // setup analysis
  //------------------------------------------------------------------------------------------------
  Analysis *ana = new Analysis;
  ana->SetUseHLT(kTRUE);
  ana->SetKeepHierarchy(kTRUE);
  ana->SetSuperModule(runLumiSel);
  ana->SetPrintScale(100);
  if (nEvents >= 0)
    ana->SetProcessNEvents(nEvents);

  //------------------------------------------------------------------------------------------------
  // organize input
  //------------------------------------------------------------------------------------------------
  Catalog *c = new Catalog(catalogDir);
  TString skimdataset = TString(dataset)+TString("/") +TString(skim);
  Dataset *d = NULL;
  TString bookstr = book;
  //if (TString(dataset).Contains("s11-h")) bookstr.ReplaceAll("local","t2mit");
  if (TString(skim).CompareTo("noskim") == 0)
    d = c->FindDataset(bookstr,dataset,fileset);
  else 
    d = c->FindDataset(bookstr,skimdataset.Data(),fileset);
  ana->AddDataset(d);
  //ana->AddFile("/mnt/hadoop/cmsprod/filefi/025/r11b-pho-n30-v1/4863E9D1-BC1C-E111-99BA-001A92810AF4.root");

  //------------------------------------------------------------------------------------------------
  // organize output
  //------------------------------------------------------------------------------------------------
  TString rootFile = TString(outputName);
  rootFile += TString("_") + TString(dataset) + TString("_") + TString(skim);
  if (TString(fileset) != TString(""))
    rootFile += TString("_") + TString(fileset);
  rootFile += TString(".root");
  ana->SetOutputName(rootFile.Data());
  //ana->SetCacheSize(64*1024*1024);
  ana->SetCacheSize(0);
  
  //------------------------------------------------------------------------------------------------
  // Say what we are doing
  //------------------------------------------------------------------------------------------------
  printf("\n==== PARAMETER SUMMARY FOR THIS JOB ====\n");
  printf("\n JSON file: %s\n  and overlap cut: %f (%s)\n",jsonFile.Data(),overlapCut,overlap);
  printf("\n Rely on Catalog: %s\n",catalogDir);
  printf("  -> Book: %s  Dataset: %s  Skim: %s  Fileset: %s <-\n",book,dataset,skim,fileset);
  printf("\n Root output: %s\n\n",rootFile.Data());  
  printf("\n========================================\n");

  //------------------------------------------------------------------------------------------------
  // run the analysis after successful initialisation
  //------------------------------------------------------------------------------------------------
  ana->Run(!gROOT->IsBatch());

  return;
}
