// $Id: runHgg2012HCP.C,v 1.2 2012/10/24 14:47:27 mingyang Exp $
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
#include "MitPhysics/Mods/interface/SeparatePileUpMod.h"

#endif

//--------------------------------------------------------------------------------------------------
void runHgg2012HCP(const char *fileset    = "0000",
		   const char *skim       = "noskim",
		   //const char *dataset = "f11--h120gg-gf-v14b-pu",
		   //const char *dataset = "r11a-pho-j16-v1",   
		   //const char *dataset = "meridiani2012-diphoj-v9",
		   //const char *dataset = "s12-h125gg-gf-v9",
		   //const char *dataset = "s12-pj40-2em-v9",
		   //const char *dataset = "s12-diphoj-3-v9",
		   //const char *dataset = "s12-zllm50-v9",
		   //const char *dataset = "f11--h121gg-gf-v14b-pu",
		   //const char *dataset = "r12a-pho-pr-v1",
		   //const char *dataset = "s12-pj40-2em-v9",
		   //const char *dataset = "s12-h125gg-vbf-v7a", 
		   const char *dataset = "s12-h125gg-vh-v7a", 
		   //const char *dataset = "s12-h125gg-gf-v7a", 
		   //const char *dataset = "r12a-pho-j13-v1", 
		   //const char *dataset = "r12b-dph-j13-v1", 
		   const char *book       = "t2mit/filefi/029",
		   const char *catalogDir = "/home/cmsprod/catalog",
		   const char *outputName = "hgg",
		   int         nEvents    = 1000)
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
  HLTMod *hltModM = new HLTMod("HLTModM");
  hltModM->AddTrigger("HLT_Mu9");
  hltModM->AddTrigger("HLT_Mu11");
  hltModM->AddTrigger("HLT_Mu15_v1");
  hltModM->SetTrigObjsName("MyHltMuonObjs");
  hltModM->SetAbortIfNotAccepted(kFALSE);

  HLTMod *hltModE = new HLTMod("HLTModE");
  hltModE->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v1",150000,161176);
  hltModE->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v2",161179,163261);
  hltModE->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v3",163262,164237);
  hltModE->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v4",165085,165888);
  hltModE->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v5",165900,166967);
  hltModE->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v6",166968,170053);
  hltModE->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6",170054,170759);
  hltModE->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7",170760,173198);
  hltModE->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v8",173199,178380);
  hltModE->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9",178381,179889);
  hltModE->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v10",179890,189999);
  hltModE->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*",190000,999999);
  hltModE->SetTrigObjsName("MyHltElecObjs");
  hltModE->SetAbortIfNotAccepted(isData);

  HLTMod *hltModES = new HLTMod("HLTModES");
  hltModES->AddTrigger("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v*",160000,189999);
  hltModES->AddTrigger("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2",160000,189999);
  hltModES->AddTrigger("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v3",160000,189999);
  hltModES->AddTrigger("HLT_Ele25_WP80_PFMT40_v*",160000,189999);
  hltModES->AddTrigger("HLT_Ele27_WP80_PFMT50_v*",160000,189999);
  hltModES->AddTrigger("HLT_Ele27_WP80_v*",190000,999999);
  hltModES->SetTrigObjsName("MyHltElecSObjs");
  hltModES->SetAbortIfNotAccepted(isData);

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
  
  hltModP->AddTrigger("HLT_Photon26_CaloIdXL_IsoXL_Photon18_CaloIdXL_IsoXL_Mass60_v*",178420,189999);
  hltModP->AddTrigger("HLT_Photon26_CaloIdXL_IsoXL_Photon18_R9IdT_Mass60_v*",178420,189999);
  hltModP->AddTrigger("HLT_Photon26_R9IdT_Photon18_CaloIdXL_IsoXL_Mass60_v*",178420,189999);
  hltModP->AddTrigger("HLT_Photon26_R9IdT_Photon18_R9IdT_Mass60_v*",178420,189999);
  hltModP->AddTrigger("HLT_Photon36_CaloIdL_IsoVL_Photon22_CaloIdL_IsoVL_v*",178420,189999);
  hltModP->AddTrigger("HLT_Photon36_CaloIdL_IsoVL_Photon22_R9Id_v*",178420,189999);
  hltModP->AddTrigger("HLT_Photon36_R9Id_Photon22_CaloIdL_IsoVL_v*",178420,189999);
  hltModP->AddTrigger("HLT_Photon36_R9Id_Photon22_R9Id_v*",178420,189999);

  hltModP->AddTrigger("HLT_Photon26_CaloId10_Iso50_Photon18_CaloId10_Iso50_Mass60_v*",190000,999999); 
  hltModP->AddTrigger("HLT_Photon26_CaloId10_Iso50_Photon18_R9Id85_Mass60_v*",190000,999999); 
  hltModP->AddTrigger("HLT_Photon26_R9Id85_OR_CaloId10_Iso50_Photon18_R9Id85_OR_CaloId10_Iso50_Mass60_v*",190000,999999); 
  hltModP->AddTrigger("HLT_Photon26_R9Id85_OR_CaloId10_Iso50_Photon18_R9Id85_OR_CaloId10_Iso50_Mass70_v*",190000,999999); 
  hltModP->AddTrigger("HLT_Photon26_R9Id85_Photon18_CaloId10_Iso50_Mass60_v*",190000,999999); 
  hltModP->AddTrigger("HLT_Photon26_R9Id85_Photon18_R9Id85_Mass60_v*",190000,999999); 
  hltModP->AddTrigger("HLT_Photon36_CaloId10_Iso50_Photon22_CaloId10_Iso50_v*",190000,999999); 
  hltModP->AddTrigger("HLT_Photon36_CaloId10_Iso50_Photon22_R9Id85_v*",190000,999999); 
  hltModP->AddTrigger("HLT_Photon36_R9Id85_OR_CaloId10_Iso50_Photon22_R9Id85_OR_CaloId10_Iso50_v*",190000,999999); 
  hltModP->AddTrigger("HLT_Photon36_R9Id85_Photon22_CaloId10_Iso50_v*",190000,999999); 
  hltModP->AddTrigger("HLT_Photon36_R9Id85_Photon22_R9Id85_v*",190000,999999); 
    
  hltModP->SetTrigObjsName("MyHltPhotObjs");
  hltModP->SetAbortIfNotAccepted(isData);
  //------------------------------------------------------------------------------------------------
  // split pfcandidates to PFPU and PFnoPU
  //------------------------------------------------------------------------------------------------
  SeparatePileUpMod* SepPUMod = new SeparatePileUpMod;
  //  SepPUMod->SetUseAllVerteces(kFALSE);
  // SepPUMod->SetVertexName("OutVtxCiC");
  SepPUMod->SetPFNoPileUpName("pfnopileupcands");
  SepPUMod->SetPFPileUpName("pfpileupcands");
  SepPUMod->SetCheckClosestZVertex(kFALSE);
  
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

  GoodPVFilterMod *goodPVFilterModE = new GoodPVFilterMod("GoodPVFilterModE");
  goodPVFilterModE->SetOutputName("GoodVertexesE");
  goodPVFilterModE->SetMinVertexNTracks(0);
  goodPVFilterModE->SetMinNDof         (4.0);
  goodPVFilterModE->SetMaxAbsZ         (24.0);
  goodPVFilterModE->SetMaxRho          (2.0);  
  goodPVFilterModE->SetIsMC(!isData);
  
  GoodPVFilterMod *goodPVFilterModES = new GoodPVFilterMod("GoodPVFilterModES");
  goodPVFilterModES->SetOutputName("GoodVertexesES");
  goodPVFilterModES->SetMinVertexNTracks(0);
  goodPVFilterModES->SetMinNDof         (4.0);
  goodPVFilterModES->SetMaxAbsZ         (24.0);
  goodPVFilterModES->SetMaxRho          (2.0);    
  

  //------------------------------------------------------------------------------------------------
  // object id and cleaning sequence
  //------------------------------------------------------------------------------------------------
  //-----------------------------------
  // Lepton Selection 
  //-----------------------------------
  ElectronIDMod* eleIdMod = new ElectronIDMod;
  eleIdMod -> SetPtMin(20);  
  eleIdMod -> SetEtaMax(2.5);
  eleIdMod -> SetApplyEcalFiducial(true);
  eleIdMod -> SetIDType("Hgg_LeptonTag_2012IdHCP");  
  eleIdMod -> SetElectronMVAWeightsSubdet0Pt10To20(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/Electrons_BDTG_NonTrigV0_Cat1.weights.xml"))));
  eleIdMod -> SetElectronMVAWeightsSubdet1Pt10To20(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/Electrons_BDTG_NonTrigV0_Cat2.weights.xml"))));  
  eleIdMod -> SetElectronMVAWeightsSubdet2Pt10To20(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/Electrons_BDTG_NonTrigV0_Cat3.weights.xml"))));  
  eleIdMod -> SetElectronMVAWeightsSubdet0Pt20ToInf(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/Electrons_BDTG_NonTrigV0_Cat4.weights.xml")))); 
  eleIdMod -> SetElectronMVAWeightsSubdet1Pt20ToInf(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/Electrons_BDTG_NonTrigV0_Cat5.weights.xml")))); 
  eleIdMod -> SetElectronMVAWeightsSubdet2Pt20ToInf(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/Electrons_BDTG_NonTrigV0_Cat6.weights.xml")))); 
  eleIdMod -> SetWhichVertex(-1);
  eleIdMod -> SetD0Cut(0.02);
  eleIdMod -> SetDZCut(0.2); //h  
  eleIdMod -> SetIsoType("PFIso_HggLeptonTag2012HCP"); //h
  eleIdMod -> SetOutputName("HggLeptonTagElectrons");
  eleIdMod -> SetRhoType(RhoUtilities::CMS_RHO_RHOKT6PFJETS);
  eleIdMod -> SetPFNoPileUpName("pfnopileupcands");
  eleIdMod -> SetInvertNExpectedHitsInnerCut(kFALSE);
  eleIdMod -> SetNExpectedHitsInnerCut(1);   
  eleIdMod -> SetApplyConversionFilterType1(kTRUE);
  
  MuonIDMod* muonIdMod = new MuonIDMod;
  // base kinematics
  muonIdMod -> SetPtMin(20.);
  muonIdMod -> SetEtaCut(2.4);
  // base ID
  muonIdMod -> SetIDType("Tight");
  muonIdMod -> SetWhichVertex(-1); // this is a 'hack'.. but hopefully good enough...
  muonIdMod -> SetD0Cut(0.2);
  muonIdMod -> SetDZCut(0.5);
  muonIdMod -> SetIsoType("PFIsoBetaPUCorrected"); //h
  muonIdMod -> SetPFIsoCut(0.2); //h
  muonIdMod -> SetOutputName("HggLeptonTagMuons");
  muonIdMod -> SetPFNoPileUpName("pfnopileupcands");
  muonIdMod -> SetPFPileUpName("pfpileupcands");
  
  PublisherMod<PFJet,Jet> *pubJet = new PublisherMod<PFJet,Jet>("JetPub");
  pubJet->SetInputName("AKt5PFJets");
  pubJet->SetOutputName("PubAKt5PFJets");
  
  PublisherMod<PFJet,Jet> *pubJetOpen = new PublisherMod<PFJet,Jet>("JetPubOpen");
  pubJetOpen->SetInputName("AKt5PFJets");
  pubJetOpen->SetOutputName("PubAKt5PFJetsOpen");  
  
  JetCorrectionMod *jetCorr = new JetCorrectionMod;
  if(isData){ 
    jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer12_V7_DATA_L1FastJet_AK5PF.txt")).Data())); 
    jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer12_V7_DATA_L2Relative_AK5PF.txt")).Data())); 
    jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer12_V7_DATA_L3Absolute_AK5PF.txt")).Data())); 
    jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer12_V7_DATA_L2L3Residual_AK5PF.txt")).Data())); 
  }
  else {
    jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer12_V7_MC_L1FastJet_AK5PF.txt")).Data())); 
    jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer12_V7_MC_L2Relative_AK5PF.txt")).Data())); 
    jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer12_V7_MC_L3Absolute_AK5PF.txt")).Data())); 
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
  
  Bool_t is25 = kFALSE;
  if (TString(book).Contains("025")) is25 = kTRUE;
  
  PhotonMvaMod *photreg = new PhotonMvaMod;
  photreg->SetRegressionVersion(3);
  photreg->SetRegressionWeights(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/gbrv3ph_52x.root")).Data()));
  photreg->SetOutputName("GoodPhotonsRegr");
  photreg->SetApplyShowerRescaling(kTRUE);
  photreg->SetIsData(isData);

  
  PhotonPairSelector         *photcic = new PhotonPairSelector("PhotonPairSelectorCiC");
  photcic->SetOutputName("GoodPhotonsCIC");
  photcic->SetOutputVtxName("OutVtxCiC");        
  photcic->SetPhotonSelType("CiCPFSelection");
  photcic->SetVertexSelType("CiCMVA2012Selection");
  photcic->SetUseSingleLegConversions(kTRUE);
  photcic->DoMCSmear(kTRUE);
  photcic->DoDataEneCorr(kTRUE);
  photcic->SetPhotonsFromBranch(kFALSE);
  photcic->SetInputPhotonsName(photreg->GetOutputName());
  
  //------------------------------------------2012 HCP--------------------------------------------------------------
  photcic->SetMCSmearFactors2012HCP(0.0099,0.0098,0.0113,0.0094,0.0195,0.0198,0.0291,0.0275,0.0325,0.0356);
  //photcic->SetMCSmearFactors2012HCP(0.01,0.011,0.0106,0.0106,0.0186,0.0196,0.0283,0.0267,0.0343,0.0345);
  photcic->AddEnCorrPerRun2012HCP(190645,190781,0.9964,0.9964,1.0020,1.0020,0.9893,1.0028,0.9871,0.9937,0.9839,0.9958);
  photcic->AddEnCorrPerRun2012HCP(190782,191042,1.0024,1.0024,1.0079,1.0079,0.9923,1.0058,0.9911,0.9977,0.9886,1.0005);
  photcic->AddEnCorrPerRun2012HCP(191043,193555,0.9935,0.9935,0.9991,0.9991,0.9861,0.9997,0.9894,0.9960,0.9864,0.9982);
  photcic->AddEnCorrPerRun2012HCP(193556,194150,0.9920,0.9920,0.9976,0.9976,0.9814,0.9951,0.9896,0.9962,0.9872,0.9990);
  photcic->AddEnCorrPerRun2012HCP(194151,194532,0.9925,0.9925,0.9981,0.9981,0.9826,0.9963,0.9914,0.9980,0.9874,0.9993);
  photcic->AddEnCorrPerRun2012HCP(194533,195113,0.9927,0.9927,0.9983,0.9983,0.9844,0.9981,0.9934,0.9999,0.9878,0.9996);
  photcic->AddEnCorrPerRun2012HCP(195114,195915,0.9929,0.9929,0.9984,0.9984,0.9838,0.9974,0.9942,1.0007,0.9878,0.9997);
  photcic->AddEnCorrPerRun2012HCP(195916,198115,0.9919,0.9919,0.9975,0.9975,0.9827,0.9964,0.9952,1.0017,0.9869,0.9987);
  photcic->AddEnCorrPerRun2012HCP(198116,199803,0.9955,0.9955,1.0011,1.0011,0.9859,0.9995,0.9893,0.9959,0.9923,1.0041);
  photcic->AddEnCorrPerRun2012HCP(199804,200048,0.9967,0.9967,1.0023,1.0023,0.9870,1.0006,0.9893,0.9959,0.9937,1.0055);
  photcic->AddEnCorrPerRun2012HCP(200049,200151,0.9980,0.9980,1.0036,1.0036,0.9877,1.0012,0.9910,0.9976,0.9980,1.0097);
  photcic->AddEnCorrPerRun2012HCP(200152,200490,0.9958,0.9958,1.0013,1.0013,0.9868,1.0004,0.9922,0.9988,0.9948,1.0065);
  photcic->AddEnCorrPerRun2012HCP(200491,200531,0.9979,0.9979,1.0035,1.0035,0.9876,1.0012,0.9915,0.9981,0.9979,1.0096);
  photcic->AddEnCorrPerRun2012HCP(200532,201656,0.9961,0.9961,1.0017,1.0017,0.9860,0.9996,0.9904,0.9970,0.9945,1.0063);
  photcic->AddEnCorrPerRun2012HCP(201657,202305,0.9969,0.9969,1.0025,1.0025,0.9866,1.0002,0.9914,0.9980,0.9999,1.0116);
  photcic->AddEnCorrPerRun2012HCP(202305,203002,0.9982,0.9982,1.0038,1.0038,0.9872,1.0008,0.9934,1.0000,1.0018,1.0135);

  //-----------------------------------------------------------------------------------------------------------------
  
  //  photcic->SetDoMCR9Scaling(kTRUE);
  //  photcic->SetMCR9Scale(1.0035, 1.0035);
  photcic->SetDoShowerShapeScaling(kTRUE); 
  photcic->SetShowerShapeType("2012ShowerShape");
  //photcic->SetDoMCErrScaling(kTRUE);
  //photcic->SetMCErrScale(1.07, 1.045); 
  //photcic->SetMCErrScale(1, 1); //ming:scale(sigE/E)
  photcic->SetJetsName(jetCorr->GetOutputName());    
  //photcic->SetRescaledBeamspotWidth(5.0);
  photcic->SetIsData(isData);
  photcic->SetApplyLeptonTag(kTRUE);
  photcic->SetLeptonTagElectronsName("HggLeptonTagElectrons");
  photcic->SetLeptonTagMuonsName("HggLeptonTagMuons");  
  photcic->Set2012HCP(kTRUE);
  
  PhotonPairSelector         *photcicnoeleveto = new PhotonPairSelector("PhotonPairSelectorCiCInvertEleVeto");
  photcicnoeleveto->SetOutputName("GoodPhotonsCICNoEleVeto");
  photcicnoeleveto->SetOutputVtxName("OutVtxCiCInvertEleVeto");      
  photcicnoeleveto->SetPhotonSelType("CiCPFSelection");
  photcicnoeleveto->SetVertexSelType("CiCMVA2012Selection");
  photcicnoeleveto->SetUseSingleLegConversions(kTRUE);
  photcicnoeleveto->DoMCSmear(kTRUE);
  photcicnoeleveto->DoDataEneCorr(kTRUE);
  photcicnoeleveto->SetPhotonsFromBranch(kFALSE);
  photcicnoeleveto->SetInputPhotonsName(photreg->GetOutputName());
  
  //------------------------------------------2012 HCP--------------------------------------------------------------
  photcicnoeleveto->SetMCSmearFactors2012HCP(0.0099,0.0098,0.0113,0.0094,0.0195,0.0198,0.0291,0.0275,0.0325,0.0356);
  //photcicnoeleveto->SetMCSmearFactors2012HCP(0.01,0.011,0.0106,0.0106,0.0186,0.0196,0.0283,0.0267,0.0343,0.0345);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(190645,190781,0.9964,0.9964,1.0020,1.0020,0.9893,1.0028,0.9871,0.9937,0.9839,0.9958);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(190782,191042,1.0024,1.0024,1.0079,1.0079,0.9923,1.0058,0.9911,0.9977,0.9886,1.0005);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(191043,193555,0.9935,0.9935,0.9991,0.9991,0.9861,0.9997,0.9894,0.9960,0.9864,0.9982);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(193556,194150,0.9920,0.9920,0.9976,0.9976,0.9814,0.9951,0.9896,0.9962,0.9872,0.9990);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(194151,194532,0.9925,0.9925,0.9981,0.9981,0.9826,0.9963,0.9914,0.9980,0.9874,0.9993);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(194533,195113,0.9927,0.9927,0.9983,0.9983,0.9844,0.9981,0.9934,0.9999,0.9878,0.9996);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(195114,195915,0.9929,0.9929,0.9984,0.9984,0.9838,0.9974,0.9942,1.0007,0.9878,0.9997);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(195916,198115,0.9919,0.9919,0.9975,0.9975,0.9827,0.9964,0.9952,1.0017,0.9869,0.9987);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(198116,199803,0.9955,0.9955,1.0011,1.0011,0.9859,0.9995,0.9893,0.9959,0.9923,1.0041);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(199804,200048,0.9967,0.9967,1.0023,1.0023,0.9870,1.0006,0.9893,0.9959,0.9937,1.0055);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(200049,200151,0.9980,0.9980,1.0036,1.0036,0.9877,1.0012,0.9910,0.9976,0.9980,1.0097);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(200152,200490,0.9958,0.9958,1.0013,1.0013,0.9868,1.0004,0.9922,0.9988,0.9948,1.0065);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(200491,200531,0.9979,0.9979,1.0035,1.0035,0.9876,1.0012,0.9915,0.9981,0.9979,1.0096);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(200532,201656,0.9961,0.9961,1.0017,1.0017,0.9860,0.9996,0.9904,0.9970,0.9945,1.0063);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(201657,202305,0.9969,0.9969,1.0025,1.0025,0.9866,1.0002,0.9914,0.9980,0.9999,1.0116);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(202305,203002,0.9982,0.9982,1.0038,1.0038,0.9872,1.0008,0.9934,1.0000,1.0018,1.0135);
  
  //-----------------------------------------------------------------------------------------------------------------
  
  //photcicnoeleveto->SetDoMCR9Scaling(kTRUE);
  //photcicnoeleveto->SetMCR9Scale(1.0035, 1.0035);
  photcicnoeleveto->SetDoShowerShapeScaling(kTRUE);
  photcicnoeleveto->SetShowerShapeType("2012ShowerShape");
  //photcicnoeleveto->SetDoMCErrScaling(kTRUE);
  //photcicnoeleveto->SetMCErrScale(1.07, 1.045);    
  //photcicnoeleveto->SetMCErrScale(1, 1);    
  photcicnoeleveto->SetApplyEleVeto(kFALSE);
  photcicnoeleveto->SetInvertElectronVeto(kTRUE);
  photcicnoeleveto->SetJetsName(jetCorr->GetOutputName());  
  //photcicnoeleveto->SetRescaledBeamspotWidth(5.0);
  photcicnoeleveto->SetIsData(isData);
  photcicnoeleveto->SetApplyLeptonTag(kTRUE);
  photcicnoeleveto->SetLeptonTagElectronsName("HggLeptonTagElectrons");
  photcicnoeleveto->SetLeptonTagMuonsName("HggLeptonTagMuons");  
  photcicnoeleveto->Set2012HCP(kTRUE);  
  
  PhotonPairSelector         *photpresel = new PhotonPairSelector("PhotonPairSelectorPresel");
  photpresel->SetOutputName("GoodPhotonsPresel");
  photpresel->SetPhotonSelType("MITPFSelection");
  photpresel->SetVertexSelType("CiCMVA2012Selection");
  photpresel->SetUseSingleLegConversions(kTRUE);
  photpresel->SetIdMVAType("2012IdMVA_globe");
  photpresel->DoMCSmear(kTRUE);
  photpresel->DoDataEneCorr(kTRUE);
  photpresel->SetPhotonsFromBranch(kFALSE);
  photpresel->SetInputPhotonsName(photreg->GetOutputName());
  //------------------------------------------2012 HCP--------------------------------------------------------------
  photpresel->SetMCSmearFactors2012HCP(0.0099,0.0098,0.0113,0.0094,0.0195,0.0198,0.0291,0.0275,0.0325,0.0356);
  //photpresel->SetMCSmearFactors2012HCP(0.01,0.011,0.0106,0.0106,0.0186,0.0196,0.0283,0.0267,0.0343,0.0345);
  photpresel->AddEnCorrPerRun2012HCP(190645,190781,0.9964,0.9964,1.0020,1.0020,0.9893,1.0028,0.9871,0.9937,0.9839,0.9958);
  photpresel->AddEnCorrPerRun2012HCP(190782,191042,1.0024,1.0024,1.0079,1.0079,0.9923,1.0058,0.9911,0.9977,0.9886,1.0005);
  photpresel->AddEnCorrPerRun2012HCP(191043,193555,0.9935,0.9935,0.9991,0.9991,0.9861,0.9997,0.9894,0.9960,0.9864,0.9982);
  photpresel->AddEnCorrPerRun2012HCP(193556,194150,0.9920,0.9920,0.9976,0.9976,0.9814,0.9951,0.9896,0.9962,0.9872,0.9990);
  photpresel->AddEnCorrPerRun2012HCP(194151,194532,0.9925,0.9925,0.9981,0.9981,0.9826,0.9963,0.9914,0.9980,0.9874,0.9993);
  photpresel->AddEnCorrPerRun2012HCP(194533,195113,0.9927,0.9927,0.9983,0.9983,0.9844,0.9981,0.9934,0.9999,0.9878,0.9996);
  photpresel->AddEnCorrPerRun2012HCP(195114,195915,0.9929,0.9929,0.9984,0.9984,0.9838,0.9974,0.9942,1.0007,0.9878,0.9997);
  photpresel->AddEnCorrPerRun2012HCP(195916,198115,0.9919,0.9919,0.9975,0.9975,0.9827,0.9964,0.9952,1.0017,0.9869,0.9987);
  photpresel->AddEnCorrPerRun2012HCP(198116,199803,0.9955,0.9955,1.0011,1.0011,0.9859,0.9995,0.9893,0.9959,0.9923,1.0041);
  photpresel->AddEnCorrPerRun2012HCP(199804,200048,0.9967,0.9967,1.0023,1.0023,0.9870,1.0006,0.9893,0.9959,0.9937,1.0055);
  photpresel->AddEnCorrPerRun2012HCP(200049,200151,0.9980,0.9980,1.0036,1.0036,0.9877,1.0012,0.9910,0.9976,0.9980,1.0097);
  photpresel->AddEnCorrPerRun2012HCP(200152,200490,0.9958,0.9958,1.0013,1.0013,0.9868,1.0004,0.9922,0.9988,0.9948,1.0065);
  photpresel->AddEnCorrPerRun2012HCP(200491,200531,0.9979,0.9979,1.0035,1.0035,0.9876,1.0012,0.9915,0.9981,0.9979,1.0096);
  photpresel->AddEnCorrPerRun2012HCP(200532,201656,0.9961,0.9961,1.0017,1.0017,0.9860,0.9996,0.9904,0.9970,0.9945,1.0063);
  photpresel->AddEnCorrPerRun2012HCP(201657,202305,0.9969,0.9969,1.0025,1.0025,0.9866,1.0002,0.9914,0.9980,0.9999,1.0116);
  photpresel->AddEnCorrPerRun2012HCP(202305,203002,0.9982,0.9982,1.0038,1.0038,0.9872,1.0008,0.9934,1.0000,1.0018,1.0135);
  //-----------------------------------------------------------------------------------------------------------------
  //photpresel->SetDoMCR9Scaling(kTRUE);
  photpresel->SetDoShowerShapeScaling(kTRUE);
  photpresel->SetShowerShapeType("2012ShowerShape");
  //photpresel->SetDoMCErrScaling(kTRUE);
  //photpresel->SetMCErrScale(1.07, 1.045); 
  //photpresel->SetMCErrScale(1, 1); 
  photpresel->SetJetsName(jetCorr->GetOutputName()); 
  //photpresel->SetRescaledBeamspotWidth(5.0);
  photpresel->SetIsData(isData);
  photpresel->SetApplyLeptonTag(kTRUE);
  photpresel->SetLeptonTagElectronsName("HggLeptonTagElectrons");
  photpresel->SetLeptonTagMuonsName("HggLeptonTagMuons");  
  photpresel->Set2012HCP(kTRUE); 

  PhotonPairSelector         *photpreselinverteleveto = new PhotonPairSelector("PhotonPairSelectorPreselInvertEleVeto");
  photpreselinverteleveto->SetOutputName("GoodPhotonsPreselInvertEleVeto");
  photpreselinverteleveto->SetOutputVtxName("OutVtxPreselInvertEleVeto");    
  photpreselinverteleveto->SetPhotonSelType("MITPFSelection");
  photpreselinverteleveto->SetIdMVAType("2012IdMVA_globe");
  photpreselinverteleveto->SetVertexSelType("CiCMVA2012Selection");
  photpreselinverteleveto->SetUseSingleLegConversions(kTRUE);
  photpreselinverteleveto->DoMCSmear(kTRUE);
  photpreselinverteleveto->DoDataEneCorr(kTRUE);
  photpreselinverteleveto->SetPhotonsFromBranch(kFALSE);
  photpreselinverteleveto->SetInputPhotonsName(photreg->GetOutputName());
  //------------------------------------------2012 HCP--------------------------------------------------------------
  photpreselinverteleveto->SetMCSmearFactors2012HCP(0.0099,0.0098,0.0113,0.0094,0.0195,0.0198,0.0291,0.0275,0.0325,0.0356);
  //photpreselinverteleveto->SetMCSmearFactors2012HCP(0.01,0.011,0.0106,0.0106,0.0186,0.0196,0.0283,0.0267,0.0343,0.0345);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(190645,190781,0.9964,0.9964,1.0020,1.0020,0.9893,1.0028,0.9871,0.9937,0.9839,0.9958);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(190782,191042,1.0024,1.0024,1.0079,1.0079,0.9923,1.0058,0.9911,0.9977,0.9886,1.0005);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(191043,193555,0.9935,0.9935,0.9991,0.9991,0.9861,0.9997,0.9894,0.9960,0.9864,0.9982);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(193556,194150,0.9920,0.9920,0.9976,0.9976,0.9814,0.9951,0.9896,0.9962,0.9872,0.9990);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(194151,194532,0.9925,0.9925,0.9981,0.9981,0.9826,0.9963,0.9914,0.9980,0.9874,0.9993);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(194533,195113,0.9927,0.9927,0.9983,0.9983,0.9844,0.9981,0.9934,0.9999,0.9878,0.9996);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(195114,195915,0.9929,0.9929,0.9984,0.9984,0.9838,0.9974,0.9942,1.0007,0.9878,0.9997);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(195916,198115,0.9919,0.9919,0.9975,0.9975,0.9827,0.9964,0.9952,1.0017,0.9869,0.9987);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(198116,199803,0.9955,0.9955,1.0011,1.0011,0.9859,0.9995,0.9893,0.9959,0.9923,1.0041);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(199804,200048,0.9967,0.9967,1.0023,1.0023,0.9870,1.0006,0.9893,0.9959,0.9937,1.0055);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(200049,200151,0.9980,0.9980,1.0036,1.0036,0.9877,1.0012,0.9910,0.9976,0.9980,1.0097);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(200152,200490,0.9958,0.9958,1.0013,1.0013,0.9868,1.0004,0.9922,0.9988,0.9948,1.0065);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(200491,200531,0.9979,0.9979,1.0035,1.0035,0.9876,1.0012,0.9915,0.9981,0.9979,1.0096);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(200532,201656,0.9961,0.9961,1.0017,1.0017,0.9860,0.9996,0.9904,0.9970,0.9945,1.0063);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(201657,202305,0.9969,0.9969,1.0025,1.0025,0.9866,1.0002,0.9914,0.9980,0.9999,1.0116);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(202305,203002,0.9982,0.9982,1.0038,1.0038,0.9872,1.0008,0.9934,1.0000,1.0018,1.0135);
  //-----------------------------------------------------------------------------------------------------------------
  photpreselinverteleveto->SetShowerShapeType("2012ShowerShape");
  photpreselinverteleveto->SetDoShowerShapeScaling(kTRUE);  
  //photpreselinverteleveto->SetDoMCErrScaling(kTRUE);
  //photpreselinverteleveto->SetMCErrScale(1.07, 1.045);  
  //photpreselinverteleveto->SetMCErrScale(1, 1);  
  photpreselinverteleveto->SetApplyEleVeto(kFALSE);
  photpreselinverteleveto->SetInvertElectronVeto(kTRUE);
  photpreselinverteleveto->SetJetsName(jetCorr->GetOutputName());  
  //photpreselinverteleveto->SetRescaledBeamspotWidth(5.0);   
  photpreselinverteleveto->SetIsData(isData);
  photpreselinverteleveto->SetApplyLeptonTag(kTRUE);
  photpreselinverteleveto->SetLeptonTagElectronsName("HggLeptonTagElectrons");
  photpreselinverteleveto->SetLeptonTagMuonsName("HggLeptonTagMuons");    
  photpreselinverteleveto->Set2012HCP(kTRUE);  
  

   PhotonPairSelector         *photpreselnosmear = new PhotonPairSelector("PhotonPairSelectorPreselNoSmear");
   photpreselnosmear->SetOutputName("GoodPhotonsPreselNoSmear");
   photpreselnosmear->SetPhotonSelType("MITPFSelection");
   photpreselnosmear->SetVertexSelType("CiCMVA2012Selection");
   photpreselnosmear->SetUseSingleLegConversions(kTRUE);  
   photpreselnosmear->SetIdMVAType("2012IdMVA_globe");
   photpreselnosmear->SetShowerShapeType("2012ShowerShape");
   photpreselnosmear->SetDoShowerShapeScaling(kTRUE);  
   photpreselnosmear->SetPhotonsFromBranch(kFALSE);
   photpreselnosmear->SetInputPhotonsName(photreg->GetOutputName());
   photpreselnosmear->SetJetsName(jetCorr->GetOutputName());  
   photpreselnosmear->SetOutputVtxName("OutVtxNoSmear");  
   photpreselnosmear->SetLeadingPtMin(30.);
   photpreselnosmear->SetTrailingPtMin(22.);
   //photpreselnosmear->SetRescaledBeamspotWidth(5.0);      
   photpreselnosmear->SetIsData(isData);  
   photpreselnosmear->SetApplyLeptonTag(kTRUE);
   photpreselnosmear->SetLeptonTagElectronsName("HggLeptonTagElectrons");
   photpreselnosmear->SetLeptonTagMuonsName("HggLeptonTagMuons");
   photpreselnosmear->Set2012HCP(kTRUE); 
   photpreselnosmear->DoMCSmear(kFALSE);
   photpreselnosmear->DoDataEneCorr(kFALSE);
   
   PhotonTreeWriter *phottreecic = new PhotonTreeWriter("PhotonTreeWriterCiC");
   phottreecic->SetPhotonsFromBranch(kFALSE);
   phottreecic->SetInputPhotonsName(photcic->GetOutputName());
   phottreecic->SetEnableJets(kTRUE);
   phottreecic->SetApplyJetId(kTRUE);
   phottreecic->SetPFJetsFromBranch(kFALSE);
   phottreecic->SetPFJetName(jetCorr->GetOutputName());
   phottreecic->SetExcludeDoublePrompt(excludedoubleprompt);
   phottreecic->SetIsData(isData);
   if (is25) phottreecic->SetEnablePFPhotons(kFALSE);
   phottreecic->SetApplyLeptonTag(kTRUE);
   phottreecic->SetLeptonTagElectronsName("HggLeptonTagElectrons");
   phottreecic->SetLeptonTagMuonsName("HggLeptonTagMuons");
   phottreecic->SetApplyVBFTag(kTRUE);
   phottreecic->SetApplyPFMetCorr(kTRUE);
   phottreecic->SetBeamspotWidth(5.0);
   
   PhotonTreeWriter *phottreecicnoeleveto = new PhotonTreeWriter("PhotonTreeWriterCiCInvertEleVeto");
   phottreecicnoeleveto->SetPhotonsFromBranch(kFALSE);
   phottreecicnoeleveto->SetInputPhotonsName(photcicnoeleveto->GetOutputName());
   phottreecicnoeleveto->SetEnableJets(kTRUE);
   phottreecicnoeleveto->SetApplyJetId(kTRUE);  
   phottreecicnoeleveto->SetPFJetsFromBranch(kFALSE);
   phottreecicnoeleveto->SetPFJetName(jetCorr->GetOutputName());
   phottreecicnoeleveto->SetApplyElectronVeto(kFALSE);
   phottreecicnoeleveto->SetExcludeDoublePrompt(excludedoubleprompt);
   phottreecicnoeleveto->SetIsData(isData);  
   if (is25) phottreecicnoeleveto->SetEnablePFPhotons(kFALSE);
   phottreecicnoeleveto->SetApplyLeptonTag(kTRUE);
   phottreecicnoeleveto->SetLeptonTagElectronsName("HggLeptonTagElectrons");
   phottreecicnoeleveto->SetLeptonTagMuonsName("HggLeptonTagMuons");
   phottreecicnoeleveto->SetApplyVBFTag(kTRUE);
   phottreecicnoeleveto->SetApplyPFMetCorr(kTRUE);
   phottreecicnoeleveto->SetBeamspotWidth(5.0);
   
   PhotonTreeWriter *phottreepresel = new PhotonTreeWriter("PhotonTreeWriterPresel");
   phottreepresel->SetPhotonsFromBranch(kFALSE);
   phottreepresel->SetInputPhotonsName(photpresel->GetOutputName());
   phottreepresel->SetEnableJets(kTRUE);
   phottreepresel->SetApplyJetId(kTRUE);    
   phottreepresel->SetPFJetsFromBranch(kFALSE);
   phottreepresel->SetPFJetName(jetCorr->GetOutputName());  
   phottreepresel->SetExcludeDoublePrompt(excludedoubleprompt);  
   phottreepresel->SetIsData(isData);  
   if (is25) phottreepresel->SetEnablePFPhotons(kFALSE);
   phottreepresel->SetApplyLeptonTag(kTRUE);
   phottreepresel->SetLeptonTagElectronsName("HggLeptonTagElectrons");
   phottreepresel->SetLeptonTagMuonsName("HggLeptonTagMuons");
   phottreepresel->SetApplyVBFTag(kTRUE);
   phottreepresel->SetApplyPFMetCorr(kTRUE);
   phottreepresel->SetBeamspotWidth(5.0);
   
  phottreepresel -> SetElectronMVAWeightsSubdet0Pt10To20(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/Electrons_BDTG_NonTrigV0_Cat1.weights.xml"))));
  phottreepresel -> SetElectronMVAWeightsSubdet1Pt10To20(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/Electrons_BDTG_NonTrigV0_Cat2.weights.xml"))));  
  phottreepresel -> SetElectronMVAWeightsSubdet2Pt10To20(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/Electrons_BDTG_NonTrigV0_Cat3.weights.xml"))));  
  phottreepresel -> SetElectronMVAWeightsSubdet0Pt20ToInf(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/Electrons_BDTG_NonTrigV0_Cat4.weights.xml")))); 
  phottreepresel -> SetElectronMVAWeightsSubdet1Pt20ToInf(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/Electrons_BDTG_NonTrigV0_Cat5.weights.xml")))); 
  phottreepresel -> SetElectronMVAWeightsSubdet2Pt20ToInf(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/Electrons_BDTG_NonTrigV0_Cat6.weights.xml")))); 



   PhotonTreeWriter *phottreepreselinverteleveto = new PhotonTreeWriter("PhotonTreeWriterPreselInvertEleVeto");
   phottreepreselinverteleveto->SetPhotonsFromBranch(kFALSE);
   phottreepreselinverteleveto->SetInputPhotonsName(photpreselinverteleveto->GetOutputName());
   phottreepreselinverteleveto->SetEnableJets(kTRUE);
   phottreepreselinverteleveto->SetApplyJetId(kTRUE);    
   phottreepreselinverteleveto->SetPFJetsFromBranch(kFALSE);
   phottreepreselinverteleveto->SetPFJetName(jetCorr->GetOutputName());  
   phottreepreselinverteleveto->SetApplyElectronVeto(kFALSE);  
   phottreepreselinverteleveto->SetExcludeDoublePrompt(excludedoubleprompt);    
   phottreepreselinverteleveto->SetIsData(isData); 
   if (is25) phottreepreselinverteleveto->SetEnablePFPhotons(kFALSE);  
   phottreepreselinverteleveto->SetApplyLeptonTag(kTRUE);
   phottreepreselinverteleveto->SetLeptonTagElectronsName("HggLeptonTagElectrons");
   phottreepreselinverteleveto->SetLeptonTagMuonsName("HggLeptonTagMuons");
   phottreepreselinverteleveto->SetApplyVBFTag(kTRUE);
   phottreepreselinverteleveto->SetApplyPFMetCorr(kTRUE);
   phottreepreselinverteleveto->SetBeamspotWidth(5.0);
  
   PhotonTreeWriter *phottreepreselnosmear = new PhotonTreeWriter("PhotonTreeWriterPreselNoSmear");
   phottreepreselnosmear->SetPhotonsFromBranch(kFALSE);
   phottreepreselnosmear->SetInputPhotonsName(photpreselnosmear->GetOutputName());
   phottreepreselnosmear->SetEnableJets(kTRUE);
   phottreepreselnosmear->SetApplyJetId(kTRUE);  
   phottreepreselnosmear->SetPFJetsFromBranch(kFALSE);
   phottreepreselnosmear->SetPFJetName(jetCorr->GetOutputName());  
   phottreepreselnosmear->SetExcludeDoublePrompt(excludedoubleprompt);  
   phottreepreselnosmear->SetIsData(isData);    
   if (is25) phottreepreselnosmear->SetEnablePFPhotons(kFALSE); 
   phottreepreselnosmear->SetApplyLeptonTag(kTRUE);
   phottreepreselnosmear->SetLeptonTagElectronsName("HggLeptonTagElectrons");
   phottreepreselnosmear->SetLeptonTagMuonsName("HggLeptonTagMuons");
   phottreepreselnosmear->SetApplyVBFTag(kTRUE);
   phottreepreselnosmear->SetApplyPFMetCorr(kTRUE);
   phottreepreselnosmear->SetBeamspotWidth(5.0); 
   
   PhotonIDMod         *photidpresel = new PhotonIDMod("PhotonIDModPresel");
   photidpresel->SetPtMin(25.0);
   photidpresel->SetOutputName("GoodPhotonsPreselid");
   photidpresel->SetIDType("MITPFSelection");
   photidpresel->SetApplyElectronVeto(kTRUE);
   photidpresel->SetIsData(isData);
   
  PhotonIDMod         *photidpreselinvert = new PhotonIDMod("PhotonIDModPreselInvert");
  photidpreselinvert->SetPtMin(25.0);
  photidpreselinvert->SetOutputName("GoodPhotonsPreselidInvert");
  photidpreselinvert->SetIDType("MITPFSelection");
  photidpreselinvert->SetApplyElectronVeto(kFALSE);
  photidpreselinvert->SetInvertElectronVeto(kTRUE);
  photidpreselinvert->SetIsData(isData);  
  
  PhotonTreeWriter *phottreesingle = new PhotonTreeWriter("PhotonTreeWriterSingle");
  phottreesingle->SetWriteDiphotonTree(kFALSE);
  phottreesingle->SetPhotonsFromBranch(kFALSE);
  phottreesingle->SetInputPhotonsName(photidpresel->GetOutputName());  
  phottreesingle->SetEnableJets(kTRUE);
  phottreesingle->SetPFJetsFromBranch(kFALSE);
  phottreesingle->SetPFJetName(jetCorr->GetOutputName());      
  phottreesingle->SetBeamspotWidth(5.0);
  phottreesingle->SetIsData(isData);
  if (is25) phottreesingle->SetEnablePFPhotons(kFALSE);  
  
  PhotonTreeWriter *phottreesingleinvert = new PhotonTreeWriter("PhotonTreeWriterSingleInvert");
  phottreesingleinvert->SetWriteDiphotonTree(kTRUE);
  phottreesingleinvert->SetPhotonsFromBranch(kFALSE);
  phottreesingleinvert->SetInputPhotonsName(photidpreselinvert->GetOutputName());  
  phottreesingleinvert->SetEnableJets(kTRUE);
  phottreesingleinvert->SetPFJetsFromBranch(kFALSE);
  phottreesingleinvert->SetPFJetName(pubJetOpen->GetOutputName());
  phottreesingleinvert->SetApplyElectronVeto(kFALSE);
  phottreesingleinvert->SetBeamspotWidth(5.0);
  phottreesingleinvert->SetIsData(isData);  
  if (is25) phottreesingleinvert->SetEnablePFPhotons(kFALSE);    
  

  //------------------------------------------------------------------------------------------------
  // making analysis chain
  //------------------------------------------------------------------------------------------------
  // this is how it always starts
  runLumiSel      ->Add(mcselmod);
    
  if (TString(dataset).Contains("-h")) {    
    mcselmod        ->Add(sysMod);
  }
  
  // high level trigger is always first
  //mcselmod         ->Add(hltModES);
  //mcselmod         ->Add(hltModE);
  
  if (!TString(dataset).Contains("meridiani")) {      
    mcselmod         ->Add(hltModP);
    hltModP         ->Add(goodPVFilterMod);
  }
  else {
    mcselmod->Add(goodPVFilterMod);
  }
  //hltModE         ->Add(goodPVFilterModE);
  //hltModES        ->Add(goodPVFilterModES);
  
  //goodPVFilterMod ->Add(muonId);
  goodPVFilterMod->Add(photreg);
  photreg->Add(pubJet);
  pubJet->Add(jetCorr);
  
  // simple object id modules
  //goodPVFilterModE -> Add(pubJetOpen);
  //pubJetOpen       -> Add(photidpreselinvert);
  //goodPVFilterModES ->Add(elecIdS);
  jetCorr          ->Add(SepPUMod); 
  SepPUMod         ->Add(eleIdMod);
  eleIdMod         ->Add(muonIdMod);
   
//   muonIdMod          ->Add(photcic);
//   muonIdMod          ->Add(photcicnoeleveto);  
//   muonIdMod          ->Add(photpreselinverteleveto);  
//   muonIdMod          ->Add(photpreselnosmear);  

  muonIdMod          ->Add(photpresel);  


//   photcic         ->Add(phottreecic);
//   photcicnoeleveto       ->Add(phottreecicnoeleveto);
//   photpreselinverteleveto    ->Add(phottreepreselinverteleveto);
//   photpreselnosmear    ->Add(phottreepreselnosmear);

  photpresel    ->Add(phottreepresel);


  //jetCorr          ->Add(photidpresel);
  //photidpresel       ->Add(phottreesingle);
  
  //photidpreselinvert       ->Add(phottreesingleinvert);  
  //elecIdS->Add(phottreeES);
  

  //TFile::SetOpenTimeout(0);
  //TFile::SetCacheFileDir("./rootfilecache",kTRUE,kTRUE);
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

  if (TString(dataset).Contains("meridiani")) {    
    ana->SetUseHLT(kFALSE);
  }
  
  //------------------------------------------------------------------------------------------------
  // organize input
  //------------------------------------------------------------------------------------------------
  Catalog *c = new Catalog(catalogDir);
  TString skimdataset = TString(dataset)+TString("/") +TString(skim);
  Dataset *d = NULL;
  TString bookstr = book;
  //if (TString(dataset).Contains("s11-h")) bookstr.ReplaceAll("local","t2mit");

//   if (TString(skim).CompareTo("noskim") == 0)
//     d = c->FindDataset(bookstr,dataset,fileset);
//   else 
//     d = c->FindDataset(bookstr,skimdataset.Data(),fileset);
//   ana->AddDataset(d);

  std::vector<TString> synchFile;

//   synchFile.push_back("26650E45-A902-E211-9523-001D09F29146.root");
//   synchFile.push_back("FC4B4004-6502-E211-B532-BCAEC518FF8A.root");
//   synchFile.push_back("063C9BE3-8F02-E211-9503-001D09F2424A.root");
//   synchFile.push_back("10946F3E-4F02-E211-ACA8-5404A63886BB.root");
//   synchFile.push_back("70EAB2D6-8F02-E211-AD24-001D09F29146.root");
//   synchFile.push_back("50F136D5-8F02-E211-A7E1-001D09F25109.root");
//   synchFile.push_back("ACF83417-6502-E211-8A06-002481E0D90C.root");
//   synchFile.push_back("50F136D5-8F02-E211-A7E1-001D09F25109.root");
//   synchFile.push_back("84D8683B-3702-E211-ACC8-BCAEC53296F3.root");
//   synchFile.push_back("64CF7800-6802-E211-AE29-003048D374CA.root");
//   synchFile.push_back("044DE8D2-8F02-E211-BEF1-00215AEDFD98.root");
//   synchFile.push_back("6A50C106-6002-E211-A8FD-003048F118E0.root");
//   synchFile.push_back("80DCA06C-5002-E211-B0FD-5404A6388692.root");
//   synchFile.push_back("1600BC48-4F02-E211-B8FA-BCAEC518FF6B.root");
//   synchFile.push_back("36E6F661-9002-E211-8ED4-001D09F253D4.root");

  synchFile.push_back("28B1B2DB-AAF5-E111-8D01-E0CB4E4408E3.root");
  synchFile.push_back("1EA891A9-77F8-E111-BD59-003048D2BB90.root");
  synchFile.push_back("38E9A015-F1F8-E111-9654-001D09F28D4A.root");
  synchFile.push_back("DC46A78E-60E0-E111-B65D-001D09F25479.root");
  synchFile.push_back("B490A799-45DD-E111-A676-001D09F2AD4D.root");
  synchFile.push_back("2E00B5D2-C9DD-E111-81B5-0025901D6268.root");

  for( unsigned int iF = 0 ; iF<synchFile.size(); ++iF) {
    ana->AddFile("root://xrootd.cmsaf.mit.edu//store/user/paus/filefi/029/r12c-dph-pr-v2/"+synchFile[iF]);   ///scratch/fabstoec/root/DA2D6D6F-BB95-E111-9564-001A92971B7E.root");
    //ana->AddFile("root://xrootd.cmsaf.mit.edu//store/user/paus/filefi/029/r12b-dph-j13-v1/"+synchFile[iF]);   ///scratch/fabstoec/root/DA2D6D6F-BB95-E111-9564-001A92971B7E.root");
  }
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

