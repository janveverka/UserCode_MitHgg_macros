// $Id: runHgg.C,v 1.4 2011/12/19 21:56:23 bendavid Exp $
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
void runHgg20125fb(const char *fileset    = "",
	    const char *skim       = "noskim",
	    //const char *dataset    = "w10-h120gg-gf-v8-pu",
           // const char *dataset    = "w10-qcd2em40-v8-pu",
            //const char *dataset    = "s11-qcd2em40-v11-pu",
            //const char *dataset    = "r11a-dph-j05-v1",          
	    //const char *dataset    = "s11-h120gg-gf-v11-pu",
	    //const char *dataset    = "s11-pj-2em20-v11-pu",
            //const char *dataset      = "s11-h130gg-gf-v11-pu",
            //const char *dataset    = "f11-pj-2em20-v14b-pu",
            //const char *dataset    = "f11--qcd-2em40-v14b-pu",
            //const char *dataset    = "f11--h120gg-gf-v14b-pu",
            //const char *dataset = "f11--h120gg-vbf-v14b-pu",
            const char *dataset    = "f11--2pibx10_25-v14b-pu",
            //const char *dataset = "s11-zeem20-powheg-v11-pu",
            //const char *dataset = "s11-wjets-v11-pu",
            //const char *dataset = "p11-zll50-v1g1-pu",
	    const char *book       = "t2mit/filefi/025",
	    const char *catalogDir = "/home/cmsprod/catalog",
	    const char *outputName = "hgg",
	    int         nEvents    = -1)
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

  TString jsonFile = TString("/home/mtouch/cms/json/") + TString(json);
  // TString jsonFile = TString("/home/mtouch/cms/json/") + TString("augmented_nov8_rereco.json");
  //TString jsonFile = TString("/home/mtouch/cms/json/augmented_nov08_rereco.json");
  printf(" json file is %s\n",json);
  Bool_t  isData   = ( (jsonFile.CompareTo("/home/mtouch/cms/json/~") != 0) );
  //  printf(" json didnt fail lol!\n");
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
  //  printf(" json didnt fail %s\n",jsonFile.Data());
  // only select on run- and lumisection numbers when valid json file present
  if ((jsonFile.CompareTo("/home/mtouch/cms/json/~") != 0) &&
      (jsonFile.CompareTo("/home/mtouch/cms/json/-") != 0)   ) {
    printf(" json didnt fail %s\n",jsonFile.Data());
    runLumiSel->AddJSONFile(jsonFile.Data());
    //runLumiSel->AddJSONFile("/home/mtouch/cms/json/augmented_nov08_rereco.json");
  }
  //    printf(" json didnt fail lol12!\n");
  if ((jsonFile.CompareTo("/home/mtouch/cms/json/-") == 0)   ) {
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
//   hltModE->AddTrigger("HLT_Photon10_L1R",132440,137028);
//   hltModE->AddTrigger("HLT_Photon15_Cleaned_L1R",138564,140401);
//   hltModE->AddTrigger("HLT_Ele15_SW_CaloEleId_L1R",141956,144114);
//   hltModE->AddTrigger("HLT_Ele17_SW_CaloEleId_L1R",144115,147145);
//   hltModE->AddTrigger("HLT_Ele17_SW_TightEleId_L1R",147146,148102);
//   hltModE->AddTrigger("HLT_Ele27_SW_TightCaloEleIdTrack_L1R_v1",147146,148102);  
//   hltModE->AddTrigger("HLT_Ele22_SW_TighterCaloIdIsol_L1R_v1",148103,159999);  
//   hltModE->AddTrigger("HLT_Ele22_SW_TighterCaloIdIsol_L1R_v2",148103,159999);  
//   hltModE->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v1",160000,999999);
//   hltModE->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v2",160000,999999);
//   hltModE->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v3",160000,999999);
//   hltModE->AddTrigger("HLT_Ele15_LW_L1R",1,1);
//   hltModE->AddTrigger("HLT_Ele17_SW_CaloEleId_L1R",1,1);
//   hltModE->AddTrigger("HLT_Ele22_SW_TighterCaloIdIsol_L1R_v2",1,1);


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
  hltModE->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v10",179890,999999);
  hltModE->SetTrigObjsName("MyHltElecObjs");
  hltModE->SetAbortIfNotAccepted(isData);

  HLTMod *hltModES = new HLTMod("HLTModES");
  hltModES->AddTrigger("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v*",160000,999999);
  hltModES->AddTrigger("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2",160000,999999);
  hltModES->AddTrigger("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v3",160000,999999);
  hltModES->AddTrigger("HLT_Ele25_WP80_PFMT40_v*",160000,999999);
  hltModES->AddTrigger("HLT_Ele27_WP80_PFMT50_v*",160000,999999);
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

  SeparatePileUpMod* SepPUModDVtx = new SeparatePileUpMod;
  //  SepPUMod->SetUseAllVerteces(kFALSE);
  //  SepPUModDVtx->SetVertexName("OutVtxCiCDVtx");
  SepPUModDVtx->SetPFNoPileUpName("pfnopileupcandsDVtx");
  SepPUModDVtx->SetPFPileUpName("pfpileupcandsDVtx");
  SepPUModDVtx->SetCheckClosestZVertex(kFALSE);

  //-----------------------------------
  // Lepton Selection    ****from Fabian's code MitHgg/macro/runHgg_LT.C
  //-----------------------------------
  ElectronIDMod* eleIdMod = new ElectronIDMod;
  eleIdMod -> SetEtMin(20.);
  eleIdMod -> SetEtaMax(2.5);
  eleIdMod -> SetApplyEcalFiducial(true);
  eleIdMod -> SetIDType("Hgg_LeptonTag_2012Id");  //was changed in the electron tool mod by Heng ** June 16 
  eleIdMod -> SetVertexName("OutVtxCiC");
  eleIdMod -> SetWhichVertex(0);
  eleIdMod -> SetD0Cut(0.02);
  eleIdMod -> SetDZCut(0.2); //h  
  eleIdMod -> SetIsoType("PFIso_HggLeptonTag2012"); //h
  eleIdMod -> SetOutputName("HggLeptonTagElectrons");
  eleIdMod -> SetRhoType(RhoUtilities::CMS_RHO_RHOKT6PFJETS);
  eleIdMod -> SetPFNoPileUpName("pfnopileupcands");

  MuonIDMod* muonIdMod = new MuonIDMod;
  // base kinematics
  muonIdMod -> SetPtMin(20.);
  muonIdMod -> SetEtaCut(2.4);
  // base ID
  muonIdMod -> SetIDType("Tight");
  muonIdMod -> SetVertexName("OutVtxCiC");
  muonIdMod -> SetWhichVertex(0); // this is a 'hack'.. but hopefully good enough...
  muonIdMod -> SetD0Cut(0.2);
  muonIdMod -> SetDZCut(0.5);
  muonIdMod -> SetIsoType("PFIsoBetaPUCorrected"); //h
  muonIdMod ->SetPFIsoCut(0.2); //h
  muonIdMod -> SetOutputName("HggLeptonTagMuons");
  muonIdMod -> SetPFNoPileUpName("pfnopileupcands");
  muonIdMod -> SetPFPileUpName("pfpileupcands");

  //****** new modules for default vertex

   ElectronIDMod* eleIdModDVtx = new ElectronIDMod;
   eleIdModDVtx -> SetEtMin(20.);
   eleIdModDVtx -> SetEtaMax(2.5);
   eleIdModDVtx -> SetApplyEcalFiducial(true);
   eleIdModDVtx -> SetIDType("Hgg_LeptonTag_2012Id");  //was changed in the electron tool mod by Heng ** June 16 
   eleIdModDVtx -> SetVertexName("OutVtxCiCDVtx");
   eleIdModDVtx -> SetWhichVertex(0);
   eleIdModDVtx -> SetD0Cut(0.02);
   eleIdModDVtx -> SetDZCut(0.2); //h  
   eleIdModDVtx -> SetIsoType("PFIso_HggLeptonTag2012"); //h
   eleIdModDVtx -> SetOutputName("HggLeptonTagElectronsDVtx");
   eleIdModDVtx -> SetRhoType(RhoUtilities::CMS_RHO_RHOKT6PFJETS);
   eleIdModDVtx -> SetPFNoPileUpName("pfnopileupcandsDVtx");

   MuonIDMod* muonIdModDVtx = new MuonIDMod;
   //base kinematics
   muonIdModDVtx -> SetPtMin(20.);
   muonIdModDVtx -> SetEtaCut(2.4);
   muonIdModDVtx -> SetIDType("Tight");
   muonIdModDVtx -> SetVertexName("OutVtxCiCDVtx");
   muonIdModDVtx -> SetWhichVertex(0);  //this is a 'hack'.. but hopefully good enough...
   muonIdModDVtx -> SetD0Cut(0.2);
   muonIdModDVtx -> SetDZCut(0.5);
   muonIdModDVtx -> SetIsoType("PFIsoBetaPUCorrected"); //h
   muonIdModDVtx -> SetPFIsoCut(0.2); //h
   muonIdModDVtx -> SetOutputName("HggLeptonTagMuonsDVtx");
   muonIdModDVtx -> SetPFNoPileUpName("pfnopileupcandsDVtx");
   muonIdModDVtx -> SetPFPileUpName("pfpileupcandsDVtx");



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
  
  PublisherMod<PFJet,Jet> *pubJetOpen = new PublisherMod<PFJet,Jet>("JetPubOpen");
  pubJetOpen->SetInputName("AKt5PFJets");
  pubJetOpen->SetOutputName("PubAKt5PFJetsOpen");  

  JetCorrectionMod *jetCorr = new JetCorrectionMod;
  if(isData){ 
    jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer12_V1_DATA_L1FastJet_AK5PF.txt")).Data())); 
    jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer12_V1_DATA_L2Relative_AK5PF.txt")).Data())); 
    jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer12_V1_DATA_L3Absolute_AK5PF.txt")).Data())); 
    jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer12_V1_DATA_L2L3Residual_AK5PF.txt")).Data())); 
  }
  else {
    jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer12_V1_MC_L1FastJet_AK5PF.txt")).Data())); 
    jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer12_V1_MC_L2Relative_AK5PF.txt")).Data())); 
    jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer12_V1_MC_L3Absolute_AK5PF.txt")).Data())); 
  }



 //  jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/START42_V12_AK5PF_L1FastJet.txt")).Data())); 
//   jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/START42_V12_AK5PF_L2Relative.txt")).Data())); 
//   jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/START42_V12_AK5PF_L3Absolute.txt")).Data())); 
//   if(isData){ 
//     jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/START42_V12_AK5PF_L2L3Residual.txt")).Data())); 
//   }
  jetCorr->SetInputName(pubJet->GetOutputName());
  jetCorr->SetCorrectedName("CorrectedJets");

  //  jetCorr->SetRhoType(RhoUtilities::MIT_RHO_RANDOM_HIGH_ETA);

  Bool_t excludedoubleprompt = kFALSE;
  if (TString(dataset).Contains("-pj")) {
    mcselmod->ExcludeProcess(18);
    mcselmod->ExcludeProcess(114);
    excludedoubleprompt = kTRUE;
    //excludedoubleprompt = kFALSE;
  }

  if (TString(dataset).Contains("-qcd2em") || TString(dataset).Contains("-qcd-2em")) {
    excludedoubleprompt = kTRUE;
    //excludedoubleprompt = kFALSE;
  }
  
  Bool_t is25 = kFALSE;
  if (TString(book).Contains("025")) is25 = kTRUE;


  PhotonMvaMod *photreg = new PhotonMvaMod;
  photreg->SetRegressionVersion(3);
  photreg->SetRegressionWeights(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/gbrv3ph_52x.root")).Data()));
  photreg->SetOutputName("GoodPhotonsRegr");
  photreg->SetIsData(isData);

    
  PhotonPairSelector         *photcic = new PhotonPairSelector("PhotonPairSelectorCiC");
  photcic->SetOutputName("GoodPhotonsCIC");
  photcic->SetOutputVtxName("OutVtxCiC");
  //  photcic->SetOutputVtxName("CiCChosenVtx");
  photcic->SetPhotonSelType("CiCPFSelection");
  //  photcic->SetPhotonSelType("CiCSelection");
  photcic->SetVertexSelType("CiCMVA2012Selection");
  //photcic->SetVertexSelType("ZeroVtxSelection");
  photcic->DoMCSmear(kTRUE);  // for synchronization only
  photcic->DoDataEneCorr(kTRUE);  //for synchronization only
  photcic->SetPhotonsFromBranch(kFALSE);
  photcic->SetInputPhotonsName(photreg->GetOutputName());
//   photcic->SetMCSmearFactors(0.0089, 0.0089, 0.0109, 0.0156, 0.0203,0.0303,0.0326,0.0318,0.0331);
//   photcic->AddEnCorrPerRun(160000,167913,0.9905,0.9905,0.9971,0.9976,1.0094,0.9994,1.0044,0.9968,1.0079);
//   photcic->AddEnCorrPerRun(170249,172619,0.9909,0.9909,0.9975,0.9994,1.0112,0.9962,1.0012,0.9962,1.0072);
//   photcic->AddEnCorrPerRun(172620,173692,0.9909,0.9909,0.9975,0.9977,1.0096,0.9963,1.0013,0.9947,1.0057);
//   photcic->AddEnCorrPerRun(175860,177139,0.9911,0.9911,0.9977,0.9990,1.0109,0.9922,0.9973,0.9967,1.0077);
//   photcic->AddEnCorrPerRun(177140,178421,0.9910,0.9910,0.9975,0.9987,1.0105,0.9921,0.9972,0.9975,1.0085);
//   photcic->AddEnCorrPerRun(178424,999999,0.9903,0.9903,0.9969,0.9976,1.0095,0.9889,0.9940,0.9976,1.0086);
//   photcic->SetMCSmearFactors(0.0067,0.0077,0.0096,0.0141,0.0196,0.0268,0.0279,0.0293,0.0301);//ming:smear(sigE/E)
//   photcic->AddEnCorrPerRun(160431,167913,0.9941,0.9941,1.0004,0.9916,1.0045,1.0033,1.0082,0.9958,1.0064);//ming:Emc/Edata
//   photcic->AddEnCorrPerRun(170000,172619,0.9954,0.9954,1.0016,0.9937,1.0066,0.9976,1.0025,0.9940,1.0046);
//   photcic->AddEnCorrPerRun(172620,173692,0.9955,0.9955,1.0017,0.9929,1.0058,0.9986,1.0035,0.9923,1.0029);
//   photcic->AddEnCorrPerRun(175830,177139,0.9958,0.9958,1.0021,0.9944,1.0073,0.9968,1.0017,0.9933,1.004);
//   photcic->AddEnCorrPerRun(177140,178421,0.9962,0.9962,1.0025,0.9946,1.0075,0.9960,1.0010,0.9944,1.005);
//   photcic->AddEnCorrPerRun(178424,180252,0.9961,0.9961,1.0024,0.9942,1.0071,0.9921,0.9970,0.9953,1.0059); 

  //2012 ICHEP promptreco 5.3/fb
  photcic->SetMCSmearFactors(0.0106,0.0108,0.0119,0.0149,0.0240,0.0375,0.0330,0.0607,0.0602);
  photcic->AddEnCorrPerRun(190450, 190781, 0.9946, 0.9946, 1.0017, 0.9976, 1.0124, 0.9932, 0.9986, 0.9859, 1.0028);
  photcic->AddEnCorrPerRun(190782, 190949, 1.0079, 1.0079, 1.0149, 0.9850, 1.0000, 1.0056, 1.0110, 0.9932, 1.0099);
  photcic->AddEnCorrPerRun(190950, 191833, 0.9974, 0.9974, 1.0044, 0.9989, 1.0137, 0.9997, 1.0051, 0.9736, 0.9907);
  photcic->AddEnCorrPerRun(191834, 193686, 0.9943, 0.9943, 1.0013, 0.9961, 1.0109, 0.9988, 1.0043, 0.9726, 0.9897);
  photcic->AddEnCorrPerRun(193746, 194210, 0.9950, 0.9950, 1.0020, 0.9948, 1.0096, 1.0004, 1.0059, 0.9875, 1.0043);
  photcic->AddEnCorrPerRun(194211, 194479, 0.9964, 0.9964, 1.0034, 0.9968, 1.0116, 0.9977, 1.0032, 0.9886, 1.0054);
  photcic->AddEnCorrPerRun(194480, 195147, 0.9974, 0.9974, 1.0044, 0.9990, 1.0138, 0.9986, 1.0041, 0.9905, 1.0073);
  photcic->AddEnCorrPerRun(195148, 195350, 0.9979, 0.9979, 1.0050, 0.9992, 1.0139, 1.0018, 1.0072, 0.9937, 1.0105);
  photcic->AddEnCorrPerRun(195378, 195530, 0.9973, 0.9973, 1.0044, 0.9985, 1.0133, 0.9954, 1.0008, 0.9856, 1.0025);
  photcic->AddEnCorrPerRun(195531, 196529, 0.9906, 0.9906, 0.9977, 0.9921, 1.0069, 1.0003, 1.0058, 1.0328, 1.0488);
  photcic->SetDoShowerShapeScaling(kTRUE); 
  photcic->SetShowerShapeType("2012ShowerShape");
  //  photcic->SetDoMCErrScaling(kTRUE);
  //  photcic->SetMCErrScale(1, 1); //ming:scale(sigE/E)
  photcic->SetJetsName(jetCorr->GetOutputName());    
  photcic->SetIsData(isData);


   PhotonPairSelector         *photcicDVtx = new PhotonPairSelector("PhotonPairSelectorCiCDefaultVtx");
   photcicDVtx->SetOutputName("GoodPhotonsCICDVtx");
   photcicDVtx->SetOutputVtxName("OutVtxCiCDVtx");
   //  photcicDVtx->SetOutputVtxName("CiCChosenVtx");
   photcicDVtx->SetPhotonSelType("CiCPFSelection");
   //  photcicDVtx->SetPhotonSelType("CiCSelection");
   photcicDVtx->SetVertexSelType("ZeroVtxSelection");
   photcicDVtx->DoMCSmear(kTRUE);// for synchronization only
   photcicDVtx->DoDataEneCorr(kTRUE);//for synchronization only
   photcicDVtx->SetPhotonsFromBranch(kFALSE);
   photcicDVtx->SetInputPhotonsName(photreg->GetOutputName());
  //2012 ICHEP promptreco 5.3/fb
   photcicDVtx->SetMCSmearFactors(0.0106,0.0108,0.0119,0.0149,0.0240,0.0375,0.0330,0.0607,0.0602);
   photcicDVtx->AddEnCorrPerRun(190450, 190781, 0.9946, 0.9946, 1.0017, 0.9976, 1.0124, 0.9932, 0.9986, 0.9859, 1.0028);
   photcicDVtx->AddEnCorrPerRun(190782, 190949, 1.0079, 1.0079, 1.0149, 0.9850, 1.0000, 1.0056, 1.0110, 0.9932, 1.0099);
   photcicDVtx->AddEnCorrPerRun(190950, 191833, 0.9974, 0.9974, 1.0044, 0.9989, 1.0137, 0.9997, 1.0051, 0.9736, 0.9907);
   photcicDVtx->AddEnCorrPerRun(191834, 193686, 0.9943, 0.9943, 1.0013, 0.9961, 1.0109, 0.9988, 1.0043, 0.9726, 0.9897);
   photcicDVtx->AddEnCorrPerRun(193746, 194210, 0.9950, 0.9950, 1.0020, 0.9948, 1.0096, 1.0004, 1.0059, 0.9875, 1.0043);
   photcicDVtx->AddEnCorrPerRun(194211, 194479, 0.9964, 0.9964, 1.0034, 0.9968, 1.0116, 0.9977, 1.0032, 0.9886, 1.0054);
   photcicDVtx->AddEnCorrPerRun(194480, 195147, 0.9974, 0.9974, 1.0044, 0.9990, 1.0138, 0.9986, 1.0041, 0.9905, 1.0073);
   photcicDVtx->AddEnCorrPerRun(195148, 195350, 0.9979, 0.9979, 1.0050, 0.9992, 1.0139, 1.0018, 1.0072, 0.9937, 1.0105);
   photcicDVtx->AddEnCorrPerRun(195378, 195530, 0.9973, 0.9973, 1.0044, 0.9985, 1.0133, 0.9954, 1.0008, 0.9856, 1.0025);
   photcicDVtx->AddEnCorrPerRun(195531, 196529, 0.9906, 0.9906, 0.9977, 0.9921, 1.0069, 1.0003, 1.0058, 1.0328, 1.0488); 
   photcicDVtx->SetDoShowerShapeScaling(kTRUE); 
   photcicDVtx->SetShowerShapeType("2012ShowerShape");
   //   photcicDVtx->SetDoMCErrScaling(kTRUE);
   //   photcicDVtx->SetMCErrScale(1, 1); //ming:scale(sigE/E)
   photcicDVtx->SetJetsName(jetCorr->GetOutputName());    
   photcicDVtx->SetIsData(isData);

  
  
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
   phottreecic->SetApplyLeptonTag(true);
   phottreecic->SetLeptonTagElectronsName("HggLeptonTagElectrons");
   phottreecic->SetLeptonTagMuonsName    ("HggLeptonTagMuons");
   phottreecic->SetApplyPFMetCorr(true);
   //  phottreecic -> SetPFNoPileUpName("pfnopileupcands");
   phottreecic-> SetPFNoPileUpName("pfnopileupcands");
   phottreecic-> SetPFPileUpName("pfpileupcands");
   phottreecic->SetDo2012LepTag(kTRUE);


   PhotonTreeWriter *phottreecicDVtx = new PhotonTreeWriter("PhotonTreeWriterCiC");
   phottreecicDVtx->SetPhotonsFromBranch(kFALSE);
   phottreecicDVtx->SetInputPhotonsName(photcicDVtx->GetOutputName());
   phottreecicDVtx->SetEnableJets(kTRUE);
   phottreecicDVtx->SetApplyJetId(kTRUE);
   phottreecicDVtx->SetPFJetsFromBranch(kFALSE);
   phottreecicDVtx->SetPFJetName(jetCorr->GetOutputName());
   phottreecicDVtx->SetExcludeDoublePrompt(excludedoubleprompt);
   phottreecicDVtx->SetIsData(isData);
   if (is25) phottreecicDVtx->SetEnablePFPhotons(kFALSE);
   phottreecicDVtx->SetApplyLeptonTag(true);
   phottreecicDVtx->SetLeptonTagElectronsName("HggLeptonTagElectronsDVtx");
   phottreecicDVtx->SetLeptonTagMuonsName    ("HggLeptonTagMuonsDVtx");
   phottreecicDVtx->SetApplyPFMetCorr(true);
   //   phottreecicDVtx -> SetPFNoPileUpName("pfnopileupcandsDVtx");
   phottreecicDVtx -> SetPFNoPileUpName("pfnopileupcandsDVtx");
   phottreecicDVtx -> SetPFPileUpName("pfpileupcandsDVtx");
   phottreecicDVtx->SetDo2012LepTag(kTRUE);



 //   PhotonTreeWriter *phottreecicDVtx = new PhotonTreeWriter("PhotonTreeWriterCiCDVtx");
//    phottreecic->SetPhotonsFromBranch(kFALSE);
//    phottreecic->SetInputPhotonsName(photcicDVtx->GetOutputName());
//    phottreecic->SetEnableJets(kTRUE);
//    phottreecic->SetApplyJetId(kTRUE);
//    phottreecic->SetPFJetsFromBranch(kFALSE);
//    phottreecic->SetPFJetName(jetCorr->GetOutputName());
//    phottreecic->SetExcludeDoublePrompt(excludedoubleprompt);
//    phottreecic->SetIsData(isData);
//    if (is25) phottreecic->SetEnablePFPhotons(kFALSE);
//    phottreecic->SetApplyLeptonTag(true);
//    phottreecic->SetLeptonTagElectronsName("HggLeptonTagElectrons");
//    phottreecic->SetLeptonTagMuonsName    ("HggLeptonTagMuons");
//    phottreecic->SetApplyPFMetCorr(true);



  PhotonIDMod         *photidcic = new PhotonIDMod("PhotonIDModPresel");
  photidcic->SetPtMin(25.0);
  photidcic->SetOutputName("GoodPhotonsPreselid");
  photidcic->SetOutputName("MITSelection");
  photidcic->SetApplyElectronVeto(kTRUE);
  photidcic->SetIsData(isData);

  PhotonTreeWriter *phottreesingle = new PhotonTreeWriter("PhotonTreeWriterSingle");
  phottreesingle->SetWriteDiphotonTree(kFALSE);
  phottreesingle->SetPhotonsFromBranch(kFALSE);
  phottreesingle->SetInputPhotonsName(photidcic->GetOutputName());  
  phottreesingle->SetIsData(isData);

  

  //------------------------------------------------------------------------------------------------
  // making analysis chain
  //------------------------------------------------------------------------------------------------
  // this is how it always starts
  runLumiSel      ->Add(mcselmod);
    
  if (TString(dataset).Contains("-h")) {    
    mcselmod        ->Add(sysMod);
  }
  
  // high level trigger is always first
  //mcselmod         ->Add(hltModE);
  //mcselmod         ->Add(hltModES);

  //hltModE         ->Add(goodPVFilterModE);
  //hltModES        ->Add(goodPVFilterModES);
  
  //goodPVFilterMod ->Add(muonId);


  // --------------------------------------------
  mcselmod         ->Add(hltModP);
  hltModP          ->Add(goodPVFilterMod);

  goodPVFilterMod  ->Add(photreg);
  photreg          ->Add(pubJet);
  pubJet           ->Add(jetCorr);

  //jetCorr->Add(eleIdMod);
  jetCorr          ->Add(photcic);
  jetCorr          ->Add(photcicDVtx);

  photcic          ->Add(SepPUMod);
  photcicDVtx      ->Add(SepPUModDVtx);

  SepPUMod         ->Add(muonIdMod);
  SepPUModDVtx     ->Add(muonIdModDVtx);

  muonIdMod         ->Add(eleIdMod);
  muonIdModDVtx     ->Add(eleIdModDVtx);

  eleIdMod       ->Add(phottreecic);   
  eleIdModDVtx   ->Add(phottreecicDVtx);
  // --------------------------------------------

  // simple object id modules
  //goodPVFilterModE ->Add(elecId);
  //goodPVFilterModES ->Add(elecIdS);

  




//   muonIdMod          ->Add(photpreselMVAVert);  
//   muonIdMod          ->Add(photpreselZeroVert);  
//   muonIdMod          ->Add(photpreselnosmear);  



//   photpreselMVAVert    ->Add(phottreepreselMVAVert);
//   photpreselZeroVert    ->Add(phottreepreselZeroVert);


//   jetCorr          ->Add(photidcic);
//   photidcic       ->Add(phottreesingle);
  
  //  elecId->Add(phottreeE);
  //elecIdS->Add(phottreeES);
 

  //TFile::SetOpenTimeout(0);
  //  TFile::SetCacheFileDir("./rootfilecache",kTRUE,kTRUE);
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
