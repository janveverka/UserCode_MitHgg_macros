// $Id: runHgg.C,v 1.2 2011/01/29 23:17:03 paus Exp $
#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TSystem.h>
#include "MitAna/DataUtil/interface/Debug.h"
#include "MitAna/Catalog/interface/Catalog.h"
#include "MitAna/TreeMod/interface/Analysis.h"
#include "MitAna/TreeMod/interface/HLTMod.h"
#include "MitAna/PhysicsMod/interface/RunLumiSelectionMod.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitPhysics/Mods/interface/GoodPVFilterMod.h"
#include "MitPhysics/Mods/interface/MuonIDMod.h"
#include "MitPhysics/Mods/interface/ElectronIDMod.h"
#include "MitPhysics/Mods/interface/ElectronCleaningMod.h"
#include "MitPhysics/Mods/interface/PhotonIDMod.h"
#include "MitPhysics/Mods/interface/PhotonCleaningMod.h"
#include "MitPhysics/Mods/interface/MergeLeptonsMod.h"
#include "MitHgg/Mods/interface/ZmmAnalysis.h"
#include "MitHgg/Mods/interface/ZeeAnalysis.h"
#include "MitHgg/Mods/interface/HggAnalysis.h"
#endif

//--------------------------------------------------------------------------------------------------
void runHgg(const char *fileset    = "0000",
	    const char *skim       = "noskim",
	    const char *dataset    = "w10-zz-z2-v8-pu",
	    //const char *dataset    = "r10b-pho-d22",
	    const char *book       = "t2mit/filefi/017",
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
    printf(" JSON file was not properly defined. EXIT!\n");
    return;
  } 
  TString jsonFile = TString("/home/cmsprod/json/") + TString(json);
  Bool_t  isData   = ( (jsonFile.CompareTo("/home/cmsprod/json/~") != 0) );

  if (gSystem->Getenv("MIT_PROD_OVERLAP")) {
    sprintf(overlap,"%s",gSystem->Getenv("MIT_PROD_OVERLAP"));
    if (EOF == sscanf(overlap,"%f",&overlapCut)) {
      printf(" Overlap was not properly defined. EXIT!\n");
      return;
    }
  }
  else {
    printf(" OVERLAP file was not properly defined. EXIT!\n");
    return;
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

  // only select on run- and lumisection numbers when valid json file present
  if ((jsonFile.CompareTo("/home/cmsprod/json/~") != 0) &&
      (jsonFile.CompareTo("/home/cmsprod/json/-") != 0)   ) {
    runLumiSel->AddJSONFile(jsonFile.Data());
  }
  if ((jsonFile.CompareTo("/home/cmsprod/json/-") == 0)   ) {
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
  hltModE->AddTrigger("HLT_Photon10_L1R");
  hltModE->AddTrigger("HLT_Ele10_LW_L1R");
  hltModE->AddTrigger("HLT_Ele15_SW_L1R");
  hltModE->AddTrigger("HLT_Ele15_SW_CaloEleId_L1R");
  hltModE->AddTrigger("HLT_Ele17_SW_CaloEleId_L1R");
  hltModE->AddTrigger("HLT_Ele17_SW_TightCaloEleId_SC8HE_L1R_v1");
  hltModE->AddTrigger("HLT_DoubleEle10_SW_L1R");
  hltModE->AddTrigger("HLT_Ele10_SW_L1R_v2");
  hltModE->AddTrigger("HLT_Ele12_SW_TighterEleId_L1R_v2");
  hltModE->AddTrigger("HLT_Ele17_SW_L1R_v2");
  hltModE->AddTrigger("HLT_Ele17_SW_Isol_L1R_v2");
  hltModE->AddTrigger("HLT_Ele17_SW_TighterEleIdIsol_L1R_v3");
  hltModE->AddTrigger("HLT_Ele17_SW_TightCaloEleId_Ele8HE_L1R_v2");
  hltModE->AddTrigger("HLT_Ele22_SW_L1R_v2");
  hltModE->AddTrigger("HLT_Ele22_SW_TighterCaloIdIsol_L1R_v2");
  hltModE->AddTrigger("HLT_Ele22_SW_TighterEleId_L1R_v3");
  hltModE->AddTrigger("HLT_Ele32_SW_TighterEleId_L1R_v2");
  hltModE->AddTrigger("HLT_DoubleEle4_SW_eeRes_L1R_v2");
  hltModE->AddTrigger("HLT_DoubleEle5_SW_Upsilon_L1R_v2");
  hltModE->AddTrigger("HLT_DoubleEle17_SW_L1R_v1");
  hltModE->AddTrigger("HLT_Photon10_Cleaned_L1R");
  hltModE->AddTrigger("HLT_Photon17_Isol_SC17HE_L1R_v1");
  hltModE->AddTrigger("HLT_Photon20_NoHE_L1R");
  hltModE->AddTrigger("HLT_Photon20_Cleaned_L1R");
  hltModE->AddTrigger("HLT_Photon20_Isol_Cleaned_L1R_v1");
  hltModE->AddTrigger("HLT_Photon22_SC22HE_L1R_v1");
  hltModE->AddTrigger("HLT_Photon30_Cleaned_L1R");
  hltModE->AddTrigger("HLT_Photon40_CaloId_Cleaned_L1R_v1");
  hltModE->AddTrigger("HLT_Photon40_Isol_Cleaned_L1R_v1");
  hltModE->AddTrigger("HLT_Photon50_Cleaned_L1R_v1");
  hltModE->AddTrigger("HLT_Photon50_NoHE_L1R");
  hltModE->AddTrigger("HLT_Photon70_Cleaned_L1R_v1");
  hltModE->AddTrigger("HLT_Photon110_NoHE_Cleaned_L1R_v1");
  hltModE->AddTrigger("HLT_DoublePhoton5_CEP_L1R_v3");
  hltModE->AddTrigger("HLT_DoublePhoton17_SingleIsol_L1R_v1");
  hltModE->AddTrigger("HLT_DoublePhoton22_L1R_v1");
  hltModE->SetTrigObjsName("MyHltElecObjs");
  hltModE->SetAbortIfNotAccepted(kFALSE);

  HLTMod *hltModP = new HLTMod("HLTModP");
  hltModP->AddTrigger("HLT_Photon10_L1R");
  hltModP->AddTrigger("HLT_Photon15_L1R");
  hltModP->AddTrigger("HLT_Photon15_Cleaned_L1R");
  hltModP->AddTrigger("HLT_Photon20_L1R");
  hltModP->AddTrigger("HLT_Photon17Isol_SC17HE_L1R");
  hltModP->AddTrigger("HLT_Ele15_SW_L1R");
  hltModP->AddTrigger("HLT_Ele15_SW_CaloEleId_L1R");
  hltModP->AddTrigger("HLT_Ele17_SW_CaloEleId_L1R");
  hltModP->AddTrigger("HLT_Ele17_SW_TightCaloEleId_SC8HE_L1R_v1");
  hltModP->AddTrigger("HLT_DoubleEle10_SW_L1R");
  hltModP->AddTrigger("HLT_Ele10_SW_L1R_v2");
  hltModP->AddTrigger("HLT_Ele12_SW_TighterEleId_L1R_v2");
  hltModP->AddTrigger("HLT_Ele17_SW_L1R_v2");
  hltModP->AddTrigger("HLT_Ele17_SW_Isol_L1R_v2");
  hltModP->AddTrigger("HLT_Ele17_SW_TighterEleIdIsol_L1R_v3");
  hltModP->AddTrigger("HLT_Ele17_SW_TightCaloEleId_Ele8HE_L1R_v2");
  hltModP->AddTrigger("HLT_Ele22_SW_L1R_v2");
  hltModP->AddTrigger("HLT_Ele22_SW_TighterCaloIdIsol_L1R_v2");
  hltModP->AddTrigger("HLT_Ele22_SW_TighterEleId_L1R_v3");
  hltModP->AddTrigger("HLT_Ele32_SW_TighterEleId_L1R_v2");
  hltModP->AddTrigger("HLT_DoubleEle4_SW_eeRes_L1R_v2");
  hltModP->AddTrigger("HLT_DoubleEle5_SW_Upsilon_L1R_v2");
  hltModP->AddTrigger("HLT_DoubleEle17_SW_L1R_v1");
  hltModP->AddTrigger("HLT_Photon10_Cleaned_L1R");
  hltModP->AddTrigger("HLT_Photon17_Isol_SC17HE_L1R_v1");
  hltModP->AddTrigger("HLT_Photon20_NoHE_L1R");
  hltModP->AddTrigger("HLT_Photon20_Cleaned_L1R");
  hltModP->AddTrigger("HLT_Photon20_Isol_Cleaned_L1R_v1");
  hltModP->AddTrigger("HLT_Photon22_SC22HE_L1R_v1");
  hltModP->AddTrigger("HLT_Photon30_Cleaned_L1R");
  hltModP->AddTrigger("HLT_Photon40_CaloId_Cleaned_L1R_v1");
  hltModP->AddTrigger("HLT_Photon40_Isol_Cleaned_L1R_v1");
  hltModP->AddTrigger("HLT_Photon50_Cleaned_L1R_v1");
  hltModP->AddTrigger("HLT_Photon50_NoHE_L1R");
  hltModP->AddTrigger("HLT_Photon70_Cleaned_L1R_v1");
  hltModP->AddTrigger("HLT_Photon110_NoHE_Cleaned_L1R_v1");
  hltModP->AddTrigger("HLT_DoublePhoton5_CEP_L1R_v3");
  hltModP->AddTrigger("HLT_DoublePhoton17_SingleIsol_L1R_v1");
  hltModP->AddTrigger("HLT_DoublePhoton22_L1R_v1");
  hltModP->SetTrigObjsName("MyHltPhotObjs");
  hltModP->SetAbortIfNotAccepted(kFALSE);
    
    
  //------------------------------------------------------------------------------------------------
  // select events with a good primary vertex
  //------------------------------------------------------------------------------------------------
  GoodPVFilterMod *goodPVFilterMod = new GoodPVFilterMod;
  goodPVFilterMod->SetMinVertexNTracks(0);
  goodPVFilterMod->SetMinNDof         (5);
  goodPVFilterMod->SetMaxAbsZ         (24.0);
  goodPVFilterMod->SetMaxRho          (2.0);

  //------------------------------------------------------------------------------------------------
  // object id and cleaning sequence
  //------------------------------------------------------------------------------------------------
  MuonIDMod *muonId = new MuonIDMod;  
  muonId->SetClassType ("Global");
  muonId->SetIDType    ("ZMuId");
  muonId->SetIsoType   ("TrackCaloSliding");
  muonId->SetApplyD0Cut(kTRUE);

  ElectronIDMod *elecId = new ElectronIDMod;
  elecId->SetIDType                    ("VBTFWorkingPoint95Id");
  elecId->SetIsoType                   ("TrackJuraSliding");
  elecId->SetApplyConversionFilterType1(kFALSE);
  elecId->SetApplyConversionFilterType2(kTRUE);
  elecId->SetChargeFilter              (kFALSE);
  elecId->SetApplyD0Cut                (kTRUE);
  elecId->SetNExpectedHitsInnerCut     (0);
  PhotonIDMod         *photId = new PhotonIDMod;

  ElectronCleaningMod *elecCleaning = new ElectronCleaningMod;
  PhotonCleaningMod   *photCleaning = new PhotonCleaningMod;

  //------------------------------------------------------------------------------------------------
  // merge modules
  //------------------------------------------------------------------------------------------------
  MergeLeptonsMod *mergeLeptonsMod = new MergeLeptonsMod;
  mergeLeptonsMod->SetMuonsName    (muonId->GetOutputName());
  mergeLeptonsMod->SetElectronsName(elecCleaning->GetOutputName());

  //------------------------------------------------------------------------------------------------
  // analyses modules
  //------------------------------------------------------------------------------------------------
  // Z -> mm analysis
  ZmmAnalysis *zmmMod = new ZmmAnalysis;
  zmmMod->SetTrigObjsName   (hltModM->GetOutputName());
  zmmMod->SetMuonName       (muonId->GetOutputName());
  zmmMod->SetMuonsFromBranch(kFALSE);
  zmmMod->SetOverlapCut(double(overlapCut));
  if (jsonFile.CompareTo("/home/cmsprod/json/~") != 0)
    zmmMod->SetIsData(kTRUE);
  else
    zmmMod->SetIsData(kFALSE);
  // Z -> ee analysis
  ZeeAnalysis *zeeMod = new ZeeAnalysis;
  zeeMod->SetTrigObjsName   (hltModE->GetOutputName());
  zeeMod->SetElecName       (elecId->GetOutputName());
  zeeMod->SetElecsFromBranch(kFALSE);
  zeeMod->SetOverlapCut     (double(overlapCut));
  if (jsonFile.CompareTo("/home/cmsprod/json/~") != 0)
    zeeMod->SetIsData(kTRUE);
  else
    zeeMod->SetIsData(kFALSE);
  HggAnalysis *hggMod = new HggAnalysis;
  hggMod->SetTrigObjsName     (hltModP->GetOutputName());
  hggMod->SetPhotonName       (photId->GetOutputName());
  hggMod->SetPhotonsFromBranch(kFALSE);
  hggMod->SetOverlapCut(double(overlapCut));
  if (jsonFile.CompareTo("/home/cmsprod/json/~") != 0)
    hggMod->SetIsData(kTRUE);
  else
    hggMod->SetIsData(kFALSE);

  //------------------------------------------------------------------------------------------------
  // making analysis chain
  //------------------------------------------------------------------------------------------------
  // this is how it always starts
  runLumiSel      ->Add(hltModM);
  // high level trigger is always first
  hltModM         ->Add(hltModE);
  hltModE         ->Add(hltModP);
  hltModP         ->Add(goodPVFilterMod);
  goodPVFilterMod ->Add(muonId);
  // simple object id modules
  muonId          ->Add(elecId);
  elecId          ->Add(photId);
  photId          ->Add(elecCleaning);
  // cleaning modules
  elecCleaning    ->Add(photCleaning);
  photCleaning    ->Add(mergeLeptonsMod);
  // lepton merging
  mergeLeptonsMod ->Add(zmmMod);
  // core analysis
  zmmMod          ->Add(zeeMod);
  zeeMod          ->Add(hggMod);

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
  if (TString(skim).CompareTo("noskim") == 0)
    d = c->FindDataset(book,dataset,fileset);
  else 
    d = c->FindDataset(book,skimdataset.Data(),fileset);
  ana->AddDataset(d);
  //ana->AddFile("root://castorcms//castor/cern.ch/user/p/paus/filler/011/s09-ttbar-7-mc3/*.root");
  //ana->AddFile("zee-skim_r10a-eg-pr-v4_noskim_0000_000.root");
  //ana->AddFile("zee-skim_p10-h110gg-gf-v26_noskim_0000_000.root");

  //------------------------------------------------------------------------------------------------
  // organize output
  //------------------------------------------------------------------------------------------------
  TString rootFile = TString(outputName);
  rootFile += TString("_") + TString(dataset) + TString("_") + TString(skim);
  if (TString(fileset) != TString(""))
    rootFile += TString("_") + TString(fileset);
  rootFile += TString(".root");
  ana->SetOutputName(rootFile.Data());
  ana->SetCacheSize(64*1024*1024);

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
