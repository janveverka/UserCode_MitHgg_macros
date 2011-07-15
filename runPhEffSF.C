// $Id: runHggAna.C,v 1.2 2011/07/08 17:55:21 fabstoec Exp $
#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TSystem.h>
#include "MitAna/DataUtil/interface/Debug.h"
#include "MitAna/Catalog/interface/Catalog.h"
#include "MitAna/TreeMod/interface/Analysis.h"
#include "MitAna/TreeMod/interface/HLTMod.h"
#include "MitAna/PhysicsMod/interface/PublisherMod.h"
#include "MitAna/PhysicsMod/interface/RunLumiSelectionMod.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitPhysics/Mods/interface/GoodPVFilterMod.h"
#include "MitPhysics/Mods/interface/PhotonIDMod.h"
#include "MitPhysics/Mods/interface/PhotonPairSelector.h"
#include "MitHgg/Mods/interface/PhotonMCDataSFMod.h"
#include "MitPhysics/Mods/interface/ElectronIDMod.h"
#endif

//--------------------------------------------------------------------------------------------------
void runPhEffSF(const char *fileset    = "0000",
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
  //TString jsonFile = TString("/home/cmsprod/json/") + TString(json);
  TString jsonFile = TString("/home/fabstoec/cms/json/") + TString(json);
  Bool_t  isData   = ( (jsonFile.CompareTo("/home/fabstoec/cms/json/~") != 0) );
  
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
  if ((jsonFile.CompareTo("/home/fabstoec/cms/json/~") != 0) &&
      (jsonFile.CompareTo("/home/fabstoec/cms/json/-") != 0)   ) {
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
  HLTMod *hltModEle = new HLTMod("HltModEle");
  hltModEle->AddTrigger("HLT_Ele32_CaloIdL_CaloIsoVL_SC17_v1");
  hltModEle->AddTrigger("HLT_Ele32_CaloIdL_CaloIsoVL_SC17_v2");
  hltModEle->AddTrigger("HLT_Ele32_CaloIdL_CaloIsoVL_SC17_v3");
  hltModEle->AddTrigger("HLT_Ele32_CaloIdL_CaloIsoVL_SC17_v4");
  hltModEle->AddTrigger("HLT_Ele32_CaloIdL_CaloIsoVL_SC17_v5");
  hltModEle->AddTrigger("HLT_Ele32_CaloIdL_CaloIsoVL_SC17_v6");
  hltModEle->AddTrigger("HLT_Ele32_CaloIdL_CaloIsoVL_SC17_v7");
  hltModEle->AddTrigger("HLT_Ele32_CaloIdL_CaloIsoVL_SC17_v8");
  hltModEle->SetTrigObjsName("MyHltPhotObjsEle");
  hltModEle->SetAbortIfNotAccepted(isData);

  //------------------------------------------------------------------------------------------------
  // select events with a good primary vertex
  //------------------------------------------------------------------------------------------------
  GoodPVFilterMod *goodPVFilterMod = new GoodPVFilterMod;
  goodPVFilterMod->SetMinVertexNTracks(0);
  goodPVFilterMod->SetMinNDof         (5);
  goodPVFilterMod->SetMaxAbsZ         (24.0);
  goodPVFilterMod->SetMaxRho          (2.0);
  goodPVFilterMod->SetIsMC            (!isData);
  
  // our MIT Photons ID...
  PhotonIDMod         *photId = new PhotonIDMod;
  photId->                SetIsoType("MITPUCorrected");
  photId->                SetApplySpikeRemoval(false);
  photId->                SetApplyPixelSeed(false);
  photId->                SetApplyElectronVetoConvRecovery(false); // no conversion, we look at electrons
  photId->                SetApplyConversionId(false);
  photId->                SetHadOverEmMax(0.05);
  photId->                SetPtMin(20.);
  photId->                SetEtaWidthEB(0.010);
  photId->                SetEtaWidthEE(0.028);
  photId->                SetAbsEtaMax(2.5);
  photId->                SetApplyR9Min(false);  
  photId->SetOutputName("MITPhotons");

  // ... and the CiC ID...
  PhotonIDMod         *photIdCiC = new PhotonIDMod;
  photIdCiC->             SetIDType("BaseLineCiC");
  photIdCiC->             SetPtMin(20.);
  photIdCiC->             SetApplyElectronVeto(false);
  photIdCiC->SetOutputName("CiCPhotons");

  // some default VBTF Electron ID
  ElectronIDMod *elecId = new ElectronIDMod;
  elecId->SetIDType                    ("VBTFWorkingPoint80Id");
  elecId->SetIsoType                   ("TrackJuraSliding");
  elecId->SetApplyConversionFilterType1(kFALSE);
  elecId->SetApplyConversionFilterType2(kTRUE);
  elecId->SetChargeFilter              (kFALSE);
  elecId->SetApplyD0Cut                (kTRUE);
  elecId->SetNExpectedHitsInnerCut     (0);
  elecId->SetVertexName                (goodPVFilterMod->GetOutputName());
  elecId->SetPtMin                     (30.);  

  PhotonMCDataSFMod *phSFMod = new PhotonMCDataSFMod;
  phSFMod->     SetIsData(isData);
  phSFMod->     SetInputPhotonsName(photId->GetOutputName());
  phSFMod->     SetPhotonsFromBranch(false);
  phSFMod->     SetPVName(goodPVFilterMod->GetOutputName());
  phSFMod->     SetPVFromBranch(false);
  phSFMod->     SetElectronName(elecId->GetOutputName());
  phSFMod->     SetEleFromBranch(false);

  phSFMod->     AddTriggerLegs("hltEle32CaloIdLCaloIsoVLSC17PixelMatchFilter","hltEle32CaloIdLCaloIsoVLSC17HEDoubleFilter");
  phSFMod->     SetTrigObjsName(hltModEle->GetOutputName());
  phSFMod->     SetTupleName("MITtuple");


  PhotonMCDataSFMod *phSFModCiC = new PhotonMCDataSFMod;
  phSFModCiC->     SetIsData(isData);
  phSFModCiC->     SetInputPhotonsName(photIdCiC->GetOutputName());
  phSFModCiC->     SetPhotonsFromBranch(false);
  phSFModCiC->     SetElectronName(elecId->GetOutputName());
  phSFModCiC->     SetEleFromBranch(false);

  phSFModCiC->     AddTriggerLegs("hltEle32CaloIdLCaloIsoVLSC17PixelMatchFilter","hltEle32CaloIdLCaloIsoVLSC17HEDoubleFilter");
  phSFModCiC->     SetTrigObjsName(hltModEle->GetOutputName());
  phSFModCiC->     SetTupleName("CiCtuple");

  //------------------------------------------------------------------------------------------------
  // making analysis chain
  //------------------------------------------------------------------------------------------------
  // this is how it always starts
  runLumiSel      ->Add(hltModEle);
  hltModEle         ->Add(goodPVFilterMod);
  goodPVFilterMod   ->Add(elecId);

  // ... the MIT flow ...
  elecId            ->Add(photId);
  photId            ->Add(phSFMod);

  // and the CiC flow...
  elecId            ->Add(photIdCiC);
  photIdCiC         ->Add(phSFModCiC);

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
