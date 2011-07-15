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
#include "MitHgg/Mods/interface/HggAnalysisMod.h"
#include "MitPhysics/Utils/interface/VertexTools.h"
#endif

//--------------------------------------------------------------------------------------------------
void runPhRes(const char *fileset    = "0000",
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
  goodPVFilterMod->SetIsMC(!isData);
  
  PhotonIDMod         *photId = new PhotonIDMod;
  photId->                SetIsoType("MITPUCorrected");
  photId->                SetApplySpikeRemoval(false);
  photId->                SetApplyPixelSeed(false);
  photId->                SetApplyElectronVetoConvRecovery(false);
  photId->                SetApplyConversionId(false);
  photId->                SetHadOverEmMax(0.05);
  photId->                SetPtMin(20.);
  photId->                SetEtaWidthEB(0.010);
  photId->                SetEtaWidthEE(0.028);
  photId->                SetAbsEtaMax(2.5);
  photId->                SetApplyR9Min(false);  

  // Pair Selectro for CiC Analysis
  PhotonPairSelector *photIdCiC = new PhotonPairSelector;
  photIdCiC->     SetIsData(isData);
  photIdCiC->     ApplyEleVeto(false);
  photIdCiC->     DoDataEneCorr(false);
  photIdCiC->     DoMCSmear(false);
  photIdCiC->     SetPhotonSelType("CiCSelection");
  photIdCiC->     SetVertexSelType("StdSelection");
  photIdCiC->     SetOutputName("CiCPhotons");
  photIdCiC->     SetTupleName("CiCTuple");

  // copy the Mod for the MIT selection
  PhotonPairSelector *photIdMIT = new PhotonPairSelector;
  photIdMIT->     SetIsData(isData);
  photIdMIT->     ApplyEleVeto(false);
  photIdMIT->     DoDataEneCorr(false);
  photIdMIT->     DoMCSmear(false);
  photIdMIT->     SetPhotonSelType("MITSelection");
  photIdMIT->     SetVertexSelType("StdSelection");
  photIdMIT->     SetInputPhotonsName(photId->GetOutputName());
  photIdMIT->     SetPhotonsFromBranch(false);
  photIdMIT->     SetOutputName("MITPhotons");
  photIdMIT->     SetPVName(goodPVFilterMod->GetOutputName());
  photIdMIT->     SetPVFromBranch(false);
  photIdMIT->     SetTupleName("MITTuple");

  // Two analysis Modules
  HggAnalysisMod *anaModCiC = new HggAnalysisMod;
  anaModCiC->SetPhotonName       (photIdCiC->GetOutputName());
  anaModCiC->SetPhotonsFromBranch(kFALSE);

  HggAnalysisMod *anaModMIT = new HggAnalysisMod;
  anaModMIT->SetPhotonName       (photIdMIT->GetOutputName());
  anaModMIT->SetPhotonsFromBranch(kFALSE);
  
  //------------------------------------------------------------------------------------------------
  // making analysis chain
  //------------------------------------------------------------------------------------------------
  // this is how it always starts
  runLumiSel      ->Add(hltModEle);

  // the MIT flow...
  hltModEle         ->Add(goodPVFilterMod);
  goodPVFilterMod ->Add(photId);
  photId          ->Add(photIdMIT);
  photIdMIT       ->Add(anaModMIT);

  // the CiC flow...
  hltModEle         ->Add(photIdCiC);
  photIdCiC       ->Add(anaModCiC);    

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
