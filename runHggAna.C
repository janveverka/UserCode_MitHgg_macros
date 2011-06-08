// $Id: runPhIso.C,v 1.4 2011/02/01 16:51:37 bendavid Exp $
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
#include "MitHgg/Mods/interface/HggAnalysisMod.h"
#endif

//--------------------------------------------------------------------------------------------------
void runHggAna(const char *fileset    = "0000",
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

  HLTMod *hltModP = new HLTMod("HLTModP");
  hltModP->AddTrigger("HLT_Photon26_IsoVL_Photon18_v1");
  hltModP->AddTrigger("HLT_Photon26_IsoVL_Photon18_v2");
  hltModP->AddTrigger("HLT_Photon26_IsoVL_Photon18_v3");
  hltModP->AddTrigger("HLT_Photon26_IsoVL_Photon18_v5");
  hltModP->AddTrigger("HLT_Photon26_IsoVL_Photon18_v6");
  hltModP->AddTrigger("HLT_Photon26_IsoVL_Photon18_v7");
  hltModP->AddTrigger("HLT_Photon26_IsoVL_Photon18_v8");

  hltModP->SetTrigObjsName("MyHltPhotObjs");
  if (isData)
    hltModP->SetAbortIfNotAccepted(kTRUE);
  else
    hltModP->SetAbortIfNotAccepted(kFALSE);
  
  //------------------------------------------------------------------------------------------------
  // select events with a good primary vertex
  //------------------------------------------------------------------------------------------------
  GoodPVFilterMod *goodPVFilterMod = new GoodPVFilterMod;
  goodPVFilterMod->SetMinVertexNTracks(0);
  goodPVFilterMod->SetMinNDof         (5);
  goodPVFilterMod->SetMaxAbsZ         (24.0);
  goodPVFilterMod->SetMaxRho          (2.0);
  goodPVFilterMod->SetVertexesName("DAPrimaryVertexes");
  
  PhotonIDMod         *photId = new PhotonIDMod;
  photId->                SetIsoType("MITPUCorrected");
  photId->                SetApplySpikeRemoval(true);
  photId->                SetApplyPixelSeed(false);
  photId->                SetApplyElectronVetoConvRecovery(true);
  photId->                SetApplyConversionId(true);
  photId->                SetHadOverEmMax(0.02);
  photId->                SetPtMin(30.);
  photId->                SetEtaWidthEB(0.010);
  photId->                SetEtaWidthEE(0.028);
  photId->                SetAbsEtaMax(2.5);
  photId->                SetApplyR9Min(false);  

  HggAnalysisMod *anaMod = new HggAnalysisMod;
  anaMod->SetTrigObjsName     (hltModP->GetOutputName());
  anaMod->SetPhotonName       (photId->GetOutputName());
  anaMod->SetPhotonsFromBranch(kFALSE);
  anaMod->SetOverlapCut(double(overlapCut));
  anaMod->SetPVName(ModNames::gkGoodVertexesName);

  if (jsonFile.CompareTo("/home/fabstoec/cms/json/~") != 0)
    anaMod->SetIsData(kTRUE);
  else
    anaMod->SetIsData(kFALSE);

  //------------------------------------------------------------------------------------------------
  // making analysis chain
  //------------------------------------------------------------------------------------------------
  // this is how it always starts
  runLumiSel      ->Add(hltModP);
  hltModP         ->Add(goodPVFilterMod);
  goodPVFilterMod ->Add(photId);
  photId          ->Add(anaMod);
  
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
