// $Id: runPhEffSF.C,v 1.1 2011/07/15 17:26:34 fabstoec Exp $
#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TSystem.h>
#include <TProfile.h>
#include <TFile.h>
#include "MitAna/DataUtil/interface/Debug.h"
#include "MitAna/Catalog/interface/Catalog.h"
#include "MitAna/TreeMod/interface/Analysis.h"
#include "MitAna/PhysicsMod/interface/PublisherMod.h"
#include "MitAna/PhysicsMod/interface/RunLumiSelectionMod.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitPhysics/Mods/interface/GoodPVFilterMod.h"
#include "MitPhysics/Mods/interface/PhotonPairSelector.h"
#include "MitPhysics/Mods/interface/PhotonTreeWriter.h"
#include "MitHgg/Mods/interface/MVASystematicsMod.h"

#endif

//--------------------------------------------------------------------------------------------------
void runMVASys(const char *fileset    = "0000",
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
  // select events with a good primary vertex
  //------------------------------------------------------------------------------------------------
  GoodPVFilterMod *goodPVFilterMod = new GoodPVFilterMod;
  goodPVFilterMod->SetMinVertexNTracks(0);
  goodPVFilterMod->SetMinNDof         (5);
  goodPVFilterMod->SetMaxAbsZ         (24.0);
  goodPVFilterMod->SetMaxRho          (2.0);
  goodPVFilterMod->SetIsMC            (!isData);
  
  MVASystematicsMod *sysMod = new MVASystematicsMod;

  TFile *fetacors = new TFile(gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/detacor.root"),"READ");
  TProfile *hetacors = (TProfile*) fetacors->Get("prof_profdatacortmp");
  hetacors->SetDirectory(0);
  fetacors->Close();

  PhotonPairSelector         *photmva = new PhotonPairSelector("PhotonPairSelectorMVA");
  photmva->SetOutputName("GoodPhotonsMVA");
  photmva->SetPhotonSelType("MVASelection");
  photmva->SetVertexSelType("CiCSelection");
  photmva->DoMCSmear(kTRUE);
  photmva->DoDataEneCorr(kTRUE);
  photmva->SetDoRegression(kTRUE);
  photmva->SetEtaCorrections(hetacors);
  photmva->SetMCSmearFactors(0.0110, 0.0121, 0.0340, 0.0285);
  photmva->AddEnCorrPerRun(160000,167913,0.00252,-0.00462,0.00379,-0.00196);
  photmva->AddEnCorrPerRun(170249,172619,-0.00006,-0.00627,0.01167,0.00544);
  photmva->AddEnCorrPerRun(172620,999999,-0.00149,-0.00691,0.02244,0.01111);
  photmva->SetLeadingPtMin(30);
  photmva->SetTrailingPtMin(22);
  photmva->SetBdtCutBarrel(-0.01);
  photmva->SetBdtCutEndcap(-0.02);
  photmva->SetIsData(isData);

  PhotonTreeWriter  * writerMIT = new PhotonTreeWriter;
  writerMIT->SetInputPhotonsName("GoodPhotonsMVA");
  writerMIT->     SetPhotonsFromBranch(false);
  writerMIT->     SetTupleName("MVAtuple");
  writerMIT->     SetIsData(isData);
  
  //------------------------------------------------------------------------------------------------
  // making analysis chain
  //------------------------------------------------------------------------------------------------
  // this is how it always starts
  runLumiSel       ->Add(goodPVFilterMod);
  goodPVFilterMod  ->Add(sysMod);
  sysMod           ->Add(photmva);
  photmva          ->Add(writerMIT);

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
