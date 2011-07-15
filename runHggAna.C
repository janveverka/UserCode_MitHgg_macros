// $Id: runHggAna.C,v 1.3 2011/07/15 17:26:34 fabstoec Exp $
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
  hltModP->AddTrigger("HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v1");
  hltModP->AddTrigger("HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v2");
  hltModP->AddTrigger("HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v3");
  hltModP->AddTrigger("HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v4");
  hltModP->AddTrigger("HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v5");
  hltModP->AddTrigger("HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v6");
  hltModP->AddTrigger("HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v7");
  hltModP->AddTrigger("HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v8");

  hltModP->AddTrigger("HLT_Photon20_R9Id_Photon18_R9Id_v1");
  hltModP->AddTrigger("HLT_Photon20_R9Id_Photon18_R9Id_v2");
  hltModP->AddTrigger("HLT_Photon20_R9Id_Photon18_R9Id_v3");
  hltModP->AddTrigger("HLT_Photon20_R9Id_Photon18_R9Id_v4");
  hltModP->AddTrigger("HLT_Photon20_R9Id_Photon18_R9Id_v5");
  hltModP->AddTrigger("HLT_Photon20_R9Id_Photon18_R9Id_v6");
  hltModP->AddTrigger("HLT_Photon20_R9Id_Photon18_R9Id_v7");
  hltModP->AddTrigger("HLT_Photon20_R9Id_Photon18_R9Id_v8");

  hltModP->SetTrigObjsName("MyHltPhotObjs");
  hltModP->SetAbortIfNotAccepted(isData);

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
  photId->                SetApplyElectronVetoConvRecovery(true);
  photId->                SetApplyConversionId(true);
  photId->                SetHadOverEmMax(0.02);
  photId->                SetPtMin(20.);
  photId->                SetEtaWidthEB(0.010);
  photId->                SetEtaWidthEE(0.028);
  photId->                SetAbsEtaMax(2.5);
  photId->                SetApplyR9Min(false);  

  // Pair Selectro for CiC Analysis
  PhotonPairSelector *photIdCiC = new PhotonPairSelector;
  photIdCiC->     SetIsData(isData);

  photIdCiC->     SetMCSmearFactors(0.0092, 0.0170,    // EB high/low R9
			 	 0.0292, 0.0289);   // EE high/low R9

  photIdCiC->     AddEnCorrPerRun  (160404, 163869,    // Run Range
				 -0.0047,  0.0025,  // EB high/low R9
				  0.0058, -0.0010); // EE high/low R9
  photIdCiC->     AddEnCorrPerRun  (165071, 165970,    // Run Range
				 -0.0007,  0.0049,  // EB high/low R9
				  0.0249,  0.0062); // EE high/low R9
  photIdCiC->     AddEnCorrPerRun  (165971, 166502,    // Run Range
				  0.0003,  0.0067,  // EB high/low R9
				  0.0376,  0.0133); // EE high/low R9
  photIdCiC->     AddEnCorrPerRun  (166503, 166861,    // Run Range
				  0.0011,  0.0063,  // EB high/low R9
				  0.0450,  0.0178); // EE high/low R9
  photIdCiC->     AddEnCorrPerRun  (166862, 999999,    // Run Range
				  0.0014,  0.0074,  // EB high/low R9
				  0.0561,  0.0273); // EE high/low R9
  photIdCiC->     SetPhotonSelType("CiCSelection");
  photIdCiC->     SetVertexSelType("CiCSelection");
  photIdCiC->     SetOutputName("CiCPhotons");
  photIdCiC->     SetTupleName("CiCtuple");

  VertexTools* vtool = VertexTools::instance(gSystem->Getenv("CMSSW_BASE"));

  // copy the Mod for the MIT selection
  PhotonPairSelector *photIdMIT = new PhotonPairSelector;
  photIdMIT->     SetIsData(isData);

  photIdMIT->     SetMCSmearFactors(0.0092, 0.0170,    // EB high/low R9
			 	 0.0292, 0.0289);   // EE high/low R9

  photIdMIT->     AddEnCorrPerRun  (160404, 163869,    // Run Range
				 -0.0047,  0.0025,  // EB high/low R9
				  0.0058, -0.0010); // EE high/low R9
  photIdMIT->     AddEnCorrPerRun  (165071, 165970,    // Run Range
				 -0.0007,  0.0049,  // EB high/low R9
				  0.0249,  0.0062); // EE high/low R9
  photIdMIT->     AddEnCorrPerRun  (165971, 166502,    // Run Range
				  0.0003,  0.0067,  // EB high/low R9
				  0.0376,  0.0133); // EE high/low R9
  photIdMIT->     AddEnCorrPerRun  (166503, 166861,    // Run Range
				  0.0011,  0.0063,  // EB high/low R9
				  0.0450,  0.0178); // EE high/low R9
  photIdMIT->     AddEnCorrPerRun  (166862, 999999,    // Run Range
				  0.0014,  0.0074,  // EB high/low R9
				  0.0561,  0.0273); // EE high/low R9
  photIdMIT->     SetPhotonSelType("MITSelection");
  photIdMIT->     SetVertexSelType("StdSelection");
  photIdMIT->     SetInputPhotonsName(photId->GetOutputName());
  photIdMIT->     SetPhotonsFromBranch(false);
  photIdMIT->     SetOutputName("MITPhotons");
  photIdMIT->     SetPVName(goodPVFilterMod->GetOutputName());
  photIdMIT->     SetPVFromBranch(false);
  photIdMIT->     SetTupleName("MITtuple");

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
  runLumiSel      ->Add(hltModP);

  // the MIT flow...
  hltModP         ->Add(goodPVFilterMod);
  goodPVFilterMod ->Add(photId);
  photId          ->Add(photIdMIT);
  photIdMIT       ->Add(anaModMIT);

  // the CiC flow...
  hltModP         ->Add(photIdCiC);
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
