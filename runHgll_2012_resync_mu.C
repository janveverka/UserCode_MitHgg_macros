// $Id: runHgll_2012_resync_mu.C,v 1.1 2012/10/29 16:39:04 khahn Exp $
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
#include "MitPhysics/Mods/interface/LeptonPairPhotonTreeWriter.h"
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
void runHgll_2012_resync_mu(const char *fileset    = "0000",
		  const char *skim       = "noskim",
		  const char *dataset    = "r11a-del-j16-v1",   
		  const char *book       = "local/filefi/025",
		  const char *catalogDir = "/home/cmsprod/catalog",
		  const char *outputName = "hgll",
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

  TString jsonFile = TString("/home/fabstoec/cms/json/") + TString(json);
  //TString jsonFile = TString("/home/auhess/cms/json/") + TString("Cert_136033-149442_7TeV_Dec22ReReco_Collisions10_JSON_v4.txt");

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


  //isData = kFALSE;
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
  runLumiSel->SetAcceptMC( kTRUE );                          // Monte Carlo events are always accepted
  
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
  //printf("\n Run lumi worked. \n\n");

  //------------------------------------------------------------------------------------------------
  // HLT information
  //------------------------------------------------------------------------------------------------
  HLTMod *hltModll = new HLTMod("HLTModll");
//     hltModll->AddTrigger("HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v*",161217, 166967);  
//     hltModll->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*",167039, 196531);
//     hltModll->AddTrigger("HLT_DoubleMu7_v*", 160431, 163869);
//     hltModll->AddTrigger("HLT_Mu13_Mu8_v*",  165088, 178380);
//    hltModll->AddTrigger("HLT_Mu17_Mu8_v*",  178420, 196531);
    hltModll->AddTrigger("HLT_Mu17_Mu8_v*",  0, 196531);
    //    hltModll->AddTrigger("HLT_Mu17_TkMu8_v*",  178420, 196531);
  //  hltModll->AddTrigger("HLT_Mu22_Mu8_v*",  178420, 196531);
  //hltModll->AddTrigger("HLT_Mu22_TkMu8_v*",  178420, 196531);
    //  hltModll->AddTrigger("HLT_IsoMu24_v*",  178420, 196531);
    //  hltModll->AddTrigger("HLT_IsoMu24_eta2p1_v*",  178420, 196531);
   hltModll->SetTrigObjsName("MyHltElecObjs");
  //  hltModll->SetAbortIfNotAccepted( false );
  hltModll->SetAbortIfNotAccepted( true );

  //------------------------------------------------------------------------------------------------
  // select events with a good primary vertex
  //------------------------------------------------------------------------------------------------
  GoodPVFilterMod *goodPVFilterMod = new GoodPVFilterMod;
  goodPVFilterMod->SetMinVertexNTracks(0);
  goodPVFilterMod->SetMinNDof         (4.0);
  goodPVFilterMod->SetMaxAbsZ         (24.0);
  goodPVFilterMod->SetMaxRho          (2.0);
  //goodPVFilterMod->SetAbortIfNotAccepted(kFALSE);
  goodPVFilterMod->SetAbortIfNotAccepted(kTRUE);
  goodPVFilterMod->SetIsMC(!isData);
    
  //------------------------------------------------------------------------------------------------
  // object id and cleaning sequence
  //------------------------------------------------------------------------------------------------
  
  // KH for sync
//   PhotonMvaMod *photreg = new PhotonMvaMod;
//   photreg->SetOutputName("GoodPhotonsRegr");
//   photreg->SetIsData(isData);  
//   photreg->SetMinNumPhotons(1);
//   photreg->SetDoPreselection(kFALSE);

  PhotonIDMod *myPhId = new PhotonIDMod;
  myPhId -> SetPtMin    (10.);
  myPhId -> SetAbsEtaMax   (2.5);
  // KH for sync
  //  myPhId -> SetInputName(photreg->GetOutputName());
  myPhId -> SetGoodElectronsFromBranch( true );
  // KH for sync
  //  myPhId -> SetPhotonsFromBranch(kFALSE);
  myPhId -> SetPhotonsFromBranch( true);

  myPhId -> SetOutputName("outputPhotons");
  myPhId -> SetIDType("TrivialSelection");

  myPhId -> DoMCSmear( true );
  //  myPhId -> DoDataEneCorr( false );

  myPhId -> SetMCSmearFactors(0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02);

  //  myPhId -> AddEnCorrPerRun(160431,167913,0.9941,0.9941,1.0004,0.9916,1.0045,1.0033,1.0082,0.9958,1.0064);//ming:Emc/Edata
  //  myPhId -> AddEnCorrPerRun(170000,172619,0.9954,0.9954,1.0016,0.9937,1.0066,0.9976,1.0025,0.9940,1.0046);
  //  myPhId -> AddEnCorrPerRun(172620,173692,0.9955,0.9955,1.0017,0.9929,1.0058,0.9986,1.0035,0.9923,1.0029);
  //  myPhId -> AddEnCorrPerRun(175830,177139,0.9958,0.9958,1.0021,0.9944,1.0073,0.9968,1.0017,0.9933,1.004);
  //  myPhId -> AddEnCorrPerRun(177140,178421,0.9962,0.9962,1.0025,0.9946,1.0075,0.9960,1.0010,0.9944,1.005);
  //  myPhId -> AddEnCorrPerRun(178424,180252,0.9961,0.9961,1.0024,0.9942,1.0071,0.9921,0.9970,0.9953,1.0059);
  //  myPhId -> SetDoMCErrScaling(kTRUE);
  //  myPhId -> SetMCErrScale(1.07, 1.045);

  myPhId -> SetIsData(isData);
  //  myPhId -> SetDoShowerShapeScaling(kTRUE);
  //  myPhId -> SetShowerShapeType("2011ShowerShape");


  PhotonIDMod *myPhIdNoSmear = new PhotonIDMod;
  myPhIdNoSmear -> SetPtMin    (10.);
  myPhIdNoSmear -> SetAbsEtaMax   (2.5);
  // KH for sync
  //  myPhIdNoSmear -> SetInputName(photreg->GetOutputName());
  myPhIdNoSmear -> SetGoodElectronsFromBranch( true );
  // KH for sync
  //  myPhIdNoSmear -> SetPhotonsFromBranch(kFALSE);
  myPhIdNoSmear -> SetPhotonsFromBranch(true);
  myPhIdNoSmear -> SetOutputName("outputPhotonsNoSmear");
  myPhIdNoSmear -> SetIDType("TrivialSelection");

  myPhIdNoSmear -> DoMCSmear( false );
  //  myPhIdNoSmear -> DoDataEneCorr( false );

  myPhIdNoSmear -> SetMCSmearFactors(0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02);



  //  myPhIdNoSmear -> AddEnCorrPerRun(160431,167913,0.9941,0.9941,1.0004,0.9916,1.0045,1.0033,1.0082,0.9958,1.0064);//ming:Emc/Edata
  //  myPhIdNoSmear -> AddEnCorrPerRun(170000,172619,0.9954,0.9954,1.0016,0.9937,1.0066,0.9976,1.0025,0.9940,1.0046);
  //  myPhIdNoSmear -> AddEnCorrPerRun(172620,173692,0.9955,0.9955,1.0017,0.9929,1.0058,0.9986,1.0035,0.9923,1.0029);
  //  myPhIdNoSmear -> AddEnCorrPerRun(175830,177139,0.9958,0.9958,1.0021,0.9944,1.0073,0.9968,1.0017,0.9933,1.004);
  //  myPhIdNoSmear -> AddEnCorrPerRun(177140,178421,0.9962,0.9962,1.0025,0.9946,1.0075,0.9960,1.0010,0.9944,1.005);
  //  myPhIdNoSmear -> AddEnCorrPerRun(178424,180252,0.9961,0.9961,1.0024,0.9942,1.0071,0.9921,0.9970,0.9953,1.0059);
  //  myPhIdNoSmear -> SetDoMCErrScaling(kTRUE);
  //  myPhIdNoSmear -> SetMCErrScale(1.07, 1.045);

  myPhIdNoSmear -> SetIsData(isData);
  //  myPhIdNoSmear -> SetDoShowerShapeScaling(kTRUE);
  //  myPhIdNoSmear -> SetShowerShapeType("2011ShowerShape");


    
  LeptonPairPhotonTreeWriter *eegtree = new LeptonPairPhotonTreeWriter;
  eegtree ->SetInputPhotonsName(myPhId->GetOutputName());
  eegtree ->SetPhotonsFromBranch( false );
  eegtree ->SetGoodElectronsFromBranch( true );
  eegtree ->SetIsData(isData);
  eegtree ->SetGoodMuonsFromBranch( true ); 
  eegtree ->SetYear(2012); 
  eegtree ->SetPhosphorDataFile("../../MyPhosphorDir/PHOSPHOR_NUMBERS_EXPFIT.txt");

  //  eegtree ->SetDoDataEleEneCorr( true );
  //  eegtree ->SetDataEleEneCorr( 1., 1. );
//   eegtree ->SetDoMCEleEneSmear( true );
//   eegtree ->SetMCEleEneSmear( 0.016, 0.028 );
  eegtree->SetTupleName("nominalTree");

  LeptonPairPhotonTreeWriter *eegtreeNoEleSmear = new LeptonPairPhotonTreeWriter;
  eegtreeNoEleSmear ->SetInputPhotonsName(myPhId->GetOutputName());
  eegtreeNoEleSmear ->SetPhotonsFromBranch( false );
  eegtreeNoEleSmear ->SetGoodElectronsFromBranch( true );
  eegtreeNoEleSmear ->SetIsData(isData);
  eegtreeNoEleSmear ->SetGoodMuonsFromBranch( true ); 
  eegtreeNoEleSmear ->SetYear(2012); 

  //  eegtreeNoEleSmear ->SetDoDataEleEneCorr( true );
  //  eegtreeNoEleSmear ->SetDataEleEneCorr( 1., 1. );
//   eegtreeNoEleSmear ->SetDoMCEleEneSmear( false );
//   eegtreeNoEleSmear ->SetMCEleEneSmear( 0.016, 0.028 );
  eegtreeNoEleSmear->SetTupleName("noEleSmearTree");

  LeptonPairPhotonTreeWriter *eegtreeNoPhSmear = new LeptonPairPhotonTreeWriter;
  eegtreeNoPhSmear ->SetInputPhotonsName(myPhIdNoSmear->GetOutputName());
  eegtreeNoPhSmear ->SetPhotonsFromBranch( false );
  eegtreeNoPhSmear ->SetGoodElectronsFromBranch( true );
  eegtreeNoPhSmear ->SetIsData(isData);
  eegtreeNoPhSmear ->SetGoodMuonsFromBranch( true ); 
  eegtreeNoPhSmear ->SetYear(2012); 

  //  eegtreeNoPhSmear ->SetDoDataEleEneCorr( true );
  //  eegtreeNoPhSmear ->SetDataEleEneCorr( 1., 1. );
//   eegtreeNoPhSmear ->SetDoMCEleEneSmear( true );
//   eegtreeNoPhSmear ->SetMCEleEneSmear( 0.016, 0.028 );
  eegtreeNoPhSmear->SetTupleName("noPhSmearTree");

  LeptonPairPhotonTreeWriter *eegtreeNoPhSmearNoEleSmear = new LeptonPairPhotonTreeWriter;
  eegtreeNoPhSmearNoEleSmear ->SetInputPhotonsName(myPhIdNoSmear->GetOutputName());
  eegtreeNoPhSmearNoEleSmear ->SetPhotonsFromBranch( false );
  eegtreeNoPhSmearNoEleSmear ->SetGoodElectronsFromBranch( true );
  eegtreeNoPhSmearNoEleSmear ->SetIsData(isData);
  eegtreeNoPhSmearNoEleSmear ->SetGoodMuonsFromBranch( true ); 
  eegtreeNoPhSmearNoEleSmear ->SetYear(2012); 
  eegtreeNoPhSmearNoEleSmear ->SetVerbose(true); 
  eegtreeNoPhSmearNoEleSmear ->SetDoElectronChannel(false); 

  //  eegtreeNoPhSmearNoEleSmear ->SetDoDataEleEneCorr( true );
  //  eegtreeNoPhSmearNoEleSmear ->SetDataEleEneCorr( 1., 1. );
//   eegtreeNoPhSmearNoEleSmear ->SetDoMCEleEneSmear( false );
//   eegtreeNoPhSmearNoEleSmear ->SetMCEleEneSmear( 0.016, 0.028 );
  eegtreeNoPhSmearNoEleSmear->SetTupleName("noPhSmearNoEleSmearTree");



  //------------------------------------------------------------------------------------------------
  // making analysis chain
  //------------------------------------------------------------------------------------------------
  // this is how it always starts

  if (TString(dataset).Contains("-h")) {    
    runLumiSel        ->Add(sysMod);
  }
  // high level trigger is always first
  runLumiSel      ->Add(hltModll);
  hltModll        ->Add(goodPVFilterMod);

  // KH replace for sync
//   goodPVFilterMod ->Add(photreg);
//   photreg->Add(myPhId);
//   photreg->Add(myPhIdNoSmear);

//  goodPVFilterMod ->Add(myPhId);
  goodPVFilterMod ->Add(myPhIdNoSmear);

  //  myPhId ->Add(eegtree);
  //  myPhId ->Add(eegtreeNoEleSmear);

  //  myPhIdNoSmear ->Add(eegtreeNoPhSmear);
  myPhIdNoSmear ->Add(eegtreeNoPhSmearNoEleSmear);


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
  //  ana->AddDataset(d);

  string fname;
  if( bookstr.Length() > 0 ) { 
    ifstream f(bookstr.Data());
    while (f >> fname) { 
      if( !(strncmp( fname.c_str(), "#", 1 ) ) ) continue; // skip commented lines
      cout << "adding inputfile : " << fname.c_str() << endl;
      //      entrymap[string(fname.c_str())] = unskimmedEntries(fname.c_str());
      //      cout << "unskimmed entries: " << entrymap[string(fname.c_str())] << endl;
      //      total_unskimmed += entrymap[string(fname.c_str())];
      ana->AddFile(fname.c_str());
    }
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
  ana->SetCacheSize(64*1024*1024);
  //ana->SetCacheSize(0);
  
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
