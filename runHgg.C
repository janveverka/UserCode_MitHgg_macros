// $Id: runHgg.C,v 1.5 2012/03/31 23:27:25 paus Exp $
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

using namespace mithep;

int  decodeEnv(char* json, char* overlap, float overlapCut);

//--------------------------------------------------------------------------------------------------
void runHgg(const char *fileset    = "",
            const char *skim       = "noskim",
            const char *dataset    = "r12a-pho-j22-v1",
            const char *book       = "local/filefi/029",
            const char *catalogDir = "/home/cmsprod/catalog",
            const char *outputName = "hgg",
            int         nEvents    = 5000)
{
  //------------------------------------------------------------------------------------------------
  // some parameters get passed through the environment
  //------------------------------------------------------------------------------------------------
  char json[1024], overlap[1024];
  float overlapCut = -1;

  if (decodeEnv(json,overlap,overlapCut) != 0)
    return;

  TString jsonFile = TString("/home/cmsprod/cms/json/") + TString(json);
  Bool_t  isData   = (jsonFile.CompareTo("/home/cmsprod/cms/json/~") != 0);

  printf("\n Initialization worked: \n\n");
  printf("   JSON   : %s (file: %s)\n",json,jsonFile.Data());
  printf("   OVERLAP: %s\n\n",overlap);
  printf("   isData : %d\n\n",isData);

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
  if ((jsonFile.CompareTo("/home/cmsprod/cms/json/~") != 0) &&
      (jsonFile.CompareTo("/home/cmsprod/cms/json/-") != 0)   ) {
    runLumiSel->AddJSONFile(jsonFile.Data());
  }
  if ((jsonFile.CompareTo("/home/cmsprod/cms/json/-") == 0)   ) {
    printf("\n WARNING -- Looking at data without JSON file: always accept.\n\n");
    runLumiSel->SetAbortIfNotAccepted(kFALSE);   // accept all events if there is no valid JSON file
  }

  printf("\n Run lumi worked. \n\n");

  //------------------------------------------------------------------------------------------------
  // HLT information
  //------------------------------------------------------------------------------------------------
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
  hltModE->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v10",179890,999999);
  hltModE->SetTrigObjsName("MyHltElecObjs");
  hltModE->SetAbortIfNotAccepted(isData);

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

  //------------------------------------------------------------------------------------------------
  // select events with a good primary vertex
  //------------------------------------------------------------------------------------------------
  GoodPVFilterMod *goodPVFilterMod = new GoodPVFilterMod;
  goodPVFilterMod->SetMinVertexNTracks  (0);
  goodPVFilterMod->SetMinNDof           (4.0);
  goodPVFilterMod->SetMaxAbsZ           (24.0);
  goodPVFilterMod->SetMaxRho            (2.0);
  goodPVFilterMod->SetAbortIfNotAccepted(kFALSE);
  goodPVFilterMod->SetIsMC              (!isData);

  GoodPVFilterMod *goodPVFilterModE = new GoodPVFilterMod("GoodPVFilterModE");
  goodPVFilterModE->SetOutputName       ("GoodVertexesE");
  goodPVFilterModE->SetMinVertexNTracks (0);
  goodPVFilterModE->SetMinNDof          (4.0);
  goodPVFilterModE->SetMaxAbsZ          (24.0);
  goodPVFilterModE->SetMaxRho           (2.0);
  goodPVFilterModE->SetIsMC             (!isData);

  //------------------------------------------------------------------------------------------------
  // object id and cleaning sequence
  //------------------------------------------------------------------------------------------------
  //ElectronIDMod *elecId = new ElectronIDMod;
  //elecId->SetVertexName           (goodPVFilterModE->GetOutputName());
  //elecId->SetIDType               ("VBTFWorkingPointLowPtId");
  //elecId->SetIsoType              ("PFIso");
  //elecId->SetNExpectedHitsInnerCut(0);
  //elecId->SetEtaMax               (999.);
  //elecId->SetChargeFilter         (kFALSE);
  //elecId->SetApplySpikeRemoval    (kFALSE);
  //elecId->SetApplyEcalFiducial    (kTRUE);
  //elecId->SetApplyEcalSeeded      (kTRUE);
  //elecId->SetPtMin                (20.0);

  PublisherMod<PFJet,Jet> *pubJet = new PublisherMod<PFJet,Jet>("JetPub");
  pubJet->SetInputName ("AKt5PFJets");
  pubJet->SetOutputName("PubAKt5PFJets");

  PublisherMod<PFJet,Jet> *pubJetOpen = new PublisherMod<PFJet,Jet>("JetPubOpen");
  pubJetOpen->SetInputName ("AKt5PFJets");
  pubJetOpen->SetOutputName("PubAKt5PFJetsOpen");

  JetCorrectionMod *jetCorr = new JetCorrectionMod;
  std::string base = std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/")).Data());

  jetCorr->AddCorrectionFromFile(base + std::string("START42_V12_AK5PF_L1FastJet.txt"));
  jetCorr->AddCorrectionFromFile(base + std::string("START42_V12_AK5PF_L2Relative.txt"));
  jetCorr->AddCorrectionFromFile(base + std::string("START42_V12_AK5PF_L3Absolute.txt"));
  if (isData)
    jetCorr->AddCorrectionFromFile(base + std::string("START42_V12_AK5PF_L2L3Residual.txt"));
  jetCorr->SetInputName    (pubJet->GetOutputName());
  jetCorr->SetCorrectedName("CorrectedJets");

  Bool_t excludeDoublePrompt = kFALSE;

  if (TString(dataset).Contains("-pj")) {
    mcselmod->ExcludeProcess(18);
    mcselmod->ExcludeProcess(114);
    excludeDoublePrompt = kTRUE;
  }
  if (TString(dataset).Contains("-qcd2em") || TString(dataset).Contains("-qcd-2em")) {
    excludeDoublePrompt = kTRUE;
  }

  PhotonMvaMod *photReg = new PhotonMvaMod;
  photReg->SetOutputName("GoodPhotonsRegr");
  photReg->SetIsData(isData);

  PhotonPairSelector *photPresel = new PhotonPairSelector("PhotonPairSelectorPresel");
  photPresel->SetOutputName       ("GoodPhotonsPresel");
  photPresel->SetPhotonSelType    ("MITSelection");
  photPresel->SetVertexSelType    ("CiCMVASelection");
  photPresel->DoMCSmear           (kTRUE);
  photPresel->DoDataEneCorr       (kTRUE);
  photPresel->SetPhotonsFromBranch(kFALSE);
  photPresel->SetInputPhotonsName (photReg->GetOutputName());
  photPresel->SetMCSmearFactors   (0.0045,0.0084,0.0109,0.0156,0.0203,0.0303,0.0326,0.0318,0.0331);
  photPresel->AddEnCorrPerRun     (160000,167913,0.9905,0.9905,0.9971,0.9976,1.0094,0.9994,1.0044,0.9968,1.0079);
  photPresel->AddEnCorrPerRun     (170249,172619,0.9909,0.9909,0.9975,0.9994,1.0112,0.9962,1.0012,0.9962,1.0072);
  photPresel->AddEnCorrPerRun     (172620,173692,0.9909,0.9909,0.9975,0.9977,1.0096,0.9963,1.0013,0.9947,1.0057);
  photPresel->AddEnCorrPerRun     (175860,177139,0.9911,0.9911,0.9977,0.9990,1.0109,0.9922,0.9973,0.9967,1.0077);
  photPresel->AddEnCorrPerRun     (177140,178421,0.9910,0.9910,0.9975,0.9987,1.0105,0.9921,0.9972,0.9975,1.0085);
  photPresel->AddEnCorrPerRun     (178424,999999,0.9903,0.9903,0.9969,0.9976,1.0095,0.9889,0.9940,0.9976,1.0086);
  photPresel->SetIsData           (isData);

  photPresel->SetJetsName         (jetCorr->GetOutputName());

  PhotonTreeWriter *photTreePresel = new PhotonTreeWriter("PhotonTreeWriterPresel");
  photTreePresel->SetPhotonsFromBranch  (kFALSE);
  photTreePresel->SetInputPhotonsName   (photPresel->GetOutputName());
  photTreePresel->SetEnableJets         (kTRUE);
  photTreePresel->SetPFJetsFromBranch   (kFALSE);
  photTreePresel->SetPFJetName          (jetCorr->GetOutputName());
  photTreePresel->SetExcludeDoublePrompt(excludeDoublePrompt);
  photTreePresel->SetIsData             (isData);

  //PhotonIDMod      *photIdCic = new PhotonIDMod("PhotonIDModPresel");
  //photIdCic->SetPtMin            (25.0);
  //photIdCic->SetOutputName       ("GoodPhotonsPreselid");
  //photIdCic->SetOutputName       ("MITSelection");
  //photIdCic->SetApplyElectronVeto(kTRUE);
  //photIdCic->SetIsData           (isData);
  //
  //PhotonTreeWriter *photTreeSingle = new PhotonTreeWriter("PhotonTreeWriterSingle");
  //photTreeSingle->SetWriteDiphotonTree(kFALSE);
  //photTreeSingle->SetPhotonsFromBranch(kFALSE);
  //photTreeSingle->SetInputPhotonsName (photIdCic->GetOutputName());
  //photTreeSingle->SetIsData           (isData);
  //
  //PhotonTreeWriter *photTreeE = new PhotonTreeWriter("PhotonTreeWriterE");
  //photTreeE->SetGoodElectronsFromBranch(kFALSE);
  //photTreeE->SetGoodElectronName       (elecId->GetOutputName());
  //photTreeE->SetLoopOnGoodElectrons    (kTRUE);
  //photTreeE->SetApplyElectronVeto      (kFALSE);
  //photTreeE->SetIsData                 (isData);

  //------------------------------------------------------------------------------------------------
  // making analysis chain
  //------------------------------------------------------------------------------------------------
  // this is how it always starts
  runLumiSel             ->Add(mcselmod);
  if (TString(dataset).Contains("-h"))
    mcselmod             ->Add(sysMod);
  // high level trigger is always first

  // something is not write with the HLT modules... need to fix
  //mcselmod               ->Add(hltModE);
  //mcselmod               ->Add(hltModP);
  //hltModP                ->Add(goodPVFilterMod);
  //hltModE                ->Add(goodPVFilterModE);

  mcselmod               ->Add(goodPVFilterMod);
  goodPVFilterMod        ->Add(photReg);
  photReg                ->Add(pubJet);
  pubJet                 ->Add(jetCorr);
  jetCorr                ->Add(photPresel);
  photPresel             ->Add(photTreePresel);

  //jetCorr                ->Add(photIdCic);
  //photIdCic              ->Add(photTreeSingle);
  //
  //mcselmod               ->Add(goodPVFilterModE);
  //goodPVFilterModE       ->Add(elecId);
  //elecId                 ->Add(photTreeE);
  
  //TFile::SetOpenTimeout(0);
  //TFile::SetCacheFileDir("./rootfilecache",kTRUE,kTRUE);
  //TFile::SetReadaheadSize(128*1024*1024);

  //------------------------------------------------------------------------------------------------
  // setup analysis
  //------------------------------------------------------------------------------------------------
  Analysis *ana = new Analysis;
  //ana->SetUseHLT(kTRUE);
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
  if (TString(skim).CompareTo("noskim") == 0)
    d = c->FindDataset(bookstr,dataset,fileset,true);
  else
    d = c->FindDataset(bookstr,skimdataset.Data(),fileset,true);
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

//--------------------------------------------------------------------------------------------------
int decodeEnv(char* json, char* overlap, float overlapCut)
{
  if (gSystem->Getenv("MIT_PROD_JSON"))
    sprintf(json,   "%s",gSystem->Getenv("MIT_PROD_JSON"));
  else {
    printf(" JSON file was not properly defined. EXIT!\n");
    return -1;
  }
  if (gSystem->Getenv("MIT_PROD_OVERLAP")) {
    sprintf(overlap,"%s",gSystem->Getenv("MIT_PROD_OVERLAP"));
    if (EOF == sscanf(overlap,"%f",&overlapCut)) {
      printf(" Overlap was not properly defined. EXIT!\n");
      return -1;
    }
  }
  else {
    printf(" OVERLAP file was not properly defined. EXIT!\n");
    return -1;
  }

  return 0;
}
