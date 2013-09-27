// Modified script dumpCiC_moriond.C from Mingiming, which in turn
// is a modified script from Josh for signal modeling

#include <stdio.h>
#include <iostream>
#include <sstream>
#include "TMath.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TCut.h"
#include "TTree.h"
#include "TNtuple.h"
#include "RooDataSet.h"
#include "RooExponential.h"
#include "RooLandau.h"
#include "RooPlot.h"
#include "RooFit.h"
#include "RooAddPdf.h"
#include "RooGaussian.h"
#include "RooCBShape.h"
#include "RooFFTConvPdf.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooHistFunc.h"
#include "RooMoment.h"
#include "RooFitResult.h"
#include "RooExtendPdf.h"
#include "RooGenericPdf.h"
#include "RooBreitWigner.h"
#include "RooBifurGauss.h"
#include "RooProdPdf.h"      
#include "TArrow.h"

#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooWorkspace.h"
#include "RooConstVar.h"
#include "TEfficiency.h"
#include "RooConstVar.h"
#include "RooAddition.h"

#include "TTreeFormula.h"
#include "TLatex.h"

#include "TString.h"

#include "RooHist.h"

#include "TMVA/Reader.h"
#include "TMVA/Tools.h"
#include "TMVA/Config.h"

#include "TRandom3.h"


// Declaration of functions.
//_____________________________________________________________________________
void dumpMvaInputs(bool debug = true,
                   TString fileName =
                      "/home/mingyang/cms/hist/hgg-2013Final8TeV/merged/"
                      "hgg-2013Final8TeV_s12-h120gg-vh-v7n_noskim.root");
template <class T>
void dumpVar(const char *name, T value, bool appendTab=true);

//_____________________________________________________________________________
void dumpMvaInputs(bool debug, TString fileName) {
 
  TFile* file = TFile::Open(fileName.Data());

  TDirectory* theDir = (TDirectory*)file->FindObjectAny("PhotonTreeWriterPresel");
  TTree* theTree = (TTree*) theDir->Get("hPhotonTree");
 
  // open the MVA files, if requested

 
  UInt_t run, lumi, evt;
  float ph1pt, ph1sceta, ph1iso1, ph1iso2, ph1iso3, ph1cov, ph1hoe, ph1r9;
  float ph2pt, ph2sceta, ph2iso1, ph2iso2, ph2iso3, ph2cov, ph2hoe, ph2r9;

  float ph1e, ph2e;

  float rho, mass;

  float ph1ecalIso3,ph1ecalIso4,ph1trackIsoSel03,ph1trackIsoWorst04;
  float ph2ecalIso3,ph2ecalIso4,ph2trackIsoSel03,ph2trackIsoWorst04;


  float ph1eerr, ph1eerrsmeared;
  float ph2eerr, ph2eerrsmeared;

  theTree->SetBranchAddress("run", &run);
  theTree->SetBranchAddress("lumi",&lumi);
  theTree->SetBranchAddress("evt", &evt);

  theTree->SetBranchAddress("mass",&mass);
  theTree->SetBranchAddress("rho",&rho);

  theTree->SetBranchAddress("ph1.pt",&ph1pt);
  theTree->SetBranchAddress("ph1.e",&ph1e);

  theTree->SetBranchAddress("ph1.eerr",&ph1eerr);
  theTree->SetBranchAddress("ph1.eerrsmeared",&ph1eerrsmeared);

  theTree->SetBranchAddress("ph1.sceta",&ph1sceta);
  theTree->SetBranchAddress("ph1.pfcic4_tIso1",&ph1iso1);
  theTree->SetBranchAddress("ph1.pfcic4_tIso2",&ph1iso2);
  theTree->SetBranchAddress("ph1.pfcic4_tIso3",&ph1iso3);
  theTree->SetBranchAddress("ph1.pfcic4_covIEtaIEta",&ph1cov);
  theTree->SetBranchAddress("ph1.pfcic4_HoE",&ph1hoe);
  theTree->SetBranchAddress("ph1.pfcic4_R9",&ph1r9);

  theTree->SetBranchAddress("ph1.pfcic4_ecalIso3",&ph1ecalIso3);
  theTree->SetBranchAddress("ph1.pfcic4_ecalIso4",&ph1ecalIso4);
  theTree->SetBranchAddress("ph1.pfcic4_trackIsoSel03",&ph1trackIsoSel03);
  theTree->SetBranchAddress("ph1.pfcic4_trackIsoWorst04",&ph1trackIsoWorst04);

  theTree->SetBranchAddress("ph2.pt",&ph2pt);
  theTree->SetBranchAddress("ph2.e",&ph2e);

  theTree->SetBranchAddress("ph2.eerr",&ph2eerr);
  theTree->SetBranchAddress("ph2.eerrsmeared",&ph2eerrsmeared);

  theTree->SetBranchAddress("ph2.sceta",&ph2sceta);
  theTree->SetBranchAddress("ph2.pfcic4_tIso1",&ph2iso1);
  theTree->SetBranchAddress("ph2.pfcic4_tIso2",&ph2iso2);
  theTree->SetBranchAddress("ph2.pfcic4_tIso3",&ph2iso3);
  theTree->SetBranchAddress("ph2.pfcic4_covIEtaIEta",&ph2cov);
  theTree->SetBranchAddress("ph2.pfcic4_HoE",&ph2hoe);
  theTree->SetBranchAddress("ph2.pfcic4_R9",&ph2r9);

  theTree->SetBranchAddress("ph2.pfcic4_ecalIso3",&ph2ecalIso3);
  theTree->SetBranchAddress("ph2.pfcic4_ecalIso4",&ph2ecalIso4);
  theTree->SetBranchAddress("ph2.pfcic4_trackIsoSel03",&ph2trackIsoSel03);
  theTree->SetBranchAddress("ph2.pfcic4_trackIsoWorst04",&ph2trackIsoWorst04);
  
  float jet1pt, jet2pt, jet1eta, jet2eta, dijetmass, zeppenfeld, dphidijetgg;
  float nnn = -99.;

  theTree->SetBranchAddress("jet1pt",&jet1pt);
  theTree->SetBranchAddress("jet2pt",&jet2pt);
  theTree->SetBranchAddress("jet1eta",&jet1eta);
  theTree->SetBranchAddress("jet2eta",&jet2eta);

  theTree->SetBranchAddress("dijetmass",&dijetmass);
  theTree->SetBranchAddress("zeppenfeld",&zeppenfeld);
  theTree->SetBranchAddress("dphidijetgg",&dphidijetgg);

  int numVtx;
  theTree->SetBranchAddress("nVtx",&numVtx);

  // presel vars
  float ecalisodr03_1,ecalisodr03_2;
  theTree->SetBranchAddress("ph1.ecalisodr03",&ecalisodr03_1);
  theTree->SetBranchAddress("ph2.ecalisodr03",&ecalisodr03_2);

  float hcalisodr03_1,hcalisodr03_2;
  theTree->SetBranchAddress("ph1.hcalisodr03",&hcalisodr03_1);
  theTree->SetBranchAddress("ph2.hcalisodr03",&hcalisodr03_2);

  float trkisohollowdr03_1,trkisohollowdr03_2;
  theTree->SetBranchAddress("ph1.trkisohollowdr03",&trkisohollowdr03_1);
  theTree->SetBranchAddress("ph2.trkisohollowdr03",&trkisohollowdr03_2);

  float hovere_1,hovere_2;
  theTree->SetBranchAddress("ph1.hovere",&hovere_1);
  theTree->SetBranchAddress("ph2.hovere",&hovere_2);

  float sieie_1,sieie_2;
  theTree->SetBranchAddress("ph1.sigietaieta",&sieie_1);
  theTree->SetBranchAddress("ph2.sigietaieta",&sieie_2);

  float idmva_ChargedIso_presel_1,idmva_ChargedIso_presel_2;
  theTree->SetBranchAddress("ph1.idmva_ChargedIso_selvtx",&idmva_ChargedIso_presel_1);
  theTree->SetBranchAddress("ph2.idmva_ChargedIso_selvtx",&idmva_ChargedIso_presel_2);


  float idmva_ChargedIso_0p2_presel_1,idmva_ChargedIso_0p2_presel_2;
  theTree->SetBranchAddress("ph1.idmva_ChargedIso_0p2_selvtx",&idmva_ChargedIso_0p2_presel_1);
  theTree->SetBranchAddress("ph2.idmva_ChargedIso_0p2_selvtx",&idmva_ChargedIso_0p2_presel_2);


//   float eleVeto_1, eleVeto_2;
//   theTree->SetBranchAddress("ph1.eleVeto",&eleVeto_1);
//   theTree->SetBranchAddress("ph2.eleVeto",&eleVeto_2);
  
  float idmva_1, idmva_2;
  theTree->SetBranchAddress("ph1.idmva",&idmva_1);
  theTree->SetBranchAddress("ph2.idmva",&idmva_2);

  float teta1, teta2;
  theTree->SetBranchAddress("ph1.eta",&teta1);
  theTree->SetBranchAddress("ph2.eta",&teta2);

  float vtxprob, ptgg, phi1, phi2, masserr, masserrwvtx, masserr_ns, masserrwvtx_ns;
  theTree->SetBranchAddress("vtxprob",&vtxprob);
  theTree->SetBranchAddress("ptgg",&ptgg);
  theTree->SetBranchAddress("ph1.phi",&phi1);
  theTree->SetBranchAddress("ph2.phi",&phi2);

  theTree->SetBranchAddress("masserrsmeared",&masserr);
  theTree->SetBranchAddress("masserrsmearedwrongvtx",&masserrwvtx);

  theTree->SetBranchAddress("masserr",&masserr_ns);
  theTree->SetBranchAddress("masserrwrongvtx",&masserrwvtx_ns);

  int vbfTag, leptonTag;
  theTree->SetBranchAddress("vbfTag",&vbfTag);
  theTree->SetBranchAddress("leptonTag",&leptonTag);

  float corrpfmet, corrpfmetphi, pfmet, pfmetphi;
  theTree->SetBranchAddress("corrpfmet",&corrpfmet);
  theTree->SetBranchAddress("corrpfmetphi",&corrpfmetphi);
  theTree->SetBranchAddress("pfmet",&pfmet);
  theTree->SetBranchAddress("pfmetphi",&pfmetphi);

  float mupt, mueta, mudr1, mudr2, mudz, mud0;
  theTree->SetBranchAddress("muonPt",&mupt);
  theTree->SetBranchAddress("muonEta",&mueta);
  theTree->SetBranchAddress("muDR1",&mudr1);
  theTree->SetBranchAddress("muDR2",&mudr2);
  theTree->SetBranchAddress("muD0",&mud0);
  theTree->SetBranchAddress("muDZ",&mudz);


  float elept, eleeta, elesceta, eledr1, eledr2, eleMass1, eleMass2;
  theTree->SetBranchAddress("elePt",&elept);
  theTree->SetBranchAddress("eleEta",&eleeta);
  theTree->SetBranchAddress("eleSCEta",&elesceta);
  theTree->SetBranchAddress("eleDR1",&eledr1);
  theTree->SetBranchAddress("eleDR2",&eledr2);
  theTree->SetBranchAddress("eleMass1",&eleMass1);
  theTree->SetBranchAddress("eleMass2",&eleMass2);

  int elemisshits;
  theTree->SetBranchAddress("eleNinnerHits",&elemisshits);


  // vertext stuff
  int vertexId1, vertexId2, vertexId3;
  float vertexMva1, vertexMva2, vertexMva3;
  //float	ptbal, ptasym, sumpt2, p2conv;
  Float_t	ptbal, ptasym, sumpt2, p2conv;
  int nleg1, nleg2;

  theTree->SetBranchAddress("vtxInd1",&vertexId1);
  theTree->SetBranchAddress("vtxInd2",&vertexId2);
  theTree->SetBranchAddress("vtxInd3",&vertexId3);

  theTree->SetBranchAddress("vtxMva1",&vertexMva1);
  theTree->SetBranchAddress("vtxMva2",&vertexMva2);
  theTree->SetBranchAddress("vtxMva3",&vertexMva3);

  theTree->SetBranchAddress("vtxBestPtbal",&ptbal);
  theTree->SetBranchAddress("vtxBestPtasym",&ptasym);
  theTree->SetBranchAddress("vtxBestSumpt2",&sumpt2);
  theTree->SetBranchAddress("vtxBestP2Conv",&p2conv);

  theTree->SetBranchAddress("vtxNleg1",&nleg1);
  theTree->SetBranchAddress("vtxNleg2",&nleg2);


  float convVtxZ1, convVtxRes1, convVtxChi1;
  float convVtxZ2, convVtxRes2, convVtxChi2;

  int convVtxIdx1, convVtxIdx2, vtxNconv;

  theTree->SetBranchAddress("vtxConv1Z",&convVtxZ1);
  theTree->SetBranchAddress("vtxConv1DZ",&convVtxRes1);
  theTree->SetBranchAddress("vtxConv1Prob",&convVtxChi1);
  theTree->SetBranchAddress("vtxConvIdx1",&convVtxIdx1);

  theTree->SetBranchAddress("vtxConv2Z",&convVtxZ2);
  theTree->SetBranchAddress("vtxConv2DZ",&convVtxRes2);
  theTree->SetBranchAddress("vtxConv2Prob",&convVtxChi2);
  theTree->SetBranchAddress("vtxConvIdx2",&convVtxIdx2);

  theTree->SetBranchAddress("vtxNconv",&vtxNconv);

  float vtxz1, vtxz2, vtxz3;
  theTree->SetBranchAddress("vtxMva1Z",&vtxz1);
  theTree->SetBranchAddress("vtxMva2Z",&vtxz2);
  theTree->SetBranchAddress("vtxMva3Z",&vtxz3);

  float eleIdMva;
  theTree->SetBranchAddress("eleIdMva",&eleIdMva);

  // additional MET tag stuff.
  float phigg, jetleadNoIDpt, jetleadNoIDphi, jetleadNoIDeta;
  float ph1scphi, ph2scphi;

  theTree->SetBranchAddress("phigg",&phigg);
  theTree->SetBranchAddress("jetleadNoIDpt",&jetleadNoIDpt);
  theTree->SetBranchAddress("jetleadNoIDphi",&jetleadNoIDphi);
  theTree->SetBranchAddress("jetleadNoIDeta",&jetleadNoIDeta);
  theTree->SetBranchAddress("ph1.scphi",&ph1scphi);
  theTree->SetBranchAddress("ph2.scphi",&ph2scphi);


  float ph1scrawe, ph1scpse, ph1e3x3, ph1e5x5, ph1e3x3seed, ph1e5x5seed;
  float ph2scrawe, ph2scpse, ph2e3x3, ph2e5x5, ph2e3x3seed, ph2e5x5seed;
  

  theTree->SetBranchAddress("ph1.scrawe",&ph1scrawe);
  theTree->SetBranchAddress("ph1.scpse",&ph1scpse);
  theTree->SetBranchAddress("ph1.e3x3",&ph1e3x3);
  theTree->SetBranchAddress("ph1.e3x3seed",&ph1e3x3seed);
  theTree->SetBranchAddress("ph1.e5x5",&ph1e5x5);
  theTree->SetBranchAddress("ph1.e5x5seed",&ph1e5x5seed);


  theTree->SetBranchAddress("ph2.scrawe",&ph2scrawe);
  theTree->SetBranchAddress("ph2.scpse",&ph2scpse);
  theTree->SetBranchAddress("ph2.e3x3",&ph2e3x3);
  theTree->SetBranchAddress("ph2.e3x3seed",&ph2e3x3seed);
  theTree->SetBranchAddress("ph2.e5x5",&ph2e5x5);
  theTree->SetBranchAddress("ph2.e5x5seed",&ph2e5x5seed);


  TRandom3 rng(0);

  int eventCounter=0;

  // Loop over the entries.
  for (int i=0; i < theTree->GetEntries(); ++i) {
   
    if (eventCounter > 9 && debug ) break;
    
    theTree->GetEntry(i);

    bool passPreselection = (mass < 100 &&
                             mass < 180 &&
                             ph1pt > mass/3. &&
                             ph2pt > 100./4);

    if (passPreselection == false) continue;

    eventCounter++;

    // Event Variables
    dumpVar("run"                    , run                    ); //  1
    dumpVar("lumi"                   , lumi                   ); //  2
    dumpVar("event"                  , evt                    ); //  3
    dumpVar("rho"                    , rho                    ); //  4

    // Leading Photon Variables
    dumpVar("pho1_ind"               , -999                   ); //  5
    dumpVar("pho1_scInd"             , -999                   ); //  6
    dumpVar("pho1_pt"                , ph1pt                  ); //  7
    dumpVar("pho1_eta"               , teta1                  ); //  8
    dumpVar("pho1_phi"               , phi1                   ); //  9
    dumpVar("pho1_e"                 , ph1e                   ); // 10
    dumpVar("pho1_eErr"              , ph1eerr                ); // 11
//     dumpVar("pho1_isConv"            , pho1_isConv            ); // 12
//     dumpVar("pho1_HoE"               , pho1_HoE               ); // 13
//     dumpVar("pho1_hcalIso03"         , pho1_hcalIso03         ); // 14
//     dumpVar("pho1_trkIso03"          , pho1_trkIso03          ); // 15
//     dumpVar("pho1_pfChargedIsoGood02", pho1_pfChargedIsoGood02); // 16
//     dumpVar("pho1_pfChargedIsoGood03", pho1_pfChargedIsoGood03); // 17
//     dumpVar("pho1_pfChargedIsoBad03" , pho1_pfChargedIsoBad03 ); // 18
//     dumpVar("pho1_pfPhotonIso03"     , pho1_pfPhotonIso03     ); // 19
//     dumpVar("pho1_pfNeutralIso03"    , pho1_pfNeutralIso03    ); // 20
//     dumpVar("pho1_sieie"             , pho1_sieie             ); // 21
//     dumpVar("pho1_sieip"             , pho1_sieip             ); // 22
//     dumpVar("pho1_etaWidth"          , pho1_etaWidth          ); // 23
//     dumpVar("pho1_phiWidth"          , pho1_phiWidth          ); // 24
//     dumpVar("pho1_r9"                , pho1_r9                ); // 25
//     dumpVar("pho1_lambdaRatio"       , pho1_lambdaRatio       ); // 26
//     dumpVar("pho1_s4Ratio"           , pho1_s4Ratio           ); // 27
//     dumpVar("pho1_scEta"             , pho1_scEta             ); // 28
//     dumpVar("pho1_ESEffSigmaRR"      , pho1_ESEffSigmaRR      ); // 29
//     dumpVar("pho1_ptOverM"           , pho1_ptOverM           ); // 30
// 
//     // Trailing Photon Variables
//     dumpVar("pho2_ind"               , pho2_ind               ); // 31
//     dumpVar("pho2_scInd"             , pho2_scInd             ); // 32
//     dumpVar("pho2_pt"                , pho2_pt                ); // 33
//     dumpVar("pho2_eta"               , pho2_eta               ); // 34
//     dumpVar("pho2_phi"               , pho2_phi               ); // 35
//     dumpVar("pho2_e"                 , pho2_e                 ); // 36
//     dumpVar("pho2_eErr"              , pho2_eErr              ); // 37
//     dumpVar("pho2_isConv"            , pho2_isConv            ); // 38
//     dumpVar("pho2_HoE"               , pho2_HoE               ); // 39
//     dumpVar("pho2_hcalIso03"         , pho2_hcalIso03         ); // 40
//     dumpVar("pho2_trkIso03"          , pho2_trkIso03          ); // 41
//     dumpVar("pho2_pfChargedIsoGood02", pho2_pfChargedIsoGood02); // 42
//     dumpVar("pho2_pfChargedIsoGood03", pho2_pfChargedIsoGood03); // 43
//     dumpVar("pho2_pfChargedIsoBad03" , pho2_pfChargedIsoBad03 ); // 44
//     dumpVar("pho2_pfPhotonIso03"     , pho2_pfPhotonIso03     ); // 45
//     dumpVar("pho2_pfNeutralIso03"    , pho2_pfNeutralIso03    ); // 46
//     dumpVar("pho2_sieie"             , pho2_sieie             ); // 47
//     dumpVar("pho2_sieip"             , pho2_sieip             ); // 48
//     dumpVar("pho2_etaWidth"          , pho2_etaWidth          ); // 49
//     dumpVar("pho2_phiWidth"          , pho2_phiWidth          ); // 50
//     dumpVar("pho2_r9"                , pho2_r9                ); // 51
//     dumpVar("pho2_lambdaRatio"       , pho2_lambdaRatio       ); // 52
//     dumpVar("pho2_s4Ratio"           , pho2_s4Ratio           ); // 53
//     dumpVar("pho2_scEta"             , pho2_scEta             ); // 54
//     dumpVar("pho2_ESEffSigmaRR"      , pho2_ESEffSigmaRR      ); // 55
//     dumpVar("pho2_ptOverM"           , pho2_ptOverM           ); // 56

    // Diphoton Variables
    dumpVar("mass"                   , mass                   ); // 57
    dumpVar("rVtxSigmaMoM"           , masserr_ns             ); // 58
    dumpVar("wVtxSigmaMoM"           , masserrwvtx_ns         ); // 59
    dumpVar("vtxIndex"               , vertexId1              ); // 60
    dumpVar("vtxProb"                , vtxprob                ); // 61
    dumpVar("ptBal"                  , ptbal                  ); // 62
    dumpVar("ptAsym"                 , ptasym                 ); // 63
    dumpVar("logSPt2"                , TMath::Log(sumpt2)     ); // 64
    dumpVar("p2Conv"                 , p2conv                 ); // 65
    dumpVar("nConv"                  , vtxNconv               ); // 66
    dumpVar("cosDPhi"                , TMath::Cos(phi1-phi2)  ); // 67

//     // Leading Jet Variables
//     dumpVar("jet1_ind"               , jet1_ind               ); // 68
//     dumpVar("jet1_pt"                , jet1_pt                ); // 69
//     dumpVar("jet1_eta"               , jet1_eta               ); // 70
// 
//     // Trailing Jet Variables
//     dumpVar("jet2_ind"               , jet2_ind               ); // 71
//     dumpVar("jet2_pt"                , jet2_pt                ); // 72
//     dumpVar("jet2_eta"               , jet2_eta               ); // 73
// 
//     // Dijet Variables
//     dumpVar("dijet_dEta"             , dijet_dEta             ); // 74
//     dumpVar("dijet_Zep"              , dijet_Zep              ); // 75
//     dumpVar("dijet_dPhi"             , dijet_dPhi             ); // 76
//     dumpVar("dijet_Mjj"              , dijet_Mjj       , false); // 77

    std::cout << std::endl;
  } // Loop over the tree entries.
  
  return;

} // void dumpMvaInputs(bool debug, TString fileName)


//_____________________________________________________________________________
template <class T>
void dumpVar(const char *name, T value, bool appendTab) {
  std::cout << name << ":" << value;
  if (appendTab) {
    std::cout << "\t";
  }
}
