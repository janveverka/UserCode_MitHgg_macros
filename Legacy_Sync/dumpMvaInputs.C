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

void dumpMvaInputs(bool debug = true, TString fileName = "/home/mingyang/cms/hist/hgg-2013Final8TeV/merged/hgg-2013Final8TeV_s12-h120gg-vh-v7n_noskim.root") {
 
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

  int convVtxIdx1, convVtxIdx2;

  theTree->SetBranchAddress("vtxConv1Z",&convVtxZ1);
  theTree->SetBranchAddress("vtxConv1DZ",&convVtxRes1);
  theTree->SetBranchAddress("vtxConv1Prob",&convVtxChi1);
  theTree->SetBranchAddress("vtxConvIdx1",&convVtxIdx1);

  theTree->SetBranchAddress("vtxConv2Z",&convVtxZ2);
  theTree->SetBranchAddress("vtxConv2DZ",&convVtxRes2);
  theTree->SetBranchAddress("vtxConv2Prob",&convVtxChi2);
  theTree->SetBranchAddress("vtxConvIdx2",&convVtxIdx2);

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

  int evCounter=0;
  
  for(int i=0; i<theTree->GetEntries(); ++i) {
   
    theTree->GetEntry(i);
    
    if( mass > 100 && mass < 180  && ph1pt > mass/3. && ph2pt > 100./4) {
      
      if(evCounter > 9 && debug ) break;

      vbfTag = 0;
      
      if (jet1pt > 30.0 &&
          jet2pt > 30.0 &&
          abs(jet1eta - jet2eta) > 3. &&
          dijetmass > 500.0 &&
          zeppenfeld < 2.5 &&
          abs(dphidijetgg) > 2.6 &&
          ph1pt > (60.0 * mass / 120.0) &&
          ph2pt > 25.) {
        vbfTag = 2;
      } else if (jet1pt > 30.0 &&
                 jet2pt > 20.0 &&
                 abs(jet1eta - jet2eta) > 3. &&
                 dijetmass > 250.0 &&
                 zeppenfeld < 2.5 &&
                 abs(dphidijetgg) > 2.6 &&
                 !(jet2pt > 30.0 && dijetmass > 500.0) &&
                 ph1pt > (60.0 * mass / 120.0) &&
                 ph2pt > 25.) {
        vbfTag = 1;
      }

      // test for evt Cat
      /*int evCat = -1;
      if ( leptonTag > 1.5 && leptonTag < 2.5 && ph2pt > mass/4) evCat = 6;
      else if ( leptonTag > 0.5 && leptonTag < 1.5 && ph2pt > mass/4) evCat = 7;
      else if ( vbfTag > 0.5 && vbfTag < 1.5) evCat = 5;
      else if ( vbfTag > 1.5 ) evCat = 4;
      else if ( TMath::Abs(ph1sceta) < 1.4442 && TMath::Abs(ph2sceta) < 1.4442 && corrpfmet > 70. && ph1pt/mass > 45./120. &&
		TMath::ACos(TMath::Cos(phigg-corrpfmetphi)) > 2.1 && 
		( jetleadNoIDpt < 50. || 
		  TMath::Sqrt( TMath::Power(ph1sceta-jetleadNoIDeta,2)+TMath::Power(TMath::ACos(TMath::Cos(ph1scphi-jetleadNoIDphi)),2)) < 0.5 || TMath::Sqrt( TMath::Power((ph2sceta-jetleadNoIDeta),2)+TMath::Power(TMath::ACos(TMath::Cos(ph2scphi-jetleadNoIDphi)),2)) < 0.5 || TMath::ACos(TMath::Cos(TMath::Abs(jetleadNoIDphi-corrpfmetphi))) < 2.7 ) && ph2pt > mass/4) evCat = 8;
      else if ( TMath::Abs(ph1sceta) < 1.5 &&  TMath::Abs(ph2sceta) < 1.5 && ph1r9 > 0.94 && ph2r9 > 0.94 && ph2pt > mass/4) evCat = 0;
      else if ( TMath::Abs(ph1sceta) < 1.5 &&  TMath::Abs(ph2sceta) < 1.5 && !(ph1r9 > 0.94 && ph2r9 > 0.94) && ph2pt > mass/4) evCat = 1;
      else if ( !(TMath::Abs(ph1sceta) < 1.5 &&  TMath::Abs(ph2sceta) < 1.5) && ph1r9 > 0.94 && ph2r9 > 0.94 && ph2pt > mass/4) evCat = 2;
      else if ( !(TMath::Abs(ph1sceta) < 1.5 &&  TMath::Abs(ph2sceta) < 1.5) && !(ph1r9 > 0.94 && ph2r9 > 0.94) && ph2pt > mass/4) evCat = 3;*/

      int evCat = -1;
      if ( leptonTag > 1.5 && leptonTag < 2.5 && ph2pt > 100./4 && ph1pt/mass > 45./120.) evCat = 6;
      else if ( leptonTag > 0.5 && leptonTag < 1.5 && ph2pt > 100./4  && ph1pt/mass > 45./120.) evCat = 7;
      else if ( vbfTag > 0.5 && vbfTag < 1.5) evCat = 5;
      else if ( vbfTag > 1.5 ) evCat = 4;
      else if ( TMath::Abs(ph1sceta) < 1.4442 && TMath::Abs(ph2sceta) < 1.4442 && corrpfmet > 70. && ph1pt/mass > 45./120. &&
		TMath::ACos(TMath::Cos(phigg-corrpfmetphi)) > 2.1 && 
		( jetleadNoIDpt < 50. || 
		  TMath::Sqrt( TMath::Power(ph1sceta-jetleadNoIDeta,2)+TMath::Power(TMath::ACos(TMath::Cos(ph1scphi-jetleadNoIDphi)),2)) < 0.5 || TMath::Sqrt( TMath::Power((ph2sceta-jetleadNoIDeta),2)+TMath::Power(TMath::ACos(TMath::Cos(ph2scphi-jetleadNoIDphi)),2)) < 0.5 || TMath::ACos(TMath::Cos(TMath::Abs(jetleadNoIDphi-corrpfmetphi))) < 2.7 ) && ph2pt > 100./4) evCat = 8;
      else if ( TMath::Abs(ph1sceta) < 1.5 &&  TMath::Abs(ph2sceta) < 1.5 && ph1r9 > 0.94 && ph2r9 > 0.94 && ph2pt > mass/4) evCat = 0;
      else if ( TMath::Abs(ph1sceta) < 1.5 &&  TMath::Abs(ph2sceta) < 1.5 && !(ph1r9 > 0.94 && ph2r9 > 0.94) && ph2pt > mass/4) evCat = 1;
      else if ( !(TMath::Abs(ph1sceta) < 1.5 &&  TMath::Abs(ph2sceta) < 1.5) && ph1r9 > 0.94 && ph2r9 > 0.94 && ph2pt > mass/4) evCat = 2;
      else if ( !(TMath::Abs(ph1sceta) < 1.5 &&  TMath::Abs(ph2sceta) < 1.5) && !(ph1r9 > 0.94 && ph2r9 > 0.94) && ph2pt > mass/4) evCat = 3;
      
      evCounter++;

      if(evCat != -1){
      //if(true){

      printf("type:0\trun:%d\tlumi:%d\tevent:%d\trho:%.4e\tr9_1:%.4e\tsceta_1:%.4e\thoe_1:%.4e\tsigieie_1:%.4e\tecaliso_1:%.4e\thcaliso_1:%.4e\ttrckiso_1:%.4e\tchpfis2_1:%.4e\tchpfis3_1:%.4e\tr9_2:%.4e\tsceta_2:%.4e\thoe_2:%.4e\tsigieie_2:%.4e\tecaliso_2:%.4e\thcaliso_2:%.4e\ttrckiso_2:%.4e\tchpfis2_2:%.4e\tchpfiso3_2:%.4e\t",	     
	     run,lumi,evt,rho,
	     ph1r9,ph1sceta,hovere_1,sieie_1,
	     ecalisodr03_1-0.012*ph1pt,
	     hcalisodr03_1-0.005*ph1pt,
	     trkisohollowdr03_1-0.002*ph1pt,
	     idmva_ChargedIso_0p2_presel_1,idmva_ChargedIso_presel_1,
	     ph2r9,ph2sceta,hovere_2,sieie_2,
	     ecalisodr03_2-0.012*ph2pt,
	     hcalisodr03_2-0.005*ph2pt,
	     trkisohollowdr03_2-0.002*ph2pt,
	     idmva_ChargedIso_0p2_presel_2,idmva_ChargedIso_presel_2);

      float bsw=4.8;

      printf("ptH:%.4e\tphoid_1:%.4e\tphoid_2:%.4e\tphoeta_1:%.4e\tphoeta_2:%.4e\tsigmrv:%.4e\tbsw:%.4e\tsigmwv:%.4e\tpt_1/m:%.4e\tpt_2/m:%.4e\tvtxprob:%.4e\tcosdphi:%.4e\tmgg:%.4e\te_1:%.4e\te_2:%.4e\teerr_1:%.4e\teerr_2:%.4e\teerrsmeared_1:%.4e\teerrsmeared_2:%.4e\tvbfevt:%d\tmuontag:%d\telectrontag:%d\tmettag:%d\tevcat:%d\t",
	     ptgg,
	     idmva_1, idmva_2,
	     teta1,teta2,
	     masserr,bsw,masserrwvtx,
	     ph1pt/mass, ph2pt/mass,
	     vtxprob,TMath::Cos(phi1-phi2),
	     mass,
	     ph1e, ph2e,
	     ph1eerr,ph2eerr,
	     ph1eerrsmeared,ph2eerrsmeared,
	     (int) (vbfTag > 0.5),
	     (int) (leptonTag > 1.5 && leptonTag < 2.5),
	     (int) (leptonTag > 0.5 && leptonTag < 1.5),
	     (int) (TMath::Abs(ph1sceta) < 1.4442 && TMath::Abs(ph2sceta) < 1.4442 && corrpfmet > 70. && ph1pt/mass > 45./120. &&
		TMath::ACos(TMath::Cos(phigg-corrpfmetphi)) > 2.1 && 
		( jetleadNoIDpt < 50. || 
		  TMath::Sqrt( TMath::Power(ph1sceta-jetleadNoIDeta,2)+TMath::Power(TMath::ACos(TMath::Cos(ph1scphi-jetleadNoIDphi)),2)) < 0.5 || TMath::Sqrt( TMath::Power((ph2sceta-jetleadNoIDeta),2)+TMath::Power(TMath::ACos(TMath::Cos(ph2scphi-jetleadNoIDphi)),2)) < 0.5 || TMath::ACos(TMath::Cos(TMath::Abs(jetleadNoIDphi-corrpfmetphi))) < 2.7 ) ), evCat);
      
      if ( ! ( TMath::Abs(ph1sceta) < 1.4442 && TMath::Abs(ph2sceta) < 1.4442 && corrpfmet > 70. && ph1pt/mass > 45./120. &&
		TMath::ACos(TMath::Cos(phigg-corrpfmetphi)) > 2.1 && 
		( jetleadNoIDpt < 50. || 
		  TMath::Sqrt( TMath::Power(ph1sceta-jetleadNoIDeta,2)+TMath::Power(TMath::ACos(TMath::Cos(ph1scphi-jetleadNoIDphi)),2)) < 0.5 || TMath::Sqrt( TMath::Power((ph2sceta-jetleadNoIDeta),2)+TMath::Power(TMath::ACos(TMath::Cos(ph2scphi-jetleadNoIDphi)),2)) < 0.5 || TMath::ACos(TMath::Cos(TMath::Abs(jetleadNoIDphi-corrpfmetphi))) < 2.7 ) ) ) {
	corrpfmet = -1.;
	corrpfmetphi = -1.;
	pfmet = -1.;
	pfmetphi = -1.;
      }

      if ( ! (leptonTag < 2.5 && leptonTag > 1.5 )) {
	mupt = -1.;
	mueta = -1.;
	mudr1=-1.;
	mudr2=-1.;
	mudz=-1.;
	mud0=-1.;
      }

      if ( !(leptonTag < 1.5 && leptonTag > 0.5) ) {
	elept = -1.;
	eleeta = -1.;
	elesceta = -1.;
	eledr1 = -1.;
	eledr2 = -1.;
	eleMass1 = -1.;
	eleMass2 = -1.;
	elemisshits = -1;
      }


     /******** printf("FileName:DUMMY.root\tvertexId1:%d\tvertexId2:%d\tvertexId3:%d\tvertexMva1:%.4f\tvertexMva2:%.4f\tvertexMva3:%.4f\tvertexdeltaz2:%.4f\tvertexdeltaz3:%.4f\tptbal:%.4f\tptasym:%.4f\tlogspt2:%.4f\tp2conv:%.4f\tconvVtxZ1:%.4f\tconvVtxdZ1:%.4f\tconvRes1:%.4f\tconvChiProb1:%.4f\tconvNtrk1:%d\tconvindex1:%d\tconvVtxZ2:%.4f\tconvVtxdZ2:%.4f\tconvRes2:%.4f\tconvChiProb2:%.4f\tconvNtrk2:%d\tconvindex2:%d\t",
	     vertexId1, vertexId2, vertexId3,
	     vertexMva1, vertexMva2,vertexMva3,
	     vtxz2-vtxz1, vtxz3-vtxz1, 
	     ptbal, ptasym, TMath::Log(sumpt2),p2conv,
	     convVtxZ1, (nleg1 > 0 ? convVtxZ1-vtxz1 : -99.), convVtxRes1, convVtxChi1,
	     nleg1, convVtxIdx1,
	     convVtxZ2, (nleg2 > 0 ? convVtxZ2-vtxz1 : -99.), convVtxRes2, convVtxChi2,
	     nleg2, convVtxIdx2
	     );*/

      printf("FileName:DUMMY.root\tvertexId1:%d\tvertexId2:%d\tvertexId3:%d\tvertexMva1:%f\tvertexMva2:%f\tvertexMva3:%f\tvertexdeltaz2:%f\tvertexdeltaz3:%f\tptbal:%f\tptasym:%5f\tlogspt2:%f\tp2conv:%f\tconvVtxZ1:%f\tconvVtxdZ1:%f\tconvRes1:%f\tconvChiProb1:%f\tconvNtrk1:%d\tconvindex1:%d\tconvVtxZ2:%f\tconvVtxdZ2:%f\tconvRes2:%f\tconvChiProb2:%.4f\tconvNtrk2:%d\tconvindex2:%d\t",
	     vertexId1, vertexId2, vertexId3,
	     vertexMva1, vertexMva2,vertexMva3,
	     vtxz2-vtxz1, vtxz3-vtxz1, 
	     ptbal, ptasym, TMath::Log(sumpt2),p2conv,
	     convVtxZ1, (nleg1 > 0 ? convVtxZ1-vtxz1 : -99.), convVtxRes1, convVtxChi1,
	     nleg1, convVtxIdx1,
	     convVtxZ2, (nleg2 > 0 ? convVtxZ2-vtxz1 : -99.), convVtxRes2, convVtxChi2,
	     nleg2, convVtxIdx2
	     );

      
      printf("\tmuind:%d\tmupt:%.4f\tmueta:%.4f\tmuiso:%.4f\tmuisopt:%.4f\tmud0:%.4f\tmudq:%.4f\tmudr1:%.4f\tmudr2:%.4f\telpt:%.4f\teleta:%.4f\telsceta:%.4f\telmva:%.4f\teliso:%.4f\telisoopt:%.4f\telAeff:%.4f\teld0:%.4f\teldZ:%.4f\telmishits:%d\teldr1:%.4f\teldr2:%.4f\telmeg1:%.4f\telmeg2:%.4f\tmetuncor:%.4f\tmetphiuncor:%.4f\tmetcor:%.4f\tmetphicor:%.4f\t",
	     (int)nnn, mupt, mueta, nnn, nnn, mud0, mudz, mudr1, mudr2,
	     elept, eleeta, elesceta, eleIdMva, nnn, nnn, nnn, nnn, nnn, elemisshits, eledr1, eledr2, eleMass1, eleMass2,
	     pfmet, pfmetphi, corrpfmet, corrpfmetphi);

      printf("\tph1scrawe:%.4f\tph1scpse:%.4f\tph1sce3x3:%.4f\tph1sce3x3seed:%.4f\tph1sce5x5:%.4f\tph1sce5x5seed:%.4f\tph2scrawe:%.4f\tph2scpse:%.4f\tph2sce3x3:%.4f\tph2sce3x3seed:%.4f\tph2sce5x5:%.4f\tph2sce5x5seed:%.4f\n",
	     ph1scrawe, ph1scpse, ph1e3x3, ph1e3x3seed, ph1e5x5, ph1e5x5seed,
	     ph2scrawe, ph2scpse, ph2e3x3, ph2e3x3seed, ph2e5x5, ph2e5x5seed);
      }
      
    }
    
  }
  
  return;
    
}
