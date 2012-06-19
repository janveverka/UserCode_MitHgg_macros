// MIT BDT Training for H->gg analysis

#include "TMVA/Reader.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TMath.h"
#include "TH1D.h"
#include "TTreeFormula.h"

double trigEff[36];

TH1D *puweights[12];
float puweight(float npu, int wset=0) {
  if (npu<0) return 1.0;
  return puweights[wset]->GetBinContent(puweights[wset]->FindFixBin(npu));
}

//pt-reweighing (but only for gf samples, use extra bool flag for now)
TH1D *ptweights = 0;
float ptweight(float genhpt, Int_t procid) {
  if (procid>0) return 1.0;
  if (genhpt<0) return 1.0;
  return ptweights->GetBinContent(ptweights->FindFixBin(genhpt));
}

//efficiency scale factors (barrel/endcap for now, but can be extended to full kinematic binning)
/*float effweight(int cat) {
  //lp11 numbers
//   if (cat==1) return 0.993*1.002;
//   else if (cat==2) return 1.011*1.011;
//   else if (cat==3) return 1.004*1.000;
//   else if (cat==4) return 1.049*0.996;
//   else return 1.0;


  //oct 18 mva numbers
//   if (cat==1) return 0.998*1.002;
//   else if (cat==2) return 1.019*1.011;
//   else if (cat==3) return 1.025*1.000;
//   else if (cat==4) return 1.068*0.996;
//   else return 1.0;

  //nov15 freeze numbers (old eveto numbers)
//   if (cat==1) return 0.985*1.002;
//   else if (cat==2) return 1.002*1.011;
//   else if (cat==3) return 1.002*1.000;
//   else if (cat==4) return 1.052*0.996;
//   else return 1.0;  

  //return 1.0;
  //ming:needs to be changed to the preselection category instead of cic category
  if (cat==1) return 0.999;
  else if (cat==2) return 0.984;
  else if (cat==3) return 1.006;
  else if (cat==4) return 1.014;
  else return 1.0;  

  }*/

//ming change
float effweight(int cat, float r9) {
  
  double scale=1.0;
  
  bool iseb = (cat==1 || cat==2);
  bool isr9 = r9>0.9;
  
  //ele veto for nov30 in CIC category (number from Doug)
  //if (cat==1) scale*=1.002;
  //else if (cat==2) scale*=1.011;
  //else if (cat==3) scale*=1.000;
  //else if (cat==4) scale*=0.996;

  //conversion safe ele veto efficiency Data/MC scale factor from 2012jan16 in CIC category (number from Doug)
  //https://hypernews.cern.ch/HyperNews/CMS/get/AUX/2012/04/16/15:22:40-54817-ElectronVeto.pdf slide 11
  if (cat==1) scale*=1.000;//100/100
  else if (cat==2) scale*=1.000;//99.7/99.7
  else if (cat==3) scale*=1.005;//99.3/98.8
  else if (cat==4) scale*=1.029;//98.1/95.3
  
  //2012jan16 presel scale factors (number from Matteo)
  //ming:are these number derived for PT>25?
  if (iseb&&isr9) scale*=0.999;
  else if (iseb&&!isr9) scale*=0.984;
  else if (!iseb&&isr9) scale*=1.006;
  else if (!iseb && !isr9) scale*=1.014;


  //2012 ICHEP presel scale factors (number from Matteo)
  //(https://hypernews.cern.ch/HyperNews/CMS/get/AUX/2012/05/28/14:42:22-22393-scale_factors.pdf)
  //(where take the idmva>-0.3 scale factor to be 1)
  //is this preselection efficiency calculated with respect to the photon with pt>25 GeV? I assume so.
  if (iseb&&isr9) scale*=0.999;
  else if (iseb&&!isr9) scale*=0.996;
  else if (!iseb&&isr9) scale*=1.007;
  else if (!iseb && !isr9) scale*=1.012; 

  return scale;
}

float trigeffweight(int cat1, int cat2, double pt1, double pt2) {
  Bool_t isb1 = cat1==1 || cat1==2;
  Bool_t isb2 = cat2==1 || cat2==2;
  
  Bool_t isr91 = cat1==1 || cat1==3;
  Bool_t isr92 = cat2==1 || cat2==3;
  
  Bool_t isb = isb1 && isb2;
  Bool_t isr9 = isr91 && isr92;
  
  //if (isb && isr9) return 1.0;
  //else if (isb && !isr9) return 0.993;
  //else if (!isb && isr9) return 1.0;
  //else if (!isb && !isr9) return 0.988;
  // else return 1.0;
  return 1.0;
  
}


// float trigeffweight(int cat1, int cat2, double pt1, double pt2) {
// 
//   Bool_t isb1 = cat1==1 || cat1==2;
//   Bool_t isb2 = cat2==1 || cat2==2;
//   
//   Bool_t isr91 = cat1==1 || cat1==3;
//   Bool_t isr92 = cat2==1 || cat2==3;
// 
//   int _cat1 = 1;
//   if ( !isb1 ) _cat1 += 4;
//   if ( !isr91 ) _cat1 += 2;
//   if ( pt1 < 40./120.*100.) _cat1 += 1;
// 
//   int _cat2 = 1;
//   if ( !isb2 ) _cat2 += 4;
//   if ( !isr92 ) _cat2 += 2;
//   if ( pt2 < 40./120.*100.) _cat2 += 1;
// 
//   if( _cat2 < _cat1) {
//     int tmp = _cat2;
//     _cat2 = _cat1;
//     _cat1 = tmp;
//   }
// 
//   int binNum = ( _cat1 - 1 ) * ( 18 - _cat1 )/2 + ( _cat2 - _cat1 );
// 
//   return trigEff[binNum];
// }

void setRelErrAlias(TTree* tree) {
  TString  corfact = "0.07";
  TString relerrsmearcor = TString::Format("0.5*sqrt( (pow(1.0+ismc*ph1.isbarrel*%s,2)*ph1.eerr*ph1.eerr+ph1.esmearing*ph1.esmearing)/ph1.e/ph1.e + (pow(1.0+ismc*ph2.isbarrel*%s,2)*ph2.eerr*ph2.eerr+ph2.esmearing*ph2.esmearing)/ph2.e/ph2.e)",corfact.Data(),corfact.Data());

  tree->SetAlias("massrelerr",relerrsmearcor.Data());
  return;
}

TH1D *getpuweights(TFile *file, TH1D *target) {
  
  TDirectory *dirmcpv = (TDirectory*)file->FindObjectAny("AnaFwkMod");
  TH1D *hnpu = (TH1D*)dirmcpv->Get("hNPU");
  TH1D *hpumc = (TH1D*)hnpu->Clone();
  
  hpumc->Sumw2();
  hpumc->Scale(1.0/hpumc->Integral(0,hpumc->GetNbinsX()+1));
  
  
  TH1D *htargettmp = new TH1D("htargettmp","", hpumc->GetNbinsX(), hpumc->GetXaxis()->GetXmin(), hpumc->GetXaxis()->GetXmax());
  htargettmp->Sumw2();
  for (int ibin = 0; ibin<=(htargettmp->GetNbinsX()+1); ++ibin) {
    htargettmp->Fill(htargettmp->GetBinCenter(ibin),target->GetBinContent(target->FindFixBin(htargettmp->GetBinCenter(ibin))));
  }
  htargettmp->Scale(1.0/htargettmp->Integral(0,htargettmp->GetNbinsX()+1));
  
  TH1D *puweights = new TH1D((*htargettmp)/(*hpumc));
    
  delete htargettmp;
  
  return puweights;
    
}

int eventCat(int cat1, int cat2) {
    Bool_t isb1 = cat1==1 || cat1==2;
  Bool_t isb2 = cat2==1 || cat2==2;
  
  Bool_t isr91 = cat1==1 || cat1==3;
  Bool_t isr92 = cat2==1 || cat2==3;
  
  Bool_t isb = isb1 && isb2;
  Bool_t isr9 = isr91 && isr92;
  
  if (isb && isr9) return 1;
  else if (isb && !isr9) return 2;
  else if (!isb && isr9) return 3;
  else if (!isb && !isr9) return 4;
  else return 1;
}

void applyTMVA_Bambu_MIT() {
  //trig eff calculated by Fabian
  trigEff[0] = 0.990  ;   
  trigEff[1] = 0.989  ;
  trigEff[2] = 0.977  ; 
  trigEff[3] = 0.977  ; 
  trigEff[4] = 0.987  ; 
  trigEff[5] = 0.991  ; 
  trigEff[6] = 0.968  ; 
  trigEff[7] = 0.973  ; 
  trigEff[8] = 0.986  ; 
  trigEff[9] = 0.975  ; 
  trigEff[10] = 0.975  ;
  trigEff[11] = 0.985  ;
  trigEff[12] = 0.988  ;
  trigEff[13] = 0.967  ;
  trigEff[14] = 0.969  ;
  trigEff[15] = 0.964  ;
  trigEff[16] = 0.964  ;
  trigEff[17] = 0.974  ;
  trigEff[18] = 0.977  ;
  trigEff[19] = 0.952  ;
  trigEff[20] = 0.958  ; 
  trigEff[21] = 0.963  ;
  trigEff[22] = 0.973  ;
  trigEff[23] = 0.974  ;
  trigEff[24] = 0.951  ;
  trigEff[25] = 0.956  ;
  trigEff[26] = 0.984  ;
  trigEff[27] = 0.987  ;
  trigEff[28] = 0.964  ;
  trigEff[29] = 0.969  ;
  trigEff[30] = 0.986  ;
  trigEff[31] = 0.968  ;
  trigEff[32] = 0.967  ;
  trigEff[33] = 0.942  ;
  trigEff[34] = 0.948  ;
  trigEff[35] = 0.949  ;
  
  //---------------------------------------------------------------
  // setting up the weight.
  //TFile *filepuest = new TFile("/home/mingyang/cms/puweight/augmented_nov08_rereco.json.68000.pileup.root","READ");//jan16 rereco
  //TFile *filepuest = new TFile("/home/mingyang/cms/puweight/latest_prompt_2012.json.68300.observed.pileup.root","READ");//2012 ICHEP
  TFile *filepuest = new TFile("/home/mingyang/cms/puweight/rereco_prompt_2012.json.69000.observed.pileup.root","READ");//2012 ICHEP
  TH1D *hpuest = (TH1D*) filepuest->Get("pileup");

  //higgs pt reweight of powheg sample to HQt; this needed for 2011 samples but not 2012 samples   
  //load pt-weights
  TFile *fileptweight = new TFile("/scratch/bendavid/root/KFactors_AllScales.root","READ");
  //ptweights= (TH1D*) fileptweight->Get("kfact125_0");  
  ptweights= (TH1D*) fileptweight->Get("kfact120_0"); 
  
  //---------------------------------------------------------------

  TString baseDir="RunLumiSelectionMod/MCProcessSelectionMod/HLTModP/GoodPVFilterMod/PhotonMvaMod/JetPub/JetCorrectionMod/PhotonPairSelectorPresel/PhotonTreeWriterPresel/";
  //TString baseDir="RunLumiSelectionMod/MCProcessSelectionMod/HLTModP/GoodPVFilterMod/PhotonMvaMod/JetPub/JetCorrectionMod/PhotonPairSelectorCiC/PhotonTreeWriterCiC/";
  
  
  //TString baseDir="RunLumiSelectionMod/MCProcessSelectionMod/HLTModP/GoodPVFilterMod/PhotonMvaMod/JetPub/JetCorrectionMod/PhotonPairSelector/PhotonTreeWriter/";

  /*TString signal_1 ="hgg-v0_f11--h125gg-gf-v14b-pu_noskim.root";
  TString signal_2 ="hgg-v0_f11--h125gg-vbf-v14b-pu_noskim.root";
  TString signal_3 ="hgg-v0_f11--h125gg-vh-v14b-pu_noskim.root";
  TString signal_4 ="hgg-v0_f11--h125gg-tt-v14b-pu_noskim.root";*/

  /*TString signal_1 ="hgg-v0_f11--h120gg-gf-v14b-pu_noskim.root";
  TString signal_2 ="hgg-v0_f11--h120gg-vbf-v14b-pu_noskim.root";
  TString signal_3 ="hgg-v0_f11--h120gg-vh-v14b-pu_noskim.root";
  TString signal_4 ="hgg-v0_f11--h120gg-tt-v14b-pu_noskim.root";

  TString bgBorn   ="hgg-v0_f11--diphoj-v14b-pu_noskim.root";
  TString bgBox_1  ="hgg-v0_f11--2pibx10_25-v14b-pu_noskim.root";
  TString bgBox_2  ="hgg-v0_f11--2pibx25_250-v14b-pu_noskim.root";
  TString bgBox_3  ="hgg-v0_f11--2pibx250-v14b-pu_noskim.root";*/

  /* TString signal_1 ="hggMvaUnblind31fb_s12-h124gg-gf-v9_noskim.root";
  TString signal_2 ="hggMvaUnblind31fb_s12-h124gg-vbf-v9_noskim.root";
  TString signal_3 ="hggMvaUnblind31fb_s12-h124gg-vh-v9_noskim.root";
  TString signal_4 ="hggMvaUnblind31fb_s12-h124gg-tt-v9_noskim.root";

  TString bgBorn   ="hggMvaUnblind31fb_s12-diphoj-v9_noskim.root";
  TString bgBox_1  ="hggMvaUnblind31fb_s12-2pibx10_25-v9_noskim.root";
  TString bgBox_2  ="hggMvaUnblind31fb_s12-2pibx25_250-v9_noskim.root";
  TString bgBox_3  ="hggMvaUnblind31fb_s12-2pibx250-v9_noskim.root";

  TString bgQCD_1    ="hggMvaUnblind31fb_s12-qcd-2em3040-v9_noskim.root";
  TString bgQCD_2    ="hggMvaUnblind31fb_s12-qcd-2em40-v9_noskim.root";
  TString bgPJ_1     ="hggMvaUnblind31fb_s12-pj20_40-2em-v9_noskim.root";
  TString bgPJ_2     ="hggMvaUnblind31fb_s12-pj40-2em-v9_noskim.root";*/

  TString signal_1 ="hggMvaTrain38fb_s12-h120gg-gf-v9_noskim.root";
  TString signal_2 ="hggMvaTrain38fb_s12-h120gg-vbf-v9_noskim.root";
  TString signal_3 ="hggMvaTrain38fb_s12-h120gg-vh-v9_noskim.root";
  TString signal_4 ="hggMvaTrain38fb_s12-h120gg-tt-v9_noskim.root";

  TString bgBorn   ="hggMvaTrain38fb_s12-diphoj-v9_noskim.root";
  TString bgBox_1  ="hggMvaTrain38fb_s12-2pibx10_25-v9_noskim.root";
  TString bgBox_2  ="hggMvaTrain38fb_s12-2pibx25_250-v9_noskim.root";
  TString bgBox_3  ="hggMvaTrain38fb_s12-2pibx250-v9_noskim.root";

  TString bgQCD_1    ="hggMvaTrain38fb_s12-qcd-2em3040-v9_noskim.root";
  TString bgQCD_2    ="hggMvaTrain38fb_s12-qcd-2em40-v9_noskim.root";
  TString bgPJ_1     ="hggMvaTrain38fb_s12-pj20_40-2em-v9_noskim.root";
  TString bgPJ_2     ="hggMvaTrain38fb_s12-pj40-2em-v9_noskim.root";


  //TString prefix = "/home/mingyang/cms/hist/hgg-v0/merged/";
  //TString prefix = "/home/mingyang/cms/hist/hggMvaUnblind31fb/merged/";
  TString prefix = "/scratch/bendavid/cms/hist/hggMvaTrain38fb/merged/";

  TFile* f_signal_1 = TFile::Open((prefix+signal_1).Data());
  TFile* f_signal_2 = TFile::Open((prefix+signal_2).Data());
  TFile* f_signal_3 = TFile::Open((prefix+signal_3).Data());
  TFile* f_signal_4 = TFile::Open((prefix+signal_4).Data());

  TFile* f_bgBorn   = TFile::Open((prefix+bgBorn).Data());
  TFile* f_bgBox_1  = TFile::Open((prefix+bgBox_1).Data());
  TFile* f_bgBox_2  = TFile::Open((prefix+bgBox_2).Data());
  TFile* f_bgBox_3  = TFile::Open((prefix+bgBox_3).Data());

  TFile* f_bgQCD_1  = TFile::Open((prefix+bgQCD_1).Data());
  TFile* f_bgQCD_2  = TFile::Open((prefix+bgQCD_2).Data());
  assert(f_bgQCD_2);
  TFile* f_bgPJ_1   = TFile::Open((prefix+bgPJ_1).Data());
  TFile* f_bgPJ_2   = TFile::Open((prefix+bgPJ_2).Data());

  printf("good0\n");


  //TFile* f_data   = TFile::Open("/home/mingyang/cms/hist/hgg-v0/merged/hgg-v0_r11-pho-j16-v1_noskim.root");
  //TFile* f_data   = TFile::Open("/scratch/bendavid/cms/hist/hgg-v0/merged/hgg-v0_f11--h125gg-gf-v14b-pu_noskim.root");
  //TFile* f_data   = TFile::Open("/home/mingyang/cms/hist/2012_ICHEP_PreApproval/merged/2012_ICHEP_PreApproval_r12_noskim.root");
  TFile* f_data   = TFile::Open("/scratch/bendavid/cms/hist/hggMvaTrain38fb/MergedPhotons2012.root");

  TTree* t_signal_1 = (TTree*) f_signal_1->Get((baseDir+TString("hPhotonTree")).Data());
  TTree* t_signal_2 = (TTree*) f_signal_2->Get((baseDir+TString("hPhotonTree")).Data());
  TTree* t_signal_3 = (TTree*) f_signal_3->Get((baseDir+TString("hPhotonTree")).Data());
  TTree* t_signal_4 = (TTree*) f_signal_4->Get((baseDir+TString("hPhotonTree")).Data());

  TTree* t_bgBorn   = (TTree*) f_bgBorn ->Get((baseDir+TString("hPhotonTree")).Data());
  TTree* t_bgBox_1  = (TTree*) f_bgBox_1->Get((baseDir+TString("hPhotonTree")).Data());
  TTree* t_bgBox_2  = (TTree*) f_bgBox_2->Get((baseDir+TString("hPhotonTree")).Data());
  TTree* t_bgBox_3  = (TTree*) f_bgBox_3->Get((baseDir+TString("hPhotonTree")).Data());

  TTree* t_bgQCD_1 = (TTree*) f_bgQCD_1->Get((baseDir+TString("hPhotonTree")).Data());
  TTree* t_bgQCD_2  = (TTree*) f_bgQCD_2->Get((baseDir+TString("hPhotonTree")).Data());
  TTree* t_bgPJ_1   = (TTree*) f_bgPJ_1 ->Get((baseDir+TString("hPhotonTree")).Data());
  TTree* t_bgPJ_2   = (TTree*) f_bgPJ_2 ->Get((baseDir+TString("hPhotonTree")).Data());

  TTree* t_data   = (TTree*) f_data ->Get((baseDir+TString("hPhotonTree")).Data());


  printf("good1\n");

  //-------------------------------------------------------
  // set up PU reweighting histograms.
  // the target is *hpuestnorm 
  //-------------------------------------------------------
  // set up PU reweighting histograms.
  // the target is *hpuestnorm 
  puweights[0] = getpuweights(f_signal_1,hpuest);

  printf("good2\n");

  puweights[1] = getpuweights(f_signal_2,hpuest);
  puweights[2] = getpuweights(f_signal_3,hpuest);
  puweights[3] = getpuweights(f_signal_4,hpuest);
  printf("good3\n");
  puweights[4] = getpuweights(f_bgBorn,hpuest);
  puweights[5] = getpuweights(f_bgBox_1,hpuest);
  puweights[6] = getpuweights(f_bgBox_2,hpuest);
  puweights[7] = getpuweights(f_bgBox_3,hpuest);
  printf("good4\n");
  puweights[8] = getpuweights(f_bgQCD_1,hpuest);
  printf("good5\n");
  puweights[9] = getpuweights(f_bgQCD_2,hpuest);
  printf("good6\n");
  puweights[10] = getpuweights(f_bgPJ_1,hpuest);
  puweights[11] = getpuweights(f_bgPJ_2,hpuest);

  printf("good3\n");

  assert(hpuest);
  assert(f_bgPJ_2);
  
  // compute the weights for the Signals
  TH1D* hNE_signal_1 = (TH1D*) f_signal_1->FindObjectAny("hDAllEvents");
  TH1D* hNE_signal_2 = (TH1D*) f_signal_2->FindObjectAny("hDAllEvents");
  TH1D* hNE_signal_3 = (TH1D*) f_signal_3->FindObjectAny("hDAllEvents");
  TH1D* hNE_signal_4 = (TH1D*) f_signal_4->FindObjectAny("hDAllEvents");

  assert(hNE_signal_1);

  //double theLumi=5089;
  //double theLumi=4800;
  
  double theLumi=5000;

  /* //lumi weights for mh=125 GeV, normalize to theLumi
  double w_signal_1 = theLumi*15.31*2.29/1000./hNE_signal_1->GetBinContent(1);
  double w_signal_2 = theLumi*1.211*2.29/1000./hNE_signal_2->GetBinContent(1);
  double w_signal_3 = theLumi*(0.5729+0.3158)*2.29/1000./hNE_signal_3->GetBinContent(1);
  double w_signal_4 = theLumi*0.08634*2.29/1000./hNE_signal_4->GetBinContent(1);*/

  //lumi weights for mh=123 GeV, normalize to theLumi
  double w_signal_1 = theLumi*20.15*2.27/1000./hNE_signal_1->GetBinContent(1);
  double w_signal_2 = theLumi*1.588*2.27/1000./hNE_signal_2->GetBinContent(1);
  double w_signal_3 = theLumi*(0.7347+0.4150)*2.27/1000./hNE_signal_3->GetBinContent(1);
  double w_signal_4 = theLumi*0.1366*2.27/1000./hNE_signal_4->GetBinContent(1);

  // FP numbers!
  bool doFP = false;
  if( doFP ) {
    w_signal_1 = theLumi*0.*1.54/100./hNE_signal_1->GetBinContent(1);
    w_signal_2 = theLumi*1.211*1.54/100./hNE_signal_2->GetBinContent(1);
    w_signal_3 = theLumi*(0.5729+0.3150)*1.54/100./hNE_signal_3->GetBinContent(1);
    w_signal_4 = theLumi*0.*1.54/100./hNE_signal_4->GetBinContent(1);
  }
  
//   // weights for mh=120 GeV
//   double w_signal_1 = theLumi*16.63*2.225/1000./hNE_signal_1->GetBinContent(1);
//   double w_signal_2 = theLumi*1.269*2.225/1000./hNE_signal_2->GetBinContent(1);
//   double w_signal_3 = theLumi*(0.6561+0.3598)*2.225/1000./hNE_signal_3->GetBinContent(1);
//   double w_signal_4 = theLumi*0.09756*2.225/1000./hNE_signal_4->GetBinContent(1);
// 
//   // FP numbers!
//   bool doFP = false;
//   if( doFP ) {
//     w_signal_1 = theLumi*0.*2.31/100./hNE_signal_1->GetBinContent(1);
//     w_signal_2 = theLumi*1.269*2.31/100./hNE_signal_2->GetBinContent(1);
//     w_signal_3 = theLumi*(0.6561+0.3598)*2.31/100./hNE_signal_3->GetBinContent(1);
//     w_signal_4 = theLumi*0.*2.31/100./hNE_signal_4->GetBinContent(1);
//   }  



  TH1D* hNE_bgBorn = (TH1D*) f_bgBorn->FindObjectAny("hDAllEvents");
  double w_bgBorn  = theLumi*1.15*81/hNE_bgBorn->GetBinContent(1);

  TH1D* hNE_bgBox_1 = (TH1D*) f_bgBox_1->FindObjectAny("hDAllEvents");
  double w_bgBox_1  = theLumi*1.15*424.8/hNE_bgBox_1->GetBinContent(1);

  TH1D* hNE_bgBox_2 = (TH1D*) f_bgBox_2->FindObjectAny("hDAllEvents");
  double w_bgBox_2  = theLumi*1.15*15.54/hNE_bgBox_2->GetBinContent(1);

  TH1D* hNE_bgBox_3 = (TH1D*) f_bgBox_3->FindObjectAny("hDAllEvents");
  double w_bgBox_3  = theLumi*1.15*0.0011805/hNE_bgBox_3->GetBinContent(1);

  TH1D* hNE_bgQCD_1 = (TH1D*) f_bgQCD_1->FindObjectAny("hDAllEvents");
  double w_bgQCD_1  = theLumi*1*(5.195e+07)*0.000235/hNE_bgQCD_1->GetBinContent(1);  
  
  TH1D* hNE_bgQCD_2 = (TH1D*) f_bgQCD_2->FindObjectAny("hDAllEvents");
  double w_bgQCD_2  = theLumi*1*(2.365e+07)*0.002175/hNE_bgQCD_2->GetBinContent(1);

  TH1D* hNE_bgPJ_1 = (TH1D*) f_bgPJ_1->FindObjectAny("hDAllEvents");
  double w_bgPJ_1  = theLumi*1.3*81930.0*0.001835/hNE_bgPJ_1->GetBinContent(1);

  TH1D* hNE_bgPJ_2 = (TH1D*) f_bgPJ_2->FindObjectAny("hDAllEvents");
  double w_bgPJ_2  = theLumi*1.3*8884.0*0.05387/hNE_bgPJ_2->GetBinContent(1);

  printf("good2\n");
  Float_t _masserr, _mass;
  Float_t _ph1pt, _ph2pt;
  Float_t _ph1eta, _ph2eta;
  Float_t _ph1phi, _ph2phi;
  Char_t _ph1cat, _ph2cat;
  Int_t _numPU;
  Float_t _relRes, _dphi, _pt1m, _pt2m;
  Float_t _pth;
  Float_t _relResWrongVtx, _vtxProb;
  Float_t _mvaval, _mvavalup, _mvavaldown;
  Float_t _mvavalidup, _mvavaliddown;
  Float_t _weight;
  Char_t  _proc;
  Int_t _evCat;
  Float_t _res, _resWrong;
  Float_t _relResUp, _relResDown;
  Float_t _relResWrongVtxUp, _relResWrongVtxDown;
  Float_t _ph1idmva, _ph2idmva;

  
  Float_t _ph1eerr, _ph1esmearing, _ph1e;
  Float_t _ph2eerr, _ph2esmearing, _ph2e;
  UChar_t _ismc;
  Bool_t _ph1isbarrel, _ph2isbarrel;
  UChar_t _ph1ispromptgen, _ph2ispromptgen;
  
  Char_t _vbftag1;
  Char_t _vbftag2;

  Float_t _jet1pt, _jet2pt;
  Float_t _jet1eta, _jet2eta;
  Float_t _jet1phi, _jet2phi;
  Float_t _jet1mass, _jet2mass;
  Float_t _dijetmass, _zeppenfeld, _dphidijetgg;  
  
  Float_t _relRes1, _relRes2;

  //add ming
  Float_t _ph1r9, _ph2r9;
  
  //-----add for sync-----
  UInt_t _run;
  UInt_t _lumi;
  UInt_t _evt;
  Float_t _ptgg;
  Float_t _deltamvtx;
  Float_t _pfmet;
  Float_t _relresraw1, _relresraw2;
  Float_t _etcorecaliso1, _etcorecaliso2;
  Float_t _etcorhcaliso1, _etcorhcaliso2;
  Float_t _etcortrkiso1, _etcortrkiso2;
  Float_t _abstrkisocic1, _abstrkisocic2;
  Float_t _hollowconetrkiso031, _hollowconetrkiso032;
  Float_t _ecaliso031, _ecaliso032;
  Float_t _hcaliso031, _hcaliso032;
  Float_t _hovere1, _hovere2;
  Float_t _sigietaieta1,_sigietaieta2;
  Float_t _vtxz, _genz;
  Float_t _rho;
  Float_t _pucorrhcalecal1, _pucorrhcalecal2;

  TMVA::Reader* tmva = new TMVA::Reader();
  printf("good3\n");
  tmva->AddVariable("masserrsmeared/mass",    &_relRes);
  tmva->AddVariable("masserrsmearedwrongvtx/mass",    &_relResWrongVtx);
  tmva->AddVariable("vtxprob",    &_vtxProb);  
  tmva->AddVariable("ph1.pt/mass"                , &_pt1m  );
  tmva->AddVariable("ph2.pt/mass"                , &_pt2m  );
  tmva->AddVariable("ph1.eta"                  , &_ph1eta);
  tmva->AddVariable("ph2.eta"                  , &_ph2eta);
  tmva->AddVariable("TMath::Cos(ph1.phi-ph2.phi)"   , &_dphi  );    
  tmva->AddVariable("ph1.idmva"   , &_ph1idmva  );    
  tmva->AddVariable("ph2.idmva"   , &_ph2idmva  );    
  printf("good4\n");
  TFile* outFile = NULL;
  //tmva->BookMVA("BDTG","/home/mingyang/cms/root/weights/HggBambu_SMDipho_Jan16_BDTG.weights.xml");
  tmva->BookMVA("BDTG","/home/mingyang/cms/root/RootFiles/hggmva_2012_jan16/diphotonMVA/weights/HggBambu_SMDipho_Jan16_BDTG.weights.xml");
  printf("good5\n");
  //outFile= new TFile("/home/mingyang/cms/root/RootFiles/hggmva_2012_jan16/Tmva_AppOutput_SMDipho_2012Jan16_sync.root","RECREATE");
  //outFile= new TFile("/home/mingyang/cms/root/RootFiles/hggmva_2012_jan16/Tmva_AppOutput_SMDipho_2012Jan16_4.8fb-1.root","RECREATE");
  //outFile= new TFile("/home/mingyang/cms/root/RootFiles/hggmva_2012_jan16/Tmva_AppOutput_SMDipho_2012Jan16_NewVBF_120.root","RECREATE");
  //outFile= new TFile("/home/mingyang/cms/root/RootFiles/hggmva_2012_jan16/diphotonMVA/Tmva_AppOutput_SMDipho_2012Jan16_CorrScaleFForMakeDatacard.root","RECREATE");
  //outFile= new TFile("/home/mingyang/cms/root/RootFiles/hggmva_2012_jan16/diphotonMVA/Tmva_AppOutput_SMDipho_2012Jan16_CorrScaleF.root","RECREATE");
  //TNtuple* outTuple = new TNtuple("MITMVAtuple","","mass:res:reswrong:vtxprob:mva:weight:proc:cat1:cat2:pt1:pt2:eta1:eta2:phi1:phi2:numPU:evcat");
  //outFile= new TFile("/home/mingyang/cms/root/RootFiles/2012_ICHEP/DiphotonMVA/oldweight/Tmva_AppOutput_SMDipho_2012ICHEP_oldweight_IdShift.root","RECREATE");
  outFile= new TFile("/home/mingyang/cms/root/RootFiles/2012_ICHEP/DiphotonMVA/oldweight/Tmva_AppOutput_SMDipho_2012ICHEP_oldweight_IdShift_38_hgg120.root","RECREATE");

  TTree* outTuple = new TTree("MITMVAtuple","");
  outTuple->Branch("mass", &_mass, "mass/F");    
  outTuple->Branch("res", &_relRes, "res/F");    
  outTuple->Branch("reswrong", &_relResWrongVtx, "reswrong/F");    
  outTuple->Branch("vtxprob", &_vtxProb, "vtxprob/F");
  outTuple->Branch("mva", &_mvaval, "mva/F");
  outTuple->Branch("mvaup", &_mvavalup, "mvaup/F");  
  outTuple->Branch("mvadown", &_mvavaldown, "mvadown/F");  
  outTuple->Branch("mvaidup", &_mvavalidup, "mvaidup/F");  
  outTuple->Branch("mvaiddown", &_mvavaliddown, "mvaiddown/F");  
  outTuple->Branch("weight", &_weight, "weight/F");
  outTuple->Branch("proc", &_proc, "proc/B");
  outTuple->Branch("cat1", &_ph1cat, "proc/B");
  outTuple->Branch("cat2", &_ph2cat, "proc/B");
  outTuple->Branch("pt1", &_ph1pt, "pt1/F");
  outTuple->Branch("pt2", &_ph2pt, "pt2/F");
  outTuple->Branch("eta1", &_ph1eta, "eta1/F");
  outTuple->Branch("eta2", &_ph2eta, "eta2/F");
  outTuple->Branch("phi1", &_ph1phi, "phi1/F");
  outTuple->Branch("phi2", &_ph2phi, "phi2/F");
  outTuple->Branch("idmva1", &_ph1idmva, "idmva1/F");
  outTuple->Branch("idmva2", &_ph2idmva, "idmva2/F");  
  outTuple->Branch("numPU", &_numPU, "numPU/I");
  outTuple->Branch("evcat", &_evCat, "evcat/I");
  //outTuple->Branch("vbftag", &_vbftag, "vbftag/B");
  outTuple->Branch("vbftag1", &_vbftag1, "vbftag1/B");
  outTuple->Branch("vbftag2", &_vbftag2, "vbftag2/B");
  outTuple->Branch("jet1pt", &_jet1pt, "jet1pt/F");//ming
  outTuple->Branch("jet2pt", &_jet2pt, "jet2pt/F");//ming
  outTuple->Branch("jet1eta", &_jet1eta, "jet1eta/F");//ming
  outTuple->Branch("jet2eta", &_jet2eta, "jet2eta/F");//ming
  outTuple->Branch("jet1phi", &_jet1phi, "jet1phi/F");
  outTuple->Branch("jet2phi", &_jet2phi, "jet2phi/F");
  outTuple->Branch("jet1mass", &_jet1mass, "jet1mass/F");
  outTuple->Branch("jet2mass", &_jet2mass, "jet2mass/F");
  outTuple->Branch("dijetmass", &_dijetmass, "dijetmass/F");
  outTuple->Branch("zeppenfeld", &_zeppenfeld, "zeppenfeld/F");
  outTuple->Branch("dphidijetgg", &_dphidijetgg, "dphidijetgg/F");  
  outTuple->Branch("relres1", &_relRes1, "relres1/F");  
  outTuple->Branch("relres2", &_relRes2, "relres2/F");  
  outTuple->Branch("resup", &_relResUp, "resup/F");  
  outTuple->Branch("resdown", &_relResDown, "resdown/F");  
  //ming add
  outTuple->Branch("r91", &_ph1r9, "r91/F");
  outTuple->Branch("r92", &_ph2r9, "r92/F");

  //-------add for sync---------
  outTuple->Branch("run", &_run, "run/i");  
  outTuple->Branch("lumi", &_lumi, "lumi/i");  
  outTuple->Branch("evt", &_evt, "evt/i");  
  outTuple->Branch("ptgg", &_ptgg, "ptgg/F");  
  outTuple->Branch("pfmet", &_pfmet, "pfmet/F");  
  outTuple->Branch("e1",&_ph1e,"e1/F");
  outTuple->Branch("e2",&_ph2e,"e2/F");
  outTuple->Branch("relresraw1",&_relresraw1,"relresraw1/F");
  outTuple->Branch("relresraw2",&_relresraw2,"relresraw2/F");  
  outTuple->Branch("etcorecaliso1",&_etcorecaliso1,"etcorecaliso1/F");
  outTuple->Branch("etcorecaliso2",&_etcorecaliso2,"etcorecaliso2/F");  
  outTuple->Branch("etcorhcaliso1",&_etcorhcaliso1,"etcorhcaliso1/F");
  outTuple->Branch("etcorhcaliso2",&_etcorhcaliso2,"etcorhcaliso2/F");    
  outTuple->Branch("etcortrkiso1",&_etcortrkiso1,"etcortrkiso1/F");
  outTuple->Branch("etcortrkiso2",&_etcortrkiso2,"etcortrkiso2/F");    
  outTuple->Branch("abstrkisocic1",&_abstrkisocic1,"abstrkisocic1/F");
  outTuple->Branch("abstrkisocic2",&_abstrkisocic2,"abstrkisocic2/F");    
  outTuple->Branch("hollowconetrkiso031",&_hollowconetrkiso031,"hollowconetrkiso031/F");
  outTuple->Branch("hollowconetrkiso032",&_hollowconetrkiso032,"hollowconetrkiso032/F");  
  outTuple->Branch("sigietaieta1",&_sigietaieta1,"sigietaieta1/F");//ming add
  outTuple->Branch("sigietaieta2",&_sigietaieta2,"sigietaieta2/F");//ming add
  outTuple->Branch("hovere1",&_hovere1,"hovere1/F");//ming add
  outTuple->Branch("hovere2",&_hovere2,"hovere2/F");//ming add
  outTuple->Branch("rho",&_rho,"rho/F");//ming add 
  outTuple->Branch("pucorrhcalecal1",& _pucorrhcalecal1,"pucorrhcalecal1/F");//ming add
  outTuple->Branch("pucorrhcalecal2",&_pucorrhcalecal2,"pucorrhcalecal2/F");//ming add
 
  outTuple->SetAutoFlush(-100000000);

  printf("good6\n");
  TTree** allTrees = new TTree*[13];
  allTrees[0] = t_signal_1;
  allTrees[1] = t_signal_2;
  allTrees[2] = t_signal_3;
  allTrees[3] = t_signal_4;
  allTrees[4] = t_bgBorn;
  allTrees[5] = t_bgBox_1;
  allTrees[6] = t_bgBox_2;
  allTrees[7] = t_bgBox_3;
  allTrees[8] = t_bgQCD_1;
  allTrees[9] = t_bgQCD_2;
  allTrees[10] = t_bgPJ_1;
  allTrees[11] = t_bgPJ_2;
  allTrees[12] = t_data;
  

  double* weights = new double[12];
  weights[0] = w_signal_1;
  weights[1] = w_signal_2;
  weights[2] = w_signal_3;
  weights[3] = w_signal_4;
  weights[4] = w_bgBorn;
  weights[5] = w_bgBox_1;
  weights[6] = w_bgBox_2;
  weights[7] = w_bgBox_3;
  weights[8] = w_bgQCD_1;
  weights[9] = w_bgQCD_2;
  weights[10] = w_bgPJ_1;
  weights[11] = w_bgPJ_2;

  printf("good7\n");

  for(int jj=0; jj<13; ++jj) {
    TTree* theTree = allTrees[jj];    
    printf("good71\n");
    theTree->SetBranchAddress("ismc"   ,     &_ismc);
    theTree->SetBranchAddress("ph1.isbarrel"   ,     &_ph1isbarrel);
    theTree->SetBranchAddress("ph2.isbarrel"   ,     &_ph2isbarrel);
    theTree->SetBranchAddress("masserrsmeared"   ,     &_res);
    theTree->SetBranchAddress("masserrsmearedwrongvtx"   ,     &_resWrong);
    theTree->SetBranchAddress("vtxprob"   ,     &_vtxProb);
    theTree->SetBranchAddress("ph1.eerr"   ,     &_ph1eerr);
    theTree->SetBranchAddress("ph2.eerr"   ,     &_ph2eerr);
    theTree->SetBranchAddress("ph1.esmearing"   ,     &_ph1esmearing);
    theTree->SetBranchAddress("ph2.esmearing"   ,     &_ph2esmearing);
    theTree->SetBranchAddress("ph1.e"   ,     &_ph1e);
    theTree->SetBranchAddress("ph2.e"   ,     &_ph2e);    

    theTree->SetBranchAddress("mass"   ,     &_mass   );
    theTree->SetBranchAddress("ph1.pt" ,     &_ph1pt  );
    theTree->SetBranchAddress("ph2.pt" ,     &_ph2pt  );
    theTree->SetBranchAddress("ph1.eta" ,  &_ph1eta  );
    theTree->SetBranchAddress("ph2.eta" ,  &_ph2eta  );
    theTree->SetBranchAddress("ph1.phi" ,  &_ph1phi  );
    theTree->SetBranchAddress("ph2.phi" ,  &_ph2phi  );
    theTree->SetBranchAddress("ph1.phcat" ,  &_ph1cat  );
    theTree->SetBranchAddress("ph2.phcat" ,  &_ph2cat  );  
    theTree->SetBranchAddress("ph1.idmva" ,  &_ph1idmva  );
    theTree->SetBranchAddress("ph2.idmva" ,  &_ph2idmva  );    
    theTree->SetBranchAddress("ph1.ispromptgen" ,  &_ph1ispromptgen  );
    theTree->SetBranchAddress("ph2.ispromptgen" ,  &_ph2ispromptgen  );    
    theTree->SetBranchAddress("numPU"      , &_numPU  );
    theTree->SetBranchAddress("genHiggspt" , &_pth    );
    
    theTree->SetBranchAddress("jet1pt" , &_jet1pt    );
    theTree->SetBranchAddress("jet2pt" , &_jet2pt    );
    theTree->SetBranchAddress("jet1eta" , &_jet1eta    );
    theTree->SetBranchAddress("jet2eta" , &_jet2eta    );  
    theTree->SetBranchAddress("jet1phi" , &_jet1phi);//ming add;  
    theTree->SetBranchAddress("jet2phi" , &_jet2phi);//ming add;  
    theTree->SetBranchAddress("jet1mass" , &_jet1mass);//ming add;  
    theTree->SetBranchAddress("jet2mass" , &_jet2mass);//ming add;  
    theTree->SetBranchAddress("dijetmass" , &_dijetmass    );
    theTree->SetBranchAddress("zeppenfeld" , &_zeppenfeld    );
    theTree->SetBranchAddress("dphidijetgg" , &_dphidijetgg    );    

    theTree->SetBranchAddress("run" , &_run    );  
    theTree->SetBranchAddress("lumi" , &_lumi    );    
    theTree->SetBranchAddress("evt" , &_evt    ); 
    theTree->SetBranchAddress("ptgg" , &_ptgg    );    
    theTree->SetBranchAddress("deltamvtx" , &_deltamvtx    );    
    

    theTree->SetBranchAddress("ph1.trkisohollowdr03" , &_hollowconetrkiso031    );    
    theTree->SetBranchAddress("ph2.trkisohollowdr03" , &_hollowconetrkiso032    );    
    theTree->SetBranchAddress("ph1.ecalisodr03" , &_ecaliso031    );    
    theTree->SetBranchAddress("ph2.ecalisodr03" , &_ecaliso032    );        
    theTree->SetBranchAddress("ph1.hcalisodr03" , &_hcaliso031    );    
    theTree->SetBranchAddress("ph2.hcalisodr03" , &_hcaliso032    );        
    theTree->SetBranchAddress("ph1.trackiso1" , &_abstrkisocic1   );    
    theTree->SetBranchAddress("ph2.trackiso1" , &_abstrkisocic2   );            
    theTree->SetBranchAddress("ph1.sigietaieta" , &_sigietaieta1   );//ming add
    theTree->SetBranchAddress("ph2.sigietaieta" , &_sigietaieta2  );//ming add   
    theTree->SetBranchAddress("ph1.hovere" , &_hovere1  );
    theTree->SetBranchAddress("ph2.hovere" , &_hovere2  );     
    theTree->SetBranchAddress("rho" , &_rho  );   
    theTree->SetBranchAddress("vtxZ" , &_vtxz);            
    theTree->SetBranchAddress("genHiggsZ" , &_genz);      
    theTree->SetBranchAddress("pfmet" , &_pfmet);    

    //ming add
    theTree->SetBranchAddress("ph1.r9" ,     &_ph1r9  );
    theTree->SetBranchAddress("ph2.r9" ,     &_ph2r9  );       
     
    TTreeFormula basecut("basecut","ph1.pt>(mass/3.0) && ph2.pt>(mass/4.0) && ph1.pt>(100.0/3.0) && ph2.pt>(100.0/4.0) && mass > 100. && ((((ph1.idmva-0.1)>-0.3) && ph1.isbarrel) || (((ph1.idmva-0.13)>-0.3) && !ph1.isbarrel)) && ((((ph2.idmva-0.1)>-0.3) && ph2.isbarrel) || (((ph2.idmva-0.13)>-0.3) && !ph2.isbarrel))",theTree);
    //TTreeFormula basecut("basecut","mass > 100.",theTree);
    //TTreeFormula vbfcut("vbfcut","jet1pt>30.0 && jet2pt>20.0 && abs(jet1eta-jet2eta)>3.5 && dijetmass>350.0 && zeppenfeld<2.5 && abs(dphidijetgg)>2.6",theTree);//ming:zeppenfeld?
    //TTreeFormula vbfcut("vbfcut","jet1pt>30.0 && jet2pt>20.0 && abs(jet1eta-jet2eta)>3.5 && dijetmass>350.0 && zeppenfeld<2.5 && abs(dphidijetgg)>2.6 && ph1.pt>(mass*55.0/120.0)",theTree);//ming:zeppenfeld?

    // TTreeFormula vbfcut1("vbfcut1","(ph1.pt>(mass/2) && ph2.pt>(25) && mass > 100. && ph1.idmva>-0.3 && ph2.idmva>-0.3) && (jet1pt>30.0 && jet2pt>20.0 && abs(jet1eta-jet2eta)>3 && dijetmass>250.0 && zeppenfeld<2.5 && abs(dphidijetgg)>2.6) && (jet2pt>30.0 && dijetmass>500.0)",theTree);
    //TTreeFormula vbfcut2("vbfcut2","(ph1.pt>(mass/2) && ph2.pt>(25) && mass > 100. && ph1.idmva>-0.3 && ph2.idmva>-0.3) && (jet1pt>30.0 && jet2pt>20.0 && abs(jet1eta-jet2eta)>3 && dijetmass>250.0 && zeppenfeld<2.5 && abs(dphidijetgg)>2.6) && !(jet2pt>30.0 && dijetmass>500.0)",theTree);

    TTreeFormula vbfcut1("vbfcut1","(ph1.pt>(mass/2) && ph2.pt>(mass/4) && mass > 100. && ((((ph1.idmva-0.1)>-0.3) && ph1.isbarrel) || (((ph1.idmva-0.13)>-0.3) && !ph1.isbarrel)) && ((((ph2.idmva-0.1)>-0.3) && ph2.isbarrel) || (((ph2.idmva-0.13)>-0.3) && !ph2.isbarrel))) && (jet1pt>30.0 && jet2pt>20.0 && abs(jet1eta-jet2eta)>3 && dijetmass>250.0 && zeppenfeld<2.5 && abs(dphidijetgg)>2.6) && (jet2pt>30.0 && dijetmass>500.0)",theTree);
    TTreeFormula vbfcut2("vbfcut2","(ph1.pt>(mass/2) && ph2.pt>(mass/4) && mass > 100. && ((((ph1.idmva-0.1)>-0.3) && ph1.isbarrel) || (((ph1.idmva-0.13)>-0.3) && !ph1.isbarrel)) && ((((ph2.idmva-0.1)>-0.3) && ph2.isbarrel) || (((ph2.idmva-0.13)>-0.3) && !ph2.isbarrel))) && (jet1pt>30.0 && jet2pt>20.0 && abs(jet1eta-jet2eta)>3 && dijetmass>250.0 && zeppenfeld<2.5 && abs(dphidijetgg)>2.6) && !(jet2pt>30.0 && dijetmass>500.0)",theTree);
    
    for(int i=1; i<=theTree->GetEntries(); ++i) {
      theTree->GetEntry(i);
      //       if( _ph1pt/_mass < 40./120. || _ph2pt/_mass < 30./120. ) continue;
      //       if( _mass < 100. ) continue;

      bool pass = basecut.EvalInstance();
      if (!pass) continue;
     
      _pt1m  =_ph1pt/_mass;
      _pt2m  =_ph2pt/_mass;
      _dphi  = TMath::Cos(_ph1phi - _ph2phi);
      
      _relRes = _res/_mass;
      _relResWrongVtx = _resWrong/_mass;
      
      double dmvtx = TMath::Sqrt(_relResWrongVtx*_relResWrongVtx - _relRes*_relRes);

      if(_ph1isbarrel){
	_ph1idmva=_ph1idmva-0.1;
      }

      if(!_ph1isbarrel){
	_ph1idmva=_ph1idmva-0.13;
      }

      if(_ph2isbarrel){
	_ph2idmva=_ph2idmva-0.1;
      }
      
      if(!_ph2isbarrel){
	_ph2idmva=_ph2idmva-0.13;
      }
      
      double id1 = _ph1idmva;
      double id2 = _ph2idmva;
      
      _ph1idmva = id1 + 0.025;
      _ph2idmva = id2 + 0.025;
      _mvavalidup = tmva->EvaluateMVA("BDTG");

      _ph1idmva = id1 - 0.025;
      _ph2idmva = id2 - 0.025;
      _mvavaliddown = tmva->EvaluateMVA("BDTG");   
      
      _ph1idmva = id1;
      _ph2idmva = id2;      
      
      _relRes1 = (1.0+_ismc*_ph1isbarrel*0.07+_ismc*(!_ph1isbarrel)*0.045)*_ph1eerr/_ph1e;
      _relRes2 = (1.0+_ismc*_ph2isbarrel*0.07+_ismc*(!_ph2isbarrel)*0.045)*_ph2eerr/_ph2e;
      
      _relResUp  = (pow(1.0+_ismc*_ph1isbarrel*0.07+_ismc*(!_ph1isbarrel)*0.045,2)*1.1*1.1*_ph1eerr*_ph1eerr+_ph1esmearing*_ph1esmearing)/_ph1e/_ph1e;//ming: 1.1 is estimated by eye from the comparison between data ans MC
      _relResUp += (pow(1.0+_ismc*_ph2isbarrel*0.07+_ismc*(!_ph2isbarrel)*0.045,2)*1.1*1.1*_ph2eerr*_ph2eerr+_ph2esmearing*_ph2esmearing)/_ph2e/_ph2e;
//       _relResUp  = _relRes1*_relRes1*1.05*1.05;
//       _relResUp += _relRes2*_relRes2*1.05*1.05;
      _relResUp  = 0.5*sqrt(_relResUp);  
      _relResWrongVtxUp = TMath::Sqrt(_relResUp*_relResUp + dmvtx*dmvtx);
      _relRes = _relResUp;
      _relResWrongVtx = _relResWrongVtxUp;
      
      _mvavalup = tmva->EvaluateMVA("BDTG");
      
      _relResDown  = (pow(1.0+_ismc*_ph1isbarrel*0.07+_ismc*(!_ph1isbarrel)*0.045,2)*0.9*0.9*_ph1eerr*_ph1eerr+_ph1esmearing*_ph1esmearing)/_ph1e/_ph1e;
      _relResDown += (pow(1.0+_ismc*_ph2isbarrel*0.07+_ismc*(!_ph2isbarrel)*0.045,2)*0.9*0.9*_ph2eerr*_ph2eerr+_ph2esmearing*_ph2esmearing)/_ph2e/_ph2e;
//       _relResDown  = _relRes1*_relRes1*0.95*0.95;
//       _relResDown += _relRes2*_relRes2*0.95*0.95;
      _relResDown  = 0.5*sqrt(_relResDown);  
      _relResWrongVtxDown = TMath::Sqrt(_relResDown*_relResDown + dmvtx*dmvtx);
      _relRes = _relResDown;
      _relResWrongVtx = _relResWrongVtxDown;
      
      _mvavaldown = tmva->EvaluateMVA("BDTG");
      
      _relRes = _res/_mass;
      _relResWrongVtx = _resWrong/_mass;
            
      _mvaval = tmva->EvaluateMVA("BDTG");
    
      _evCat=eventCat(_ph1cat,_ph2cat);
      _proc = jj;
      
      //_vbftag = vbfcut.EvalInstance();
      _vbftag1 = vbfcut1.EvalInstance();
      _vbftag2 = vbfcut2.EvalInstance();
     
      if( (_proc==4) && ((_evt>=270000 && _evt<=390000) || (_evt>=480000 && _evt<=510000) || (_evt>=1170000 && _evt<=1200000) || (_evt>=1230000 && _evt<=1260000)  || (_evt>=1890000 && _evt<=1920000)  || (_evt>=2670000 && _evt<=2730000) || (_evt>=3090000 && _evt<=3210000) || (_evt>=3270000 && _evt<=3300000) || (_evt>=3510000 && _evt<=3600000) || (_evt>=4890000 && _evt<=5010000) || (_evt>=5340000 && _evt<=5370000) || (_evt>=5430000 && _evt<=5500000))) continue;
    
      //-----add for syn----
      _etcorecaliso1 = _ecaliso031 - 0.012*_ph1pt;
      _etcorecaliso2 = _ecaliso032 - 0.012*_ph2pt;
      _etcorhcaliso1 = _hcaliso031 - 0.005*_ph1pt;
      _etcorhcaliso2 = _hcaliso032 - 0.005*_ph2pt;
      _etcortrkiso1 = _hollowconetrkiso031 - 0.002*_ph1pt;
      _etcortrkiso2 = _hollowconetrkiso032 - 0.002*_ph2pt;
      _pucorrhcalecal1= _hcaliso031+ _ecaliso031-0.17*_rho;
      _pucorrhcalecal2= _hcaliso032+_ecaliso032-0.17*_rho;
      _relresraw1 = _ph1eerr/_ph1e;
      _relresraw2 = _ph2eerr/_ph2e;  
     
      
      // compute the weight
      /*if (jj<11) {
        _weight = weights[jj] * puweight(_numPU,jj) * ptweight(_pth,jj) * effweight(_ph1cat) *effweight(_ph2cat) * trigeffweight(_ph1cat,_ph2cat,_ph1pt, _ph2pt);
        if (!_ph1ispromptgen && !_ph2ispromptgen) _weight*=1.0/1.3;
      }
      else _weight = 1.0;*/

      //ming change
      /*if (jj<11) {
        //_weight = weights[jj] * puweight(_numPU,jj) * ptweight(_pth,jj) * effweight(_ph1cat,_ph1r9) *effweight(_ph2cat,_ph2r9) * trigeffweight(_ph1cat,_ph2cat,_ph1pt, _ph2pt);
	_weight = weights[jj] * puweight(_numPU,jj) * effweight(_ph1cat,_ph1r9) *effweight(_ph2cat,_ph2r9) * trigeffweight(_ph1cat,_ph2cat,_ph1pt, _ph2pt);
        if (!_ph1ispromptgen && !_ph2ispromptgen) _weight*=1.0/1.3;
	}*/

      if (jj<12) {
	_weight = weights[jj] * puweight(_numPU,jj);
      }
      else _weight = 1.0;
      
      outTuple->Fill();
      
    }
    printf("good8\n");
  }
    
  outTuple->Write();
  outFile->Close();

  return;
}
