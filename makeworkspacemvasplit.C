//----------------------------------------------------------------------------------------------------------------
//ming: this macro produce three RooWorkspace which will be used as the input for background and signal modeling
//RooWorkspace wdata contains: 
//1) RooDataSet for data for each category passing preselection, within kinematic acceptance and idmva1>-0.3 && idmva2>-0.3 
//2) lumi
//RooWorkspace wmclow and wmclow contains: 
//1) RooDataSet for mc at each mass point passing preselection, within kinematic acceptance and idmva1>-0.3 && idmva2>-0.3 cateorized in mva category, production mechanism and right or wrong vertex
//1) weighted by lumi,pu,higgs pt,preselection scale factor and electron veto scale factor, and trigger efficency
//2) cross section times BR for each production mechanism at each mass
//3) lumi

#include <sstream>
#include <iostream>
#include "TMath.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TCut.h"
#include "TTree.h"
#include "TNtuple.h"
#include "RooRealVar.h"
#include "RooArgSet.h"
#include "RooArgList.h"
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
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooWorkspace.h"
#include "RooConstVar.h"
#include "TEfficiency.h"
#include "TMVA/Reader.h"

using namespace RooFit;

Float_t *g_mvavars;
TMVA::Reader *g_reader;

Double_t getmvaval(Float_t var0, Float_t var1, Float_t var2, Float_t var3, Float_t var4, Float_t var5, Float_t var6, Float_t var7, Float_t var8, Float_t var9) {
  g_mvavars[0] = var0;
  g_mvavars[1] = var1;
  g_mvavars[2] = var2;
  g_mvavars[3] = var3;
  g_mvavars[4] = var4;
  g_mvavars[5] = var5;
  g_mvavars[6] = var6;
  g_mvavars[7] = var7;
  g_mvavars[8] = var8;
  g_mvavars[9] = var9;
  //return g_reader->EvaluateMVA("Ystar");
  return g_reader->EvaluateMVA("BDTG");
}

//pileup reweighting
TH1D *puweights[4];//ming:the reason for various puweights[i] is that sometime we might want different puweights for different MC sample in case different MC has different pu distribution
float puweight(float npu, int wset=0) {
  if (npu<0) return 1.0;
  return puweights[wset]->GetBinContent(puweights[wset]->FindFixBin(npu));
}

void setpuweights(TFile *file, TH1D *target, int wset=0) {//ming:set puweights for each process; might not be necessary when the pu distributions for different mc are the same
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
  
  puweights[wset] = new TH1D((*htargettmp)/(*hpumc));
    
  delete htargettmp;
  
}

//pt-reweighing (but only for gf samples, use extra bool flag for now);reweight the higgs pt distribution for gf
//ming: the gf pt distribution from powheg is wrong and needs to be corrected with weights; the pt distribution for other higgs mechanism is fine and needs no correction
TH1D *ptweights = 0;
float ptweight(float genhpt, Int_t procid) {
  if (procid>0) return 1.0;
  if (genhpt<0) return 1.0;
  
 // if (genhpt>200.0) {
    //double wgt = ptweights->GetBinContent(ptweights->FindFixBin(genhpt));
    //printf("genhpt = %5f, wgt = %5f\n",genhpt,wgt);
 // }
  
  return ptweights->GetBinContent(ptweights->FindFixBin(genhpt));
}

//efficiency scale factors (barrel/endcap for now, but can be extended to full kinematic binning)
float effweight(int cat, float r9) {
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
  
  //mva id  - ele veto only
  // if (cat==1) return 1.002;
  //   else if (cat==2) return 1.011;
  //   else if (cat==3) return 1.000;
  //   else if (cat==4) return 0.996;
  //   else return 1.0;  
  
  double scale=1.0;
  
  bool iseb = (cat==1 || cat==2);
  bool isr9 = r9>0.9;
  
  //mva id  - ele veto only
  //if (cat==1) scale*=1.002;
  //else if (cat==2) scale*=1.011;
  //else if (cat==3) scale*=1.000;
  //else if (cat==4) scale*=0.996;
  
  //feb2 presel scale factors
  //if (iseb&&isr9) scale*=0.99206;
  //else if (iseb&&!isr9) scale*=0.99741;
  //else if (!iseb&&isr9) scale*=1.01633;
  //else if (!iseb && !isr9) scale*=1.02007;
  
  //ele veto (from nov30) in CIC category (number from Doug)
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
  
  return scale;
}

//nov15 freeze
//ming:trigeff from Vladmir 
//ming:these trigeff are for four eta, r9 categories considering full trigger sets
float trigeffweight(int cat1, int cat2, double pt1, double pt2) {
  Bool_t isb1 = cat1==1 || cat1==2;
  Bool_t isb2 = cat2==1 || cat2==2;
  
  Bool_t isr91 = cat1==1 || cat1==3;
  Bool_t isr92 = cat2==1 || cat2==3;
  
  Bool_t isb = isb1 && isb2;
  Bool_t isr9 = isr91 && isr92;
  
  if (isb && isr9) return 1.0;//both barrel and both high r9
  else if (isb && !isr9) return 0.993;//both barrel not both high r9
  else if (!isb && isr9) return 1.0;//not both barrel but both high r9
  else if (!isb && !isr9) return 0.988;//not both barrel not both high r9
  else return 1.0;
}

//oct18 mva stuff
double trigEff[36];
double trigEffErr[36];
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

// float vtxweight(float dz) {
//   if (
// }

//append data from file to RooDataSet adding a column which has the weight by given cross section
//(Note that the weight is not enabled at this stage, only made available for subsequent use)
void append(RooDataSet &data, TFile *infile, TCut sel, Double_t xsec, Int_t procid, TString modname) {//ming:what's the purpose of this?

    TDirectory *hdirfwk = (TDirectory*) infile->FindObjectAny("AnaFwkMod");
    const TH1D *hDAllEvents = (TH1D*)hdirfwk->Get("hDAllEvents");
  
  
    TDirectory *hdir = (TDirectory*) infile->FindObjectAny(modname);
    TTree *hdata = (TTree*)hdir->Get("hPhotonTree");  
    

    RooRealVar xsecweight("xsecweight","xsecweight",xsec/(double)hDAllEvents->GetEntries());
    RooRealVar procidx("procidx","procidx",procid);
    
    
    RooArgSet varlistsmall = *data.get();
    varlistsmall.remove(xsecweight,kFALSE,kTRUE);//ming:?
    varlistsmall.remove(procidx,kFALSE,kTRUE);//ming:?
    
    
    RooDataSet newdata("newdata","newdata",varlistsmall,RooFit::Import(*hdata),RooFit::Cut(sel)); 
    newdata.addColumn(xsecweight);
    newdata.addColumn(procidx);
   
    
    data.append(newdata);//ming:append newdata to data  
}


void makeworkspacemvasplit() {
  //---set plot style---
  gROOT->Macro("MitStyle.C");
  gStyle->SetErrorX(0);
  gStyle->SetOptStat(1110); 
  gROOT->ForceStyle();  

  //---change directory---
  //gSystem->cd("/home/mingyang/cms/root/RootFiles/hggmva_2012_jan16/RooWorkSpace");
  gSystem->cd("./test/");

  //---set overall parameters---
  bool doff = false;//ming:choose whether to do fermiphobic or not
  bool dott = true;//ming:choose to do tt or not
  float totallumi = 5089;
  TString modname = "PhotonTreeWriterPresel";

  //---define selection cuts---
  //overall cut
  TCut loosesel = "mass>=100 && mass<=180 && ph1.pt>(mass/3.0) && ph2.pt>(mass/4.0) && ph1.pt>(100.0/3.0) && ph2.pt>(100.0/4.0)";
  TCut presel = "ph1.idmva>-0.3 && ph2.idmva>-0.3";
  TCut vbfcut = "(jet1pt>30.0 && jet2pt>20.0 && abs(jet1eta-jet2eta)>3.5 && dijetmass>350.0 && zeppenfeld<2.5 && abs(dphidijetgg)>2.6) && ph1.pt>(55.0*mass/120.0)";
  
  //category cut: 2012jan16 mva && vbf
  TCut cat0 = loosesel && presel && !vbfcut && "bdt>= 0.89";
  TCut cat1 = loosesel && presel && !vbfcut && "bdt>= 0.74 && bdt<0.89";
  //TCut cat2 = loosesel && presel && !vbfcut && "bdt>= 0.545 && bdt<0.74";
  //TCut cat3 = loosesel && presel && !vbfcut && "bdt>= 0.05 && bdt<0.545"; 
  TCut cat2 = loosesel && presel && !vbfcut && "bdt>= 0.55 && bdt<0.74";
  TCut cat3 = loosesel && presel && !vbfcut && "bdt>= 0.05 && bdt<0.55";   
  TCut cat4 = loosesel && presel &&  vbfcut && "bdt>= 0.05";

  TString CMEnergy = TString("_8TeV");
  std::vector<TString> catnames;  
  catnames.push_back("cat0"+CMEnergy);
  catnames.push_back("cat1"+CMEnergy);
  catnames.push_back("cat2"+CMEnergy);
  catnames.push_back("cat3"+CMEnergy);
  catnames.push_back("cat4"+CMEnergy);

  /*std::vector<TString> catnames;
  std::vector<TCut> catcuts;
  catnames.push_back("cat0");
  catnames.push_back("cat1");
  catnames.push_back("cat2");
  catnames.push_back("cat3");
  catnames.push_back("cat4");*/

  std::vector<TCut> catcuts;
  catcuts.push_back(cat0);
  catcuts.push_back(cat1);
  catcuts.push_back(cat2);
  catcuts.push_back(cat3);
  catcuts.push_back(cat4);
    
  //---load the weights---
  //trigEff
  //ming:the trigEff[i] listed below are computed by Fabian which consider different combinations of eta, pt and r9 but only for a subset of triggers
  //ming:we have never used these number yet but use the one computed by Vladmir for only four eta and r9 categories at this moment
  trigEff[0] = 0.990298;  trigEffErr[0] = 0.00109912;
  trigEff[1] = 0.988554;  trigEffErr[1] = 0.000531549;
  trigEff[2] = 0.977048;  trigEffErr[2] = 0.000820424;
  trigEff[3] = 0.976983;  trigEffErr[3] = 0.000997418;
  trigEff[4] = 0.987159;  trigEffErr[4] = 0.000655609;
  trigEff[5] = 0.990696;  trigEffErr[5] = 0.000517864;
  trigEff[6] = 0.967779;  trigEffErr[6] = 0.00107305;
  trigEff[7] = 0.972583;  trigEffErr[7] = 0.000897231;
  trigEff[8] = 0.986463;  trigEffErr[8] = 0.000383804;
  trigEff[9] = 0.975357;  trigEffErr[9] = 0.000610688;
  trigEff[10] = 0.974622;  trigEffErr[10] = 0.000726944;
  trigEff[11] = 0.985374;  trigEffErr[11] = 0.00047033;
  trigEff[12] = 0.987551;  trigEffErr[12] = 0.000359073;
  trigEff[13] = 0.96588;  trigEffErr[13] = 0.000842646;
  trigEff[14] = 0.969466;  trigEffErr[14] = 0.000663178;
  trigEff[15] = 0.964202;  trigEffErr[15] = 0.000841842;
  trigEff[16] = 0.964274;  trigEffErr[16] = 0.0008088;
  trigEff[17] = 0.973583;  trigEffErr[17] = 0.000583466;
  trigEff[18] = 0.976898;  trigEffErr[18] = 0.000532823;
  trigEff[19] = 0.951719;  trigEffErr[19] = 0.00093695;
  trigEff[20] = 0.958458;  trigEffErr[20] = 0.000848359;
  trigEff[21] = 0.963204;  trigEffErr[21] = 0.00111653;
  trigEff[22] = 0.973444;  trigEffErr[22] = 0.00070631;
  trigEff[23] = 0.974054;  trigEffErr[23] = 0.000792509;
  trigEff[24] = 0.951071;  trigEffErr[24] = 0.000952073;
  trigEff[25] = 0.955716;  trigEffErr[25] = 0.000738544;
  trigEff[26] = 0.983985;  trigEffErr[26] = 0.00070867;
  trigEff[27] = 0.987473;  trigEffErr[27] = 0.000503261;
  trigEff[28] = 0.964332;  trigEffErr[28] = 0.000669253;
  trigEff[29] = 0.969182;  trigEffErr[29] = 0.000861804;
  trigEff[30] = 0.98573;  trigEffErr[30] = 0.000555197;
  trigEff[31] = 0.967892;  trigEffErr[31] = 0.000810588;
  trigEff[32] = 0.967382;  trigEffErr[32] = 0.00090159;
  trigEff[33] = 0.942403;  trigEffErr[33] = 0.00136035;
  trigEff[34] = 0.947787;  trigEffErr[34] = 0.00128003;
  trigEff[35] = 0.948598;  trigEffErr[35] = 0.000990943;
  
  //pileup weights
  //TFile *filepuest = new TFile("/scratch/bendavid/root/puweightsNov13/2011_0100_73500.pileup.root","READ");
  //TFile *filepuest = new TFile("/home/mingyang/cms/puweight/augmented_nov08_rereco.json.68000.pileup.root","READ");
  TFile *filepuest = new TFile("/home/mingyang/cms/puweight/latest_prompt_2012.json.68300.observed.pileup.root","READ");
  TH1D *hpuest = (TH1D*) filepuest->Get("pileup");
  
  //pt-weights
  TFile *fileptweight = new TFile("/scratch/bendavid/root/KFactors_AllScales.root","READ");//ming:?
  
  //---set up reader---
  //get new reader
  g_reader = new TMVA::Reader();
  //get input variables
  g_mvavars = new Float_t[10];
  g_reader->AddVariable("masserrsmeared/mass",    &g_mvavars[0]);
  g_reader->AddVariable("masserrsmearedwrongvtx/mass",    &g_mvavars[1]);
  g_reader->AddVariable("vtxprob",    &g_mvavars[2]);
  g_reader->AddVariable("ph1.pt/mass"                , &g_mvavars[3]);
  g_reader->AddVariable("ph2.pt/mass"                , &g_mvavars[4]);
  g_reader->AddVariable("ph1.eta"                  , &g_mvavars[5]);
  g_reader->AddVariable("ph2.eta"                  , &g_mvavars[6]);
  g_reader->AddVariable("TMath::Cos(ph1.phi-ph2.phi)"   , &g_mvavars[7]);
  g_reader->AddVariable("ph1.idmva"                , &g_mvavars[8]);
  g_reader->AddVariable("ph2.idmva"                , &g_mvavars[9]);
  //get weight file
  //g_reader->BookMVA("Ystar","/scratch/bendavid/root/hggmvaJan2/HggBambu_SM_Jan2_BDTG.weights.xml");
  g_reader->BookMVA("BDTG","/home/mingyang/cms/root/RootFiles/hggmva_2012_jan16/diphotonMVA/weights/HggBambu_SMDipho_Jan16_BDTG.weights.xml");  

  //---define RooRealVar---
  //variables already in the input ntuple
  RooRealVar hmass("mass","m_{#gamma#gamma}",100.0,180.0,"GeV");
  hmass.setBins(320);
  RooRealVar masserrsmeared("masserrsmeared","#sigma_m#gamma#gamma",0.0);  
  RooRealVar masserrsmearedwrongvtx("masserrsmearedwrongvtx","#sigma_m#gamma#gamma",0.0);   
  RooRealVar vtxprob("vtxprob","",0.0);  
  RooRealVar ph1Cat("ph1.phcat","",0.0);
  RooRealVar ph2Cat("ph2.phcat","",0.0);
  RooRealVar ptgg("ptgg","",0.0);
  RooRealVar higgspt("genHiggspt","",0.0);
  RooRealVar higgsZ("genHiggsZ","",0.0);
  RooRealVar vtxZ("vtxZ","",0.0);
  RooRealVar numPU("numPU","",0.0);
  RooRealVar ismc("ismc","",0.0);
  RooRealVar ph1isbarrel("ph1.isbarrel","",0.0);
  RooRealVar ph2isbarrel("ph2.isbarrel","",0.0);
  RooRealVar ph1sceta("ph1.sceta","",0.0);
  RooRealVar ph2sceta("ph2.sceta","",0.0);  
  RooRealVar ph1scphi("ph1.scphi","",0.0);
  RooRealVar ph2scphi("ph2.scphi","",0.0);
  RooRealVar ph1eta("ph1.eta","",0.0);
  RooRealVar ph2eta("ph2.eta","",0.0);  
  RooRealVar ph1phi("ph1.phi","",0.0);
  RooRealVar ph2phi("ph2.phi","",0.0);    
  RooRealVar ph1e("ph1.e","",0.0);
  RooRealVar ph2e("ph2.e","",0.0);  
  RooRealVar ph1eerr("ph1.eerr","",0.0);
  RooRealVar ph2eerr("ph2.eerr","",0.0);  
  RooRealVar ph1esmearing("ph1.esmearing","",0.0);
  RooRealVar ph2esmearing("ph2.esmearing","",0.0);  
  RooRealVar ph1pt("ph1.pt","",0.0);
  RooRealVar ph2pt("ph2.pt","",0.0);  
  RooRealVar ph1idmva("ph1.idmva","",0.0);
  RooRealVar ph2idmva("ph2.idmva","",0.0);   
  RooRealVar ph1r9("ph1.r9","",0.0);
  RooRealVar ph2r9("ph2.r9","",0.0);    
  RooRealVar jet1pt("jet1pt","",0.0);
  RooRealVar jet2pt("jet2pt","",0.0);
  RooRealVar jet1eta("jet1eta","",0.0);
  RooRealVar jet2eta("jet2eta","",0.0);  
  RooRealVar dijetmass("dijetmass","",0.0);  
  RooRealVar zeppenfeld("zeppenfeld","",0.0);  
  RooRealVar dphidijetgg("dphidijetgg","",0.0);  
  RooRealVar pfmet("pfmet","",0.0);  

  RooArgSet varlist;
  varlist.add(hmass);
  varlist.add(masserrsmeared);  
  varlist.add(masserrsmearedwrongvtx);  
  varlist.add(vtxprob);    
  varlist.add(ph1Cat);
  varlist.add(ph2Cat);
  varlist.add(ptgg);
  varlist.add(higgspt);
  varlist.add(higgsZ);
  varlist.add(vtxZ);
  varlist.add(numPU);
  //varlist.add(ismc);
  //varlist.add(ph1isbarrel);
  //varlist.add(ph2isbarrel);
  varlist.add(ph1sceta);
  varlist.add(ph2sceta);  
  varlist.add(ph1scphi);
  varlist.add(ph2scphi);    
  varlist.add(ph1eta);
  varlist.add(ph2eta);  
  varlist.add(ph1phi);
  varlist.add(ph2phi);    
  varlist.add(ph1e);
  varlist.add(ph2e);      
  varlist.add(ph1eerr);
  varlist.add(ph2eerr);        
  varlist.add(ph1esmearing);
  varlist.add(ph2esmearing);          
  varlist.add(ph1pt);
  varlist.add(ph2pt);
  varlist.add(ph1idmva);
  varlist.add(ph2idmva);
  varlist.add(ph1r9);
  varlist.add(ph2r9);  
  varlist.add(jet1pt);
  varlist.add(jet2pt);
  varlist.add(jet1eta);
  varlist.add(jet2eta);
  varlist.add(dijetmass);
  varlist.add(zeppenfeld);
  varlist.add(dphidijetgg);
  varlist.add(pfmet);
  
  //extra variables to be added which not included in the input ntuple yet
  RooRealVar xsecweight("xsecweight","xsecweight",1.0);
  RooRealVar procidx("procidx","",0.0);  
  
  RooArgSet varlistw = varlist;
  varlistw.add(xsecweight);
  varlistw.add(procidx);
  
  //RooFormulaVar newmassf computed from other variables
  RooFormulaVar newmassf("CMS_hgg_mass","m_{#gamma#gamma}","mass",RooArgList(hmass));//ming:the reason to add this newmassf is to change the name of the mass varibale to match the name for the combination datacards
  RooFormulaVar bdtf("bdt","","getmvaval(masserrsmeared/mass,masserrsmearedwrongvtx/mass,vtxprob,ph1.pt/mass,ph2.pt/mass,ph1.eta,ph2.eta,TMath::Cos(ph1.phi-ph2.phi),ph1.idmva,ph2.idmva)",varlist);
  
  //---make the RooDataSet contain data---
  //get the input tree
  //TFile *datafile = new TFile("/home/bendavid/cms/hist/hgg-v0/MergedPho2011.root","READ");
  TFile *datafile = new TFile("/home/mingyang/cms/hist/hgg-v0/merged/hgg-v0_r11-pho-j16-v1_noskim.root","READ");//ming:data
  TDirectory *datadir = (TDirectory*) datafile->FindObjectAny(modname);
  TTree *hdata = (TTree*)datadir->Get("hPhotonTree");  
  //make RooDataSet from the input tree passing loosesel and only selecting branches within varlist 
  RooDataSet calldata("calldata","",varlist,RooFit::Import(*hdata),RooFit::Cut(loosesel));
  //add new variables bdtf and newmassf
  RooRealVar *bdt = (RooRealVar*)calldata.addColumn(bdtf);
  //the step below on newmass is to set the range for fitting
  RooRealVar *newmass = (RooRealVar*)calldata.addColumn(newmassf);
  newmass->setRange(100.0,180.0);
  newmass->setUnit("GeV");
  newmass->setBins(320);
  //make some control plot: the bdt distribution of calldata
  RooPlot *dplot = bdt->frame(-1.0,1.0,100);
  calldata.plotOn(dplot);
  dplot->Draw();
  //return;
  
  //---make RooWorkspace---
  RooWorkspace *w = 0;
  RooWorkspace *wdata = new RooWorkspace("cms_hgg_workspace_data","") ;
  RooWorkspace *wmclow = new RooWorkspace("cms_hgg_workspace_mclow","") ;//ming:we split mc to wmclow and wmchigh to decrease the size of the workspace to not hit the internal roofit limit
  RooWorkspace *wmchigh = new RooWorkspace("cms_hgg_workspace_mchigh","") ;

  //import IntLumi to workspace
  RooRealVar *IntLumi = new RooRealVar("IntLumi","",totallumi);
  IntLumi->setConstant();
  wdata->import(*IntLumi);
  wmclow->import(*IntLumi);
  
  //import cdata to wdata
  for (int icat=0; icat<catcuts.size(); ++icat) {
    //ming: make cdata from calldata with category cut 
    RooDataSet cdata(TString::Format("data_mass_%s",catnames.at(icat).Data()),"",*calldata.get(),RooFit::Import(calldata),RooFit::Cut(catcuts.at(icat)));
    wdata->import(cdata);   
  }
  
  //define mass points and cross sections
  //ming: we choose 110,115,120,130,140,150 at this moment, and interpolate 125 from 120 and 130; 115 is choosen here to be the test of interpolation between 110 and 120
  //ming: for the next step, we could choose the mass points using 5 GeV steps instead
  std::vector<int> mhs;//mass points
  std::vector<double> gfxsecs;//gf cross sections with h->gg BR multiplied
  std::vector<double> vbfxsecs;//vbf cross sections with h->gg BR multiplied 
  std::vector<double> vtthxsecs;//vh plus tth cross sections with h->gg BR multiplied;we use vtthxsecs and vhxsecs here instead of vh and tth is because some historical reason: the vh and tth sample used to be merged
  std::vector<double> vhxsecs;//vh cross sections with h->gg BR multiplied
  
  mhs.push_back(110);//ming:the cross section numbers listed below are not all right
  gfxsecs.push_back(0.03908);//cross section*BR=19.84pb*0.00197=0.03908 pb
  vbfxsecs.push_back(0.00275);
  vtthxsecs.push_back(0.00290);
  vhxsecs.push_back(0.00265);
  //ffxsecs.push_back(0.08323692+0.08023015);//ming: this should be the fermiphobic cross section
  
  mhs.push_back(115);
  gfxsecs.push_back(0.03862);
  vbfxsecs.push_back(0.00284);
  vtthxsecs.push_back(0.00272);  
  vhxsecs.push_back(0.00248);
  //ffxsecs.push_back(0.0481518+0.042125595);
  
  mhs.push_back(120);
  gfxsecs.push_back(0.03742);
  vbfxsecs.push_back(0.00286);
  vtthxsecs.push_back(0.00251);
  vhxsecs.push_back(0.00229);
  //ffxsecs.push_back(0.02927583+0.023436813);
  
  mhs.push_back(130);
  gfxsecs.push_back(0.03191);
  vbfxsecs.push_back(0.00261);
  vtthxsecs.push_back(0.00193);  
  vhxsecs.push_back(0.00176);
  //ffxsecs.push_back(0.01224394+0.008260946);
    
  mhs.push_back(140);
  gfxsecs.push_back(0.02353);
  vbfxsecs.push_back(0.00204);
  vtthxsecs.push_back(0.00129);    
  vhxsecs.push_back(0.00117);
  //ffxsecs.push_back(0.005656604+0.0032417933);
  
  mhs.push_back(150);
  gfxsecs.push_back(0.01428);
  vbfxsecs.push_back(0.001307912);
  vtthxsecs.push_back(0.0007073224);    
  vhxsecs.push_back(0.000641104);  
  
  std::vector<double> effaccs;
  
  //loop over mass points
  for (UInt_t i=0; i<mhs.size(); ++i) {
    //split workspace for high mass and low mass
    int masspoint = mhs.at(i);
    if (masspoint<=120) w = wmclow;
    else w = wmchigh;
    
    //load mass dependent pt weights
    ptweights= (TH1D*) fileptweight->Get(TString::Format("kfact%d_0",masspoint));
    
    //get input 
    //TString samplestring = "/home/bendavid/cms/hist/hgg-v0/merged/hgg-v0_f11--h%dgg-%s-v14b-pu_noskim.root";
    TString samplestring = "/home/mingyang/cms/hist/hgg-v0/merged/hgg-v0_f11--h%dgg-%s-v14b-pu_noskim.root";
    TString vhsamplestring = samplestring;
    if (masspoint==110 || masspoint==130 || masspoint==140) vhsamplestring.ReplaceAll("f11--","f11-"); 
    TString ttsamplestring = samplestring;
    
    TFile *fgf = new TFile(TString::Format(samplestring,masspoint,"gf"),"READ");
    TFile *fvbf = new TFile(TString::Format(samplestring,masspoint,"vbf"),"READ");
    TFile *fvh = new TFile(TString::Format(vhsamplestring,masspoint,"vh"),"READ");
    TFile *ftt = new TFile(TString::Format(ttsamplestring,masspoint,"tt"),"READ");
    
    setpuweights(fgf,hpuest,0);
    setpuweights(fvbf,hpuest,1);
    setpuweights(fvh,hpuest,2);
    if (dott) setpuweights(ftt,hpuest,3);
    
    //define aggregate weight, so far using xsec, pileup, pt-reweighting and efficiency scale factors
    RooArgList weightvarlist(*IntLumi,xsecweight,numPU,higgspt,procidx,ph1Cat,ph2Cat,ph1pt,ph2pt);
    weightvarlist.add(ph1r9);
    weightvarlist.add(ph2r9);
    
    RooFormulaVar totweight("totweight","totweight","IntLumi*xsecweight*puweight(numPU,procidx)*ptweight(genHiggspt,procidx)*effweight(ph1.phcat,ph1.r9)*effweight(ph2.phcat,ph2.r9)*trigeffweight(ph1.phcat,ph2.phcat,ph1.pt,ph2.pt)",weightvarlist);
    
    //import cross section for each mass point to workspace
    double totalxsec = 0.0;
    if (!doff) totalxsec += gfxsecs.at(i);
    totalxsec += vbfxsecs.at(i);
    totalxsec += vhxsecs.at(i);
    if (dott && !doff) totalxsec += vtthxsecs.at(i)-vhxsecs.at(i);//ming:what does vtthxsecs mean?

    RooRealVar *XSBR = new RooRealVar(TString::Format("XSBR_%d",masspoint),"",totalxsec);//total cross section
    XSBR->setConstant();
    w->import(*XSBR);  

    RooRealVar *XSBR_gf = new RooRealVar(TString::Format("XSBR_ggh_%d",masspoint),"",gfxsecs.at(i));//gf cross section
    XSBR_gf->setConstant();
    w->import(*XSBR_gf);    

    RooRealVar *XSBR_vbf = new RooRealVar(TString::Format("XSBR_vbf_%d",masspoint),"",vbfxsecs.at(i));//vbf cross section
    XSBR_vbf->setConstant();
    w->import(*XSBR_vbf);    
    
    RooRealVar *XSBR_vh = new RooRealVar(TString::Format("XSBR_wzh_%d",masspoint),"",vhxsecs.at(i));//vh cross section
    XSBR_vh->setConstant();
    w->import(*XSBR_vh);    
    
    RooRealVar *XSBR_tt = new RooRealVar(TString::Format("XSBR_tth_%d",masspoint),"",vtthxsecs.at(i)-vhxsecs.at(i));//tt cross section
    XSBR_tt->setConstant();
    w->import(*XSBR_tt);        
    
    TCut rightvtx = "abs(genHiggsZ-vtxZ)<1.0";
    TCut wrongvtx = "abs(genHiggsZ-vtxZ)>=1.0";
    TCut rightcut = loosesel && rightvtx;
    TCut wrongcut = loosesel && wrongvtx;
    TCut allcut = loosesel;
    
    //make RooDataSet for mc with right vertex, wrong vertex and all vertex
    //right vertex
    RooDataSet mcsigdata(TString::Format("mcsigdata_m%d",masspoint),"",varlistw);    
    if (!doff) append(mcsigdata,fgf,rightcut,gfxsecs.at(i), 0, modname);//append fgf to mcsigdata with xsecweight=gfxsecs.at(i)/hDAllEvents->GetEntries() and procidx=0 added
    append(mcsigdata,fvbf,rightcut,vbfxsecs.at(i), 1, modname);
    append(mcsigdata,fvh,rightcut,vhxsecs.at(i), 2, modname);  
    if (dott && !doff) append(mcsigdata,ftt,rightcut,vtthxsecs.at(i)-vhxsecs.at(i), 3, modname);  
    mcsigdata.addColumn(totweight);
    mcsigdata.addColumn(bdtf);
    mcsigdata.addColumn(newmassf);
    
    //wrong vertex
    RooDataSet mcsigwrongdata(TString::Format("mcsigwrongdata_m%d",masspoint),"",varlistw);    
    if (!doff) append(mcsigwrongdata,fgf,wrongcut,gfxsecs.at(i), 0, modname);
    append(mcsigwrongdata,fvbf,wrongcut,vbfxsecs.at(i), 1, modname);
    append(mcsigwrongdata,fvh,wrongcut,vhxsecs.at(i), 2, modname);  
    if (dott && !doff) append(mcsigwrongdata,ftt,wrongcut,vtthxsecs.at(i)-vhxsecs.at(i), 3, modname);  
    mcsigwrongdata.addColumn(totweight);
    mcsigwrongdata.addColumn(bdtf);
    mcsigwrongdata.addColumn(newmassf);

    //combined
    RooDataSet mcsigalldata(TString::Format("mcsigalldata_m%d",masspoint),"",varlistw);    
    if (!doff) append(mcsigalldata,fgf,allcut,gfxsecs.at(i), 0, modname);
    append(mcsigalldata,fvbf,allcut,vbfxsecs.at(i), 1, modname);
    append(mcsigalldata,fvh,allcut,vhxsecs.at(i), 2, modname);
    if (dott && !doff) append(mcsigalldata,ftt,allcut,vtthxsecs.at(i)-vhxsecs.at(i), 3, modname);  
    mcsigalldata.addColumn(totweight);    
    mcsigalldata.addColumn(bdtf);    
    mcsigalldata.addColumn(newmassf);
      
    //make control plot  
    if (masspoint==120) {
      new TCanvas;
      RooPlot *dplota = bdt->frame(-1.0,1.0,100);
      calldata.plotOn(dplota,Rescale(1.0/calldata.sumEntries()));
      mcsigalldata.plotOn(dplota,MarkerColor(kRed),Rescale(1.0/mcsigalldata.sumEntries()));
      dplota->Draw();
      //return;
    }
    
    //loop over all categories
    for (UInt_t icat=0; icat<catcuts.size(); ++icat) {

      TString catname = catnames.at(icat);
 
      //---import RooAbsData for each category each mechanism for right vertex, wrong vertex and all vertex

      //import RooAbsData for each category each mechanism for right vertex
      RooDataSet *mcsigwdata = new RooDataSet(TString::Format("sig_mass_m%d_rv_%s",masspoint,catname.Data()),"",*mcsigdata.get(),RooFit::Import(mcsigdata),RooFit::Cut(catcuts.at(icat)),RooFit::WeightVar("totweight"));//e.g. sig_mass_m110_rv_cat0 
      RooAbsData *mcsigwdataggh = mcsigwdata->reduce("procidx==0");
      mcsigwdataggh->SetName(TString::Format("sig_ggh_mass_m%d_rv_%s",masspoint,catname.Data()));//e.g. sig_ggh_mass_m110_rv_cat0
      RooAbsData *mcsigwdatavbf = mcsigwdata->reduce("procidx==1");
      mcsigwdatavbf->SetName(TString::Format("sig_vbf_mass_m%d_rv_%s",masspoint,catname.Data()));
      RooAbsData *mcsigwdatawzh = mcsigwdata->reduce("procidx==2");
      mcsigwdatawzh->SetName(TString::Format("sig_wzh_mass_m%d_rv_%s",masspoint,catname.Data()));
      RooAbsData *mcsigwdatatth = mcsigwdata->reduce("procidx==3");
      mcsigwdatatth->SetName(TString::Format("sig_tth_mass_m%d_rv_%s",masspoint,catname.Data()));      
      w->import(*mcsigwdata);
      w->import(*mcsigwdataggh);
      w->import(*mcsigwdatavbf);
      w->import(*mcsigwdatawzh);
      w->import(*mcsigwdatatth);
      
      //import RooAbsData for each category each mechanism for wrong vertex
      RooDataSet *mcsigwrongwdata = new RooDataSet(TString::Format("sig_mass_m%d_wv_%s",masspoint,catname.Data()),"",*mcsigwrongdata.get(),RooFit::Import(mcsigwrongdata),RooFit::Cut(catcuts.at(icat)),RooFit::WeightVar("totweight"));     
      RooAbsData *mcsigwrongwdataggh = mcsigwrongwdata->reduce("procidx==0");
      mcsigwrongwdataggh->SetName(TString::Format("sig_ggh_mass_m%d_wv_%s",masspoint,catname.Data()));
      RooAbsData *mcsigwrongwdatavbf = mcsigwrongwdata->reduce("procidx==1");
      mcsigwrongwdatavbf->SetName(TString::Format("sig_vbf_mass_m%d_wv_%s",masspoint,catname.Data()));
      RooAbsData *mcsigwrongwdatawzh = mcsigwrongwdata->reduce("procidx==2");
      mcsigwrongwdatawzh->SetName(TString::Format("sig_wzh_mass_m%d_wv_%s",masspoint,catname.Data()));
      RooAbsData *mcsigwrongwdatatth = mcsigwrongwdata->reduce("procidx==3");
      mcsigwrongwdatatth->SetName(TString::Format("sig_tth_mass_m%d_wv_%s",masspoint,catname.Data()));      
      w->import(*mcsigwrongwdata);
      w->import(*mcsigwrongwdataggh);
      w->import(*mcsigwrongwdatavbf);
      w->import(*mcsigwrongwdatawzh);
      w->import(*mcsigwrongwdatatth);
      
      //import RooAbsData for each category each mechanism for all vertex
      RooDataSet *mcsigallwdata = new RooDataSet(TString::Format("sig_mass_m%d_%s",masspoint,catname.Data()),"",*mcsigalldata.get(),RooFit::Import(mcsigalldata),RooFit::Cut(catcuts.at(icat)),RooFit::WeightVar("totweight"));     
      RooAbsData *mcsigallwdataggh = mcsigallwdata->reduce("procidx==0");
      mcsigallwdataggh->SetName(TString::Format("sig_ggh_mass_m%d_%s",masspoint,catname.Data()));
      RooAbsData *mcsigallwdatavbf = mcsigallwdata->reduce("procidx==1");
      mcsigallwdatavbf->SetName(TString::Format("sig_vbf_mass_m%d_%s",masspoint,catname.Data()));
      RooAbsData *mcsigallwdatawzh = mcsigallwdata->reduce("procidx==2");
      mcsigallwdatawzh->SetName(TString::Format("sig_wzh_mass_m%d_%s",masspoint,catname.Data()));
      RooAbsData *mcsigallwdatatth = mcsigallwdata->reduce("procidx==3");
      mcsigallwdatatth->SetName(TString::Format("sig_tth_mass_m%d_%s",masspoint,catname.Data()));      
      w->import(*mcsigallwdata);
      w->import(*mcsigallwdataggh);
      w->import(*mcsigallwdatavbf);
      w->import(*mcsigallwdatawzh);
      w->import(*mcsigallwdatatth);

      //---get some control information to print---
      //compute acceptance*efficiency and right vertex fraction
      double eaccnum = mcsigwdata->sumEntries()+mcsigwrongwdata->sumEntries();//total number of events passing the accpetance and selection in this category at this mass
      double eaccden = totalxsec*totallumi;//total number of events produced at this mass
      //double eaccden = (vbfxsecs.at(i)+vhxsecs.at(i))*1e6;

      double eacc = eaccnum/eaccden;//acceptance*efficiency
      double eaccerrlo = TEfficiency::ClopperPearson(Int_t(eaccden), Int_t(eaccnum), 0.683, kFALSE) - eacc;//ming:error on the eacc?
      double eaccerrhi = TEfficiency::ClopperPearson(Int_t(eaccden), Int_t(eaccnum), 0.683, kTRUE) - eacc;//ming:erro on the eacc?
      
      printf("effacc = %5e\n",eacc);
      
      //return;
      
      effaccs.push_back(eacc);
      
      //double fright = mcsigwdata->sumEntries()/(mcsigwdata->sumEntries()+mcsigwrongwdata->sumEntries());
      double frightnum = mcsigwdata->sumEntries();//number of events with right vertex in this category
      double frightden = mcsigwdata->sumEntries()+mcsigwrongwdata->sumEntries();//total number of events in this category
      double fright = frightnum/frightden;//fraction of events with right vtx
      double frighterrlo = TEfficiency::ClopperPearson(Int_t(frightden), Int_t(frightnum), 0.683, kFALSE) - fright;
      double frighterrhi = TEfficiency::ClopperPearson(Int_t(frightden), Int_t(frightnum), 0.683, kTRUE) - fright;
    }

    //---print some control information---
    double toteffacc = 0.0;
    for (int icat=0; icat<catnames.size(); ++icat) {
      printf("%s: %5f\n",catnames.at(icat).Data(),effaccs.at(icat));
      toteffacc+=effaccs.at(icat);
    }
    printf("total effacc = %5f\n",toteffacc);
    
  }
  
  //save everything to file with RooWorkspace
  //w->import(fullsigpdf);
  //w->import(*combhvtxslides[0]);
  //w->Print();
  //w->writeToFile("CMS-HGG-mchigh.root") ;
  //wdata->writeToFile("CMS-HGG-data-2012jan16_test.root") ;
  //wmclow->writeToFile("CMS-HGG-mclow-2012jan16_test.root") ;
  //wmchigh->writeToFile("CMS-HGG-mchigh-2012jan16_test.root") ;
  //   TFile *foutfile = new TFile("CMS-HGG.root","RECREATE");
  //   w->Write();
  //   
  //   foutfile->Close();

  TString dataFile = TString("CMS-HGG-data-2012jan16"+CMEnergy+".root");
  TString mclowFile = TString("CMS-HGG-mclow-2012jan16"+CMEnergy+".root");
  TString mchighFile = TString("CMS-HGG-mchigh-2012jan16"+CMEnergy+".root"); 
  wdata->writeToFile(dataFile) ;
  wmclow->writeToFile(mclowFile) ;
  wmchigh->writeToFile(mchighFile) ; 
  
  return; 
  
}
