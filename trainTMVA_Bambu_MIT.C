#include <cstdlib>
#include <iostream> 
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TPluginManager.h"

//#include "TMVAGui.C"

#if not defined(__CINT__) || defined(__MAKECINT__)
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"

#include "TH1D.h"


#endif

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


TH1D *puweights[11];
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
float effweight(int cat) {
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

  //2012 Jan16 freeze numbers(ming: accroding to Matteo's email
  //scale factor for the effieicnecy of the single photon
  //if (cat==1) return 0.999;
  //else if (cat==2) return 0.984;
  //else if (cat==3) return 1.006;
  //else if (cat==4) return 1.014;
  //else return 1.0;  

  //2012 ICHEP Open Box
  //scale factor for the effieicnecy of the single photon
  if (cat==1) return 1.0;
  else if (cat==2) return 1.0;
  else if (cat==3) return 1.0;
  else if (cat==4) return 1.0;
  else return 1.0;  
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
  //else return 1.0;

  if (isb && isr9) return 1.0;
  else if (isb && !isr9) return 1.0;
  else if (!isb && isr9) return 1.0;
  else if (!isb && !isr9) return 1.0;
  else return 1.0;
  
}

void trainTMVA_Bambu_MIT() {
  
  //-------------set up the weights-------------------
  //pu target
  //TFile *filepuest = new TFile("/scratch/bendavid/root/puweightsNov13/2011_0100_73500.pileup.root","READ");
  //TFile *filepuest = new TFile("/home/mingyang/cms/puweight/augmented_nov08_rereco.json.68000.pileup.root","READ");
  //TH1D *hpuest = (TH1D*) filepuest->Get("pileup");
  
  //load higgs pt-weights
  //TFile *fileptweight = new TFile("/scratch/bendavid/root/KFactors_AllScales.root","READ");
  //TFile *fileptweight = new TFile("./root/KFactors_AllScales.root","READ");
  //ptweights= (TH1D*) fileptweight->Get("kfact125_0");//ming: how about higgs 123  

  //-------------get input tree-------------------
  TString prefix = "/home/mingyang/cms/root/RootFiles/2012_ICHEP/DiphotonMVA/";
  TString rootfile ="diphomva_training_dijet.root";
  
  std::cout<<"  loading files... "<<std::endl;

  TFile* f_train_file = TFile::Open((prefix+rootfile).Data());

  std::cout<<"  ... done"<<std::endl;

  TTree* t_signal_1 = (TTree*) f_train_file->Get((TString("gluglu_H_gg_124_pu2012")).Data());
  TTree* t_signal_2 = (TTree*) f_train_file->Get((TString("vbf_H_gg_124_pu2012")).Data());
  TTree* t_signal_3 = (TTree*) f_train_file->Get((TString("wz_H_gg_124_pu2012")).Data());
  TTree* t_signal_4 = (TTree*) f_train_file->Get((TString("tt_H_gg_124_pu2012")).Data());

  TTree* t_bgBorn   = (TTree*) f_train_file->Get((TString("diphotonjets")).Data());
  TTree* t_bgBox_1  = (TTree*) f_train_file->Get((TString("box_10")).Data());
  TTree* t_bgBox_2  = (TTree*) f_train_file->Get((TString("box_25")).Data());
  TTree* t_bgBox_3  = (TTree*) f_train_file->Get((TString("box_250")).Data());

  TTree* t_bgQCD30_pf = (TTree*) f_train_file->Get((TString("qcd_pf_30")).Data());
  TTree* t_bgQCD30_ff = (TTree*) f_train_file->Get((TString("qcd_ff_30")).Data()); 
  TTree* t_bgQCD40_pf = (TTree*) f_train_file->Get((TString("qcd_pf_40")).Data());
  TTree* t_bgQCD40_ff = (TTree*) f_train_file->Get((TString("qcd_ff_40")).Data()); 
  
  TTree* t_bgGJET20_pf = (TTree*) f_train_file->Get((TString("gjet_pf_20")).Data());
  TTree* t_bgGJET40_pf = (TTree*) f_train_file->Get((TString("gjet_pf_40")).Data());

  //-------------add new branch procidx to input tree-------------------
  float iproc = 0.;
  TBranch* procB_signal_1 = t_signal_1->Branch("procidx", &iproc, "procidx/F");//define new branch
  for(int i=1; i<=t_signal_1->GetEntries(); ++i)
    procB_signal_1->Fill();//fill the new branch

  iproc = 1.;
  TBranch* procB_signal_2 = t_signal_2->Branch("procidx", &iproc, "procidx/F");
  for(int i=1; i<=t_signal_2->GetEntries(); ++i)
    procB_signal_2->Fill();

  iproc = 2.;
  TBranch* procB_signal_3 = t_signal_3->Branch("procidx", &iproc, "procidx/F");
  for(int i=1; i<=t_signal_3->GetEntries(); ++i)
    procB_signal_3->Fill();

  iproc = 3.;
  TBranch* procB_signal_4 = t_signal_4->Branch("procidx", &iproc, "procidx/F");
  for(int i=1; i<=t_signal_4->GetEntries(); ++i)
    procB_signal_4->Fill();
  
  iproc = 4.;
  TBranch* procB_bgBorn = t_bgBorn->Branch("procidx", &iproc, "procidx/F");
  for(int i=1; i<=t_bgBorn->GetEntries(); ++i)
    procB_bgBorn->Fill();

  iproc = 5.;
  TBranch* procB_bgBox_1 = t_bgBox_1->Branch("procidx", &iproc, "procidx/F");
  for(int i=1; i<=t_bgBox_1->GetEntries(); ++i)
    procB_bgBox_1->Fill();

  iproc = 6.;
  TBranch* procB_bgBox_2 = t_bgBox_2->Branch("procidx", &iproc, "procidx/F");
  for(int i=1; i<=t_bgBox_2->GetEntries(); ++i)
    procB_bgBox_2->Fill();

  iproc = 7.;
  TBranch* procB_bgBox_3 = t_bgBox_3->Branch("procidx", &iproc, "procidx/F");
  for(int i=1; i<=t_bgBox_3->GetEntries(); ++i)
    procB_bgBox_3->Fill();

  iproc = 8.;
  TBranch* procB_bgQCD30_pf = t_bgQCD30_pf->Branch("procidx", &iproc, "procidx/F");
  for(int i=1; i<=t_bgQCD30_pf->GetEntries(); ++i)
    procB_bgQCD30_pf->Fill();

  iproc = 9.;
  TBranch* procB_bgQCD30_ff = t_bgQCD30_ff->Branch("procidx", &iproc, "procidx/F");
  for(int i=1; i<=t_bgQCD30_ff->GetEntries(); ++i)
    procB_bgQCD30_ff->Fill();

  iproc = 10.;
  TBranch* procB_bgQCD40_pf = t_bgQCD40_pf->Branch("procidx", &iproc, "procidx/F");
  for(int i=1; i<=t_bgQCD40_pf->GetEntries(); ++i)
    procB_bgQCD40_pf->Fill();
  
  iproc = 11.;
  TBranch* procB_bgQCD40_ff = t_bgQCD40_ff->Branch("procidx", &iproc, "procidx/F");
  for(int i=1; i<=t_bgQCD40_ff->GetEntries(); ++i)
    procB_bgQCD40_ff->Fill();
   
  iproc = 12.;
  TBranch* procB_bgGJET20_pf = t_bgGJET20_pf->Branch("procidx", &iproc, "procidx/F");
  for(int i=1; i<=t_bgGJET20_pf->GetEntries(); ++i)
    procB_bgGJET20_pf->Fill();

  iproc = 13.;
  TBranch* procB_bgGJET40_pf = t_bgGJET40_pf->Branch("procidx", &iproc, "procidx/F");
  for(int i=1; i<=t_bgGJET40_pf->GetEntries(); ++i)
    procB_bgGJET40_pf->Fill();

  //----------------------Train MVA-------------------------------------------
  TFile* tmvaOutput = new TFile("/home/mingyang/cms/root/RootFiles/2012_ICHEP/DiphotonMVA/Tmva_TrainOutput_SMDipho_2012ICHEP_DijetTagCorrSmear_.root","RECREATE");
  TMVA::Factory* tmva = new TMVA::Factory("Hgg_SMDipho_2012ICHEP_DijetTagCorrSmear",tmvaOutput);
  
  //add input signal tree
  tmva->AddSignalTree(t_signal_1);
  tmva->AddSignalTree(t_signal_2);
  tmva->AddSignalTree(t_signal_3);
  tmva->AddSignalTree(t_signal_4); 
  
  //add input background tree
  tmva->AddBackgroundTree(t_bgBorn);
  tmva->AddBackgroundTree(t_bgBox_1);
  tmva->AddBackgroundTree(t_bgBox_2);
  tmva->AddBackgroundTree(t_bgBox_3);

  tmva->AddBackgroundTree(t_bgQCD30_pf);
  tmva->AddBackgroundTree(t_bgQCD30_ff);
  tmva->AddBackgroundTree(t_bgQCD40_pf);
  tmva->AddBackgroundTree(t_bgQCD40_ff);

  tmva->AddBackgroundTree(t_bgGJET20_pf);  
  tmva->AddBackgroundTree(t_bgGJET40_pf);

  //add training variable
  tmva->AddVariable("sigmaMrvoM",    'D');
  tmva->AddVariable("sigmaMwvoM",    'D');
  tmva->AddVariable("vtxprob",    'D');  
  tmva->AddVariable("ptoM1" ,  'D');
  tmva->AddVariable("ptoM2" ,  'D');
  tmva->AddVariable("eta1" ,  'D');
  tmva->AddVariable("eta2" ,  'D');
  tmva->AddVariable("cosphi" ,  'D');
  tmva->AddVariable("idmva1" ,  'D');
  tmva->AddVariable("idmva2" ,  'D');

  //add signal weight: 1/lumi * 1/resolution
  tmva->SetSignalWeightExpression("xsec_weight*(vtxprob/sigmaMrvoM + (1.0-vtxprob)/sigmaMwvoM)"); 
  //add background weight: 1/lumi 
  tmva->SetBackgroundWeightExpression("xsec_weight");
  
  TCut selectionCut = "ptoM1 > 40./120 && ptoM2 > 30./120. && mass < 180. && mass > 100. && idmva1>-0.3 && idmva2>-0.3";
 
  TCut diphojCut = "!( (procidx==4) && ((event>=270000 && event<=390000) || (event>=480000 && event<=510000) || (event>=1170000 && event<=1200000) || (event>=1230000 && event<=1260000)  || (event>=1890000 && event<=1920000)  || (event>=2670000 && event<=2730000) || (event>=3090000 && event<=3210000) || (event>=3270000 && event<=3300000) || (event>=3510000 && event<=3600000) || (event>=4890000 && event<=5010000) || (event>=5340000 && event<=5370000) || (event>=5430000 && event<=5500000)))";
 
  TCut totalCut = selectionCut && diphojCut;

  tmva->PrepareTrainingAndTestTree(totalCut,"nTrain_Signal=0:nTrain_Background=0:nTest_Signal=20:nTest_Background=20");
  tmva->BookMethod( TMVA::Types::kBDT, "BDTG" ,
		    "!H:!V:NTrees=2000::BoostType=Grad:Shrinkage=0.1:UseBaggedGrad=F:nCuts=2000:MaxDepth=3:NNodesMax=100000:UseYesNoLeaf=F:nEventsMin=10000:");
  
  tmva->TrainAllMethods();
  tmva->TestAllMethods();
  tmva->EvaluateAllMethods();
  
  return;
}
