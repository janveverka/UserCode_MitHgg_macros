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

// MIT BDT Training for H->gg analysis

/*double trigEff[36];
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
trigEff[35] = 0.949  ;*/

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
  if (cat==1) return 0.999;
  else if (cat==2) return 0.984;
  else if (cat==3) return 1.006;
  else if (cat==4) return 1.014;
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

//ming: this setRelErrAlias seems not be used
void setRelErrAlias(TTree* tree) {
  TString  corfact = "0.07";
  TString relerrsmearcor = TString::Format("0.5*sqrt( (pow(1.0+ismc*ph1.isbarrel*%s,2)*ph1.eerr*ph1.eerr+ph1.esmearing*ph1.esmearing)/ph1.e/ph1.e + (pow(1.0+ismc*ph2.isbarrel*%s,2)*ph2.eerr*ph2.eerr+ph2.esmearing*ph2.esmearing)/ph2.e/ph2.e)",corfact.Data(),corfact.Data());

  tree->SetAlias("massrelerr",relerrsmearcor.Data());
  return;
}

void trainTMVA_Bambu_MIT() {
  
  //-------------set up the weights-------------------
  //pu target
  //TFile *filepuest = new TFile("/scratch/bendavid/root/puweightsNov13/2011_0100_73500.pileup.root","READ");
  TFile *filepuest = new TFile("/home/mingyang/cms/puweight/augmented_nov08_rereco.json.68000.pileup.root","READ");
  TH1D *hpuest = (TH1D*) filepuest->Get("pileup");

  
  //load higgs pt-weights
  TFile *fileptweight = new TFile("/scratch/bendavid/root/KFactors_AllScales.root","READ");
  //TFile *fileptweight = new TFile("./root/KFactors_AllScales.root","READ");
  ptweights= (TH1D*) fileptweight->Get("kfact125_0");//ming: how about higgs 123  

  //-------------get input tree-------------------
  //TString prefix = "rfio:///castor/cern.ch/user/f/fabstoec/BambuYstar_Oct14/";  
  TString prefix = "/home/mingyang/cms/hist/hgg-v0/merged/";

  TString signal_1 ="hgg-v0_f11--h123gg-gf-v14b-pu_noskim.root";
  TString signal_2 ="hgg-v0_f11--h123gg-vbf-v14b-pu_noskim.root";
  TString signal_3 ="hgg-v0_f11--h123gg-vh-v14b-pu_noskim.root";
  TString signal_4 ="hgg-v0_f11--h123gg-tt-v14b-pu_noskim.root";

  TString bgBorn   ="hgg-v0_f11--diphoj-v14b-pu_noskim.root";
  TString bgBox_1  ="hgg-v0_f11--2pibx10_25-v14b-pu_noskim.root";
  TString bgBox_2  ="hgg-v0_f11--2pibx25_250-v14b-pu_noskim.root";
  TString bgBox_3  ="hgg-v0_f11--2pibx250-v14b-pu_noskim.root";

  TString bgQCDl   ="hgg-v0_f11--qcd-2em3040-v14b-pu_noskim.root";
  TString bgQCD    ="hgg-v0_f11--qcd-2em40-v14b-pu_noskim.root";
  TString bgPJ     ="hgg-v0_f11-pj-2em20-v14b-pu_noskim.root";

  std::cout<<"  loading files... "<<std::endl;

  TFile* f_signal_1 = TFile::Open((prefix+signal_1).Data());
  TFile* f_signal_2 = TFile::Open((prefix+signal_2).Data());
  TFile* f_signal_3 = TFile::Open((prefix+signal_3).Data());
  TFile* f_signal_4 = TFile::Open((prefix+signal_4).Data());

  TFile* f_bgBorn   = TFile::Open((prefix+bgBorn).Data());
  TFile* f_bgBox_1  = TFile::Open((prefix+bgBox_1 ).Data());
  TFile* f_bgBox_2  = TFile::Open((prefix+bgBox_2 ).Data());
  TFile* f_bgBox_3  = TFile::Open((prefix+bgBox_3 ).Data());

  TFile* f_bgQCDl = TFile::Open((prefix+bgQCDl ).Data());
  TFile* f_bgQCD  = TFile::Open((prefix+bgQCD ).Data());
  TFile* f_bgPJ   = TFile::Open((prefix+bgPJ ).Data());

  std::cout<<"  ... done"<<std::endl;

  TString baseDir="RunLumiSelectionMod/MCProcessSelectionMod/HLTModP/GoodPVFilterMod/PhotonMvaMod/JetPub/JetCorrectionMod/PhotonPairSelectorPresel/PhotonTreeWriterPresel/";

  TTree* t_signal_1 = (TTree*) f_signal_1->Get((baseDir+TString("hPhotonTree")).Data());
  TTree* t_signal_2 = (TTree*) f_signal_2->Get((baseDir+TString("hPhotonTree")).Data());
  TTree* t_signal_3 = (TTree*) f_signal_3->Get((baseDir+TString("hPhotonTree")).Data());
  TTree* t_signal_4 = (TTree*) f_signal_4->Get((baseDir+TString("hPhotonTree")).Data());

  TTree* t_bgBorn   = (TTree*) f_bgBorn ->Get((baseDir+TString("hPhotonTree")).Data());
  TTree* t_bgBox_1  = (TTree*) f_bgBox_1->Get((baseDir+TString("hPhotonTree")).Data());
  TTree* t_bgBox_2  = (TTree*) f_bgBox_2->Get((baseDir+TString("hPhotonTree")).Data());
  TTree* t_bgBox_3  = (TTree*) f_bgBox_3->Get((baseDir+TString("hPhotonTree")).Data());

  TTree* t_bgQCDl = (TTree*) f_bgQCDl->Get((baseDir+TString("hPhotonTree")).Data());
  TTree* t_bgQCD  = (TTree*) f_bgQCD->Get((baseDir+TString("hPhotonTree")).Data());
  TTree* t_bgPJ   = (TTree*) f_bgPJ ->Get((baseDir+TString("hPhotonTree")).Data());

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
  TBranch* procB_bgQCDl = t_bgQCDl->Branch("procidx", &iproc, "procidx/F");
  for(int i=1; i<=t_bgQCDl->GetEntries(); ++i)
    procB_bgQCDl->Fill();
  iproc = 9.;
  TBranch* procB_bgQCD = t_bgQCD->Branch("procidx", &iproc, "procidx/F");
  for(int i=1; i<=t_bgQCD->GetEntries(); ++i)
    procB_bgQCD->Fill();  
  iproc = 10.;
  TBranch* procB_bgPJ = t_bgPJ->Branch("procidx", &iproc, "procidx/F");
  for(int i=1; i<=t_bgPJ->GetEntries(); ++i)
    procB_bgPJ->Fill();

  //--------------get PU weights-----------------------------------------
  puweights[0] = getpuweights(f_signal_1,hpuest);
  puweights[1] = getpuweights(f_signal_2,hpuest);
  puweights[2] = getpuweights(f_signal_3,hpuest);
  puweights[3] = getpuweights(f_signal_4,hpuest);

  puweights[4] = getpuweights(f_bgBorn,hpuest);
  puweights[5] = getpuweights(f_bgBox_1,hpuest);
  puweights[6] = getpuweights(f_bgBox_2,hpuest);
  puweights[7] = getpuweights(f_bgBox_3,hpuest);
  
  puweights[8] = getpuweights(f_bgQCDl,hpuest);
  puweights[9] = getpuweights(f_bgQCD,hpuest);
  puweights[10] = getpuweights(f_bgPJ,hpuest);
   
  //--------------get 1/Lumi_MC------------------------------------------
  TH1D* hNE_signal_1 = (TH1D*) f_signal_1->FindObjectAny("hDAllEvents");//ming: where's it?
  TH1D* hNE_signal_2 = (TH1D*) f_signal_2->FindObjectAny("hDAllEvents");
  TH1D* hNE_signal_3 = (TH1D*) f_signal_3->FindObjectAny("hDAllEvents");
  TH1D* hNE_signal_4 = (TH1D*) f_signal_4->FindObjectAny("hDAllEvents");

  //double w_signal_1 = 0.03742/hNE_signal_1->GetBinContent(1);//ming:pb
  //double w_signal_2 = 0.00286/hNE_signal_2->GetBinContent(1);
  //double w_signal_3 = 0.00229/hNE_signal_3->GetBinContent(1);
  //double w_signal_4 = 0.00022/hNE_signal_4->GetBinContent(1);
  
  double w_signal_1 = 0.03607/hNE_signal_1->GetBinContent(1);//ming:pb
  double w_signal_2 = 0.00281/hNE_signal_2->GetBinContent(1);
  double w_signal_3 = 0.00214/hNE_signal_3->GetBinContent(1);
  double w_signal_4 = 0.00021/hNE_signal_4->GetBinContent(1);

  TH1D* hNE_bgBorn = (TH1D*) f_bgBorn->FindObjectAny("hDAllEvents");
  double w_bgBorn  = 1.15*154.7/hNE_bgBorn->GetBinContent(1);

  TH1D* hNE_bgBox_1 = (TH1D*) f_bgBox_1->FindObjectAny("hDAllEvents");
  double w_bgBox_1  = 1.3*358.2/hNE_bgBox_1->GetBinContent(1);

  TH1D* hNE_bgBox_2 = (TH1D*) f_bgBox_2->FindObjectAny("hDAllEvents");
  double w_bgBox_2  = 1.3*12.37/hNE_bgBox_2->GetBinContent(1);

  TH1D* hNE_bgBox_3 = (TH1D*) f_bgBox_3->FindObjectAny("hDAllEvents");
  double w_bgBox_3  = 1.3*0.00021/hNE_bgBox_3->GetBinContent(1);

  TH1D* hNE_bgQCDl = (TH1D*) f_bgQCDl->FindObjectAny("hDAllEvents");
  double w_bgQCDl  = 1.3*10868.0/hNE_bgQCDl->GetBinContent(1);  
  
  TH1D* hNE_bgQCD = (TH1D*) f_bgQCD->FindObjectAny("hDAllEvents");
  double w_bgQCD  = 1.3*1.87E+07*0.0023/hNE_bgQCD->GetBinContent(1);//ming:0.0023 should be the double EM enriched filter efficiency

  TH1D* hNE_bgPJ = (TH1D*) f_bgPJ->FindObjectAny("hDAllEvents");
  double w_bgPJ  = 1.3*77100.*0.0065/hNE_bgPJ->GetBinContent(1);
 
  //----------------------Train MVA-------------------------------------------
  //TFile* tmvaOutput = new TFile("/scratch/bendavid/root/hggmvaJan2/tmvaOutputfpJan2.root","RECREATE");
  TFile* tmvaOutput = new TFile("/home/mingyang/cms/root/RootFiles/hggmva_2012_jan16/Tmva_TrainOutput_SMDipho_2012Jan16.root","RECREATE");//root file to store training information
  //TMVA::Factory* tmva = new TMVA::Factory("HggBambu_FP_Jan2",tmvaOutput);
  TMVA::Factory* tmva = new TMVA::Factory("HggBambu_SMDipho_Jan16",tmvaOutput);
  
  //add input signal tree
  tmva->AddSignalTree(t_signal_1,w_signal_1);//each input sample is normalized to unit luminosity
  tmva->AddSignalTree(t_signal_2,w_signal_2);
  tmva->AddSignalTree(t_signal_3,w_signal_3);
  tmva->AddSignalTree(t_signal_4,w_signal_4); 
  
  //add input background tree
  tmva->AddBackgroundTree(t_bgBorn,w_bgBorn);
  tmva->AddBackgroundTree(t_bgBox_1,w_bgBox_1);
  tmva->AddBackgroundTree(t_bgBox_2,w_bgBox_2);
  tmva->AddBackgroundTree(t_bgBox_3,w_bgBox_3);
  tmva->AddBackgroundTree(t_bgQCDl,w_bgQCDl);  
  tmva->AddBackgroundTree(t_bgQCD,w_bgQCD);
  tmva->AddBackgroundTree(t_bgPJ,w_bgPJ);

  //add training variable
  tmva->AddVariable("masserrsmeared/mass",    'F');
  tmva->AddVariable("masserrsmearedwrongvtx/mass",    'F');
  tmva->AddVariable("vtxprob",    'F');  
  tmva->AddVariable("ph1.pt/mass" ,  'F');
  tmva->AddVariable("ph2.pt/mass" ,  'F');
  tmva->AddVariable("ph1.eta" ,  'F');
  tmva->AddVariable("ph2.eta" ,  'F');
  tmva->AddVariable("TMath::Cos(ph1.phi-ph2.phi)" ,  'F');
  tmva->AddVariable("ph1.idmva" ,  'F');
  tmva->AddVariable("ph2.idmva" ,  'F');

  //tmva->AddVariable("TMath::Sqrt(ph1.pt^2+ph2.pt^2+2.*ph1.pt*ph2.pt*TMath::Cos(ph1.scphi-ph2.scphi))/mass" ,  'F');

  //add signal weight: 1/resolution*puweight*ptweight*presel eff scale factor* trigger eff
  tmva->SetSignalWeightExpression("(vtxprob*mass/masserrsmeared + (1.0-vtxprob)*mass/masserrsmearedwrongvtx)*puweight(numPU,procidx)*ptweight(genHiggspt,procidx)*effweight(ph1.phcat)*effweight(ph2.phcat)*trigeffweight(ph1.phcat,ph2.phcat,ph1.pt,ph2.pt)"); 
  //add background weight: it seems not necessary to include ptweight here...
  //ming: apply no scale factor to fake fake; not sure how the scale factors determined though
  tmva->SetBackgroundWeightExpression("((!ph1.ispromptgen && !ph2.ispromptgen)*1.0/1.3 + !(!ph1.ispromptgen && !ph2.ispromptgen)*1.0)*puweight(numPU,procidx)*ptweight(genHiggspt,procidx)*effweight(ph1.phcat)*effweight(ph2.phcat)*trigeffweight(ph1.phcat,ph2.phcat,ph1.pt,ph2.pt)");
  
  TCut selectionCut = "ph1.pt/mass > 40./120 && ph2.pt/mass > 30./120. && mass < 180. && mass > 100. && ph1.idmva>-0.3 && ph2.idmva>-0.3";
  //dec9id
  //tmva->PrepareTrainingAndTestTree(selectionCut,"nTrain_Signal=151800:nTrain_Background=89700:");
  //tmva->PrepareTrainingAndTestTree(selectionCut,"nTrain_Signal=95400:nTrain_Background=89700:");
  //tmva->PrepareTrainingAndTestTree(selectionCut,"SplitMode=Alternate");
  tmva->PrepareTrainingAndTestTree(selectionCut,"nTrain_Signal=0:nTrain_Background=0:nTest_Signal=10:nTest_Background=10");
  tmva->BookMethod( TMVA::Types::kBDT, "BDTG" ,
		    "!H:!V:NTrees=2000::BoostType=Grad:Shrinkage=0.1:UseBaggedGrad=F:nCuts=2000:MaxDepth=3:NNodesMax=100000:UseYesNoLeaf=F:nEventsMin=10000:");
  
  //tmva->BookMethod(TMVA::Types::kLD, "LD", "H:!V:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10");
  
  tmva->TrainAllMethods();
  tmva->TestAllMethods();
  tmva->EvaluateAllMethods();
  
  return;
}
