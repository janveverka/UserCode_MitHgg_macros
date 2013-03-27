// modified script from Josh for signal modeling

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

#include "TGraphAsymmErrors.h"


//---------------------------------------------------------------------------------------------------------------------------
// Include the header file with the model constans

#include "modelConstants_8TeV.h"
using namespace RooFit;

// --------------------------------------------------------------------------------------
// Weighting functions needed
// 1. PU weighting
TH1D* puweights;
TGraphAsymmErrors* vtxWeights_pass;
TGraphAsymmErrors* vtxWeights_fail;

// Graphs to potentially cottect IDMVA variables
TGraph* scaleIDMVA_EB_hR9;
TGraph* scaleIDMVA_EB_lR9;
TGraph* scaleIDMVA_EE_hR9;
TGraph* scaleIDMVA_EE_lR9;

//Graphs tom potentially correct SigEoE
TGraph* scaleSigEoE_0_10;
TGraph* scaleSigEoE_10_15;
TGraph* scaleSigEoE_15_20;
TGraph* scaleSigEoE_20_25;

float vtxWeight(float ptgg, float vtxz, float genHiggsZ ) {
  
  float diffVtx = TMath::Abs(vtxz - genHiggsZ);

  if (diffVtx < 0.1) {
    // this is the correct guy
    if (!vtxWeights_pass) return 1.0;
    return vtxWeights_pass->Eval(ptgg);
  }

  if (!vtxWeights_fail) return 1.0;
  return vtxWeights_fail->Eval(ptgg);

}


// BeamSpot reweighting
float bsweight( float vtxz, float genHiggsZ ) {

  float diffVtx = TMath::Abs(vtxz - genHiggsZ);

  if (diffVtx < 0.1) return 1.;
  
  float newBSmean1  = 9.9391e-02;
  float newBSmean2  = 1.8902e-01;
  float newBSnorm1  = 5.3210e+00;
  float newBSnorm2  = 4.1813e+01;
  float newBSsigma1 = 9.7530e-01;
  float newBSsigma2 = 7.0811e+00;

  float oldBSmean1  = 7.2055e-02;
  float oldBSmean2  = 4.9986e-01;
  float oldBSnorm1  = 3.5411e+00;
  float oldBSnorm2  = 4.0258e+01;
  float oldBSsigma1 = 7.9678e-01;
  float oldBSsigma2 = 8.5356e+00;

  float newBSgaus1 = newBSnorm1*TMath::Exp(-0.5*TMath::Power( (diffVtx - newBSmean1)/newBSsigma1 , 2));
  float newBSgaus2 = newBSnorm2*TMath::Exp(-0.5*TMath::Power( (diffVtx - newBSmean2)/newBSsigma2 , 2));
  float oldBSgaus1 = oldBSnorm1*TMath::Exp(-0.5*TMath::Power( (diffVtx - oldBSmean1)/oldBSsigma1 , 2));
  float oldBSgaus2 = oldBSnorm2*TMath::Exp(-0.5*TMath::Power( (diffVtx - oldBSmean2)/oldBSsigma2 , 2));

  float reweight = 1.1235*(newBSgaus1 + newBSgaus2)/(oldBSgaus1 + oldBSgaus2);

  return reweight;

}


//float puweight(Double_t npu, int wset=0) {
float puweight(Double_t npu) {

  //return 1.;

  if( !puweights ) return 1.;
  return puweights->GetBinContent(puweights->FindFixBin(npu));
}


void setpuweights(TFile *file, TH1D *target) {//ming:set puweights for each process; might not be necessary when the pu distributions for different mc are the same
  TDirectory *dirmcpv = (TDirectory*)file->FindObjectAny("AnaFwkMod");
  TH1D       *hnpu  = (TH1D*)dirmcpv->Get("hNPU");
//   TH1D       *hpumc = new TH1D("thNPU","",60,0.,60.);
//   for( int i=0; i<=61; ++i)
//     hpumc->SetBinContent(i,hnpu->GetBinContent(i));
  
  TH1D       *hpumc = (TH1D*)hnpu->Clone();
  
  hpumc->Sumw2();
  hpumc->Scale(1.0/hpumc->Integral(0,hpumc->GetNbinsX()+1));
  
  TH1D *htargettmp = new TH1D("htargettmp","", hpumc->GetNbinsX(), hpumc->GetXaxis()->GetXmin(), hpumc->GetXaxis()->GetXmax());
  htargettmp->Sumw2();
  for (int ibin = 0; ibin<=(htargettmp->GetNbinsX()+1); ++ibin) {
    htargettmp->Fill(htargettmp->GetBinCenter(ibin),target->GetBinContent(target->FindFixBin(htargettmp->GetBinCenter(ibin))));
  }
  htargettmp->Scale(1.0/htargettmp->Integral(0,htargettmp->GetNbinsX()+1));
  
  if( puweights )
    delete puweights;  // delete old guy, not needed anymore..

  puweights = new TH1D((*htargettmp)/(*hpumc));

  delete htargettmp;
  return;
}


// 2. PT weights
TH1D *ptweights = 0;
float ptweight(float genhpt) {

  if (genhpt<0 || !ptweights) return 1.0;

  //return 1.0;
  return ptweights->GetBinContent(ptweights->FindFixBin(genhpt));
}


// set to hold all the per-photon efficiency scale factors
// only PHEFFSCALESET s that are turned ON will be considered
// map holds the name of the set, than a pointer to another map, holding pairs of <etamin,etamax> and a vector with entries: #R9 boundaries, (R9boundaries), (scalefactors)
// so length of vector must be (1 + #R9boundaries + #R9boundaries+1) = 4 in the case of one R9 boundary
std::map<TString, std::map< std::pair<double,double>, std::vector<double>* >*> phEffScaleSets;

float effweightCard( float sceta, float r9 ) {

  double scale = 1.;
  
  // loop over all scale factor sets
  for( std::map<TString, std::map< std::pair<double,double>, std::vector<double>* >*>::iterator it = phEffScaleSets.begin(); it != phEffScaleSets.end(); ++it ){
    std::map< std::pair<double, double>, std::vector<double>* >* thisSetMap = it->second;

    for( std::map< std::pair<double, double>, std::vector<double>* >::iterator itEta = thisSetMap->begin(); itEta != thisSetMap->end(); ++itEta ) {

      std::pair<double, double> etaBound = itEta->first;
      std::vector<double>*      values   = itEta->second;

      if ( sceta > (float) etaBound.first && sceta < (float) etaBound.second ) { // photon is in this eta-bin
	if ( r9 < (float) values->at(0) ) scale *= values->at(1);
	else                              scale *= values->at(2);
      }      
    }
  }
  
  return scale;
}

float trigeffscale() {  

  return (float) 0.9968;

}

bool ApplySysAsFriend(TTree* intree, TTree* subtree, std::map<TString,TString>* nuisMap, TString name, TH2D** histos_SYS = NULL) {
  

  std::vector<TTreeFormula*> inputforms_tot;
  inputforms_tot.push_back(new TTreeFormula("hpt","hpt",intree));
  inputforms_tot.push_back(new TTreeFormula("hy","hy",intree));

  std::vector<TTreeFormula*> inputforms_sub;
  inputforms_sub.push_back(new TTreeFormula("genHiggspt","genHiggspt",subtree));

  TTree *friendtree_tot = new TTree();
  TTree *friendtree_sub = new TTree();
  
  friendtree_tot->SetName(TString::Format("mvasystree_tot_%s",name.Data()).Data());
  friendtree_sub->SetName(TString::Format("mvasystree_sub_%s",name.Data()).Data());
  

  std::cout<<" adding friends for "<<nuisMap->size()<<" nuissances."<<std::endl;

  std::vector<float> _nuisUp   (nuisMap->size());
  std::vector<float> _nuisDo   (nuisMap->size());

  int nuisCounter=0;
  for(std::map<TString,TString>::iterator nIt = nuisMap->begin(); nIt != nuisMap->end(); ++nIt ) {

    TString   upName = TString::Format("%s_up",(nIt->first).Data());
    TString   doName = TString::Format("%s_do",(nIt->first).Data());
  
    friendtree_tot->Branch(upName.Data(),  &(_nuisUp[nuisCounter]),TString::Format("%s/F",upName.Data()).Data());
    friendtree_tot->Branch(doName.Data(),  &(_nuisDo[nuisCounter]),TString::Format("%s/F",doName.Data()).Data());

    friendtree_sub->Branch(upName.Data(),  &(_nuisUp[nuisCounter]),TString::Format("%s/F",upName.Data()).Data());
    friendtree_sub->Branch(doName.Data(),  &(_nuisDo[nuisCounter]),TString::Format("%s/F",doName.Data()).Data());
    nuisCounter++;
  }
  
  Float_t  _tot_pt = 0.;
  Float_t  _tot_y  = 0.;
  Float_t  _genPt  = 0.;

  int currenttree = -1;
  Long64_t subTreeEntry = 0;
  subtree->LoadTree(subTreeEntry);
  int thistree = subtree->GetTreeNumber();
  bool newtree = currenttree!=thistree;
  currenttree = thistree;
  if (newtree) inputforms_sub[0]->Notify();
  _genPt     = inputforms_sub[0]->EvalInstance();
  
  for (Long64_t iev=0; iev<intree->GetEntries(); ++iev) {
    if (iev%100000==0) printf("%i\n",int(iev));
    intree->LoadTree(iev);
    thistree = intree->GetTreeNumber();
    newtree  = currenttree!=thistree;
    currenttree = thistree;
    
    if (newtree) {
      inputforms_tot[0]->Notify();
      inputforms_tot[1]->Notify();
    }
    _tot_pt = inputforms_tot[0]->EvalInstance();
    _tot_y  = inputforms_tot[1]->EvalInstance();
    
    double nomVal = (histos_SYS ? ( _tot_pt < 300. ? histos_SYS[0]->Interpolate(TMath::Abs(_tot_y),_tot_pt) : 1.) : 1.);

    nuisCounter=0;
    for(std::map<TString,TString>::iterator nIt = nuisMap->begin(); nIt != nuisMap->end(); ++nIt ) {
      bool isGFscale = !( (nIt->first).CompareTo("CMS_hgg_n_sc_gf") );
      
      float _nuis_up        = ( histos_SYS ? ( _tot_pt < 300. ? histos_SYS[2*nuisCounter+1]->Interpolate(TMath::Abs(_tot_y),_tot_pt) / nomVal : 1. ) : 1.);
      if( isGFscale && histos_SYS ) _nuis_up /= 10.;
      float _nuis_do        = ( histos_SYS ? ( _tot_pt < 300. ? histos_SYS[2*nuisCounter+2]->Interpolate(TMath::Abs(_tot_y),_tot_pt) / nomVal : 1. ) : 1.);
      
      _nuisUp[nuisCounter] = _nuis_up;
      _nuisDo[nuisCounter] = _nuis_do;
      
      nuisCounter++;
    }

    friendtree_tot->Fill();

    if( _genPt == _tot_pt ) {
      friendtree_sub->Fill();      
      subTreeEntry++;
      subtree->LoadTree(subTreeEntry);
      thistree = subtree->GetTreeNumber();
      newtree = currenttree!=thistree;
      currenttree = thistree;
      if (newtree) inputforms_sub[0]->Notify();
      _genPt     = inputforms_sub[0]->EvalInstance();
    }

  }
  
  intree ->AddFriend(friendtree_tot);
  subtree->AddFriend(friendtree_sub);

  return (subTreeEntry == subtree->GetEntries() );

}


TTree *ApplyAsFriend(TTree *intree, TString tmvaweights, const std::vector<std::string> &vars, std::string targetname, bool scaleIDMVA = false, bool scaleSigEoE = false) {
  
  int nvars = vars.size();
    
  //initialize TTreeFormulas to read variables from TTree
  std::vector<TTreeFormula*> inputforms;
  for (std::vector<std::string>::const_iterator it = vars.begin(); 
      it != vars.end(); ++it) {
    inputforms.push_back(new TTreeFormula(it->c_str(),it->c_str(),intree));
  }
  
  Float_t target = 0.;
  Float_t *vals = new Float_t[nvars+12];
  
  inputforms.push_back(new TTreeFormula("ph1.sceta","ph1.sceta",intree));
  inputforms.push_back(new TTreeFormula("ph2.sceta","ph2.sceta",intree));

  inputforms.push_back(new TTreeFormula("ph1.r9","ph1.r9",intree));
  inputforms.push_back(new TTreeFormula("ph2.r9","ph2.r9",intree));

  inputforms.push_back(new TTreeFormula("ph1.e","ph1.e",intree));
  inputforms.push_back(new TTreeFormula("ph1.eerr","ph1.eerr",intree));
  inputforms.push_back(new TTreeFormula("ph1.eerrsmeared","ph1.eerrsmeared",intree));

  inputforms.push_back(new TTreeFormula("ph2.e","ph2.e",intree));
  inputforms.push_back(new TTreeFormula("ph2.eerr","ph2.eerr",intree));
  inputforms.push_back(new TTreeFormula("ph2.eerrsmeared","ph2.eerrsmeared",intree));

  inputforms.push_back(new TTreeFormula("mass","mass",intree));
  inputforms.push_back(new TTreeFormula("deltamvtx","deltamvtx",intree));


  //initialize tmva reader
  TMVA::Reader* tmva = new TMVA::Reader();
  for (unsigned int ivar=0; ivar<vars.size(); ++ivar) {
    tmva->AddVariable(vars.at(ivar),&vals[ivar]);
  }  
  tmva->BookMVA("BDTG",tmvaweights);  
  
  //initialize new friend tree
  TTree *friendtree = new TTree();
  friendtree->SetName(targetname.c_str());
  friendtree->Branch(targetname.c_str(),&target,TString::Format("%s/F",targetname.c_str()));
  
  int currenttree = -1;
  for (Long64_t iev=0; iev<intree->GetEntries(); ++iev) {
    if (iev%100000==0) printf("%i\n",int(iev));
    intree->LoadTree(iev);
    int thistree = intree->GetTreeNumber();
    bool newtree = currenttree!=thistree;
    currenttree = thistree;
    
    for (int i=0; i<nvars+12; ++i) {
      if (newtree) inputforms[i]->Notify();
      vals[i] = inputforms[i]->EvalInstance();
    }
    
    if (scaleIDMVA) {
      bool ph1_isHR9 = ( vals[nvars] > 0.94 );
      bool ph2_isHR9 = ( vals[nvars+1] > 0.94 );
      bool ph1_isEB  = (TMath::Abs(vals[nvars+2]) < 1.5);
      bool ph2_isEB  = (TMath::Abs(vals[nvars+3]) < 1.5);
      
      if(ph1_isEB) {
	if(ph1_isHR9)
	  vals[nvars-2] = scaleIDMVA_EB_hR9->Eval( (double) vals[nvars-2] );
	else
	  vals[nvars-2] = scaleIDMVA_EB_lR9->Eval( (double) vals[nvars-2] );
      } else {
	if(ph1_isHR9)
	  vals[nvars-2] = scaleIDMVA_EE_hR9->Eval( (double) vals[nvars-2] );
	else
	  vals[nvars-2] = scaleIDMVA_EE_lR9->Eval( (double) vals[nvars-2] );
      }
      
      if(ph2_isEB) {
	if(ph2_isHR9)
	  vals[nvars-1] = scaleIDMVA_EB_hR9->Eval( (double) vals[nvars-1] );
	else
	  vals[nvars-1] = scaleIDMVA_EB_lR9->Eval( (double) vals[nvars-1] );
      } else {
	if(ph2_isHR9)
	  vals[nvars-1] = scaleIDMVA_EE_hR9->Eval( (double) vals[nvars-1] );
	else
	  vals[nvars-1] = scaleIDMVA_EE_lR9->Eval( (double) vals[nvars-1] );
      }
      
    }
    
    if (scaleSigEoE) {
      
      double _mass            = vals[nvars+10];
      double _deltamvtx       = vals[nvars+11];
      
      double _ph1_e           = vals[nvars+4];
      double _ph1_eerr        = vals[nvars+5];
      double _ph1_eerrsmaered = vals[nvars+6];
      
      double _ph2_e           = vals[nvars+7];
      double _ph2_eerr        = vals[nvars+8];
      double _ph2_eerrsmaered = vals[nvars+9];
      
      double _ph1_smear = TMath::Power(_ph1_eerrsmaered,2) - TMath::Power(_ph1_eerr,2);
      double _ph2_smear = TMath::Power(_ph2_eerrsmaered,2) - TMath::Power(_ph2_eerr,2);
      
      double _ph1_sigEoE = 0.; //TMath::Exp( TMath::Log(_ph1_eerr/_ph1_e) );
      double _ph2_sigEoE = 0.; //TMath::Exp( TMath::Log(_ph2_eerr/_ph2_e) );
      
      if( TMath::Abs(vals[nvars]) < 1.0 ) {
	_ph1_sigEoE = TMath::Exp( scaleSigEoE_0_10 ->Eval( TMath::Log(_ph1_eerr/_ph1_e) ) );	
      } else if(  TMath::Abs(vals[nvars]) < 1.5 ) {
	_ph1_sigEoE = TMath::Exp( scaleSigEoE_10_15->Eval( TMath::Log(_ph1_eerr/_ph1_e) ) );	
      } else if(  TMath::Abs(vals[nvars]) < 2.0 ) {
	_ph1_sigEoE = TMath::Exp( scaleSigEoE_15_20->Eval( TMath::Log(_ph1_eerr/_ph1_e) ) );	
      } else {
	_ph1_sigEoE = TMath::Exp( scaleSigEoE_20_25->Eval( TMath::Log(_ph1_eerr/_ph1_e) ) );	
      }
      
      if( TMath::Abs(vals[nvars+1]) < 1.0 ) {
	_ph2_sigEoE = TMath::Exp( scaleSigEoE_0_10 ->Eval( TMath::Log(_ph2_eerr/_ph2_e) ) );	
      } else if(  TMath::Abs(vals[nvars+1]) < 1.5 ) {
	_ph2_sigEoE = TMath::Exp( scaleSigEoE_10_15->Eval( TMath::Log(_ph2_eerr/_ph2_e) ) );	
      } else if(  TMath::Abs(vals[nvars+1]) < 2.0 ) {
	_ph2_sigEoE = TMath::Exp( scaleSigEoE_15_20->Eval( TMath::Log(_ph2_eerr/_ph2_e) ) );	
      } else {
	_ph2_sigEoE = TMath::Exp( scaleSigEoE_20_25->Eval( TMath::Log(_ph2_eerr/_ph2_e) ) );	
      }
      
      double _ph1_eerr_mod = _ph1_sigEoE * _ph1_e;
      double _ph2_eerr_mod = _ph2_sigEoE * _ph2_e;
      
      double _ph1_eerrsmeared_mod = TMath::Sqrt( TMath::Power(_ph1_eerr_mod,2) + _ph1_smear );
      double _ph2_eerrsmeared_mod = TMath::Sqrt( TMath::Power(_ph2_eerr_mod,2) + _ph2_smear );
      
      double _ph1_sigEoE_mod = _ph1_eerrsmeared_mod/_ph1_e;
      double _ph2_sigEoE_mod = _ph2_eerrsmeared_mod/_ph2_e;
      
      double _masserrsmeared         = 0.5*_mass*TMath::Sqrt( TMath::Power( _ph1_sigEoE_mod, 2) + TMath::Power( _ph2_sigEoE_mod, 2) );
      double _masserrsmearedwrongvtx = TMath::Sqrt(_masserrsmeared*_masserrsmeared + _deltamvtx*_deltamvtx);
      
      vals[0] = _masserrsmeared/_mass;
      vals[1] = _masserrsmearedwrongvtx/_mass;
    } 

    target = tmva->EvaluateMVA("BDTG");    
    friendtree->Fill();    
  }
  
  //clear TMVA reader
  delete tmva;
  
  //clear TTreeFormulas
  for (std::vector<TTreeFormula*>::const_iterator it = inputforms.begin(); 
        it != inputforms.end(); ++it) {
      delete *it;
  }
  
  delete[] vals;
  
  intree->AddFriend(friendtree);
  return friendtree;
  
}

//---------------------------------------------------------------------------------------------------------------------------
// Fwd helper function declaration... implementation at the end
bool readWeightsFromConfigCard(TString fileName, 
			       TH1D*& puhisto, TFile*& ptweight, TString& vtxWeightFile);

bool readParmsAndCatsFromConfigCard( TString fileName,
				     bool& computeMVAvar,
				     bool& correctIDMVAvar,
				     bool& correctSigEoEvar,
				     TString& mvaWeightFile,
				     TString& mvaDefFile,
				     TString& idmvaCorrFile,
				     TString& sigeoeCorrFile,
				     TString& projectDir,
				     TString& modname,
				     TString& treename,
				     double& massmax, double& massmin,
				     double&theLumi,
				     std::map<TString,TString>& auxCats,
				     std::map<TString,TCut>& anaCats,
				     TCut& baseCut
				     );



bool readConfigCard(TString configCardName, 
		    std::map<TString,TString>& procMCfileMap,
		    std::map<TString,bool>&    procOn,
		    std::map<int,TString>&     procIdxMap,			
		    int& numProcs, int& numCats,
		    std::vector<TString>& catNames,
		    std::vector<TString>& catDesc);

RooDataSet *makedset(TString name, TTree *tree, TCut cut, RooRealVar *hmass, TString xVarName, RooRealVar *weight) {
 
  RooDataSet *dset = new RooDataSet(name,"",RooArgSet(*hmass,*weight),weight->GetName());
  RooRealVar *vset = (RooRealVar*)dset->get()->find(hmass->GetName());
  
//   std::cout<<" -----------------------------------------"<<std::endl;
//   std::cout<<" Creating dataset with name = "<<name<<" :"<<std::endl;
//   std::cout<<"   "<<cut<<std::endl;
//   std::cout<<" -----------------------------------------"<<std::endl;

  tree->SetEstimate(tree->GetEntries());

  TH1D* dummy = new TH1D(TString::Format("dummy_%s",name.Data()),"",10000,0.,10000.);

  Int_t nev = tree->Draw(TString::Format("%s>>%s",xVarName.Data(),dummy->GetName()).Data(),cut,"goff");
  
  //   float target = 0.;
  //   tree->SetBranchAddress(xVarName.Data(), &target);

  //  std::cout<<"  got estimate for #events = "<<nev<<std::endl;

//   TString massrecomp = " TMath::Sqrt(ph1.e*ph2.e*(1-gencostheta)*2.) ";
//   Int_t nev = tree->Draw(massrecomp.Data(),cut,"goff");

  double *vals    = tree->GetV1();
  double *weights = tree->GetW();
  
  for (int iev=0; iev<nev; ++iev) {
    vset->setVal(vals[iev]);
    dset->add(*vset,weights[iev]);
  }
  
  //  std::cout<<"  ds creation done "<<std::endl;


  // remove the dummy from above...
  delete dummy;
  return dset;
  
  
}

//---------------------------------------------------------------------------------------------------------------------------
#define DOBDT


//void createSignalModels(TString configCardName="mettag.config", modeltype model = STANDARDMODEL) {
void fitmvasysCard_HGG(
		       TString configCardName="/home/fabstoec/cms/root/templateHGG_8TeV_Moriond.config") {

  //************************************************************************************
  // SETUP PART
  // -----------------------------------------------------
  std::map<TString,double*> processCrossSectionMap;
  initProcessCSArrayMap(processCrossSectionMap);

  TH1::AddDirectory(false);


  // generate a map for the guys, so we can acces by names
  std::map<TString,double*>* modelParmsMap = new std::map<TString,double*>();
  modelParmsMap->insert(std::pair<TString,double*>(TString("smbr"),smbr));
  modelParmsMap->insert(std::pair<TString,double*>(TString("ffbr") ,ffbr));
  modelParmsMap->insert(std::pair<TString,double*>(TString("sm4br"),sm4br));
  modelParmsMap->insert(std::pair<TString,double*>(TString("gghxsec"),gghxsec));
  modelParmsMap->insert(std::pair<TString,double*>(TString("vbfxsec"),vbfxsec));
  modelParmsMap->insert(std::pair<TString,double*>(TString("wzhxsec"),wzhxsec));

  modelParmsMap->insert(std::pair<TString,double*>(TString("tthxsec"),tthxsec));
  modelParmsMap->insert(std::pair<TString,double*>(TString("sm4gghxsec"),sm4gghxsec));

  //----------------------------------------------------------------------------
  // sum the cross-sections properly
  std::vector<TString>                      sumRuleNames;
  std::vector<std::vector<const char*>* >   sumRuleComps;

  std::vector<const char*> myVec;
  
  std::vector<const char*> sumNames_SM;
  sumNames_SM.push_back("wzhxsec");
  sumNames_SM.push_back("gghxsec");
  sumNames_SM.push_back("vbfxsec");
  sumNames_SM.push_back("tthxsec");
  sumRuleNames.push_back("smxsec");
  sumRuleComps.push_back(&sumNames_SM);
  
  std::vector<const char*> sumNames_FF;
  sumNames_FF.push_back("wzhxsec");
  sumNames_FF.push_back("vbfxsec");
  sumRuleNames.push_back("ffxsec");
  sumRuleComps.push_back(&sumNames_FF);
  
  std::vector<const char*> sumNames_SM4;
  sumNames_SM4.push_back("wzhxsec");
  sumNames_SM4.push_back("sm4gghxsec");
  sumNames_SM4.push_back("vbfxsec");
  sumNames_SM4.push_back("tthxsec");
  sumRuleNames.push_back("sm4xsec");
  sumRuleComps.push_back(&sumNames_SM4);
  //----------------------------------------------------------------------------

  std::map<TString,std::vector<const char*>*> procxseclist_list;

  std::vector<const char*> smlist;
  smlist.push_back("gghxsec");
  smlist.push_back("vbfxsec");
  smlist.push_back("wzhxsec");
  smlist.push_back("tthxsec");
  procxseclist_list.insert(std::pair<TString,std::vector<const char*>*>("sm",&smlist));

  std::vector<const char*> sm4list;
  sm4list.push_back("sm4gghxsec");
  sm4list.push_back("vbfxsec");
  sm4list.push_back("wzhxsec");
  sm4list.push_back("tthxsec");
  procxseclist_list.insert(std::pair<TString,std::vector<const char*>*>("sm4",&sm4list));

  std::vector<const char*> fflist;
  fflist.push_back("vbfxsec");
  fflist.push_back("wzhxsec");
  procxseclist_list.insert(std::pair<TString,std::vector<const char*>*>("ff",&fflist));

  
  //************************************************************************************

  // number of processes, Categories, Models and mass points
  int numProcs  = -1;
  int numCats   = -1;

  // process names and toggle process ON/OFF
  //std::vector<TString> procnames;
  //std::vector<bool>    procon;

  std::map<int,TString>     procIdxMap;
  std::map<TString,TString> procMCfileMap;
  std::map<TString,bool>    procon;
  
  double massmax = -1.;
  double massmin = -1.;
  
  // auxiliary vectro for base cat names...
  std::vector<TString> catnamesbase;
  std::vector<TString> catdesc;
  
  TString    projectDir;
  TString    modname;
  TString    treename;
  
  std::map<TString,TString> auxCatMap;
  std::map<TString,TCut>    anaCatMap;  
  TCut theBaseCut;
  
  bool    computeMVAvar;
  bool    correctIDMVAvar;
  bool    correctSigEoEvar;

  TString mvaWeightFile;
  TString mvaDefFile;

  TString idmvaCorrFile;
  TString sigeoeCorrFile;

  bool status = readConfigCard(configCardName, 
			       procMCfileMap,
			       procon,
			       procIdxMap,
			       numProcs, numCats,
			       catnamesbase,
			       catdesc);


  if(!status && true ) {
    std::cerr<<" ERROR when readin input card "<<configCardName<<"."<<std::endl;
    return;
  }

  // 1. Set up the PU and pt-weights from input files
  TH1D*  hpuest       = NULL;
  TFile* fileptweight = NULL;
  TString vtxWeightFileName = "";

  status = readWeightsFromConfigCard( configCardName       , 
				      hpuest              ,
				      fileptweight        ,
				      vtxWeightFileName   );
  
  if( !status ) {
    std::cerr<<" ERROR: Could not read weight information from file "<<configCardName.Data()<<std::endl;
    return;
  }

  if( ! fileptweight )
    std::cerr<<" [INFO] Higgs pt-rewighting turned OFF."<<std::endl;

  if( !hpuest )
    std::cerr<<" [INFO] PU-rewighting switched OFF."<<std::endl;

  // --------------------------------------------------------------------------------------------------
  // set up the Vtx weights
  
  if ( vtxWeightFileName.CompareTo("") ) {
    //TFile* vtxWeightFile = TFile::Open("/home/fabstoec/cms/cmssw/029/CMSSW_5_3_2_patch4/src/UserCode/HiggsAnalysis/HiggsTo2photons/h2gglobe/Macros/vertex_reweighing_mva_HCP2012_unblind.root");
    TFile* vtxWeightFile = TFile::Open(vtxWeightFileName.Data());
    TString graphName_pass = "ratioVertex_cat0_pass";
    TString graphName_fail = "ratioVertex_cat0_fail";
    
    vtxWeights_pass = (TGraphAsymmErrors*) vtxWeightFile->Get(graphName_pass.Data());
    vtxWeights_fail = (TGraphAsymmErrors*) vtxWeightFile->Get(graphName_fail.Data());
    
    if (!vtxWeights_pass || !vtxWeights_fail) {
      std::cerr<<" WARNING: Could not load vtx weights."<<std::endl;
      return;
    }
  } else {
    vtxWeights_pass = NULL;
    vtxWeights_fail = NULL;
    std::cout<<" [INFO] Vtx Reweigthing turned OFF."<<std::endl;
  }

  // ---------------------------------------------------------------------------------------------------------

  double lumi = 1.;

  status = readParmsAndCatsFromConfigCard( configCardName      ,
					   computeMVAvar       ,
					   correctIDMVAvar     ,
					   correctSigEoEvar    ,
					   mvaWeightFile       ,
					   mvaDefFile          ,
					   idmvaCorrFile       ,
					   sigeoeCorrFile      ,
					   projectDir          ,
					   modname             ,
					   treename            ,
					   massmax, massmin    ,
					   lumi,
					   auxCatMap           ,
					   anaCatMap           ,
					   theBaseCut
					   );
  
  
  if( !status ) {
    std::cerr<<" ERROR: Could not read Category information from file "<<configCardName.Data()<<std::endl;
    return;
  }

  // ------------------------------------------------------------------------------------------
  // set-up IDMVA scale functions, if required
  TFile* idmvaScaleFile = NULL;
  if( correctIDMVAvar ) {
    idmvaScaleFile = TFile::Open(idmvaCorrFile.Data());
    scaleIDMVA_EB_hR9 = (TGraph*) idmvaScaleFile->Get( "idmvascale_ABCD_hR9_EB" );
    scaleIDMVA_EB_lR9 = (TGraph*) idmvaScaleFile->Get( "idmvascale_ABCD_lR9_EB" );
    scaleIDMVA_EE_hR9 = (TGraph*) idmvaScaleFile->Get( "idmvascale_ABCD_hR9_EE" );
    scaleIDMVA_EE_lR9 = (TGraph*) idmvaScaleFile->Get( "idmvascale_ABCD_lR9_EE" );
  }

  TFile* sigEoEscaleFile = NULL;
  if( correctSigEoEvar ) {
    //sigEoEscaleFile = TFile::Open("/home/fabstoec/cms/root/PhotonIDMVA_new/Moriond13_phSigEoE.root");
    sigEoEscaleFile = TFile::Open(sigeoeCorrFile.Data());

    scaleSigEoE_0_10  = (TGraph*) sigEoEscaleFile->Get( "sigeoescale_ABCD_0_10"  );
    scaleSigEoE_10_15 = (TGraph*) sigEoEscaleFile->Get( "sigeoescale_ABCD_10_15" );
    scaleSigEoE_15_20 = (TGraph*) sigEoEscaleFile->Get( "sigeoescale_ABCD_15_20" );
    scaleSigEoE_20_25 = (TGraph*) sigEoEscaleFile->Get( "sigeoescale_ABCD_20_25" );
  }


  // ==================================================================================

  // everything is set up.... starting real work...

  // ------------------------------------------------------------------
  // General ROOT setup
  gROOT->Macro("MitStyle.C");
  gStyle->SetErrorX(0); 
  gStyle->SetOptStat(0);
  gROOT->ForceStyle();  
  gStyle->SetOptStat(1110);

  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(0);

  gStyle->SetLabelColor(1, "XYZ");
  gStyle->SetLabelFont(42, "XYZ");
  gStyle->SetLabelOffset(0.007, "XYZ");
  gStyle->SetLabelSize(0.035, "XYZ");

  gStyle->SetTitleColor(1, "XYZ");
  gStyle->SetTitleFont(42, "XYZ");
  gStyle->SetTitleOffset(1., "XZ");
  gStyle->SetTitleOffset(1.5, "Y");
  gStyle->SetTitleSize(0.04, "XYZ");

  RooMsgService::instance().getStream(1).removeTopic(RooFit::Caching);
  RooMsgService::instance().getStream(0).removeTopic(RooFit::Minimization);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::Minimization);
  RooMsgService::instance().getStream(0).removeTopic(RooFit::Plotting);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::Plotting);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::Fitting);
  RooMsgService::instance().getStream(0).removeTopic(RooFit::Eval);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::Eval);
  // ------------------------------------------------------------------

  RooRealVar*   hmass = new RooRealVar("CMS_hgg_mass","m_{#gamma#gamma}",massmin,massmax,"GeV");
  hmass->setRange(massmin,massmax);
  hmass->setBins( (int) (1.0*(massmax-massmin)) );
  hmass->SetTitle("m_{#gamma#gamma}");
  hmass->setUnit("GeV");  

  RooRealVar*   hpt = new RooRealVar("CMS_hpt","p_{T}",0.,10000.,"GeV");

  RooRealVar mnom("MH","m_{h}",120.,massmin,massmax,"GeV");
  mnom.setConstant();
  
  //---RooRealVar to load dataset or used for calculation---
  RooRealVar *weight = new RooRealVar("weight","",1.0);//this seems to be useless variable
  weight->removeRange();

  // open the MVA files, if requested
  TFile* tmvaOutput = NULL;
  TString weights   = "";
  std::vector<std::string> *varlist = NULL;

  if( computeMVAvar ) {
    tmvaOutput = TFile::Open(mvaDefFile.Data(),"READ");//root file to store training information  
    varlist = (std::vector<std::string>*)tmvaOutput->Get("varlist");
    weights = mvaWeightFile;
  }  
  
  // loop over all processes that are on
  TFile *friendtmp = TFile::Open("friendtmp.root","RECREATE");
  

  // do a prcname<->MVAsys map
  std::map<TString,TString> procToSysMap;
  //procToSysMap.insert(std::pair<TString,TString>("ggh","theorySystematics/GF_mvasys.root"));
  //procToSysMap.insert(std::pair<TString,TString>("vbf","theorySystematics/VBF_mvasys.root"));
  procToSysMap.insert(std::pair<TString,TString>("ggh","GF"));
  procToSysMap.insert(std::pair<TString,TString>("vbf","VBF"));

  // add nuissances
  std::map<TString,TString> nuissancesGF;
  std::map<TString,TString> nuissancesVBF;

  std::vector<TString> allNuisances;
  allNuisances.push_back("CMS_hgg_n_sc_gf");
  allNuisances.push_back("CMS_hgg_n_sc_vbf");

  nuissancesGF.insert(std::pair<TString,TString>("CMS_hgg_n_sc_gf","up:down"));
  nuissancesVBF.insert(std::pair<TString,TString>("CMS_hgg_n_sc_vbf","up:down"));
  for( int iPdf=1;iPdf <=26; ++iPdf ) {
    nuissancesGF.insert(std::pair<TString,TString>(TString::Format("CMS_hgg_n_pdf_%d",iPdf),TString::Format("PDF_%d:PDF_%d",2*(iPdf-1)+1,2*(iPdf-1)+2)));
    nuissancesVBF.insert(std::pair<TString,TString>(TString::Format("CMS_hgg_n_pdf_%d",iPdf),TString::Format("PDF_%d:PDF_%d",2*(iPdf-1)+1,2*(iPdf-1)+2)));

    allNuisances.push_back( TString::Format("CMS_hgg_n_pdf_%d",iPdf)  );

  }
  
  std::map<TString,std::map<TString,TString>*> procToNuisMap;
  procToNuisMap.insert(std::pair<TString,std::map<TString,TString>*>("ggh",&nuissancesGF));
  procToNuisMap.insert(std::pair<TString,std::map<TString,TString>*>("vbf",&nuissancesVBF));


  std::map<TString, std::map<TString, std::map<TString,TString>*>*> allNuisMap;
  for( UInt_t iNuis = 0; iNuis < allNuisances.size(); ++iNuis ) {
    std::map<TString, std::map<TString, std::map<TString,TString>*>*>::iterator it = allNuisMap.find(allNuisances[iNuis]);
    if ( it == allNuisMap.end() ) {
      allNuisMap.insert(std::pair<TString, std::map<TString, std::map<TString,TString>*>*>( allNuisances[iNuis], new std::map<TString, std::map<TString,TString>*>() ) );
      it = allNuisMap.find(allNuisances[iNuis]);
    } 
    for (UInt_t iCat=0; iCat<catnamesbase.size(); ++iCat)
      it->second->insert(std::pair<TString, std::map<TString,TString>*>(catnamesbase[iCat], new std::map<TString,TString>()));
  }


  TH2D*** histos_SYS = new TH2D**[numProcs];

  TFile** theFile    = new TFile*     [numProcs];
  TDirectory** hdir  = new TDirectory*[numProcs];
  TTree** theTree    = new TTree*     [numProcs];
  TTree** tallevts   = new TTree*     [numProcs];

  TFile** sysFiles   = new TFile*     [numProcs];

  // full Cat names including the process name

  // keep alll the datasets... really?
  RooDataSet*** datasets = new RooDataSet**[numProcs];
  
  TCut theWeight       = "puweight(numPU)*bsweight(vtxZ,genHiggsZ)*trigeffscale()*effweightCard(ph1.sceta,ph1.r9)*effweightCard(ph2.sceta,ph2.r9)*vtxWeight(ptgg,vtxZ,genHiggsZ)";
  TCut masscutReal     =  theBaseCut && TCut(TString::Format("mass > %f && mass < %f",massmin,massmax).Data());    

  for(int iProc = 0; iProc < numProcs; ++iProc) {
    
    // get the proc name from the map
    TString procname = procIdxMap.find(iProc)->second;
    if( !(procon.find(procname)->second) ) continue;
    
    std::vector<TString> catnames;
    for(int iCat = 0; iCat < numCats; ++iCat)
      catnames.push_back(TString::Format("%s_%s",catnamesbase[iCat].Data(),procname.Data()));
    
    std::cout<<"  #Categories = "<<numCats<<std::endl;
    
    
    // load the systematics ROOT file
    std::map<TString,TString>::iterator fIt = procToSysMap.find(procname);
    if ( fIt == procToSysMap.end() ) {
      std::cout<<" WARNING: No sysmetaics file for process "<<procname<<"."<<std::endl;
      std::cout<<"          Skippoing process and assign -/- for all systematics."<<std::endl;
      continue;
    }
    
    // load all the tress for this process...
    double masspoint = 120.;  // THE MASS...!
      
    // load the correct tree from the correct file...
    std::map<TString,TString>::iterator tempIt = procMCfileMap.find(procname);
    if( tempIt == procMCfileMap.end() ) {
      std::cout<< " ERROR *** "<<std::endl;
      return;
    }
    
   TString samplestring = tempIt->second;      

    theFile[iProc] = TFile::Open(TString::Format(samplestring, (int) masspoint),"READ");
    if (!theFile[iProc]) {
      std::cerr<<" ERROR: Could not open file "<<TString::Format(samplestring, (int) masspoint).Data()<<"."<<std::endl;
      return;
    }
    
    // load the Tree
    hdir    [iProc] = (TDirectory*) theFile[iProc]->FindObjectAny(modname);
    theTree [iProc] = (TTree*)      hdir   [iProc]->Get(treename.Data());
    tallevts[iProc] = (TTree*)      theFile[iProc]->FindObjectAny("hMVAtuple");

    if( computeMVAvar )
      ApplyAsFriend(theTree[iProc],weights,*varlist,"bdt",correctIDMVAvar,correctSigEoEvar);
    //ApplyAsFriend(theTree[iProc],weights,*varlist,"bdt");

      
    sysFiles[iProc] = TFile::Open(TString::Format("theorySystematics/%s_mvasys.root",fIt->second.Data()),"READ");
    if ( !sysFiles[iProc] ) {
      std::cerr<<" ERROR: Could not open systematics file "<<TString::Format("theorySystematics/%s_mvasys.root",fIt->second.Data())<<"."<<std::endl;
      return;
    }

    // loop over all the uissances for this process
    std::map<TString,std::map<TString,TString>*>::iterator nmIt =  procToNuisMap.find(procname);
    if ( nmIt == procToNuisMap.end() ) {
      std::cerr<<" ERROR: COuld not find nuissance map for process "<<procname<<"."<<std::endl;
      return;
    }
    
    std::map<TString,TString>* nuisMap = nmIt->second;
    if (!nuisMap ){
      std::cerr<<" ERROR: Could not find nuissance list for process "<<procname<<"."<<std::endl;
      return;
    }

    // load all the histograms needed
    int numNuis = nuisMap->size();
    histos_SYS[iProc] = new TH2D*[2*numNuis+1];    
    std::cout<<" Nuissance list successfully loaded for process "<<procname<<"."<<std::endl;
    std::cout<<"         Name of histogram fro central scales : "<<TString::Format("%s_cent",(fIt->second).Data())<<std::endl;
    histos_SYS[iProc][0] = (TH2D*) ((TH2D*) sysFiles[iProc]->FindObjectAny(TString::Format("%s_cent",(fIt->second).Data())))->Clone();

    int nuisCounter=0;
    for(std::map<TString,TString>::iterator nIt = nuisMap->begin(); nIt != nuisMap->end(); ++nIt ) {
      Ssiz_t split =  (nIt->second).First(':');      
      TString upVar = (nIt->second)(0,split);
      TString doVar = (nIt->second)(split+1,(nIt->second).Length()-split-1);
      histos_SYS[iProc][2*nuisCounter+1] = (TH2D*) ((TH2D*)sysFiles[iProc]->FindObjectAny(TString::Format("%s_%s",(fIt->second).Data(),upVar.Data())))->Clone(); 
      histos_SYS[iProc][2*nuisCounter+2] = (TH2D*) ((TH2D*)sysFiles[iProc]->FindObjectAny(TString::Format("%s_%s",(fIt->second).Data(),doVar.Data())))->Clone(); 


//       histos_SYS[2*nuisCounter+1]->Draw();
//       histos_SYS[2*nuisCounter+2]->Draw();

      nuisCounter++;

    }

    bool sysFine = ApplySysAsFriend(tallevts[iProc],theTree[iProc],nuisMap,procname,histos_SYS[iProc]);    
    //bool sysFine = ApplySysAsFriend(tallevts[iProc],theTree[iProc],nuisMap,procname,NULL);

    if( !sysFine ) {
      std::cerr<<" ERROR: Systematics could not be assigned."<<std::endl;
      return;
    }

    std::cout<<"  Systematic weight assiged as friend."<<std::endl;


    // set the PU weight correctly...
    setpuweights(theFile[iProc], hpuest);


    datasets[iProc] = new RooDataSet*[catnames.size()*(nuisMap->size()*4 + 2 )];

    // =========== EFFECTIVE SMAERING ==================
    // keep all datasets, also the un-smeared ones.

    int dsCounter=0;

    for (UInt_t iCat=0; iCat<catnames.size(); ++iCat) {      
      std::cout<<"  TESTING CAT: "<<catnamesbase.at(iCat)<<std::endl;

      // first get the nominal eff x acc fro this category
            
      TString dataAllNameTot=TString::Format("sig_%s_mass_m%d_%s_tot",procname.Data(),(int) masspoint,catnamesbase.at(iCat).Data());
      TString dataAllName   =TString::Format("sig_%s_mass_m%d_%s",    procname.Data(),(int) masspoint,catnamesbase.at(iCat).Data());

      std::map<TString,TCut>::iterator catIt = anaCatMap.find(catnamesbase[iCat]);
      if( catIt == anaCatMap.end() ) {
	std::cerr<<" ERROR: Cannot find TCut String for ANACAT with name "<<catnamesbase[iCat]<<"."<<std::endl;
	return;
      }
      

//       RooDataSet *mcsigallwdata           = makedset(dataAllName   ,theTree[iProc], theWeight*(masscutReal  && catIt->second), hmass,"mass", weight);
//       RooDataSet *mcsigallwdataTot        = makedset(dataAllNameTot,tallevts[iProc], "(hpt>0.)",                               hpt,  "hpt" , weight);

      datasets[iProc][dsCounter]        = makedset(dataAllName   ,theTree[iProc], theWeight*(masscutReal  && catIt->second), hmass,"mass", weight);
      dsCounter++;
      datasets[iProc][dsCounter]        = makedset(dataAllNameTot,tallevts[iProc], "(hpt>0.)",                               hpt,  "hpt" , weight);
      dsCounter++;

      double eaccnum = ( datasets[iProc][dsCounter-2]->sumEntries()    );
      double eaccden = ( datasets[iProc][dsCounter-1]->sumEntries() );
      double eacc = eaccnum/eaccden;
      //printf("eacc = %5f, eaccnum = %5f, eaccden = %5f\n",eacc,eaccnum,eaccden);
//       double eaccerrlo = TEfficiency::ClopperPearson(Int_t(eaccden), Int_t(eaccnum), 0.683, kFALSE) - eacc;
//       double eaccerrhi = TEfficiency::ClopperPearson(Int_t(eaccden), Int_t(eaccnum), 0.683, kTRUE) - eacc;


      double nomEff = eacc;

      for(std::map<TString,TString>::iterator nIt = nuisMap->begin(); nIt != nuisMap->end() ; ++nIt ) {

	// firs we do the up variation of this nuissance
	TString weightName = TString::Format("%s_up",(nIt->first).Data());
	
	dataAllNameTot=TString::Format("ds_%s_%s_%s_tot_up",procname.Data(),catnamesbase.at(iCat).Data(),(nIt->first).Data());
	dataAllName   =TString::Format("ds_%s_%s_%s_up",    procname.Data(),catnamesbase.at(iCat).Data(),(nIt->first).Data());

 	datasets[iProc][dsCounter]        = makedset(dataAllName   ,theTree [iProc], TCut(weightName)*theWeight*(masscutReal  && catIt->second), hmass, "mass", weight);
	dsCounter++;
	datasets[iProc][dsCounter]        = makedset(dataAllNameTot,tallevts[iProc], TCut(weightName)*TCut("(hpt>0.)"),                          hpt,    "hpt", weight);
	dsCounter++;

	eaccnum = ( datasets[iProc][dsCounter-2]->sumEntries()    );
	eaccden = ( datasets[iProc][dsCounter-1]->sumEntries() );
	eacc = eaccnum/eaccden;
	//printf("eacc UP = %5f, eaccnum = %5f, eaccden = %5f\n",eacc,eaccnum,eaccden);
// 	double eaccerrlo = TEfficiency::ClopperPearson(Int_t(eaccden), Int_t(eaccnum), 0.683, kFALSE) - eacc;
// 	double eaccerrhi = TEfficiency::ClopperPearson(Int_t(eaccden), Int_t(eaccnum), 0.683, kTRUE) - eacc;


	double upRatio = eacc/nomEff;

	// firs we do the up variation of this nuissance
	weightName = TString::Format("%s_do",(nIt->first).Data());

	dataAllNameTot=TString::Format("ds_%s_%s_%s_tot_do",procname.Data(),catnamesbase.at(iCat).Data(),(nIt->first).Data());
	dataAllName   =TString::Format("ds_%s_%s_%s_do",    procname.Data(),catnamesbase.at(iCat).Data(),(nIt->first).Data());
		
	datasets[iProc][dsCounter]        = makedset(dataAllName   ,theTree [iProc], TCut(weightName)*theWeight*(masscutReal  && catIt->second), hmass, "mass", weight);
	dsCounter++;
	datasets[iProc][dsCounter]        = makedset(dataAllNameTot,tallevts[iProc], TCut(weightName)*TCut("(hpt>0.)"),                          hpt,    "hpt", weight);
	dsCounter++;

	eaccnum = ( datasets[iProc][dsCounter-2]->sumEntries()    );
	eaccden = ( datasets[iProc][dsCounter-1]->sumEntries() );
	eacc = eaccnum/eaccden;
	//printf("eacc DO = %5f, eaccnum = %5f, eaccden = %5f\n",eacc,eaccnum,eaccden);
// 	eaccerrlo = TEfficiency::ClopperPearson(Int_t(eaccden), Int_t(eaccnum), 0.683, kFALSE) - eacc;
// 	eaccerrhi = TEfficiency::ClopperPearson(Int_t(eaccden), Int_t(eaccnum), 0.683, kTRUE) - eacc;

	double downRatio = eacc/nomEff;

	printf("%-20s\t%.4f/%.4f\n",(nIt->first).Data(),upRatio,downRatio);
	allNuisMap.find((nIt->first).Data())->second->find(catnamesbase.at(iCat))->second->insert(std::pair<TString,TString>(procname,TString::Format("%.4f/%.4f", (nomEff > 0. ? upRatio: 1.), (nomEff > 0. ? downRatio: 1.))));
	
      }
    }

  }
  
  std::cout<<" 8TeV analysis done."<<std::endl;
  
  return;
  
}




bool readConfigCard(TString configCardName, 
		    std::map<TString,TString>& procMCfileMap,
		    std::map<TString,bool>&    procOn,
		    std::map<int,TString>&     procIdxMap,
		    int& numProcs, int& numCats,
		    std::vector<TString>& catNames,
		    std::vector<TString>& catDesc) {
  
  FILE* configFile = fopen(configCardName.Data(),"r");
  if ( !configFile ) {
    std::cerr<<" Inputfile "<<configCardName<<" not found."<<std::endl;
    return false;
  }
  
  char line[1000];
  
  std::cout<<" Reading model from file "<<configCardName<<"..."<<std::endl;
    
  int   whichCat = -1 ;
  
  char dummyString[100];
  float dummyFloat;
  int dummy1;
  //# ------------------------------------------------------------------------------------------------------------------------------------------------
  //# Section on the Porcesses/Events
  //#	num proc	num photon cats		num cats(auxiliary)	num cats (analysis)	num cats (trigger)	numModels	numMasses
  //# ------------------------------------------------------------------------------------------------------------------------------------------------
  //INIT	4		4			7			4			16			1		5
  while (fgets(line,1000,configFile)) {
    if ( !sscanf(&line[0], "#") ) {   // not a document line      
      if ( sscanf(line, "INIT %d %d %d", &numProcs, &numCats, &dummy1) ) { // this is the INIT line	
	if( numProcs < 1 || numCats < 1 ) {
	  std::cerr<<" Error in config file "<<configCardName<<" : Number of processes and categories must be positive."<<std::endl;
	  return false;
	}
	
	catNames.resize(numCats);
	catDesc.resize(numCats);
	
	// # analysis categpry definitions
	// #	idx	name		eff. smearing	BG-model/Order	TCuts *use only AUXCATs from above*		StDesc(String descriptor) (*for plots and such*)
	// # ------------------------------------------------------------------------------------------------------------------------------------------------
	// ANACAT	0	cat0		0.005432	Bern/5		" basecut && !vbfcut && baseline0 "		StDesc(Baseline Cat 1)
	// ANACAT	1	cat1		0.005196	Bern/5		" basecut && !vbfcut && baseline1 "		StDesc(Baseline Cat 2)
	// ANACAT	2	cat2		0.007464	Bern/5		" basecut && !vbfcut && baseline2 "		StDesc(Baseline Cat 3)
	// ANACAT	3	cat3		0.012978	Bern/5		" basecut && !vbfcut && baseline3 "		StDesc(Baseline Cat 4)
	// ANACAT	4	cat4		0.008722	Bern/3		" masscut && vbfcut "				StDesc(Baseline VBVA Cat)
      } else if (  sscanf(line, "ANACAT %d %s %f ", &whichCat, dummyString, &dummyFloat) ) { // this is the SMEARING line
	if(whichCat < 0 || whichCat >= numCats) {
	  std::cerr<<" Error in config file "<<configCardName<<" : Setting smearing for Cat "<<whichCat<<" with only "<<numCats<<" total Cats."<<std::endl;
	  return false;
	}
	catNames[whichCat]  = TString(dummyString);
	std::string theLine = line;
	size_t fPos2 = theLine.rfind("Desc(");
	theLine.erase(0,fPos2+5);
	fPos2 = theLine.find_last_of(")");
	theLine.erase(fPos2);
	catDesc[whichCat] = TString(theLine);
      } else {
	
	// additional setup
	int theProc = -1;
	char name[30];
	char onoff[3];
    
	char MCfileName[200];

	//int testmass = 0;

	if ( sscanf(line,"PROC %d %s %s file:%s",&theProc, name, onoff, MCfileName) ) {
	  //# ------------------------------------------------------------------------------------------------------------------------------------------------
	  //PROC    0	ggh	OFF			file:/scratch/fabstoec/cms/hist/hgg-7TeV-janReReco/merged/hgg-7TeV-janReReco_f11--h%dgg-gf-v14b-pu_noskim.root
	  std::cout<<" adding proc with idx = "<<theProc<<"  and name "<<name<<std::endl;
	  procMCfileMap.insert(std::pair<TString,TString>(TString(name)      ,TString(MCfileName)));
	  procIdxMap   .insert(std::pair<int,TString>    (theProc            ,TString(name)));
	  procOn       .insert(std::pair<TString,bool>   (TString(name)      ,strcmp(onoff,"OFF")));

	}
      }
    }
  }
  
  fclose(configFile);
  return true;
  
}




bool readParmsAndCatsFromConfigCard( TString fileName,
				     bool&    computeMVAvar,
				     bool&    correctIDMVAvar,
				     bool&    correctSigEoEvar,
				     TString& mvaWeightFile,
				     TString& mvaDefFile,
				     TString& idmvaCorrFile,
				     TString& sigeoeCorrFile,
				     TString& projectDir,
				     TString& modname,
				     TString& treename,
				     double& massmax, double& massmin,
				     double&theLumi,
				     std::map<TString,TString>& auxCats,
				     std::map<TString,TCut>& anaCats,
				     TCut& baseCut
				     ) {

  FILE* configFile = fopen(fileName.Data(),"r");
  if ( !configFile ) {
    std::cerr<<" Inputfile "<<fileName<<" not found."<<std::endl;
    return false;
  }
  
  char line[1000];
  computeMVAvar = false;
  std::cout<<" Reading paramater from file "<<fileName<<"...";
  
  while (fgets(line,1000,configFile)) {

    if(line[0] == '#') continue;

    char name[400];
    int catIdx = -1;
    char catName[400];
    float massval = -1.;
    if ( sscanf(line,"PROJECTDIR %s",name ) ) projectDir = TString(name);
    else if  ( sscanf(line,"MODNAME %s",name) ) modname = TString(name);
    else if  ( sscanf(line,"TREENAME %s",name) ) treename = TString(name);
    else if( sscanf(line,"LUMI %f",&massval) ) theLumi = (double) massval;
    else if( sscanf(line,"MINMSS %f",&massval) ) massmin = (double) massval;
    else if( sscanf(line,"MAXMSS %f",&massval) ) massmax = (double) massval;
    else if  ( sscanf(line,"AUXCAT %d %s",&catIdx, catName) ) {
      // parsing the ctegory definition like:
      // AUXCAT	0	masscut		" mass>100.0 && mass<180. "
      std::string totLine = line;
      int startCat = totLine.find_first_of("\"");
      int endCat   = totLine.find_last_of("\"");
      std::string catLine = totLine.substr(startCat+1, endCat-startCat-1);
      auxCats.insert(std::pair<TString,TString>(TString(catName),TString(catLine.c_str())));
    }
    else if  ( sscanf(line,"COMPUTEMVA ON %s",name) ) { mvaWeightFile = TString(name); computeMVAvar = true; }
    else if  ( sscanf(line,"CORRECTIDMVA ON %s", name) ) { idmvaCorrFile = TString(name); correctIDMVAvar = true; }
    else if  ( sscanf(line,"CORRECTSIGEOE ON %s", name) ) { sigeoeCorrFile = TString(name); correctSigEoEvar = true; }
    else if  ( sscanf(line,"MVADEFFILE %s",name) ) { mvaDefFile = TString(name);}
    else if  ( sscanf(line,"ANACAT %d %s",&catIdx, catName) ) {
      // ANACAT	0	cat0		" basecut && !vbfcut && baseline0 "		0.005432		Bern/5
      std::string totLine = line;
      int startCat = totLine.find_first_of("\"");
      int endCat   = totLine.find_last_of("\"");
      std::string catLine = totLine.substr(startCat+1, endCat-startCat-1);
      
      // erase starting empties...
      size_t fPos = catLine.find_first_not_of(" ");
      catLine.erase(0,fPos);
      
      // test string for all 
      std::string theCutLine = "";
      while ( catLine.size() > 0 ) {
	fPos = catLine.find_first_not_of(" ()|&!");
	if( fPos != std::string::npos ) {
	  theCutLine.append(catLine.substr(0,fPos));	
	  catLine.erase(0,fPos);
	}
	fPos = catLine.find_first_of(" ()|&!|");
	std::string label;
	if( fPos != std::string::npos ) {
	  label = catLine.substr(0,fPos);
	  catLine.erase(0,fPos);
	  // and erase starting enpties...
	  fPos = catLine.find_first_not_of(" ");
	  catLine.erase(0,fPos);
	} else {
	  label = catLine;
	  catLine = "";
	}
	std::map<TString,TString>::iterator it = auxCats.find(TString(label));
	if ( it == auxCats.end() ){
	  std::cerr<<" ERROR: Could not find AUXCAT with name > "<<label<<" < in list."<<std::endl;
	  return false;
	}
	theCutLine.append(" ( ");
	theCutLine.append( (it->second).Data() );
	theCutLine.append(" ) ");
      }
      anaCats.insert(std::pair<TString,TCut>(TString(catName),TCut(theCutLine.c_str())));
    }
    else if  ( sscanf(line,"BASECAT %s",name) ) {
      std::string totLine = line;
      int startCat = totLine.find_first_of("\"");
      int endCat   = totLine.find_last_of("\"");
      std::string catLine = totLine.substr(startCat+1, endCat-startCat-1);
      
      // erase starting empties...
      size_t fPos = catLine.find_first_not_of(" ");
      catLine.erase(0,fPos);
      
      // test string for all 
      std::string theCutLine = "";
      while ( catLine.size() > 0 ) {
	fPos = catLine.find_first_not_of(" ()|&!");
	if( fPos != std::string::npos ) {
	  theCutLine.append(catLine.substr(0,fPos));	
	  catLine.erase(0,fPos);
	}
	fPos = catLine.find_first_of(" ()|&!|");
	std::string label;
	if( fPos != std::string::npos ) {
	  label = catLine.substr(0,fPos);
	  catLine.erase(0,fPos);
	  // and erase starting enpties...
	  fPos = catLine.find_first_not_of(" ");
	  catLine.erase(0,fPos);
	} else {
	  label = catLine;
	  catLine = "";
	}
	std::map<TString,TString>::iterator it = auxCats.find(TString(label));
	if ( it == auxCats.end() ){
	  std::cerr<<" ERROR: Could not find AUXCAT with name > "<<label<<" < in list."<<std::endl;
	  return false;
	}
	theCutLine.append(" ( ");
	theCutLine.append( (it->second).Data() );
	theCutLine.append(" ) ");
      }
      baseCut = TCut(theCutLine.c_str());
    }
  }
  
  std::cout<<" done"<<std::endl;
  
  fclose(configFile);
  return true;
}


bool readWeightsFromConfigCard(TString fileName, 
			       TH1D*& puweight, TFile*& ptfile, TString& vtxWeightFileName) {
  
  FILE* configFile = fopen(fileName.Data(),"r");
  if ( !configFile ) {
    std::cerr<<" Inputfile "<<fileName<<" not found."<<std::endl;
    return false;
  }
  
  char line[1000];
  
  std::cout<<" Reading weight information from file "<<fileName<<"...";
  
  while (fgets(line,1000,configFile)) {
    // test number of photon categories

    char weightFileName[150];
    
    char onoff[3];
    char name[400];

    if(line[0] == '#') continue;

    if( sscanf(line,"PUREWEIGHFILE ON %s", weightFileName) ) {
      //   //Load pileup weights
      TFile *filepuest = new TFile(weightFileName);
      if( !filepuest ) {
	std::cerr<<" ERROR: Could not open PU file with name "<<weightFileName<<"."<<std::endl;
	fclose(configFile);
	return false;
      }
      TH1D *hpuest = (TH1D*) filepuest->Get("pileup");
      if( !hpuest ) {
	std::cerr<<" ERROR: Could not read histogram <pileup> from PU file with name "<<weightFileName<<"."<<std::endl;
	fclose(configFile);
	return false;
      }

      puweight = hpuest;

    } else if( sscanf(line,"PTREWEIGHFILE ON %s", weightFileName) ) {      
      //load pt-weights
      ptfile = new TFile(weightFileName,"READ");
      if( !ptfile ) {
	std::cerr<<" ERROR: Could not open PU file with name "<<weightFileName<<"."<<std::endl;
	fclose(configFile);
	return false;
      }

    } else if( sscanf(line,"VTXWEIGHTFILE ON %s", weightFileName) ) {
      vtxWeightFileName=TString(weightFileName);
      std::cout<<std::endl<<" [INFO] Vtx reweight fileName = "<<vtxWeightFileName<<std::endl;
    } else if( sscanf(line,"PHEFFSCALESET %s %s {", onoff, name) ) {
      // check if set is switched ON
      bool setIsOn = strcmp(onoff,"OFF");
      float minEta, maxEta, R9split, sf_low, sf_high;
      std::map<std::pair<double,double>, std::vector<double>* >* thisSetsMap = NULL;
      if( setIsOn ) {
	thisSetsMap = new std::map<std::pair<double,double>,std::vector<double>* >();      
	std::cout<<" [INFO] Creating photon efficiency scale-factors set with name < "<<TString(name)<<" > :"<<std::endl;
      } else
	std::cout<<" [WARN] Photon efficiency scale-factors set with name < "<<TString(name)<<" >"<<" present in Card but turned OFF."<<std::endl;
      
      // get next lines until block closes '}'
      while (fgets(line,1000,configFile)) {
	if(line[0] == '}') break;
	if( setIsOn ) {
	  if( sscanf(line," %f %f %f SF(%f) SF(%f)", &minEta, &maxEta, &R9split, &sf_low, &sf_high) ) {
	    std::vector<double>* valVec = new std::vector<double>();
	    valVec->push_back(R9split);
	    valVec->push_back(sf_low);
	    valVec->push_back(sf_high);
	    thisSetsMap->insert( std::pair< std::pair<double, double>, std::vector<double>* > ( std::pair<double,double>(minEta,maxEta), valVec ) );
	    std::cout<<" [INFO]                    Eta-Range [ "<<minEta<<" , "<<maxEta<<" ]  (R9 < "<<R9split<<"):"<<sf_low<<"  (R9 >= "<<R9split<<"):"<<sf_high<<std::endl;
	  } else {	    
	    std::cerr<<" [ERROR] Photon eff. scale factor line not well formated:"<<std::endl<<"   "<<
	      line<<std::endl;
	    return false;
	  }
	  phEffScaleSets.insert( std::pair<TString, std::map< std::pair<double,double>, std::vector<double>* >*> (TString(name), thisSetsMap) );
	}	
      }

//     } else if( sscanf(line,"TRIGCAT %d %d %d %f", &theTrigCatIdx, &thePhCatLabel, &thePhCat2Label, &trigEff) ) {
//       if ( trigWeights->size() == 0 || theTrigCatIdx >= (int) trigWeights->size() ) {
// 	std::cerr<<" ERROR: Trigger Categories not properly set up: NCat = "<<trigWeights->size()<<" and asking for index idx = "<<theTrigCatIdx<<"."<<std::endl;
// 	fclose(configFile);
// 	return false;
//       }
//       (*trigWeights)[theTrigCatIdx] = trigEff;
//       int sumCats = 100*thePhCatLabel+thePhCat2Label;
//       trigArray->insert(std::pair<int,int>(sumCats,theTrigCatIdx));
//     }
    }
  }
  
  std::cout<<std::endl;
  fclose(configFile);
  
  return true;
  
}
