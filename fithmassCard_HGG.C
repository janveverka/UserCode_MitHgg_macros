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


void computeEffSigma( RooAbsData* data, RooRealVar* hmass, double nomMass, double& _lowBound, double& _upBound, double& _CI, double binSize1 = 0.1, double binSize2 = 0.005) {

  double binSize = binSize1;
  double lowBound = nomMass-5.;
  double upBound  = nomMass-5.;
      
  double allNorm = data->sumEntries();
  double theVal = 0.;
  
  double theInterval = 1000;
  
  double theLowBound = 0.;
  double theUpBound  = 0.;
  
  while (lowBound < nomMass) {

    //std::cout<<"  testing low bound = "<<lowBound<<"   nomMass = "<<nomMass<<std::endl;
    
    upBound = lowBound;
    double tmpVal  = 0.;
    while ( tmpVal < (0.6823 * allNorm) && ( upBound < (nomMass+5.)) ) {
      upBound += binSize;
      hmass->setRange("intrange",lowBound,upBound);
      tmpVal = data->sumEntries("1","intrange");
    }
    
    if( tmpVal > (0.6823 * allNorm) ) {
      
      if ( (upBound - lowBound) < theInterval ){
	theLowBound = lowBound;
	theUpBound  = upBound;
	theInterval = (upBound - lowBound);
	theVal      = tmpVal;
      }
      
      lowBound += binSize;
    } else {
      //std::cout<<" Sigma eff: stopping loop, upper bound reached limit: "<<lowBound<<"  -  "<<upBound<<"  with "<<tmpVal/allNorm<<std::endl;
      break;
    }
  }
  
//   std::cout<<"  found decent first bining:"<<std::endl;
//   std::cout<<"  "<<theLowBound<<"  -  "<<theUpBound<<"  with "<<theVal/allNorm<<std::endl;

  // second round with finer binning
  double stopTestLow = theLowBound+3*binSize;
  double startTestUp = theUpBound -5*binSize;
  lowBound = theLowBound-10*binSize;
  binSize = binSize2;
  
  while (lowBound < stopTestLow ) {
    
    //std::cout<<"  testing low bound = "<<lowBound<<std::endl;
    upBound = startTestUp;
    double tmpVal  = 0.;
    while ( tmpVal < (0.6823 * allNorm) && ( upBound < (nomMass+5.) ) ) {
      upBound += binSize;
      
      hmass->setRange("intrange",lowBound,upBound);
      tmpVal = data->sumEntries("1","intrange");
    }
      
    //std::cout<<"      up = "<<upBound<<"  val = "<<tmpVal/allNorm<<"  range = "<<upBound-lowBound<<std::endl;
      
    if( tmpVal > (0.6823 * allNorm) ) {
      
      if ( (upBound - lowBound) < theInterval ){
	//std::cout<<"  found improved binning with "<<lowBound<<"  -  "<<upBound<<"  CL = "<<tmpVal/allNorm<<std::endl;
	theLowBound = lowBound;
	theUpBound  = upBound;
	theInterval = (upBound - lowBound);
	theVal      = tmpVal;
      }
      
      lowBound += binSize;
    } else {
      //std::cout<<" Sigma eff: stopping loop, upper bound reached limit: "<<lowBound<<"  -  "<<upBound<<"  with "<<tmpVal/allNorm<<std::endl;
      break;
    }
  }

  _lowBound = theLowBound;
  _upBound  = theUpBound;
  _CI       = theVal/allNorm;

}


void computeEffSigmaBinned( TH1D* data, double nomMass, double& _lowBound, double& _upBound, double& _CI, double binSize1 = 0.1, double binSize2 = 0.005) {

  double binSize = binSize1;
  double lowBound = nomMass-5.;
  double upBound  = nomMass-5.;
      
  double allNorm = data->Integral();
  double theVal = 0.;
  
  double theInterval = 1000;
  
  double theLowBound = 0.;
  double theUpBound  = 0.;
  
  while (lowBound < nomMass) {

    //std::cout<<"  testing low bound = "<<lowBound<<"   nomMass = "<<nomMass<<std::endl;
    
    upBound = lowBound;
    double tmpVal  = 0.;
    while ( tmpVal < (0.6823 * allNorm) && ( upBound < (nomMass+5.)) ) {
      upBound += binSize;
      tmpVal = data->Integral(data->FindBin(lowBound),data->FindBin(upBound));
    }
    
    if( tmpVal > (0.6823 * allNorm) ) {
      
      if ( (upBound - lowBound) < theInterval ){
	theLowBound = lowBound;
	theUpBound  = upBound;
	theInterval = (upBound - lowBound);
	theVal      = tmpVal;
      }
      
      lowBound += binSize;
    } else {
      //std::cout<<" Sigma eff: stopping loop, upper bound reached limit: "<<lowBound<<"  -  "<<upBound<<"  with "<<tmpVal/allNorm<<std::endl;
      break;
    }
  }
  
//   std::cout<<"  found decent first bining:"<<std::endl;
//   std::cout<<"  "<<theLowBound<<"  -  "<<theUpBound<<"  with "<<theVal/allNorm<<std::endl;

  // second round with finer binning
  double stopTestLow = theLowBound+3*binSize;
  double startTestUp = theUpBound -5*binSize;
  lowBound = theLowBound-10*binSize;
  binSize = binSize2;
  
  while (lowBound < stopTestLow ) {
    
    //std::cout<<"  testing low bound = "<<lowBound<<std::endl;
    upBound = startTestUp;
    double tmpVal  = 0.;
    while ( tmpVal < (0.6823 * allNorm) && ( upBound < (nomMass+5.) ) ) {
      upBound += binSize;
      tmpVal = data->Integral(data->FindBin(lowBound),data->FindBin(upBound));
    }
      
    //std::cout<<"      up = "<<upBound<<"  val = "<<tmpVal/allNorm<<"  range = "<<upBound-lowBound<<std::endl;
      
    if( tmpVal > (0.6823 * allNorm) ) {
      
      if ( (upBound - lowBound) < theInterval ){
	//std::cout<<"  found improved binning with "<<lowBound<<"  -  "<<upBound<<"  CL = "<<tmpVal/allNorm<<std::endl;
	theLowBound = lowBound;
	theUpBound  = upBound;
	theInterval = (upBound - lowBound);
	theVal      = tmpVal;
      }
      
      lowBound += binSize;
    } else {
      //std::cout<<" Sigma eff: stopping loop, upper bound reached limit: "<<lowBound<<"  -  "<<upBound<<"  with "<<tmpVal/allNorm<<std::endl;
      break;
    }
  }

  _lowBound = theLowBound;
  _upBound  = theUpBound;
  _CI       = theVal/allNorm;

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

// float effweight( float sceta, float r9 ) {

//   double scale = 1.;

//   bool iseb   = TMath::Abs(sceta) < 1.5;
//   bool isr9   = (r9 > 0.90);
//   bool isr9_2 = (r9 > 0.94);


//   if( false ) { // old pre-selection...
  
//     // electron veto scale-factros
//     if      (  iseb &&  isr9_2 ) scale *= 0.998;
//     else if (  iseb && !isr9_2 ) scale *= 0.984;
//     else if ( !iseb &&  isr9_2 ) scale *= 0.992;
//     else if ( !iseb && !isr9_2 ) scale *= 0.961;

//     // pre-selection
//     if      (  iseb &&  isr9 ) scale *=  1.003;
//     else if (  iseb && !isr9 ) scale *=  0.996;
//     else if ( !iseb &&  isr9 ) scale *=  0.997;
//     else if ( !iseb && !isr9 ) scale *=  1.009;

//     //BDT cut scale factors
//     if      (  iseb &&  isr9_2 ) scale *= 1.0003;  
//     else if (  iseb && !isr9_2 ) scale *= 1.0002; 
//     else if ( !iseb &&  isr9_2 ) scale *= 0.9980; 
//     else if ( !iseb && !isr9_2 ) scale *= 0.9999; 

//     //   // electron veto scale-factros
//     //   if      (  iseb &&  isr9_2 ) scale *= 0.998;
//     //   else if (  iseb && !isr9_2 ) scale *= 0.984;
//     //   else if ( !iseb &&  isr9_2 ) scale *= 0.992;
//     //   else if ( !iseb && !isr9_2 ) scale *= 0.961;

//     //   // pre-selection * BDT cut factors
//     //   if      (  iseb &&  isr9 ) scale *= 0.998 * 1.0001;
//     //   else if (  iseb && !isr9 ) scale *= 1.002 * 0.9995;
//     //   else if ( !iseb &&  isr9 ) scale *= 1.008 * 0.9998;
//     //   else if ( !iseb && !isr9 ) scale *= 0.996 * 0.9968;

//   } else {
  
//     // electron veto scale-factros
//     if      (  iseb &&  isr9_2 ) scale *= 0.998;
//     else if (  iseb && !isr9_2 ) scale *= 0.984;
//     else if ( !iseb &&  isr9_2 ) scale *= 0.992;
//     else if ( !iseb && !isr9_2 ) scale *= 0.961;

//     // pre-selection
//     if      (  iseb &&  isr9 ) scale *=  1.000;
//     else if (  iseb && !isr9 ) scale *=  0.984;
//     else if ( !iseb &&  isr9 ) scale *=  1.003;
//     else if ( !iseb && !isr9 ) scale *=  0.997;

//     //BDT cut scale factors
//     if      (  iseb &&  isr9_2 ) scale *= 1.0000;  
//     else if (  iseb && !isr9_2 ) scale *= 0.9987; 
//     else if ( !iseb &&  isr9_2 ) scale *= 1.0000; 
//     else if ( !iseb && !isr9_2 ) scale *= 0.9978; 

//   }

//   return scale;

// }

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
      
      //std::cout<<"   -> "<< _masserrsmeared/_mass<<"       "<<_masserrsmearedwrongvtx/_mass<<std::endl;
      
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
bool validateInput(int iProc, int iCat, int nProcs, int nCats);

bool readWeightsFromConfigCard(TString fileName, 
			       TH1D*& puhisto, TFile*& ptweight, TString& vtxWeightFile);

bool readConfigCard(TString configCardName, 
		    std::map<TString,TString>& procMCfileMap,
		    std::map<TString,bool>&    procOn,
		    std::map<int,TString>&     procIdxMap,			
		    int& numProcs, int& numCats, int& numMasses,
		    bool*& catIsOn,
		    std::vector<double>& mhs,
		    std::vector<int>   & mhs_fit,
		    std::vector<double>& smearingv,
		    std::vector<TString>& catNames,
		    std::vector<TString>& catDesc,
		    std::map<TString,TString>& catToSmearMap,
		    std::map<TString,float>*& startVals,
		    std::map<TString,int>*& startValsFixed
		    );



bool readConfigCardNuissances(TString configCardName, int numCats,
			      std::vector<RooAbsReal*>& nsigcat,
			      std::vector<RooAbsReal*>& nsigcat_2nd,
			      std::vector<RooAbsReal*>& nsigcatModel,
			      std::vector<RooAbsReal*>& nuissances,
			      std::vector<RooAbsReal*>& finalnorm,
			      std::vector<RooAbsReal*>& finalnorm_2nd,
			      std::vector<RooAbsReal*>& finalnormModel,
			      std::vector<TString>      catnames);


RooAbsPdf* generateMultiGaussian(RooRealVar* mass, 
				 RooRealVar* nomMass,
				 TString procName, TString quali,
				 std::map<TString,RooRealVar*>& fitparms,
				 std::map<TString,float>&       startvals,
				 std::map<TString,int>&         valsfixed,
				 bool rightVtx = true  );

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
				     TString& nosmearmodname,
				     TString& nosmeartreename,
				     double& massmax, double& massmin,
				     double& fitmassmax, double& fitmassmin,
				     double& lumi,
				     std::map<TString,TString>& auxCats,
				     std::map<TString,TCut>& anaCats,
				     TCut& theBaseCut
				     );

bool resetStartValues(std::map<TString,RooRealVar*>& fitparms,
		      std::map<TString,float>&       startvals);


RooDataSet *makedset(TString name, TTree *tree, TCut cut, RooRealVar *hmass, RooRealVar *weight, double nomMass, bool& empty, double moreWeights = 1., TH1D* hist = NULL) {
 
  RooDataSet *dset = new RooDataSet(name,"",RooArgSet(*hmass,*weight),weight->GetName());
  RooRealVar *vset = (RooRealVar*)dset->get()->find(hmass->GetName());
  
//   std::cout<<" -----------------------------------------"<<std::endl;
//   std::cout<<" Creating dataset with name = "<<name<<" :"<<std::endl;
//   std::cout<<"   "<<cut<<std::endl;
//   std::cout<<" -----------------------------------------"<<std::endl;

  tree->SetEstimate(tree->GetEntries());
  Int_t nev = tree->Draw("mass",cut,"goff");

//   TString massrecomp = " TMath::Sqrt(ph1.e*ph2.e*(1-gencostheta)*2.) ";
//   Int_t nev = tree->Draw(massrecomp.Data(),cut,"goff");

  double *vals = tree->GetV1();
  double *weights = tree->GetW();
  
  for (int iev=0; iev<nev; ++iev) {
    vset->setVal(vals[iev]);
    dset->add(*vset,weights[iev]*moreWeights);

    if ( hist )
      hist->Fill(vals[iev], weights[iev]*moreWeights);
  }
  
  empty = false;

  double dummyWeight = 0.00000001;
  // if there's no event, add one dummy event at the nominal mass with small weight...
  if (nev == 0) {
    empty = true;
    std::cout<<" WARNING: dataset with name "<<name<<" has not events. The Cuts are: "<<std::endl;
    std::cout<<"         Cut = "<<cut<<std::endl;
    std::cout<<"   adding a dummy Event with weight "<<dummyWeight<<" @mass = "<<nomMass<<std::endl;
    vset->setVal(nomMass);
    dset->add(*vset,dummyWeight);
    vset->setVal(nomMass-2.);
    dset->add(*vset,dummyWeight/2.);
    vset->setVal(nomMass+2.);
    dset->add(*vset,dummyWeight/2.);
  }

  return dset;
  
}

// ============================================================================================
#define DOBDT


//void createSignalModels(TString configCardName="mettag.config", modeltype model = STANDARDMODEL) {
void fithmassCard_HGG(bool fitonly        = false, 
		      bool doeffsigma     = true,
		      bool doeffsmear     = true,
		      bool doSecondSignal = false,
		      TString baseString  = "moriond13",
#ifdef DOBDT
		      //TString configCardName="/home/fabstoec/cms/root/templateHGG_8TeV_Moriond.config") {
		      TString configCardName="templateHGG_8TeV_Moriond.config") {
#else
                      TString configCardName="/home/fabstoec/cms/root/templateHGG_8TeV_CiC_Moriond.config") {
#endif

  //************************************************************************************
  // SETUP PART
  // -----------------------------------------------------
  std::map<TString,double*> processCrossSectionMap;
  initProcessCSArrayMap(processCrossSectionMap);

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

  std::vector<const char*>* procxseclist_names = NULL;
  
  procxseclist_names = &smlist;
  //procxseclist_names = &fflist;
  
  //   switch(model) {
  //   case STANDARDMODEL:
  //     procxseclist_names = &smlist;
  //     break;
  //   case STANDARDMODEL4:
  //     procxseclist_names = &sm4list;
  //     break;
  //   case FERMIOPHOBIC:
  //     procxseclist_names = &fflist;
  //     break;
  //   default:
  //     std::cerr<<" Model "<<model<<" not implemented."<<std::endl;
  //     return;
  //   }
  
  //************************************************************************************

  // number of processes, Categories, Models and mass points
  int numProcs  = -1;
  int numCats   = -1;
  int numMasses = -1;

  // decide if category is switched ON/OFF
  bool* catIsOn = NULL;

  // process names and toggle process ON/OFF
  //std::vector<TString> procnames;
  //std::vector<bool>    procon;

  std::map<int,TString>     procIdxMap;
  std::map<TString,TString> procMCfileMap;
  std::map<TString,bool>    procon;
  
  // smearing values FIX-ME: compute on the fly...
  std::vector<double> smearingv;

  // define the masses we're testing to interpolate the models between
  std::vector<double> mhs;
  std::vector<int>    mhs_fit;
  double massmax = -1.;
  double massmin = -1.;

  double fitmassmax = -1.;
  double fitmassmin = -1.;
  

  // auxiliary vectro for base cat names...
  std::vector<TString> catnamesbase;
  std::vector<TString> catdesc;


  TString    projectDir;
  TString    modname;
  TString    treename;
  TString    nosmearmodname;
  TString    nosmeartreename;


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
  
  std::map<TString,float>      * startvals = NULL; //new std::map<TString,float>      [numCats];
  std::map<TString,int  >      * valsfixed = NULL; //new std::map<TString,float>      [numCats];



  std::map<TString,TString> catToSmearMap;
//   catToSmearMap.insert(std::pair<TString,TString>("hgg_8TeV_hcp_mtag","cat4"));
//   catToSmearMap.insert(std::pair<TString,TString>("hgg_8TeV_hcp_etag","cat4"));
//   catToSmearMap.insert(std::pair<TString,TString>("hgg_8TeV_hcp_dijetloose","cat4"));
//   catToSmearMap.insert(std::pair<TString,TString>("hgg_8TeV_hcp_dijettight","cat4"));
//   catToSmearMap.insert(std::pair<TString,TString>("hgg_8TeV_hcp_mettag","cat4"));
//   catToSmearMap.insert(std::pair<TString,TString>("hgg_8TeV_hcp_bdt0","cat0"));
//   catToSmearMap.insert(std::pair<TString,TString>("hgg_8TeV_hcp_bdt1","cat1"));
//   catToSmearMap.insert(std::pair<TString,TString>("hgg_8TeV_hcp_bdt2","cat2"));
//   catToSmearMap.insert(std::pair<TString,TString>("hgg_8TeV_hcp_bdt3","cat3"));

  bool status = readConfigCard(configCardName, 
			       procMCfileMap,
			       procon,
			       procIdxMap,
			       numProcs, numCats, numMasses,
			       catIsOn,
			       mhs,
			       mhs_fit,
			       smearingv,
			       catnamesbase,
			       catdesc,
			       catToSmearMap,
			       startvals,
			       valsfixed);


  if(!status) {
    std::cerr<<" ERROR when readin input card "<<configCardName<<"."<<std::endl;
    return;
  }

  // compute how many mh we're actuallt fitting
  int numFitMasses    = 0;
  double minFitMass   = 1000.;
  double maxFitMass   = 0.;
  double fitBinSize   = 0.;

  bool lastWasFit = false;

  for( unsigned int iM =0; iM < mhs.size(); ++iM ) {
    if( mhs_fit[iM] ) {
      if( lastWasFit ) fitBinSize = mhs[iM] - mhs[iM-1];
      if( mhs[iM] < minFitMass ) minFitMass = mhs[iM];
      if( mhs[iM] > maxFitMass ) maxFitMass = mhs[iM];
      numFitMasses++;
      lastWasFit = true;
    } else
      lastWasFit = false;
  }
  
  typedef std::map<double,double>        massValMap;
  typedef std::map<TString,massValMap*>  catValMap;
  typedef std::map<TString,catValMap*>   procValMap;

  typedef std::pair<double,double>       massValPair;
  typedef std::pair<TString,massValMap*> catValPair;
  typedef std::pair<TString,catValMap*>  procValPair;

  // keep track of eff x acc and Nevents per mass/cat/proc
  procValMap effaccMap;
  procValMap neventsMap;

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
    TFile* vtxWeightFile = TFile::Open("/home/fabstoec/cms/cmssw/029/CMSSW_5_3_2_patch4/src/UserCode/HiggsAnalysis/HiggsTo2photons/h2gglobe/Macros/vertex_reweighing_mva_HCP2012_unblind.root");
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
					   nosmearmodname             ,
					   nosmeartreename            ,
					   massmax, massmin    ,
					   fitmassmax, fitmassmin    ,
					   lumi,
					   auxCatMap           ,
					   anaCatMap           ,
					   theBaseCut
					   );
  

  if( !status ) {
    std::cerr<<" ERROR: Could not read Category information from file "<<configCardName.Data()<<std::endl;
    return;
  }


  std::cout<<" Checking for directory "<<projectDir+TString("/model")<<"..."<<std::endl;
  if ( gSystem->mkdir(TString::Format("%s/model",projectDir.Data()), true) < 0 )
    std::cerr<<" WARNING: Directory "<<TString::Format("%s/model",projectDir.Data()) <<" does not exists. Creating it..."<<std::endl;
    
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
    sigEoEscaleFile = TFile::Open("/home/fabstoec/cms/root/PhotonIDMVA_new/Moriond13_phSigEoE.root");

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

  // output workspace ... again: name configuirable ? FIX-ME
  RooWorkspace* wOut = new RooWorkspace("wsig","");
  
  RooRealVar*   hmass = new RooRealVar("CMS_hgg_mass","m_{#gamma#gamma}",massmin,massmax,"GeV");
  hmass->setRange(massmin,massmax);
  hmass->setBins( (int) (1.0*(massmax-massmin)) );
  hmass->SetTitle("m_{#gamma#gamma}");
  hmass->setUnit("GeV");  

  // This is a second hmass, used creating the signal models...
  RooRealVar*   hmass2 = new RooRealVar("model_mass","m_{#gamma#gamma}",massmin,massmax,"GeV");
  hmass2->SetTitle("m_{#gamma#gamma}");
  hmass2->setUnit("GeV");  

  RooRealVar mnom("MH","m_{h}",110.,massmin,massmax,"GeV");
  mnom.setConstant();
  
  RooRealVar mnom2("MH_BG","m_{h}",110.,massmin,massmax,"GeV");
  mnom2.setVal(125.);
  mnom2.setConstant();

  RooRealVar the2ndMu("CMS_hgg_2ndMu","2ndMu",-100,100);
  the2ndMu.setVal(1.0);
  the2ndMu.setConstant();

  //---RooRealVar to load dataset or used for calculation---
  RooRealVar *weight = new RooRealVar("weight","",1.0);//this seems to be useless variable
  weight->removeRange();

  // maps for the histFuncs and the Additions
  std::map<TString,RooHistFunc*> histFuncMap;
  std::map<TString,RooHistFunc*> histFuncMap2nd;

  std::map<TString,RooAddition*> addFuncMap;
  std::map<TString,RooAddition*> addFuncMap2nd;
  
  int numHists = modelParmsMap->size();
  TH1D** parmHists     = new TH1D*       [numHists];

  RooDataHist** parmDataHists = new RooDataHist*[numHists];
  RooHistFunc** parmHistFuncs = new RooHistFunc*[numHists];

  // stuff for second Higgs model (to be used as BG)
  RooDataHist** parmDataHists2nd = new RooDataHist*[numHists];
  RooHistFunc** parmHistFuncs2nd = new RooHistFunc*[numHists];

  int iPair = 0;
  for(std::map<TString,double*>::iterator pair=modelParmsMap->begin(); pair!=modelParmsMap->end(); ++pair, ++iPair) {
    TString name = pair->first;
    double* data = pair->second;
    //parmHists[iPair] = new TH1D(histName.str().c_str(),"",numsmpoints,smmasses[0]-0.25,smmasses[numsmpoints-1]+0.25);
    parmHists[iPair] = new TH1D(TString::Format("h%ss_%s",name.Data(),baseString.Data()).Data(),"",numsmpoints,smmasses[0]-0.25,smmasses[numsmpoints-1]+0.25);

    for(int ipoint=0; ipoint<numsmpoints; ++ipoint)
      parmHists[iPair]->Fill(smmasses[ipoint],data[ipoint]);

    parmDataHists[iPair] = new RooDataHist(TString::Format("d%ss_%s",name.Data(),baseString.Data()).Data(),"",RooArgList(mnom),parmHists[iPair]);
    parmHistFuncs[iPair] = new RooHistFunc(TString::Format("f%ss_%s",name.Data(),baseString.Data()).Data(),"",RooArgList(mnom),(*parmDataHists[iPair]),1);
    histFuncMap.insert(std::pair<TString,RooHistFunc*>(name,parmHistFuncs[iPair]));

    if( doSecondSignal ) {
      parmDataHists2nd[iPair] = new RooDataHist(TString::Format("d%ss_%s_2nd",name.Data(),baseString.Data()).Data(),"",RooArgList(mnom2),parmHists[iPair]);      
      parmHistFuncs2nd[iPair] = new RooHistFunc(TString::Format("f%ss_%s_2nd",name.Data(),baseString.Data()).Data(),"",RooArgList(mnom2),(*parmDataHists2nd[iPair]),1);
      histFuncMap2nd.insert(std::pair<TString,RooHistFunc*>(name,parmHistFuncs2nd[iPair]));
    }
  }
  

  // create the summed up functions
  for(unsigned int iSum=0; iSum<sumRuleNames.size(); iSum++) {
    std::vector<const char*>* comps = sumRuleComps[iSum];
    TString sumName                 = sumRuleNames[iSum];
    RooArgList sumList;
    RooArgList sumList2nd;
    for(unsigned int iComp=0; iComp<comps->size(); ++iComp) {
      RooHistFunc* tempFunc = NULL;

      if( histFuncMap.find((TString)((*comps)[iComp])) != histFuncMap.end() )
	tempFunc = histFuncMap.find(  (TString)((*comps)[iComp]) )->second;
      if( !tempFunc ) {
	std::cerr<<" Cannot find histfunc with name "<<TString( comps->at(iComp) )<<"."<<std::endl;
	return;
      }
      sumList.add( *tempFunc );

      if( doSecondSignal ) {
	if( histFuncMap2nd.find((TString)((*comps)[iComp])) != histFuncMap2nd.end() )
	  tempFunc = histFuncMap2nd.find(  (TString)((*comps)[iComp]) )->second;
	if( !tempFunc ) {
	  std::cerr<<" Cannot find histfunc with name "<<TString( comps->at(iComp) )<<"."<<std::endl;
	  return;
	}
	sumList2nd.add( *tempFunc );
      }
    }
    RooAddition* tempAdd = new RooAddition(TString::Format("f%ss_%s",sumName.Data(),baseString.Data()).Data(),"",sumList);
    addFuncMap.insert(std::pair<TString,RooAddition*>(sumName,tempAdd));

    if( doSecondSignal ) {
      RooAddition* tempAdd2nd = new RooAddition(TString::Format("f%ss_%s_2nd",sumName.Data(),baseString.Data()).Data(),"",sumList);
      addFuncMap2nd.insert(std::pair<TString,RooAddition*>(sumName,tempAdd2nd));
    }
  }
  

  std::vector<RooAbsReal*> procxseclist;
  for(unsigned int iComp=0; iComp < procxseclist_names->size(); ++iComp) {

//     stringstream pSS;
//     pSS <<  TString(procxseclist_names->at(iComp));

    TString theName = procxseclist_names->at(iComp);
    
    RooAbsReal* tempAbs = NULL;
    //if( addFuncMap.find((TString) pSS.str().c_str()) != addFuncMap.end() ) {
    if( addFuncMap.find( theName ) != addFuncMap.end() ) {
      tempAbs = addFuncMap.find( theName )->second;
    } else if ( histFuncMap.find( theName ) != histFuncMap.end() ) {
      tempAbs = histFuncMap.find( theName )->second;
    }
    if( !tempAbs ) {
      std::cerr<<" Cannot find addfunc with name "<<theName<<"."<<std::endl;
      return;
    }
    procxseclist.push_back(tempAbs);
  }
  
  std::vector<TTree*> trees;
  std::vector<TTree*> treesNS;

  std::vector<double> ntots;
  std::vector<TFile*> files;

  // open the MVA files, if requested
  TFile* tmvaOutput = NULL;
  TString weights   = "";
  std::vector<std::string> *varlist = NULL;

  if( computeMVAvar ) {
    tmvaOutput = new TFile(mvaDefFile.Data(),"READ");//root file to store training information  
    varlist = (std::vector<std::string>*)tmvaOutput->Get("varlist");
    weights = mvaWeightFile;
  }  
  
  // ========================= BEGIN effective smearing ============================
  // keep a RooDataSet for each mass and each Cat
  RooDataSet*** effSmearData    = new RooDataSet**[numCats+1];
  RooDataSet*** effSmearDataNS  = new RooDataSet**[numCats];

  for(int iCat = 0; iCat < numCats + 1; ++iCat) {
    
    effSmearData   [iCat] = new RooDataSet*[mhs.size()];
    if( iCat < numCats )
      effSmearDataNS [iCat] = new RooDataSet*[mhs.size()];
    
    for( unsigned int iMass=0; iMass < mhs.size();  ++iMass ) {
      effSmearData   [iCat][iMass] = new RooDataSet(TString::Format("effsmearset_%d_%d",   (int)mhs[iMass],iCat).Data(),"",RooArgSet(*hmass2,*weight),weight->GetName());
      if( iCat < numCats )
	effSmearDataNS [iCat][iMass] = new RooDataSet(TString::Format("effsmearsetNS_%d_%d", (int)mhs[iMass],iCat).Data(),"",RooArgSet(*hmass2,*weight),weight->GetName());
    }
  }


  RooAddPdf** effSmearPdf    = new RooAddPdf*[numCats+1];
  RooArgList* addPdfList     = new RooArgList[numCats+1];
  RooArgList* addPdfCoefList = new RooArgList[numCats+1];
  // ========================= END   effective smearing ============================


  std::map<TString,RooDataSet*> keepAllDataMap;


  TH1D*** effSmearHists   = new TH1D**[mhs.size()];
  TH1D*** effSmearHistsNS = new TH1D**[mhs.size()];

  //TCanvas** tmpCan = new TCanvas*[numCats];

  // loop over all processes that are on
  TFile *friendtmp = new TFile("friendtmp.root","RECREATE");
  friendtmp->ls();
  for(int iProc = 0; iProc < numProcs; ++iProc) {

    // get the proc name from the map
    TString procname = procIdxMap.find(iProc)->second;
    if( !(procon.find(procname)->second) ) continue;


    // add an entry for the maps to store eff x acc and nevents
    effaccMap.insert(procValPair(procname, new catValMap()));
        
    // load all the tress for this process...
    trees.resize(0);
    ntots.resize(0);
    for (UInt_t i=0; i<mhs.size(); ++i) {
      double masspoint = mhs.at(i);
      

      // load the correct tree from the correct file...
      std::map<TString,TString>::iterator tempIt = procMCfileMap.find(procname);
      if( tempIt == procMCfileMap.end() ) {
	std::cout<< " ERROR *** "<<std::endl;
	return;
      }
      
      TString samplestring = tempIt->second;      

      std::cout<<" Opening file for mass "<<masspoint<<": "<<TString::Format(samplestring, (int) masspoint)<<std::endl;
      
      //TFile *theFile = new TFile(TString::Format(samplestring, (int) masspoint),"READ");
      TFile *theFile = TFile::Open(TString::Format(samplestring, (int) masspoint),"READ");
      if (!theFile) {
	std::cerr<<" ERROR: Could not open file "<<TString::Format(samplestring, (int) masspoint).Data()<<"."<<std::endl;
	return;
      }
      
      effSmearHists  [i] = new TH1D*[numCats];
      effSmearHistsNS[i] = new TH1D*[numCats];
      for(int iCat=0; iCat < numCats; ++iCat) {
	//effSmearHists  [iCat] = new TH1D(TString::Format("effsmearHist_%s",catnamesbase.at(iCat).Data()).Data(),"",(int) ((massmax-massmin)/0.001), massmin, massmax);
	effSmearHistsNS[i][iCat] = new TH1D(TString::Format("effsmearHistNS_%d_%s",i,catnamesbase.at(iCat).Data()).Data(),"",(int) ((massmax-massmin)/0.001), massmin, massmax);
      }


      // load the Tree
      TDirectory *hdir = (TDirectory*) theFile->FindObjectAny(modname);
      TH1D *hallevts   = (TH1D*)       theFile->Get("AnaFwkMod/hDAllEvents");
      ntots.push_back(hallevts->GetSumOfWeights());
      
      TTree *tree    = (TTree*)hdir->Get(treename.Data());
      TTree *treeNS  = NULL;
      if(doeffsmear) {
	TDirectory *nsdir = (TDirectory*) theFile->FindObjectAny( nosmearmodname );
	treeNS = (TTree*)nsdir->FindObjectAny(nosmeartreename.Data());  // FIX-ME
      }

      if( computeMVAvar ) {
	ApplyAsFriend(tree,weights,*varlist,"bdt",correctIDMVAvar,correctSigEoEvar);
	if(doeffsmear)
	  ApplyAsFriend(treeNS,weights,*varlist,"bdt",correctIDMVAvar,correctSigEoEvar);
      }

      trees.push_back( tree );
      if(doeffsmear)
	treesNS.push_back( treeNS );

      files.push_back(theFile);

    }

    TString projProcDir = TString::Format("%s/%s",projectDir.Data(), procname.Data());
    std::cout<<" Checking for directory "<<projProcDir<<"..."<<std::endl;
    
    if ( gSystem->mkdir( projProcDir, true) == 0 )
      std::cerr<<" WARNING: Directory "<<projProcDir<<" does not exists. Creating it..."<<std::endl;
    
    // full Cat names including the process name
    std::vector<TString> catnames;
    for(int iCat = 0; iCat < numCats; ++iCat) {
      stringstream pSS;
      pSS << catnamesbase[iCat].Data() << "_" << procname;
      catnames.push_back(pSS.str().c_str());
    }
    
    // collection for the fitparameters
    std::map<TString,RooRealVar*>* fitparms  = new std::map<TString,RooRealVar*>[numCats];
    
    // generate the RooAbsPdf s for the right and wrong Vtx hypothesis    
    RooAbsPdf** combh_rv      = new RooAbsPdf*[numCats];    // right hypothesis
    RooAbsPdf** combh_wv      = new RooAbsPdf*[numCats];    // right hypothesis
    
    std::cout<<"  #Categories = "<<numCats<<std::endl;
    
    for(int iCat=0; iCat < numCats; ++iCat) {
      
      combh_rv[iCat]      = generateMultiGaussian(hmass2,&mnom,
						  procname,catnamesbase[iCat],fitparms[iCat],startvals[iProc*numCats+iCat],valsfixed[iProc*numCats+iCat]);
      combh_wv[iCat]      = generateMultiGaussian(hmass2,&mnom,
						  procname,catnamesbase[iCat],fitparms[iCat],startvals[iProc*numCats+iCat],valsfixed[iProc*numCats+iCat], false );
      if (!combh_rv[iCat] || !combh_wv[iCat]) return;
      
    }

    //parameters calculated from mc events rather than fitting: fracright and effacc
    RooRealVar** effacc      = new RooRealVar*[numCats];
    RooRealVar** effaccModel = new RooRealVar*[numCats];   // needed for mdeling, mass cuts not the same...
    
    // right/wrong Vtx fraction
    RooRealVar** fracright      = new RooRealVar*[numCats];
    RooRealVar** fracrightModel = new RooRealVar*[numCats];

    // generate the combined PDFs
    RooAddPdf** combhvtx   = new RooAddPdf* [numCats];

    for(int iCat = 0; iCat < numCats; ++iCat) {
      
      TString fracrightName   = TString::Format("fracright_%s",catnames[iCat].Data());
      fracright[iCat]         = new RooRealVar(fracrightName.Data(),"fracright",   0.9,0.0,1.0);
      
      fitparms[iCat].insert(std::pair<TString,RooRealVar*>(fracrightName,fracright[iCat]));
      startvals[iProc*numCats+iCat].insert(std::pair<TString,float>(fracrightName,0.9));
      
      fracrightName   = TString::Format("fracrightModel_%s",catnames[iCat].Data());
      fracrightModel[iCat]    = new RooRealVar(fracrightName.Data(),"fracrightModel",   0.9,0.0,1.0);            
      fitparms[iCat].insert(std::pair<TString,RooRealVar*>(fracrightName,fracright[iCat]));
      startvals[iProc*numCats+iCat].insert(std::pair<TString,float>(fracrightName,0.9));

      TString effaccName   = TString::Format("effacc_%s",catnames[iCat].Data());
      effacc[iCat]         = new RooRealVar(effaccName.Data(),"effacc",   0.9,0.0,1.0);
      
      fitparms[iCat].insert(std::pair<TString,RooRealVar*>(effaccName,effacc[iCat]));
      startvals[iProc*numCats+iCat].insert(std::pair<TString,float>(effaccName,0.9));
      
      effaccName   = TString::Format("effaccModel_%s",catnames[iCat].Data());
      effaccModel[iCat]    = new RooRealVar(effaccName.Data(),"effaccModel",   0.9,0.0,1.0);            
      fitparms[iCat].insert(std::pair<TString,RooRealVar*>(effaccName,effacc[iCat]));
      startvals[iProc*numCats+iCat].insert(std::pair<TString,float>(effaccName,0.9));

      // This is the combined PDF  (FIX-ME:  Possible name problem ???)
      combhvtx[iCat] = new RooAddPdf(TString::Format("combhvtx_%s_%d",procname.Data(),iCat).Data(),"combhvtx",RooArgList( (*combh_rv[iCat]),(*combh_wv[iCat]) ), RooArgList( (*fracright[iCat]) ));
      
    }
    
    //define histograms to keep track of fit parameters for each category each mass in the step of 10 GeV
    TH1F ***fitparmhists = new TH1F**[numCats];
    for (UInt_t iCat=0; iCat<catnames.size(); ++iCat) {
      fitparmhists[iCat] = new TH1F*[fitparms[iCat].size()];
      unsigned int iparm = 0;
      for(std::map<TString,RooRealVar*>::iterator pair = fitparms[iCat].begin(); pair != fitparms[iCat].end(); ++pair) {	
	TString histname = TString("hist")+ pair->first + catnames.at(iCat);
	//fitparmhists[iCat][iparm] = new TH1F(histname,histname,(int) numMasses,mhs[0]-(mhs[1]-mhs[0])/2.,mhs[numMasses-1]+(mhs[1]-mhs[0])/2.);
	//fitparmhists[iCat][iparm] = new TH1F(histname,histname,(int) numMasses,mhs[0]-2.5,mhs[mhs.size()-1]+2.5);
	fitparmhists[iCat][iparm] = new TH1F(histname,histname,(int) numFitMasses,minFitMass-fitBinSize/2.,maxFitMass+fitBinSize/2.);
	iparm++;
      }
    }
    
    //---define 2d hist to fill the fit status for each mass and each category---
    TH2F *histfitstatus_rv = new TH2F("histfitstatus_rv","histfitstatus_rv",(int) numMasses,mhs[0]-5.,mhs[numMasses-1]+5.,catnames.size(),-0.5,catnames.size()-0.5);//0 fit; 1 not fit
    TH2F *histfitstatus_wv = new TH2F("histfitstatus_wv","histfitstatus_wv",(int) numMasses,mhs[0]-5.,mhs[numMasses-1]+5.,catnames.size(),-0.5,catnames.size()-0.5);//0 fit; 1 not fit
    

    //---add canvas for plot---
    TCanvas* rightCan = new TCanvas("rightCan","right vertex",numMasses*250,numCats*250); 
    rightCan->Divide(numMasses,numCats);
    TCanvas* wrongCan = new TCanvas("wrongCan","wrong vertex",numMasses*250,numCats*250); 
    wrongCan->Divide(numMasses,numCats);
//     TCanvas* allCan = new TCanvas("allCan","combined",numMasses*250,numCats*250); 
//     allCan->Divide(numMasses, numCats);
    
    TCanvas* dummy = new TCanvas();
    
    // =========== EFFECTIVE SMAERING ==================
    // keep all datasets, also the un-smeared ones.
    for (UInt_t iCat=0; iCat<catnames.size(); ++iCat) {
      
      if ( !catIsOn[iProc*numCats+iCat] ) continue;  // don't do fits in this Cat
      
      // keep track of eff x acc /nevents
      effaccMap.find(procname)->second->insert(catValPair(catnamesbase[iCat],new massValMap()));
      //neventsMap.find(procname)->second->insert(catname[iCat],new massValMap());

      int fitMass = -1;

      for (UInt_t i=0; i<mhs.size(); ++i) {
	std::cout<<" MH = "<<mhs[i]<<std::endl;
	
	double masspoint = mhs.at(i);

	//mass points string
	std::stringstream numstringstr;
	numstringstr<<masspoint;
	TString numstring(numstringstr.str());
	//get dataset for this category for this mass point
	
	TString dataAllName;
	
	dataAllName=TString::Format("sig_%s_mass_m",procname.Data())+numstring+TString("_")+catnamesbase.at(iCat);
	
	// set the PU weight correctly...
	if( hpuest )
	  setpuweights(files[i], hpuest);
	
	// find the TCut string for the category we're in
	
	std::map<TString,TCut>::iterator catIt = anaCatMap.find(catnamesbase[iCat]);
	if( catIt == anaCatMap.end() ) {
	  std::cerr<<" ERROR: Cannot find TCut String for ANACAT with name "<<catnamesbase[iCat]<<"."<<std::endl;
	  return;
	}


        //TCut theWeight = "puweight(numPU)*bsweight(vtxZ,genHiggsZ)*trigeffscale()*effweight(ph1.sceta,ph1.r9)*effweight(ph2.sceta,ph2.r9)*vtxWeight(ptgg,vtxZ,genHiggsZ)";
        TCut theWeight = "puweight(numPU)*bsweight(vtxZ,genHiggsZ)*trigeffscale()*effweightCard(ph1.sceta,ph1.r9)*effweightCard(ph2.sceta,ph2.r9)*vtxWeight(ptgg,vtxZ,genHiggsZ)";
        //TCut theWeight = "puweight(numPU)*trigeffscale()*effweight(ph1.sceta,ph1.r9)*effweight(ph2.sceta,ph2.r9)";
	
 	TCut masscutModel   =  theBaseCut && TCut(TString::Format("mass > %f && mass < %f",massmin,massmax).Data());
 	TCut masscutReal    =  theBaseCut && TCut(TString::Format("mass > %f && mass < %f",massmin,massmax).Data());

	
	TCut rightVtx = "(TMath::Abs(genHiggsZ-vtxZ)<1.0)";
	TCut wrongVtx = !rightVtx;


	bool isDataSetEmpty = false;

	RooDataSet *mcsigallwdata_rv        = makedset(dataAllName+TString("_rv"),   trees[i], theWeight*(masscutReal  && catIt->second && rightVtx), hmass, weight, mhs[i], isDataSetEmpty);
	RooDataSet *mcsigallwdata_wv        = makedset(dataAllName+TString("_wv"),   trees[i], theWeight*(masscutReal  && catIt->second && wrongVtx), hmass, weight, mhs[i], isDataSetEmpty);
	RooDataSet *mcsigallwdata           = makedset(dataAllName,                  trees[i], theWeight*(masscutReal  && catIt->second            ), hmass, weight, mhs[i], isDataSetEmpty);
	
	dataAllName=TString::Format("sig_%s_massModel_m",procname.Data())+numstring+TString("_")+catnamesbase.at(iCat);

	RooDataSet *mcsigallwdataModel_rv   = makedset(dataAllName+TString("_rv"),   trees[i], theWeight*(masscutModel && catIt->second && rightVtx), hmass2, weight, mhs[i], isDataSetEmpty);
	RooDataSet *mcsigallwdataModel_wv   = makedset(dataAllName+TString("_wv"),   trees[i], theWeight*(masscutModel && catIt->second && wrongVtx), hmass2, weight, mhs[i], isDataSetEmpty);
	RooDataSet *mcsigallwdataModel      = makedset(dataAllName,                  trees[i], theWeight*(masscutModel && catIt->second            ), hmass2, weight, mhs[i], isDataSetEmpty);


	if( !mcsigallwdata_rv ||
	    !mcsigallwdata_wv ||
	    !mcsigallwdata ||
	    !mcsigallwdataModel_rv ||
	    !mcsigallwdataModel_wv ||
	    !mcsigallwdataModel  ) {
	  std::cerr<<" ERROR: Could not create some Dataset from ROOT file."<<std::endl;
	  return;
	}
	  

	//set higgs mass and fit range
	mnom.setVal(mhs.at(i));
	hmass2->setRange("plotrange",mhs.at(i)-15.,mhs.at(i)+15.);
	
	// reset the fir parameters to the starting values
	if( ! resetStartValues(fitparms[iCat],startvals[iProc*numCats+iCat]) ) {
	  std::cerr<<" ERROR: Could not reset parameters to starting values."<<std::endl;
	  return;
	}
	
	//fit
	RooFitResult *fitres_rv = 0;
	RooFitResult *fitres_wv = 0;
	
	RooAbsPdf* rightpdf = combh_rv[iCat];
	RooAbsPdf* wrongpdf = combh_wv[iCat];


	// perfomr thr fit only if the mass point is a fit mass guy...

	if ( mhs_fit[i] == 0 ) {
	  printf("FIT SKIPPING:  proc:%d mass point:%d cat:%d\n",iProc,(int) mhs.at(i),iCat); 

	  double eaccnumModel = ( isDataSetEmpty ? 0. : mcsigallwdataModel_rv->sumEntries()+mcsigallwdataModel_wv->sumEntries());
	  double eaccden = ntots[i];
	  double eacc = eaccnumModel/eaccden;
	  //double eaccerrlo = TEfficiency::ClopperPearson(Int_t(eaccden), Int_t(eaccnumModel), 0.683, kFALSE) - eacc;
	  //double eaccerrhi = TEfficiency::ClopperPearson(Int_t(eaccden), Int_t(eaccnumModel), 0.683, kTRUE) - eacc;
	  printf("eacc = %5f, eaccnum = %5f, eaccden = %5f\n",eacc,eaccnumModel,eaccden);

	  // keep track....
	  effaccMap.find(procname)->second->find(catnamesbase[iCat])->second->insert(massValPair(mhs[i],eacc));

	} else {
	  fitMass++;
	  printf("FITSTART:  proc:%d mass point:%d cat:%d\n",iProc,(int) mhs.at(i),iCat); 
	  
	  hmass2->setRange("fitrange",mhs.at(i)-fitmassmin,mhs.at(i)+fitmassmax);
	  
	  fitres_rv      = rightpdf     ->fitTo(*mcsigallwdataModel_rv,     Strategy(0),Minimizer("Minuit2",""),Minos(kFALSE),SumW2Error(kFALSE), Save(kTRUE),NumCPU(8),RooFit::Range("fitrange"));
	  printf("FITRES:  proc:%d mass point:%d cat:%d fit status right vertex:%d \n",iProc,(int) mhs.at(i),iCat,fitres_rv->status()); 
	  
	  if( !fitres_rv ) {
	    std::cerr<<" ERROR: Did not get any fitresult."<<std::endl;
	    return;
	  }
	  
	  fitres_wv      = wrongpdf     ->fitTo(*mcsigallwdataModel_wv,     Strategy(0),Minimizer("Minuit2",""),Minos(kFALSE),SumW2Error(kFALSE), Save(kTRUE),NumCPU(8),RooFit::Range("fitrange"));
	  printf("FITRES:  proc:%d mass point:%d cat:%d fit status right vertex:%d \n",iProc,(int) mhs.at(i),iCat,fitres_wv->status()); 
	  
	  if( !fitres_wv ) {
	    std::cerr<<" ERROR: Did not get any fitresult."<<std::endl;
	    return;
	  }
	  
	  
	  // ============================================
	  for( int iVtx =0; iVtx < 2; ++iVtx ) {
	    
	    TString vtxString = ( iVtx ? "wv" : "rv" );
	    
	    // sort the Gaussians by POSITION !
	    std::vector< std::map<TString,RooRealVar*>::iterator > allGs;
	    std::vector< std::map<TString,RooRealVar*>::iterator > allGs_s;
	    std::vector< std::map<TString,RooRealVar*>::iterator > allGs_f;
	    std::vector< double                                  > allMeans;
	    
	    /// ======== need to compute the correct coefficents for recurive adding
	    int gCounter=1;
	    std::map<TString,double> realCoefs;
	    std::map<TString,RooRealVar*>::iterator gIt_f = fitparms[iCat].find( TString("f1_") + vtxString+TString("_")+ procname +TString("_")+catnamesbase[iCat]);
	    double prodLastCoefs = 1.;
	    while ( gIt_f != fitparms[iCat].end() ) {
	      double thisCoef = gIt_f->second->getVal();
	      double realCoef = prodLastCoefs*thisCoef;
	      realCoefs.insert(std::pair<TString,double>(gIt_f->first,realCoef));
	      prodLastCoefs = prodLastCoefs*(1-thisCoef);
	      gCounter++;
	      gIt_f = fitparms[iCat].find( TString::Format("f%d_",gCounter) + vtxString+TString("_")+ procname +TString("_")+catnamesbase[iCat]);
	    }
	    
	    std::map<TString,RooRealVar*>::iterator gIt   = fitparms[iCat].find( TString("dmG1_") + vtxString+TString("_") + procname +TString("_")+catnamesbase[iCat]);
	    std::map<TString,RooRealVar*>::iterator gIt_s = fitparms[iCat].find( TString("sigmaG1_") + vtxString+TString("_")+ procname +TString("_")+catnamesbase[iCat]);
	    gIt_f = fitparms[iCat].find( TString("f1_") + vtxString+TString("_")+ procname +TString("_")+catnamesbase[iCat]);
	    gCounter=1;
	    // if there's no 'f1' must be one component Gaussian, we can skip all the sorting...
	    if ( gIt_f == fitparms[iCat].end() ) {
	      // BUT: Set the sigma to the absolute vakue never-the-less
	      gIt_s->second->setVal( TMath::Abs(gIt_s->second->getVal()) );
	    
	      // otherwise we sort...
	    } else {
	      while ( gIt != fitparms[iCat].end() ) {
		allGs   .push_back(gIt  );
		allGs_s .push_back(gIt_s);
		allGs_f .push_back(gIt_f);  // careful! The last componeten push_es_back ( fitparms[iCat].end() ) !!!
		//allMeans.push_back( (gIt->second)->getVal() );
		allMeans.push_back( TMath::Abs((gIt_s->second)->getVal()) );
		gCounter++;
		gIt   = fitparms[iCat].find( TString::Format("dmG%d_",gCounter) + vtxString+TString("_")+ procname +TString("_")+catnamesbase[iCat]);
		gIt_s = fitparms[iCat].find( TString::Format("sigmaG%d_",gCounter) + vtxString+TString("_")+ procname +TString("_")+catnamesbase[iCat]);
		gIt_f = fitparms[iCat].find( TString::Format("f%d_",gCounter) + vtxString+TString("_")+ procname +TString("_")+catnamesbase[iCat]);
	      }
	      
	      // loop over them and sort
	      std::vector<double> sortMean;
	      std::vector<double> sortSigma;
	      std::vector<double> sortFrac;
	      std::vector<double> realFrac;
	      
	      //double sumFrac = 0.;
	      
	      while ( allGs.size() > 0 ) {
		
		double       minVal = 100000.;
		unsigned int minPos = 0;   
		std::vector< std::map<TString,RooRealVar*>::iterator >::iterator 
		  minG   = allGs.begin();   
		unsigned int iG     = 0;
		for( std::vector< std::map<TString,RooRealVar*>::iterator >::iterator t_it = allGs.begin(); t_it != allGs.end(); ++t_it ) {
		  if (allMeans[iG] < minVal) {
		    minG   = t_it;
		    minPos = iG;
		    minVal = allMeans[iG];
		  }
		  iG++;
		}
		sortMean .push_back( (allGs[minPos]->second)->getVal() );
		sortSigma.push_back(  minVal );
		//sortMean .push_back(  minVal );
		//sortSigma.push_back( (allGs_s[minPos]->second)->getVal() );
		// !!! CAUTION: fro the last component (careful, does not have to be last in map!), push_back -1
		if ( allGs_f[minPos] == fitparms[iCat].end() ) {
		  sortFrac .push_back( prodLastCoefs );		
		  //sortFrac .push_back( -1. );		
		} else {
		  sortFrac .push_back( realCoefs.find(allGs_f[minPos]->first)->second );
		  //sortFrac .push_back( (allGs_f[minPos]->second)->getVal() );
		  //sumFrac += (allGs_f[minPos]->second)->getVal();
		}
		
		std::vector< std::map<TString,RooRealVar*>::iterator >::iterator it_s    = allGs_s.begin();
		std::vector< std::map<TString,RooRealVar*>::iterator >::iterator it_f    = allGs_f.begin();
		std::vector< double >::iterator                                  it_mean = allMeans.begin();
		
		for(unsigned int iPos=0; iPos<minPos; ++iPos) {
		  it_s++;
		  it_f++;
		  it_mean++;
		}
		
		allGs   .erase( minG );
		allGs_s .erase( it_s );
		allGs_f .erase( it_f );
		allMeans.erase( it_mean );
	      }
	      
	      
	      // now we have the sorted values, assign them back to the RooRealVars
	      gCounter=0;
	      gIt   = fitparms[iCat].find( TString("dmG1_") + vtxString+TString("_")+ procname +TString("_")+catnamesbase[iCat]);
	      gIt_s = fitparms[iCat].find( TString("sigmaG1_") + vtxString+TString("_")+ procname +TString("_")+catnamesbase[iCat]);
	      gIt_f = fitparms[iCat].find( TString("f1_") + vtxString+TString("_")+ procname +TString("_")+catnamesbase[iCat]);
	      
	      prodLastCoefs = 1.;
	      while ( gIt != fitparms[iCat].end() ) {
		(gIt->second)  ->setVal( sortMean[gCounter] );
		(gIt_s->second)->setVal( TMath::Abs(sortSigma[gCounter]) );
		if( gIt_f != fitparms[iCat].end() ) {
		  double thisCoef = sortFrac[gCounter]/prodLastCoefs;
		  gIt_f->second->setVal( thisCoef );
		  prodLastCoefs *= (1-thisCoef);
		  // 		  if ( sortFrac[gCounter] > 0. )
		  // 		    (gIt_f->second)->setVal( sortFrac[gCounter] );
		  // 		  else
		  // 		    (gIt_f->second)->setVal( 1. - sumFrac );   // this guy was the last in the list and had no fraction, so assign 1.-Sum(fractions)
		  
		}
		gCounter++;
		gIt   = fitparms[iCat].find( TString::Format("dmG%d_",gCounter+1) + vtxString+TString("_")+ procname +TString("_")+catnamesbase[iCat]);
		gIt_s = fitparms[iCat].find( TString::Format("sigmaG%d_",gCounter+1) + vtxString+TString("_")+ procname +TString("_")+catnamesbase[iCat]);
		gIt_f = fitparms[iCat].find( TString::Format("f%d_",gCounter+1) + vtxString+TString("_")+ procname +TString("_")+catnamesbase[iCat]);
	      }
	    }
	  }
	
	  double eaccnum = ( isDataSetEmpty ? 0. : mcsigallwdata_rv->sumEntries()+mcsigallwdata_wv->sumEntries() );
	  double eaccden = ntots[i];
	  double eacc = eaccnum/eaccden;
	  //printf("eacc = %5f, eaccnum = %5f, eaccden = %5f\n",eacc,eaccnum,eaccden);
	  double eaccerrlo = TEfficiency::ClopperPearson(Int_t(eaccden), Int_t(eaccnum), 0.683, kFALSE) - eacc;
	  double eaccerrhi = TEfficiency::ClopperPearson(Int_t(eaccden), Int_t(eaccnum), 0.683, kTRUE) - eacc;
	  
	  effacc[iCat]->setVal(eacc);
	  effacc[iCat]->setAsymError(eaccerrlo,eaccerrhi);
	  
	  
	  //compute the fraction of right vertex
	  //double fright = mcsigwdata->sumEntries()/(mcsigwdata->sumEntries()+mcsigwrongwdata->sumEntries());
	  double frightnum = mcsigallwdata_rv->sumEntries();
	  double frightden = mcsigallwdata_rv->sumEntries()+mcsigallwdata_wv->sumEntries();
	  double fright = frightnum/frightden;
	  double frighterrlo = TEfficiency::ClopperPearson(Int_t(frightden), Int_t(frightnum), 0.683, kFALSE) - fright;
	  double frighterrhi = TEfficiency::ClopperPearson(Int_t(frightden), Int_t(frightnum), 0.683, kTRUE) - fright;
	  
	  fracright[iCat]->setVal(fright);
	  fracright[iCat]->setAsymError(frighterrlo,frighterrhi);
	  
	  // ================== BEGIN effective smearing =========================
	  double eaccnumModel = ( isDataSetEmpty ? 0. : mcsigallwdataModel_rv->sumEntries()+mcsigallwdataModel_wv->sumEntries());
	  eacc = eaccnumModel/eaccden;
	  eaccerrlo = TEfficiency::ClopperPearson(Int_t(eaccden), Int_t(eaccnumModel), 0.683, kFALSE) - eacc;
	  eaccerrhi = TEfficiency::ClopperPearson(Int_t(eaccden), Int_t(eaccnumModel), 0.683, kTRUE) - eacc;
	  printf("eacc = %5f, eaccnum = %5f, eaccden = %5f\n",eacc,eaccnumModel,eaccden);
	  effaccModel[iCat]->setVal(eacc);
	  effaccModel[iCat]->setAsymError(eaccerrlo,eaccerrhi);
	  
	  // keep track....
	  effaccMap.find(procname)->second->find(catnamesbase[iCat])->second->insert(massValPair(mhs[i],eacc));


	  double frightnumModel = mcsigallwdataModel_rv->sumEntries();
	  frightden = mcsigallwdataModel_rv->sumEntries()+mcsigallwdataModel_wv->sumEntries();
	  fright = frightnumModel/frightden;
	  printf("fright = %5f, frightnum = %5f, frightden = %5f\n",fright,frightnumModel,frightden);
	  frighterrlo = TEfficiency::ClopperPearson(Int_t(frightden), Int_t(frightnumModel), 0.683, kFALSE) - fright;
	  frighterrhi = TEfficiency::ClopperPearson(Int_t(frightden), Int_t(frightnumModel), 0.683, kTRUE) - fright;
	  fracrightModel[iCat]->setVal(fright);
	  fracrightModel[iCat]->setAsymError(frighterrlo,frighterrhi);
	  // ================== END effective smearing =========================
	  printf("eff done\n");
	  
	  for( int iVtx = 0; iVtx < 2; ++iVtx ) {
	    TCanvas *chfit = new TCanvas;
	    TString plotname    = catnames.at(iCat) + (iVtx ? TString("_wv_") : TString("_rv_")) + numstring + TString(".pdf");//e.g. rightvtx110cat0_ggh.pdf      
	    TString plotnameLog = catnames.at(iCat) + (iVtx ? TString("_wv_") : TString("_rv_")) + numstring + TString("_log.pdf");//e.g. rightvtx110cat0_ggh.pdf      
	    RooPlot *hplot = hmass2->frame(Bins(60),Range("plotrange"));
	    
	    RooDataSet* theData = ( iVtx ? mcsigallwdataModel_wv : mcsigallwdataModel_rv );
	    RooAbsPdf* theModel = ( iVtx ? wrongpdf : rightpdf );
	    
	    theData->plotOn(hplot,RooFit::MarkerStyle(25),RooFit::MarkerSize(1.2));	
	    theModel->plotOn(hplot,RooFit::LineColor(kBlue),Range("plotrange"),NormRange("fitrange"));  
	    hplot->SetTitle("");
	    hplot->Draw();
	    if( !fitonly ) {
	      chfit->SaveAs( TString::Format("%s/%s/%s",projectDir.Data(),procname.Data(),plotname.Data() ));
	    }
	    
	    if( !iVtx ) {
	      rightCan->cd(iCat*numMasses+fitMass+1);
	      hplot->Draw(); 
	    } else {
	      wrongCan->cd(iCat*numMasses+fitMass+1);
	      hplot->Draw(); 
	    }
	    
	    TString plotnameallOverview = TString("allvtx_") + (iVtx ? TString("wv") : TString("rv") ) +TString(".pdf");  
	    if ( !fitonly  ) {
	      if( !iVtx )
		rightCan->SaveAs( TString::Format("%s/%s/%s",projectDir.Data(),procname.Data(),plotnameallOverview.Data() ) ); 
	      else
		wrongCan->SaveAs( TString::Format("%s/%s/%s",projectDir.Data(),procname.Data(),plotnameallOverview.Data() ) ); 
	    }
	  }
	  
	  printf ("Sum of weights: %5f \n",mcsigallwdata->sumEntries());  
	  
	  unsigned int iParm=0;
	  for(std::map<TString,RooRealVar*>::iterator pair = fitparms[iCat].begin(); pair != fitparms[iCat].end(); ++pair) {
	    fitparmhists[iCat][iParm] -> Fill(mhs.at(i),pair->second->getVal());
	    fitparmhists[iCat][iParm] -> SetBinError(fitparmhists[iCat][iParm]->FindFixBin(mhs.at(i)),pair->second->getError());
	    if (pair->second->hasAsymError()) {
	      double avgerror = (pair->second->getErrorLo() + pair->second->getErrorHi())/2.0;
	      fitparmhists[iCat][iParm]->SetBinError(fitparmhists[iCat][iParm]->FindFixBin(mhs.at(i)),avgerror); 
	    }
	    iParm++;
	  }
	  histfitstatus_rv->Fill(mhs.at(i),iCat,fitres_rv->status());
	  histfitstatus_wv->Fill(mhs.at(i),iCat,fitres_wv->status());
	  
	  dummy->cd();
	}

	// ================== BEGIN effective smearing =========================
	// need to know the cross-section x BR for normalization	
	RooAbsReal* theProcCSList = procxseclist[iProc];
	if ( histFuncMap.find("smbr")== histFuncMap.end()) {
	  std::cout<<" can't find smbr "<<std::endl;
	  return;
	}
	
	RooAbsReal* theBR = histFuncMap.find("smbr")->second;
	mnom.setVal(mhs[i]);
	
	dataAllName=TString::Format("sig_%s_massModelNom_m",procname.Data())+numstring+TString("_")+catnamesbase.at(iCat);
	RooDataSet *mcsigallwdataModel_effSmear = NULL;
	bool dummyBool = true;
	mcsigallwdataModel_effSmear   = makedset(dataAllName,   trees[i],   theWeight*(masscutModel && catIt->second), hmass2, weight, mhs[i], dummyBool, theProcCSList->getVal()*theBR->getVal()/ntots[i]*lumi*1000.);
	  
	// =========================================
	// we got it, push it back for later reference, using a map by name FIX-ME: Should use nonModel for eff x acc numbers at the end... (does not matter for Hgg)
	keepAllDataMap.insert(std::pair<TString,RooDataSet*> ( dataAllName, mcsigallwdataModel_effSmear) );
	
	RooDataSet *mcsigallwdataModel_noSmear    = NULL;
	if(doeffsmear) {
	  dataAllName=TString::Format("sig_%s_massModelNoSmear_m",procname.Data())+numstring+TString("_")+catnamesbase.at(iCat);
	  
	  mcsigallwdataModel_noSmear    = makedset(dataAllName,   treesNS[i], theWeight*(masscutModel && catIt->second), hmass2, weight, mhs[i], dummyBool, theProcCSList->getVal()*theBR->getVal()/ntots[i]*lumi*1000.,effSmearHistsNS[i][iCat]);
	}
	
	// append this DataSet to the per-mass/per-cat datasets
	effSmearData   [iCat][i]->append(*mcsigallwdataModel_effSmear);	
	// add all processes & Cats to gether
	effSmearData[numCats][i]->append(*mcsigallwdataModel_effSmear);	
	
	if( doeffsmear )
	  effSmearDataNS [iCat][i]->append(*mcsigallwdataModel_noSmear);	
	
	// ================== END effective smearing =========================
	std::cout<<"    Nev = "<<mcsigallwdataModel_effSmear->sumEntries()<<"   ( Ntot = "<<ntots[i]<<"  L = "<<1000.*lumi<<"  CS = "<<theProcCSList->getVal()<<"  BR = "<<theBR->getVal()<<" )"<<std::endl;
      }
    }

    if( fitonly ) return;

    // ---------------------------------------------------------------
    // prepare RooDataSet and RooHistFunc arrays
    RooDataHist*** fitparmdatas    = new RooDataHist**[numCats];
    RooDataHist*** fitparmdatas2nd = new RooDataHist**[numCats];
    
    std::map<TString,RooHistFunc*>* fitparmfuncs    = new std::map<TString,RooHistFunc*>[numCats];
    std::map<TString,RooHistFunc*>* fitparmfuncs2nd = new std::map<TString,RooHistFunc*>[numCats];
    
    for(int iCat=0; iCat<numCats; ++iCat) {

      fitparmdatas[iCat] = new RooDataHist*[fitparms[iCat].size()];
      
      if( doSecondSignal )
	fitparmdatas2nd[iCat] = new RooDataHist*[fitparms[iCat].size()];

      // loop over all fitparams
      // -----------------------------------------
      unsigned int iParm=0;
      for(std::map<TString,RooRealVar*>::iterator pair = fitparms[iCat].begin(); pair != fitparms[iCat].end(); ++pair) {
	TString dataname=TString("data") + pair->first;// + catnames.at(iCat);
	TString funcname=TString("func") + pair->first;// + catnames.at(iCat);
	fitparmdatas[iCat][iParm] = new RooDataHist(dataname,dataname,RooArgList(mnom),fitparmhists[iCat][iParm]);
	fitparmfuncs[iCat].insert(std::pair<TString,RooHistFunc*>(pair->first, new RooHistFunc(funcname,funcname,RooArgList(mnom),*fitparmdatas[iCat][iParm],1)));
	
	if( doSecondSignal ) {	  
	  TString dataname2nd=TString("data2nd") + pair->first;// + catnames.at(iCat);
	  TString funcname2nd=TString("func2nd") + pair->first;// + catnames.at(iCat);
	  fitparmdatas2nd[iCat][iParm] = new RooDataHist(dataname2nd,dataname2nd,RooArgList(mnom2),fitparmhists[iCat][iParm]);
	  fitparmfuncs2nd[iCat].insert(std::pair<TString,RooHistFunc*>(pair->first, new RooHistFunc(funcname2nd,funcname2nd,RooArgList(mnom2),*fitparmdatas2nd[iCat][iParm],1)));
	}

	// do the plottiong here if requests
	if( !fitonly ) {
	  TString plotname = TString("func_") + pair->first + TString(".pdf");//e.g funcsigma1cat0_ggh.pdf
	  TCanvas *cfunctest = new TCanvas;
	  //  fitparmhists[iCat*nparms + iparm]->Draw();
	  RooPlot *hploteffacc = mnom.frame(Bins(100),Range(mhs[0]-5.,mhs[mhs.size()-1]+5.));
	  fitparmdatas[iCat][iParm]->plotOn(hploteffacc);  
	  //  fitparmfuncs[iCat][iparm]->plotOn(hploteffacc,RooFit::LineColor(kBlue));  
	  fitparmfuncs[iCat].find(pair->first)->second->plotOn(hploteffacc,RooFit::LineColor(kBlue));  
	  hploteffacc->SetTitle("");
	  hploteffacc->GetYaxis()->SetTitle(pair->first.Data());
	  hploteffacc->Draw(); 

	  cfunctest->SaveAs(TString::Format("%s/%s/%s",projectDir.Data(),procname.Data(),plotname.Data() ) );
	  delete cfunctest;
	}
	iParm++;

      }
    }

    // ---------------------------------------------------------------
    //---draw fit status---
    TCanvas *cfitstatus_rv = new TCanvas;
    histfitstatus_rv->Draw("COL");
    if( !fitonly )
      cfitstatus_rv->SaveAs( TString::Format("%s/%s/%s",projectDir.Data(),procname.Data(),"fitstatus_rv.pdf") );

    TCanvas *cfitstatus_wv = new TCanvas;
    histfitstatus_wv->Draw("COL");
    if( !fitonly )
      cfitstatus_wv->SaveAs( TString::Format("%s/%s/%s",projectDir.Data(),procname.Data(),"fitstatus_wv.pdf") );
    
    //---define final pdfs in each category---
   //---nuissance---

    //the nuissance and the systematics are redone in combine
    RooRealVar nuissancedeltafracright("CMS_hgg_nuissancedeltafracright","",1.0,0.1,10.0);
    nuissancedeltafracright.setConstant();

    RooConstVar   **smears                   = new RooConstVar*[catnames.size()];
    RooRealVar    **nuissancedeltasmears     = new RooRealVar*[catnames.size()];
    RooRealVar    **nuissancedeltams     = new RooRealVar*[catnames.size()];
    RooFormulaVar **smearmods            = new RooFormulaVar*[catnames.size()];

    RooFormulaVar*** meanslidesG_rv  = new RooFormulaVar**[catnames.size()];
    RooFormulaVar*** sigmaslidesG_rv  = new RooFormulaVar**[catnames.size()];
    RooFormulaVar*** fracslides_rv   = new RooFormulaVar**[catnames.size()];

    RooFormulaVar*** meanslidesG_wv  = new RooFormulaVar**[catnames.size()];
    RooFormulaVar*** sigmaslidesG_wv  = new RooFormulaVar**[catnames.size()];
    RooFormulaVar*** fracslides_wv   = new RooFormulaVar**[catnames.size()];


    RooFormulaVar*** meanslidesG_rv_2nd  = new RooFormulaVar**[catnames.size()];
    RooFormulaVar*** sigmaslidesG_rv_2nd  = new RooFormulaVar**[catnames.size()];
    RooFormulaVar*** fracslides_rv_2nd   = new RooFormulaVar**[catnames.size()];

    RooFormulaVar*** meanslidesG_wv_2nd  = new RooFormulaVar**[catnames.size()];
    RooFormulaVar*** sigmaslidesG_wv_2nd  = new RooFormulaVar**[catnames.size()];
    RooFormulaVar*** fracslides_wv_2nd   = new RooFormulaVar**[catnames.size()];


    RooGaussian*** gslides_rv  = new RooGaussian**[catnames.size()];
    RooGaussian*** gslides_rv_2nd  = new RooGaussian**[catnames.size()];
    RooGaussian*** gslidesModel_rv  = new RooGaussian**[catnames.size()];

    RooGaussian*** gslides_wv  = new RooGaussian**[catnames.size()];
    RooGaussian*** gslides_wv_2nd  = new RooGaussian**[catnames.size()];
    RooGaussian*** gslidesModel_wv  = new RooGaussian**[catnames.size()];


    RooAddPdf**   combhslides_rv      = new RooAddPdf*[catnames.size()];    
    RooAddPdf**   combhslides_rv_2nd      = new RooAddPdf*[catnames.size()];    
    RooAddPdf**   combhslidesModel_rv      = new RooAddPdf*[catnames.size()];  

    RooAddPdf**   combhslides_wv      = new RooAddPdf*[catnames.size()];    
    RooAddPdf**   combhslides_wv_2nd      = new RooAddPdf*[catnames.size()];    
    RooAddPdf**   combhslidesModel_wv      = new RooAddPdf*[catnames.size()];  

    RooFormulaVar** fracrightmodslides      = new RooFormulaVar*[catnames.size()];  
    RooFormulaVar** fracrightmodslides_2nd      = new RooFormulaVar*[catnames.size()];  
    RooFormulaVar** fracrightmodslidesModel = new RooFormulaVar*[catnames.size()];  


    RooAddPdf  **combhvtxslides      = new RooAddPdf*[catnames.size()];
    RooAddPdf  **combhvtxslides_2nd      = new RooAddPdf*[catnames.size()];
    RooAddPdf  **combhvtxslidesModel = new RooAddPdf*[catnames.size()];

    RooAbsPdf  **finalpdfslides      = new RooAbsPdf*[catnames.size()];
    RooAbsPdf  **finalpdfslides_2nd      = new RooAbsPdf*[catnames.size()];
    RooAbsPdf  **finalpdfslidesModel = new RooAbsPdf*[catnames.size()];

    std::vector<RooAbsReal*> finalnormslides;
    finalnormslides.resize(catnames.size());

    std::vector<RooAbsReal*> finalnormslides_2nd;
    finalnormslides_2nd.resize(catnames.size());


    // ======== BEGIN ==============
    std::vector<RooAbsReal*> finalnormslidesModel;
    finalnormslidesModel.resize(catnames.size());
    // =============================

    RooRealVar* globalScale = new RooRealVar("CMS_hgg_globalscale","",0.0,-0.004717,0.004717);
    globalScale->setConstant();


    for(int iCat=0; iCat < numCats; ++iCat) {
      
      smears[iCat] = new RooConstVar(TString("smear_")+catnames.at(iCat),"",smearingv.at(iCat));
      
      //nuissancedeltasmears[iCat] = new RooRealVar(TString("CMS_hgg_nuissancedeltasmear")+catnamesbase.at(iCat),"",0.0, -smearingv.at(iCat),smearingv.at(iCat));
      nuissancedeltasmears[iCat] = new RooRealVar(TString::Format("CMS_hgg_nuissancedeltasmear%s",(catToSmearMap.find(catnamesbase.at(iCat))->second).Data()),"",0.0, -smearingv.at(iCat),smearingv.at(iCat));
      nuissancedeltasmears[iCat]->setConstant();
      //nuissancedeltams[iCat] = new RooRealVar(TString("CMS_hgg_nuissancedeltam")+catnamesbase.at(iCat),"",0.0,-5.0,5.0);
      nuissancedeltams[iCat] = new RooRealVar(TString::Format("CMS_hgg_nuissancedeltam%s",(catToSmearMap.find(catnamesbase.at(iCat))->second).Data()),"",0.0,-5.0,5.0);
      nuissancedeltams[iCat]->setConstant();
					      
      smearmods[iCat] = new RooFormulaVar(TString("smearmod")+catnames.at(iCat),"","@0*(@1 + @2)",RooArgList(mnom,*smears[iCat],*nuissancedeltasmears[iCat]));

      // ------------------------------------------------------------------------------------------------------------------
      // Right & Wrong Vertex Models
      for( int iVtx = 0; iVtx < 2; ++iVtx ) {

	TString vtxString = ( iVtx ? "wv" : "rv" );
	
	RooArgList compListModel;
	RooArgList compList;
	RooArgList fracList;

	RooArgList compList_2nd;
	RooArgList fracList_2nd;

	// first counter how many guys we got (little ugly....)
	int gCounter=1;
	std::map<TString,RooRealVar*>::iterator gIt   = fitparms[iCat].find( TString("dmG1_") + vtxString +TString("_") + procname +TString("_")+catnamesbase[iCat]);
	std::map<TString,RooRealVar*>::iterator gIt_s = fitparms[iCat].find( TString("sigmaG1_") + vtxString +TString("_") + procname +TString("_")+catnamesbase[iCat]);
	std::map<TString,RooRealVar*>::iterator gIt_f = fitparms[iCat].find( TString("f1_") + vtxString +TString("_") + procname +TString("_")+catnamesbase[iCat]);

	while ( gIt != fitparms[iCat].end() ) {
	  gCounter++;
	  gIt   = fitparms[iCat].find( TString::Format("dmG%d_",gCounter) + vtxString +TString("_") + procname +TString("_")+catnamesbase[iCat]);
	  gIt_s = fitparms[iCat].find( TString::Format("sigmaG%d_",gCounter) + vtxString +TString("_") + procname +TString("_")+catnamesbase[iCat]);
	  gIt_f = fitparms[iCat].find( TString::Format("f%d_",gCounter) + vtxString +TString("_") + procname +TString("_")+catnamesbase[iCat]);
	}
	


	// create the correct size variable arays
	if ( iVtx ) {
	  meanslidesG_wv[iCat]  = new RooFormulaVar*[gCounter-1];
	  sigmaslidesG_wv[iCat] = new RooFormulaVar*[gCounter-1];
	  fracslides_wv[iCat]   = new RooFormulaVar*[gCounter-1];
	  gslides_wv[iCat]      = new RooGaussian*[gCounter-1];


	  if (doSecondSignal) {
	    meanslidesG_wv_2nd[iCat]  = new RooFormulaVar*[gCounter-1];
	    sigmaslidesG_wv_2nd[iCat] = new RooFormulaVar*[gCounter-1];
	    fracslides_wv_2nd[iCat]   = new RooFormulaVar*[gCounter-1];
	    gslides_wv_2nd[iCat]      = new RooGaussian*[gCounter-1];
	  }

	  // ============= BEGIN effective smearing ==============
	  gslidesModel_wv[iCat]      = new RooGaussian*[gCounter-1];
	  // =====================================================
	} else {
	  meanslidesG_rv[iCat]  = new RooFormulaVar*[gCounter-1];
	  sigmaslidesG_rv[iCat] = new RooFormulaVar*[gCounter-1];
	  fracslides_rv[iCat]   = new RooFormulaVar*[gCounter-1];
	  gslides_rv[iCat]      = new RooGaussian*[gCounter-1];

	  if (doSecondSignal) {
	    meanslidesG_rv_2nd[iCat]  = new RooFormulaVar*[gCounter-1];
	    sigmaslidesG_rv_2nd[iCat] = new RooFormulaVar*[gCounter-1];
	    fracslides_rv_2nd[iCat]   = new RooFormulaVar*[gCounter-1];
	    gslides_rv_2nd[iCat]      = new RooGaussian*[gCounter-1];
	  }

	  // ============= BEGIN effective smearing ==============
	  gslidesModel_rv[iCat]      = new RooGaussian*[gCounter-1];
	  // =====================================================
	}

	// loop over all gaussians and produce 'slides'
	gCounter=1;
	gIt   = fitparms[iCat].find( TString("dmG1_") + vtxString +TString("_") + procname +TString("_")+catnamesbase[iCat]);
	
	TString slideName;
	TString parName;
	std::map<TString,RooHistFunc*>::iterator pair;
	
	while ( gIt != fitparms[iCat].end() ) {
	  
	  slideName=TString::Format("meanslideG%d%s%s",gCounter,catnames.at(iCat).Data(),vtxString.Data());
	  parName  =TString::Format("dmG%d_%s_%s_%s",gCounter,vtxString.Data(),procname.Data(),catnamesbase[iCat].Data());
	  pair = fitparmfuncs[iCat].find( parName );
	  if( pair == fitparmfuncs[iCat].end() ) {
	    std::cerr<<" ERROR: Could not find RooHistFunction for parameter "<<parName<<"."<<std::endl;
	    return;
	  }
	  if ( iVtx )
	    meanslidesG_wv[iCat][gCounter-1] = new RooFormulaVar(slideName.Data(),"","@0 + @1 +@0*@2+@0*@3",RooArgList(mnom,*(pair->second),*nuissancedeltams[iCat],*globalScale));
	  else
	    meanslidesG_rv[iCat][gCounter-1] = new RooFormulaVar(slideName.Data(),"","@0 + @1 +@0*@2+@0*@3",RooArgList(mnom,*(pair->second),*nuissancedeltams[iCat],*globalScale));	  
	  
	  slideName=TString::Format("sigmaslideG%d%s%s",gCounter,catnames.at(iCat).Data(),vtxString.Data());
	  parName  =TString::Format("sigmaG%d_%s_%s_%s",gCounter,vtxString.Data(),procname.Data(),catnamesbase[iCat].Data());
	  pair = fitparmfuncs[iCat].find( parName );
	  if( pair == fitparmfuncs[iCat].end() ) {
	    std::cerr<<" ERROR: Could not find RooHistFunction for parameter "<<parName<<"."<<std::endl;
	    return;
	  }
	  
	  if ( iVtx ) 
	    sigmaslidesG_wv[iCat][gCounter-1] = new RooFormulaVar(slideName.Data(),"","TMath::Max(0.01,sqrt(@0*@0-@3*@3*@2*@2 +@1*@1))",RooArgList(*(pair->second),*smearmods[iCat],*smears[iCat],mnom));//why set min 0.01?	
	  else
	    sigmaslidesG_rv[iCat][gCounter-1] = new RooFormulaVar(slideName.Data(),"","TMath::Max(0.01,sqrt(@0*@0-@3*@3*@2*@2 +@1*@1))",RooArgList(*(pair->second),*smearmods[iCat],*smears[iCat],mnom));//why set min 0.01?	
	  
	  slideName=TString::Format("fracslide%d%s%s",gCounter,catnames.at(iCat).Data(),vtxString.Data());
	  parName  =TString::Format("f%d_%s_%s_%s",gCounter,vtxString.Data(),procname.Data(),catnamesbase[iCat].Data());
	  pair = fitparmfuncs[iCat].find( parName );
	  if( pair == fitparmfuncs[iCat].end() ) {
	    std::cerr<<" WARNING: Could not find RooHistFunction for parameter "<<parName<<"."<<std::endl;
	    std::cout<<"          Must be last component."<<std::endl;
	  } else {
	    if( iVtx ) {
	      fracslides_wv[iCat][gCounter-1] = new RooFormulaVar(slideName.Data(),"","@0",RooArgList(*(pair->second)));      	  
	      fracList.add(* (fracslides_wv[iCat][gCounter-1]) );
	    } else {
	      fracslides_rv[iCat][gCounter-1] = new RooFormulaVar(slideName.Data(),"","@0",RooArgList(*(pair->second)));      	  
	      fracList.add(* (fracslides_rv[iCat][gCounter-1]) );
	    }
	  }
	  
	  slideName=TString::Format("gslide%d%s%s",gCounter,catnames.at(iCat).Data(),vtxString.Data());
	  if( iVtx )
	    gslides_wv[iCat][gCounter-1] = new RooGaussian(slideName.Data(),"",*hmass,*(meanslidesG_wv[iCat][gCounter-1]),*(sigmaslidesG_wv[iCat][gCounter-1]));
	  else
	    gslides_rv[iCat][gCounter-1] = new RooGaussian(slideName.Data(),"",*hmass,*(meanslidesG_rv[iCat][gCounter-1]),*(sigmaslidesG_rv[iCat][gCounter-1]));
	  
	  // ============== BEGIN effective smearing ======================
	  slideName=TString::Format("gslideModel%d%s%s",gCounter,catnames.at(iCat).Data(),vtxString.Data());
	  if( iVtx ) {
	    gslidesModel_wv[iCat][gCounter-1] = new RooGaussian(slideName.Data(),"",*hmass2,*(meanslidesG_wv[iCat][gCounter-1]),*(sigmaslidesG_wv[iCat][gCounter-1]));
	    compListModel.add(* (gslidesModel_wv[iCat][gCounter-1]) );
	  } else {
	    gslidesModel_rv[iCat][gCounter-1] = new RooGaussian(slideName.Data(),"",*hmass2,*(meanslidesG_rv[iCat][gCounter-1]),*(sigmaslidesG_rv[iCat][gCounter-1]));
	    compListModel.add(* (gslidesModel_rv[iCat][gCounter-1]) );
	  }
	  
	  // ==============================================================
	  
	  compList.add( (iVtx ? ( *(gslides_wv[iCat][gCounter-1]) ): ( *(gslides_rv[iCat][gCounter-1]) )));
	  
	  // ========================================================================
	  // models fro second Higgs as BG
	  if( doSecondSignal ) {
	    
	    slideName=TString::Format("meanslideG%d%s%s_2nd",gCounter,catnames.at(iCat).Data(),vtxString.Data());
	    parName  =TString::Format("dmG%d_%s_%s_%s",gCounter,vtxString.Data(),procname.Data(),catnamesbase[iCat].Data());
	    pair = fitparmfuncs2nd[iCat].find( parName );
	    if( pair == fitparmfuncs2nd[iCat].end() ) {
	      std::cerr<<" ERROR: Could not find RooHistFunction for parameter "<<parName<<"."<<std::endl;
	      return;
	    }
	    if ( iVtx )
	      meanslidesG_wv_2nd[iCat][gCounter-1] = new RooFormulaVar(slideName.Data(),"","@0 + @1 +@0*@2+@0*@3",RooArgList(mnom2,*(pair->second),*nuissancedeltams[iCat],*globalScale));
	    else
	      meanslidesG_rv_2nd[iCat][gCounter-1] = new RooFormulaVar(slideName.Data(),"","@0 + @1 +@0*@2+@0*@3",RooArgList(mnom2,*(pair->second),*nuissancedeltams[iCat],*globalScale));	  
	    
	    slideName=TString::Format("sigmaslideG%d%s%s_2nd",gCounter,catnames.at(iCat).Data(),vtxString.Data());
	    parName  =TString::Format("sigmaG%d_%s_%s_%s",gCounter,vtxString.Data(),procname.Data(),catnamesbase[iCat].Data());
	    pair = fitparmfuncs2nd[iCat].find( parName );
	    if( pair == fitparmfuncs2nd[iCat].end() ) {
	      std::cerr<<" ERROR: Could not find RooHistFunction for parameter "<<parName<<"."<<std::endl;
	      return;
	    }
	    
	    if ( iVtx ) 
	      sigmaslidesG_wv_2nd[iCat][gCounter-1] = new RooFormulaVar(slideName.Data(),"","TMath::Max(0.01,sqrt(@0*@0-@3*@3*@2*@2 +@1*@1))",RooArgList(*(pair->second),*smearmods[iCat],*smears[iCat],mnom2));//why set min 0.01?	
	    else
	      sigmaslidesG_rv_2nd[iCat][gCounter-1] = new RooFormulaVar(slideName.Data(),"","TMath::Max(0.01,sqrt(@0*@0-@3*@3*@2*@2 +@1*@1))",RooArgList(*(pair->second),*smearmods[iCat],*smears[iCat],mnom2));//why set min 0.01?	
	    
	    slideName=TString::Format("fracslide%d%s%s_2nd",gCounter,catnames.at(iCat).Data(),vtxString.Data());
	    parName  =TString::Format("f%d_%s_%s_%s",gCounter,vtxString.Data(),procname.Data(),catnamesbase[iCat].Data());
	    pair = fitparmfuncs2nd[iCat].find( parName );
	    if( pair == fitparmfuncs2nd[iCat].end() ) {
	      std::cerr<<" WARNING: Could not find RooHistFunction for parameter "<<parName<<"."<<std::endl;
	      std::cout<<"          Must be last component."<<std::endl;
	    } else {
	      if( iVtx ) {
		fracslides_wv_2nd[iCat][gCounter-1] = new RooFormulaVar(slideName.Data(),"","@0",RooArgList(*(pair->second)));      	  
		fracList_2nd.add(* (fracslides_wv_2nd[iCat][gCounter-1]) );
	      } else {
		fracslides_rv_2nd[iCat][gCounter-1] = new RooFormulaVar(slideName.Data(),"","@0",RooArgList(*(pair->second)));      	  
		fracList_2nd.add(* (fracslides_rv_2nd[iCat][gCounter-1]) );
	      }
	    }
	    
	    slideName=TString::Format("gslide%d%s%s_2nd",gCounter,catnames.at(iCat).Data(),vtxString.Data());
	    if( iVtx )
	      gslides_wv_2nd[iCat][gCounter-1] = new RooGaussian(slideName.Data(),"",*hmass,*(meanslidesG_wv_2nd[iCat][gCounter-1]),*(sigmaslidesG_wv_2nd[iCat][gCounter-1]));
	    else
	      gslides_rv_2nd[iCat][gCounter-1] = new RooGaussian(slideName.Data(),"",*hmass,*(meanslidesG_rv_2nd[iCat][gCounter-1]),*(sigmaslidesG_rv_2nd[iCat][gCounter-1]));
	    
	    compList_2nd.add( (iVtx ? ( *(gslides_wv_2nd[iCat][gCounter-1]) ): ( *(gslides_rv_2nd[iCat][gCounter-1]) )));
	  }
	  
	  // go to next Gaussian...
	  gCounter++;
	  gIt   = fitparms[iCat].find( TString::Format("dmG%d_",gCounter) + vtxString +TString("_") + procname +TString("_")+catnamesbase[iCat]);
	}
	
	if ( iVtx  ) {
	  combhslides_wv[iCat]      = new RooAddPdf(TString("combhslide")+catnames.at(iCat)+TString("_wv"),"",compList,fracList,true);
	  combhslidesModel_wv[iCat] = new RooAddPdf(TString("combhslideModel")+catnames.at(iCat)+TString("_wv"),"",compListModel,fracList,true);
	} else {
	  combhslides_rv[iCat]      = new RooAddPdf(TString("combhslide")+catnames.at(iCat)+TString("_rv"),"",compList,fracList,true);
	  combhslidesModel_rv[iCat] = new RooAddPdf(TString("combhslideModel")+catnames.at(iCat)+TString("_rv"),"",compListModel,fracList,true);
	}
	
	if( doSecondSignal ) {
	  if ( iVtx  )
	    combhslides_wv_2nd[iCat]      = new RooAddPdf(TString("combhslide")+catnames.at(iCat)+TString("_wv_2nd"),"",compList_2nd,fracList_2nd,true);
	  else
	    combhslides_rv_2nd[iCat]      = new RooAddPdf(TString("combhslide")+catnames.at(iCat)+TString("_rv_2nd"),"",compList_2nd,fracList_2nd,true);
	}
      }

      // add right & wrong together
      
      TString fracName = TString::Format("fracright_%s",catnames[iCat].Data());
      std::map<TString,RooHistFunc*>::iterator pair = fitparmfuncs[iCat].find(fracName);
      if( pair == fitparmfuncs[iCat].end() ) {
	std::cerr<<" ERROR: Could not find RooHistFunction for parameter "<<fracName<<"."<<std::endl;
	return;
      }      
      fracrightmodslides[iCat] = new RooFormulaVar(TString("fracrightmodslide")+catnames.at(iCat),"","TMath::Min(@0*@1,1.0)",RooArgList(nuissancedeltafracright,* (pair->second) ));      
      
      fracName = TString::Format("fracrightModel_%s",catnames[iCat].Data());
      pair = fitparmfuncs[iCat].find(fracName);
      if( pair == fitparmfuncs[iCat].end() ) {
	std::cerr<<" ERROR: Could not find RooHistFunction for parameter "<<fracName<<"."<<std::endl;
	return;
      }
      
      fracrightmodslidesModel[iCat] = new RooFormulaVar(TString("fracrightmodslideModel")+catnames.at(iCat),"","TMath::Min(@0*@1,1.0)",RooArgList(nuissancedeltafracright,* (pair->second) ));      
      
      combhvtxslides     [iCat] = new RooAddPdf(TString("combhvtxslide")     +catnames.at(iCat),"",RooArgList(*combhslides_rv     [iCat],*combhslides_wv     [iCat]),RooArgList(*fracrightmodslides     [iCat]));
      combhvtxslidesModel[iCat] = new RooAddPdf(TString("combhvtxslideModel")+catnames.at(iCat),"",RooArgList(*combhslidesModel_rv[iCat],*combhslidesModel_wv[iCat]),RooArgList(*fracrightmodslidesModel[iCat]));
      
      
      pair = fitparmfuncs[iCat].find(TString::Format("effacc_%s",catnames[iCat].Data()));
      
      if( pair == fitparmfuncs[iCat].end() ) {
	std::cerr<<" ERROR: Could not find RooHistFunction for parameter "<<TString::Format("effacc_%s",catnames[iCat].Data())<<"."<<std::endl;
	return;
      }
      finalnormslides[iCat] = pair->second;
      finalpdfslides[iCat]  = combhvtxslides[iCat];
      
      // ============== BEGIN effective smearing ======================
      pair = fitparmfuncs[iCat].find(TString::Format("effaccModel_%s",catnames[iCat].Data()));
      
      if( pair == fitparmfuncs[iCat].end() ) {
	std::cerr<<" ERROR: Could not find RooHistFunction for parameter "<<TString::Format("effaccModel_%s",catnames[iCat].Data())<<"."<<std::endl;
	return;
      }
      finalnormslidesModel[iCat] = pair->second;
      finalpdfslidesModel[iCat]  = combhvtxslidesModel[iCat];
      
      // ============== END effective smearing ======================

      if( doSecondSignal ) {

	fracName = TString::Format("fracright_%s",catnames[iCat].Data());
	std::map<TString,RooHistFunc*>::iterator pair2nd = fitparmfuncs2nd[iCat].find(fracName);
	if( pair2nd == fitparmfuncs2nd[iCat].end() ) {
	  std::cerr<<" ERROR: Could not find RooHistFunction for parameter "<<fracName<<"."<<std::endl;
	  return;
	}      
	fracrightmodslides_2nd[iCat] = new RooFormulaVar(TString("fracrightmodslide2nd")+catnames.at(iCat),"","TMath::Min(@0*@1,1.0)",RooArgList(nuissancedeltafracright,* (pair2nd->second) ));      
	
	combhvtxslides_2nd     [iCat] = new RooAddPdf(TString("combhvtxslide2nd")     +catnames.at(iCat),"",RooArgList(*combhslides_rv_2nd     [iCat],*combhslides_wv_2nd     [iCat]),RooArgList(*fracrightmodslides_2nd     [iCat]));

	
	pair2nd = fitparmfuncs2nd[iCat].find(TString::Format("effacc_%s",catnames[iCat].Data()));
      
	if( pair2nd == fitparmfuncs2nd[iCat].end() ) {
	  std::cerr<<" ERROR: Could not find RooHistFunction for parameter "<<TString::Format("effacc_%s",catnames[iCat].Data())<<"."<<std::endl;
	  return;
	}
	finalnormslides_2nd[iCat] = pair2nd->second;
	finalpdfslides_2nd[iCat]  = combhvtxslides_2nd[iCat];
      }
      // ========== END SECOND MODEL =============
    }

    // FIX-ME: Add back interpolation test...
    
    // open again the config file and read the normalization buisiness...
    // ... Catgeory normalizations
    std::vector<RooAbsReal*> nsigcats;
    std::vector<RooAbsReal*> nsigcats_2nd;
    std::vector<RooAbsReal*> nsigcatsModel; /// ======== BEGIN =========
    std::vector<RooAbsReal*> addNuissance;
    
    std::cout<<" ======================================================== "<<std::endl;

    status = readConfigCardNuissances(configCardName, numCats,
				      nsigcats,
				      nsigcats_2nd,
				      nsigcatsModel,
				      addNuissance,
				      finalnormslides,
				      finalnormslides_2nd,
				      finalnormslidesModel,
				      catnames);

    
    std::cout<<" ======================================================== "<<std::endl;

    if(!status) {
      std::cerr<<" ERROR when readin input card "<<configCardName<<"."<<std::endl;
      return;
    }
    
    RooArgList addnorm;
    RooArgList addnormff;
    
    RooArgList addpdfs;
    RooArgList addpdfsff;
    
    RooArgList addcoeffs;  
    RooArgList addcoeffsff;  
    
    std::vector<RooAbsReal*> combnorms;
    std::vector<RooAbsReal*> combnormsff;

    RooArgList addnorm_2nd;
    RooArgList addnormff_2nd;
    
    RooArgList addpdfs_2nd;
    RooArgList addpdfsff_2nd;
    
    RooArgList addcoeffs_2nd;  
    RooArgList addcoeffsff_2nd;  
    
    std::vector<RooAbsReal*> combnorms_2nd;
    std::vector<RooAbsReal*> combnormsff_2nd;
    
    //final pdfs
    for (UInt_t icat=0; icat<catnames.size(); ++icat) {

      RooAbsPdf *hggpdfsmabs = finalpdfslides[icat];
      hggpdfsmabs->SetName(TString("hggpdfsmabs_")+catnames.at(icat));
      
      RooAbsPdf *hggpdfsmrel = (RooAbsPdf*)finalpdfslides[icat]->Clone(TString::Format("hggpdfsmrel_%s",catnames.at(icat).Data()));//a function of mnom
      // ========= BEGIN ==========
      RooAbsPdf *hggpdfsmrelModel = (RooAbsPdf*)finalpdfslidesModel[icat]->Clone(TString::Format("hggpdfsmrelModel_%s",catnames.at(icat).Data()));//a function of mnom
      // ==========================
      
      RooAbsReal* tempAbs = addFuncMap.find("smxsec")->second;
      if (!tempAbs) {
	std::cerr<<" ERROR: Could not find addFunc with name smxsec."<<std::endl;
	return;
      }
      
      if ( !(procxseclist[iProc]) ) {
	std::cout<< "Missing proclist "<<std::endl;
	return;
      }
      if ( !nsigcats[icat] ) {
	std::cout<< "Missing nsigcats "<<std::endl;
	return;
      }

      if ( !nsigcatsModel[icat] ) {
	std::cout<< "Missing nsigcats "<<std::endl;
	return;
      }
      
      RooFormulaVar *nsigsmabs = new RooFormulaVar(TString::Format("hggpdfsmabs_%s_norm",catnames.at(icat).Data()),"","@0*@1/@2",RooArgList(*(procxseclist[iProc]),*nsigcats[icat],*tempAbs));

      if ( histFuncMap.find("smbr")== histFuncMap.end()) {
	std::cout<<" can't find smbr "<<std::endl;
	return;
      }
      
      tempAbs = histFuncMap.find("smbr")->second;
      if (!tempAbs) {
	std::cerr<<" ERROR: Could not find histFunc with name smbr."<<std::endl;
	return;
      }
      
      RooFormulaVar *nsigsmrel      = new RooFormulaVar(TString::Format("hggpdfsmrel_%s_norm",catnames.at(icat).Data()),"","@0*@1*@2",RooArgList( *tempAbs, *procxseclist[iProc],*nsigcats[icat]));

      RooExtendPdf *sigpdfsmabs = new RooExtendPdf(TString::Format("sigpdfsmabs%s",catnames.at(icat).Data()),"",*hggpdfsmabs,*nsigsmabs);
      RooExtendPdf *sigpdfsmrel = new RooExtendPdf(TString::Format("sigpdfsmrel%s",catnames.at(icat).Data()),"",*hggpdfsmrel,*nsigsmrel);

      // ========    BEGIN ============
      RooFormulaVar *nsigsmrelModel = new RooFormulaVar(TString::Format("hggpdfsmrelModel_%s_norm",catnames.at(icat).Data()),"","@0*@1*@2",RooArgList( *tempAbs, *procxseclist[iProc],*nsigcatsModel[icat]));

      RooExtendPdf *sigpdfsmrelModel = new RooExtendPdf(TString::Format("sigpdfsmrelModel%s",catnames.at(icat).Data()),"",*hggpdfsmrelModel,*nsigsmrelModel);

      addPdfCoefList[icat].add(*nsigsmrelModel);
      addPdfList[icat].add(*sigpdfsmrelModel);

      addPdfCoefList[numCats].add( * (RooFormulaVar*) (nsigsmrelModel->clone(TString::Format("hggpdfsmrelModel_clone_%s_norm",catnames.at(icat).Data()))));
      addPdfList[numCats].add( * (RooExtendPdf*)  (sigpdfsmrelModel->clone(TString::Format("sigpdfsmrelModelclone%s",catnames.at(icat).Data()))));
      
      // ===============================
      
      if(!fitonly) {
	wOut->import(*sigpdfsmabs,RecycleConflictNodes());
	wOut->import(*sigpdfsmrel,RecycleConflictNodes());
	//wOut->import(*sigpdfsmrelModel,RecycleConflictNodes());
      }

      addnorm.add(*nsigsmrel);
      combnorms.push_back(nsigsmrel);
    }

    //--------------------------------------------------------------------------------
    if( doSecondSignal ) {
      
      //final pdfs
      for (UInt_t icat=0; icat<catnames.size(); ++icat) {
      
	RooAbsPdf *hggpdfsmabs = finalpdfslides_2nd[icat];
	hggpdfsmabs->SetName(TString::Format("hggpdfsmabs_%s_2nd",catnames.at(icat).Data()).Data());
	
	RooAbsPdf *hggpdfsmrel = (RooAbsPdf*)finalpdfslides_2nd[icat]->Clone(TString::Format("hggpdfsmrel_%s_2nd",catnames.at(icat).Data()));//a function of mnom
	
	RooAbsReal* tempAbs = addFuncMap2nd.find("smxsec")->second;
	if (!tempAbs) {
	  std::cerr<<" ERROR: Could not find addFunc with name smxsec."<<std::endl;
	  return;
	}	
	if ( !(procxseclist[iProc]) ) {
	  std::cout<< "Missing proclist "<<std::endl;
	  return;
	}
	if ( !nsigcats_2nd[icat] ) {
	  std::cout<< "Missing nsigcats "<<std::endl;
	  return;
	}
      
	RooFormulaVar *nsigsmabs = new RooFormulaVar(TString::Format("hggpdfsmabs_%s_2nd_norm",catnames.at(icat).Data()),"","@0*@1/@2*@3",RooArgList(*(procxseclist[iProc]),*nsigcats_2nd[icat],*tempAbs,the2ndMu));
	
	if ( histFuncMap2nd.find("smbr")== histFuncMap2nd.end()) {
	  std::cout<<" can't find smbr "<<std::endl;
	  return;
	}
	
	tempAbs = histFuncMap2nd.find("smbr")->second;
	if (!tempAbs) {
	  std::cerr<<" ERROR: Could not find histFunc with name smbr."<<std::endl;
	  return;
	}
	
	RooFormulaVar *nsigsmrel      = new RooFormulaVar(TString::Format("hggpdfsmrel_%s_2nd_norm",catnames.at(icat).Data()),"","@0*@1*@2*@3",RooArgList( *tempAbs, *procxseclist[iProc],*nsigcats_2nd[icat],the2ndMu));
	
	RooExtendPdf *sigpdfsmabs = new RooExtendPdf(TString::Format("sigpdfsmabs%s_2nd",catnames.at(icat).Data()),"",*hggpdfsmabs,*nsigsmabs);
	RooExtendPdf *sigpdfsmrel = new RooExtendPdf(TString::Format("sigpdfsmrel%s_2nd",catnames.at(icat).Data()),"",*hggpdfsmrel,*nsigsmrel);
	
	if(!fitonly) {
	  wOut->import(*sigpdfsmabs,RecycleConflictNodes());
	  wOut->import(*sigpdfsmrel,RecycleConflictNodes());
	}
      }
      //--------------------------------------------------------------------------------
    }
  } // go to next process...
  
  if( !fitonly )
    wOut->writeToFile( TString::Format("%s/model/ubersignalmodel.root",projectDir.Data() ) );
  
  // =================   BEGIN effective smearing.... ===================
  TCanvas***  effSmearCan  = new TCanvas**[numCats+1];
  TLegend***  legend       = new TLegend**[numCats+1];
  TLatex***   catdesc_l    = new TLatex** [numCats+1];
  TArrow***   fwhmArrow    = new TArrow** [numCats+1];
  TLatex***   fwhmText     = new TLatex** [numCats+1];
  RooCurve*** nomcurve     = new RooCurve**[numCats+1];
  RooHist***  datacurve    = new RooHist** [numCats+1];
  RooPlot***  hplot        = new RooPlot**[numCats+1];
  TGraph***   effSigmaGraph = new TGraph**[numCats+1];

  if( doeffsigma ) {
    
    TH1D* dummyData = new TH1D("dummyData","",1,0.,10.);
    dummyData->SetMarkerStyle(25);
    dummyData->SetMarkerSize(1.2);

    TLatex* cmsprel = new TLatex(0.15,0.95,"CMS preliminary");
    cmsprel->SetTextAlign(12);
    cmsprel->SetTextSize(0.04);
    cmsprel->SetNDC();

    //for( int iCat = 0; iCat < numCats + 1; ++iCat) {
    for( int iCat = 0; iCat < numCats; ++iCat) {
      
      bool onForAllProcs = true;  // Cat must be ON for all PROCS, otherwise dont' do the signal model
      if ( iCat < numCats ) {
	for(int iProc =0; iProc < numProcs; ++iProc)
	  onForAllProcs = onForAllProcs && catIsOn[iProc*numCats+iCat];
      }
      if ( !onForAllProcs ) continue;  // don't do fits in this Cat
      
      effSmearCan[iCat]   = new TCanvas*[mhs.size()];
      legend     [iCat]   = new TLegend*[mhs.size()];
      catdesc_l  [iCat]   = new TLatex* [mhs.size()];
      fwhmArrow  [iCat]   = new TArrow* [mhs.size()];
      fwhmText   [iCat]   = new TLatex* [mhs.size()];
      nomcurve   [iCat]   = new RooCurve*[mhs.size()];
      datacurve  [iCat]   = new RooHist*[mhs.size()];
      hplot      [iCat]   = new RooPlot*[mhs.size()];
      effSigmaGraph[iCat] = new TGraph*[mhs.size()];

      effSmearPdf[iCat] = new RooAddPdf(TString::Format("effsmearpdf_%d",iCat).Data(),"",addPdfList[iCat],addPdfCoefList[iCat]);
      
      if (iCat < numCats)
	std::cout<<"  ----------------- CATEGORY "<<catnamesbase[iCat]<<" --------------------- "<<std::endl;
      else
	std::cout<<"  ----------------------- ALL CATEGORIES  --------------------- "<<std::endl;
      
      for(unsigned int iMass = 0; iMass < mhs.size(); ++iMass) {
	
	// only do this for the test mass
	//if( mhs_fit[iMass] > 0 ) continue;
	
	effSmearCan[iCat][iMass] = new TCanvas();
	catdesc_l  [iCat][iMass] = new TLatex(0.95,0.95,( iCat < numCats ? catdesc[iCat].Data() : "All Classes") );
	catdesc_l  [iCat][iMass]->SetNDC();
	catdesc_l  [iCat][iMass]->SetTextAlign(32);
	catdesc_l  [iCat][iMass]->SetTextSize(0.04);

	effSmearCan[iCat][iMass] ->cd();
	

	//std::cout<<"  setting mass to "<<mhs[iMass]<<std::endl;
	mnom.setVal(mhs[iMass]);
	
	//hmass2->setRange("plotrange",mhs.at(iMass)-12.,mhs.at(iMass)+8.);
	hmass2->setRange("plotrange",mhs.at(iMass)-15.,mhs.at(iMass)+10.);
	hplot[iCat][iMass] = hmass2->frame(Bins(50),Range("plotrange"));
	hmass2->setRange("fitrange",mhs.at(iMass)-fitmassmin,mhs.at(iMass)+fitmassmax);
	effSmearData[iCat][iMass]->plotOn(hplot[iCat][iMass],RooFit::MarkerStyle(25),RooFit::MarkerSize(1.2));
	effSmearPdf [iCat]       ->plotOn(hplot[iCat][iMass],NormRange("fitrange"));
	hplot[iCat][iMass]->SetMinimum(1e-5);
	
	nomcurve [iCat][iMass] = new RooCurve(*hplot[iCat][iMass]->getCurve(hplot[iCat][iMass]->nameOf(1)));
	datacurve[iCat][iMass] = new RooHist(*hplot[iCat][iMass]->getHist (hplot[iCat][iMass]->nameOf(0)));
	
	hplot[iCat][iMass]->SetTitle("");
	//hplot[iCat][iMass]->Draw();
	
	// find the maximum...
	
	double lowBound = -1.;
	double upBound  = -1.;
	
	double theVal   = -1.;
	double binSize = 0.001;
	
	bool useBinnedForEffSigma = true;

	if( iCat < numCats ) {

	if ( useBinnedForEffSigma ) {	  	  
	  effSmearHists[iMass][iCat] = (TH1D*) effSmearPdf[iCat]->createHistogram(TString::Format("effhist_%d_%s",iMass,(iCat < numCats ? catnamesbase[iCat].Data() : "combined" )).Data(), *hmass2, RooFit::Binning((int) ((massmax-massmin)/binSize)));
	  computeEffSigmaBinned( effSmearHists[iMass][iCat], mhs.at(iMass), lowBound, upBound, theVal, 0.1, binSize);
	} else
	  computeEffSigma( effSmearData[iCat][iMass], hmass2, mhs.at(iMass), lowBound, upBound, theVal, 0.1, binSize);

	}

	// cdreate the TGraph for the eff-sigma thingy...
	int numPoints = (upBound-lowBound) / binSize + 3;
	double* effSigmaShapeX = new double[numPoints];
	double* effSigmaShapeY = new double[numPoints];
	
	effSigmaShapeX[0] = lowBound;
	effSigmaShapeY[0] = 0.;
	
	for(int iP =0; iP < numPoints - 2; ++iP){
	  
	  effSigmaShapeX[iP+1] = lowBound + iP*binSize;
	  effSigmaShapeY[iP+1] = nomcurve[iCat][iMass]->interpolate(lowBound + iP*binSize);
	}
	
	effSigmaShapeX[numPoints-1] = upBound;
	effSigmaShapeY[numPoints-1] = 0.;
	
	
	effSigmaGraph[iCat][iMass] = new TGraph(numPoints,effSigmaShapeX,effSigmaShapeY);
	
	effSigmaGraph[iCat][iMass]->SetFillColor(18);
	effSigmaGraph[iCat][iMass]->SetLineColor(kGray);
	
	hplot[iCat][iMass]->addObject(effSigmaGraph[iCat][iMass],"FL");
	hplot[iCat][iMass]->addPlotable(datacurve[iCat][iMass],"P");
	hplot[iCat][iMass]->addObject(nomcurve[iCat][iMass]);
	
	hplot[iCat][iMass]->Draw();
	

	legend[iCat][iMass] = new TLegend(0.18,0.65,0.49,0.9);
	legend[iCat][iMass]->SetFillColor(0);
	legend[iCat][iMass]->SetFillStyle(0);
	legend[iCat][iMass]->AddEntry(dummyData,"Simulation","PLE");
	legend[iCat][iMass]->AddEntry(nomcurve[iCat][iMass],"Parametric Model","L");
	legend[iCat][iMass]->AddEntry(effSigmaGraph[iCat][iMass],TString::Format("#sigma_{eff} = %.2f GeV",(upBound-lowBound)/2.).Data(),"F");
	legend[iCat][iMass]->Draw();
	
	cmsprel->Draw();
	catdesc_l[iCat][iMass]->Draw();
	
	// =============== DO FWHM COMPUTATION =======================
	// 1. find the maximum
	double thePos    = mhs.at(iMass) - 5.;
	double maxVal    = -10.;
	double maxPos    = -10.;
	while ( true ) {
	  //hmass2->setVal(thePos);
	  //double theVal = effSmearPdf[iCat]->getVal();
	  double theValNom = nomcurve[iCat][iMass]->interpolate(thePos);
	  //if ( nomcurve[iCat][iMass]->interpolate(thePos) > maxVal ) {
	  if ( theValNom > maxVal ) {
	    maxPos = thePos;
	    //maxVal = nomcurve[iCat][iMass]->interpolate(thePos);
	    maxVal = theValNom;
	  } else break; // we're going down again...
	  thePos += binSize;
	}
	// 2. ... so half-max-value is:
	double HMValue = maxVal/2.;
	
	// 3. find the left end of the arrow
	double minDiff = 10000000.;
	double leftEnd = mhs.at(iMass) - 5.;
	thePos    = mhs.at(iMass) - 5.;

	while ( thePos < maxPos ) {
	  //hmass2->setVal(thePos);
	  //double theVal = effSmearPdf[iCat]->getVal();
	  double theValNom = nomcurve[iCat][iMass]->interpolate(thePos);
	  //double theDiff = TMath::Abs(nomcurve[iCat][iMass]->interpolate(thePos) - HMValue);
	  double theDiff = TMath::Abs(theValNom - HMValue);
	  if (theDiff < minDiff) {
	    minDiff = theDiff;
	    leftEnd = thePos;
	  } else break; // gone to far...
	  thePos += binSize;
	}
	
	// 4. find the right end of the arrow
	minDiff = 100.;
	double rightEnd = mhs.at(iMass);
	thePos    = maxPos;
	
	while ( thePos < ( mhs.at(iMass) + 5.) ) {
	  //hmass2->setVal(thePos);
	  //double theVal = effSmearPdf[iCat]->getVal();
	  double theValNom = nomcurve[iCat][iMass]->interpolate(thePos);
	  //double theDiff = TMath::Abs(nomcurve[iCat][iMass]->interpolate(thePos) - HMValue);
	  double theDiff = TMath::Abs(theValNom - HMValue);
	  if (theDiff < minDiff) {
	    minDiff = theDiff;
	    rightEnd = thePos;
	  } else break; // getting worse... so break;
	  thePos += binSize;
	}
	
	// 5. plot the arrow
	fwhmArrow[iCat][iMass] = new TArrow(leftEnd,HMValue,rightEnd,HMValue,0.03,"<>");
	fwhmArrow[iCat][iMass]->Draw();

	fwhmText [iCat][iMass] = new TLatex(leftEnd-2.,HMValue,TString::Format("FWHM = %.2f GeV",(rightEnd-leftEnd)).Data());
	fwhmText [iCat][iMass]->SetTextAlign(32);
	fwhmText [iCat][iMass]->SetTextSize(0.035);
	fwhmText [iCat][iMass]->Draw();
	
	std::cout<<"              Found range as [ "<<lowBound<<" - "<<upBound<<" ] with CI = "<<theVal*100.<<" % --> eff. sigma = "<<(upBound-lowBound)/2.<<"  ( "<<(upBound-lowBound)/2./mhs.at(iMass)*100<<" % )"<<std::endl;
	std::cout<<"                        FWHM   "<<(rightEnd-leftEnd)<<" GeV  ( "<<(rightEnd-leftEnd)*100./mhs.at(iMass)<<" % )"<<std::endl;
	double numEventsTot = effSmearData   [iCat][iMass]->sumEntries();
	std::cout<<"     #Events ( L = "<<lumi<<" fb-1 ) = "<< numEventsTot << std::endl;
	double debugAll = 0.;
	if ( iCat < numCats ) {
	  for( int iProc =0; iProc < numProcs; ++iProc ) {
	    TString procname   = procIdxMap.find(iProc)->second;
	    if( !(procon.find(procname)->second) ) continue;
	    TString dataAllName= TString::Format("sig_%s_massModelNom_m%d_%s",procname.Data(),(int) mhs.at(iMass),catnamesbase[iCat].Data());	
	    std::map<TString,RooDataSet*>::iterator tmpDataIt = keepAllDataMap.find(dataAllName);
	    if (tmpDataIt == keepAllDataMap.end() ) {
	      std::cerr<<"  ERROR: Could not find Dataset with name "<<dataAllName<<"."<<std::endl;
	      return;
	    }
	    std::cout<<"                    (  PROC = "<<procname<<" :  "<<(tmpDataIt->second)->sumEntries()/numEventsTot*100.<<" % )"<<std::endl;
	    debugAll += (tmpDataIt->second)->sumEntries();
	  }
	}
	//std::cout<<"   DEBUG  "<<debugAll<<"      "<<numEventsTot<<std::endl;
	
	
	//if( doeffsmear &&  (mhs[iMass] > 122.5 && mhs[iMass] < 127.5) && iCat < numCats ) {
	if( doeffsmear ) {
	  
	  double effSigma1 = (upBound-lowBound)/2.;
	  
	  // compute the no-smear effective width...
	  if ( useBinnedForEffSigma ) {
	    computeEffSigmaBinned( effSmearHistsNS[iMass][iCat], mhs.at(iMass), lowBound, upBound, theVal, 0.1, binSize);
	  } else
	    computeEffSigma( effSmearDataNS[iCat][iMass], hmass2, mhs.at(iMass), lowBound, upBound, theVal, 0.1, binSize);
	  
	  std::cout<<"  NO    SMEAR Found range as [ "<<lowBound<<" - "<<upBound<<" ] with CI = "<<theVal*100.<<" % --> eff. sigma = "<<(upBound-lowBound)/2.<<"  ( "<<(upBound-lowBound)/2./mhs.at(iMass)*100<<" % )"<<std::endl;	
	  double effSigma2 = (upBound-lowBound)/2.;
	  
	  std::cout<<" [ "<<mhs.at(iMass)<<" ]   REL. EFF. SMEARING    = "<<TMath::Sqrt(TMath::Max(0.0,effSigma1*effSigma1-effSigma2*effSigma2))/mhs.at(iMass)<<std::endl;

// 	  if ( mhs_fit[iMass] == 0 )
// 	    smears[iCat] -> setVal(TMath::Sqrt(TMath::Max(0.0,effSigma1*effSigma1-effSigma2*effSigma2))/mhs.at(iMass));
	  
	}
	
	if( !fitonly )
	  effSmearCan[iCat][iMass]->SaveAs(TString::Format("%s/model/sigmodel_%d_%s.pdf",projectDir.Data(),(int) mhs.at(iMass),( iCat < numCats ? catnamesbase[iCat].Data() : "combined" )));
      }
    }
  }

  std::cout<<" 8TeV analysis done."<<std::endl;
  
  // print out the summary
  for (int iCat = 0; iCat < numCats; ++iCat) {
    std::cout<<" Results for Cat: "<<catnamesbase[iCat]<<std::endl;
    for(int iProc = 0; iProc < numProcs; ++iProc) {
      TString procname = procIdxMap.find(iProc)->second;
      std::cout<<"   proc = "<<procname<<std::endl;
      
      // check if proc was on for this cat..
      procValMap::iterator it = effaccMap.find(procname);
      if ( it == effaccMap.end() ) {
	std::cout<<"     was OFF"<<std::endl;
	continue;
      }
      
      catValMap::iterator it2 = it->second->find(catnamesbase[iCat]);
      if (it2 == it->second->end() ) {
	std::cout<<"     was OFF"<<std::endl;
	continue;
      }
      
      for(unsigned int iMass = 0; iMass < mhs.size(); ++iMass) 
	std::cout<<"        mass: eff x acc = "<<mhs[iMass]<<"    "<<it2->second->find(mhs[iMass])->second<<std::endl;
    }
  }
  
  return;
  
}


RooAbsPdf* generateMultiGaussian(RooRealVar* mass, 
				 RooRealVar* nomMass,
				 TString procName, TString quali,
				 std::map<TString,RooRealVar*>& fitparms,
				 std::map<TString,float>&       startvals,
				 std::map<TString,int>&         valsfixed,
				 bool rightVtx) {

  RooArgList compList;
  RooArgList fracList;

  TString vtxString = "rv";
  if( !rightVtx )
    vtxString = "wv";

  // now loop over the Gaussians and add them
  TString mName = "dmG1"+vtxString;
  TString sName = "sigmaG1"+vtxString;
  TString fName = "f1"+vtxString;
  std::map<TString,float>::iterator g_mean   = startvals.find(mName);
  std::map<TString,float>::iterator g_sigma  = startvals.find(sName);
  std::map<TString,float>::iterator g_f      = startvals.find(fName);
  int gCounter = 1;
  TString parName;
  bool isLast = false;
  while ( g_mean  != startvals.end() ) {
    if ( g_sigma == startvals.end() || g_f == startvals.end() ) {
      std::cout<<( ( g_sigma == startvals.end() || isLast )? " ERROR" : "WARNING")<<": Need a complete Gaus information.  Sigma ? "<<( g_sigma == startvals.end() )<<"  Frac ? "<<( g_f == startvals.end() )<<std::endl;
      if (g_sigma == startvals.end() || isLast )
	return NULL;
      else {
	isLast = true;
	std::cout<<"   Must be last Gauss. Test next round..."<<std::endl;
      }
    }    

    parName = TString::Format("dmG%d_%s_%s_%s",gCounter,vtxString.Data(),procName.Data(),quali.Data());
    RooRealVar* dmG = new RooRealVar( parName.Data(),"",-5.,5.);
    dmG -> removeRange();
    if ( valsfixed.find(mName)->second )
      dmG->setConstant();
    fitparms.insert (std::pair<TString, RooRealVar*>(parName,dmG));
    double tempVal = g_mean->second;
    startvals.erase(mName);
    startvals.insert(std::pair<TString,float>( parName,tempVal));

    parName = TString::Format("sigmaG%d_%s_%s_%s",gCounter,vtxString.Data(),procName.Data(),quali.Data());
    RooRealVar* sigmaG = new RooRealVar( parName.Data(),"",-5.,5.);
    sigmaG -> removeRange();
    if ( valsfixed.find(sName)->second )
      sigmaG->setConstant();
    fitparms.insert (std::pair<TString, RooRealVar*>(parName,sigmaG));
    tempVal = g_sigma->second;
    startvals.erase(sName);
    startvals.insert(std::pair<TString,float>( parName,tempVal));

    if( !(g_f == startvals.end()) ) {
      parName = TString::Format("f%d_%s_%s_%s",gCounter,vtxString.Data(),procName.Data(),quali.Data());
      RooRealVar* f = new RooRealVar( parName.Data(),"",0.,1.);
      if ( valsfixed.find(fName)->second )
	f->setConstant();
      fitparms.insert (std::pair<TString, RooRealVar*>(parName,f));
      tempVal = g_f->second;
      startvals.erase(fName);
      startvals.insert(std::pair<TString,float>( parName,tempVal));
      fracList.add(*f);
    }

    parName = TString::Format("meanG%d_%s_%s_%s",gCounter,vtxString.Data(),procName.Data(),quali.Data());
    RooFormulaVar* meanG = new RooFormulaVar(parName.Data(),"","@0+@1",RooArgList(*nomMass,*dmG));
    
    // construct the Gaussian & add it to the list
    parName = TString::Format("g%d_%s_%s_%s",gCounter,vtxString.Data(),procName.Data(),quali.Data());
    RooGaussian* g = new RooGaussian(parName.Data(),"",*mass,*meanG,*sigmaG);

    compList.add(*g);

    gCounter++;
    mName = TString::Format("dmG%d%s",gCounter,vtxString.Data());
    sName = TString::Format("sigmaG%d%s",gCounter,vtxString.Data());
    fName = TString::Format("f%d%s",gCounter,vtxString.Data());

    g_mean   = startvals.find(mName);
    g_sigma  = startvals.find(sName);
    g_f      = startvals.find(fName);

  }

  parName = TString::Format("combh_%s_%s_%s",vtxString.Data(),procName.Data(),quali.Data());
  RooAddPdf* comb = new RooAddPdf(parName.Data(),"",compList,fracList,true);
  //RooAddPdf* comb = new RooAddPdf(parName.Data(),"",compList,fracList);
  RooAbsPdf* returnPdf = (RooAbsPdf*) comb;

  return returnPdf;
}


bool validateInput(int iProc, int iCat, int nProcs, int nCats) {

  if(iProc >= nProcs || iCat >= nCats ) {
    std::cerr<<" Error in Settings: "<<iProc<<"  "<<iCat<<"  "<<nProcs<<"  "<<nCats<<std::endl;
    return false;
  }

  return true;
}





bool readConfigCard(TString configCardName, 
		    std::map<TString,TString>& procMCfileMap,
		    std::map<TString,bool>&    procOn,
		    std::map<int,TString>&     procIdxMap,
		    int& numProcs, int& numCats, int& numMasses,
		    bool*& catIsOn,
		    std::vector<double>& mhs,
		    std::vector<int>& mhs_fit,
		    std::vector<double>& smearingv,
		    std::vector<TString>& catNames,
		    std::vector<TString>& catDesc,
		    std::map<TString,TString>& catToSmearMap,
		    std::map<TString,float>*& startvals,
		    std::map<TString,int>*& valsfixed) { 

  
  FILE* configFile = fopen(configCardName.Data(),"r");
  if ( !configFile ) {
    std::cerr<<" Inputfile "<<configCardName<<" not found."<<std::endl;
    return false;
  }
  
  mhs.resize(0);
  
  char line[1000];
  
  std::cout<<" [INFO] Reading model from file < "<<configCardName<<" >."<<std::endl;
  
  int   iProc = -1;
  int   iCat  = -1;
  int   massIdx = -1;
  
  int   whichCat = -1 ;
  float smearing= -1.;
  float massval  = -1.;
  
  char rightStart[100];
  
  char dummyString[100];
  
  //# ------------------------------------------------------------------------------------------------------------------------------------------------
  //# Section on the Porcesses/Events
  //#	num proc	num photon cats		num cats(auxiliary)	num cats (analysis)	num cats (trigger)	numModels	numMasses
  //# ------------------------------------------------------------------------------------------------------------------------------------------------
  //INIT	4		4			7			4			16			1		5
  while (fgets(line,1000,configFile)) {
    if ( !sscanf(&line[0], "#") ) {   // not a document line      
      if ( sscanf(line, "INIT %d %d %d", &numProcs, &numCats, &numMasses) ) { // this is the INIT line	
	if( numProcs < 1 || numCats < 1  ){
	  std::cerr<<" Error in config file "<<configCardName<<" : Number of processes and categories must be positive."<<std::endl;
	  return false;
	}

	catIsOn = new bool[numProcs*numCats];

	smearingv.resize(numCats);
	catNames.resize(numCats);
	catDesc.resize(numCats);

	//_startVals       = new float*[numProcs*numCats];
	startvals        = new std::map<TString,float>[numProcs*numCats];
	valsfixed        = new std::map<TString,int>[numProcs*numCats];

	//std::cout<<" startVals has size = "<<numProcs*numCats<<std::endl;

	// # ------------------------------------------------------------------------------------------------------------------------------------------------
	// # analysis categpry definitions
	// #	idx	name		eff. smearing	BG-model/Order	TCuts *use only AUXCATs from above*		StDesc(String descriptor) (*for plots and such*)
	// # ------------------------------------------------------------------------------------------------------------------------------------------------
	// ANACAT	0	cat0		0.005432	Bern/5		" basecut && !vbfcut && baseline0 "		StDesc(Baseline Cat 1)
	// ANACAT	1	cat1		0.005196	Bern/5		" basecut && !vbfcut && baseline1 "		StDesc(Baseline Cat 2)
	// ANACAT	2	cat2		0.007464	Bern/5		" basecut && !vbfcut && baseline2 "		StDesc(Baseline Cat 3)
	// ANACAT	3	cat3		0.012978	Bern/5		" basecut && !vbfcut && baseline3 "		StDesc(Baseline Cat 4)
	// ANACAT	4	cat4		0.008722	Bern/3		" masscut && vbfcut "				StDesc(Baseline VBVA Cat)
      } else if (  sscanf(line, "ANACAT %d %s %f ", &whichCat, dummyString, &smearing) ) { // this is the SMEARING line
	if(whichCat < 0 || whichCat >= numCats) {
	  std::cerr<<" Error in config file "<<configCardName<<" : Setting smearing for Cat "<<whichCat<<" with only "<<numCats<<" total Cats."<<std::endl;
	  return false;
	}
	smearingv[whichCat] = (double) smearing;
	catNames[whichCat]  = TString(dummyString);
	std::string theLine = line;
	size_t fPos2 = theLine.rfind("Desc(");
	theLine.erase(0,fPos2+5);
	fPos2 = theLine.find_last_of(")");
	theLine.erase(fPos2);
	catDesc[whichCat] = TString(theLine);
	
	if ( true ) {
	  // find smearing/scalinf cat. This defines how smearing/scaling is correlated among cats (also with other analyses)
	  theLine = line;
	  fPos2 = theLine.rfind("Smear(");
	  theLine.erase(0,fPos2+6);
	  fPos2 = theLine.find_first_of(")");
	  theLine.erase(fPos2);	  
	  catToSmearMap.insert(std::pair<TString,TString>(catNames[whichCat],TString(theLine)));

	}
      } else if  (  
		  sscanf(line, "OFF %d %d ( %s )", &iProc, &iCat, rightStart) ||
		  sscanf(line, "ON  %d %d ( %s )", &iProc, &iCat, rightStart)
		  ) {
	
	if( iProc < 0 || iCat < 0 ) break;
	if (!validateInput(iProc,iCat,numProcs,numCats)) break;
	if ( sscanf(line, "OFF %s", dummyString) )
	  catIsOn[iProc*numCats+iCat] = false;
	else
	  catIsOn[iProc*numCats+iCat] = true;
	
	std::cout<<" Cat #"<<iCat<<" for Proc #"<<iProc<<" is "<<(catIsOn[iProc*numCats+iCat] ? "on" : "off")<<std::endl;

	// need to have TWO line, first is right Vtx, second wrniog Vtx hypothesis
	for( int iVtx = 0; iVtx < 2; ++iVtx ) {
	  std::string parStr(line);
	  size_t startPos = parStr.find_first_of("(");
	  size_t endPos   = parStr.find_first_of(")");
	  parStr = parStr.substr(startPos+1, endPos-startPos-1);
	  size_t lastPos = parStr.find_first_not_of(" ");
	  while ( parStr.find_first_not_of(" ") != std::string::npos ) {
	    // string not empty...
	    size_t pos = parStr.find_first_of(",");
	    std::string thisPar = parStr.substr(lastPos,pos);
	    size_t idx = thisPar.find_first_of(":");

	    // check if the parameter should be fixed
	    bool isFixed = ( thisPar.find_first_of("X") != std::string::npos );

	    float parValue = -1.;
	    char  parName[10];
	    //sscanf(thisPar.substr(0,idx).c_str(),"%s",&parName,&parValue);
	    sscanf(thisPar.substr(0,idx).c_str(),"%s",parName);
	    sscanf( ( isFixed ? thisPar.substr(idx+1,thisPar.find_first_of("X")-idx-1).c_str() : thisPar.substr(idx+1).c_str() ),"%f",&parValue);
	    TString parNameString = TString(parName) + (iVtx ? TString("wv"): TString("rv"));
	    startvals[iProc*numCats+iCat].insert(std::pair<TString,float>(parNameString,parValue));
	    valsfixed[iProc*numCats+iCat].insert(std::pair<TString,int>  (parNameString,(int) isFixed));
	    parStr.erase(lastPos,pos);
	    lastPos = parStr.find_first_not_of(" ");
	  }
	  if( !iVtx ) {
	    // read the next line... must have correct shape
	    fgets(line,1000,configFile);
	    while ( sscanf(&line[0], "#") ) fgets(line,1000,configFile); 
	    if (!sscanf(line, "%s ( %s )", dummyString, rightStart))
	      std::cerr<<"    Line for starting values has wrong format: "<<line<<std::endl;
	  }
	}
	
      } else {
	
	// additional setup
	int theProc = -1;
	char name[30];

	char onoff[3];
	//char directory[100];
	//char procCSLabel[400];
    
	char MCfileName[200];

	//int testmass = 0;

	if ( sscanf(line,"PROC %d %s %s file:%s",&theProc, name, onoff, MCfileName) ) {
	  //# ------------------------------------------------------------------------------------------------------------------------------------------------
	  //PROC    0	ggh	OFF			file:/scratch/fabstoec/cms/hist/hgg-7TeV-janReReco/merged/hgg-7TeV-janReReco_f11--h%dgg-gf-v14b-pu_noskim.root
	  std::cout<<" adding proc with idx = "<<theProc<<"  and name "<<name<<std::endl;
	  procMCfileMap.insert(std::pair<TString,TString>(TString(name)      ,TString(MCfileName)));
	  procIdxMap   .insert(std::pair<int,TString>    (theProc            ,TString(name)));
	  procOn       .insert(std::pair<TString,bool>   (TString(name)      ,strcmp(onoff,"OFF")));

// 	if( sscanf(line,"PROC %d %s %s",&theProc,&name,&onoff) ) {
// 	  procNames[theProc]=TString(name).Strip();
// 	  procOn[theProc]=strcmp(onoff,"OFF");
	}
	else if ( sscanf(line,"MASS %d %f",&massIdx, &massval) ) {
	  if( massIdx > (int) mhs.size() && false ) {
	    std::cerr<<" ERROR: Ordering in MASS section wrong. Must be idx continuous."<<std::endl;
	    return false;
	  }
	  mhs    .push_back( (double) massval );
	  mhs_fit.push_back( ( massIdx < 0 ? 0 : 1 ) );  // use onyl for testing/eff. smear etc.
	}
	//} else {
	//std::cerr<<" Input line "<<line<<" not recopgnized."<<std::endl;
	//return false;
	//}
      }
    }
  }
  
  fclose(configFile);
  return true;
  
}



bool readConfigCardNuissances(TString configCardName, int numCats,
			      std::vector<RooAbsReal*>& nsigcat,
			      std::vector<RooAbsReal*>& nsigcat_2nd,
			      std::vector<RooAbsReal*>& nsigcatModel,
			      std::vector<RooAbsReal*>& nuissances,
			      std::vector<RooAbsReal*>& finalnorm,
			      std::vector<RooAbsReal*>& finalnorm_2nd,
			      std::vector<RooAbsReal*>& finalnormModel,
			      std::vector<TString> catnames) {
  
  nsigcat.resize(numCats);
  nsigcat_2nd.resize(numCats);
  nsigcatModel.resize(numCats);
  nuissances.resize(0);
  
  //TString cardName = TString("../")+configCardName;  
  TString cardName = configCardName;  
  FILE* configFile = fopen(cardName.Data(),"r");
  if ( !configFile ) {
    std::cerr<<" Inputfile "<<configCardName<<" not found."<<std::endl;
    return false;
  }
  
  char line[1000];
  
  std::cout<<" Reading model from file "<<configCardName<<"..."<<std::endl;
  
  while (fgets(line,1000,configFile)) {
    // must be some nuissance or normalization stuff...
    int indNuis = -1;
    float startVal = -1.;
    float minVal = -1.;
    float maxVal = -1.;
    char name[30];
    if( sscanf(line,"NUIS %d %s %f %f %f",&indNuis, name, &startVal, &minVal, &maxVal) ) {
      if (indNuis > numCats ) return false;
      if (indNuis != (int) nuissances.size()) return false;
      nuissances.push_back( new RooRealVar(name,"",startVal,minVal,maxVal) );
      ((RooRealVar*) nuissances[nuissances.size()-1])->setConstant();
    } else {
      // normalization stuff...
      int catInd = -1;
      int nParms = -1;
      char formula[30];
      char parList[50];
      if( sscanf(line,"NSIG %d %d %s %s",&catInd,&nParms, formula, parList) ) {
	// need to transform the parList...
	RooArgList theList;
	RooArgList theList_2nd;
	RooArgList theListModel;
	if( !strcmp(formula,"DEFAULT") ) {	 
	  theList.add( *(finalnorm[catInd]) );
	  theList_2nd.add( *(finalnorm_2nd[catInd]) );
	  theListModel.add( *(finalnormModel[catInd]) );
	  nsigcat     [catInd] = new RooFormulaVar(TString::Format("nsig%s",catnames.at(catInd).Data()),"","@0",theList);
	  nsigcat_2nd [catInd] = new RooFormulaVar(TString::Format("nsig2nd%s",catnames.at(catInd).Data()),"","@0",theList_2nd);
	  nsigcatModel[catInd] = new RooFormulaVar(TString::Format("nsigModel%s",catnames.at(catInd).Data()),"","@0",theListModel);
	} else {	
	  std::string parStr(parList);
	  size_t lastPos = 0;
	  while ( parStr.find_last_of(",") != std::string::npos ) {
	    size_t pos = parStr.find_first_of(",");
	    std::string thisPar = parStr.substr(lastPos,pos);
	    int idx = -1;
	    if ( sscanf(thisPar.c_str(),"NUIS(%d)",&idx) ) {
	      theList.add(* (nuissances[idx]) );
	      theList_2nd.add(* (nuissances[idx]) );
	      theListModel.add(* (nuissances[idx]) );
	    } else if ( sscanf(thisPar.c_str(),"NOMINAL(%d)",&idx) ) {
	      theList.add(* (finalnorm[idx]) );
	      theList_2nd.add(* (finalnorm_2nd[idx]) );
	      theListModel.add(* (finalnormModel[idx]) );
	    } else {
	      std::cerr<<" ERROR: Cannot recognize Par Type "<<thisPar.c_str()<<"."<<std::endl;
	      return false;
	    }
	    parStr.replace(pos,1,"X");
	    lastPos = pos+1;
	  }
	  // last parameter
	  std::string lastPar = parStr.substr(lastPos);
	  int idx = -1;
	  if ( sscanf(lastPar.c_str(),"NUIS(%d)",&idx) ) {
	    theList.add(* (nuissances[idx]) );
	    theList_2nd.add(* (nuissances[idx]) );
	    theListModel.add(* (nuissances[idx]) );
	  } else if ( sscanf(lastPar.c_str(),"NOMINAL(%d)",&idx) ) {
	    theList.add(* (finalnorm[idx]) );
	    theList_2nd.add(* (finalnorm_2nd[idx]) );
	    theListModel.add(* (finalnormModel[idx]) );
	  }  else {
	    std::cerr<<" ERROR: Cannot recognize Par Type "<<lastPar.c_str()<<"."<<std::endl;
	    return false;
	  }
	  RooFormulaVar* tempForm = new RooFormulaVar(TString::Format("nsig%s",catnames.at(catInd).Data()),"",formula,theList);
	  nsigcat[catInd] = tempForm;
	  RooFormulaVar* tempForm2nd = new RooFormulaVar(TString::Format("nsig2nd%s",catnames.at(catInd).Data()),"",formula,theList);
	  nsigcat_2nd[catInd] = tempForm2nd;
	  RooFormulaVar* tempFormModel = new RooFormulaVar(TString::Format("nsigModel%s",catnames.at(catInd).Data()),"",formula,theListModel);
	  nsigcatModel[catInd] = tempFormModel;
	}
      }   
    }
  }

  fclose(configFile);
  return true;
}


bool resetStartValues(std::map<TString,RooRealVar*>& fitparms,
		      std::map<TString,float>&       startvals) {
  
  for(std::map<TString,RooRealVar*>::iterator pair = fitparms.begin(); pair != fitparms.end(); ++pair) {
    std::map<TString,float>::iterator pair2 = startvals.find(pair->first);
    if( pair2 != startvals.end() )
      pair->second->setVal(pair2->second);
    else
      return false;
  }

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
				     TString& nosmearmodname,
				     TString& nosmeartreename,
				     double& massmax, double& massmin,
				     double& fitmassmax, double& fitmassmin,
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
  computeMVAvar    = false;
  correctIDMVAvar  = false;
  correctSigEoEvar = false;

  std::cout<<" Reading paramater from file "<<fileName<<"...";
  
  while (fgets(line,1000,configFile)) {

    if(line[0] == '#') continue;

    char name[400];
    int catIdx = -1;
    char catName[400];
    //char theCat[400];
    float massval = -1.;
    if ( sscanf(line,"PROJECTDIR %s", name ) ) projectDir = TString(name);
    else if  ( sscanf(line,"MODNAME %s", name) ) modname = TString(name);
    else if  ( sscanf(line,"NOSMEARMOD %s", name) ) nosmearmodname = TString(name);
    else if  ( sscanf(line,"TREENAME %s", name) ) treename = TString(name);
    else if  ( sscanf(line,"NOSMEARTREE %s", name) ) nosmeartreename = TString(name);
    else if( sscanf(line,"LUMI %f",&massval) ) theLumi = (double) massval;
    else if( sscanf(line,"MINMSS %f",&massval) ) massmin = (double) massval;
    else if( sscanf(line,"MAXMSS %f",&massval) ) massmax = (double) massval;
    else if( sscanf(line,"FITMSSMIN %f",&massval) ) fitmassmin = (double) massval;
    else if( sscanf(line,"FITMSSMAX %f",&massval) ) fitmassmax = (double) massval;
    else if  ( sscanf(line,"AUXCAT %d %s",&catIdx,  catName) ) {
      // parsing the ctegory definition like:
      // AUXCAT	0	masscut		" mass>100.0 && mass<180. "
      std::string totLine = line;
      int startCat = totLine.find_first_of("\"");
      int endCat   = totLine.find_last_of("\"");
      std::string catLine = totLine.substr(startCat+1, endCat-startCat-1);
      auxCats.insert(std::pair<TString,TString>(TString(catName),TString(catLine.c_str())));
    }
    else if  ( sscanf(line,"COMPUTEMVA ON %s", name) ) { mvaWeightFile = TString(name); computeMVAvar = true; }
    else if  ( sscanf(line,"CORRECTIDMVA ON %s", name) ) { idmvaCorrFile = TString(name); correctIDMVAvar = true; }
    else if  ( sscanf(line,"CORRECTSIGEOE ON %s", name) ) { sigeoeCorrFile = TString(name); correctSigEoEvar = true; }
    else if  ( sscanf(line,"MVADEFFILE %s", name) ) { mvaDefFile = TString(name);}
    else if  ( sscanf(line,"ANACAT %d %s",&catIdx, catName) ) {
      // ANACAT	0	cat0		" basecut && !vbfcut && baseline0 "		0.005432		Bern/5
      std::string totLine = line;
      int startCat = totLine.find_first_of("\"");
      int endCat   = totLine.find_last_of("\"");
      std::string catLine = totLine.substr(startCat+1, endCat-startCat-1);
      
      // erase starting empties...
      //unsigned int fPos = catLine.find_first_not_of(" ");
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
    else if  ( sscanf(line,"BASECAT %s", name) ) {
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
