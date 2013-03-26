#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
#include <iostream>
#include "TH1.h"
#include "TFile.h"
#include "TLatex.h"
#include "TPostScript.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TF1.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TRandom.h"
#include "TNtuple.h"
#include "TMath.h"
#include "TTimeStamp.h"

#include "Math/DistFuncMathCore.h"

#include <sstream>


#include "RooWorkspace.h"
#include "RooAbsPdf.h"
#include "RooRealIntegral.h"
#include "RooGaussian.h"
#include "RooAddPdf.h"
#include "RooAddition.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooSimultaneous.h"
#include "RooRealVar.h"
#include "RooConstVar.h"
#include "RooCategory.h"
#include "RooFitResult.h"
#include "RooFFTConvPdf.h"
#include "RooPlot.h"
#include "RooMinimizer.h"
#include "RooNLLVar.h"
#include "RooCachedPdf.h"
#include "TTimeStamp.h"
#include "RooExponential.h"
#include "RooGenericPdf.h"
#include "RooChebychev.h"
#include "RooBernstein.h"
#include "RooRandom.h"
#include "RooPolynomial.h"
#include "TSystem.h"
#include "RooExtendPdf.h"
#include "RooChi2Var.h"

#endif

// ===========================================================================================

#include "FstModels.h"
#include "FstModelDef.h"

// ===========================================================================================

double doFit_Chi2(int fitMode       = 1,
		  int fitOrder      = 2,
		  TString catName   = "hgg_8TeV_2013moriond_bdt0",
		  bool doSave       = false,
		  bool verbose      = false,
		  bool blinded      = false,
		  // MODIFY ME PLEASE !!!
		  TString modelFile = "/home/fabstoec/Hgg_BGBiasStudy/MoriondHggWorkspaces/databkg/bkgdatawithfit.root"
		  ) {


  // MODIFY ME PLEASE!!!
  TString dataWSname = "wbkg";

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetLabelColor(1, "XYZ");
  gStyle->SetLabelFont(42, "XYZ");
  gStyle->SetLabelOffset(0.007, "XYZ");
  gStyle->SetLabelSize(0.035, "XYZ");

  gStyle->SetTitleColor(1, "XYZ");
  gStyle->SetTitleFont(42, "XYZ");
  gStyle->SetTitleOffset(1., "XYZ");
  gStyle->SetTitleSize(0.04, "XYZ");

  gStyle->SetPalette(1);

  int dummyNDF = 0;
  TString modelName = "";
  
  fstHgg_getNameAndNDF(fitMode,dummyNDF, modelName, dummyNDF);
  
  TString labelRight      = TString::Format("%d%s",fitOrder,modelName.Data());
  TString labelRightFinal = TString::Format("  Model = %d%s",fitOrder,modelName.Data());
  
  // load the Data workspace
  TFile* dataFile = TFile::Open(modelFile.Data());
  
  RooWorkspace* wdata = (RooWorkspace*) dataFile->Get(dataWSname.Data());
  
  // ======================================================
  // Shut TF** up

  if ( !verbose ) {
    
    RooMsgService::instance().setSilentMode(true);
    RooMsgService::instance().getStream(0).removeTopic(RooFit::Caching);
    RooMsgService::instance().getStream(1).removeTopic(RooFit::Caching);
    RooMsgService::instance().getStream(0).removeTopic(RooFit::Minimization);
    RooMsgService::instance().getStream(1).removeTopic(RooFit::Minimization);
    RooMsgService::instance().getStream(0).removeTopic(RooFit::Plotting);
    RooMsgService::instance().getStream(1).removeTopic(RooFit::Plotting);
    RooMsgService::instance().getStream(0).removeTopic(RooFit::Fitting);
    RooMsgService::instance().getStream(1).removeTopic(RooFit::Fitting);
    RooMsgService::instance().getStream(0).removeTopic(RooFit::Eval);
    RooMsgService::instance().getStream(1).removeTopic(RooFit::Eval);
    RooMsgService::instance().getStream(0).removeTopic(RooFit::Integration);
    RooMsgService::instance().getStream(1).removeTopic(RooFit::Integration);
    
    RooMsgService::instance().setStreamStatus(0,false);
    RooMsgService::instance().setStreamStatus(1,false);
    
  }
  // -----------------------------------------------------
  // read in the Workspaces  

  // read the data (binned for quickness...)
  RooDataSet*   dataCat       = NULL;
  RooDataHist*  dataCatHist   = NULL;
  
  // MODIFY ME PLEASE: make sure this is THE MASS name...
  RooRealVar* mass     = wdata->var("CMS_hgg_mass");

  // ================================
  // MODIFY ME PLEASE !!!!
  dataCat     = (RooDataSet*) wdata->data(TString::Format("data_%s",catName.Data()));
  dataCatHist = dataCat->binnedClone(TString::Format("datahist_%s",catName.Data()));
  // ================================  
  
  
  fstModel* bgfitcat = NULL;

  bgfitcat = fstHgg_getModel(mass,fitMode,fitOrder,0,"fit");
  
  if ( !bgfitcat ) {
    std::cout<<"  Fitting Mode = "<<fitMode<<" not implemented."<<std::endl;
    return 0.;
  }  
  
  RooAbsData* blindData = NULL;
  if( blinded )
    blindData = dataCat->reduce(RooFit::Cut("CMS_hgg_mass < 110. || CMS_hgg_mass > 150."), RooFit::Name("blinded"));
  else
    blindData = dataCat;
  
  RooRealVar*    toynorm = new RooRealVar("toynorm","",0.,20000.);
  toynorm->setVal(dataCat->sumEntries()); // goodn starting value
  toynorm->removeRange();
  
  // produce the truth model here...
  
  RooExtendPdf* toypdf = new RooExtendPdf("toypdf","",*(bgfitcat->getPdf()),*toynorm);
  
  RooChi2Var*   toyChi2 = (RooChi2Var*) toypdf->createChi2(*dataCatHist,RooFit::DataError(RooAbsData::Poisson));
  RooMinimizer* toyMinim2 = new RooMinimizer(*toyChi2);
  toyMinim2->fit("m");
  
  RooNLLVar*    toyNll = (RooNLLVar*) toypdf->createNLL(*dataCat,RooFit::Extended(true));
  RooMinimizer* toyMinim = new RooMinimizer(*toyNll);
  toyMinim->fit("m");
  
  double theChi2Val = toyChi2->getVal();
  double theNllVal  = toyNll->getVal();
  
  double band_1sigma[160];
  double band_2sigma[160];
  double band_X[160];
  
  TGraph* g_1sigma = NULL;
  TGraph* g_2sigma = NULL;
  
  RooPlot* frame = mass->frame(RooFit::Bins(80));
  //RooPlot* frame = mass->frame();
  dataCat->plotOn(frame,RooFit::Invisible());
  
  if (doSave) {
    
    // compute the bands in 1 GeV steps
    RooRealVar* nlim = new RooRealVar("nlim","",0.,1000.);
    nlim->removeRange();
    
    nlim->setVal(dataCat->sumEntries()/10.);
    
    for( int iBin = 0; iBin < 80; ++iBin ) {
      
      band_X[iBin]       = 100.+(double)iBin + 0.5;
      band_X[159 - iBin] = 100.+(double)iBin + 0.5;
      
      mass->setRange("nlimrange",100.+(double)iBin, 101+(double)iBin);
      RooExtendPdf epdf("epdf","",*(bgfitcat->getPdf()),*nlim,"nlimrange");
      
      if(toyNll)   delete toyNll;
      if(toyMinim) delete toyMinim;
      
      toyNll = (RooNLLVar*) epdf.createNLL(*dataCatHist,RooFit::Extended(true));
      toyMinim = new RooMinimizer(*toyNll);
      toyMinim->fit("m");
      toyMinim->minos();
      
      band_1sigma[iBin]       = nlim->getVal() +    nlim->getAsymErrorLo();
      band_2sigma[iBin]       = nlim->getVal() + 2.*nlim->getAsymErrorLo();
      
      band_1sigma[159 - iBin] = nlim->getVal() +    nlim->getAsymErrorHi();
      band_2sigma[159 - iBin] = nlim->getVal() + 2.*nlim->getAsymErrorHi();
      
    }
    
    g_1sigma = new TGraph(160,band_X,band_1sigma);
    g_2sigma = new TGraph(160,band_X,band_2sigma);
    
    g_1sigma->SetFillColor(kYellow);
    g_2sigma->SetFillColor(kGreen);
    
    frame->addObject( g_2sigma, "F");
    frame->addObject( g_1sigma, "F");
    
    //     g_2sigma->Draw("F");
    //     g_1sigma->Draw("F");
    
  }
  
  toypdf->plotOn(frame,RooFit::LineColor(kRed));
  blindData->plotOn(frame);
  frame->SetMinimum(1e-5);
  
  
  TCanvas* can = new TCanvas();
  can->cd();
  frame->Draw();
  
  TLatex* textModel=new TLatex(0.9,0.95,labelRightFinal.Data());
  textModel->SetNDC();
  textModel->SetTextAlign(33);
  textModel->SetTextFont(42);
  textModel->SetTextSize(0.04);// dflt=28
  textModel->Draw();
  
  TLatex* textCat=new TLatex(0.1,0.95,catName.Data());
  textCat->SetNDC();
  textCat->SetTextAlign(13);
  textCat->SetTextFont(42);
  textCat->SetTextSize(0.04);// dflt=28
  textCat->Draw();
  
  
  if( doSave ) 
    can->SaveAs(TString::Format("plots/truth_%s_%s.pdf",catName.Data(),labelRight.Data()).Data());
  
  std::cout<<TString::Format("\tC:\t%s\tM:\t%s\tNLL = %.2f",catName.Data(),labelRight.Data(),theNllVal).Data()<<std::endl;
  std::cout<<TString::Format("\tC:\t%s\tM:\t%s\tChi2/NDF = %.2f/%d (%.2f)",catName.Data(),labelRight.Data(),theChi2Val,mass->getBins(),theChi2Val/mass->getBins()).Data()<<std::endl;
  
  return theNllVal;
}

