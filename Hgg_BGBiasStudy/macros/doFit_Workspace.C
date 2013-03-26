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

#endif

#include "FstModels.h"
#include "FstModelDef.h"


// ===========================================================================================

double doFit_Workspace(int genMode       = 1,   // mode fot the GEN funtion
		       int genOrder      = 1,
		       int fitMode       = 1,
		       int fitOrder      = 1,
		       TString catName   = "hgg_8TeV_2013moriond_bdt0",
		       bool plotBias     = true,
		       bool verbose      = false,
		       int numMass       = 41,
		       // MODIFY ME PLEASE !!!
		       TString modelDir = "/home/fabstoec/Hgg_BGBiasStudy/MoriondHggWorkspaces"
		       ) {
  
  // load the Data workspace
  // MODIFY ME PLEASE !!!
  TString dataFileName = modelDir + TString("/databkg/bkgdatawithfit.root");
  TFile* dataFile = TFile::Open(dataFileName.Data());

  // MODIFY ME PLEASE !!!
  RooWorkspace* wdata = (RooWorkspace*) dataFile->Get("wbkg");

  // load the signal WS
  // MODIFY ME PLEASE !!!
  TString sigFileName = modelDir + TString("/model/ubersignalmodel.root");
  TFile* sigFile = TFile::Open(sigFileName.Data());
  // MODIFY ME PLEASE !!!
  RooWorkspace* wsig = (RooWorkspace*) sigFile->Get("wsig");
  
  // ==============================================

  Double_t _nbg_fit    [numMass];
  Double_t _nbg_fit_err[numMass];
  Double_t _nbg        [numMass];
  Double_t _nbg_err    [numMass];
  
  Double_t _bias       [numMass];

  // ======================================================
  // Shut TF** up
  
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

  // -----------------------------------------------------
  // read in the Workspaces  

  // read the data (binned for quickness...)
  RooAbsData*  dataCat   = NULL;
  dataCat = (RooDataHist*) wdata->data(TString::Format("databinned_%s",catName.Data()));

  // open signal worksapces...
  std::vector<TString> procNames;
  procNames.push_back("ggh");
  procNames.push_back("vbf");
  procNames.push_back("wzh");
  procNames.push_back("tth");
  // open all the worksapces  
  
  RooAbsPdf*     sigpdfcat       = NULL;
  RooRealVar* MH       = NULL;  // hypothesis mass  
  RooRealVar* mass     = NULL;  // di-photon mass

  double* leftBound  = new double[numMass];
  double* rightBound = new double[numMass];

  RooFitResult* bgfitres  = NULL;

  if (verbose)
    std::cout<<" CheckingCat "<<catName<<std::endl;
  
  RooArgList compList;
  // ============================================================
  // MODIFY ME PLEASE (maybe a lot to do... sorry..., until the next '=======')
  for( unsigned int iProc = 0; iProc < procNames.size(); ++iProc ){
    if ( !MH   ) MH   = wsig->var("MH");
    if ( !mass ) mass = wsig->var("CMS_hgg_mass");
    
    // get the normalization Variable
    RooFormulaVar* absNormVar = (RooFormulaVar*) wsig->function(TString::Format("hggpdfsmrel_%s_%s_norm",catName.Data(),procNames[iProc].Data()).Data());
    // get the PDF
    RooAbsPdf*    sigPdfRel = wsig->pdf(TString::Format("hggpdfsmrel_%s_%s",catName.Data(),procNames[iProc].Data()).Data());
    RooExtendPdf* sigPdf    = new RooExtendPdf(TString::Format("hggpdfsmrel_%s_%s_ext",catName.Data(),procNames[iProc].Data()).Data(),"",*sigPdfRel,*absNormVar);      
    compList.add(*sigPdf);
  }
  
  sigpdfcat = new RooAddPdf(TString::Format("hggpdfsmrel_tmp_%s",catName.Data()).Data(),"",compList);
  // ============================================================  

  // compute FWHM for each mass
  for(int iMass = 0; iMass < numMass; ++iMass) {
    double theMass = 110. + iMass * (150. - 110.)/ (double) ( numMass - 1 );
    MH->setVal(theMass);
    
    mass->setRange("plotrange",theMass-15., theMass+10.);
    mass->setBins(100,"plotrange");
    
    
    RooPlot* frame1 = mass->frame(RooFit::Range("plotrange"),RooFit::Bins(100));
    sigpdfcat    ->plotOn(frame1);
    RooCurve* nomcurve = new RooCurve(*frame1->getCurve(frame1->nameOf(0)));
    
    // =============== DO FWHM COMPUTATION =======================
    // 1. find the maximum
    double thePos    = theMass - 8.;
    double maxVal    = -10.;
    double maxPos    = -10.;
    double binSize  = 0.001;
    while ( thePos < theMass + 8. ) {
      double theVal = nomcurve->interpolate(thePos);
      if ( theVal > maxVal ) {
	maxPos = thePos;
	maxVal = theVal;
      } //else break; // we're going down again...
      thePos += binSize;
    }
    // 2. ... so half-max-value is:
    double HMValue = maxVal/2.;
    
    // 3. find the left end of the arrow
    double minDiff = 10000000.;
    double leftEnd = theMass - 8.;
    thePos    = theMass - 8.;
    
    while ( thePos < maxPos ) {
      double theVal = nomcurve->interpolate(thePos);
      double theDiff = TMath::Abs(theVal - HMValue);
      if (theDiff < minDiff) {
	minDiff = theDiff;
	leftEnd = thePos;
      } else break; // gone to far...
      thePos += binSize;
    }
    
    // 4. find the right end of the arrow
    minDiff = 100.;
    double rightEnd = theMass;
    thePos          = maxPos;
    
    while ( thePos < theMass + 8. ) {
      double theVal = nomcurve->interpolate(thePos);
      double theDiff = TMath::Abs(theVal - HMValue);
      if (theDiff < minDiff) {
	minDiff = theDiff;
	rightEnd = thePos;
      } else break; // getting worse... so break;
      thePos += binSize;
    }
    
    leftBound[iMass] = leftEnd;
    rightBound[iMass] = rightEnd;
    
    if( verbose ) {
      std::cout<<" FWHM @"<<theMass<<":  ["<<leftEnd<<" - "<<rightEnd<<"]"<<"  -> "<<rightEnd-leftEnd<<std::endl;
    }
  }
  
  // ============= IMPORTANT ==================
  // chose which model to use here
  
  fstModel* bgpdfcat = NULL;
  fstModel* bgfitcat = NULL;
  
  bgpdfcat = fstHgg_getModel(mass,genMode,genOrder,0,"truth");      
  bgfitcat = fstHgg_getModel(mass,fitMode,fitOrder,0,"fit");      
    
  if( !bgpdfcat || !bgfitcat ) {
    std::cerr<<"ERROR: Could not generate some model.."<<std::endl;
    return -100.;
  }

  RooAbsData* dataHCat = NULL;

  RooRealVar*    toynorm = new RooRealVar("toynorm","",0.,20000.);
  toynorm->setVal(dataCat->sumEntries()); // goodn starting value
  toynorm->removeRange();
  
  
  // produce the truth model here...
  RooExtendPdf* toypdf = new RooExtendPdf("toypdf","",*(bgpdfcat->getPdf()),*toynorm);
  RooNLLVar* toyNll = (RooNLLVar*) toypdf->createNLL(*dataCat,RooFit::Extended(true));
  RooMinimizer* toyMinim = new RooMinimizer(*toyNll);
  toyMinim->fit("m");
  
  double* numBG     = new double[numMass];
  double* numBGerr  = new double[numMass];
  RooRealVar* nbglim = new RooRealVar("nbglim","",0.,1000.);
  nbglim->removeRange();
  RooMinimizer* bgminim = NULL;
  RooNLLVar*    bgnll   = NULL;

  // store max bias & max-bias-mass
  double maxBias     = 0.;
  double maxBiasMass = 0.;


  for( int iMass = 0; iMass < numMass; ++iMass ) {
    
    if ( bgfitres ) delete bgfitres;
    
    mass->setRange("errrange",leftBound[iMass],rightBound[iMass]);
    if (iMass == 0 )
      nbglim->setVal( dataCat->sumEntries()/10. );
    else
      nbglim->setVal( nbglim->getVal()*0.8 );
    RooExtendPdf ebgpdf("ebgpdf","",*(bgpdfcat->getPdf()),*nbglim,"errrange");
    
    if( bgnll ) delete bgnll;
    if( bgminim ) delete bgminim;
    
    bgnll = (RooNLLVar*) ebgpdf.createNLL(*dataCat,RooFit::Extended(true));
    bgminim = new RooMinimizer(*bgnll);
    bgminim->fit("m");
    bgfitres = bgminim->lastMinuitFit();
    if ( !bgfitres->status() ) {
      bgminim->minos();
      bgfitres = bgminim->lastMinuitFit();
    }
    if ( bgfitres->status() ) {
      std::cerr<<" ERROR: Fit for truth model not converged. Status = "<<bgfitres->status()<<std::endl;
      return -1.;
    }
    numBG   [iMass] = nbglim->getVal();

    if ( nbglim->hasAsymError() ) {
      numBGerr[iMass] = ( nbglim->getAsymErrorHi() - nbglim->getAsymErrorLo() )/2.;
      if (nbglim->getAsymErrorHi() == 0)
 	numBGerr[iMass] = - nbglim->getAsymErrorLo();
      if (nbglim->getAsymErrorLo() == 0)
 	numBGerr[iMass] =  nbglim->getAsymErrorHi();
    } else
      numBGerr[iMass] = nbglim->getError();
    
    double theMass = 110. + iMass * (150. - 110.)/ (double) ( numMass - 1 );
    
    if(verbose)
      std::cout<<" TRUTH N_bg (mh="<<theMass<<") : "<<numBG[iMass]<<" +- "<<numBGerr[iMass]<<std::endl;    
  }

  if( bgnll ) delete bgnll;
  if( bgminim ) delete bgminim;
  
  // fit also the BG only fit model, to have decent starting values
  delete toyMinim;
  RooNLLVar* bgonlyNll = (RooNLLVar*) bgfitcat->getPdf()->createNLL(*dataCat);
  toyMinim = new RooMinimizer(*bgonlyNll);
  toyMinim->fit("m");
  bgfitcat->storeValues();

  // produce the total S+B Pdf
  RooRealVar*    snorm  = new RooRealVar("snorm","",0.,1000.);
  snorm->removeRange();
  RooRealVar*    bgnorm = new RooRealVar("bgnorm","",0.,20000.);
  bgnorm->setVal(dataCat->sumEntries()); // goodn starting value
  bgnorm->removeRange();
    
  RooNLLVar*    theNll = NULL;
  RooMinimizer* minim  = NULL;

  RooRealVar* nlim = new RooRealVar("nlim","",0.,1000.);
  nlim->removeRange();  

  if( dataHCat ) delete dataHCat;      
  dataHCat    = toypdf->generateBinned(*mass,RooFit::Name("toydata"),RooFit::NumEvents(dataCat->sumEntries()),RooFit::Asimov());
  
  // ---------------------------------------------------
  // fit the BG only to get good starting values
  if (theNll) delete theNll;
  if (minim)  delete minim;
  
  theNll = (RooNLLVar*) bgfitcat->getPdf()->createNLL(*dataHCat);
  minim  = new RooMinimizer(*theNll);
  minim->fit("m");
  //bgfitcat->getPdf()->fitTo(*dataHCat);
  bgfitcat->storeValues();
  // ---------------------------------------------------
  
  nlim->setVal( toynorm->getVal()/20.);
  // loop over the mass points
  for(int iMass = 0; iMass < numMass; ++iMass){
    
    if (iMass > 0) nlim->setVal( nlim->getVal() *0.8 ); 
    
    // we have the toy data, now fit with the BG only model (get differntial BG uncertanioty
    double theMass = 110. + iMass * (150. - 110.)/ (double) ( numMass - 1 );
    
    MH->setVal(theMass);
    mass->setRange("errrange",leftBound[iMass],rightBound[iMass]);
    RooExtendPdf epdf("epdf","",*(bgfitcat->getPdf()),*nlim,"errrange");
    
    if (theNll) delete theNll;
    if (minim)  delete minim;
    theNll = NULL;
    minim  = NULL;

    theNll = (RooNLLVar*) epdf.createNLL(*dataHCat,RooFit::Extended(true));
    
    // reset the BG value to the initial fit to the toy
    double nlimStart = nlim->getVal();
    do {
      if ( bgfitres ) delete bgfitres;
      //if ( theNll   ) delete theNll;
      if ( minim    ) delete minim;
      nlimStart *= 0.95;
      nlim->setVal(nlimStart);
      nlim->setConstant(false);
      bgfitcat ->resetValues();
      minim  = new RooMinimizer(*theNll);
      minim->fit("m");
      bgfitres = minim->lastMinuitFit();
      if ( !bgfitres->status() ) {
	minim->minos();
	bgfitres = minim->lastMinuitFit();
      }
      
    } while (nlim->getVal() < 0. && !bgfitres->status() );
    
    _nbg    [iMass] = numBG[iMass];
    _nbg_err[iMass] = numBGerr[iMass];
    
    _nbg_fit    [iMass] = nlim->getVal();
    
    bool hasAsym = false;
    
    if ( nlim->hasAsymError() ) {
      
      hasAsym = true;
      
      _nbg_fit_err[iMass] = ( nlim->getAsymErrorHi() - nlim->getAsymErrorLo() )/2.;
      if (nlim->getAsymErrorHi() == 0) {
	//std::cout<<" Error Hi = 0 !"<<std::endl;
	_nbg_fit_err[iMass] = - nlim->getAsymErrorLo();
      }
      if (nlim->getAsymErrorLo() == 0) {
	//std::cout<<" Error Lo = 0 !"<<std::endl;
	_nbg_fit_err[iMass] =  nlim->getAsymErrorHi();
      }
    }
    else
      _nbg_fit_err[iMass] = nlim->getError();
    

    _bias[iMass] = (_nbg_fit[iMass]-_nbg[iMass])/_nbg_fit_err[iMass];

    if( TMath::Abs(_bias[iMass]) > TMath::Abs(maxBias) ) {
      maxBias     = _bias[iMass];
      maxBiasMass = theMass;
    }

    if( verbose ) {
      std::cout<<"            Bias (@"<<theMass<<")= "<<_bias[iMass]<<std::endl; 
      std::cout<<"                  N_TRUTH  ["<<leftBound[iMass]<<" - "<<rightBound[iMass]<<"] = "<<_nbg[iMass]<<" +- "<<_nbg_err[iMass]<<std::endl;
      std::cout<<"                  N_FIT    ["<<leftBound[iMass]<<" - "<<rightBound[iMass]<<"] = "<<_nbg_fit[iMass]<<" +- "<<_nbg_fit_err[iMass]<<"  asym ? "<<hasAsym<<std::endl;
      if( nlim->hasAsymError() )
	std::cout<<"                  N_FIT    ["<<leftBound[iMass]<<" - "<<rightBound[iMass]<<"] = "<<_nbg_fit[iMass]<<" + "<<nlim->getAsymErrorHi()<<" - "<<nlim->getAsymErrorLo()<<std::endl;
    }

    if (bgfitres->status() ) {
      std::cout<<" *** Skipping toy #"<<std::endl;
      std::cout<<"     Reason: BG did not vonverge for m = "<<theMass<<" with NBG = "<<nlim->getVal()<<"  status = "<<bgfitres->status()<<std::endl;
      
      std::cout<<"            Bias (@"<<theMass<<")= "<<_bias[iMass]<<std::endl;
      std::cout<<"                  N_TRUTH  ["<<leftBound[iMass]<<" - "<<rightBound[iMass]<<"] = "<<_nbg[iMass]<<" +- "<<_nbg_err[iMass]<<std::endl;
      std::cout<<"                  N_FIT    ["<<leftBound[iMass]<<" - "<<rightBound[iMass]<<"] = "<<_nbg_fit[iMass]<<" +- "<<_nbg_fit_err[iMass]<<std::endl;
      
      break; // skip this toy...
    }    
  }
  
  std::cout<<" Max Bias @mh = "<<maxBiasMass<<" with value "<<maxBias<<std::endl;

  // plot the entire story...
  if( plotBias ) {
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

    
    stringstream yaxisTitle;
    stringstream labelRight;
    stringstream labelRight2;
    stringstream labelRightTotal;
    
    labelRight <<genOrder;
    
    labelRightTotal <<"  Truth = ";
    
    switch(genMode) {
    case 1:
      labelRight << "Exp";
      break;
    case 2:
      labelRight << "Pow";
      break;
    case 3:
      labelRight << "Pol";
      break;
    case 4:
      labelRight << "Lau";
      break;
    case 5:
      labelRight << "Ber";
      break;
    case 6:
      labelRight << "Pol(+Gauss)";
      break;
    case 7:
      labelRight.str("");
      labelRight << "Truth = Exp(-"<<genOrder<<"Pol)";
      break;
    case 8:
      labelRight.str("");
      labelRight << "Truth = Pow*Exp(-"<<genOrder<<"Pol)";
      break;
    case 9:
      labelRight << "Greg";
      break;
    default:
      std::cout<<" Mode not implemenyed."<<std::endl;
      return -1.;
    }
    
    labelRightTotal<<labelRight.str().c_str()<< "   Fit = ";
    
    labelRight2 << fitOrder;
    
    switch(fitMode) {
    case 1:
      labelRight2 << "Exp";
      break;
    case 2:
      labelRight2 << "Pow";
      break;
    case 3:
      labelRight2 << "Pol";
      break;
    case 4:
      labelRight2 << "Lau";
      break;
    case 5:
      labelRight2 << "Ber";
      break;
    case 6:
      labelRight2 << "Pol(+Gauss)";
      break;
    case 7:
      labelRight2 << "Exp(-"<<fitOrder<<"Pol)";
      break;
    case 8:
      labelRight2 << "Pow*Exp(-"<<fitOrder<<"Pol)";
      break;
    case 9:
      labelRight2 << "Greg";
      break;
    default:
      std::cout<<" Mode not implemenyed."<<std::endl;
      return -1.;
    }
    
    labelRightTotal<<labelRight2.str().c_str();

    yaxisTitle << "Median( (N_{fit}^{FWHM}-N_{truth}^{FWHM})/#Delta N_{fit}^{FWHM} )";
  
    Double_t x_mass[2*numMass];
    double massStepSize = (150.-110.)/(double) (numMass-1);
  
    for( int iMass = 0; iMass < numMass; ++iMass ) {
      x_mass[            iMass] = 110. + (double) iMass * massStepSize;
      x_mass[2*numMass-1-iMass] = 110. + (double) iMass * massStepSize;
    }

    Double_t* mean = new Double_t[numMass];
  
    double theValue = TMath::Abs(maxBias);
    theValue = TMath::Max(theValue, 0.2); // plot until 0.2 at least

    TGraph* g_mean    = NULL;

    for(int iMass = 0; iMass < numMass; ++iMass)
      mean[iMass] = _bias[iMass];

    g_mean    = new TGraph(numMass,x_mass,mean);
    
    g_mean->SetLineColor(kRed);
    g_mean->SetLineWidth(2.);
    
    TH1D* dummy = new TH1D( TString::Format("dummy_%d%d%d%d",genMode,genOrder,fitMode,fitOrder).Data() ,"",1,105.,155.);
    dummy->SetBinContent(1,0.);
    dummy->GetYaxis()->SetLimits(-1.1*theValue, 1.1*theValue);
    dummy->SetMaximum( 1.1*theValue);
    dummy->SetMinimum(-1.1*theValue);
    
    dummy->SetLineStyle(kDashed);

    dummy->GetXaxis()->SetTitle("m_{H}   [GeV]");  
    dummy->GetYaxis()->SetTitle(yaxisTitle.str().c_str());

    // 5-times-smaller region
    double fiveX[4];
    double fiveY[4];
    
    fiveX[0] = 100;
    fiveX[1] = 100;
    fiveX[2] = 160;
    fiveX[3] = 160;
    
    fiveY[0] = -0.2;
    fiveY[1] =  0.2;
    fiveY[2] =  0.2;
    fiveY[3] = -0.2;
    
    TGraph* fiveG = new TGraph(4,fiveX,fiveY);
    fiveG->SetFillColor(kBlue-10);
    
    TCanvas* can = new TCanvas();
    dummy ->Draw();
    
    fiveG->Draw("F");
    dummy->Draw("AXIS SAME");
    
    g_mean   ->Draw("LP");
    
    dummy->Draw("SAME");


    TLatex* textModel=new TLatex(0.9,0.95,labelRightTotal.str().c_str());
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

    TLatex* textAsimov=new TLatex(0.15,0.85,"ASIMOV TOY");
    textAsimov->SetNDC();
    textAsimov->SetTextAlign(13);
    textAsimov->SetTextFont(42);
    textAsimov->SetTextSize(0.04);// dflt=28
    textAsimov->Draw();

    can->SaveAs(TString::Format("plots/biasAsimov_%s_%s_vs_%s.pdf",catName.Data(),labelRight.str().c_str(),labelRight2.str().c_str()).Data());

  }

  return maxBias;
}

