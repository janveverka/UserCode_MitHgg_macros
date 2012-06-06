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

#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooWorkspace.h"
#include "RooConstVar.h"
#include "TEfficiency.h"
#include "RooConstVar.h"
#include "RooAddition.h"

#include "TString.h"

//---------------------------------------------------------------------------------------------------------------------------
// Include the header file with the model constans
#include "modelConstants_8TeV.h"

using namespace RooFit;

//---------------------------------------------------------------------------------------------------------------------------
// Fwd helper function declaration... implementation at the end
bool validateInput(int iProc, int iCat, int nProcs, int nCats);

bool readConfigCard(TString configCardName, 
		    std::vector<TString>& procNames, std::vector<bool>& procOn,
		    TString& inDir, TString& outDir, 
		    TString& wsPrefix,
		    int& numProcs, int& numCats, int& numModels, int& numMasses,
		    bool*& catIsOn,
		    double& maxmass, double& minmass,
		    std::vector<double>& mhs,
		    std::vector<double>& smearingv,
		    std::vector<TString>& catNames,
		    int*& _right, int*& _wrong,
		    float**& _meanStart, float**& _sigmaStart, float**& _fracStart,
		    float**& _meanStartWrong, float**& _sigmaStartWrong, float**& _fracStartWrong);

bool readConfigCardNuissances(TString configCardName, int numCats,
			      std::vector<RooAbsReal*>& nsigcat,
			      std::vector<RooAbsReal*>& nuissances,
			      std::vector<RooAbsReal*>& finalnorm,
			      std::vector<TString>      catnames);

RooAbsPdf* generateMultiGaussian(RooRealVar* mass, 
				 RooRealVar* nomMass,
				 int numComp,
				 RooGaussian**& g,
				 TString procName, TString quali,
				 std::map<TString,RooRealVar*>& fitparms,
				 std::map<TString,float>&       startvals,
				 float* startMean, float* startSigma, float* startFrac, bool fix=false);

bool resetStartValues(std::map<TString,RooRealVar*>& fitparms,
		      std::map<TString,float>&       startvals);

RooDataSet *cwdset(RooDataSet *indata, RooRealVar *mvar, RooRealVar *wvar, RooRealVar *pidvar, TString name, Double_t weightscale, int filterproc);
void *appendcwd(RooDataSet *outdata, RooDataSet *indata, RooRealVar *mvar, RooRealVar *wvar, Double_t weightscale);


//---------------------------------------------------------------------------------------------------------------------------


//void createSignalModels(TString configCardName="mettag.config", modeltype model = STANDARDMODEL) {
void fithmassCard(TString configCardName="template.config", bool fitonly = false) {
  
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

  //   modelParmsMap->insert(std::pair<TString,double*>(TString("zhxsec"),zhxsec));
  //   modelParmsMap->insert(std::pair<TString,double*>(TString("whxsec"),whxsec));

  modelParmsMap->insert(std::pair<TString,double*>(TString("wzhxsec"),wzhxsec));

  modelParmsMap->insert(std::pair<TString,double*>(TString("tthxsec"),tthxsec));
  modelParmsMap->insert(std::pair<TString,double*>(TString("sm4gghxsec"),sm4gghxsec));

  //----------------------------------------------------------------------------
  // sum the cross-sections properly
  std::vector<TString>                      sumRuleNames;
  std::vector<std::vector<const char*>* >   sumRuleComps;

  std::vector<const char*> myVec;
  
//   std::vector<const char*> sumNames_WZH;
//   sumNames_WZH.push_back("whxsec");
//   sumNames_WZH.push_back("zhxsec");
//   sumRuleNames.push_back("wzhxsec");
//   sumRuleComps.push_back(&sumNames_WZH);

  std::vector<const char*> sumNames_SM;
  //sumNames_SM.push_back("whxsec");
  //sumNames_SM.push_back("zhxsec");
  sumNames_SM.push_back("wzhxsec");
  sumNames_SM.push_back("gghxsec");
  sumNames_SM.push_back("vbfxsec");
  sumNames_SM.push_back("tthxsec");
  sumRuleNames.push_back("smxsec");
  sumRuleComps.push_back(&sumNames_SM);
  
  std::vector<const char*> sumNames_FF;
  //sumNames_FF.push_back("whxsec");
  //sumNames_FF.push_back("zhxsec");
  sumNames_FF.push_back("wzhxsec");
  sumNames_FF.push_back("vbfxsec");
  sumRuleNames.push_back("ffxsec");
  sumRuleComps.push_back(&sumNames_FF);
  
  std::vector<const char*> sumNames_SM4;
  //sumNames_SM4.push_back("whxsec");
  //sumNames_SM4.push_back("zhxsec");
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

  int numProcs  = -1;
  int numCats   = -1;
  int numModels = -1;
  int numMasses = -1;

  int* _right = NULL;
  int* _wrong = NULL;

  bool* catIsOn = NULL;

  float** _meanStartRight = NULL;
  float** _sigmaStartRight = NULL;
  float** _fracStartRight = NULL;

  float** _meanStartWrong = NULL;
  float** _sigmaStartWrong = NULL;
  float** _fracStartWrong = NULL;


  std::vector<TString> procnames;
  std::vector<bool>    procon;
  TString inputDir;
  TString outputDir;
  TString wsPrefix;

  std::vector<double> smearingv;

  // define the masses we're testing to interpolate the models between
  std::vector<double> mhs;
//   for(int iMass=0; iMass<5; ++iMass)
//     mhs.push_back(110. + iMass * 10.);
  //  setting up the common stuff
  double massmax = -1.;
  double massmin = -1.;

  std::vector<TString> catnamesbase;

  bool status = readConfigCard(configCardName, 
			       procnames,
			       procon,
			       inputDir,
			       outputDir,
			       wsPrefix,
			       numProcs, numCats, numModels, numMasses,
			       catIsOn,
			       massmax, massmin, mhs,
			       smearingv,
			       catnamesbase,
			       _right, _wrong,
			       _meanStartRight, _sigmaStartRight, _fracStartRight,
			       _meanStartWrong, _sigmaStartWrong, _fracStartWrong);

  if(!status) {
    std::cerr<<" ERROR when readin input card "<<configCardName<<"."<<std::endl;
    return;
  }

  // ------------------------------------------------------------------
  // General ROOT setup
  gROOT->Macro("MitStyle.C");
  gStyle->SetErrorX(0); 
  gStyle->SetOptStat(0);
  gROOT->ForceStyle();  
  gStyle->SetOptStat(1110);

  RooMsgService::instance().getStream(1).removeTopic(RooFit::Caching);
  RooMsgService::instance().getStream(0).removeTopic(RooFit::Minimization);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::Minimization);
  RooMsgService::instance().getStream(0).removeTopic(RooFit::Plotting);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::Plotting);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::Fitting);
  RooMsgService::instance().getStream(0).removeTopic(RooFit::Eval);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::Eval);
  // ------------------------------------------------------------------

  // loading the Workspaces... should possibly configurable... ? FIX-ME
  TFile* fdata  = new TFile(inputDir+TString::Format("/%s-mclow.root",wsPrefix.Data()));
  TFile* fdata2 = new TFile(inputDir+TString::Format("/%s-mchigh.root",wsPrefix.Data()));
  
  if( !fdata || !fdata2 ) {
    std::cerr<<" ERROR: Could not load input file "<<TString(inputDir+( !fdata ? TString::Format("%s-mclow.root",wsPrefix.Data()) : TString::Format("%s-mchigh.root",wsPrefix.Data()) ))<<"."<<std::endl;
    return;
  }
  
  // get the workspaces from the file ... configurable ? FIX-ME...
  RooWorkspace* win1 = (RooWorkspace*)fdata->Get("cms_hgg_workspace_mclow");
  RooWorkspace* win2 = (RooWorkspace*)fdata2->Get("cms_hgg_workspace_mchigh");  
  RooWorkspace *win = win1;

  if( !win1 || !win2 ) {
    std::cerr<<" ERROR: Could not load workspace "<<( !win1 ? TString("cms_hgg_workspace_mclow") : TString("cms_hgg_workspace_mchigh") )<<"."<<std::endl;
    return;
  }
  

  // output workspace ... again: name configuirable ? FIX-ME
  RooWorkspace* wOut = new RooWorkspace("wsig","");

  RooRealVar *hmass = win->var("CMS_hgg_mass");
  hmass->setRange(massmin,massmax);
  hmass->setBins( (int) (4.0*(massmax-massmin)) );
  hmass->SetTitle("m_{#gamma#gamma}");
  hmass->setUnit("GeV");  
  RooRealVar mnom("MH","m_{h}",110.,massmin,massmax,"GeV");
  mnom.setConstant();
  
//   std::vector<TString> catnamesbase;
//   for(int iCat=0; iCat < numCats; ++iCat) {
//     stringstream pSS;
//     pSS << "cat"<<iCat;
//     catnamesbase.push_back(pSS.str().c_str());
//   }


  //---RooRealVar to load dataset or used for calculation---
  RooRealVar *weight = new RooRealVar("weight","",1.0);//this seems to be useless variable
  weight->removeRange();
  RooRealVar *procidx = win->var("procidx");  
  RooRealVar *IntLumi = win->var("IntLumi");
  
  // maps for the histFuncs and the Additions
  std::map<TString,RooHistFunc*> histFuncMap;
  std::map<TString,RooAddition*> addFuncMap;
  
  int numHists = modelParmsMap->size();
  TH1D** parmHists     = new TH1D*       [numHists];
  RooDataHist** parmDataHists = new RooDataHist*[numHists];
  RooHistFunc** parmHistFuncs = new RooHistFunc*[numHists];
  int iPair = 0;
  for(std::map<TString,double*>::iterator pair=modelParmsMap->begin(); pair!=modelParmsMap->end(); ++pair, ++iPair) {
    TString name = pair->first;
    double* data = pair->second;
    stringstream histName;
    histName << "h" <<name.Data()<<"s";
    parmHists[iPair] = new TH1D(histName.str().c_str(),"",numsmpoints,smmasses[0]-0.25,smmasses[numsmpoints-1]+0.25);
    for(int ipoint=0; ipoint<numsmpoints; ++ipoint)
      parmHists[iPair]->Fill(smmasses[ipoint],data[ipoint]);
    histName.str("");
    histName << "d"<<name.Data()<<"s";
    parmDataHists[iPair] = new RooDataHist(histName.str().c_str(),"",RooArgList(mnom),parmHists[iPair]);
    histName.str("");
    histName << "f"<<name.Data()<<"s";
    parmHistFuncs[iPair] = new RooHistFunc(histName.str().c_str(),"",RooArgList(mnom),(*parmDataHists[iPair]),1);
    histFuncMap.insert(std::pair<TString,RooHistFunc*>(name,parmHistFuncs[iPair]));
  }
  
  // create the summed up functions
  for(unsigned int iSum=0; iSum<sumRuleNames.size(); iSum++) {
    std::vector<const char*>* comps = sumRuleComps[iSum];
    TString sumName                 = sumRuleNames[iSum];
    stringstream addName;
    addName << "f" << sumName.Data() << "s";
    RooArgList sumList;
    for(unsigned int iComp=0; iComp<comps->size(); ++iComp) {
      RooHistFunc* tempFunc = histFuncMap.find(  (TString)((*comps)[iComp]) )->second;
      if(! tempFunc) {
	std::cerr<<" Cannot find histfunc with name "<<TString( comps->at(iComp) )<<"."<<std::endl;
	return;
      }
      sumList.add( *tempFunc );
    }
    RooAddition* tempAdd = new RooAddition(addName.str().c_str(),"",sumList);
    addFuncMap.insert(std::pair<TString,RooAddition*>(sumName,tempAdd));
  }
  
  std::vector<RooAbsReal*> procxseclist;
  for(unsigned int iComp=0; iComp < procxseclist_names->size(); ++iComp) {
    stringstream pSS;
    pSS <<  TString(procxseclist_names->at(iComp));
    RooAbsReal* tempAbs = addFuncMap.find((TString) pSS.str().c_str())->second;
    if( !tempAbs ) tempAbs = histFuncMap.find((TString) pSS.str().c_str())->second;
    if( !tempAbs ) {
      std::cerr<<" Cannot find addfunc with name "<<pSS.str().c_str()<<"."<<std::endl;
      return;
    }
    procxseclist.push_back(tempAbs);
  }

  // loop over all processes that are on
  for(int iProc = 0; iProc < numProcs; ++iProc) {
    if( !procon[iProc] ) continue;
    
    // change to output directory
    int filterproc = iProc;
    std::cout<<" Changing directory to "<<TString::Format(outputDir+TString("/%s"),procnames.at(filterproc).Data())<<std::endl;
    if ( !gSystem->cd(TString::Format(outputDir+TString("/%s"),procnames.at(filterproc).Data())) ) {
      std::cerr<<" ERROR: Directory "<<TString::Format(outputDir+TString("/%s"),procnames.at(filterproc).Data())<<" does not exists."<<std::endl;
      std::cerr<<"        Skipping process "<<procnames[filterproc]<<"."<<std::endl;
      continue;
    }


    // full Cat names including the process name
    std::vector<TString> catnames;
    for(int iCat = 0; iCat < numCats; ++iCat) {
      stringstream pSS;
      pSS << catnamesbase[iCat].Data() << "_" << procnames[iProc];
      catnames.push_back(pSS.str().c_str());
    }


    // collection for the fitparameters
    std::map<TString,RooRealVar*>* fitparms  = new std::map<TString,RooRealVar*>[numCats];
    std::map<TString,float>      * startvals = new std::map<TString,float>      [numCats];

    // generate the RooAbsPdf s for the right and wrong Vtx hypothesis    
    RooAbsPdf** combh       = new RooAbsPdf*[numCats];    // right hypothesis
    RooAbsPdf** combhwrong  = new RooAbsPdf*[numCats];    // ... well ...

    RooGaussian*** gRight   = new RooGaussian**[numCats];
    RooGaussian*** gWrong   = new RooGaussian**[numCats];

    for(int iCat=0; iCat < numCats; ++iCat) {
      combh[iCat]      = generateMultiGaussian(hmass,&mnom,_right[iProc*numCats+iCat],gRight[iCat],procnames[iProc],"rVtx",fitparms[iCat],startvals[iCat],
					       _meanStartRight[iProc*numCats+iCat],_sigmaStartRight[iProc*numCats+iCat],_fracStartRight[iProc*numCats+iCat]);
      combhwrong[iCat] = generateMultiGaussian(hmass,&mnom,_wrong[iProc*numCats+iCat],gWrong[iCat],procnames[iProc],"wVtx",fitparms[iCat],startvals[iCat],
					       _meanStartWrong[iProc*numCats+iCat],_sigmaStartWrong[iProc*numCats+iCat],_fracStartWrong[iProc*numCats+iCat], false);
    }

    //parameters calculated from mc events rather than fitting: fracright and effacc
    RooRealVar** fracright = new RooRealVar*[numCats];
    RooRealVar** effacc    = new RooRealVar*[numCats];

    // generate the combined PDFs
    RooAddPdf** combhvtx   = new RooAddPdf* [numCats];

    for(int iCat = 0; iCat < numCats; ++iCat) {
      stringstream pSS;
      pSS << "fracright_" << procnames[iProc] << "_" << iCat;      
      fracright[iCat] = new RooRealVar(pSS.str().c_str(),"fracright",0.9,0.0,1.0);

      fitparms[iCat].insert(std::pair<TString,RooRealVar*>(TString(pSS.str().c_str()),fracright[iCat]));
      startvals[iCat].insert(std::pair<TString,float>(TString(pSS.str().c_str()),0.9));

      pSS.str("");
      pSS << "effacc_" << procnames[iProc] << "_" << iCat;      
      effacc[iCat]    = new RooRealVar(pSS.str().c_str(),"effacc",   0.9,0.0,1.0);
      
      fitparms[iCat].insert(std::pair<TString,RooRealVar*>(TString(pSS.str().c_str()),effacc[iCat]));
      startvals[iCat].insert(std::pair<TString,float>(TString(pSS.str().c_str()),0.9));
      
      pSS.str("");
      pSS << "combhvtx_" << procnames[iProc] << "_" << iCat;            
      combhvtx[iCat] = new RooAddPdf(pSS.str().c_str(),"combhvtx",RooArgList( (*combh[iCat]),(*combhwrong[iCat]) ), RooArgList( (*fracright[iCat]) ));

    }

    //define histograms to keep track of fit parameters for each category each mass in the step of 10 GeV
    TH1F ***fitparmhists = new TH1F**[numCats];
    for (UInt_t iCat=0; iCat<catnames.size(); ++iCat) {
      fitparmhists[iCat] = new TH1F*[fitparms[iCat].size()];
      unsigned int iparm = 0;
      for(std::map<TString,RooRealVar*>::iterator pair = fitparms[iCat].begin(); pair != fitparms[iCat].end(); ++pair) {	
	TString histname = TString("hist")+ pair->first + catnames.at(iCat);
	fitparmhists[iCat][iparm] = new TH1F(histname,histname,(int) numMasses,mhs[0]-5.,mhs[numMasses-1]+5.);
	iparm++;
      }
    }
    
    //---define 2d hist to fill the fit status for each mass and each category---
    TH2F *histfitstatus = new TH2F("histfitstatus","histfitstatus",(int) numMasses,mhs[0]-5.,mhs[numMasses-1]+5.,catnames.size(),-0.5,catnames.size()-0.5);//0 fit; 1 not fit
    TH2F *histfitstatuswrong = new TH2F("histfitstatuswrong","histfitstatuswrong",(int) numMasses,mhs[0]-5.,mhs[numMasses-1]+5.,catnames.size(),-0.5,catnames.size()-0.5);//0 fit; 1 not fit
    
    //---defined the datasets for test---
    std::vector<RooDataSet*> testdsets;
    std::vector<RooDataSet*> testdsetswrong;
    std::vector<RooDataSet*> testdsetsall;
    
    //---add canvas for plot---
    TCanvas* rightCan = new TCanvas("rightCan","right vertex",numMasses*250,numCats*250); 
    rightCan->Divide(numMasses,numCats);
    TCanvas* wrongCan = new TCanvas("wrongCan","wrong vertex",numMasses*250,numCats*250); 
    wrongCan->Divide(numMasses,numCats);
    TCanvas* allCan = new TCanvas("allCan","combined",numMasses*250,numCats*250); 
    allCan->Divide(numMasses, numCats);
    
    TCanvas* dummy = new TCanvas();
    
    for (UInt_t iCat=0; iCat<catnames.size(); ++iCat) {

      if ( !catIsOn[iProc*numCats+iCat] ) continue;  // don't do fits in this Cat

      for (UInt_t i=0; i<mhs.size(); ++i) {
	//workspace

	std::cout<<" MH = "<<mhs[i]<<std::endl;

	if (mhs[i]<=120) win = win1;
	else win = win2;
	//mass points string
	std::stringstream numstringstr;
	numstringstr<<mhs.at(i);
	TString numstring(numstringstr.str());
	TString vbfstring = numstring;
	//get dataset for this category for this mass point
	
	TString dataRightName;
	TString dataWrongName;
	TString dataAllName;
	
	if(filterproc!=4){
	  dataRightName=TString::Format("sig_%s_mass_m",procnames[filterproc].Data())+numstring+TString("_rv_")+catnamesbase.at(iCat);
	  dataWrongName=TString::Format("sig_%s_mass_m",procnames[filterproc].Data())+numstring+TString("_wv_")+catnamesbase.at(iCat);
	  dataAllName=TString::Format("sig_%s_mass_m",procnames[filterproc].Data())+numstring+TString("_")+catnamesbase.at(iCat);
	}
	
	if(filterproc==4){
	  dataRightName=TString("sig_mass_m")+numstring+TString("_rv_")+catnamesbase.at(iCat);
	  dataWrongName=TString("sig_mass_m")+numstring+TString("_wv_")+catnamesbase.at(iCat);
	  dataAllName=TString("sig_mass_m")+numstring+TString("_")+catnamesbase.at(iCat);  
	}
      
	//right vertex
	RooDataSet *mcsigdata = (RooDataSet*)win->data(dataRightName.Data());

// 	RooPlot* frame = hmass->frame(); 
// 	mcsigdata->plotOn(frame);
// 	frame->Draw();
// 	return;

	if( !mcsigdata ) {
	  std::cerr<<" Could not find dataset "<<dataRightName.Data()<<" in workspace."<<std::endl;
	  return;
	}
	const double weightscale = static_cast<double>(mcsigdata->numEntries())/mcsigdata->sumEntries();
	RooDataSet *mcsigwdata = cwdset(mcsigdata,hmass,weight,procidx,TString("mcsigwdata") + numstring+catnames.at(iCat),weightscale,filterproc);
	
	//wrong vertex
	RooDataSet *mcsigwrongdata = (RooDataSet*)win->data(dataWrongName.Data());
	if( !mcsigwrongdata ) {
	  std::cerr<<" Could not find dataset "<<dataWrongName.Data()<<" in workspace."<<std::endl;
	  return;
	}
	RooDataSet *mcsigwrongwdata = cwdset(mcsigwrongdata,hmass,weight,procidx,TString("mcsigwrongwdata") + numstring+catnames.at(iCat),weightscale,filterproc);
	
	//combined
	RooDataSet *mcsigalldata = (RooDataSet*)win->data(dataAllName.Data());
	if( !mcsigalldata ) {
	  std::cerr<<" Could not find dataset "<<dataAllName.Data()<<" in workspace."<<std::endl;
	  return;
	}
	RooDataSet *mcsigallwdata = cwdset(mcsigalldata,hmass,weight,procidx,TString("mcsigallwdata") + numstring+catnames.at(iCat),weightscale,filterproc);
	
	//printf("sigdata = %i, sigwrong = %i, sigall = %i\n",mcsigdata!=0, mcsigwrongdata!=0,mcsigalldata!=0);
	printf("right = %5f, wrong = %5f, all = %5f\n",mcsigdata->sumEntries(),mcsigwrongdata->sumEntries(), mcsigalldata->sumEntries());
	printf("right = %5f, wrong = %5f, all = %5f\n",mcsigwdata->sumEntries(),mcsigwrongwdata->sumEntries(), mcsigallwdata->sumEntries());

	//set higgs mass and fit range

	mnom.setVal(mhs.at(i));

	hmass->setRange("higgsrange",massmin,massmax);
	hmass->setRange("plotrange",TMath::Max(massmin,mhs.at(i)-30.0),TMath::Min(massmax,mhs.at(i)+20.0));
	//hmass->setRange("plotrange",100.,160.);

	// reset the fir parameters to the starting values
	if( ! resetStartValues(fitparms[iCat],startvals[iCat]) ) {
	  std::cerr<<" ERROR: Could not reset parameters to starting values."<<std::endl;
	  return;
	}
	
	//fit
	RooFitResult *fitres = 0;
	RooFitResult *fitreswrong = 0;
	
	RooAbsPdf* rightpdf = combh[iCat];
	RooAbsPdf* wrongpdf = combhwrong[iCat];
	RooAbsPdf* allpdf   = combhvtx[iCat];

	printf("FITSTART:  proc:%d mass point:%d cat:%d\n",filterproc,(int) mhs.at(i),iCat); 

	//mnom.setVal(110.);
	fitres      = rightpdf     ->fitTo(*mcsigwdata,     Strategy(0),Minimizer("Minuit2",""),Minos(kFALSE),SumW2Error(kFALSE), Save(kTRUE),NumCPU(8));
	fitreswrong = wrongpdf     ->fitTo(*mcsigwrongwdata,Strategy(0),Minimizer("Minuit2",""),Minos(kFALSE),SumW2Error(kFALSE), Save(kTRUE),NumCPU(8));    

	//fitres      = rightpdf     ->fitTo(*mcsigwdata,     SumW2Error(kFALSE), Save(kTRUE) );
	//fitreswrong = wrongpdf     ->fitTo(*mcsigwrongwdata,SumW2Error(kFALSE), Save(kTRUE) );    

	//mnom.setVal(mhs.at(i));
	//std::cout<<" *** nom set to "<<mhs.at(i)<<std::endl;

	printf("FITRES:  proc:%d mass point:%d cat:%d fit status right vertex:%d wrong vertex:%d \n",filterproc,mhs.at(i),iCat,fitres->status(),fitreswrong->status()); 

	if( !fitres || !fitreswrong ) {
	  std::cerr<<" ERROR: Did not get any fitresult."<<std::endl;
	  return;
	}

	//compute acceptance*efficiency 
	RooRealVar *mtotalxsec = win->var(TString::Format("XSBR_%s_%d",procnames[iProc].Data(),(int) mhs.at(i)));
	printf("lumi = %5f, xsec = %5f\n",IntLumi->getVal(),mtotalxsec->getVal());
	double eaccnum = mcsigwdata->sumEntries()+mcsigwrongwdata->sumEntries();
	double eaccden = IntLumi->getVal()*mtotalxsec->getVal()*weightscale;//reover the weightscale
	double eacc = eaccnum/eaccden;
	printf("eacc = %5f, eaccnum = %5f, eaccden = %5f\n",eacc,eaccnum,eaccden);
	double eaccerrlo = TEfficiency::ClopperPearson(Int_t(eaccden), Int_t(eaccnum), 0.683, kFALSE) - eacc;
	double eaccerrhi = TEfficiency::ClopperPearson(Int_t(eaccden), Int_t(eaccnum), 0.683, kTRUE) - eacc;
	printf("eff done\n");
	
	//compute the fraction of right vertex
	//double fright = mcsigwdata->sumEntries()/(mcsigwdata->sumEntries()+mcsigwrongwdata->sumEntries());
	double frightnum = mcsigwdata->sumEntries();
	double frightden = mcsigwdata->sumEntries()+mcsigwrongwdata->sumEntries();
	double fright = frightnum/frightden;
	double frighterrlo = TEfficiency::ClopperPearson(Int_t(frightden), Int_t(frightnum), 0.683, kFALSE) - fright;
	double frighterrhi = TEfficiency::ClopperPearson(Int_t(frightden), Int_t(frightnum), 0.683, kTRUE) - fright;
	
	fracright[iCat]->setVal(fright);
	fracright[iCat]->setAsymError(frighterrlo,frighterrhi);
	
	effacc[iCat]->setVal(eacc);
	effacc[iCat]->setAsymError(eaccerrlo,eaccerrhi);
      
	//TCanvas *chfit = new TCanvas;
	TString plotname = TString("rightvtx") + numstring + catnames.at(iCat) + TString(".eps");//e.g. rightvtx110cat0_ggh.eps      
	RooPlot *hplot = hmass->frame(Bins(100),Range("plotrange"));
	

	mcsigwdata->plotOn(hplot);

	unsigned int* colors = new unsigned int[3];
	colors[0] = kOrange;
	colors[1] = kMagenta;
	colors[2] = kRed;

	for(int iComp=0; iComp < _right[iProc*numCats+iCat]; ++iComp)
	  rightpdf->plotOn(hplot,Components( *(gRight[iCat][iComp]) ),LineColor(colors[iComp]),Range("higgsrange"),NormRange("higgsrange"));
	
	rightpdf->plotOn(hplot,RooFit::LineColor(kBlue),Range("higgsrange"),NormRange("higgsrange"));  
	hplot->SetTitle("");
	hplot->Draw();  
	//chfit->SaveAs(plotname);
	
	//TCanvas *chfitwrong = new TCanvas;
	TString plotnamewrong = TString("wrongvtx") + numstring + catnames.at(iCat) + TString(".eps");            
	RooPlot *hplotwrong = hmass->frame(Bins(40),Range("plotrange"));
	mcsigwrongwdata->plotOn(hplotwrong);

	for(int iComp=0; iComp < _wrong[iProc*numCats+iCat]; ++iComp)
	  wrongpdf->plotOn(hplotwrong,Components( *(gWrong[iCat][iComp]) ),LineColor(colors[iComp]),Range("higgsrange"),NormRange("higgsrange"));

	wrongpdf->plotOn(hplotwrong,RooFit::LineColor(kBlue),Range("higgsrange"),NormRange("higgsrange"));  
	hplotwrong->SetTitle("");
	hplotwrong->Draw();        
	//chfitwrong->SaveAs(plotnamewrong);
	
	//TCanvas *chfitall = new TCanvas;
	TString plotnameall = TString("allvtx") + numstring + catnames.at(iCat) + TString(".eps");                  
	RooPlot *hplotall = hmass->frame(Bins(100),Range("plotrange"));
	mcsigallwdata->plotOn(hplotall);
	allpdf->plotOn(hplotall,RooFit::LineColor(kBlue),Range("higgsrange"),NormRange("higgsrange"));  
	hplotall->SetTitle("");
	hplotall->Draw();   
	//chfitall->SaveAs(plotnameall);
	
	TString plotnameOverview = TString("rightvtx") +  TString(".eps");
	TString plotnamewrongOverview = TString("wrongvtx") + TString(".eps");  
	TString plotnameallOverview = TString("allvtx") + TString(".eps");  
	
	rightCan->cd(iCat*numMasses+i+1);
	hplot->Draw();  
	wrongCan->cd(iCat*numMasses+i+1);
	hplotwrong->Draw(); 
	allCan->cd(iCat*numMasses+i+1);
	hplotall->Draw(); 

	rightCan->SaveAs(plotnameOverview);
	wrongCan->SaveAs(plotnamewrongOverview);
	allCan->SaveAs(plotnameallOverview); 
	printf ("Sum of weights: right = %5f, wrong = %5f, all = %5f, right+wrong = %5f\n",mcsigwdata->sumEntries(),mcsigwrongwdata->sumEntries(), mcsigallwdata->sumEntries(),mcsigwdata->sumEntries()+mcsigwrongwdata->sumEntries());  
	printf ("Number of entries: right = %5f, wrong = %5f, all = %5f, right+wrong = %5f\n",mcsigdata->sumEntries(),mcsigwrongdata->sumEntries(), mcsigalldata->sumEntries(),mcsigdata->sumEntries()+mcsigwrongdata->sumEntries());  
	//printf("data weights = %5e\n",mcsigwdata->sumEntries());
	//printf("mass = %i, status = %i, statuswrong = %i\n",mhs.at(i),fitres->status(),fitreswrong->status());

	// -----------------------------------------

	// need to sort the width parameters...
	double* sortedSigma = new double[_right[iProc*numCats+iCat]];
	for(int iComp=0; iComp < _right[iProc*numCats+iCat]; ++iComp) {
	  stringstream pSS;
	  pSS << "sigma"<<iComp<<"_"<<procnames[iProc]<<"_rVtx";
	  std::map<TString,RooRealVar*>::iterator pair = fitparms[iCat].find(TString(pSS.str().c_str()));
	  if( pair == fitparms[iCat].end() ) {
	    std::cerr<<" ERROR: Could not find fir parameter with name "<<pSS.str().c_str()<<std::endl;
	    return;
	  }
	  int insertPosition = iComp;
	  // also take the absolute values...
	  double thisVal = TMath::Abs(pair->second->getVal());
	  for(int jComp = 0; jComp < iComp; ++jComp) {
	    if(sortedSigma[jComp] < thisVal) {
	      for(int kComp = iComp-1; kComp >= jComp; --kComp)
		sortedSigma[kComp+1] = sortedSigma[kComp];
	      insertPosition = jComp;
	      break;
	    }
	  }
	  sortedSigma[insertPosition] = thisVal;
	}
	for(int iComp=0; iComp < _right[iProc*numCats+iCat]; ++iComp) {
	  stringstream pSS;
	  pSS << "sigma"<<iComp<<"_"<<procnames[iProc]<<"_rVtx";
	  std::map<TString,RooRealVar*>::iterator pair = fitparms[iCat].find(TString(pSS.str().c_str()));
	  if( pair == fitparms[iCat].end() ) {
	    std::cerr<<" ERROR: Could not find fir parameter with name "<<pSS.str().c_str()<<std::endl;
	    return;
	  }
	  pair->second->setVal(sortedSigma[iComp]);
	}

	for(int iComp=0; iComp < _wrong[iProc*numCats+iCat]; ++iComp) {
	  stringstream pSS;
	  pSS << "sigma"<<iComp<<"_"<<procnames[iProc]<<"_wVtx";
	  std::map<TString,RooRealVar*>::iterator pair = fitparms[iCat].find(TString(pSS.str().c_str()));
	  if( pair == fitparms[iCat].end() ) {
	    std::cerr<<" ERROR: Could not find fir parameter with name "<<pSS.str().c_str()<<std::endl;
	    return;
	  }
	  int insertPosition = iComp;
	  // also take the absolute values...
	  double thisVal = TMath::Abs(pair->second->getVal());
	  for(int jComp = 0; jComp < iComp; ++jComp) {
	    if(sortedSigma[jComp] < thisVal) {
	      for(int kComp = iComp-1; kComp >= jComp; --kComp)
		sortedSigma[kComp+1] = sortedSigma[kComp];
	      insertPosition = jComp;
	      break;
	    }
	  }
	  sortedSigma[insertPosition] = thisVal;
	}
	for(int iComp=0; iComp < _wrong[iProc*numCats+iCat]; ++iComp) {
	  stringstream pSS;
	  pSS << "sigma"<<iComp<<"_"<<procnames[iProc]<<"_wVtx";
	  std::map<TString,RooRealVar*>::iterator pair = fitparms[iCat].find(TString(pSS.str().c_str()));
	  if( pair == fitparms[iCat].end() ) {
	    std::cerr<<" ERROR: Could not find fir parameter with name "<<pSS.str().c_str()<<std::endl;
	    return;
	  }
	  pair->second->setVal(sortedSigma[iComp]);
	}

	

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
	histfitstatus->Fill(mhs.at(i),iCat,fitres->status());
        histfitstatuswrong->Fill(mhs.at(i),iCat,fitreswrong->status());

	dummy->cd();
      }
    }

    // ---------------------------------------------------------------
    // prepare RooDataSet and RooHistFunc arrays
    RooDataHist*** fitparmdatas = new RooDataHist**[numCats];
    
    std::map<TString,RooHistFunc*>* fitparmfuncs = new std::map<TString,RooHistFunc*>[numCats];
    
    for(int iCat=0; iCat<numCats; ++iCat) {

      fitparmdatas[iCat] = new RooDataHist*[fitparms[iCat].size()];
      //fitparmfuncs[iCat] = new RooHistFunc*[fitparms[iCat].size()];

      // loop over all fitparams
      // -----------------------------------------
      unsigned int iParm=0;
      for(std::map<TString,RooRealVar*>::iterator pair = fitparms[iCat].begin(); pair != fitparms[iCat].end(); ++pair) {
	TString dataname=TString("data") + pair->first + catnames.at(iCat);
	TString funcname=TString("func") + pair->first + catnames.at(iCat);
	fitparmdatas[iCat][iParm] = new RooDataHist(dataname,dataname,RooArgList(mnom),fitparmhists[iCat][iParm]);
	fitparmfuncs[iCat].insert(std::pair<TString,RooHistFunc*>(pair->first, new RooHistFunc(funcname,funcname,RooArgList(mnom),*fitparmdatas[iCat][iParm],1)));
	
	// do the plottiong here if requests
	if( !fitonly ) {
	  TString plotname = TString("func") + pair->first + catnames.at(iCat) + TString(".eps");//e.g funcsigma1cat0_ggh.eps
	  TCanvas *cfunctest = new TCanvas;
	  //  fitparmhists[iCat*nparms + iparm]->Draw();
	  RooPlot *hploteffacc = mnom.frame(Bins(100),Range(105,155));
	  fitparmdatas[iCat][iParm]->plotOn(hploteffacc);  
	  //  fitparmfuncs[iCat][iparm]->plotOn(hploteffacc,RooFit::LineColor(kBlue));  
	  fitparmfuncs[iCat].find(pair->first)->second->plotOn(hploteffacc,RooFit::LineColor(kBlue));  
	  hploteffacc->SetTitle("");
	  hploteffacc->GetYaxis()->SetTitle(pair->first.Data());
	  hploteffacc->Draw(); 
	  cfunctest->SaveAs(plotname);
	  delete cfunctest;
	}
	iParm++;

      }
    }
    // ---------------------------------------------------------------
    //---draw fit status---
    TCanvas *cfitstatus = new TCanvas;
    histfitstatus->Draw("COL");
    cfitstatus->SaveAs("fitstatus.eps");
    
    TCanvas *cfitstatuswrong = new TCanvas;
    histfitstatuswrong->Draw("COL");
    cfitstatuswrong->SaveAs("fitstatuswrong.eps");
    
    //---nuissance---
    //the nuissance and the systematics are redone in combine
    RooRealVar nuissancedeltafracright("CMS_hgg_nuissancedeltafracright","",1.0,0.1,10.0);
    nuissancedeltafracright.setConstant();
    
    //---define final pdfs in each category---
    RooConstVar   **smears               = new RooConstVar*[catnames.size()];
    RooRealVar    **nuissancedeltasmears = new RooRealVar*[catnames.size()];
    RooRealVar    **nuissancedeltams     = new RooRealVar*[catnames.size()];
    RooFormulaVar **smearmods            = new RooFormulaVar*[catnames.size()];

    RooFormulaVar*** meanslides   = new RooFormulaVar**[catnames.size()];
    RooFormulaVar*** sigmaslides  = new RooFormulaVar**[catnames.size()];
    RooFormulaVar*** wmeanslides  = new RooFormulaVar**[catnames.size()];
    RooFormulaVar*** wsigmaslides = new RooFormulaVar**[catnames.size()];
    
    RooGaussian*** gslides  = new RooGaussian**[catnames.size()];
    RooGaussian*** wgslides = new RooGaussian**[catnames.size()];

    RooAddPdf**   combhslides      = new RooAddPdf*[catnames.size()];
    RooAddPdf**   combhwrongslides = new RooAddPdf*[catnames.size()];

    RooFormulaVar** fracrightmodslides = new RooFormulaVar*[catnames.size()];  

    RooAddPdf  **combhvtxslides = new RooAddPdf*[catnames.size()];
    RooAbsPdf  **finalpdfslides = new RooAbsPdf*[catnames.size()];
    std::vector<RooAbsReal*> finalnormslides;
    finalnormslides.resize(catnames.size());


    for(int iCat=0; iCat < numCats; ++iCat) {
      smears[iCat] = new RooConstVar(TString("smear")+catnames.at(iCat),"",smearingv.at(iCat));
      nuissancedeltasmears[iCat] = new RooRealVar(TString("CMS_hgg_nuissancedeltasmear")+catnamesbase.at(iCat),"",0.0, -smearingv.at(iCat),smearingv.at(iCat));
      nuissancedeltasmears[iCat]->setConstant();
      nuissancedeltams[iCat] = new RooRealVar(TString("CMS_hgg_nuissancedeltam")+catnamesbase.at(iCat),"",0.0,-5.0,5.0);
      nuissancedeltams[iCat]->setConstant();
      smearmods[iCat] = new RooFormulaVar(TString("smearmod")+catnames.at(iCat),"","@0*(@1 + @2)",RooArgList(mnom,*smears[iCat],*nuissancedeltasmears[iCat]));

      // ------------------------------------------------------------------------------------------------------------------
      // Right Vertex Models
      meanslides [iCat] = new RooFormulaVar*[_right[iProc*numCats+iCat]];
      sigmaslides[iCat] = new RooFormulaVar*[_right[iProc*numCats+iCat]];
      gslides    [iCat] = new RooGaussian  *[_right[iProc*numCats+iCat]];
      
      RooArgList compListRight;
      RooArgList fracListRight;

      for(int iComp=0; iComp < _right[iProc*numCats+iCat]; ++iComp) {
	stringstream pSS;
	pSS<<"mean"<<(iComp+1)<<"slide"+catnames.at(iCat);
	stringstream pSSmean;
	pSSmean << "dm"<<iComp<<"_"<<procnames[iProc]<<"_rVtx"; // careful !!! this must match above definition!
	std::map<TString,RooHistFunc*>::iterator pair = fitparmfuncs[iCat].find(TString(pSSmean.str().c_str()));
	if( pair == fitparmfuncs[iCat].end() ) {
	  std::cerr<<" ERROR: Could not find RooHistFunction for parameter "<<pSSmean.str().c_str()<<"."<<std::endl;
	  return;
	}	
	meanslides [iCat][iComp] = new RooFormulaVar(pSS.str().c_str(),"","@0 + @1 + @0*@2",RooArgList(mnom,*(pair->second),*nuissancedeltams[iCat]));
	pSS.str("");
	pSSmean.str("");
	pSS<<"sigma"<<(iComp+1)<<"slide"+catnames.at(iCat);	
	pSSmean << "sigma"<<iComp<<"_"<<procnames[iProc]<<"_rVtx"; // careful !!! this must match above definition!
	pair = fitparmfuncs[iCat].find(TString(pSSmean.str().c_str()));
	if( pair == fitparmfuncs[iCat].end() ) {
	  std::cerr<<" ERROR: Could not find RooHistFunction for parameter "<<pSSmean.str().c_str()<<"."<<std::endl;
	  return;
	}	
	sigmaslides[iCat][iComp] = new RooFormulaVar(pSS.str().c_str(),"","TMath::Max(0.01,sqrt(@0*@0-@3*@3*@2*@2 +@1*@1))",RooArgList(*(pair->second),*smearmods[iCat],*smears[iCat],mnom));//why set min 0.01?	
	pSS.str("");
	pSS<<"g"<<(iComp+1)<<"slide" << catnames.at(iCat);
	gslides[iCat][iComp] = new RooGaussian(pSS.str().c_str(),"",*hmass,*(meanslides[iCat][iComp]),*(sigmaslides[iCat][iComp]));
	compListRight.add(* (gslides[iCat][iComp]) );
	
	
	if(iComp > 0) {
	  // collect the fraction parameters
	  pSSmean.str("");
	  pSSmean << "f"<<(iComp-1)<<"_"<<procnames[iProc]<<"_rVtx"; // careful !!! this must match above definition!
	  pair = fitparmfuncs[iCat].find(TString(pSSmean.str().c_str()));
	  if( pair == fitparmfuncs[iCat].end() ) {
	    std::cerr<<" ERROR: Could not find RooHistFunction for parameter "<<pSSmean.str().c_str()<<"."<<std::endl;
	    return;
	  }
	  fracListRight.add( *(pair->second) );
	}
      }

      combhslides[iCat] = new RooAddPdf(TString("combhslide")+catnames.at(iCat),"",compListRight,fracListRight);

      // ------------------------------------------------------------------------------------------------------------------
      // Wrong Vertex Models
      wmeanslides [iCat] = new RooFormulaVar*[_wrong[iProc*numCats+iCat]];
      wsigmaslides[iCat] = new RooFormulaVar*[_wrong[iProc*numCats+iCat]];
      wgslides    [iCat] = new RooGaussian  *[_wrong[iProc*numCats+iCat]];
      
      RooArgList compListWrong;
      RooArgList fracListWrong;

      for(int iComp=0; iComp < _wrong[iProc*numCats+iCat]; ++iComp) {
	stringstream pSS;
	pSS<<"wmean"<<(iComp+1)<<"slide"+catnames.at(iCat);
	stringstream pSSmean;
	pSSmean << "dm"<<iComp<<"_"<<procnames[iProc]<<"_wVtx"; // careful !!! this must match above definition!
	std::map<TString,RooHistFunc*>::iterator pair = fitparmfuncs[iCat].find(TString(pSSmean.str().c_str()));
	if( pair == fitparmfuncs[iCat].end() ) {
	  std::cerr<<" ERROR: Could not find RooHistFunction for parameter "<<pSSmean.str().c_str()<<"."<<std::endl;
	  return;
	}	
	wmeanslides [iCat][iComp] = new RooFormulaVar(pSS.str().c_str(),"","@0 + @1 + @0*@2",RooArgList(mnom,*(pair->second),*nuissancedeltams[iCat]));
	pSS.str("");
	pSSmean.str("");
	pSS<<"wsigma"<<(iComp+1)<<"slide"+catnames.at(iCat);	
	pSSmean << "sigma"<<iComp<<"_"<<procnames[iProc]<<"_wVtx"; // careful !!! this must match above definition!
	pair = fitparmfuncs[iCat].find(TString(pSSmean.str().c_str()));
	if( pair == fitparmfuncs[iCat].end() ) {
	  std::cerr<<" ERROR: Could not find RooHistFunction for parameter "<<pSSmean.str().c_str()<<"."<<std::endl;
	  return;
	}	
	wsigmaslides[iCat][iComp] = new RooFormulaVar(pSS.str().c_str(),"","TMath::Max(0.01,sqrt(@0*@0-@3*@3*@2*@2 +@1*@1))",RooArgList(*(pair->second),*smearmods[iCat],*smears[iCat],mnom));//why set min 0.01?	
	pSS.str("");
	pSS<<"g"<<(iComp+1)<<"slide" << catnames.at(iCat);
	wgslides[iCat][iComp] = new RooGaussian(pSS.str().c_str(),"",*hmass,* (wmeanslides[iCat][iComp]),*(wsigmaslides[iCat][iComp]));
	compListWrong.add(* (wgslides[iCat][iComp]) );
	
	
	if(iComp > 0) {
	  // collect the fraction parameters
	  pSSmean.str("");
	  pSSmean << "f"<<(iComp-1)<<"_"<<procnames[iProc]<<"_wVtx"; // careful !!! this must match above definition!
	  pair = fitparmfuncs[iCat].find(TString(pSSmean.str().c_str()));
	  if( pair == fitparmfuncs[iCat].end() ) {
	    std::cerr<<" ERROR: Could not find RooHistFunction for parameter "<<pSSmean.str().c_str()<<"."<<std::endl;
	    return;
	  }
	  fracListWrong.add( *(pair->second) );
	}
      }

      combhwrongslides[iCat] = new RooAddPdf(TString("combhwrongslide")+catnames.at(iCat),"",compListWrong,fracListWrong);

      // add right & wrong together
      stringstream pSS;
      pSS << "fracright_" << procnames[iProc] << "_" << iCat;
      std::map<TString,RooHistFunc*>::iterator pair = fitparmfuncs[iCat].find(TString(pSS.str().c_str()));
      if( pair == fitparmfuncs[iCat].end() ) {
	std::cerr<<" ERROR: Could not find RooHistFunction for parameter "<<pSS.str().c_str()<<"."<<std::endl;
	return;
      }
      
      fracrightmodslides[iCat] = new RooFormulaVar(TString("fracrightmodslide")+catnames.at(iCat),"","TMath::Min(@0*@1,1.0)",RooArgList(nuissancedeltafracright,* (pair->second) ));      
      combhvtxslides[iCat] = new RooAddPdf(TString("combhvtxslide")+catnames.at(iCat),"",RooArgList(*combhslides[iCat],*combhwrongslides[iCat]),RooArgList(*fracrightmodslides[iCat]),kTRUE);

      pSS.str("");
      pSS << "effacc_" << procnames[iProc] << "_" << iCat;
      pair = fitparmfuncs[iCat].find(TString(pSS.str().c_str()));
      
      if( pair == fitparmfuncs[iCat].end() ) {
	std::cerr<<" ERROR: Could not find RooHistFunction for parameter "<<pSS.str().c_str()<<"."<<std::endl;
	return;
      }
      finalnormslides[iCat] = pair->second;

      //std::cout<" ********* "<<<finalnormslides[iCat]->getVal()<<std::endl;

      // This is smoewhat dummy, but to keep annotation....
      finalpdfslides[iCat]  = combhvtxslides[iCat];
    }
    
    // FIX-ME: Add back interpolation test...

    // open again the config file and read the normalization buisiness...
    // ... Catgeory normalizations
    std::vector<RooAbsReal*> nsigcats;
    std::vector<RooAbsReal*> addNuissance;
    
    std::cout<<" ======================================================== "<<std::endl;

    status = readConfigCardNuissances(configCardName, numCats,
				      nsigcats,
				      addNuissance,
				      finalnormslides,
				      catnames);
    
    std::cout<<" ======================================================== "<<std::endl;

    if(!status) {
      std::cerr<<" ERROR when readin input card "<<configCardName<<"."<<std::endl;
      return;
    }
    
    //    return;

//     // ************************ COPY PASTE ONE-TO-ONE ***********************************
//     //r9 migration effect on effacc in categories
//     RooRealVar nuissancedeltar9fracbarrel("CMS_hgg_nuissancedeltar9fracbarrel","",1.0,0.1,10.0);
//     nuissancedeltar9fracbarrel.setConstant();
//     RooRealVar nuissancedeltar9fracmixed("CMS_hgg_nuissancedeltar9fracmixed","",1.0,0.1,10.0);
//     nuissancedeltar9fracmixed.setConstant();  
    
//     std::vector<RooAbsReal*> nsigcats;
//     RooFormulaVar nsigcat0(TString::Format("nsig%s",catnames.at(0).Data()),"","@0*@1",RooArgList(nuissancedeltar9fracbarrel,*finalnormslides[0]));
//     RooFormulaVar nsigcat1(TString::Format("nsig%s",catnames.at(1).Data()),"","(1.0-@0)*@1 + @2",RooArgList(nuissancedeltar9fracbarrel,*finalnormslides[0],*finalnormslides[1]));
//     RooFormulaVar nsigcat2(TString::Format("nsig%s",catnames.at(2).Data()),"","@0*@1",RooArgList(nuissancedeltar9fracmixed,*finalnormslides[2]));
//     RooFormulaVar nsigcat3(TString::Format("nsig%s",catnames.at(3).Data()),"","(1.0-@0)*@1 + @2",RooArgList(nuissancedeltar9fracmixed,*finalnormslides[2],*finalnormslides[3]));
//     RooFormulaVar nsigcat4(TString::Format("nsig%s",catnames.at(4).Data()),"","@0",RooArgList(*finalnormslides[4]));
//     //   RooFormulaVar nsigcat5(TString::Format("nsig%s",catnames.at(5).Data()),"","@0",RooArgList(*finalnormslides[5]));
//     //   RooFormulaVar nsigcat6(TString::Format("nsig%s",catnames.at(6).Data()),"","@0",RooArgList(*finalnormslides[6]));
    
//     nsigcats.push_back(&nsigcat0);
//     nsigcats.push_back(&nsigcat1);  
//     nsigcats.push_back(&nsigcat2);
//     nsigcats.push_back(&nsigcat3);
//     nsigcats.push_back(&nsigcat4);  
//     //   nsigcats.push_back(&nsigcat5);  
//     //   nsigcats.push_back(&nsigcat6);  
  
    RooArgList addnorm;
    RooArgList addnormff;
    
    RooArgList addpdfs;
    RooArgList addpdfsff;
    
    RooArgList addcoeffs;  
    RooArgList addcoeffsff;  
    
    std::vector<RooAbsReal*> combnorms;
    std::vector<RooAbsReal*> combnormsff;
    
    //final pdfs
    for (UInt_t icat=0; icat<catnames.size(); ++icat) {
                  
      RooAbsPdf *hggpdfsmabs = finalpdfslides[icat];
      hggpdfsmabs->SetName(TString("hggpdfsmabs_")+catnames.at(icat));
      RooAbsPdf *hggpdfsmrel = (RooAbsPdf*)finalpdfslides[icat]->Clone(TString::Format("hggpdfsmrel_%s",catnames.at(icat).Data()));//a function of mnom
      
//       RooAbsPdf *hggpdfffabs = finalpdfslides[icat];
//       hggpdfffabs->SetName(TString("hggpdfffabs_")+catnames.at(icat));
//       RooAbsPdf *hggpdfffrel = (RooAbsPdf*)finalpdfslides[icat]->Clone(TString::Format("hggpdfffrel_%s",catnames.at(icat).Data()));


//       RooAbsPdf *hggpdfsm4abs = (RooAbsPdf*)finalpdfslides[icat]->Clone(TString::Format("hggpdfsm4abs_%s",catnames.at(icat).Data()));    
//       RooAbsPdf *hggpdfsm4rel = (RooAbsPdf*)finalpdfslides[icat]->Clone(TString::Format("hggpdfsm4rel_%s",catnames.at(icat).Data()));



// =====================================================================

      RooAbsReal* tempAbs = addFuncMap.find("smxsec")->second;
      if (!tempAbs) {
	std::cerr<<" ERROR: Could not find addFunc with name smxsec."<<std::endl;
	return;
      }
      RooFormulaVar *nsigsmabs = new RooFormulaVar(TString::Format("hggpdfsmabs_%s_norm",catnames.at(icat).Data()),"","@0*@1/@2",RooArgList(*(procxseclist[iProc]),*nsigcats[icat],*tempAbs));
      
      tempAbs = histFuncMap.find("smbr")->second;
      if (!tempAbs) {
	std::cerr<<" ERROR: Could not find histFunc with name smbr."<<std::endl;
	return;
      }
      
      RooFormulaVar *nsigsmrel = new RooFormulaVar(TString::Format("hggpdfsmrel_%s_norm",catnames.at(icat).Data()),"","@0*@1*@2",RooArgList(*tempAbs,*procxseclist[iProc],*nsigcats[icat]));
      

      for(std::map<TString,RooAddition*>::iterator it = addFuncMap.begin(); it != addFuncMap.end(); ++it) {
	std::cout<<it->first<<std::endl;
      }

// =====================================================================

//       RooAbsReal* tempAbs = addFuncMap.find("ffxsec")->second;

//       if (!tempAbs) {
// 	std::cerr<<" ERROR: Could not find addFunc with name ffxsec."<<std::endl;
// 	return;
//       }

      
//       RooFormulaVar *nsigffabs = new RooFormulaVar(TString::Format("hggpdfffabs_%s_norm",catnames.at(icat).Data()),"","@0*@1/@2",RooArgList(*(procxseclist[iProc]),*nsigcats[icat],*tempAbs));
      
//       std::cout<<nsigffabs->getVal()<<std::endl;
      
//       tempAbs = histFuncMap.find("ffbr")->second;
//       if (!tempAbs) {
// 	std::cerr<<" ERROR: Could not find histFunc with name ffbr."<<std::endl;
// 	return;
//       }
      
//        RooFormulaVar *nsigffrel = new RooFormulaVar(TString::Format("hggpdfffrel_%s_norm",catnames.at(icat).Data()),"","@0*@1*@2",RooArgList(*tempAbs,*procxseclist[iProc],*nsigcats[icat]));


// =====================================================================
      
//       RooFormulaVar *nsigffabs = new RooFormulaVar(TString::Format("hggpdfffabs_%s_norm",catnames.at(icat).Data()),"","@0*@1/@2",RooArgList(*(procxseclist[iProc]),*nsigcats[icat],*fffxsecs));
//       RooFormulaVar *nsigffrel = new RooFormulaVar(TString::Format("hggpdfffrel_%s_norm",catnames.at(icat).Data()),"","@0*@1*@2",RooArgList(*fffbrs,*procxsecs,*nsigcats[icat]));    

//       RooFormulaVar *nsigsm4abs = new RooFormulaVar(TString::Format("hggpdfsm4abs_%s_norm",catnames.at(icat).Data()),"","@0*@1/@2",RooArgList(*(procxseclist[iProc]),*nsigcats[icat],*fsm4xsecs));
//       RooFormulaVar *nsigsm4rel = new RooFormulaVar(TString::Format("hggpdfsm4rel_%s_norm",catnames.at(icat).Data()),"","@0*@1*@2",RooArgList(*fsm4brs,*procxsecssm4,*nsigcats[icat]));


      RooExtendPdf *sigpdfsmabs = new RooExtendPdf(TString::Format("sigpdfsmabs%s",catnames.at(icat).Data()),"",*hggpdfsmabs,*nsigsmabs);
      RooExtendPdf *sigpdfsmrel = new RooExtendPdf(TString::Format("sigpdfsmrel%s",catnames.at(icat).Data()),"",*hggpdfsmrel,*nsigsmrel);
      //RooExtendPdf *sigpdfffabs = new RooExtendPdf(TString::Format("sigpdfffabs%s",catnames.at(icat).Data()),"",*hggpdfffabs,*nsigffabs);
      //RooExtendPdf *sigpdfffrel = new RooExtendPdf(TString::Format("sigpdfffrel%s",catnames.at(icat).Data()),"",*hggpdfffrel,*nsigffrel);

//       RooExtendPdf *sigpdfsm4abs = new RooExtendPdf(TString::Format("sigpdfsm4abs%s",catnames.at(icat).Data()),"",*hggpdfsm4abs,*nsigsm4abs);
//       RooExtendPdf *sigpdfsm4rel = new RooExtendPdf(TString::Format("sigpdfsm4rel%s",catnames.at(icat).Data()),"",*hggpdfsm4rel,*nsigsm4rel);

      if(!fitonly) {
	wOut->import(*sigpdfsmabs,RecycleConflictNodes());
	wOut->import(*sigpdfsmrel,RecycleConflictNodes());
	//wOut->import(*sigpdfffabs,RecycleConflictNodes());
	//wOut->import(*sigpdfffrel,RecycleConflictNodes());    
      }

//       w->import(*sigpdfsm4abs,RecycleConflictNodes());
//       w->import(*sigpdfsm4rel,RecycleConflictNodes());    

      addnorm.add(*nsigsmrel);
      combnorms.push_back(nsigsmrel);
      
      //addnormff.add(*nsigffrel);
      //combnormsff.push_back(nsigffrel);      
    }
    
    //save everything to file with RooWorkspace
  
    //wOut->Print();
    if(!fitonly)
      wOut->writeToFile("ubersignalmodel.root") ;
  
    gSystem->cd("../../");

  } // go to next process...

  return;

}

// -----------------------------------------------------------------------------------------------------------------
RooAbsPdf* generateMultiGaussian(RooRealVar* mass, 
				 RooRealVar* nomMass,
				 int numComp,
				 RooGaussian**& g,
				 TString procName, TString quali,
				 std::map<TString,RooRealVar*>& fitparms,
				 std::map<TString,float>&       startvals,
				 float* startMean, float* startSigma, float* startFrac, bool fix) {

  // generate one Gaussian per component
  RooRealVar** dm    = new RooRealVar*[numComp];
  RooRealVar** sigma = new RooRealVar*[numComp];
  RooRealVar** f = NULL;

  RooFormulaVar** mean = new RooFormulaVar*[numComp];

  //RooGaussian** g = new RooGaussian*[numComp];
  g = new RooGaussian*[numComp];

  RooAddPdf* comb = NULL;
  
  RooAbsPdf* returnPdf = NULL;

  if(numComp > 1) f = new RooRealVar*[numComp-1];

  stringstream pSS;
  for(int iComp=0; iComp < numComp; ++iComp) {
    pSS.str("");
    pSS << "dm"<<iComp<<"_"<<procName<<"_"<<quali;
    dm[iComp] = new RooRealVar(pSS.str().c_str(),"",-1.,1.);
    dm[iComp] -> removeRange();
    if( fix ) dm[iComp]->setConstant();
    fitparms.insert(std::pair<TString, RooRealVar*>(pSS.str().c_str(),dm[iComp]));
    startvals.insert(std::pair<TString,float>(pSS.str().c_str(),startMean[iComp]));
    pSS.str("");
    pSS << "mean"<<iComp<<"_"<<procName<<"_"<<quali;
    mean[iComp] = new RooFormulaVar(pSS.str().c_str(),"","@0+@1",RooArgList(*nomMass,*dm[iComp]));
    pSS.str("");
    pSS << "sigma"<<iComp<<"_"<<procName<<"_"<<quali;
    sigma[iComp] = new RooRealVar(pSS.str().c_str(),"",-1.,1.);
    sigma[iComp] -> removeRange();
    //if( fix ) sigma[iComp]->setConstant();
    fitparms.insert(std::pair<TString,RooRealVar*>(pSS.str(),sigma[iComp]));
    startvals.insert(std::pair<TString,float>(pSS.str().c_str(),startSigma[iComp]));

    pSS.str("");
    pSS << "g"<<iComp<<"_"<<procName<<"_"<<quali;
    g[iComp] = new RooGaussian(pSS.str().c_str(),"",*mass,*mean[iComp],*sigma[iComp]);

    if(iComp > 0) {
      pSS.str("");
      pSS << "f"<<iComp-1<<"_"<<procName<<"_"<<quali;
      f[iComp-1] = new RooRealVar(pSS.str().c_str(),"",0.,1.);
      fitparms.insert(std::pair<TString,RooRealVar*>(pSS.str(),f[iComp-1]));
      startvals.insert(std::pair<TString,float>(pSS.str().c_str(),startFrac[iComp-1]));
      //f[iComp-1] -> removeRange();
    }      
  }
   
  RooArgList gausList;
  RooArgList fracList;
  if( numComp > 1) {
    for(int iComp = 0; iComp < numComp; ++iComp) {
      gausList.add( (*g[iComp]) );
      if(iComp > 0) fracList.add( (*f[iComp-1]) );
    }
    pSS.str("");
    pSS << "combh_"<<procName<<"_"<<quali;
    comb = new RooAddPdf(pSS.str().c_str(),"",gausList,fracList,true);
    returnPdf = (RooAbsPdf*) comb;
  } else
    returnPdf = (RooAbsPdf*) g[0];
  
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
		    std::vector<TString>& procNames, std::vector<bool>& procOn,
		    TString& inputDir, TString& outputDir,
		    TString& wsPrefix,
		    int& numProcs, int& numCats, int& numModels, int& numMasses,
		    bool*& catIsOn,
		    double& massmax, double& massmin,
		    std::vector<double>& mhs,
		    std::vector<double>& smearingv,
		    std::vector<TString>& catNames,
		    int*& _right, int*& _wrong,
		    float**& _meanStartRight, float**& _sigmaStartRight, float**& _fracStartRight,
		    float**& _meanStartWrong, float**& _sigmaStartWrong, float**& _fracStartWrong) {

  
  FILE* configFile = fopen(configCardName.Data(),"r");
  if ( !configFile ) {
    std::cerr<<" Inputfile "<<configCardName<<" not found."<<std::endl;
    return false;
  }
  
  mhs.resize(0);
  
  char line[100];
  
  std::cout<<" Reading model from file "<<configCardName<<"..."<<std::endl;
  
  int   iProc = -1;
  int   iCat  = -1;
  int   massIdx = -1;
  int   nRight = -1;
  int   nWrong = -1;
  
  int   whichCat = -1 ;
  float smearing = -1.;
  float massval  = -1.;
  
  char rightStart[100];
  char wrongStart[100];
  
  char dummyString[100];
  
  int dummy1, dummy2, dummy3;
  //# ------------------------------------------------------------------------------------------------------------------------------------------------
  //# Section on the Porcesses/Events
  //#	num proc	num photon cats		num cats(auxiliary)	num cats (analysis)	num cats (trigger)	numModels	numMasses
  //# ------------------------------------------------------------------------------------------------------------------------------------------------
  //INIT	4		4			7			4			16			1		5
  while (fgets(line,100,configFile)) {
    nRight = -1;
    nWrong = -1;
    if ( !sscanf(&line[0], "#") ) {   // not a document line      
      if ( sscanf(line, "INIT %d %d %d %d %d %d %d", &numProcs, &dummy1, &dummy2, &numCats, &dummy3, &numModels, &numMasses) ) { // this is the INIT line	
	if( numProcs < 1 || numCats < 1 ) {
	  std::cerr<<" Error in config file "<<configCardName<<" : Number of processes and categories must be positive."<<std::endl;
	  return false;
	}
	if( _right || _wrong ) {
	  std::cerr<<" Error in config file "<<configCardName<<" : Double initialization."<<std::endl;
	  return false;
	}
	_right = new int[numProcs*numCats];
	_wrong = new int[numProcs*numCats];

	catIsOn = new bool[numProcs*numCats];

	smearingv.resize(numCats);
	catNames.resize(numCats);

	_meanStartRight  = new float*[numProcs*numCats];
	_sigmaStartRight = new float*[numProcs*numCats];
	_fracStartRight  = new float*[numProcs*numCats];

	_meanStartWrong  = new float*[numProcs*numCats];
	_sigmaStartWrong = new float*[numProcs*numCats];
	_fracStartWrong  = new float*[numProcs*numCats];

	procNames.resize(numProcs);
	procOn.resize(numProcs);
	// # ------------------------------------------------------------------------------------------------------------------------------------------------
	// # analysis categpry definitions
	// #	idx	name		eff. smearing	BG-model/Order	TCuts *use only AUXCATs from above*		StDesc(String descriptor) (*for plots and such*)
	// # ------------------------------------------------------------------------------------------------------------------------------------------------
	// ANACAT	0	cat0		0.005432	Bern/5		" basecut && !vbfcut && baseline0 "		StDesc(Baseline Cat 1)
	// ANACAT	1	cat1		0.005196	Bern/5		" basecut && !vbfcut && baseline1 "		StDesc(Baseline Cat 2)
	// ANACAT	2	cat2		0.007464	Bern/5		" basecut && !vbfcut && baseline2 "		StDesc(Baseline Cat 3)
	// ANACAT	3	cat3		0.012978	Bern/5		" basecut && !vbfcut && baseline3 "		StDesc(Baseline Cat 4)
	// ANACAT	4	cat4		0.008722	Bern/3		" masscut && vbfcut "				StDesc(Baseline VBF Cat)
      } else if (  sscanf(line, "ANACAT %d %s %f", &whichCat, &dummyString, &smearing) ) { // this is the SMEARING line
	if(whichCat < 0 || whichCat >= numCats) {
	  std::cerr<<" Error in config file "<<configCardName<<" : Setting smearing for Cat "<<whichCat<<" with only "<<numCats<<" total Cats."<<std::endl;
	  return false;
	}
	smearingv[whichCat] = (double) smearing;
	catNames[whichCat]  = TString(dummyString);
      } else if  (  sscanf(line, "OFF %d %d %d ( %s ) %d ( %s )", &iProc, &iCat, &nRight, &rightStart, &nWrong, &wrongStart ) ) {
	if( iProc < 0 || iCat < 0 || nRight < 1 || nWrong < 1 ) break;
	if (!validateInput(iProc,iCat,numProcs,numCats)) break;
	_right[iProc*numCats+iCat] = nRight;
	_wrong[iProc*numCats+iCat] = nWrong;

	catIsOn[iProc*numCats+iCat] = false;

	// read the parameters
	if( nRight > 3 || nWrong > 3 ) {
	  std::cerr<<" To many Gaussians: "<<nRight<<"  "<<nWrong<<std::endl;
	  return false;
	}

	switch(nRight) {
	case 1:
	  _meanStartRight[iProc*numCats+iCat]  = new float[1];
	  _sigmaStartRight[iProc*numCats+iCat] = new float[1];
	  sscanf(rightStart, "m:%fs:%f", &(_meanStartRight[iProc*numCats+iCat][0]), &(_sigmaStartRight[iProc*numCats+iCat][0]));
	  break;
	case 2:
	  _meanStartRight[iProc*numCats+iCat]  = new float[2];
	  _sigmaStartRight[iProc*numCats+iCat] = new float[2];
	  _fracStartRight[iProc*numCats+iCat]  = new float[1];
	  sscanf(rightStart, "m:%fm:%fs:%fs:%ff:%f", &(_meanStartRight[iProc*numCats+iCat][0]), &(_meanStartRight[iProc*numCats+iCat][1]), &(_sigmaStartRight[iProc*numCats+iCat][0]), &(_sigmaStartRight[iProc*numCats+iCat][1]), &(_fracStartRight[iProc*numCats+iCat][0]));
	  break;
	case 3:
	  _meanStartRight[iProc*numCats+iCat]  = new float[3];
	  _sigmaStartRight[iProc*numCats+iCat] = new float[3];
	  _fracStartRight[iProc*numCats+iCat]  = new float[2];
	  sscanf(rightStart, "m:%fm:%fm:%fs:%fs:%fs:%ff:%ff:%f", &(_meanStartRight[iProc*numCats+iCat][0]), &(_meanStartRight[iProc*numCats+iCat][1]), &(_meanStartRight[iProc*numCats+iCat][2]), &(_sigmaStartRight[iProc*numCats+iCat][0]), &(_sigmaStartRight[iProc*numCats+iCat][1]), &(_sigmaStartRight[iProc*numCats+iCat][2]), &(_fracStartRight[iProc*numCats+iCat][0]), &(_fracStartRight[iProc*numCats+iCat][1]));
	  break;
	}

	switch(nWrong) {
	case 1:
	  _meanStartWrong[iProc*numCats+iCat]  = new float[1];
	  _sigmaStartWrong[iProc*numCats+iCat] = new float[1];
	  sscanf(wrongStart, "m:%fs:%f", &(_meanStartWrong[iProc*numCats+iCat][0]), &(_sigmaStartWrong[iProc*numCats+iCat][0]));
	  break;
	case 2:
	  _meanStartWrong[iProc*numCats+iCat]  = new float[2];
	  _sigmaStartWrong[iProc*numCats+iCat] = new float[2];
	  _fracStartWrong[iProc*numCats+iCat]  = new float[1];
	  sscanf(wrongStart, "m:%fm:%fs:%fs:%ff:%f", &(_meanStartWrong[iProc*numCats+iCat][0]), &(_meanStartWrong[iProc*numCats+iCat][1]), &(_sigmaStartWrong[iProc*numCats+iCat][0]), &(_sigmaStartWrong[iProc*numCats+iCat][1]), &(_fracStartWrong[iProc*numCats+iCat][0]));
	  break;
	case 3:
	  _meanStartWrong[iProc*numCats+iCat]  = new float[3];
	  _sigmaStartWrong[iProc*numCats+iCat] = new float[3];
	  _fracStartWrong[iProc*numCats+iCat]  = new float[2];
	  sscanf(wrongStart, "m:%fm:%fm:%fs:%fs:%fs:%ff:%ff:%f", &(_meanStartWrong[iProc*numCats+iCat][0]), &(_meanStartWrong[iProc*numCats+iCat][1]), &(_meanStartWrong[iProc*numCats+iCat][2]), &(_sigmaStartWrong[iProc*numCats+iCat][0]), &(_sigmaStartWrong[iProc*numCats+iCat][1]), &(_sigmaStartWrong[iProc*numCats+iCat][2]), &(_fracStartWrong[iProc*numCats+iCat][0]), &(_fracStartWrong[iProc*numCats+iCat][1]));
	  break;
	}

      } else if (  sscanf(line, "ON %d %d %d ( %s ) %d ( %s )", &iProc, &iCat, &nRight, &rightStart, &nWrong, &wrongStart ) ) {
	if( iProc < 0 || iCat < 0 || nRight < 1 || nWrong < 1 ) break;
	if (!validateInput(iProc,iCat,numProcs,numCats)) break;
	_right[iProc*numCats+iCat] = nRight;
	_wrong[iProc*numCats+iCat] = nWrong;

	catIsOn[iProc*numCats+iCat] = true;

	// read the parameters
	if( nRight > 3 || nWrong > 3 ) {
	  std::cerr<<" To many Gaussians: "<<nRight<<"  "<<nWrong<<std::endl;
	  return false;
	}

	switch(nRight) {
	case 1:
	  _meanStartRight[iProc*numCats+iCat]  = new float[1];
	  _sigmaStartRight[iProc*numCats+iCat] = new float[1];
	  sscanf(rightStart, "m:%fs:%f", &(_meanStartRight[iProc*numCats+iCat][0]), &(_sigmaStartRight[iProc*numCats+iCat][0]));
	  break;
	case 2:
	  _meanStartRight[iProc*numCats+iCat]  = new float[2];
	  _sigmaStartRight[iProc*numCats+iCat] = new float[2];
	  _fracStartRight[iProc*numCats+iCat]  = new float[1];
	  sscanf(rightStart, "m:%fm:%fs:%fs:%ff:%f", &(_meanStartRight[iProc*numCats+iCat][0]), &(_meanStartRight[iProc*numCats+iCat][1]), &(_sigmaStartRight[iProc*numCats+iCat][0]), &(_sigmaStartRight[iProc*numCats+iCat][1]), &(_fracStartRight[iProc*numCats+iCat][0]));
	  break;
	case 3:
	  _meanStartRight[iProc*numCats+iCat]  = new float[3];
	  _sigmaStartRight[iProc*numCats+iCat] = new float[3];
	  _fracStartRight[iProc*numCats+iCat]  = new float[2];
	  sscanf(rightStart, "m:%fm:%fm:%fs:%fs:%fs:%ff:%ff:%f", &(_meanStartRight[iProc*numCats+iCat][0]), &(_meanStartRight[iProc*numCats+iCat][1]), &(_meanStartRight[iProc*numCats+iCat][2]), &(_sigmaStartRight[iProc*numCats+iCat][0]), &(_sigmaStartRight[iProc*numCats+iCat][1]), &(_sigmaStartRight[iProc*numCats+iCat][2]), &(_fracStartRight[iProc*numCats+iCat][0]), &(_fracStartRight[iProc*numCats+iCat][1]));
	  break;
	}

	switch(nWrong) {
	case 1:
	  _meanStartWrong[iProc*numCats+iCat]  = new float[1];
	  _sigmaStartWrong[iProc*numCats+iCat] = new float[1];
	  sscanf(wrongStart, "m:%fs:%f", &(_meanStartWrong[iProc*numCats+iCat][0]), &(_sigmaStartWrong[iProc*numCats+iCat][0]));
	  break;
	case 2:
	  _meanStartWrong[iProc*numCats+iCat]  = new float[2];
	  _sigmaStartWrong[iProc*numCats+iCat] = new float[2];
	  _fracStartWrong[iProc*numCats+iCat]  = new float[1];
	  sscanf(wrongStart, "m:%fm:%fs:%fs:%ff:%f", &(_meanStartWrong[iProc*numCats+iCat][0]), &(_meanStartWrong[iProc*numCats+iCat][1]), &(_sigmaStartWrong[iProc*numCats+iCat][0]), &(_sigmaStartWrong[iProc*numCats+iCat][1]), &(_fracStartWrong[iProc*numCats+iCat][0]));
	  break;
	case 3:
	  _meanStartWrong[iProc*numCats+iCat]  = new float[3];
	  _sigmaStartWrong[iProc*numCats+iCat] = new float[3];
	  _fracStartWrong[iProc*numCats+iCat]  = new float[2];
	  sscanf(wrongStart, "m:%fm:%fm:%fs:%fs:%fs:%ff:%ff:%f", &(_meanStartWrong[iProc*numCats+iCat][0]), &(_meanStartWrong[iProc*numCats+iCat][1]), &(_meanStartWrong[iProc*numCats+iCat][2]), &(_sigmaStartWrong[iProc*numCats+iCat][0]), &(_sigmaStartWrong[iProc*numCats+iCat][1]), &(_sigmaStartWrong[iProc*numCats+iCat][2]), &(_fracStartWrong[iProc*numCats+iCat][0]), &(_fracStartWrong[iProc*numCats+iCat][1]));
	  break;
	}


      } else if (  sscanf(line, "OFF ~ %d %d ( %s ) %d ( %s )", &iCat, &nRight, &rightStart, &nWrong, &wrongStart) ) {
	if( iProc < 0 || iCat < 0 || nRight < 1 || nWrong < 1 ) break;
	if (!validateInput(iProc,iCat,numProcs,numCats)) break;
	_right[iProc*numCats+iCat] = nRight;
	_wrong[iProc*numCats+iCat] = nWrong;

	catIsOn[iProc*numCats+iCat] = false;

	// read the parameters
	if( nRight > 3 || nWrong > 3 ) {
	  std::cerr<<" To many Gaussians: "<<nRight<<"  "<<nWrong<<std::endl;
	  return false;
	}

	switch(nRight) {
	case 1:
	  _meanStartRight[iProc*numCats+iCat]  = new float[1];
	  _sigmaStartRight[iProc*numCats+iCat] = new float[1];
	  sscanf(rightStart, "m:%fs:%f", &(_meanStartRight[iProc*numCats+iCat][0]), &(_sigmaStartRight[iProc*numCats+iCat][0]));
	  break;
	case 2:
	  _meanStartRight[iProc*numCats+iCat]  = new float[2];
	  _sigmaStartRight[iProc*numCats+iCat] = new float[2];
	  _fracStartRight[iProc*numCats+iCat]  = new float[1];
	  sscanf(rightStart, "m:%fm:%fs:%fs:%ff:%f", &(_meanStartRight[iProc*numCats+iCat][0]), &(_meanStartRight[iProc*numCats+iCat][1]), &(_sigmaStartRight[iProc*numCats+iCat][0]), &(_sigmaStartRight[iProc*numCats+iCat][1]), &(_fracStartRight[iProc*numCats+iCat][0]));
	  break;
	case 3:
	  _meanStartRight[iProc*numCats+iCat]  = new float[3];
	  _sigmaStartRight[iProc*numCats+iCat] = new float[3];
	  _fracStartRight[iProc*numCats+iCat]  = new float[2];
	  sscanf(rightStart, "m:%fm:%fm:%fs:%fs:%fs:%ff:%ff:%f", &(_meanStartRight[iProc*numCats+iCat][0]), &(_meanStartRight[iProc*numCats+iCat][1]), &(_meanStartRight[iProc*numCats+iCat][2]), &(_sigmaStartRight[iProc*numCats+iCat][0]), &(_sigmaStartRight[iProc*numCats+iCat][1]), &(_sigmaStartRight[iProc*numCats+iCat][2]), &(_fracStartRight[iProc*numCats+iCat][0]), &(_fracStartRight[iProc*numCats+iCat][1]));
	  break;
	}

	switch(nWrong) {
	case 1:
	  _meanStartWrong[iProc*numCats+iCat]  = new float[1];
	  _sigmaStartWrong[iProc*numCats+iCat] = new float[1];
	  sscanf(wrongStart, "m:%fs:%f", &(_meanStartWrong[iProc*numCats+iCat][0]), &(_sigmaStartWrong[iProc*numCats+iCat][0]));
	  break;
	case 2:
	  _meanStartWrong[iProc*numCats+iCat]  = new float[2];
	  _sigmaStartWrong[iProc*numCats+iCat] = new float[2];
	  _fracStartWrong[iProc*numCats+iCat]  = new float[1];
	  sscanf(wrongStart, "m:%fm:%fs:%fs:%ff:%f", &(_meanStartWrong[iProc*numCats+iCat][0]), &(_meanStartWrong[iProc*numCats+iCat][1]), &(_sigmaStartWrong[iProc*numCats+iCat][0]), &(_sigmaStartWrong[iProc*numCats+iCat][1]), &(_fracStartWrong[iProc*numCats+iCat][0]));
	  break;
	case 3:
	  _meanStartWrong[iProc*numCats+iCat]  = new float[3];
	  _sigmaStartWrong[iProc*numCats+iCat] = new float[3];
	  _fracStartWrong[iProc*numCats+iCat]  = new float[2];
	  sscanf(wrongStart, "m:%fm:%fm:%fs:%fs:%fs:%ff:%ff:%f", &(_meanStartWrong[iProc*numCats+iCat][0]), &(_meanStartWrong[iProc*numCats+iCat][1]), &(_meanStartWrong[iProc*numCats+iCat][2]), &(_sigmaStartWrong[iProc*numCats+iCat][0]), &(_sigmaStartWrong[iProc*numCats+iCat][1]), &(_sigmaStartWrong[iProc*numCats+iCat][2]), &(_fracStartWrong[iProc*numCats+iCat][0]), &(_fracStartWrong[iProc*numCats+iCat][1]));
	  break;
	}
      } else if (  sscanf(line, "ON ~ %d %d ( %s ) %d ( %s )", &iCat, &nRight, &rightStart, &nWrong, &wrongStart) ) {
	if( iProc < 0 || iCat < 0 || nRight < 1 || nWrong < 1 ) break;
	if (!validateInput(iProc,iCat,numProcs,numCats)) break;
	_right[iProc*numCats+iCat] = nRight;
	_wrong[iProc*numCats+iCat] = nWrong;

	catIsOn[iProc*numCats+iCat] = true;

	// read the parameters
	if( nRight > 3 || nWrong > 3 ) {
	  std::cerr<<" To many Gaussians: "<<nRight<<"  "<<nWrong<<std::endl;
	  return false;
	}

	switch(nRight) {
	case 1:
	  _meanStartRight[iProc*numCats+iCat]  = new float[1];
	  _sigmaStartRight[iProc*numCats+iCat] = new float[1];
	  sscanf(rightStart, "m:%fs:%f", &(_meanStartRight[iProc*numCats+iCat][0]), &(_sigmaStartRight[iProc*numCats+iCat][0]));
	  break;
	case 2:
	  _meanStartRight[iProc*numCats+iCat]  = new float[2];
	  _sigmaStartRight[iProc*numCats+iCat] = new float[2];
	  _fracStartRight[iProc*numCats+iCat]  = new float[1];
	  sscanf(rightStart, "m:%fm:%fs:%fs:%ff:%f", &(_meanStartRight[iProc*numCats+iCat][0]), &(_meanStartRight[iProc*numCats+iCat][1]), &(_sigmaStartRight[iProc*numCats+iCat][0]), &(_sigmaStartRight[iProc*numCats+iCat][1]), &(_fracStartRight[iProc*numCats+iCat][0]));
	  break;
	case 3:
	  _meanStartRight[iProc*numCats+iCat]  = new float[3];
	  _sigmaStartRight[iProc*numCats+iCat] = new float[3];
	  _fracStartRight[iProc*numCats+iCat]  = new float[2];
	  sscanf(rightStart, "m:%fm:%fm:%fs:%fs:%fs:%ff:%ff:%f", &(_meanStartRight[iProc*numCats+iCat][0]), &(_meanStartRight[iProc*numCats+iCat][1]), &(_meanStartRight[iProc*numCats+iCat][2]), &(_sigmaStartRight[iProc*numCats+iCat][0]), &(_sigmaStartRight[iProc*numCats+iCat][1]), &(_sigmaStartRight[iProc*numCats+iCat][2]), &(_fracStartRight[iProc*numCats+iCat][0]), &(_fracStartRight[iProc*numCats+iCat][1]));
	  break;
	}

	switch(nWrong) {
	case 1:
	  _meanStartWrong[iProc*numCats+iCat]  = new float[1];
	  _sigmaStartWrong[iProc*numCats+iCat] = new float[1];
	  sscanf(wrongStart, "m:%fs:%f", &(_meanStartWrong[iProc*numCats+iCat][0]), &(_sigmaStartWrong[iProc*numCats+iCat][0]));
	  break;
	case 2:
	  _meanStartWrong[iProc*numCats+iCat]  = new float[2];
	  _sigmaStartWrong[iProc*numCats+iCat] = new float[2];
	  _fracStartWrong[iProc*numCats+iCat]  = new float[1];
	  sscanf(wrongStart, "m:%fm:%fs:%fs:%ff:%f", &(_meanStartWrong[iProc*numCats+iCat][0]), &(_meanStartWrong[iProc*numCats+iCat][1]), &(_sigmaStartWrong[iProc*numCats+iCat][0]), &(_sigmaStartWrong[iProc*numCats+iCat][1]), &(_fracStartWrong[iProc*numCats+iCat][0]));
	  break;
	case 3:
	  _meanStartWrong[iProc*numCats+iCat]  = new float[3];
	  _sigmaStartWrong[iProc*numCats+iCat] = new float[3];
	  _fracStartWrong[iProc*numCats+iCat]  = new float[2];
	  sscanf(wrongStart, "m:%fm:%fm:%fs:%fs:%fs:%ff:%ff:%f", &(_meanStartWrong[iProc*numCats+iCat][0]), &(_meanStartWrong[iProc*numCats+iCat][1]), &(_meanStartWrong[iProc*numCats+iCat][2]), &(_sigmaStartWrong[iProc*numCats+iCat][0]), &(_sigmaStartWrong[iProc*numCats+iCat][1]), &(_sigmaStartWrong[iProc*numCats+iCat][2]), &(_fracStartWrong[iProc*numCats+iCat][0]), &(_fracStartWrong[iProc*numCats+iCat][1]));
	  break;
	}
      } else {
	
	// additional setup
	int theProc = -1;
	char name[30];
	char onoff[3];
	char directory[100];
	if( sscanf(line,"PROC %d %s %s",&theProc,&name,&onoff) ) {
	  procNames[theProc]=TString(name).Strip();
	  procOn[theProc]=std::strcmp(onoff,"OFF");
	} else if  ( sscanf(line,"WS1PREFIX %s",&name) ) wsPrefix = TString(name);
	else if  ( sscanf(line,"PROJECTDIR %s",&directory) ) {
	  inputDir=TString(directory).Strip();
	  outputDir=TString(directory).Strip();	  
	}
	else if( sscanf(line,"MINMSS %f",&massval) ) massmin = (double) massval;
	else if( sscanf(line,"MAXMSS %f",&massval) ) massmax = (double) massval;
	else if ( sscanf(line,"MASS %d %f",&massIdx, &massval) ) {
	  if( massIdx > mhs.size() ) {
	    std::cerr<<" ERROR: Ordering in MASS section wrong. Must be idx continuous."<<std::endl;
	    return false;
	  }
	  mhs.push_back( (double) massval );
	}
	//} else {
	//std::cerr<<" Input line "<<line<<" not recopgnized."<<std::endl;
	//return false;
	//}
      }
    }
  }
  
  std::cout<<" done."<<std::endl;
  
  fclose(configFile);
  if(!_right || !_wrong) return false;
  return true;
  
}

bool readConfigCardNuissances(TString configCardName, int numCats,
			      std::vector<RooAbsReal*>& nsigcat,
			      std::vector<RooAbsReal*>& nuissances,
			      std::vector<RooAbsReal*>& finalnorm,
			      std::vector<TString> catnames) {
  
  nsigcat.resize(numCats);
  nuissances.resize(0);
  
  TString cardName = TString("../../")+configCardName;  
  FILE* configFile = fopen(cardName.Data(),"r");
  if ( !configFile ) {
    std::cerr<<" Inputfile "<<configCardName<<" not found."<<std::endl;
    return false;
  }
  
  char line[100];
  
  std::cout<<" Reading model from file "<<configCardName<<"..."<<std::endl;
  
  while (fgets(line,100,configFile)) {
    // must be some nuissance or normalization stuff...
    int indNuis = -1;
    float startVal = -1.;
    float minVal = -1.;
    float maxVal = -1.;
    char name[30];
    if( sscanf(line,"NUIS %d %s %f %f %f",&indNuis, &name, &startVal, &minVal, &maxVal) ) {
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
      if( sscanf(line,"NSIG %d %d %s %s",&catInd,&nParms,&formula, &parList) ) {
	// need to transform the parList...
	RooArgList theList;
	if( !strcmp(formula,"DEFAULT") ) {	 
	  theList.add( *(finalnorm[catInd]) );
	  nsigcat[catInd] = new RooFormulaVar(TString::Format("nsig%s",catnames.at(catInd).Data()),"","@0",theList);
	} else {	
	  std::string parStr(parList);
	  size_t lastPos = 0;
	  while ( parStr.find_last_of(",") != std::string::npos ) {
	    size_t pos = parStr.find_first_of(",");
	    std::string thisPar = parStr.substr(lastPos,pos);
	    int idx = -1;
	    if ( sscanf(thisPar.c_str(),"NUIS(%d)",&idx) )
	      theList.add(* (nuissances[idx]) );
	    else if ( sscanf(thisPar.c_str(),"NOMINAL(%d)",&idx) ) 
	      theList.add(* (finalnorm[idx]) );
	    else {
	      std::cerr<<" ERROR: Cannot recognize Par Type "<<thisPar.c_str()<<"."<<std::endl;
	      return false;
	    }
	    parStr.replace(pos,1,"X");
	    lastPos = pos+1;
	  }
	  // last parameter
	  std::string lastPar = parStr.substr(lastPos);
	  int idx = -1;
	  if ( sscanf(lastPar.c_str(),"NUIS(%d)",&idx) )
	    theList.add(* (nuissances[idx]) );
	  else if ( sscanf(lastPar.c_str(),"NOMINAL(%d)",&idx) )
	    theList.add(* (finalnorm[idx]) );
	  else {
	    std::cerr<<" ERROR: Cannot recognize Par Type "<<lastPar.c_str()<<"."<<std::endl;
	    return false;
	  }
	  RooFormulaVar* tempForm = new RooFormulaVar(TString::Format("nsig%s",catnames.at(catInd).Data()),"",formula,theList);
	  nsigcat[catInd] = tempForm;
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


//---function to change the weights of the dataset---
//we rescale the weight by applying an overall scale factor so that the weights on average are roughly one; we do this because miniut doesn't like event weights far from 1  
RooDataSet *cwdset(RooDataSet *indata, RooRealVar *mvar, RooRealVar *wvar, RooRealVar *pidvar, TString name, Double_t weightscale, int filterproc) {
  RooDataSet *outdata = new RooDataSet(name,"",RooArgList(*mvar,*wvar),wvar->GetName());
  for (Int_t ient=0; ient<indata->numEntries(); ++ient) {
    const RooArgSet *ent = indata->get(ient);
    mvar->setVal(static_cast<RooAbsReal*>(ent->find(mvar->GetName()))->getVal());
    outdata->add(*mvar,weightscale*indata->weight());
  }
  return outdata;
}

void *appendcwd(RooDataSet *outdata, RooDataSet *indata, RooRealVar *mvar, RooRealVar *wvar, Double_t weightscale) {
  //RooDataSet *outdata = new RooDataSet(name,"",RooArgList(*mvar,*wvar),wvar->GetName());
  for (Int_t ient=0; ient<indata->numEntries(); ++ient) {
    const RooArgSet *ent = indata->get(ient);
    mvar->setVal(static_cast<RooAbsReal*>(ent->find(mvar->GetName()))->getVal());
    outdata->add(*mvar,weightscale*indata->weight());
  }
  
  return outdata;
  
}
