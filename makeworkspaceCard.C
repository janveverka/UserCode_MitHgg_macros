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

#include "modelConstants_8TeV.h"

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
    
  return g_reader->EvaluateMVA("Ystar");
  
}

// --------------------------------------------------------------------------------------
// Weighting functions needed
// 1. PU weighting
std::vector<TH1D*> puweights;

float puweight(float npu, int wset=0) {
  if ( npu  < 0                ) return 1.0;
  if ( wset >= (int) puweights.size() ) return 1.0;  
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
  
  if( wset >= puweights.size() ) puweights.resize(wset+1);
  puweights[wset] = new TH1D((*htargettmp)/(*hpumc));

  delete htargettmp;
  
}


// void setpuweights(TFile *file, TH1D *target, int wset=0) {
//   TDirectory *dirmcpv = (TDirectory*)file->FindObjectAny("AnaFwkMod");
//   TH1D *hnpu  = (TH1D*)dirmcpv -> Get("hNPU");
//   TH1D *hpumc = (TH1D*)hnpu    -> Clone();
//   hpumc->Sumw2();
//   hpumc->Scale(1.0/hpumc->GetSumOfWeights());  
//   if( wset >= puweights.size() ) puweights.resize(wset+1);
//   puweights[wset] = new TH1D((*target)/(*hpumc));
// }

// 2. PT weights
TH1D *ptweights = 0;
float ptweight(float genhpt, Int_t procid) {
  if (procid>0) return 1.0;
  if (genhpt<0) return 1.0;

  return 1.0;

  return ptweights->GetBinContent(ptweights->FindFixBin(genhpt));
}

// 3. Effcicnecy weights
std::vector<float> theEffWeights;
std::map<int,int> effCatToArrayMap;
float effweight( int cat ) {
  std::map<int,int>::iterator arrayInd = effCatToArrayMap.find( cat );
  if ( arrayInd == effCatToArrayMap.end() ) return 1.0;
  return theEffWeights[arrayInd->second];
}

std::vector<float> theTrigEffWeights;
std::map<int,int> trigCatToArrayMap;
float trigeffweight( int cat1, int cat2 ) {  
  int sumCat = 100*cat1+cat2;
  std::map<int,int>::iterator arrayInd = trigCatToArrayMap.find(sumCat);
  if ( arrayInd == trigCatToArrayMap.end() ) return 1.0;
  return theTrigEffWeights[arrayInd->second];
}
// --------------------------------------------------------------------------------------

//append data from file to RooDataSet adding a column which has the weight by given cross section
//(Note that the weight is not enabled at this stage, only made available for subsequent use)
void append(RooDataSet &data, TFile *infile, TCut sel, Double_t xsec, Int_t procid, TString modname);

bool readWeightsFromConfigCard(TString fileName, 
			       std::vector<float>* effWeight, std::vector<float>* trigWeights,
			       std::map<int,int> * effArray , std::map<int,int> * trigArray,
			       TH1D*& puhisto, TFile*& ptweight);

bool readParmsAndCatsFromConfigCard( TString fileName,
				     float& totallumi,
				     bool& computeMVAvar,
				     TString& mvaWeightFile,
				     TString& projectDir,
				     TString& datafilename,
				     TString& modname,
				     TString& treename,
				     TString& wsPrefix,
				     double& massmax, double& massmin,
				     std::map<TString,TString>& auxCats,
				     std::map<TString,TCut>& anaCats,
				     TCut& theBaseCut
				     );

bool readModelProcInfoFromCard(TString fileName,
			       TString& brLabel,
			       std::map<TString,TString>& procCSmap,
			       std::map<int,TString>&     procIdxMap,
			       std::map<TString,TString>& procMCfileMap,
			       std::vector<double>& mhs
			       );
// --------------------------------------------------------------------------------------

void makeworkspaceCard(TString inputCardName="template.config") {     
  
  gROOT->Macro("MitStyle.C");
  gStyle->SetErrorX(0); 
  gStyle->SetOptStat(0);
  gROOT->ForceStyle();  
  

  // 1. Set up the PU and pt-weights from input files
  TH1D*  hpuestnorm   = NULL;
  TFile* fileptweight = NULL;
  bool readStatus = readWeightsFromConfigCard( inputCardName       , 
					       &theEffWeights      ,
					       &theTrigEffWeights  ,
					       &effCatToArrayMap   ,
					       &trigCatToArrayMap  ,
					       hpuestnorm          ,
					       fileptweight         );
  
  
  if( !readStatus ) {
    std::cerr<<" ERROR: Could not read weight information from file "<<inputCardName.Data()<<std::endl;
    return;
  }

  if( ! fileptweight ) {
    std::cerr<<" ERROR: No pt-rewighting file loaded. Check config card for parameter PTREWEIGHFILE."<<std::endl;
    return;
  }

  if( !hpuestnorm ) {
    std::cerr<<" ERROR: No PU-rewighting histogram loaded. Check config card for parameter PUREWEIGHFILE."<<std::endl;
    return;
  }
  
  // ---------------------------------------------------------------------------------------------------------
  float      totallumi = -1.;

  double     massmax = -1.;
  double     massmin = -1.;

  TString    projectDir;
  TString    datafilename;
  TString    modname;
  TString    treename;
  TString    wsPrefix;


  std::map<TString,TString> auxCatMap;
  std::map<TString,TCut> anaCatMap;
  
  TCut theBaseCut;

  bool    computeMVAvar;
  TString mvaWeightFile;

  readStatus = readParmsAndCatsFromConfigCard( inputCardName       ,
					       totallumi           ,
					       computeMVAvar       ,
					       mvaWeightFile       ,
					       projectDir          ,
					       datafilename        ,
					       modname             ,
					       treename            ,
					       wsPrefix            ,
					       massmax, massmin,
					       auxCatMap           ,
					       anaCatMap           ,
					       theBaseCut
					       );
  

  if( !readStatus ) {
    std::cerr<<" ERROR: Could not read Category information from file "<<inputCardName.Data()<<std::endl;
    return;
  }
  // ---------------------------------------------------------------------------------------------------------  
  
  TString brLabel;
  std::map<TString,TString> procCSmap;
  std::map<int,TString>     procIdxMap;
  std::map<TString,TString> procMCfileMap;
  std::vector<double> mhs;

  readStatus = readModelProcInfoFromCard( inputCardName    ,				 
					  brLabel          ,
					  procCSmap        ,
					  procIdxMap       ,
					  procMCfileMap    ,
					  mhs
					  );

  if( !readStatus ) {
    std::cerr<<" ERROR: Could not read Model information from file "<<inputCardName.Data()<<std::endl;
    return;
  }
  // ---------------------------------------------------------------------------------------------------------  
  
  // test modle information....
  // initialize the map mass->index
  std::map<double,int> massArrayMap;
  initMassArrayMap(numsmpoints, smmasses, massArrayMap);
  std::map<TString,double*> processCrossSectionMap;
  initProcessCSArrayMap(processCrossSectionMap);

  std::map<TString,double*>::iterator it = processCrossSectionMap.find(brLabel);
  if ( it == processCrossSectionMap.end() ) {
    std::cerr<<" ERROR: Could not find BR array with name > "<<brLabel<<" < in file modelConstants.h."<<std::endl;
    return;
  }
  for( std::map<TString,TString>::iterator it2 = procCSmap.begin(); it2 != procCSmap.end(); ++it2) {
    it = processCrossSectionMap.find(it2->second);
    if ( it == processCrossSectionMap.end() ) {
      std::cerr<<" ERROR: Could not find CrossSection array with name > "<<it2->second<<" < in file modelConstants.h."<<std::endl;
      return;
    }
  }
  
  bool doff = false;
  if ( !gSystem->cd(projectDir.Data()) ) {
    std::cerr<<" ERROR: Could not change to project Directory > "<<projectDir.Data()<<" <. Directory does not exist."<<std::endl;
    return;
  }
  
  gStyle->SetOptStat(1110);    
  TCut puweight = "puweight(numPU)";
  
  //define categories
  std::vector<TString> catnames;
  std::vector<TCut>     catcuts;
  
  //define variables
  RooRealVar hmass("mass","m_{#gamma#gamma}",massmin,massmax,"GeV");
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
  RooRealVar corrpfmet("corrpfmet","",0.0);  
  RooRealVar leptonTag("leptonTag","",0.0);  

  int numBins = (int) (massmax - massmin) / 4;  
  hmass.setBins(numBins);


  // if required, set up the DiPhoton MVA
  if( computeMVAvar ) {

    g_reader = new TMVA::Reader();
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
  
    g_reader->BookMVA("Ystar",mvaWeightFile.Data());
  }

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
  varlist.add(corrpfmet);
  varlist.add(leptonTag);
  varlist.add(ptgg);  
  
  //extra variables for reweighting
  RooRealVar xsecweight("xsecweight","xsecweight",1.0);
  RooRealVar procidx("procidx","",0.0);  
  
  RooArgSet varlistw = varlist;
  varlistw.add(xsecweight);
  varlistw.add(procidx);
  
  RooFormulaVar newmassf("CMS_hgg_mass","m_{#gamma#gamma}","mass",RooArgList(hmass));
  RooFormulaVar* bdtf = NULL;
  if ( computeMVAvar )
    bdtf = new RooFormulaVar("bdt","","getmvaval(masserrsmeared/mass,masserrsmearedwrongvtx/mass,vtxprob,ph1.pt/mass,ph2.pt/mass,ph1.eta,ph2.eta,TMath::Cos(ph1.phi-ph2.phi),ph1.idmva,ph2.idmva)",varlist);
  
  TFile *datafile = new TFile(datafilename.Data(),"READ");
  if ( !datafile ) {
    std::cerr<<" ERROR: Could not open datafile with name > "<<datafilename.Data()<<" <."<<std::endl;
    return;
  }
  TDirectory *datadir = (TDirectory*) datafile->FindObjectAny(modname);
  if ( !datadir ) {
    std::cerr<<" ERROR: Could not find directory > "<<modname.Data()<<" < in datafile with name > "<<datafilename.Data()<<" <."<<std::endl;
    return;
  }
  TTree *hdata = (TTree*)datadir->Get(treename.Data());  
  if ( !hdata ) {
    std::cerr<<" ERROR: Could not find Tree > "<<treename.Data()<<" < in directory > "<<modname.Data()<<" in datafile with name > "<<datafilename.Data()<<" <."<<std::endl;
    return;
  }
  
  RooDataSet calldata("calldata","",varlist,RooFit::Import(*hdata),RooFit::Cut(theBaseCut)); 
  RooRealVar *newmass = (RooRealVar*)calldata.addColumn(newmassf);

  RooRealVar *bdt = NULL;
  if ( computeMVAvar )
    bdt = (RooRealVar*)calldata.addColumn(*bdtf);
  
  newmass->setRange(100.0,180.0);
  newmass->setUnit("GeV");
  newmass->setBins(320);
  
  RooWorkspace *w = 0;
  RooWorkspace *wdata = new RooWorkspace("cms_hgg_workspace_data","") ;
  RooWorkspace *wmclow = new RooWorkspace("cms_hgg_workspace_mclow","") ;
  RooWorkspace *wmchigh = new RooWorkspace("cms_hgg_workspace_mchigh","") ;

  for( std::map<TString,TCut>::iterator ttIt = anaCatMap.begin(); ttIt != anaCatMap.end(); ++ttIt) {
    TString catname = ttIt->first; //catnames.at(icat);
    RooDataSet cdata(TString::Format("data_mass_%s",catname.Data()),"",*calldata.get(),RooFit::Import(calldata),RooFit::Cut( ttIt->second )); 

    wdata->import(cdata);
  }
  
  RooRealVar *IntLumi = new RooRealVar("IntLumi","",totallumi);
  IntLumi->setConstant();

  wdata->import(*IntLumi);  
  wmclow->import(*IntLumi);  
  wmchigh->import(*IntLumi);  
  
  for (UInt_t i=0; i<mhs.size(); ++i) {
    
    //     int masspoint = mhs.at(i);  
    double masspoint = mhs.at(i);
    
    if ( masspoint < 122.5 ) w = wmclow;
    else w = wmchigh;
    
    double gfmasspoint = masspoint;
    if (masspoint==150.) gfmasspoint = 145.;
    
    //load correct pt weights
    ptweights= (TH1D*) fileptweight->Get(TString::Format("kfact%d_0",(int)masspoint));
    if ( !ptweights ) {
      std::cerr<<" ERROR: No pt-weight histogram for mass point "<<(int) masspoint<<" in Kfactors file."<<std::endl;
      return;
    }

    double totalxsec = 0.0;

    TCut rightvtx = "abs(genHiggsZ-vtxZ)<1.0";
    TCut wrongvtx = "abs(genHiggsZ-vtxZ)>=1.0";

    TCut rightcut =  theBaseCut && rightvtx;
    TCut wrongcut =  theBaseCut && wrongvtx;
    TCut allcut   =  theBaseCut;
    
    //right vertex
    RooDataSet mcsigdata(TString::Format("mcsigdata_m%d",(int) masspoint),"",varlistw);
    //wrong vertex
    RooDataSet mcsigwrongdata(TString::Format("mcsigwrongdata_m%d",(int) masspoint),"",varlistw);    
    //combined
    RooDataSet mcsigalldata(TString::Format("mcsigalldata_m%d",(int) masspoint),"",varlistw);

    //define aggregate weight, so far using xsec, pileup, pt-reweighting and efficiency scale factors
    RooArgList weightvarlist(*IntLumi,xsecweight,numPU,higgspt,procidx,ph1Cat,ph2Cat,ph1pt,ph2pt);
    weightvarlist.add(ph1r9);
    weightvarlist.add(ph2r9);
    
    RooFormulaVar totweight("totweight","totweight","IntLumi*xsecweight*puweight(numPU,procidx)*ptweight(genHiggspt,procidx)*effweight(ph1.phcat)*effweight(ph2.phcat)*trigeffweight(ph1.phcat,ph2.phcat)",weightvarlist);

    
    for( std::map<int,TString>::iterator tIt = procIdxMap.begin(); tIt != procIdxMap.end(); ++tIt ) {
      int theProcIdx   = tIt->first;
      TString procName = tIt->second;

      std::map<TString,TString>::iterator tempIt = procMCfileMap.find(procName);
      if( tempIt == procMCfileMap.end() ) {
	std::cout<< " ERROR *** "<<std::endl;
	return;
      }

      TString samplestring = tempIt->second;
      
//       if ( !strcmp( (procName).Data(), "wzh" ) )
// 	if ( (int) masspoint==110 || (int) masspoint==130 || (int) masspoint==140) samplestring.ReplaceAll("f11--","f11-");
      
      //TFile *theFile = new TFile(TString::Format(samplestring, (int) masspoint),"READ");
      int masspointDummy=110;
      TFile *theFile = new TFile(TString::Format(samplestring, (int) masspointDummy),"READ");
      if (!theFile) {
	std::cerr<<" ERROR: Could not open file "<<TString::Format(samplestring, (int) masspoint).Data()<<"."<<std::endl;
	return;
      }


      setpuweights(theFile,hpuestnorm,theProcIdx);

      std::map<double,int>::iterator maIdx = massArrayMap.find(masspoint);
      if ( maIdx == massArrayMap.end() ) {
	std::cerr<<" ERROR: Could not find cross-section for mass "<<masspoint<<"."<<std::endl;
	return;
      }
      
      std::map<TString,TString>::iterator tempIt3 = procCSmap.find(procName);
      if ( tempIt3 == procCSmap.end() ) {
	std::cout<<" ERROR: Cound not find XSlist for process "<<procName<<std::endl;
	return;
      }
      
      std::map<TString,double*>::iterator tempIt2 = processCrossSectionMap.find(tempIt3->second);
      if ( tempIt2 == processCrossSectionMap.end() ) {
	std::cout<<" ERROR: Cound not find XS for process "<<tempIt3->second<<std::endl;
	return;
      }

      double theXSBR = (tempIt2->second)[maIdx->second];
      totalxsec += theXSBR;

      RooRealVar *XSBR_proc = new RooRealVar(TString::Format("XSBR_%s_%d",procName.Data(), (int)masspoint),"",theXSBR);
      XSBR_proc->setConstant();
      w->import(*XSBR_proc); 
      
      append(mcsigdata     ,theFile,rightcut,theXSBR, theProcIdx, modname);
      append(mcsigwrongdata,theFile,wrongcut,theXSBR, theProcIdx, modname);
      append(mcsigalldata  ,theFile,allcut  ,theXSBR, theProcIdx, modname);
    }

    mcsigdata.addColumn(totweight);
    mcsigdata.addColumn(newmassf);
    if( computeMVAvar )
      mcsigdata.addColumn(*bdtf);

    mcsigwrongdata.addColumn(totweight);
    mcsigwrongdata.addColumn(newmassf);
    if( computeMVAvar )
      mcsigwrongdata.addColumn(*bdtf);    

    mcsigalldata.addColumn(totweight);    
    mcsigalldata.addColumn(newmassf);
    if( computeMVAvar )
      mcsigalldata.addColumn(*bdtf);

    RooRealVar *XSBR = new RooRealVar(TString::Format("XSBR_%d",(int) masspoint),"",totalxsec);
    XSBR->setConstant();
    w->import(*XSBR);  
              
    //loop over Categories
    //for (UInt_t icat=0; icat<catcuts.size(); ++icat) {
    for(std::map<TString,TCut>::iterator catIt = anaCatMap.begin(); catIt != anaCatMap.end(); ++catIt) {

      TString catname = catIt->first;
      RooDataSet *mcsigwdata      = new RooDataSet(TString::Format("sig_mass_m%d_rv_%s",(int) masspoint,catname.Data()),"",*mcsigdata.get()     ,RooFit::Import(mcsigdata)     ,RooFit::Cut(catIt->second),RooFit::WeightVar("totweight"));
      RooDataSet *mcsigwrongwdata = new RooDataSet(TString::Format("sig_mass_m%d_wv_%s",(int) masspoint,catname.Data()),"",*mcsigwrongdata.get(),RooFit::Import(mcsigwrongdata),RooFit::Cut(catIt->second),RooFit::WeightVar("totweight"));
      RooDataSet *mcsigallwdata   = new RooDataSet(TString::Format("sig_mass_m%d_%s",(int) masspoint,catname.Data()),"",*mcsigalldata.get()  ,RooFit::Import(mcsigalldata)  ,RooFit::Cut(catIt->second),RooFit::WeightVar("totweight"));

      for( std::map<int,TString>::iterator tIt = procIdxMap.begin(); tIt != procIdxMap.end(); ++tIt ) {

	int procCounter = tIt->first;

	RooAbsData *mcsigwdata_proc = mcsigwdata->reduce(TCut(TString::Format("procidx==%d",procCounter)));       
	mcsigwdata_proc->SetName(TString::Format("sig_%s_mass_m%d_rv_%s",(tIt->second).Data(),(int) masspoint,catname.Data()));
	
	RooAbsData *mcsigwrongwdata_proc = mcsigwrongwdata->reduce(TCut(TString::Format("procidx==%d",procCounter)));       
	mcsigwrongwdata_proc->SetName(TString::Format("sig_%s_mass_m%d_wv_%s",(tIt->second).Data(),(int) masspoint,catname.Data()));
	
	RooAbsData *mcsigallwdata_proc = mcsigallwdata->reduce(TCut(TString::Format("procidx==%d",procCounter)));       
	mcsigallwdata_proc->SetName(TString::Format("sig_%s_mass_m%d_%s",(tIt->second).Data(),(int) masspoint,catname.Data()));

	w->import(*mcsigwdata_proc);
	w->import(*mcsigwrongwdata_proc);
	w->import(*mcsigallwdata_proc);
      }

      w->import(*mcsigwdata);
      w->import(*mcsigwrongwdata);
      w->import(*mcsigallwdata);	

    }
  }

  wdata->writeToFile(TString::Format("%s-data.root",wsPrefix.Data())) ;
  wmclow->writeToFile(TString::Format("%s-mclow.root",wsPrefix.Data())) ;
  wmchigh->writeToFile(TString::Format("%s-mchigh.root",wsPrefix.Data())) ;
  
  // go back to original directory
  gSystem->cd( ".." );
  return;
  
}


// -------------------------------------------------------------------------------------------------------------------------
// helper-function implementations
void append(RooDataSet &data, TFile *infile, TCut sel, Double_t xsec, Int_t procid, TString modname) {

  TDirectory *hdirfwk = (TDirectory*) infile->FindObjectAny("AnaFwkMod");
  const TH1D *hDAllEvents = (TH1D*)hdirfwk->Get("hDAllEvents");
  
  TDirectory *hdir = (TDirectory*) infile->FindObjectAny(modname);
  TTree *hdata = (TTree*)hdir->Get("hPhotonTree");  
    
  RooRealVar xsecweight("xsecweight","xsecweight",xsec/(double)hDAllEvents->GetEntries());
  RooRealVar procidx("procidx","procidx",procid);
  
  RooArgSet varlistsmall = *data.get();
  varlistsmall.remove(xsecweight,kFALSE,kTRUE);
  varlistsmall.remove(procidx,kFALSE,kTRUE);
  
  RooDataSet newdata("newdata","newdata",varlistsmall,RooFit::Import(*hdata),RooFit::Cut(sel)); 
  newdata.addColumn(xsecweight);
  newdata.addColumn(procidx);
  
  data.append(newdata);        
  return;
}

bool readWeightsFromConfigCard(TString fileName, 
			       std::vector<float>* effWeights, std::vector<float>* trigWeights,
			       std::map<int,int> * effArray  , std::map<int,int> * trigArray,
			       TH1D*& puweight, TFile*& ptfile) {
  
  FILE* configFile = fopen(fileName.Data(),"r");
  if ( !configFile ) {
    std::cerr<<" Inputfile "<<fileName<<" not found."<<std::endl;
    return false;
  }
  
  char line[200];
  
  std::cout<<" Reading weight information from file "<<fileName<<"...";
  
  while (fgets(line,200,configFile)) {
    // test number of photon categories
    int numProcs = -1;
    int numPhCats   = -1;
    int numAuxCats  = -1;
    int numAnaCats  = -1;
    int numTrigCats = -1;
    float effScale  = -1.;
    float vetoScale = -1.;
    float trigEff   = -1.;
    int thePhCatIdx     = -1;
    int thePhCatLabel   = -1;
    int thePhCat2Label  = -1;
    int theTrigCatIdx   = -1;

    char puWeightFileName[150];
    char ptWeightFileName[150];
    
    if(line[0] == '#') continue;

    if( sscanf(line,"PUREWEIGHFILE %s",&puWeightFileName) ) {
      //   //Load pileup weights
      TFile *filepuest = new TFile(puWeightFileName);
      if( !filepuest ) {
	std::cerr<<" ERROR: Could not open PU file with name "<<puWeightFileName<<"."<<std::endl;
	fclose(configFile);
	return false;
      }
      TH1D *hpuest = (TH1D*) filepuest->Get("pileup");
      if( !hpuest ) {
	std::cerr<<" ERROR: Could not read histogram <pileup> from PU file with name "<<puWeightFileName<<"."<<std::endl;
	fclose(configFile);
	return false;
      }

      puweight = hpuest;

//       puweight = new TH1D("hNPU", "hNPU", 51, -0.5, 50.5);
//       for (int i=0; i<51; ++i) {
// 	puweight->Fill(i,hpuest->GetBinContent(hpuest->GetXaxis()->FindFixBin(i)));
//       }  
//       puweight->Sumw2();
//       puweight->Scale(1.0/puweight->GetSumOfWeights());

    } else if( sscanf(line,"PTREWEIGHFILE %s",&ptWeightFileName) ) {      
      //load pt-weights
      ptfile = new TFile(ptWeightFileName,"READ");
      if( !ptfile ) {
	std::cerr<<" ERROR: Could not open PU file with name "<<ptWeightFileName<<"."<<std::endl;
	fclose(configFile);
	return false;
      }
    } else if( sscanf(line,"INIT %d %d %d %d %d", &numProcs, &numPhCats, &numAuxCats, &numAnaCats, &numTrigCats) ) {
      if (numPhCats*numPhCats != numTrigCats) {
	std::cerr<<" ERROR: Number of trigger Cats must be Number of photon cats squared for now."<<std::endl;
	fclose(configFile);
	return false;
      }
      effWeights ->resize(numPhCats);
      trigWeights->resize(numTrigCats);
      effArray   ->clear();
      trigArray  ->clear();
    } else if( sscanf(line,"PHOCAT %d %d %f %f", &thePhCatIdx, &thePhCatLabel, &effScale, &vetoScale) ) {
      if ( effWeights->size() == 0 || thePhCatIdx >= (int) effWeights->size() ) {
	std::cerr<<" ERROR: Photon Categories not properly set up: NCat = "<<effWeights->size()<<" and asking for index idx = "<<thePhCatIdx<<"."<<std::endl;
	fclose(configFile);
	return false;
      }
      (*effWeights)[thePhCatIdx] = effScale*vetoScale;
      effArray->insert(std::pair<int,int>(thePhCatLabel,thePhCatIdx));
    } else if( sscanf(line,"TRIGCAT %d %d %d %f", &theTrigCatIdx, &thePhCatLabel, &thePhCat2Label, &trigEff) ) {
      if ( trigWeights->size() == 0 || theTrigCatIdx >= (int) trigWeights->size() ) {
	std::cerr<<" ERROR: Trigger Categories not properly set up: NCat = "<<trigWeights->size()<<" and asking for index idx = "<<theTrigCatIdx<<"."<<std::endl;
	fclose(configFile);
	return false;
      }
      (*trigWeights)[theTrigCatIdx] = trigEff;
      int sumCats = 100*thePhCatLabel+thePhCat2Label;
      trigArray->insert(std::pair<int,int>(sumCats,theTrigCatIdx));
    }
  }

  std::cout<<" done."<<std::endl;
  fclose(configFile);
  
  return true;
  
}

bool readModelProcInfoFromCard(TString fileName,
			       TString& brLabel,
			       std::map<TString,TString>& procCSmap,
			       std::map<int,TString>&     procIdxMap,
			       std::map<TString,TString>& procMCfileMap,
			       std::vector<double>& mhs
			       ) {
  
  FILE* configFile = fopen(fileName.Data(),"r");
  if ( !configFile ) {
    std::cerr<<" Inputfile "<<fileName<<" not found."<<std::endl;
    return false;
  }
  
  char line[400];
  
  std::cout<<" Reading Model and Process information from file "<<fileName<<"...";
  
  mhs.resize(0);

  while (fgets(line,400,configFile)) {
    if(line[0] == '#') continue;    
    // we only care for model 0 in this macro...
    char modelName[100];
    char brName[100];
    char procName[100];
    char procCSLabel[100];

    char MCfileName[200];

    float mass = -100.;
    int modelIdx = -1;
    int procIdx  = -1;
    int massIdx  = -1;

    if ( sscanf(line,"MODEL %d %s %s",&modelIdx, &modelName, &brName) ) {
      if(modelIdx == 0) brLabel = TString(brName);
    } else if ( sscanf(line,"PROC %d %s %s file:%s",&procIdx, &procName, &procCSLabel, &MCfileName) ) {
      //# ------------------------------------------------------------------------------------------------------------------------------------------------
      //PROC    0	ggh	OFF			file:/scratch/fabstoec/cms/hist/hgg-7TeV-janReReco/merged/hgg-7TeV-janReReco_f11--h%dgg-gf-v14b-pu_noskim.root
      std::cout<<" adding poroc with idx = "<<procIdx<<"  and name "<<procName<<std::endl;
      procMCfileMap.insert(std::pair<TString,TString>(TString(procName),TString(MCfileName)));
      procIdxMap   .insert(std::pair<int,TString>    (procIdx          ,TString(procName)));
    } else if ( sscanf(line,"MODPROS %d %d %s",&modelIdx, &procIdx, &procCSLabel) ) {
      if( procIdxMap.find(procIdx) ==  procIdxMap.end() ) {
	std::cerr<<" ERROR: Cannot access process with idx = "<<procIdx<<". Process not defined in PROC section."<<std::endl;
	return false;
      }      
      procCSmap.insert(std::pair<TString,TString>(   procIdxMap.find(procIdx)->second, TString(procCSLabel)  ));
      std::cout<<"inserting proc with idx = "<<procIdx<<"  name = "<<procIdxMap.find(procIdx)->second<<"   and label = "<<TString(procCSLabel)<<std::endl;
    } else if ( sscanf(line,"MASS %d %f",&massIdx, &mass) ) {
      if( massIdx > mhs.size() ) {
	std::cerr<<" ERROR: Ordering in MASS section wrong. Must be idx continuous."<<std::endl;
	return false;
      }
      mhs.push_back( (double) mass );
    }
  }
  
  std::cout<<" done."<<std::endl;
  fclose(configFile);
  return true;      
}

bool readParmsAndCatsFromConfigCard( TString fileName,
				     float& totallumi,
				     bool&    computeMVAvar,
				     TString& mvaWeightFile,
				     TString& projectDir,
				     TString& datafilename,
				     TString& modname,
				     TString& treename,
				     TString& wsPrefix,
				     double& massmax, double& massmin,
				     std::map<TString,TString>& auxCats,
				     std::map<TString,TCut>& anaCats,
				     TCut& baseCut
				     ) {

  FILE* configFile = fopen(fileName.Data(),"r");
  if ( !configFile ) {
    std::cerr<<" Inputfile "<<fileName<<" not found."<<std::endl;
    return false;
  }
  
  char line[200];
  computeMVAvar = false;
  std::cout<<" Reading paramater from file "<<fileName<<"...";
  
  while (fgets(line,200,configFile)) {

    if(line[0] == '#') continue;

    char name[200];
    int catIdx = -1;
    char catName[200];
    char theCat[200];
    float lumi    = -1.;
    float massval = -1.;
    if ( sscanf(line,"PROJECTDIR %s",&name ) ) projectDir = TString(name);
    else if  ( sscanf(line,"LUMI %f",&lumi) ) totallumi = lumi;
    else if  ( sscanf(line,"DATAFILE %s",&name) ) datafilename = TString(name);
    else if  ( sscanf(line,"MODNAME %s",&name) ) modname = TString(name);
    else if  ( sscanf(line,"TREENAME %s",&name) ) treename = TString(name);
    else if  ( sscanf(line,"WS1PREFIX %s",&name) ) wsPrefix = TString(name);
    else if( sscanf(line,"MINMSS %f",&massval) ) massmin = (double) massval;
    else if( sscanf(line,"MAXMSS %f",&massval) ) massmax = (double) massval;
    else if  ( sscanf(line,"AUXCAT %d %s",&catIdx, &catName) ) {
      // parsing the ctegory definition like:
      // AUXCAT	0	masscut		" mass>100.0 && mass<180. "
      std::string totLine = line;
      int startCat = totLine.find_first_of("\"");
      int endCat   = totLine.find_last_of("\"");
      std::string catLine = totLine.substr(startCat+1, endCat-startCat-1);
      auxCats.insert(std::pair<TString,TString>(TString(catName),TString(catLine.c_str())));
    }
    else if  ( sscanf(line,"COMPUTEMVA ON %s",&name) ) { mvaWeightFile = TString(name); computeMVAvar = true; }
    else if  ( sscanf(line,"ANACAT %d %s",&catIdx, &catName) ) {
      // ANACAT	0	cat0		" basecut && !vbfcut && baseline0 "		0.005432		Bern/5
      std::string totLine = line;
      int startCat = totLine.find_first_of("\"");
      int endCat   = totLine.find_last_of("\"");
      std::string catLine = totLine.substr(startCat+1, endCat-startCat-1);
      
      // erase starting empties...
      int fPos = catLine.find_first_not_of(" ");
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
    else if  ( sscanf(line,"BASECAT %s",&name) ) {
      std::string totLine = line;
      int startCat = totLine.find_first_of("\"");
      int endCat   = totLine.find_last_of("\"");
      std::string catLine = totLine.substr(startCat+1, endCat-startCat-1);
      
      // erase starting empties...
      int fPos = catLine.find_first_not_of(" ");
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
