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
#include "RooChebychev.h"
#include "RooBernstein.h"
#include "TLatex.h"
#include "TGraphAsymmErrors.h"
#include "RooMinimizer.h"
#include "RooStats/RooStatsUtils.h"
#include "RooAddition.h"


#include "TTreeFormula.h"

#include "TString.h"

#include "TMVA/Reader.h"
#include "TMVA/Tools.h"
#include "TMVA/Config.h"

using namespace RooFit;

#include "FstModels.h"     // some helpfer functions for different models (well... I cut away all but Bernstein, but can be brought back...)

// RooDataSet *makedset(TString name, TTree *tree, TCut cut, RooRealVar *hmass, RooRealVar* sbweight = NULL, Double_t tWeight = 1.,
// 		     std::map<int,std::map<int,std::vector<int>*>*>* eventMap = NULL,
// 		     int* numSkipped = NULL) {
  
//   RooDataSet *dset = NULL;
//   if( sbweight )
//     dset = new RooDataSet(name,"",RooArgSet(*hmass,*sbweight),"sbweight");
//   else
//     dset = new RooDataSet(name,"",RooArgSet(*hmass));
  
//   RooRealVar *vset = (RooRealVar*)dset->get()->find(hmass->GetName());
  

//   std::cout<<"  ======================================= "<<std::endl;
//   std::cout<<"  "<<cut<<std::endl;

//   tree->SetEstimate(tree->GetEntries());

//   Int_t nev = tree->Draw("mass:run:lumi:evt",cut,"goff");    
//   std::cout<<"  Nev = "<<nev<<std::endl;

//   double *vals  = tree->GetV1();
//   double *runs  = tree->GetV2();
//   double *lumis = tree->GetV3();
//   double *evts  = tree->GetV4();
//   //double *weights = tree->GetW();  // this is Data, has no weights.
  
//   for (int iev=0; iev<nev; ++iev) {
//     vset->setVal(vals[iev]);

//     bool noskip = true;  // if false, we skip event
//     if ( eventMap ) {  // remove events in the map
      
//       noskip = ( eventMap->find( (int) runs[iev] ) == eventMap->end() || 
// 		 eventMap->find( (int) runs[iev] ) -> second->find( (int) lumis[iev] ) == eventMap->find( (int) runs[iev] ) ->second ->end() );
//       if ( !noskip ){
// 	noskip = true;
// 	for(std::vector<int>::iterator it = eventMap->find( (int) runs[iev] ) -> second->find( (int) lumis[iev] ) -> second ->begin(); 
// 	    it != eventMap->find( (int) runs[iev] ) -> second->find( (int) lumis[iev] ) -> second ->end(); ++it ){
// 	  if (evts[iev] == (*it)) {
// 	    noskip = false;
// 	    if(numSkipped) (*numSkipped)++;
// 	    std::cout<<" BAD ECAL Laser List: Skipping event "<<runs[iev]<<"  "<<lumis[iev]<<"  "<<evts[iev]<<std::endl;
// 	    break;
// 	  }
// 	}
//       }
//     }
      
//     if (noskip)
//       dset->add(RooArgSet(*vset),tWeight);
//   }
  
//   return dset;
  
// }

RooDataSet *makedset(TString name, TTree **trees, int numTrees, TCut cut, RooRealVar *hmass, RooRealVar* sbweight = NULL, Double_t tWeight = 1.,
		     int catCounter = -1,
		     bool dumpEvents = false,
		     std::map<int,std::map<int,std::vector<int>*>*>* eventMap = NULL,
		     int* numSkipped = NULL) {
  
  RooDataSet *dset = NULL;
  if( sbweight )
    dset = new RooDataSet(name,"",RooArgSet(*hmass,*sbweight),"sbweight");
  else
    dset = new RooDataSet(name,"",RooArgSet(*hmass));
  
  RooRealVar *vset = (RooRealVar*)dset->get()->find(hmass->GetName());
  

  std::cout<<"  ======================================= "<<std::endl;
  std::cout<<"  "<<cut<<std::endl;

  for(int iTree = 0; iTree < numTrees; ++iTree ) {
    TTree* tree = trees[iTree];

    tree->SetEstimate(tree->GetEntries());

    Int_t nev = tree->Draw("mass:run:lumi:evt",cut,"goff");    
    std::cout<<"  Nev = "<<nev<<std::endl;

    double *vals  = NULL;
    double *runs  = NULL;
    double *lumis = NULL;
    double *evts  = NULL;
    //double *weights = tree->GetW();  // this is Data, has no weights.

    double* rhos     = NULL;
    double* nvtx     = NULL;
    double* bdt      = NULL;
    
    double* ph1pt    = NULL;
    double* ph1sceta = NULL;
    double* ph1idmva = NULL;
    double* ph1r9    = NULL;
    
    double* ph2pt    = NULL;
    double* ph2sceta = NULL;
    double* ph2idmva = NULL;
    double* ph2r9    = NULL;
    

    double* cosDphi  = NULL;
    double* sigMoM   = NULL;
    double* sigMoM_wv= NULL;
    double* pvtx     = NULL;

  
    if( dumpEvents ) {
      vals  = new double[nev];
      runs  = new double[nev];
      lumis = new double[nev];
      evts  = new double[nev];

      double* tmp1 = tree->GetV1();
      double* tmp2 = tree->GetV2();
      double* tmp3 = tree->GetV3();
      double* tmp4 = tree->GetV4();

      for (int iev=0; iev<nev; ++iev) {
	vals[iev] = tmp1[iev];
	runs[iev] = tmp2[iev];
	lumis[iev] = tmp3[iev];
	evts[iev] = tmp4[iev];
      }


      rhos     = new double[nev];
      nvtx     = new double[nev];
      bdt      = new double[nev];
      
      ph1pt    = new double[nev];
      ph1sceta = new double[nev];
      ph1idmva = new double[nev];
      ph1r9    = new double[nev];
      
      ph2pt    = new double[nev];
      ph2sceta = new double[nev];
      ph2idmva = new double[nev];
      ph2r9    = new double[nev];

      cosDphi  = new double[nev];
      sigMoM   = new double[nev];
      sigMoM_wv= new double[nev];
      pvtx     = new double[nev];



      nev = tree->Draw("ph1.pt:ph1.sceta:ph1.idmva:bdt",cut,"goff");
      tmp1   = tree->GetV1();
      tmp2   = tree->GetV2();
      tmp3   = tree->GetV3();
      tmp4   = tree->GetV4();
      for (int iev=0; iev<nev; ++iev) {
	ph1pt[iev] = tmp1[iev];
	ph1sceta[iev] = tmp2[iev];
	ph1idmva[iev] = tmp3[iev];
	bdt[iev]      = tmp4[iev];
      }

      nev = tree->Draw("ph2.pt:ph2.sceta:ph2.idmva",cut,"goff");
      tmp1   = tree->GetV1();
      tmp2   = tree->GetV2();
      tmp3   = tree->GetV3();
      for (int iev=0; iev<nev; ++iev) {
	ph2pt[iev] = tmp1[iev];
	ph2sceta[iev] = tmp2[iev];
	ph2idmva[iev] = tmp3[iev];
      }
      
      nev = tree->Draw("ph1.r9:ph2.r9:rho:nVtx",cut,"goff");
      tmp1   = tree->GetV1();
      tmp2   = tree->GetV2();
      tmp3   = tree->GetV3();
      tmp4   = tree->GetV3();
      for (int iev=0; iev<nev; ++iev) {
	ph1r9[iev] = tmp1[iev];
	ph2r9[iev] = tmp2[iev];
	rhos[iev]  = tmp3[iev];
	nvtx[iev]  = tmp4[iev];
      }

      nev = tree->Draw("masserrsmeared/mass:masserrsmearedwrongvtx/mass:TMath::Cos(ph1.scphi-ph2.scphi):vtxprob",cut,"goff");
      tmp1   = tree->GetV1();
      tmp2   = tree->GetV2();
      tmp3   = tree->GetV3();
      tmp4   = tree->GetV4();
      for (int iev=0; iev<nev; ++iev) {
	ph1r9[iev] = tmp1[iev];

	sigMoM   [iev]=tmp1[iev];
	sigMoM_wv[iev]=tmp2[iev];
	cosDphi  [iev]=tmp3[iev];
	pvtx     [iev]=tmp4[iev];

      }
    } else {

      vals  = tree->GetV1();
      runs  = tree->GetV2();
      lumis = tree->GetV3();
      evts  = tree->GetV4();
      
    }

    for (int iev=0; iev<nev; ++iev) {
      vset->setVal(vals[iev]);
      
      bool noskip = true;  // if false, we skip event
      if ( eventMap ) {  // remove events in the map
	
	noskip = ( eventMap->find( (int) runs[iev] ) == eventMap->end() || 
		   eventMap->find( (int) runs[iev] ) -> second->find( (int) lumis[iev] ) == eventMap->find( (int) runs[iev] ) ->second ->end() );
	if ( !noskip ){
	  noskip = true;
	  for(std::vector<int>::iterator it = eventMap->find( (int) runs[iev] ) -> second->find( (int) lumis[iev] ) -> second ->begin(); 
	      it != eventMap->find( (int) runs[iev] ) -> second->find( (int) lumis[iev] ) -> second ->end(); ++it ){
	    if (evts[iev] == (*it)) {
	      noskip = false;
	      if(numSkipped) (*numSkipped)++;
	      std::cout<<"  Skipping even "<<runs[iev]<<"  "<<lumis[iev]<<"  "<<evts[iev]<<std::endl;
	      break;
	    }
	  }
	}
      }
      
      if (noskip) {
	dset->add(RooArgSet(*vset),tWeight);
	
	if( dumpEvents ) {

	  int cicClass = 0;
	  if( TMath::Abs(ph1sceta[iev]) > 1.5 || TMath::Abs(ph2sceta[iev]) > 1.5 ) cicClass+=2;
	  if( ph1r9[iev] < 0.94 || ph2r9[iev] < 0.94 ) cicClass++;


 	  printf("type:0\trun:%d\tlumi:%d\tevent:%d\tnvtx:%d\trho:%.4e\tmgg:%.4e\tr9_1:%.4e\tsceta_1:%.4e\tidmva_1:%.4e\tptom_1:%.4e\tr9_2:%.4e\tsceta_2:%.4e\tidmva_2:%.4e\tptom_2:%.4e\tbdt:%.4e\tciccat:%d\tevcat:%d\tcosdphi:%.4e\tsigmom_rv:%.4e\tsigmom_wv:%.4e\tvtxprob:%.4e\n",
 		 (int) runs[iev], (int) lumis[iev],(unsigned int) evts[iev],
		 (int) nvtx[iev], rhos[iev],vals[iev],		 
 		 ph1r9[iev],ph1sceta[iev],ph1idmva[iev],ph1pt[iev]/vals[iev],
 		 ph2r9[iev],ph2sceta[iev],ph2idmva[iev],ph2pt[iev]/vals[iev],
 		 bdt[iev],cicClass,catCounter,
		 cosDphi  [iev],
		 sigMoM   [iev],
		 sigMoM_wv[iev],
		 pvtx     [iev]
		 );		 
	  
	  // 	  printf("type:0\trun:%d\tlumi:%d\tevent:%d\tmass:%.4e\tph1_pt:%.4e\t%s\n",
	  // 		 (int) runs[iev], (int) lumis[iev],(int) evts[iev], vals[iev],
// 		 ph1pt[iev],
// 		 name.Data());		 
	  
	  
	}
      }
    }
  }
  
  return dset;
  
}


TTree *ApplyAsFriend(TTree *intree, TString tmvaweights, const std::vector<std::string> &vars, std::string targetname)
{
  
  int nvars = vars.size();
    
  //initialize TTreeFormulas to read variables from TTree
  std::vector<TTreeFormula*> inputforms;
  for (std::vector<std::string>::const_iterator it = vars.begin(); 
      it != vars.end(); ++it) {
    inputforms.push_back(new TTreeFormula(it->c_str(),it->c_str(),intree));
  }
  
  Float_t target = 0.;
  Float_t *vals = new Float_t[nvars];
  
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
    
    for (int i=0; i<nvars; ++i) {
      if (newtree) inputforms[i]->Notify();
      vals[i] = inputforms[i]->EvalInstance();
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

// ----------------------------------------------------------------------------------------------------
bool readFromConfigCard( TString fileName,
			 float& totallumi,
			 bool& computeMVAvar,
			 TString& mvaWeightFile,
			 TString& mvaDefFile,
			 TString& projectDir,
			 //TString& datafilename,
			 std::vector<std::string>& datafilenames,
			 std::vector<TString>& procNameList,
			 TString& modname,
			 TString& treename,
			 TString& wsPrefix,
			 double& massmax, double& massmin,
			 std::map<TString,TString>& auxCats,
			 std::map<TString,TCut>   & anaCats,
			 std::map<TString,TString>& catDesc ,		
			 std::map<TString,int>    & polOrder,
			 TCut& theBaseCut,
			 double& theCMenergy
			 );
// ----------------------------------------------------------------------------------------------------

void createECALFilterMap(std::string fileName, std::map<int,std::map<int,std::vector<int>* >* >* map) {

  FILE* file = fopen(fileName.c_str(),"r");
  char line[100];
  int run, lumi, event;

  int lastrun  = -1;
  int lastlumi = -1;
  
  while( fgets(line,100,file) ) {
    if( sscanf(line,"%d %d %d",&run, &lumi, &event) ){            
      if(run != lastrun)
	map->insert(std::pair<int, std::map<int, std::vector<int>*>*>(run, new std::map<int,std::vector<int>*>()));
      if(lumi != lastlumi)
	map->find(run)->second->insert(std::pair<int,std::vector<int>*>(lumi, new std::vector<int>()));

      map->find(run)->second->find(lumi)->second->push_back(event);
      
      lastrun  = run;
      lastlumi = lumi;
    }
  }
  
  return;
}


#define DOBDT

void fitbkgdataCard_HGG(bool dobands  = false, 
#ifdef DOBDT
			TString configCard="templateHGG_8TeV_Moriond.config", 
#else
			TString configCard="templateHGG_8TeV_Moriond_CiC.config", 
#endif
			bool dosignal     = true, // plot the signal model (needs to be present)
			double signalMass = 125.,
			double theMuVal   = 1.,
			bool blinded  = false,  // blind the data in the plots?
			bool binned   = true,   // use binned data for all fits
			bool plotDerivative = false,
			bool dumpEventsOnly = false,
			bool verbose  = false  ) {
  
  
  // ======== ECAL LASER FILTER ===========
  std::map<int, std::map<int,std::vector<int>*>*> ecalLaserMap;
  std::string filterFileName = "ECALLaserCorrectionFilteredEvents.txt";
  createECALFilterMap(filterFileName, &ecalLaserMap);

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

  // =============================================================
  // load signal WS if required...
  TFile* sigWSfile = NULL;
  RooWorkspace* wsig = NULL;
  
  // ---------------------------------------------------------------------------------------------------------
  float      totallumi = -1.;

  double     massmax = -1.;
  double     massmin = -1.;

  TString    projectDir;
  //TString    datafilename;
  std::vector<std::string>    datafilenames;
  std::vector<TString> procnamelist;
  TString    modname;
  TString    treename;
  TString    wsPrefix;

  std::map<TString,TString> auxCatMap;
  std::map<TString,TCut>    anaCatMap;
  std::map<TString,TString> catdesc;
  std::map<TString,int>     polorder;
  
  std::map<TString,fstBernModel*>    catModelMap;

  TCut theBaseCut;

  bool    computeMVAvar;
  TString mvaWeightFile;
  TString mvaDefFile;

  double theCMenergy = -1.;

  bool readStatus = readFromConfigCard( configCard       ,
					totallumi           ,
					computeMVAvar       ,
					mvaWeightFile       ,
					mvaDefFile          ,
					projectDir          ,
					//datafilename       ,
					datafilenames       ,
					procnamelist,
					modname             ,
					treename            ,
					wsPrefix            ,
					massmax, massmin,
					auxCatMap           ,
					anaCatMap           ,
					catdesc,
					polorder,
					theBaseCut          ,
					theCMenergy
					);
  
  if( !readStatus ) {
    std::cerr<<" ERROR: Could not read from card > "<<configCard.Data()<<" <."<<std::endl;
    return;
  }

  // ======================================================================
  // extracting signal models from WS if requested...

  RooRealVar* hmass = NULL;

  if( dosignal ) {
    sigWSfile = TFile::Open(TString::Format("%s/model/ubersignalmodel.root", projectDir.Data()) ) ;  
    if( sigWSfile ) {
      wsig = (RooWorkspace*) sigWSfile->Get("wsig");
      hmass = wsig->var("CMS_hgg_mass");
      hmass->setBins(4*(int)(massmax-massmin));
    } else {
      std::cout<<" [WARNING] Could not find file with signal workspace <"<<TString::Format("%s/model/ubsersignal.root", projectDir.Data())<<"."<<std::endl;
      std::cout<<"           Will not overlay signal model in plots. Work continues..."<<std::endl;
      dosignal = false; 
    }
  }

  if ( !dosignal ) {
    hmass = new RooRealVar("CMS_hgg_mass","m_{#gamma#gamma}",massmin,massmax,"GeV");
    hmass->setRange(massmin,massmax);
    hmass->setBins(4*(int)(massmax-massmin));
  }
  
  std::cout<<" Checking for directory: "<<TString::Format("%s/databkg/",projectDir.Data())<<" ... "<<std::endl;
  
  if( gSystem->mkdir(TString::Format("%s/databkg/",projectDir.Data())) == 0 )
    std::cerr<<" WARNING: Could not find directory to "<<TString::Format("%s/databkg/",projectDir.Data()).Data()<<". Creating it."<<std::endl;
  
  hmass->setRange("fitrange",massmin,massmax);

  RooRealVar *sbweight = new RooRealVar("sbweight","",0.);
  
  int numFiles = datafilenames.size();

  TFile** datafiles = new TFile*[numFiles];
  TTree** hdata     = new TTree*[numFiles];

  for (int iFile = 0; iFile<numFiles; ++iFile) {
    datafiles[iFile] = TFile::Open(datafilenames[iFile].c_str(),"READ");
    if ( !datafiles[iFile] ) {
      std::cerr<<" ERROR: Could not open datafile with name > "<<datafilenames[iFile]<<" <."<<std::endl;
      return;
    }    

    TDirectory *datadir = (TDirectory*) datafiles[iFile]->FindObjectAny(modname);
    if ( !datadir ) {
      std::cerr<<" ERROR: Could not find directory > "<<modname.Data()<<" < in datafile with name > "<<datafilenames[iFile]<<" <."<<std::endl;
      return;
    }
    
    hdata[iFile] = (TTree*)datadir->Get(treename.Data());  
    if ( !hdata[iFile] ) {
      std::cerr<<" ERROR: Could not find Tree > "<<treename.Data()<<" < in directory > "<<modname.Data()<<" in datafile with name > "<<datafilenames[iFile]<<" <."<<std::endl;
      return;
    }
  }


//   TFile *datafile = new TFile(datafilename.Data(),"READ");
//   if ( !datafile ) {
//     std::cerr<<" ERROR: Could not open datafile with name > "<<datafilename.Data()<<" <."<<std::endl;
//     return;
//   }

//   TDirectory *datadir = (TDirectory*) datafile->FindObjectAny(modname);
//   if ( !datadir ) {
//     std::cerr<<" ERROR: Could not find directory > "<<modname.Data()<<" < in datafile with name > "<<datafilename.Data()<<" <."<<std::endl;
//     return;
//   }

//   TTree *hdata = (TTree*)datadir->Get(treename.Data());  
//   if ( !hdata ) {
//     std::cerr<<" ERROR: Could not find Tree > "<<treename.Data()<<" < in directory > "<<modname.Data()<<" in datafile with name > "<<datafilename.Data()<<" <."<<std::endl;
//     return;
//   }
  
  // open the MVA files, if requested
  TFile* tmvaOutput = NULL;
  TString weights   = "";
  std::vector<std::string> *varlist = NULL;
  TFile *friendtmp = NULL;

  if( computeMVAvar ) {
    friendtmp  = new TFile("friendtmp.root","RECREATE");
    tmvaOutput = new TFile(mvaDefFile.Data(),"READ");//root file to store training information  
    varlist = (std::vector<std::string>*)tmvaOutput->Get("varlist");
    weights = mvaWeightFile;
    for( int iFile=0; iFile<numFiles; ++iFile)
      ApplyAsFriend(hdata[iFile],weights,*varlist,"bdt");
  }    

  // ----------------------------------------------------------------------
  // some auxiliray vectro (don't know the meaning of all of them ... yet...
  std::vector<RooAbsData*> data_vec;
  std::vector<RooAbsData*> data_vecplot;

  std::vector<RooAbsData*> databinned_vec;

  std::vector<RooAbsPdf*>  pdfShape_vec;   // vector to store the NOT-EXTENDED PDFs (aka pdfshape)
  std::vector<RooAbsPdf*>  pdf_vec;        // vector to store the EXTENDED PDFs
  
  std::vector<RooAbsReal*> normu_vec;      // this holds the normalization vars for each Cat (needed in bands for combined cat)

  std::map<TString,int>    anaCatIndexMap;
  std::map<TString,int>    normVecIndexMap;

  RooArgList               normList;       // list of range-limityed normalizations (needed for error bands on combined category)

  // ----------------------------------------------------------------------
  // define output works
  RooWorkspace *wOut = new RooWorkspace("wbkg","wbkg") ;
  
  // util;ities for the combined fit
  RooCategory     finalcat  ("finalcat",  "finalcat") ;  
  
  for(std::map<TString,TCut>::iterator it = anaCatMap.begin(); it != anaCatMap.end(); ++it) {
    TString catname = it->first;
    finalcat.defineType(catname);
  }

  RooSimultaneous fullbkgpdf("fullbkgpdf","fullbkgpdf",finalcat);

  RooDataSet      datacomb      ("datacomb",      "datacomb",      RooArgList(*hmass,*sbweight,finalcat),"sbweight") ;
  RooDataSet      datacombplot  ("datacombplot",  "datacombplot",  RooArgList(*hmass,*sbweight,finalcat),"sbweight") ;

  // map to get the total data binned...
  std::map<std::string,RooDataHist*> binnedDataMap;
  
  RooDataSet *datacombcat     = new RooDataSet("data_combcat",""    ,RooArgList(*hmass,*sbweight),"sbweight") ;
  RooDataSet *datacombcatplot = new RooDataSet("data_combcatplot","",RooArgList(*hmass,*sbweight),"sbweight") ;
  
  // add the 'combcat' to the list...if more than one cat
  if( anaCatMap.size() > 1 ) {
    std::cout<<"  Inserting Combined Event Class into Analysis Class Map."<<std::endl;
    anaCatMap.insert(std::pair<TString,TCut>("combcat",TCut("")));
    catdesc.insert(std::pair<TString,TString>("combcat","Combined"));
  }
  
  int catCounter = -1;
  TString skipCatName = "";
  
  TCut masscut   =  theBaseCut && TCut(TString::Format("mass > %f && mass < %f",massmin,massmax).Data());

  // =================================================
  // stuff needed for the signal model
  std::map<TString, RooExtendPdf*> catSignalPdfMap;
  std::map<TString, double       > catSignalNeventsMap;
  RooConstVar* theLumi = new RooConstVar("theLumi","",totallumi*1000.);
  RooConstVar* theMu   = new RooConstVar("theMu","",theMuVal);  
  RooFormulaVar* snorm = new RooFormulaVar("snorm","","@0*@1",RooArgList(*theMu,*theLumi));
  RooRealVar* MH = NULL;
  // =================================================

  std::map<TString,int*> skippedEvents;
  for(std::map<TString,TCut>::iterator it = anaCatMap.begin(); it != anaCatMap.end(); ++it) {
    
    TString catname = it->first;
    
    skippedEvents.insert(std::pair<TString,int*>(catname, new int(0)));
    
    if ( !catname.CompareTo("combcat") ) continue; // NEED TO TAKE COMBCAT AT THE END !!!
    catCounter++;
    
    // =====================================================================
    // load signal PDF for this event class if required....
    if( dosignal ) {
      RooArgList compList;
      if ( !MH   ) MH   = wsig->var("MH");
      MH->setVal(signalMass);
      for( unsigned int iProc = 0; iProc < procnamelist.size(); ++iProc ){
	
	// get the normalization Variable
	RooFormulaVar* absNormVar = (RooFormulaVar*) wsig->function(TString::Format("hggpdfsmrel_%s_%s_norm",catname.Data(),procnamelist[iProc].Data()).Data());
	// get the PDF
	RooAbsPdf*    sigPdfRel = wsig->pdf(TString::Format("hggpdfsmrel_%s_%s",catname.Data(),procnamelist[iProc].Data()).Data());
	RooExtendPdf* sigPdf    = new RooExtendPdf(TString::Format("hggpdfsmrel_%s_%s_ext",catname.Data(),procnamelist[iProc].Data()).Data(),"",*sigPdfRel,*absNormVar);      
	compList.add(*sigPdf);
      }
      RooAddPdf* sigpdfcat       = new RooAddPdf(TString::Format("hggpdfsmrel_tmp_%s",catname.Data()).Data(),"",compList);
      RooExtendPdf* sigpdfcatExt = new RooExtendPdf(TString::Format("hggpdfsmrel_tmp_ext_%s",catname.Data()).Data(),"",*sigpdfcat,*snorm);
      catSignalPdfMap    .insert( std::pair<TString,RooExtendPdf*> ( catname, sigpdfcatExt ) );
      catSignalNeventsMap.insert( std::pair<TString,double       > ( catname, sigpdfcatExt->expectedEvents(*hmass) ) );
      
      //std::cout<<" NUMBER OF EXPECTED EVENT FOR CLASS "<<catname<<" : "<<sigpdfcatExt->expectedEvents(*hmass)<<std::endl;

    }
    // =====================================================================

    std::cout<<"  Creating dataset with name = "<<it->first<<std::endl;

    // check if we're in a sub-cat or the comb-cat

    
    RooDataSet* data     = makedset(TString("data_")+catname, hdata, numFiles, masscut && it->second, hmass, sbweight, 1., catCounter, dumpEventsOnly, &ecalLaserMap, skippedEvents.find(catname)->second);
    
    if( dumpEventsOnly ) continue;

    RooDataSet* dataplot = (RooDataSet*) data->reduce( "CMS_hgg_mass < 110. || CMS_hgg_mass > 150." );
    
    // append the data to the combined data...
    RooDataSet *datacat     = new RooDataSet(TString("datacat")+catname,    "",*hmass,Index(finalcat),Import(catname,*data)) ;
    RooDataSet *datacatplot = new RooDataSet(TString("datacatplot")+catname,"",*hmass,Index(finalcat),Import(catname,*dataplot)) ;

    if( !datacat ) {
      std::cout<<" ERROR: Could not load data "<<TString("datacat")+catname<<"."<<std::endl;
      return;
    }
    
    datacomb    .append(*datacat);
    datacombplot.append(*datacatplot);
    
    datacombcat    ->append(*data);
    datacombcatplot->append(*dataplot);
    
    // normalization for this category
    RooRealVar *nbkg = new RooRealVar(TString::Format("CMS_hgg_%s_bkgshape_norm",catname.Data()),"",(double) datacat->sumEntries(),0.0,250e3);
    
    // we keep track of the normalization vars only for N-1 cats, naming convetnions hystoric...
    if( ( catCounter < (int) anaCatMap.size() - 2 ) ) {
      RooRealVar* cbkg = new RooRealVar(TString::Format("cbkg%s",catname.Data()),"",0.0,0.0,1e3);
      cbkg->removeRange();

      normu_vec.push_back(cbkg);
      std::cout<<" *** PUSHING BACK in normu_vec: "<<catname<<" ... size of normu_vec is now: "<<normu_vec.size()<<std::endl;
      
      normVecIndexMap.insert(std::pair<TString,int>(catname,catCounter));
      normList.add(*cbkg);
    } else {
      skipCatName = catname;
    }
    
    /// generate the Bernstrin polynomial (FIX-ME: add possibility ro create other models...)
    fstBernModel* theBGmodel = new fstBernModel(hmass, polorder.find(catname)->second, catCounter, catname, "hgg");            // using my dedicated class...
    
    catModelMap.insert( std::pair<TString,fstBernModel*>(it->first,theBGmodel) );
    std::cout<<" Added BG model for catname "<<it->first<<std::endl;


    std::cout<<" model name is "<<theBGmodel->getPdf()->GetName()<<" with order "<<polorder.find(catname)->second<<"."<<std::endl;

    RooAbsPdf*    bkgshape   = theBGmodel->getPdf();                                              // the BG shape
    RooAbsPdf*    bkgpdf     = new RooExtendPdf(TString("bkgpdf")+catname,"",*bkgshape,*nbkg);    // the extended PDF
    
    // add the extedned PDF to the RooSimultaneous holding all models...
    fullbkgpdf.addPdf(*bkgpdf,catname);
    // store the NON-EXTENDED PDF for usgae to compute the error bands later..
    pdfShape_vec.push_back(bkgshape);
    pdf_vec     .push_back(bkgpdf);
    data_vec    .push_back(data);      
    //data_vecplot.push_back( data->reduce( TString::Format("CMS_hzg_mass_%s < 120. || CMS_hzg_mass_%s > 140.", (doElectrons ? "ee" : "mm"), (doElectrons ? "ee" : "mm")).Data() ) );
    data_vecplot.push_back( data->reduce( "CMS_hgg_mass < 110. || CMS_hgg_mass > 150." ) );

    // store the index for alter reference
    anaCatIndexMap.insert(std::pair<TString,int>(catname,catCounter));
    
    // generate the binned dataset (to be put into the workspace... just in case...)
    RooDataHist *databinned = new RooDataHist(TString("databinned_")+catname,"",*hmass,*data);
    databinned_vec    .push_back(databinned);

    binnedDataMap.insert( std::pair<std::string,RooDataHist*> ( catname.Data(), databinned ) );

    wOut->import(*data);
    wOut->import(*databinned);
  }

  if (dumpEventsOnly) return;

  // check if we have a combined set (alwways true if more than 1 Cat)
  std::map<TString,TCut>::iterator it = anaCatMap.find("combcat");
  if (it != anaCatMap.end()) {
    
    catCounter++;
    std::cout<<" *********** SHOULD NEVER BE HERE ************* "<<std::endl;
    TString catname = it->first;
    RooDataSet* data = datacombcat;   // we're looking at the last cat (by construction the combination)
    RooDataSet* dataplot = datacombcatplot ;
    
    data_vec.push_back( data );
    data_vecplot.push_back( dataplot );
    
    // sum up all the cts PDFs for combined PDF
    RooArgList subpdfs;
    for (unsigned int ipdf=0; ipdf<pdf_vec.size(); ++ipdf) {
      subpdfs.add(*pdf_vec.at(ipdf));
    }
    RooAddPdf* bkgpdf = new RooAddPdf(TString("bkgpdf")+catname,"",subpdfs);
    pdfShape_vec.push_back(bkgpdf);      
    pdf_vec     .push_back(bkgpdf);  // I don't think this is really needed though....

    anaCatIndexMap.insert(std::pair<TString,int>(catname,catCounter));

    // generate the binned dataset (to be put into the workspace... just in case...)
    RooDataHist *databinned = new RooDataHist(TString("databinned_")+catname,"",*hmass,*data);
    
    wOut->import(*data);
    wOut->import(*databinned);
  }
  
  // generate the combined binned dataset from the map
  RooDataHist* datacomb_binned = new RooDataHist("datacomb_binned","",*hmass, finalcat, binnedDataMap);
  
  RooFitResult *fullbkgfitres = NULL;
  
  if( binned ) {
    fullbkgpdf.fitTo(*datacomb_binned,Strategy(1),Minos(kFALSE),Save(kTRUE));
    fullbkgpdf.fitTo(*datacomb_binned,Strategy(1),Minos(kFALSE),Save(kTRUE));
    fullbkgfitres = fullbkgpdf.fitTo(*datacomb_binned,Strategy(2),Minos(kFALSE),Save(kTRUE));
  } else {
    fullbkgpdf.fitTo(datacomb,Strategy(1),Minos(kFALSE),Save(kTRUE));
    fullbkgpdf.fitTo(datacomb,Strategy(1),Minos(kFALSE),Save(kTRUE));
    fullbkgfitres = fullbkgpdf.fitTo(datacomb,Strategy(2),Minos(kFALSE),Save(kTRUE));
  }

  // in principle we're done now, so store the results in the output workspace
  wOut->import(datacomb);  
  wOut->import(*datacomb_binned);  
  wOut->import(fullbkgpdf);
  //wOut->import(*fullbkgfitres);

  if( verbose ) wOut->Print();

  wOut->writeToFile( TString::Format("%s/databkg/bkgdatawithfit.root", projectDir.Data()) ) ;  
  
  
  // --------------------------------------------------------------------------------------------
  // Now comesd the plotting
  // chage the Statistics style...
  gStyle->SetOptStat(1110);
  
  // we want to plot in 1GeV bins (apparently...)
  UInt_t nbins = (UInt_t) (massmax-massmin);
  
  // here we'll store the curves for the bands...
  std::vector<RooCurve*> fitcurves;
  fitcurves.resize(anaCatMap.size());

  // loop again over the cats
  TCanvas **canbkg = new TCanvas*[anaCatMap.size()];
  TCanvas **canbkg_der  = new TCanvas*[anaCatMap.size()];
  TCanvas **canbkg_der2 = new TCanvas*[anaCatMap.size()];

  RooPlot** plot   = new RooPlot*[anaCatMap.size()];
  //RooPlot** plot_der   = new RooPlot*[anaCatMap.size()];

  TLatex** lat  = new TLatex*[anaCatMap.size()];
  TLatex** lat2 = new TLatex*[anaCatMap.size()];

  TLatex** derlat  = new TLatex*[anaCatMap.size()];
  TLatex** derlat2 = new TLatex*[anaCatMap.size()];

  bool doCombCatLast = true;

  
  TLatex* topText = new TLatex(0.15,0.94,  "CMS preliminary");
  topText->SetNDC();
  topText->SetTextAlign(11);
  topText->SetTextFont(42);
  topText->SetTextSize(0.04);

  TLatex* topText2 = new TLatex(0.95,0.94,    TString::Format("#sqrt{s} = %.1f TeV, L = %.2f fb^{-1}",theCMenergy,totallumi).Data());
  topText2->SetNDC();
  topText2->SetTextAlign(31);
  topText2->SetTextFont(42);
  topText2->SetTextSize(0.04);

  
  for(std::map<TString,TCut>::iterator itNew = anaCatMap.begin(); itNew != anaCatMap.end(); ) {
    TString catname = itNew->first;

    std::cout<<catname<<std::endl;

    // combined last... so skip it if we hit it the first time
    if ( !catname.CompareTo("combcat") && doCombCatLast ) {
      itNew++;
      std::cout<<" Trying to skip combcat..."<<std::endl;
      if ( itNew == anaCatMap.end() )  { // 'combcat' is per chance the last one...
	doCombCatLast = false;
	itNew = anaCatMap.find("combcat");
	std::cout<<" ...  combcat is last by chance."<<std::endl;
      }
      continue;
    }

    int catIdx = anaCatIndexMap.find(catname)->second;
    
    // plot the data and the fit 
    canbkg[catIdx] = new TCanvas;

    plot  [catIdx] = hmass->frame(Bins(nbins),Range("fitrange"));

    // first plot the data invisibly... and put the fitted BG model on top...
    data_vec    .at(catIdx)->plotOn(plot[catIdx],RooFit::LineColor(kWhite),MarkerColor(kWhite),Invisible());

    if( dosignal && catname.CompareTo("combcat") ) {
      catSignalPdfMap.find(catname)->second->plotOn(plot[catIdx],RooFit::DrawOption("CFC"),RooFit::FillColor(18),RooFit::Normalization(catSignalNeventsMap.find(catname)->second, RooAbsReal::NumEvent));
      catSignalPdfMap.find(catname)->second->plotOn(plot[catIdx],RooFit::Normalization(catSignalNeventsMap.find(catname)->second, RooAbsReal::NumEvent));
    }

    pdfShape_vec.at(catIdx)->plotOn(plot[catIdx],RooFit::LineColor(kRed),Range("fitrange"),NormRange("fitrange"));
    if( blinded )
      data_vecplot.at(catIdx)->plotOn(plot[catIdx]);
    
    // if toggled on, plot also the Data visibly
    if( !blinded ) {
      data_vec.at(catIdx)->plotOn(plot[catIdx]);
    }

    
    // some cosmetics..
    plot[catIdx]->SetTitle("");      
    plot[catIdx]->SetMinimum(1e-5);
    plot[catIdx]->SetMaximum(1.40*plot[catIdx]->GetMaximum());
    plot[catIdx]->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
    plot[catIdx]->Draw();       


    // legend....
    int numEntries = 1;
    if(  dosignal &&  dobands ) numEntries = 4;
    if( !dosignal &&  dobands ) numEntries = 3;
    if(  dosignal && !dobands ) numEntries = 2;

    bool iscombcat       = ( !catname.CompareTo("combcat") );


    TLegend *legmc = new TLegend(0.68,0.80-(numEntries*0.06),0.97,0.90);
    legmc->AddEntry(plot[catIdx]->getObject( ( (dosignal && !iscombcat) ? 4 : 2 ) ),"Data","LPE");
    legmc->AddEntry(plot[catIdx]->getObject( ( (dosignal && !iscombcat) ? 3 : 1 ) ),"Bkg Model","L");
    
    // this part computes the 1/2-sigma bands.    
    TGraphAsymmErrors *onesigma = NULL;
    TGraphAsymmErrors *twosigma = NULL;
    
    TGraph* firstDer = NULL;
    TGraph* secondDer = NULL;
    
    RooAddition* sumcatsnm1 = NULL;
    
    if ( dobands ) { //&& icat == (catnames.size() - 1) ) {
      //if ( dobands && !iscombcat ) { //&& icat == (catnames.size() - 1) ) {

      onesigma = new TGraphAsymmErrors();
      twosigma = new TGraphAsymmErrors();
      
      // get the PDF for this cat from the vector
      RooAbsPdf *thisPdf = pdfShape_vec.at(catIdx); 

      // get the nominal fir curve
      RooCurve *nomcurve = dynamic_cast<RooCurve*>(plot[catIdx]->getObject( (dosignal && !iscombcat) ? 3 : 1));
      fitcurves[catIdx] = nomcurve;

      //RooAbsData *datanorm = ( iscombcat ? &datacomb : data_vec.at(catIdx) );
      RooAbsData *datanorm = ( iscombcat ? (binned ? (RooAbsData*) datacomb_binned : (RooAbsData*) &datacomb) : (binned ? databinned_vec.at(catIdx) : data_vec.at(catIdx) ));

      // this si the nornmalization in the 'sliding-window' (i.e. per 'test-bin')
      RooRealVar *nlim = new RooRealVar(TString::Format("nlim%s",catname.Data()),"",0.0,0.0,10.0);
      nlim->removeRange();
      
      if( iscombcat ) {
	// ----------- HISTORIC NAMING  ----------------------------------------
	sumcatsnm1 = new RooAddition("sumcatsnm1","",normList);   // summing all normalizations epect the last Cat
	// this is the normlization of the last Cat
	RooFormulaVar *nlast = new RooFormulaVar("nlast","","TMath::Max(0.1,@0-@1)",RooArgList(*nlim,*sumcatsnm1));
	// ... and adding it ot the list of norms
	normu_vec.push_back(nlast);
	normVecIndexMap.insert(std::pair<TString,int>(skipCatName,catCounter-1));
      }
      
      for (int i=1; i<(plot[catIdx]->GetXaxis()->GetNbins()+1); ++i) {
	
	// this defines the 'binning' we use for the error bands
	double lowedge = plot[catIdx]->GetXaxis()->GetBinLowEdge(i);
	double upedge = plot[catIdx]->GetXaxis()->GetBinUpEdge(i);
	double center = plot[catIdx]->GetXaxis()->GetBinCenter(i);
	
	// get the nominal value at the center of the bin
	double nombkg = nomcurve->interpolate(center);
	nlim->setVal(nombkg);
	hmass->setRange("errRange",lowedge,upedge);
	
	// this is the new extended PDF whith the normalization restricted to the bin-area
	RooAbsPdf *extLimPdf = NULL;
	if( iscombcat ) {
	  extLimPdf = new RooSimultaneous("epdf","",finalcat);
	  // loop over the cats and generate temporary extended PDFs
	  
	  for(std::map<TString,TCut>::iterator it2 = anaCatMap.begin(); it2 != anaCatMap.end(); ++it2) {
	    if ( !(it2->first).CompareTo("combcat") ) continue;
	    int catIdxNew  = anaCatIndexMap .find( it2->first )->second;

	    int normIdx    = normVecIndexMap.find( it2->first )->second;

	    std::cout<<" combining classes .. "<<it2->first<<"   with idx = "<<catIdxNew<<"   and normIdx = "<<normIdx<<std::endl;
	    

	    RooRealVar *rvar = dynamic_cast<RooRealVar*>(normu_vec.at( normIdx ));
	    if (rvar) rvar->setVal(fitcurves.at( catIdxNew )->interpolate(center));
	    RooExtendPdf *ecpdf = new RooExtendPdf(TString::Format("ecpdf%s",(it2->first).Data()),"",*pdfShape_vec.at( catIdxNew ),*normu_vec.at( normIdx ),"errRange");
	    static_cast<RooSimultaneous*>(extLimPdf)->addPdf(*ecpdf,(it2->first));
	  }
	} else
	  extLimPdf = new RooExtendPdf("extLimPdf","",*thisPdf,*nlim,"errRange");
	
	RooAbsReal *nll = extLimPdf->createNLL(*datanorm,Extended(),NumCPU(1));
	RooMinimizer minim(*nll);
	minim.setStrategy(0);
	//double clone = 1.0 - 2.0*RooStats::SignificanceToPValue(1.0);
	double cltwo = 1.0 - 2.0*RooStats::SignificanceToPValue(2.0);
	
	if (iscombcat) minim.setStrategy(2);
	
	minim.migrad();
	
	if (!iscombcat) { 
	  minim.minos(*nlim);
	}
	else {
	  minim.hesse();
	  nlim->removeAsymError();
	}
	
	if( verbose ) 
	  printf("errlo = %5f, errhi = %5f\n",nlim->getErrorLo(),nlim->getErrorHi());
	
	onesigma->SetPoint(i-1,center,nombkg);
	onesigma->SetPointError(i-1,0.,0.,-nlim->getErrorLo(),nlim->getErrorHi());
	
	// to get the 2-sigma bands...
	minim.setErrorLevel(0.5*pow(ROOT::Math::normal_quantile(1-0.5*(1-cltwo),1.0), 2)); // the 0.5 is because qmu is -2*NLL
	// eventually if cl = 0.95 this is the usual 1.92!      
	
	if (!iscombcat) { 
	  minim.migrad();
	  minim.minos(*nlim);
	}
	else {
	  nlim->setError(2.0*nlim->getError());
	  nlim->removeAsymError();          
	}
	
	twosigma->SetPoint(i-1,center,nombkg);
	twosigma->SetPointError(i-1,0.,0.,-nlim->getErrorLo(),nlim->getErrorHi());      
	
	// for memory clean-up
	delete nll;
	delete extLimPdf;
      }
      
      hmass->setRange("errRange",massmin,massmax);
      
      if( verbose )
	onesigma->Print("V");
      
      // plot[catIdx] the error bands
      twosigma->SetLineColor(kGreen);
      twosigma->SetFillColor(kGreen);
      twosigma->SetMarkerColor(kGreen);
      twosigma->Draw("L3 SAME");     
      
      onesigma->SetLineColor(kYellow);
      onesigma->SetFillColor(kYellow);
      onesigma->SetMarkerColor(kYellow);
      onesigma->Draw("L3 SAME");
      
      plot[catIdx]->Draw("SAME");
      
      // and add the error bands to the legend
      legmc->AddEntry(onesigma,"#pm1 #sigma","F");  
      legmc->AddEntry(twosigma,"#pm2 #sigma","F");  
    }
    
    if( dosignal && !iscombcat )
      legmc->AddEntry(plot[catIdx]->getObject( 1 ),"Signal Model","F");

    // rest of the legend ....
    legmc->SetBorderSize(0);
    legmc->SetFillStyle(0);
    legmc->Draw();   

    lat[catIdx]  = new TLatex(massmin+3.,0.9*plot[catIdx]->GetMaximum(),TString::Format("#scale[0.7]{#splitline{CMS preliminary}{#sqrt{s} = %.1f TeV L = %.2f fb^{-1}}}",theCMenergy,totallumi));
    lat2[catIdx] = new TLatex(massmin+3.,0.75*plot[catIdx]->GetMaximum(),catdesc.find(catname)->second);
    
    lat[catIdx] ->Draw();
    lat2[catIdx]->Draw();
    
    // -------------------------------------------------------    
    // save canvas in different formats
    canbkg[catIdx]->SaveAs( TString::Format("%s/databkg/databkg%s.pdf",projectDir.Data(),catname.Data()) );
//     canbkg[catIdx]->SaveAs(TString("databkg") + catname + TString(".eps"));
//     canbkg[catIdx]->SaveAs(TString("databkg") + catname + TString(".gif"));
//     canbkg[catIdx]->SaveAs(TString("databkg") + catname + TString(".root"));              
    
    if (plotDerivative) {


      derlat [catIdx] = new TLatex(0.1,0.95,TString::Format("#scale[0.7]{#splitline{CMS preliminary}{#sqrt{s} = %.1f TeV L = %.2f fb^{-1}}}",theCMenergy,totallumi));
      derlat [catIdx]->SetNDC();
      
      derlat2[catIdx] = new TLatex(0.25,0.85,catdesc.find(catname)->second);
      derlat2[catIdx]->SetNDC();

      
      firstDer = new TGraph();

      for (int i=1; i<(plot[catIdx]->GetXaxis()->GetNbins()+1); ++i) {
	
	// this defines the 'binning' we use for the error bands
	double lowedge = plot[catIdx]->GetXaxis()->GetBinLowEdge(i);
	double upedge = plot[catIdx]->GetXaxis()->GetBinUpEdge(i);
	double center = plot[catIdx]->GetXaxis()->GetBinCenter(i);
	
	RooCurve *nomcurve = dynamic_cast<RooCurve*>(plot[catIdx]->getObject( (dosignal && !iscombcat) ? 3 : 1));

	// get the nominal value at the center of the bin
	double deri = (nomcurve->interpolate(upedge) - nomcurve->interpolate(lowedge))/(upedge-lowedge); 
	firstDer->SetPoint(i-1,center,deri);
      }

      // plot the data and the fit 

      canbkg_der[catIdx] = new TCanvas;
      canbkg_der[catIdx]->cd();

      firstDer->SetLineColor(kBlue);
      firstDer->SetLineWidth(2);
      firstDer->Draw("AL");

      firstDer->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
      firstDer->GetYaxis()->SetTitle("d^{2}N/dm_{#gamma#gamma}^{2} [1/GeV^{2}]");

      //derlat [catIdx]->Draw();
      topText->Draw();
      topText2->Draw();
      derlat2[catIdx]->Draw();
      
      //canbkg_der[catIdx]->SaveAs(TString("databkg_der_") + catname + TString(".pdf"));
      canbkg_der[catIdx]->SaveAs(TString::Format("%s/databkg/databkg_der_%s.pdf",projectDir.Data(),catname.Data()));

      secondDer = new TGraph();

      for (int i=1; i<(plot[catIdx]->GetXaxis()->GetNbins()+1); ++i) {
	
	// this defines the 'binning' we use for the error bands
	double lowedge = plot[catIdx]->GetXaxis()->GetBinLowEdge(i);
	double upedge = plot[catIdx]->GetXaxis()->GetBinUpEdge(i);
	double center = plot[catIdx]->GetXaxis()->GetBinCenter(i);
	
	// get the nominal value at the center of the bin
	double deri = (firstDer->Eval(upedge) - firstDer->Eval(lowedge))/(upedge-lowedge); 
	secondDer->SetPoint(i-1,center,deri);
      }

      // plot the data and the fit 
      canbkg_der2[catIdx] = new TCanvas;

      secondDer->SetLineColor(kRed);
      secondDer->SetLineWidth(2);
      secondDer->Draw("AL");

      secondDer->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
      secondDer->GetYaxis()->SetTitle("d^{3}N/dm_{#gamma#gamma}^{3} [1/GeV^{3}]");

      topText->Draw();
      topText2->Draw();
      derlat2[catIdx]->Draw();

      //canbkg_der2[catIdx]->SaveAs(TString("databkg_der2_") + catname + TString(".pdf"));
      canbkg_der2[catIdx]->SaveAs(TString::Format("%s/databkg/databkg_der2_%s.pdf",projectDir.Data(),catname.Data()));

    }
    
    if( !doCombCatLast ) break;

    itNew++;    
    if (itNew == anaCatMap.end() && doCombCatLast ) {
      itNew = anaCatMap.find("combcat");  // do this at the end...
      doCombCatLast = false;
    }
  }
  
  
  if( verbose ) {
    printf("IntLumi = %5f\n",totallumi);
    printf("ndata:\n");
    for(std::map<TString,TCut>::iterator itNew = anaCatMap.begin(); itNew != anaCatMap.end(); ++itNew) {
      TString catname = itNew->first;
      printf("%s N = %i  (skipped events = %d)\n",catname.Data(),data_vec.at( anaCatIndexMap.find(catname)->second )->numEntries(), *(skippedEvents.find(catname)->second));      
    }   
    printf("\n");
  } 
  

  return;
  
}


// ----------------------------------------------------------------------------------------
bool readFromConfigCard( TString fileName,
			 float& totallumi,
			 bool&    computeMVAvar,
			 TString& mvaWeightFile,
			 TString& mvaDefFile,
			 TString& projectDir,
			 //TString& datafilename,
			 std::vector<std::string>& datafilenames,
			 std::vector<TString>& procNameList,
			 TString& modname,
			 TString& treename,
			 TString& wsPrefix,
			 double& massmax, double& massmin,
			 std::map<TString,TString>& auxCats,
			 std::map<TString,TCut>   & anaCats,
			 std::map<TString,TString>& catDesc,
			 std::map<TString,int>    & polOrders,
			 TCut& baseCut,
			 double& theCMenergy
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

    char name[1000];
    char smearcat[1000];
    int catIdx = -1;
    char catName[400];
    char theCat[400];
    float lumi    = -1.;
    float massval = -1.;
    float smear   = -1.;
    int polOrder  = 0;
    char onoff[3];

    if ( sscanf(line,"PROJECTDIR %s", name ) ) projectDir = TString(name);
    else if  ( sscanf(line,"LUMI %f", &lumi) ) totallumi = lumi;
    //else if  ( sscanf(line,"DATAFILE %s",&name) ) datafilename = TString(name);
    else if  ( sscanf(line,"DATAFILE %s", name) ) datafilenames.push_back(name);
    else if  ( sscanf(line,"MODNAME %s", name) ) modname = TString(name);
    else if  ( sscanf(line,"TREENAME %s", name) ) treename = TString(name);
    else if  ( sscanf(line,"WS1PREFIX %s", name) ) wsPrefix = TString(name);
    else if( sscanf(line,"MINMSS %f",&massval) ) massmin = (double) massval;
    else if( sscanf(line,"MAXMSS %f",&massval) ) massmax = (double) massval;
    else if( sscanf(line,"CMENERGY %f",&massval) ) theCMenergy = (double) massval;
    else if  ( sscanf(line,"AUXCAT %d %s",&catIdx, catName) ) {
      // parsing the ctegory definition like:
      // AUXCAT	0	masscut		" mass>100.0 && mass<180. "
      std::string totLine = line;
      int startCat = totLine.find_first_of("\"");
      int endCat   = totLine.find_last_of("\"");
      std::string catLine = totLine.substr(startCat+1, endCat-startCat-1);
      auxCats.insert(std::pair<TString,TString>(TString(catName),TString(catLine.c_str())));
    }
    else if  ( sscanf(line,"COMPUTEMVA ON %s",  name) ) { mvaWeightFile = TString(name); computeMVAvar = true; }
    else if  ( sscanf(line,"MVADEFFILE %s",  name) ) { mvaDefFile = TString(name);}
    else if  ( sscanf(line,"ANACAT %d %s %f %s Bern/%d",&catIdx, name, &smear, smearcat, &polOrder) ) {
      // ANACAT	0	hzgcat0		0.005432	Bern/5		" basecut "					StDesc(H#rightarrow Zg#rightarrow 2mu#gamma)
      
      polOrders.insert(std::pair<TString,int>(TString(name),polOrder));
      
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

	fPos = catLine.find_first_of(" ()|&!");
	std::string label;
	if( fPos != std::string::npos ) {
	  label = catLine.substr(0,fPos);
	  catLine.erase(0,fPos);
	  // and erase starting entries...
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
      anaCats.insert(std::pair<TString,TCut>(TString(name),TCut(theCutLine.c_str())));

      std::string theLine = line;
      size_t fPos2 = theLine.rfind("Desc(");
      if( fPos2 == std::string::npos ) {
	std::cerr<<" ERROR: ANACAT "<<catIdx<<" with name "<<name<<" has no string descriptor (StDexc(...))."<<std::endl;
	return false;
      }
      theLine.erase(0,fPos2+5);
      fPos2 = theLine.find_last_of(")");
      if( fPos2 == std::string::npos ) {
	std::cerr<<" ERROR: ANACAT "<<catIdx<<" with name "<<name<<" has no valid string descriptor (StDexc(...))."<<std::endl;
	return false;
      }
      theLine.erase(fPos2);
      catDesc.insert(std::pair<TString,TString>(TString(name),TString(theLine)));
    } else if  ( sscanf(line,"BASECAT %s", name) ) {
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
    } else {

      // additional setup
      int theProc = -1;
      char name[30];
      
      char onoff[3];
      char MCfileName[200];
      
      
      if ( sscanf(line,"PROC %d %s %s file:%s",&theProc, name, onoff, MCfileName) ) {
	//# ------------------------------------------------------------------------------------------------------------------------------------------------
	//PROC    0	ggh	OFF			file:/scratch/fabstoec/cms/hist/hgg-7TeV-janReReco/merged/hgg-7TeV-janReReco_f11--h%dgg-gf-v14b-pu_noskim.root
	procNameList.push_back( TString(name) );
      }
    }
  }
  
  std::cout<<" done"<<std::endl;
  
  fclose(configFile);
  return true;
}
