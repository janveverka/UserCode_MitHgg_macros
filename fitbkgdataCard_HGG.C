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

RooDataSet *makedset(TString name, TTree *tree, TCut cut, RooRealVar *hmass, RooRealVar* sbweight = NULL, Double_t tWeight = 1.,
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

  tree->SetEstimate(tree->GetEntries());

  Int_t nev = tree->Draw("mass:run:lumi:evt",cut,"goff");    
  std::cout<<"  Nev = "<<nev<<std::endl;

  double *vals  = tree->GetV1();
  double *runs  = tree->GetV2();
  double *lumis = tree->GetV3();
  double *evts  = tree->GetV4();
  //double *weights = tree->GetW();  // this is Data, has no weights.
  
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
      
    if (noskip)
      dset->add(RooArgSet(*vset),tWeight);
  }
  
  return dset;
  
}

TTree* ApplyVBFMVAAsFriend(TTree *intree, std::string targetname, int mass) {

  
  TMVA::Reader* fReader = new TMVA::Reader( "!Color:!Silent:Error" );       
  
  TString Weights = TString("/home/fabstoec/cms/cmssw/029/CMSSW_5_3_2_patch4/src/MitPhysics/data/TMVA_vbf_6var_mjj100_diphopt_phopt_BDTG.weights.xml");
  
  float _jet1pt, _jet2pt, _deltajeteta, _dijetmass, _zeppenfeld, _dphidijetgg, _diphoptOverdiphomass, _pho1ptOverdiphomass, _pho2ptOverdiphomass;
  
  float _jet1eta, _jet2eta, _ph1pt, _ph2pt, _mass, _ptgg;
  

  // input variables
  intree->SetBranchAddress("jet1pt", &_jet1pt);
  intree->SetBranchAddress("jet2pt", &_jet2pt);
  intree->SetBranchAddress("jet1eta", &_jet1eta);
  intree->SetBranchAddress("jet2eta", &_jet2eta);
  intree->SetBranchAddress("ph1.pt", &_ph1pt);
  intree->SetBranchAddress("ph2.pt", &_ph2pt);
  intree->SetBranchAddress("dijetmass", &_dijetmass);
  intree->SetBranchAddress("dphidijetgg", &_dphidijetgg);
  intree->SetBranchAddress("zeppenfeld", &_zeppenfeld);
  intree->SetBranchAddress("mass", &_mass);
  intree->SetBranchAddress("ptgg", &_ptgg);

  // TMVA input variables
  fReader->AddVariable("jet1pt",&_jet1pt);
  fReader->AddVariable("jet2pt",&_jet2pt);
  fReader->AddVariable("abs(jet1eta-jet2eta)",&_deltajeteta);
  fReader->AddVariable("mj1j2",&_dijetmass);
  fReader->AddVariable("zepp",&_zeppenfeld);
  fReader->AddVariable("dphi",&_dphidijetgg);
  fReader->AddVariable("diphopt/diphoM",&_diphoptOverdiphomass);
  fReader->AddVariable("pho1pt/diphoM",&_pho1ptOverdiphomass);
  fReader->AddVariable("pho2pt/diphoM",&_pho2ptOverdiphomass);

  fReader->BookMVA("BDT method",Weights);

  assert(fReader);
  
  Float_t target = 0.;
  
  //initialize new friend tree
  //TTree *friendtree = new TTree(TString::Format("mvatree_%s_%d",targetname.c_str(),mass).Data(),"");
  TTree *friendtree = new TTree();
  friendtree->SetName(TString::Format("mvatree_%s_%d",targetname.c_str(),mass).Data());
  friendtree->Branch(targetname.c_str(),&target,TString::Format("%s/F",targetname.c_str()));
  
  int currenttree = -1;
  for (Long64_t iev=0; iev<intree->GetEntries(); ++iev) {
    if (iev%100000==0) printf("%i\n",int(iev));
    intree->GetEntry(iev);
    
    target = -99.;
    // assign varibles
    _deltajeteta = TMath::Abs(_jet1eta-_jet2eta);
    _diphoptOverdiphomass = _ptgg/_mass;
    _pho1ptOverdiphomass  = _ph1pt/_mass;
    _pho2ptOverdiphomass  = _ph2pt/_mass;

    if( (_pho1ptOverdiphomass > 40/120) && (_pho2ptOverdiphomass > 30/120) && (_jet1pt > 30) && (_jet2pt > 20) && (_dijetmass > 250) ) {
      target = fReader->EvaluateMVA("BDT method");
    }
    
    friendtree->Fill();
    
  }
  
  //clear TMVA reader
  delete fReader;
    
  intree->AddFriend(friendtree);
  return friendtree;
  
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
			 TString& datafilename,
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

void fitbkgdataCard_HGG(bool dobands  = true, 
#ifdef DOBDT
			TString configCard="templateHGG_8TeV.config", 
#else
			TString configCard="templateHGG_8TeV_CiC.config", 
#endif
			bool dosignal = false, // plot the signal model (needs to be present)
			bool blinded  = false,  // blind the data in the plots?
			bool verbose  = true  ) {
  
  
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
  std::map<TString,TCut>    anaCatMap;
  std::map<TString,TString> catdesc;
  std::map<TString,int>     polorder;
  
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
					datafilename        ,
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
    
  if( !gSystem->cd(TString::Format("%s/databkg/",projectDir.Data())) ) {
    std::cerr<<" ERROR: Could not change directory to "<<TString::Format("%s/databkg/",projectDir.Data()).Data()<<"."<<std::endl;
    return;
  }
  
  RooRealVar* hmass = new RooRealVar("CMS_hgg_mass","m_{#gamma#gamma}",massmin,massmax,"GeV");
  hmass->setRange(massmin,massmax);
  hmass->setBins(4*(int)(massmax-massmin));
  hmass->setRange("fitrange",massmin,massmax);

  RooRealVar *sbweight = new RooRealVar("sbweight","",0.);
  
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
    ApplyAsFriend(hdata,weights,*varlist,"bdt");
//     ApplyVBFMVAAsFriend(hdata,"vbfmva",125);

  }  



  // ----------------------------------------------------------------------
  // some auxiliray vectro (don't know the meaning of all of them ... yet...
  std::vector<RooAbsData*> data_vec;
  std::vector<RooAbsData*> data_vecplot;
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

  RooDataSet *datacombcat     = new RooDataSet("data_combcat",""    ,RooArgList(*hmass,*sbweight),"sbweight") ;
  RooDataSet *datacombcatplot = new RooDataSet("data_combcatplot","",RooArgList(*hmass,*sbweight),"sbweight") ;
  
  // add the 'combcat' to the list...if more than one cat
  if( anaCatMap.size() > 1 ) {
    anaCatMap.insert(std::pair<TString,TCut>("combcat",TCut("")));
    catdesc.insert(std::pair<TString,TString>("combcat","Combined"));
  }

  int catCounter = -1;
  TString skipCatName = "";

  TCut masscut   =  theBaseCut && TCut(TString::Format("mass > %f && mass < %f",massmin,massmax).Data());

  std::map<TString,int*> skippedEvents;
  for(std::map<TString,TCut>::iterator it = anaCatMap.begin(); it != anaCatMap.end(); ++it) {

    TString catname = it->first;

    skippedEvents.insert(std::pair<TString,int*>(catname, new int(0)));

    if ( !catname.CompareTo("combcat") ) continue; // NEED TO TAKE COMBCAT AT THE END !!!
    catCounter++;

    std::cout<<"  Creating dataset with name = "<<it->first<<std::endl;

    // check if we're in a sub-cat or the comb-cat

    //RooDataSet* data     = makedset(TString("data_")+catname, hdata, masscut && it->second, hmass, sbweight, 1., &ecalLaserMap, skippedEvents.find(catname)->second);
    RooDataSet* data     = makedset(TString("data_")+catname, hdata, masscut && it->second, hmass, sbweight, 1.);
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
    
    // we keep track of the normalizario vars only for N-1 cats, naming convetnions hystoric...
    if( ( catCounter < (int) anaCatMap.size() - 2 ) ) {
      RooRealVar* cbkg = new RooRealVar(TString::Format("cbkg%s",catname.Data()),"",0.0,0.0,1e3);
      cbkg->removeRange();
      normu_vec.push_back(cbkg);
      normVecIndexMap.insert(std::pair<TString,int>(catname,catCounter));
      normList.add(*cbkg);
    } else {
      skipCatName = catname;
    }
    
    /// generate the Bernstrin polynomial (FIX-ME: add possibility ro create other models...)
    fstBernModel* theBGmodel = new fstBernModel(hmass, polorder.find(catname)->second, catCounter, catname, "hgg");            // using my dedicated class...
    
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
    
    wOut->import(*data);
    wOut->import(*databinned);
  }

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

  // fit the RooSimultaneous to the combined dataset -> (we could also fit each cat separately)
  //RooDataHist* datacombHist = dataComb->binnedClone();
  fullbkgpdf.fitTo(datacomb,Strategy(1),Minos(kFALSE),Save(kTRUE));
  fullbkgpdf.fitTo(datacomb,Strategy(1),Minos(kFALSE),Save(kTRUE));
  RooFitResult *fullbkgfitres = fullbkgpdf.fitTo(datacomb,Strategy(2),Minos(kFALSE),Save(kTRUE));

  // in principle we're done now, so store the results in the output workspace
  wOut->import(datacomb);  
  wOut->import(fullbkgpdf);
  wOut->import(*fullbkgfitres);

  if( verbose ) wOut->Print();

  wOut->writeToFile("bkgdatawithfit.root") ;  
  
  
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
  RooPlot** plot   = new RooPlot*[anaCatMap.size()];

  TLatex** lat  = new TLatex*[anaCatMap.size()];
  TLatex** lat2 = new TLatex*[anaCatMap.size()];

  bool doCombCatLast = true;

  for(std::map<TString,TCut>::iterator it = anaCatMap.begin(); it != anaCatMap.end(); ) {
    TString catname = it->first;

    // combined last... so skip it if we hit it the first time
    if ( !catname.CompareTo("combcat") && doCombCatLast ) {
      it++;
      if ( it == anaCatMap.end() )  { // 'combcat' is per chance the last one...
	doCombCatLast = false;
	it = anaCatMap.find("combcat");
      }
      continue;
    }
    
    int catIdx = anaCatIndexMap.find(catname)->second;
    
    // plot the data and the fit 
    canbkg[catIdx] = new TCanvas;
    plot  [catIdx] = hmass->frame(Bins(nbins),Range("fitrange"));

    // first plot the data invisibly... and put the fitted BG model on top...
    data_vec    .at(catIdx)->plotOn(plot[catIdx],RooFit::LineColor(kWhite),MarkerColor(kWhite),Invisible());
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
    plot[catIdx]->GetXaxis()->SetTitle("m_{#gamma#gamma} (GeV/c^{2})");
    plot[catIdx]->Draw();       


    // legend....
    TLegend *legmc = new TLegend(0.68,0.70,0.97,0.90);
    legmc->AddEntry(plot[catIdx]->getObject(2),"Data","LPE");
    legmc->AddEntry(plot[catIdx]->getObject(1),"Bkg Model","L");
    
    // this part computes the 1/2-sigma bands.    
    TGraphAsymmErrors *onesigma = NULL;
    TGraphAsymmErrors *twosigma = NULL;


    RooAddition* sumcatsnm1 = NULL;

    bool iscombcat       = ( !catname.CompareTo("combcat") );
    if ( dobands && !iscombcat ) { //&& icat == (catnames.size() - 1) ) {
    //if ( dobands ) { //&& icat == (catnames.size() - 1) ) {

      onesigma = new TGraphAsymmErrors();
      twosigma = new TGraphAsymmErrors();

      // get the PDF for this cat from the vector
      RooAbsPdf *thisPdf = pdfShape_vec.at(catIdx); 

      // get the nominal fir curve
      RooCurve *nomcurve = dynamic_cast<RooCurve*>(plot[catIdx]->getObject(1));
      fitcurves[catIdx] = nomcurve;

      RooAbsData *datanorm = ( iscombcat ? &datacomb : data_vec.at(catIdx) );

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
	normVecIndexMap.insert(std::pair<TString,int>(skipCatName,catCounter));
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
	    int catIdx  = anaCatIndexMap.find( it2->first )->second;
	    int normIdx = normVecIndexMap.find( it2->first )->second;
	    RooRealVar *rvar = dynamic_cast<RooRealVar*>(normu_vec.at( normIdx ));
	    if (rvar) rvar->setVal(fitcurves.at( catIdx )->interpolate(center));
	    RooExtendPdf *ecpdf = new RooExtendPdf(TString::Format("ecpdf%s",(it2->first).Data()),"",*pdfShape_vec.at( catIdx ),*normu_vec.at( normIdx ),"errRange");
	    static_cast<RooSimultaneous*>(extLimPdf)->addPdf(*ecpdf,(it2->first));
	  }
	} else
	  extLimPdf = new RooExtendPdf("extLimPdf","",*thisPdf,*nlim,"errRange");
	
	RooAbsReal *nll = extLimPdf->createNLL(*datanorm,Extended(),NumCPU(1));
	RooMinimizer minim(*nll);
	minim.setStrategy(0);
	double clone = 1.0 - 2.0*RooStats::SignificanceToPValue(1.0);
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
    canbkg[catIdx]->SaveAs(TString("databkg") + catname + TString(".pdf"));
    canbkg[catIdx]->SaveAs(TString("databkg") + catname + TString(".eps"));
    canbkg[catIdx]->SaveAs(TString("databkg") + catname + TString(".gif"));
    canbkg[catIdx]->SaveAs(TString("databkg") + catname + TString(".root"));              

    if( !doCombCatLast ) break;

    it++;
    if (it == anaCatMap.end() && doCombCatLast ) {
      it = anaCatMap.find("combcat");  // do this at the end...
      doCombCatLast = false;
    }
  }
  

  if( verbose ) {
    printf("IntLumi = %5f\n",totallumi);
    printf("ndata:\n");
    for(std::map<TString,TCut>::iterator it = anaCatMap.begin(); it != anaCatMap.end(); ++it) {
      TString catname = it->first;
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
			 TString& datafilename,
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

    if ( sscanf(line,"PROJECTDIR %s",&name ) ) projectDir = TString(name);
    else if  ( sscanf(line,"LUMI %f",&lumi) ) totallumi = lumi;
    else if  ( sscanf(line,"DATAFILE %s",&name) ) datafilename = TString(name);
    else if  ( sscanf(line,"MODNAME %s",&name) ) modname = TString(name);
    else if  ( sscanf(line,"TREENAME %s",&name) ) treename = TString(name);
    else if  ( sscanf(line,"WS1PREFIX %s",&name) ) wsPrefix = TString(name);
    else if( sscanf(line,"MINMSS %f",&massval) ) massmin = (double) massval;
    else if( sscanf(line,"MAXMSS %f",&massval) ) massmax = (double) massval;
    else if( sscanf(line,"CMENERGY %f",&massval) ) theCMenergy = (double) massval;
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
    else if  ( sscanf(line,"MVADEFFILE %s",&name) ) { mvaDefFile = TString(name);}
    else if  ( sscanf(line,"ANACAT %d %s %f %s Bern/%d",&catIdx, &name, &smear, &smearcat, &polOrder) ) {
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
    } else if  ( sscanf(line,"BASECAT %s",&name) ) {
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
