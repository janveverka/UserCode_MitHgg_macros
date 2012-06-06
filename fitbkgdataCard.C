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

using namespace RooFit;

#include "FstModels.h"     // some helpfer functions for different models (well... I cut away all but Bernstein, but can be brought back...)

// ----------------------------------------------------------------------------------------------------
bool readFromConfigCard(TString cardName,
			TString& projectDir,
			std::vector<TString>& catNames,
			std::vector<TString>& catDesc ,		
			std::vector<int>    & polOrder,
			double& massmin,
			double& massmax,
			double& theCMenergy
			);


// ----------------------------------------------------------------------------------------------------

void fitbkgdataCard(TString configCard="template.config", 
		    bool dobands  = true,  // create baerror bands for BG models
		    bool dosignal = false, // plot the signal model (needs to be present)
		    bool blinded  = true,  // blind the data in the plots?
		    bool verbose  = true  ) {
  
  gROOT->Macro("MitStyle.C");
  gStyle->SetErrorX(0); 
  gStyle->SetOptStat(0);
  gROOT->ForceStyle();  
  
  TString projectDir;

  std::vector<TString> catdesc;
  std::vector<TString> catnames;  
  std::vector<int>     polorder;

  double massmin = -1.;
  double massmax = -1.;

  double theCMenergy = -1.;

  bool readStatus = readFromConfigCard( configCard,
					projectDir,
					catnames,
					catdesc,
					polorder,
					massmin,
					massmax,
					theCMenergy
					);
  
  if( !readStatus ) {
    std::cerr<<" ERROR: Could not read from card > "<<configCard.Data()<<" <."<<std::endl;
    return;
  }
  
  TFile *fdata = new TFile(TString::Format("%s/CMS-HGG-data.root",projectDir.Data()),"READ");
  if( !fdata ) {
    std::cerr<<" ERROR: Could not open file "<<projectDir.Data()<<"/CMS-HGG-data.root."<<std::endl;
    return;
  }
  
  if( !gSystem->cd(TString::Format("%s/databkg/",projectDir.Data())) ) {
    std::cerr<<" ERROR: Could not change directory to "<<TString::Format("%s/databkg/",projectDir.Data()).Data()<<"."<<std::endl;
    return;
  }
  
  // ----------------------------------------------------------------------
  // load the input workspace....
  RooWorkspace* win = (RooWorkspace*)fdata->Get("cms_hgg_workspace_data");
  if( !win ) {
    std::cerr<<" ERROR: Could not load workspace > cms_hgg_workspace_data < from file > "<<TString::Format("%s/CMS-HGG-data.root",projectDir.Data()).Data()<<" <."<<std::endl;
    return;
  }

  RooRealVar *intLumi = win->var("IntLumi");
  RooRealVar *hmass   = win->var("CMS_hgg_mass");
  if( !intLumi || !hmass ) {
    std::cerr<<" ERROR: Could not load needed variables > IntLumi < or > CMS_hgg_mass < forom input workspace."<<std::endl;
    return;
  }

  //win->Print();

  hmass->setRange(massmin,massmax);
  hmass->setBins(4*(int)(massmax-massmin));
  hmass->SetTitle("m_{#gamma#gamma}");
  hmass->setUnit("GeV");
  hmass->setRange("fitrange",massmin,massmax);

  hmass->setRange("blind1",100.,110.);
  hmass->setRange("blind2",150.,180.);
  
  // ----------------------------------------------------------------------
  // some auxiliray vectro (don't know the meaning of all of them ... yet...
  std::vector<RooAbsData*> data_vec;
  std::vector<RooAbsPdf*>  pdfShape_vec;   // vector to store the NOT-EXTENDED PDFs (aka pdfshape)
  std::vector<RooAbsPdf*>  pdf_vec;        // vector to store the EXTENDED PDFs
  
  std::vector<RooAbsReal*> normu_vec;      // this holds the normalization vars for each Cat (needed in bands for combined cat)

  RooArgList               normList;       // list of range-limityed normalizations (needed for error bands on combined category)

  //std::vector<RooRealVar*> coeffv;
  //std::vector<RooAbsReal*> normu_vecv; // ???

  // ----------------------------------------------------------------------
  // define output works
  RooWorkspace *wOut = new RooWorkspace("wbkg","wbkg") ;
  
  // util;ities for the combined fit
  RooCategory     finalcat  ("finalcat",  "finalcat") ;  
  RooSimultaneous fullbkgpdf("fullbkgpdf","fullbkgpdf",finalcat);
  RooDataSet      datacomb  ("datacomb",  "datacomb",  RooArgList(*hmass,finalcat)) ;

  RooDataSet *datacombcat = new RooDataSet("data_combcat","",RooArgList(*hmass)) ;
  
  // add the 'combcat' to the list...if more than one cat
  if( catnames.size() > 1 ) {
    catnames.push_back("combcat");    
    catdesc.push_back("Combined");
  }
  
  for (UInt_t icat=0; icat<catnames.size(); ++icat) {
    TString catname = catnames.at(icat);
    finalcat.defineType(catname);
    
    // check if we're in a sub-cat or the comb-cat
    RooDataSet *data   = NULL;
    RooDataSet *inData = NULL;
    if( icat < (catnames.size() - 1) || catnames.size() == 1) { // this is NOT the last cat (which is by construction the combination)
      inData = (RooDataSet*)win->data(TString("data_mass_")+catname);
      if( !inData ) {
	std::cerr<<" ERROR: Could not find dataset > data_mass_"<<catname.Data()<<" < in input workspace."<<std::endl;
	return;
      }
      data = new RooDataSet(TString("data_")+catname,"",*hmass,Import(*inData));  // copy the dataset (why?)
      
      // append the data to the combined data...
      RooDataSet *datacat = new RooDataSet(TString("datacat")+catname,"",*hmass,Index(finalcat),Import(catname,*data)) ;
      datacomb.append(*datacat);
      datacombcat->append(*data);
      
      // normalization for this category
      RooRealVar *nbkg = new RooRealVar(TString::Format("CMS_hgg_%s_bkgshape_norm",catname.Data()),"",800.0,0.0,25e3);
      
      // we keep track of the normalizario vars only for N-1 cats, naming convetnions hystoric...
      if( catnames.size() > 2 && icat < (catnames.size() - 2) ) {
	RooRealVar* cbkg = new RooRealVar(TString::Format("cbkg%s",catname.Data()),"",0.0,0.0,1e3);
	cbkg->removeRange();
	normu_vec.push_back(cbkg);
	normList.add(*cbkg);
      }
      
      /// generate the Bernstrin polynomial (FIX-ME: add possibility ro create other models...)
      fstBernModel* theBGmodel = new fstBernModel(hmass, polorder[icat], icat, catname);            // using my dedicated class...
      
      std::cout<<" model name is "<<theBGmodel->getPdf()->GetName()<<std::endl;

      RooAbsPdf*    bkgshape   = theBGmodel->getPdf();                                              // the BG shape
      RooAbsPdf*    bkgpdf     = new RooExtendPdf(TString("bkgpdf")+catname,"",*bkgshape,*nbkg);    // the extended PDF
      
      // add the extedned PDF to the RooSimultaneous holding all models...
      fullbkgpdf.addPdf(*bkgpdf,catname);
      // store the NON-EXTENDED PDF for usgae to compute the error bands later..
      pdfShape_vec.push_back(bkgshape);
      pdf_vec     .push_back(bkgpdf);
      data_vec    .push_back(data);
      
    } else {
      data = datacombcat;   // we're looking at the last cat (by construction the combination)
      data_vec.push_back(data);
      
      // sum up all the cts PDFs for combined PDF
      RooArgList subpdfs;
      for (int ipdf=0; ipdf<pdf_vec.size(); ++ipdf) {
	subpdfs.add(*pdf_vec.at(ipdf));
      }
      RooAddPdf* bkgpdf = new RooAddPdf(TString("bkgpdf")+catname,"",subpdfs);
      pdfShape_vec.push_back(bkgpdf);      
      pdf_vec     .push_back(bkgpdf);  // I don't think this is really needed though....
    }
    
    // generate the binned dataset (to be put into the workspace... just in case...)
    RooDataHist *databinned = new RooDataHist(TString("databinned_")+catname,"",*hmass,*data);
    
    wOut->import(*data);
    wOut->import(*databinned);

  }
  
  std::cout<<" ***************** "<<std::endl;

  // fit the RooSimultaneous to the combined dataset -> (we could also fit each cat separately)
  fullbkgpdf.fitTo(datacomb,Strategy(1),Minos(kFALSE),Save(kTRUE));
  RooFitResult *fullbkgfitres = fullbkgpdf.fitTo(datacomb,Strategy(2),Minos(kFALSE),Save(kTRUE));
  
  // in principle we're done now, so store the results in the output workspace
  wOut->import(datacomb);  
  wOut->import(fullbkgpdf);
  wOut->import(*fullbkgfitres);

  std::cout<<" ***************** "<<std::endl;
  

  if( verbose ) wOut->Print();

  
  std::cout<<" ***************** "<<std::endl;

  wOut->writeToFile("bkgdatawithfit.root") ;  
  
  if( verbose ) {
    printf("IntLumi = %5f\n",intLumi->getVal());
    printf("ndata:\n");
    for (UInt_t icat=0; icat<catnames.size(); ++icat) {    
      printf("%i ",data_vec.at(icat)->numEntries());      
    }   
    printf("\n");
  } 
  
  // --------------------------------------------------------------------------------------------
  // Now comesd the plotting
  // chage the Statistics style...
  gStyle->SetOptStat(1110);
  
  // we want to plot in 1GeV bins (apparently...)
  UInt_t nbins = (UInt_t) (massmax-massmin);
  
  // here we'll store the curves for the bands...
  std::vector<RooCurve*> fitcurves;
  
  // loop again over the cats
  TCanvas **canbkg = new TCanvas*[catnames.size()];
  RooPlot** plot   = new RooPlot*[catnames.size()];

  TLatex** lat  = new TLatex*[catnames.size()];
  TLatex** lat2 = new TLatex*[catnames.size()];

  std::cout<<"  beofre plotting..."<<std::endl;
  

  for (UInt_t icat=0; icat<catnames.size(); ++icat) {
    TString catname = catnames.at(icat);
    

    std::cout<<" trying to plot #"<<icat<<std::endl;

    // plot the data and the fit 
    canbkg[icat] = new TCanvas;
    plot  [icat] = hmass->frame(Bins(nbins),Range("fitrange"));
    
    std::cout<<" trying to plot #"<<icat<<std::endl;

    // first plot the data invisibly... and put the fitted BG model on top...
    data_vec    .at(icat)->plotOn(plot[icat],RooFit::LineColor(kWhite),MarkerColor(kWhite),Invisible());
    pdfShape_vec.at(icat)->plotOn(plot[icat],RooFit::LineColor(kRed),Range("fitrange"),NormRange("fitrange"));
    
    std::cout<<" trying to plot #"<<icat<<std::endl;


    // if toggled on, plot also the Data visibly
    if( !blinded ) {
      data_vec.at(icat)->plotOn(plot[icat]);
    }
   
    std::cout<<" trying to plot #"<<icat<<std::endl;

    // some cosmetics...
    plot[icat]->SetTitle("");      
    plot[icat]->SetMinimum(0.0);
    plot[icat]->SetMaximum(1.40*plot[icat]->GetMaximum());
    plot[icat]->GetXaxis()->SetTitle("m_{#gamma#gamma} (GeV/c^{2})");
    plot[icat]->Draw();       
            

    std::cout<<" trying to plot #"<<icat<<std::endl;

    // legend....
    TLegend *legmc = new TLegend(0.68,0.70,0.97,0.90);
    legmc->AddEntry(plot[icat]->getObject(2),"Data","LPE");
    legmc->AddEntry(plot[icat]->getObject(1),"Bkg Model","L");
    
    // this part computes the 1/2-sigma bands.    
    TGraphAsymmErrors *onesigma = NULL;
    TGraphAsymmErrors *twosigma = NULL;
    
    std::cout<<" trying ***  to plot #"<<icat<<std::endl;

    RooAddition* sumcatsnm1 = NULL;

    if ( dobands ) { //&& icat == (catnames.size() - 1) ) {

      onesigma = new TGraphAsymmErrors();
      twosigma = new TGraphAsymmErrors();

      // get the PDF for this cat from the vector
      RooAbsPdf *thisPdf = pdfShape_vec.at(icat); 

      // get the nominal fir curve
      RooCurve *nomcurve = dynamic_cast<RooCurve*>(plot[icat]->getObject(1));
      fitcurves.push_back(nomcurve);

      bool iscombcat       = ( icat == (catnames.size() - 1) && catnames.size() > 1);
      RooAbsData *datanorm = ( iscombcat ? &datacomb : data_vec.at(icat) );

      // this si the nornmalization in the 'sliding-window' (i.e. per 'test-bin')
      RooRealVar *nlim = new RooRealVar(TString::Format("nlim%s",catnames.at(icat).Data()),"",0.0,0.0,10.0);
      nlim->removeRange();

      if( iscombcat ) {
	// ----------- HISTORIC NAMING  ----------------------------------------
	sumcatsnm1 = new RooAddition("sumcatsnm1","",normList);   // summing all normalizations epect the last Cat
	// this is the normlization of the last Cat
	RooFormulaVar *nlast = new RooFormulaVar("nlast","","TMath::Max(0.1,@0-@1)",RooArgList(*nlim,*sumcatsnm1));
	// ... and adding it ot the list of norms
	normu_vec.push_back(nlast);
      }

      //if (icat == 1 && catnames.size() == 2) continue; // only 1 cat, so don't need combination

      for (int i=1; i<(plot[icat]->GetXaxis()->GetNbins()+1); ++i) {
	
	// this defines the 'binning' we use for the error bands
        double lowedge = plot[icat]->GetXaxis()->GetBinLowEdge(i);
        double upedge = plot[icat]->GetXaxis()->GetBinUpEdge(i);
        double center = plot[icat]->GetXaxis()->GetBinCenter(i);
        
	// get the nominal value at the center of the bin
        double nombkg = nomcurve->interpolate(center);
        nlim->setVal(nombkg);
        hmass->setRange("errRange",lowedge,upedge);

	// this is the new extended PDF whith the normalization restricted to the bin-area
        RooAbsPdf *extLimPdf = NULL;
	if( iscombcat ) {
	  extLimPdf = new RooSimultaneous("epdf","",finalcat);
	  // loop over the cats and generate temporary extended PDFs
	  for (int jcat=0; jcat<(catnames.size()-1); ++jcat) {
            RooRealVar *rvar = dynamic_cast<RooRealVar*>(normu_vec.at(jcat));
            if (rvar) rvar->setVal(fitcurves.at(jcat)->interpolate(center));
            RooExtendPdf *ecpdf = new RooExtendPdf(TString::Format("ecpdf%s",catnames.at(jcat).Data()),"",*pdfShape_vec.at(jcat),*normu_vec.at(jcat),"errRange");
            static_cast<RooSimultaneous*>(extLimPdf)->addPdf(*ecpdf,catnames.at(jcat));
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
      
      // plot[icat] the error bands
      twosigma->SetLineColor(kGreen);
      twosigma->SetFillColor(kGreen);
      twosigma->SetMarkerColor(kGreen);
      twosigma->Draw("L3 SAME");     
      
      onesigma->SetLineColor(kYellow);
      onesigma->SetFillColor(kYellow);
      onesigma->SetMarkerColor(kYellow);
      onesigma->Draw("L3 SAME");
      
      plot[icat]->Draw("SAME");
    
      // and add the error bands to the legend
      legmc->AddEntry(onesigma,"#pm1 #sigma","F");  
      legmc->AddEntry(twosigma,"#pm2 #sigma","F");  
    }
    
    std::cout<<" trying ***2  to plot #"<<icat<<std::endl;

    // rest of the legend ....
    legmc->SetBorderSize(0);
    legmc->SetFillStyle(0);
    legmc->Draw();   

    lat[icat]  = new TLatex(103.0,0.9*plot[icat]->GetMaximum(),TString::Format("#scale[0.7]{#splitline{CMS preliminary}{#sqrt{s} = %.1f TeV L = %.2f fb^{-1}}}",theCMenergy,intLumi->getVal()));
    lat2[icat] = new TLatex(103.0,0.75*plot[icat]->GetMaximum(),catdesc.at(icat));

    lat[icat] ->Draw();
    lat2[icat]->Draw();
    
    // -------------------------------------------------------    
    // save canvas in different formats
    canbkg[icat]->SaveAs(TString("databkg") + catname + TString(".pdf"));
    canbkg[icat]->SaveAs(TString("databkg") + catname + TString(".eps"));
    canbkg[icat]->SaveAs(TString("databkg") + catname + TString(".root"));              
  }
  
  return;
  
}


// ----------------------------------------------------------------------------------------

bool readFromConfigCard(TString fileName,
			TString& projectDir,
			std::vector<TString>& catNames,
			std::vector<TString>& catDesc,
			std::vector<int>    & polOrders,
			double& massmin,
			double& massmax,
			double& theCMenergy    
			) {
  
  catNames.resize(0);
  catDesc.resize(0);
  polOrders.resize(0);

  FILE* configFile = fopen(fileName.Data(),"r");
  if ( !configFile ) {
    std::cerr<<" Inputfile "<<fileName<<" not found."<<std::endl;
    return false;
  }
  
  char line[200];
  
  std::cout<<" Reading weight information from file "<<fileName<<"...";
  
  while (fgets(line,200,configFile)) {
    char name[200];
    int catIdx = -1;
    float smear = -1.;
    int polOrder = -1;
    float massval = -1.;
    if( sscanf(line,"PROJECTDIR %s",&name) ) projectDir = TString(name);
    else if( sscanf(line,"MINMSS %f",&massval) ) massmin = (double) massval;
    else if( sscanf(line,"MAXMSS %f",&massval) ) massmax = (double) massval;
    else if( sscanf(line,"CMENERGY %f",&massval) ) theCMenergy = (double) massval;
    else if( sscanf(line,"ANACAT %d %s %f Bern/%d",&catIdx, &name, &smear, &polOrder) ) {
      if( catIdx != (int) catNames.size() ) {
	std::cerr<<" ERROR: ANACAT categories not in order."<<std::endl;
	return false;
      }
      catNames.push_back(TString(name));
      polOrders.push_back(polOrder);

      std::string theLine = line;
      size_t fPos = theLine.rfind("Desc(");
      if( fPos == std::string::npos ) {
	std::cerr<<" ERROR: ANACAT "<<catIdx<<" with name "<<name<<" has no string descriptor (StDexc(...))."<<std::endl;
	return false;
      }
      theLine.erase(0,fPos+5);
      fPos = theLine.find_last_of(")");
      if( fPos == std::string::npos ) {
	std::cerr<<" ERROR: ANACAT "<<catIdx<<" with name "<<name<<" has no valid string descriptor (StDexc(...))."<<std::endl;
	return false;
      }
      theLine.erase(fPos);
      catDesc.push_back(TString(theLine));
    }
  }
  
  fclose(configFile);
  std::cout<<" done."<<std::endl;

  return true;
}

