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
#include "RooConstVar.h"
#include "RooMinimizer.h"
#include "TGaxis.h"
#include "TLatex.h"
#include "TArrow.h"

using namespace RooFit;

//effsigma function from Chris
Double_t effSigma(TH1 * hist)
{

  TAxis *xaxis = hist->GetXaxis();
  Int_t nb = xaxis->GetNbins();
  if(nb < 10) {
    cout << "effsigma: Not a valid histo. nbins = " << nb << endl;
    return 0.;
  }
  
  Double_t bwid = xaxis->GetBinWidth(1);
  if(bwid == 0) {
    cout << "effsigma: Not a valid histo. bwid = " << bwid << endl;
    return 0.;
  }
  Double_t xmax = xaxis->GetXmax();
  Double_t xmin = xaxis->GetXmin();
  Double_t ave = hist->GetMean();
  Double_t rms = hist->GetRMS();

  Double_t total=0.;
  for(Int_t i=0; i<nb+2; i++) {
    total+=hist->GetBinContent(i);
  }
//   if(total < 100.) {
//     cout << "effsigma: Too few entries " << total << endl;
//     return 0.;
//   }
  Int_t ierr=0;
  Int_t ismin=999;
  
  Double_t rlim=0.683*total;
  Int_t nrms=rms/(bwid);    // Set scan size to +/- rms
  if(nrms > nb/10) nrms=nb/10; // Could be tuned...

  Double_t widmin=9999999.;
  for(Int_t iscan=-nrms;iscan<nrms+1;iscan++) { // Scan window centre
    Int_t ibm=(ave-xmin)/bwid+1+iscan;
    Double_t x=(ibm-0.5)*bwid+xmin;
    Double_t xj=x;
    Double_t xk=x;
    Int_t jbm=ibm;
    Int_t kbm=ibm;
    Double_t bin=hist->GetBinContent(ibm);
    total=bin;
    for(Int_t j=1;j<nb;j++){
      if(jbm < nb) {
        jbm++;
        xj+=bwid;
        bin=hist->GetBinContent(jbm);
        total+=bin;
        if(total > rlim) break;
      }
      else ierr=1;
      if(kbm > 0) {
        kbm--;
        xk-=bwid;
        bin=hist->GetBinContent(kbm);
        total+=bin;
        if(total > rlim) break;
      }
      else ierr=1;
    }
    Double_t dxf=(total-rlim)*bwid/bin;
    Double_t wid=(xj-xk+bwid-dxf)*0.5;
    if(wid < widmin) {
      widmin=wid;
      ismin=iscan;
    }   
  }
  if(ismin == nrms || ismin == -nrms) ierr=3;
  if(ierr != 0) cout << "effsigma: Error of type " << ierr << endl;
  
  return widmin;
  
}


//pileup reweighting
TH1D *puweights = 0;
float puweight(float npu) {
  if (npu<0) return 1.0;
  return puweights->GetBinContent(puweights->FindFixBin(npu));
}

//pt-reweighing (but only for gf samples, use extra bool flag for now)
TH1D *ptweights = 0;
float ptweight(float genhpt, bool isgfsample) {
  if (!isgfsample) return 1.0;
  if (genhpt<0) return 1.0;
  return ptweights->GetBinContent(ptweights->FindFixBin(genhpt));
}

//efficiency scale factors (barrel/endcap for now, but can be extended to full kinematic binning)
float effweight(bool iseb1, bool iseb2) {
  const float ebscale = 0.992;
  const float eescale = 0.951;
  
  float effw = 1.0;
  
  if (iseb1) effw *= ebscale;
  else effw *= eescale;
  
  if (iseb2) effw *= ebscale;
  else effw *= eescale;  
  
  return effw;
}

//append data from file to RooDataSet adding a column which has the weight by given cross section
//(Note that the weight is not enabled at this stage, only made available for subsequent use)
void append(RooDataSet &data, TFile *infile, TCut sel, Double_t xsec, Bool_t isgfsample) {

    TDirectory *hdirfwk = (TDirectory*) infile->FindObjectAny("AnaFwkMod");
    const TH1D *hDAllEvents = (TH1D*)hdirfwk->Get("hDAllEvents");
  
  
    TDirectory *hdir = (TDirectory*) infile->FindObjectAny("HGGMod");
    TTree *hdata = (TTree*)hdir->Get("hHggNtuple");  
    

    RooRealVar xsecweight("xsecweight","xsecweight",1e6*xsec/(double)hDAllEvents->GetEntries());
    RooRealVar isgf("isgf","isgf",(double)isgfsample);
    
    
    RooArgSet varlistsmall = *data.get();
    varlistsmall.remove(xsecweight,kFALSE,kTRUE);
    varlistsmall.remove(isgf,kFALSE,kTRUE);
    
    
    RooDataSet newdata("newdata","newdata",varlistsmall,RooFit::Import(*hdata),RooFit::Cut(sel)); 
    newdata.addColumn(xsecweight);
    newdata.addColumn(isgf);
    
    
    data.append(newdata);
    
    
}

RooAbsReal *gcdf=0;
RooRealVar *gcdfvar=0;
float cdfinterval(float center, float width) {
  gcdfvar->setVal(center+width);
  float cdfhigh = gcdf->getVal();
  gcdfvar->setVal(center-width);
  float cdflow = gcdf->getVal();

  printf("cdfinterval = %5f\n",cdfhigh-cdflow);
  
  return (cdfhigh-cdflow);
}

RooDataSet *cwdset(RooDataSet *indata, RooRealVar *mvar, RooRealVar *wvar, TString name, Double_t weightscale) {
  
    RooDataSet *outdata = new RooDataSet(name,"",RooArgList(*mvar,*wvar),wvar->GetName());
    for (Int_t ient=0; ient<indata->numEntries(); ++ient) {
      const RooArgSet *ent = indata->get(ient);
      mvar->setVal(static_cast<RooAbsReal*>(ent->find(mvar->GetName()))->getVal());
      outdata->add(*mvar,weightscale*indata->weight());
    }
    
    return outdata;
    
}


void plothmasssplit() { 
  
  
  gROOT->Macro("MitStyle.C");
  gStyle->SetErrorX(0); 
  gStyle->SetOptStat(0);
  gROOT->ForceStyle();  
  
  
  //Load pileup weights
  //TFile *filepuest = new TFile("/home/bendavid/cms/root/pudists/puweightsDec22.root","READ");
  //TFile *filepuest = new TFile("/home/bendavid/cms/root/pudists/puweightsGoodCollisions2011.root","READ");
  //TFile *filepuest = new TFile("/home/bendavid/cms/root/pudists/puweightsCombined219.root","READ");
  TFile *filepuest = new TFile("/home/bendavid/cms/root/pudists/PuWeightsMay19Combined.root","READ");
  puweights= (TH1D*) filepuest->Get("pileup");

  //load pt-weights
  TFile *fileptweight = new TFile("/home/bendavid/cms/root/pudists/Kfactors.root","READ");
  
  
  gSystem->cd("./bambuOutputFeb16smmvavbf/signalplots");
  gStyle->SetOptStat(1110);
  
 
  

    TFile *fdata = new TFile("/scratch/bendavid/root/bambuOutputFeb16smmvavbf/CMS-HGG-mclow.root","READ");

    
    TFile *fextraggh = new TFile("/scratch/bendavid/root/bambuOutputFeb16smmvavbf/ggh/extra.root","READ");    
    TFile *fextravbf = new TFile("/scratch/bendavid/root/bambuOutputFeb16smmvavbf/vbf/extra.root","READ");    
    TFile *fextrawzh = new TFile("/scratch/bendavid/root/bambuOutputFeb16smmvavbf/wzh/extra.root","READ");    
    TFile *fextratth = new TFile("/scratch/bendavid/root/bambuOutputFeb16smmvavbf/tth/extra.root","READ");    
    
    
    TFile *fpdfsggh = new TFile("/scratch/bendavid/root/bambuOutputFeb16smmvavbf/ggh/ubersignalmodel.root","READ");
    TFile *fpdfsvbf = new TFile("/scratch/bendavid/root/bambuOutputFeb16smmvavbf/vbf/ubersignalmodel.root","READ");
    TFile *fpdfswzh = new TFile("/scratch/bendavid/root/bambuOutputFeb16smmvavbf/wzh/ubersignalmodel.root","READ");
    TFile *fpdfstth = new TFile("/scratch/bendavid/root/bambuOutputFeb16smmvavbf/tth/ubersignalmodel.root","READ");
   
   
    RooWorkspace* win = (RooWorkspace*)fdata->Get("cms_hgg_workspace_mclow");
    
    RooWorkspace* wextraggh = (RooWorkspace*)fextraggh->Get("wextra");
    RooWorkspace* wextravbf = (RooWorkspace*)fextravbf->Get("wextra");
    RooWorkspace* wextrawzh = (RooWorkspace*)fextrawzh->Get("wextra");
    RooWorkspace* wextratth = (RooWorkspace*)fextratth->Get("wextra");

    RooWorkspace* wsigggh = (RooWorkspace*)fpdfsggh->Get("wsig");
    RooWorkspace* wsigvbf = (RooWorkspace*)fpdfsvbf->Get("wsig");    
    RooWorkspace* wsigwzh = (RooWorkspace*)fpdfswzh->Get("wsig");
    RooWorkspace* wsigtth = (RooWorkspace*)fpdfstth->Get("wsig");    
    
    RooWorkspace *wcomb = new RooWorkspace("wcomb","");
    
    //win->Print();
//     wsig->Print();
//       wextra->Print();
//      return;
    
    RooRealVar *hmass = win->var("CMS_hgg_mass");
    //RooRealVar *hmass = win->var("mass");  
    hmass->setRange(100,180);
    hmass->setBins(320);
    hmass->SetTitle("m_{#gamma#gamma}");
    hmass->setUnit("GeV/c^{2}");
  
    RooRealVar *mnom = wsigggh->var("MH");
    
  std::vector<TString> catdesc;
  
  //Feb2 smmvavbf
    catdesc.push_back("#scale[0.5]{BDT >= 0.89}");
    catdesc.push_back("#scale[0.5]{0.72 <= BDT < 0.89}");
    catdesc.push_back("#scale[0.5]{0.55 <= BDT < 0.72}");
    catdesc.push_back("#scale[0.5]{0.05 <= BDT < 0.55}");
     catdesc.push_back("#scale[0.5]{BDT >= 0.05 VBF Tag}");
    catdesc.push_back("#scale[0.5]{All Categories Combined}");    


  //define categories
  std::vector<TString> catnames;
  std::vector<TCut> catcuts;
  
  catnames.push_back("cat0");
  catnames.push_back("cat1");
  catnames.push_back("cat2");
  catnames.push_back("cat3");
  catnames.push_back("cat4");
  
  catnames.push_back("combcat");



  std::vector<RooRealVar*> fitparms;
  std::vector<double> fitparmsinit;

  

  RooRealVar *weight = new RooRealVar("weight","",1.0);
  weight->removeRange();
  const double weightscale = 1.0;

  double sumofweights=0.0;

  std::vector<double> effaccs;
  std::vector<double> sumweights;
  std::vector<double> numevents;
  
  std::vector<double> effsigmas;
  
  std::vector<double> effsigmashist;
  std::vector<double> effsigmapshist;
  
  
  const float testmass = 120.0;
    
  TString numstring="120";
  
  //RooRealVar *mtotalxsec = win->var(TString("XSBR_")+numstring);
  
  hmass->setBins(6000);

  //return;
  RooArgList catpdfs;  
  for (UInt_t icat=0; icat<catnames.size(); ++icat) {

    printf("begin loop\n");

    hmass->setRange("higgsrange",TMath::Max(100.0,testmass-20.0),TMath::Min(160.0,testmass+15.0));
    mnom->setVal(testmass);  

  
    RooDataSet *mcsigalldata = 0;
    
    
    printf("combcat?\n");
    if (catnames.at(icat)=="combcat") {
      printf("using combined dataset");
      mcsigalldata = (RooDataSet*)wextraggh->data(TString("sig_mass_m")+numstring+TString("_")+catnames.at(icat)+TString("_ggh"));
      
      RooDataSet *mcsigalldatavbf = (RooDataSet*)wextravbf->data(TString("sig_mass_m")+numstring+TString("_")+catnames.at(icat)+TString("_vbf"));
      RooDataSet *mcsigalldatawzh = (RooDataSet*)wextrawzh->data(TString("sig_mass_m")+numstring+TString("_")+catnames.at(icat)+TString("_wzh"));
      RooDataSet *mcsigalldatatth = (RooDataSet*)wextratth->data(TString("sig_mass_m")+numstring+TString("_")+catnames.at(icat)+TString("_tth"));
      
      printf("mcsigalldata = %p, mcsigalldatavbf = %p",(void*)mcsigalldata,(void*)mcsigalldatavbf);
      
      mcsigalldata->append(*mcsigalldatavbf);
      mcsigalldata->append(*mcsigalldatawzh);
      mcsigalldata->append(*mcsigalldatatth);

    }
    else {
      RooDataSet *mcsigalldataggh = (RooDataSet*)win->data(TString("sig_ggh_mass_m")+numstring+TString("_")+catnames.at(icat));
      RooDataSet *mcsigalldatavbf = (RooDataSet*)win->data(TString("sig_vbf_mass_m")+numstring+TString("_")+catnames.at(icat));
      RooDataSet *mcsigalldatawzh = (RooDataSet*)win->data(TString("sig_wzh_mass_m")+numstring+TString("_")+catnames.at(icat));
      RooDataSet *mcsigalldatatth = (RooDataSet*)win->data(TString("sig_tth_mass_m")+numstring+TString("_")+catnames.at(icat));
      
      mcsigalldata = mcsigalldataggh;
      mcsigalldataggh->append(*mcsigalldatavbf);
      mcsigalldataggh->append(*mcsigalldatawzh);
      mcsigalldataggh->append(*mcsigalldatatth);

    }
    
    TH1 *tmphist = mcsigalldata->createHistogram("tmphist",*hmass);
    double effsigmahist = effSigma(tmphist);
    delete tmphist;
        
    effsigmashist.push_back(effsigmahist);
    
    sumweights.push_back(mcsigalldata->sumEntries());
    numevents.push_back(mcsigalldata->numEntries());
    
    //RooDataSet *mcsigallwdata = cwdset(mcsigalldata,hmass,weight,TString("mcsigallwdata") + numstring+catnames.at(icat),weightscale);
    //RooAbsPdf  *fullpdf = wsig->pdf(TString("hggpdf_")+catnames.at(icat));
    //RooAbsReal *fullnorm = wsig->function(TString("hggpdf_")+catnames.at(icat) + TString("_norm"));
    
    printf("getting pdfs\n");
    
    RooAddPdf *fullpdf = 0;
    
    if (catnames.at(icat)=="combcat") {
      
      fullpdf = new RooAddPdf(TString::Format("fullpdf_%s",catnames.at(icat).Data()),"",catpdfs);
        
      
    }
    else {
      RooAbsPdf *fullpdfggh = wsigggh->pdf(TString::Format("sigpdfsmrel%s_ggh",catnames.at(icat).Data()));
      RooAbsPdf *fullpdfvbf = wsigvbf->pdf(TString::Format("sigpdfsmrel%s_vbf",catnames.at(icat).Data()));
      RooAbsPdf *fullpdfwzh = wsigwzh->pdf(TString::Format("sigpdfsmrel%s_wzh",catnames.at(icat).Data()));
      RooAbsPdf *fullpdftth = wsigtth->pdf(TString::Format("sigpdfsmrel%s_tth",catnames.at(icat).Data()));
      
      printf("fullpdfggh = %p\n",(void*)fullpdfggh);
      
      fullpdf = new RooAddPdf(TString::Format("fullpdf_%s",catnames.at(icat).Data()),"",RooArgList(*fullpdfggh,*fullpdfvbf,*fullpdfwzh,*fullpdftth));
      
      catpdfs.add(*fullpdf);
    }
    //effaccs.push_back(fullnorm->getVal());

    TH1 *tmphistp = fullpdf->createHistogram("tmphistp",*hmass);
    double effsigmaphist = effSigma(tmphistp);
    effsigmapshist.push_back(effsigmaphist);
    delete tmphistp;
    
    wcomb->import(*fullpdf,RecycleConflictNodes());
    
    printf("mcsigalldata = %i, fullpdf = %i\n",mcsigalldata!=0, fullpdf!=0);
    
    RooGaussian testgaus("testgaus","",*hmass,RooConst(120.),RooConst(1.0));

    RooAbsPdf *wrapdf = fullpdf;
    
    
    //calculate effective sigma and corresponding interval
    RooAbsReal *cdf = fullpdf->createCdf(*hmass, ScanAllCdf());
    //RooAbsReal *cdf = wrapdf.createCdf(*hmass);
    
    new TCanvas;
    RooPlot *cplot = hmass->frame(Bins(100),Range("higgsrange"));
    cdf->plotOn(cplot);
    wrapdf->plotOn(cplot);
    cplot->Draw();
    
  //  return;
    
    
    gcdf = cdf;
    gcdfvar = hmass;

    hmass->setVal(160.0);
    double maxcdf = cdf->getVal();
    float center = testmass-10.0;
    float minwidth = 999.0;
    float mlmin = 0.0;
    float mhmin = 0.0;
    float step=0.01;
    for (int i=0; i<2000; ++i) {
      float mlow = center+i*0.01;
      hmass->setVal(mlow);
      float cdflo = cdf->getVal();
      for (int j=i; j<2000; ++j) {
        float mhigh = center+j*0.01;
        hmass->setVal(mhigh);
        float cdfhi = cdf->getVal();
        if ( (cdfhi-cdflo)>0.683 ) {
          if ( (mhigh-mlow)<minwidth) {
            minwidth = mhigh-mlow;
            mlmin = mlow;
            mhmin = mhigh;
          }
          break;
        }
      }
      
    }
    float sigmaeff = minwidth/2.0;
    effsigmas.push_back(sigmaeff);   
   

    RooRealVar *clonedmass = (RooRealVar*)wrapdf->getDependents(&RooArgSet(*hmass))->find("CMS_hgg_mass");

    
    

    
    printf("sigmaeff = %5f\n",sigmaeff);

    
    sumofweights+=mcsigalldata->sumEntries();
    TGaxis::SetMaxDigits(2);
    
    //TH1 *hmcall = mcsigalldata->createHistogram(TString::Format("hmcall_%s",catnames.at(icat).Data()),*hmass,Binning(100,100.0,135.0));
    
    Int_t numbins = 100;
    TCanvas *ccompint = new TCanvas;
    TString plotname = TString("effsigma") + catnames.at(icat) + TString(".eps");
    RooPlot *hplotcompint = hmass->frame(Bins(numbins),Range("higgsrange"));
    mcsigalldata->plotOn(hplotcompint,Invisible());
    //fullpdf->plotOn(hplotcompint,DrawOption("F"),LineColor(kGray),FillStyle(3344),FillColor(kGray),NormRange("higgsrange"),Range(mlmin,mhmin),VLines());
    fullpdf->plotOn(hplotcompint,DrawOption("F"),LineColor(kGray),FillStyle(1001),FillColor(19),NormRange("higgsrange"),Range(mlmin,mhmin),VLines());      
    fullpdf->plotOn(hplotcompint,LineWidth(1.5),LineColor(kGray),NormRange("higgsrange"),Range(mlmin,mhmin),VLines());      
    //mcsigalldata->plotOn(hplotcompint,MarkerStyle(),LineWidth(2));
    mcsigalldata->plotOn(hplotcompint,MarkerStyle(25));
    //hmcall->Draw("Same");
    fullpdf->plotOn(hplotcompint,RooFit::LineColor(kBlue),NormRange("higgsrange"),Range("higgsrange"));
    if (0) {
      mnom->setVal(testmass+5.0);
      fullpdf->plotOn(hplotcompint,RooFit::LineColor(kOrange),NormRange("higgsrange"),Range("higgsrange"));
      mnom->setVal(testmass);
    }
    
    //hplotcompint->SetMaximum(1.0*hplotcompint->GetMaximum());
    hplotcompint->SetTitle("");
    hplotcompint->GetXaxis()->SetNoExponent(kTRUE);
    hplotcompint->GetXaxis()->SetTitle("m_{#gamma#gamma} (GeV/c^{2})");
    //hplotcompint->GetYaxis()->SetTitle("m_{#gamma#gamma} (GeV/c^{2})");
    hplotcompint->Draw();    
    
    TLegend *legmc = new TLegend(0.18,0.6,0.50,0.9);  
    legmc->AddEntry(hplotcompint->getObject(3),"Simulation","LPE");
    //legmc->AddEntry(hmcall,"Smeared MC","LE");
    legmc->AddEntry(hplotcompint->getObject(4),"Parametric Model","L");
    legmc->AddEntry(hplotcompint->getObject(1),TString::Format("#sigma_{eff} = %.2f GeV/c^{2}",sigmaeff),"LF");
    legmc->SetBorderSize(0);
    legmc->SetFillStyle(0);
    legmc->Draw();          
    
    TLatex lat(123.0,0.9*hplotcompint->GetMaximum(),"#scale[0.7]{#splitline{CMS preliminary}{Simulation}}");
    lat.Draw();
    
    TLatex lat2(125.5,0.75*hplotcompint->GetMaximum(),catdesc.at(icat));
    lat2.Draw();    
    
    //double plotnorm = mcsigalldata->sumEntries()*35.0/(double)numbins;
    double plotnorm = hplotcompint->getFitRangeNEvt()*hplotcompint->getFitRangeBinW();
    
    //double fwhmpos = plotnorm*phalf;
    RooCurve *mcurve = dynamic_cast<RooCurve*>(hplotcompint->getObject(4));
    
    double pmax = mcurve->getYAxisMax();
    double phalf = pmax/2.0;
    double fwhmpos = phalf;    
    
    
    double mmax = 0.0;
    double mleft = 0.0;
    double mright = 0.0;
    double dpminleft = 1e6;
    double dpminright = 1e6;
    double pmaxtmp = 0.0;
    for (int i=0; i<2000; ++i) {
      double mlow = center+i*0.01;
      
      double pval = mcurve->interpolate(mlow);
      //printf("pval = %5e\n",pval);
      if (pval>pmaxtmp) {
        pmaxtmp = pval;
        mmax = mlow;
      }
    }
    
    printf("pmax = %5f\n",pmax);
    
    
    for (int i=0; i<2000; ++i) {
      double mlow = center+i*0.01;

      double pval = mcurve->interpolate(mlow);
      double dp = TMath::Abs(pval-phalf);
      
      if (dp<dpminleft && mlow<mmax) {
        dpminleft = dp;
        mleft = mlow;
      }
      
      if (dp<dpminright && mlow>mmax) {
        dpminright = dp;
        mright = mlow;
      }      
      
    }    
    double fwhm = mright-mleft;    
    printf("mleft = %5f, mright = %5f, phalf = %5f\n",mleft,mright,phalf);
    
    //TLine fwhmarrow(mleft,100e-3,mright,100e-3);
    //TLine fwhmarrow(0.1,0.1,0.2,0.2,0.05,"<>");
    printf("plotnorm = %5f\n",plotnorm);
    TArrow fwhmarrow;
    fwhmarrow.SetLineWidth(2);
    fwhmarrow.DrawArrow(mleft,fwhmpos,mright,fwhmpos,0.025,"<>");
    TLatex latfwhm(mleft-14.0,fwhmpos,TString::Format("#scale[0.6]{FWHM = %2g GeV/c^{2}}",fwhm));
    latfwhm.Draw("");    
    
    ccompint->SaveAs(TString("effsigma") + catnames.at(icat) + TString(".eps"));
    ccompint->SaveAs(TString("effsigma") + catnames.at(icat) + TString(".pdf"));
    
    
    //ccompint->SaveAs(TString("effsigma") + catnames.at(icat) + TString(".png"));
    ccompint->SaveAs(TString("effsigma") + catnames.at(icat) + TString(".root"));
    
    //return;
    
    
    //return;

    
  }
  
  wcomb->writeToFile("combsignal.root") ;

  
  for (int i=0; i<effsigmashist.size(); ++i) {
    printf("effsigma = %5f, effsigmahist = %5f, effsigmaphist = %5f\n", effsigmas.at(i),   effsigmashist.at(i), effsigmapshist.at(i));

  }
  
  
}