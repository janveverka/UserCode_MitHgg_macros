#include <cstdlib>
#include <iostream> 
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TPluginManager.h"

//#include "TMVAGui.C"

#include "TRandom.h"

#include "GBRTrainer.h"
#include "GBRForest.h"
#include "GBRApply.h"
#include "TCut.h"
#include "TH2D.h"
#include "TChain.h"
#include "TObjArray.h"
#include "TTreeFormula.h"
#include "TCanvas.h"
#include "Cintex/Cintex.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraphAsymmErrors.h"
#include "TLine.h"
#include "THStack.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooCBShape.h"
#include "RooGaussian.h"
#include "RooProdPdf.h"
#include "RooAddPdf.h"
#include "RooMinuit.h"
#include "RooNLLVar.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TH1.h"
#include "RooPoisson.h"
#include "RooConstVar.h"
#include "RooDataHist.h"
#include "RooNDKeysPdf.h"
#include "TRandom.h"
#include "TMath.h"
#include "THnSparse.h"
#include "TFile.h"
#include "TH2D.h"
#include "TTree.h"
#include "TBranchElement.h"
#include "TTreeFormula.h"
#include "TROOT.h"
#include "TMVA/Tools.h"
#include "TMVA/Factory.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TSystem.h"
#include "TF1.h"
#include "THnSparse.h"
#include "RooFFTConvPdf.h"
#include "RooBreitWigner.h"
#include "RooFit.h"

using namespace RooFit;


float rndmassapply(float vtxprob, float sigmaright, float sigmawrong) {
  double rndu = gRandom->Uniform();
  return rndu >= vtxprob ? fabs(sigmaright*gRandom->Gaus(0.,1.)) : fabs(sigmawrong*gRandom->Gaus(0.,1.));
}

TH1D *getpuweights(TFile *file, TH1D *target) {
  
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
  
  TH1D *puweights = new TH1D((*htargettmp)/(*hpumc));
    
  delete htargettmp;
  
  return puweights;
    
}

double getweight(TFile *file, double xsec) {
 
  TDirectory *dir = (TDirectory*)file->FindObjectAny("AnaFwkMod");
  TH1D *hallevts = (TH1D*)dir->Get("hDAllEvents");
  
  return xsec/hallevts->GetSumOfWeights();
  
}


TH1D *puweights[50];
float puweight(float npu, int wset=0) {
  if (npu<0) return 1.0;
  return puweights[wset]->GetBinContent(puweights[wset]->FindFixBin(npu));
}

float xsecweights[50];
float xsecweight(int procidx=0) {
  return xsecweights[procidx];
}

void initweights(TChain *chain, TH1D *target, float *xsecs, float lumi) {
 
  TObjArray *files = chain->GetListOfFiles();
  for (int i=0; i<files->GetEntries(); ++i) {    
    TFile *file = TFile::Open(files->At(i)->GetTitle(),"READ");
    
    puweights[i] = getpuweights(file,target);
    xsecweights[i] = getweight(file,lumi*xsecs[i]);
    
    file->Close();    
  } 
  
  chain->SetAlias("procidx","This->GetTreeNumber()");
  
}


void fillcbhists(TString name, TString var, TCut cut, RooAbsPdf &zpdf,RooRealVar &mass, RooRealVar &sigma, RooRealVar &m0, RooRealVar *vset, int i, TH1F *hresmc, TH1F* hresd, TH1F *hress, TH1F *hmmc, TH1F* hmd, TTree *tmc, TTree *tdata, RooDataSet &tmpdset) {
      double resscale = 1.0;
      double smearscale = 2.0/sqrt(2.0)/91.1876;
  

      
      tmpdset.reset();
      Int_t nwe = tmc->Draw(var, cut);
      double *vals = tmc->GetV1();
      double *weights = tmc->GetW();
      for (UInt_t j=0; j<nwe; ++j) {
        vset->setVal(vals[j]);
        if (vset->getVal()>65. && vset->getVal()<115)
          tmpdset.add(*vset,weights[j]);
      }
      double sigmamc = -99.;
      double sigmamcerr = -99.;      
      if (tmpdset.sumEntries()>10) {
        zpdf.fitTo(tmpdset,Strategy(1),NumCPU(16));
        sigmamc = TMath::Abs(sigma.getVal())*resscale;
        sigmamcerr = sigma.getError()*resscale;
        hresmc->SetBinContent(i,sigmamc);
        hresmc->SetBinError(i,sigmamcerr);
        hmmc->SetBinContent(i,m0.getVal());
        hmmc->SetBinError(i,m0.getError());
        
        if (i==6) {
          TCanvas *cfit = new TCanvas;
          RooPlot *plot = mass.frame(100);
          tmpdset.plotOn(plot);
          zpdf.plotOn(plot);
          plot->SetTitle("");
          plot->Draw();
          if (var.Contains("sqrt")) cfit->SaveAs(TString::Format("%sstdmcfit.eps",name.Data()));
          else cfit->SaveAs(TString::Format("%s%smcfit.eps",name.Data(),var.Data()));
          new TCanvas;
        }
        
      }    
      
      tmpdset.reset();
      Int_t nwed = tdata->Draw(var, cut);
      double *valsd = tdata->GetV1();
      double *weightsd = tdata->GetW();
      for (UInt_t j=0; j<nwed; ++j) {
        vset->setVal(valsd[j]);
        if (vset->getVal()>65. && vset->getVal()<115)
          tmpdset.add(*vset,weightsd[j]);
      }
      double sigmad = -99.;
      double sigmaderr = -99.;      
      if (tmpdset.sumEntries()>10) {
        zpdf.fitTo(tmpdset,Strategy(1),NumCPU(16));
        sigmad = TMath::Abs(sigma.getVal())*resscale;
        sigmaderr = sigma.getError()*resscale;
        hresd->SetBinContent(i,sigmad);
        hresd->SetBinError(i,sigmaderr);      
        hmd->SetBinContent(i,m0.getVal());
        hmd->SetBinError(i,m0.getError());        
        
        if (i==6) {
          TCanvas *cfit = new TCanvas;
          RooPlot *plot = mass.frame(100);
          tmpdset.plotOn(plot);
          zpdf.plotOn(plot);
          plot->SetTitle("");
          plot->Draw();
          if (var.Contains("sqrt")) cfit->SaveAs(TString::Format("%sstddatafit.eps",name.Data()));
          else cfit->SaveAs(TString::Format("%s%sdatafit.eps",name.Data(),var.Data()));
          new TCanvas;
          
          printf("%s\n",TString(cut).Data());
          
        }
      }    
      if (sigmamc>0.0 && sigmad>0.0) {
        double smear = sqrt(sigmad*sigmad-sigmamc*sigmamc);
        double smearerr = (1.0/smear)*sqrt(sigmad*sigmad*sigmaderr*sigmaderr + sigmamc*sigmamc*sigmamcerr*sigmamcerr);
        if (smear>0.0 && smear<100.0) {
          hress->SetBinContent(i,smear*smearscale);
          hress->SetBinError(i,smearerr*smearscale);
        }
      }
}

void cbfits(TString name, TCut cut, TTree *tmc, TTree *tdata, TString xvar1, TString xvar2, int nbins, double xmin, double xmax, TString lform, TString labelx, double startsigma=1.0, bool dodcor=false) {
  const double zmass = 91.1876;
  const double zwidth = 2.4952;
  
  RooRealVar mass("mass","m_{ee}",65.0,115.0,"GeV");
  mass.setBins(10e3,"cache");
  
  //RooRealVar m0("m0","",91.0,85.0,95.0);
  RooRealVar m0("m0","",0.0,-10.0,10.0);
  m0.removeRange();
  
  RooRealVar sigma("sigma","",startsigma,0.3,5.5);
  //sigma.removeRange();
  
  RooRealVar alpha("alpha","",1.0,-1000.0,1000.0);
  alpha.removeRange();
  
  RooRealVar ncb("ncb","",10.0,-1000.0,1000.0);
  ncb.removeRange();
  
  RooCBShape cbpdf("cbpdf","",mass,m0,sigma,alpha,ncb);
  
  RooBreitWigner bwpdf("bwpdf","",mass,RooConst(zmass),RooConst(zwidth));
  
  RooFFTConvPdf zpdf("zpdf","",mass,bwpdf,cbpdf);   

  
  if (1) {
  

    
    new TCanvas;
    //TCut zcutr = evtcut && "(ph1.elept>25.0 && ph2.elept>25.0 && massmvacorele>60 && massmvacorele<120)";
    //TCut zcutr = evtcut && "(ph1.elept>25.0 && ph2.elept>25.0 && ph1.isbarrel && ph2.isbarrel && massmvacorele>60 && massmvacorele<120)";
    //TCut zcutr = evtcut && "(ph1.elept>25.0 && ph2.elept>25.0 && abs(ph1.sceta)<0.8 && abs(ph2.sceta)<0.8 && massmvacorele>60 && massmvacorele<120)";
    //TCut zcutr = evtcut && "(ph1.elept>25.0 && ph2.elept>25.0 && abs(ph1.sceta)>1.25 && abs(ph2.sceta)>1.25 && ph1.isbarrel && ph2.isbarrel && massmvacorele>60 && massmvacorele<120)";
    //TCut zcutr = evtcut && "(ph1.elept>25.0 && ph2.elept>25.0 && abs(ph1.sceta)<0.8 && abs(ph2.sceta)>1.3 && ph1.emvacorerr/ph1.emvacor<0.015 && ph2.isbarrel && massmvacorele>60 && massmvacorele<120)";
    //TCut zcutr = evtcut && "(ph1.elept>25.0 && ph2.elept>25.0 && abs(ph1.sceta)<1.5 && massmvacorele>60 && massmvacorele<120)";
//     TH2F *herrscatter = new TH2F("herrscatter","",50,0.0,0.05,50,0.0,0.05);
//     hmcele->Draw("ph1.emvacorerr/ph1.emvacor:ph2.emvacorerr/ph2.emvacor>>herrscatter",zcutr);
//     herrscatter->Draw("COL");
//     return;
    
    Int_t nztotmc = tmc->Draw("mass",cut);
    Int_t nztotdat = tdata->Draw("mass",cut);
    
    tmc->SetEstimate(nztotmc);
    tdata->SetEstimate(nztotdat);
    
    RooDataSet tmpdset("tmpdset","",mass);
    RooRealVar *vset = (RooRealVar*)tmpdset.get()->find(mass.GetName());
    
    RooArgSet *pset = zpdf.getParameters(tmpdset);
    RooArgSet *psetinit = (RooArgSet*)pset->snapshot(kTRUE);  
    
    
    
//     const int numresbins = 25;
//     const float resmax = 0.05;
//     //const double resscale = 2.0/sqrt(2.0)/zmass;
//     //const double resscale =  2.0/sqrt(2.0)/zmass;     
//     const double resscale = 1.0/zmass;
    
    TH1F *hresmc = new TH1F("hresmc","",nbins,xmin,xmax);
    TH1F *hresd = new TH1F("hresd","",nbins,xmin,xmax);
    TH1F *hress = new TH1F("hress","",nbins,xmin,xmax);
    TH1F *hmmc = new TH1F("hmmc","",nbins,xmin,xmax);
    TH1F *hmd = new TH1F("hmd","",nbins,xmin,xmax);
    
    TH1F *hrescormc = new TH1F("hrescormc","",nbins,xmin,xmax);
    TH1F *hrescord = new TH1F("hrescord","",nbins,xmin,xmax);
    TH1F *hrescors = new TH1F("hrescors","",nbins,xmin,xmax);
    TH1F *hmcormc = new TH1F("hmcormc","",nbins,xmin,xmax);
    TH1F *hmcord = new TH1F("hmcord","",nbins,xmin,xmax);
    
    TH1F *hresdcormc = new TH1F("hresdcormc","",nbins,xmin,xmax);
    TH1F *hresdcord = new TH1F("hresdcord","",nbins,xmin,xmax);
    TH1F *hresdcors = new TH1F("hresdcors","",nbins,xmin,xmax);
    TH1F *hmdcormc = new TH1F("hmdcormc","",nbins,xmin,xmax);
    TH1F *hmdcord = new TH1F("hmdcord","",nbins,xmin,xmax);
    
    for (int i=1; i<(hresmc->GetNbinsX()+1); ++i) {
      float siglow = hresmc->GetBinLowEdge(i);
      float sighigh = hresmc->GetBinLowEdge(i+1);
      TCut sigcut1(TString::Format("(%s>=%f && %s<%f)",xvar1.Data(),siglow,xvar1.Data(),sighigh));
      //TCut sigcut2(TString::Format("(%s>=%f && %s<%f)",xvar2.Data(),siglow,xvar2.Data(),sighigh));
      //TCut sigcut(TString::Format("(ph1.emvacorerr/ph1.emvacor>=%f && ph2.emvacorerr/ph2.emvacor>=%f && ph1.emvacorerr/ph1.emvacor<%f && ph2.emvacorerr/ph2.emvacor<%f)",siglow,siglow,sighigh,sighigh));
      //TCut sigcut(TString::Format("(ph1.emvacorerr/ph1.emvacor>=%f && ph2.emvacorerr/ph2.emvacor>=%f && ph1.emvacorerr/ph1.emvacor<%f && ph2.emvacorerr/ph2.emvacor<%f)",siglow,siglow,sighigh,sighigh));
      //TCut sigcut(TString::Format("(abs(ph1.sceta)>=%f && abs(ph2.sceta)>=%f && abs(ph1.sceta)<%f && abs(ph2.sceta)<%f)",siglow,siglow,sighigh,sighigh));
      //TCut sigcut(TString::Format("(ph2.emvacorerr/ph2.emvacor>=%f && ph2.emvacorerr/ph2.emvacor<%f)",siglow,sighigh));
      //TCut sigcut(TString::Format("(abs(ph2.sceta)>=%f && abs(ph2.sceta)<%f)",siglow,sighigh));
      TCut zcutrsig = cut && sigcut1;
      //fillcbhists(name,"sqrt(2.0*ph1.sce*ph2.sce*(1.0-costhetaele))",zcutrsig,zpdf,mass,sigma,m0,vset,i,hresmc,hresd,hress,hmmc,hmd, tmc,tdata,tmpdset);
      //fillcbhists(name,"massmvacorele",zcutrsig,zpdf,mass,sigma,m0,vset, i,hrescormc,hrescord,hrescors,hmcormc,hmcord,tmc,tdata,tmpdset);
      fillcbhists(name,"mass",zcutrsig,zpdf,mass,sigma,m0,vset,i,hresmc,hresd,hress,hmmc,hmd, tmc,tdata,tmpdset);
      if (dodcor) fillcbhists(name,"massmvacoreledcor",zcutrsig,zpdf,mass,sigma,m0, vset, i,hresdcormc,hresdcord,hresdcors,hmdcormc,hmdcord,tmc,tdata,tmpdset);

      
    }
    //new TCanvas;
    TF1 *fdiagres = new TF1("fdiag",lform,xmin,xmax);
    fdiagres->SetLineColor(8);
    //fdiagres->Draw("SAME"); 
    
    TF1 *fflat = new TF1("fdiag","0",xmin,xmax);
    fflat->SetLineColor(8);
    
    
    TCanvas *cresmc = new TCanvas;
    hresmc->SetMinimum(0.0);
    hresmc->SetMaximum(1.6*hresmc->GetMaximum());
    hresmc->SetLineColor(kMagenta);
    hresmc->SetMarkerColor(kMagenta);
    hrescormc->SetLineColor(kBlue);
    hrescormc->SetMarkerColor(kBlue);
    
    hresd->SetMinimum(0.0);
    hresd->SetMaximum(1.6*hresd->GetMaximum());    
    hrescord->SetLineColor(kRed);
    hrescord->SetMarkerColor(kRed);
    hresdcord->SetLineColor(8);
    hresdcord->SetMarkerColor(8);    
    hresd->GetXaxis()->SetTitle(labelx);
    hresd->GetYaxis()->SetTitle("#sigma_{CB} (GeV)");    
    
    hresmc->GetXaxis()->SetTitle(labelx);
    hresmc->GetYaxis()->SetTitle("#sigma_{CB} (GeV)");
    hresmc->SetMinimum(0.0);
    hresd->Draw("E");
    //hrescord->Draw("ESAME");
    if (dodcor) hresdcord->Draw("ESAME");
    fdiagres->Draw("SAME");    
    hresmc->Draw("ESAME");
    //hrescormc->Draw("ESAME");
    fdiagres->Draw("SAME");
    
    TLegend *lege = new TLegend(0.45,0.6,0.95,0.9);
    lege->AddEntry(hresd,"Z->ee Data","LPE");
    //lege->AddEntry(hrescord,"Z->ee Data (Raw+Regression)","LPE");
    if (dodcor) lege->AddEntry(hresdcord,"Z->ee Data (Raw+Regression+Res)","LPE");    
    lege->AddEntry(hresmc,"Z->ee MC","LPE");
    //lege->AddEntry(hrescormc,"Z->ee MC (Raw+Regression)","LPE");
    lege->SetBorderSize(0);
    lege->SetFillStyle(0);
    lege->Draw();   
    cresmc->SaveAs(TString::Format("%sresmc.eps",name.Data()));
    
    //hresmc->Draw("E");
    TCanvas *cresd = new TCanvas;

    hresd->Draw("E");
    //hrescord->Draw("ESAME");
    if (dodcor) hresdcord->Draw("ESAME");
    fdiagres->Draw("SAME");

    TLegend *legd = new TLegend(0.45,0.6,0.95,0.9);
    legd->AddEntry(hresd,"Z->ee Data","LPE");
    //legd->AddEntry(hrescord,"Z->ee Data","LPE");
    if (dodcor)legd->AddEntry(hresdcord,"Z->ee Data (Raw+Regression+Res)","LPE");
    legd->SetBorderSize(0);
    legd->SetFillStyle(0);
    legd->Draw();      
    cresd->SaveAs(TString::Format("%sresd.eps",name.Data()));
    
    TCanvas *cress = new TCanvas;
//    hress
    hrescors->SetLineColor(kRed);
    hrescors->SetMarkerColor(kRed);
    hresdcors->SetLineColor(8);
    hresdcors->SetMarkerColor(8);
    
    hress->GetXaxis()->SetTitle(labelx);
    hress->GetYaxis()->SetTitle("Single Photon Relative Smearing");
    hress->SetMinimum(0.0);
    hress->SetMaximum(1.6*hress->GetMaximum());
    hress->Draw("E");
    //hrescors->Draw("ESAME");
    //hresdcors->Draw("ESAME");
    legd->Draw();
    cress->SaveAs(TString::Format("%sress.eps",name.Data()));
    
    
    TCanvas *cmmc = new TCanvas;
    hmmc->SetMinimum(0.0);
    hmmc->SetMaximum(1.6*hmmc->GetMaximum());
    hmmc->SetLineColor(kMagenta);
    hmmc->SetMarkerColor(kMagenta);
    hmcormc->SetLineColor(kBlue);
    hmcormc->SetMarkerColor(kBlue);
    
    hmd->SetMinimum(-3.0);
    hmd->SetMaximum(7.0);    
    hmcord->SetLineColor(kRed);
    hmcord->SetMarkerColor(kRed);
    hmdcord->SetLineColor(8);
    hmdcord->SetMarkerColor(8);    
    hmd->GetXaxis()->SetTitle(labelx);
    hmd->GetYaxis()->SetTitle("#Delta m (GeV)");    
    
    hmmc->GetXaxis()->SetTitle(labelx);
    hmmc->GetYaxis()->SetTitle("#Delta m (GeV)");
    hmmc->SetMinimum(0.0);
    hmd->Draw("E");
    //hmcord->Draw("ESAME");
    //hmdcord->Draw("ESAME");
//    fdiagm->Draw("SAME");    
    hmmc->Draw("ESAME");
    //hmcormc->Draw("ESAME");
    fflat->Draw("SAME");
    
    lege->Draw();   
    cmmc->SaveAs(TString::Format("%smmc.eps",name.Data()));

    
//     hresd->Draw("E");
//     hresmc->Draw("ESame");
   
//     cresd->SaveAs("rescompare.eps");
//     
//     TCanvas *cress = new TCanvas;
//     hress->SetMinimum(0.0);
//     hress->SetMaximum(0.03);
//     hress->Draw("E");
//     cress->SaveAs("ressmear.eps");
    //return;
  }  
}


void plotstack(TString name, TTree *tdata, TTree *tmc, TString draw, TCut cut, TString xlabel, TString ylabel, Double_t xmin, Double_t xmax, UInt_t nbins=100, bool rescale=false, bool dolog=false, bool legleft=false, TString drawalt1="", TString drawalt2="", TString altlabel="", TCut altcut1=TCut(), TCut altcut2=TCut()) {
  

  
//   TCut promptpromptcut = "proc>3 && proc<11 && (ispromptgen1&&ispromptgen2)";
//   TCut promptfakecut =   "proc>3 && proc<11 && (ispromptgen1&&!ispromptgen2)";
//   TCut fakefakecut =   "proc>3 && proc<11 && (!ispromptgen1)";
//   TCut dycut = "proc==11";
//   TCut datacut = "proc==12";

//   TCut promptpromptcut = "proc>3 && proc<11 && (ispromptgen1&&ispromptgen2)";
//   TCut promptfakecut =   "proc>3 && proc<=11 && ( (ispromptgen1&&!ispromptgen2) || (!ispromptgen1&&ispromptgen2))";
//   TCut fakefakecut =   "proc>3 && proc<=11 && (!ispromptgen1&&!ispromptgen2)";
//   TCut dycut = "proc==11 && (ispromptgen1&&ispromptgen2)";
//   TCut datacut = "proc==12";  
  
  //TCut weightc = "(weight)"; 
  

  TH1D *hdata = new TH1D(name+"data","",nbins,xmin,xmax);
 
  
  TH1D *hbkg = new TH1D(name+"bkg","",nbins,xmin,xmax);
  TH1D *hbkgalt1 = new TH1D(name+"bkgalt1","",nbins,xmin,xmax);
  TH1D *hbkgalt2 = new TH1D(name+"bkgalt2","",nbins,xmin,xmax);
  

  hbkg->Sumw2();
  hbkgalt1->Sumw2();
  hbkgalt2->Sumw2(); 
  hdata->Sumw2();
   
  tdata->Draw(draw+">>"+name+"data",cut,"goff");

  tmc->Draw(draw+">>"+name+"bkg",cut,"goff");  
  
  
  bool dosig = false;
  //bool dosig = false;
  


  
  
  bool doalt = drawalt1!="" && drawalt2!="";
  if (doalt) {    
    TCut acut1 = cut;
    TCut acut2 = cut;
    
    if (altcut1!=TCut()) acut1 = altcut1;
    if (altcut2!=TCut()) acut2 = altcut2;
    
    tmc->Draw(drawalt1+">>"+name+"bkgalt1",acut1,"goff");  
    tmc->Draw(drawalt2+">>"+name+"bkgalt2",acut2,"goff");  
  }
  
  Double_t mcsum = hbkg->GetSumOfWeights();
  Double_t datasum = hdata->GetSumOfWeights();
  Double_t scale = datasum/mcsum;
  
  if (rescale) {
    hbkg->Scale(scale);
    hbkgalt1->Scale(scale);    
    hbkgalt2->Scale(scale);    
  }
 
  //if (dosig) hsig->Scale(datasum/hsig->GetSumOfWeights());
 
  UInt_t numbins = hbkg->GetXaxis()->GetNbins();
  TGraphAsymmErrors *bkgalt = new TGraphAsymmErrors;
  TGraphAsymmErrors *bkgaltratio = new TGraphAsymmErrors;
  //bkgalt->SetFillColor(k
  for (int i=1; i<(numbins+1); ++i) {
    double nominal = hbkg->GetBinContent(i); 
    double alt1 = hbkgalt1->GetBinContent(i);
    double alt2 = hbkgalt2->GetBinContent(i);
    
    double xerr = 0.5*hbkg->GetXaxis()->GetBinWidth(i);
    double uperr = TMath::Max(0.0,TMath::Max(alt1,alt2)-nominal );
    double downerr = TMath::Min(0.0,TMath::Min(alt1,alt2)-nominal );
  
    double staterr = hbkg->GetBinError(i);
    
//     uperr = sqrt(uperr*uperr+staterr*staterr);
//     downerr = -sqrt(downerr*downerr+staterr*staterr);
    
    //printf("nominal=%5f, uperr=%5f, downerr=%5f, alt1=%5f, alt2=%5f\n",nominal,uperr,downerr,alt1,alt2);
    
    bkgalt->SetPoint(i-1,hbkg->GetXaxis()->GetBinCenter(i),nominal);
    bkgaltratio->SetPoint(i-1,hbkg->GetXaxis()->GetBinCenter(i),1.0);    
    
    bkgalt->SetPointError(i-1,xerr,xerr,TMath::Abs(downerr),TMath::Abs(uperr));
    
    double up = TMath::Max(nominal,TMath::Max(alt1,alt2));
    double down = TMath::Min(nominal,TMath::Min(alt1,alt2));
    
    bkgaltratio->SetPointError(i-1,xerr,xerr,1.0-down/nominal,up/nominal-1.0);
    
    //bkgaltratio->SetPointError(i-1,0,0,TMath::Abs(downerr)/nominal,TMath::Abs(uperr)/nominal);

//     bkgalt->SetPointError(i-1,0,0,5000.0,5000.0);
    
    bkgalt->SetFillStyle(3244);
    bkgalt->SetFillColor(kRed);
    bkgaltratio->SetFillStyle(3244);
    bkgaltratio->SetFillColor(kRed);    

    
  }
  
  

  hbkg->SetFillStyle(1001);
  hbkg->SetFillColor(kBlue-7);


  
  if (hdata->GetMaximum()>hbkg->GetMaximum()) {
    hbkg->SetMaximum(1.1*hbkg->GetMaximum());
  }
     
  
  TH1D *hratio = new TH1D( (*hdata)/(*hbkg) );
  hratio->GetYaxis()->SetTitle("Data/MC");
  
  TCanvas *c1 = new TCanvas;

  bool doratio=true;
  
  TPad*c1npad1 = 0;
  TPad*c1npad2 = 0;
  if (doratio) {
    //TCanvas* c1 = new TCanvas("1","Plots1",5,30,900,600); // create new canvas
    c1->Clear();
    c1npad1 = new TPad("c1npad1","",0,0.2,1,1); // create pad for the plots
    c1npad2 = new TPad("c1npad2","",0,0,1,0.2); // create pad for the ratio
    c1->cd(); // set working area to main area of canvas
    c1npad1->Draw(); // draw larger pad in main area
    c1npad1->cd(); // change working area to inside main pad
  }
  
  hbkg->Draw("HIST");
  hbkg->GetXaxis()->SetTitle(xlabel);
  hbkg->GetYaxis()->SetTitle(ylabel);    
  if (doalt) bkgalt->Draw("2 SAME");
  
  hdata->Draw("PESAME");
  
  
  
  TLegend *leg = 0;
  if (legleft) leg = new TLegend(0.20,0.7,0.47,0.9);
  else leg = new TLegend(0.65,0.7,0.92,0.9);
  
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);  
  
  
  leg->AddEntry(hdata,"Data","LPE");
  leg->AddEntry(hbkg,"MC","LF");
  
  if (doalt) leg->AddEntry(bkgalt,altlabel,"F");
  
  leg->Draw();
    
  if (doratio) {
    c1->cd(); // go back to working area of whole canvas
    c1npad2->Draw(); // draw smaller pad in main canvas
    c1npad2->cd(); // set working area to smaller pad
    hratio->GetYaxis()->SetTitle("Ratio");
    hratio->GetYaxis()->SetTitleSize(0.12);
    hratio->GetYaxis()->SetTitleOffset(0.3);
    hratio->GetYaxis()->SetLabelSize(0.1);
    hratio->GetXaxis()->SetLabelSize(0.1);
    //hratio->GetYaxis()->SetRangeUser(0.8,1.2);
    hratio->GetYaxis()->SetNdivisions(4,kTRUE);
    hratio->SetTickLength(0.01,"Y");
    //hratio->SetMarkerStyle(7);  
    if (hratio->GetMaximum()>1.5) hratio->SetMaximum(1.5);
    if (hratio->GetMinimum()<0.8) hratio->SetMinimum(0.8);
    hratio->Draw("PE");
    if (doalt) {
      bkgaltratio->Draw("2 Same");
      hratio->Draw("PESame");
    }
    TLine *line = new TLine(xmin,1.,xmax,1.);
    line->Draw();
    
    if (draw=="mva") {
      hratio->SetMaximum(1.15);
      hratio->SetMinimum(0.85);
      
      TLine *line1 = new TLine(0.05,hratio->GetMinimum(),0.05,hratio->GetMaximum());
      TLine *line2 = new TLine(0.55,hratio->GetMinimum(),0.55,hratio->GetMaximum());
      TLine *line3 = new TLine(0.72,hratio->GetMinimum(),0.72,hratio->GetMaximum());
      TLine *line4 = new TLine(0.89,hratio->GetMinimum(),0.89,hratio->GetMaximum());
      
      line1->SetLineWidth(2); line1->SetLineStyle(9);
      line2->SetLineWidth(2); line2->SetLineStyle(9);
      line3->SetLineWidth(2); line3->SetLineStyle(9);
      line4->SetLineWidth(2); line4->SetLineStyle(9);
      
      line1->Draw();
      line2->Draw();
      line3->Draw();
      line4->Draw(); 
    }    
    
    c1->cd(); // go back to main canvas
  }
  
//   else {
//     hstack->Draw();
//     hstack->GetXaxis()->SetTitle(xlabel);
//     hstack->GetYaxis()->SetTitle(ylabel);    
//     if (doalt) bkgaltratio->Draw("4 SAME");
//     hdata->Draw("PESAME");
//   }
  

  
  //TLatex lat(1.5,0.9*hstack->GetMaximum(),"#scale[0.7]{#splitline{CMS preliminary}{#sqrt{s} = 7 TeV L = 4.76 fb^{-1}}}");
  TLatex lat(1.5,0.9*hbkg->GetMaximum(),"#scale[0.7]{#splitline{CMS preliminary}{#sqrt{s} = 8 TeV L = 3.1 fb^{-1}}}");
  lat.Draw();  
  
  printf("data entries = %5f\n",hdata->GetSumOfWeights());
  
  if (dolog) {
    if (doratio) {
      c1npad1->SetLogy();
    }
    else {
      c1->SetLogy();
    }
    hbkg->SetMinimum(0.5);
  }
  
  c1->SaveAs(name+".eps");
  c1->SaveAs(name+".pdf");
  
  //hpromptprompt->Draw();
  
}

void plotstack(TString name, TTree *tree, TString draw, TCut cut, TString xlabel, TString ylabel, Double_t xmin, Double_t xmax, UInt_t nbins=100, bool splitsingle = false, bool rescale=false, bool dolog=false, bool legleft=false, TString drawalt1="", TString drawalt2="", TString altlabel="", TCut altcut1=TCut(), TCut altcut2=TCut()) {
  
  TCut promptpromptcut = "procidx>3 && procidx<12 && (ph1.ispromptgen && ph2.ispromptgen)";
  TCut promptfakecut =   "procidx>3 && procidx<12 && ( (ph1.ispromptgen&&!ph2.ispromptgen) || (!ph1.ispromptgen&&ph2.ispromptgen))";
  TCut fakefakecut =   "procidx>3 && procidx<12 && (!ph1.ispromptgen&&!ph2.ispromptgen)";
  TCut dycut = "procidx==12";
  TCut bkgcut = "procidx>3 && procidx<13";
  TCut sigcut = "procidx<=3";
  TCut datacut = "procidx==13";

  if (splitsingle) {
    promptpromptcut = "procidx>3 && procidx<12 && (ph1.ispromptgen)";
    promptfakecut =   "0==1";
    fakefakecut =   "procidx>3 && procidx<12 && (!ph1.ispromptgen)";       
  }
  
//   TCut promptpromptcut = "proc>3 && proc<11 && (ispromptgen1&&ispromptgen2)";
//   TCut promptfakecut =   "proc>3 && proc<11 && (ispromptgen1&&!ispromptgen2)";
//   TCut fakefakecut =   "proc>3 && proc<11 && (!ispromptgen1)";
//   TCut dycut = "proc==11";
//   TCut datacut = "proc==12";

//   TCut promptpromptcut = "proc>3 && proc<11 && (ispromptgen1&&ispromptgen2)";
//   TCut promptfakecut =   "proc>3 && proc<=11 && ( (ispromptgen1&&!ispromptgen2) || (!ispromptgen1&&ispromptgen2))";
//   TCut fakefakecut =   "proc>3 && proc<=11 && (!ispromptgen1&&!ispromptgen2)";
//   TCut dycut = "proc==11 && (ispromptgen1&&ispromptgen2)";
//   TCut datacut = "proc==12";  
  
  //TCut weightc = "(weight)"; 
  
  THStack *hstack = new THStack(name,"");
  
  TH1D *hpromptprompt = new TH1D(name+"promptprompt","",nbins,xmin,xmax);
  TH1D *hpromptfake = new TH1D(name+"promptfake","",nbins,xmin,xmax);
  TH1D *hfakefake = new TH1D(name+"fakefake","",nbins,xmin,xmax);
  TH1D *hdy = new TH1D(name+"dy","",nbins,xmin,xmax);
  TH1D *hdata = new TH1D(name+"data","",nbins,xmin,xmax);
 
  TH1D *hsig = new TH1D(name+"sig","",nbins,xmin,xmax);

  
  TH1D *hbkg = new TH1D(name+"bkg","",nbins,xmin,xmax);
  TH1D *hbkgalt1 = new TH1D(name+"bkgalt1","",nbins,xmin,xmax);
  TH1D *hbkgalt2 = new TH1D(name+"bkgalt2","",nbins,xmin,xmax);
  
  hpromptprompt->Sumw2();
  hpromptfake->Sumw2();
  hfakefake->Sumw2();
  hdy->Sumw2();
  hbkg->Sumw2();
  hsig->Sumw2();
  hbkgalt1->Sumw2();
  hbkgalt2->Sumw2(); 
  hdata->Sumw2();
  
//   tree->Draw(draw+">>"+name+"promptprompt",weightc*(promptpromptcut&&cut),"goff");
//   tree->Draw(draw+">>"+name+"promptfake",weightc*(promptfakecut&&cut),"goff");
//   tree->Draw(draw+">>"+name+"fakefake",weightc*(fakefakecut&&cut),"goff");
//   tree->Draw(draw+">>"+name+"dy",weightc*(dycut&&cut),"goff");
//   tree->Draw(draw+">>"+name+"data",(datacut&&cut),"goff");
// 
//   tree->Draw(draw+">>"+name+"bkg",weightc*(bkgcut&&cut),"goff");
  
  
  tree->Draw(draw+">>"+name+"promptprompt",promptpromptcut*cut,"goff");
  tree->Draw(draw+">>"+name+"promptfake",promptfakecut*cut,"goff");
  tree->Draw(draw+">>"+name+"fakefake",fakefakecut*cut,"goff");
  tree->Draw(draw+">>"+name+"dy",dycut*cut,"goff");
  tree->Draw(draw+">>"+name+"data",datacut*cut,"goff");

  tree->Draw(draw+">>"+name+"bkg",bkgcut*cut,"goff");  
  
  
  bool dosig = false;
  //bool dosig = false;
  
  if (dosig) {
    tree->Draw(draw+">>"+name+"sig",sigcut*cut,"goff");  
  }

  
  
  bool doalt = drawalt1!="" && drawalt2!="";
  if (doalt) {    
    TCut acut1 = bkgcut*cut;
    TCut acut2 = bkgcut*cut;
    
    if (altcut1!=TCut()) acut1 = bkgcut*altcut1;
    if (altcut2!=TCut()) acut2 = bkgcut*altcut2;
    
    tree->Draw(drawalt1+">>"+name+"bkgalt1",acut1,"goff");  
    tree->Draw(drawalt2+">>"+name+"bkgalt2",acut2,"goff");  
  }
  
  Double_t mcsum = hpromptprompt->GetSumOfWeights() + hpromptfake->GetSumOfWeights() + hfakefake->GetSumOfWeights() + hdy->GetSumOfWeights();
  Double_t datasum = hdata->GetSumOfWeights();
  Double_t scale = datasum/mcsum;
  
  if (rescale) {
    hpromptprompt->Scale(scale);
    hpromptfake->Scale(scale);
    hfakefake->Scale(scale);
    hdy->Scale(scale);
    hbkg->Scale(scale);
    hbkgalt1->Scale(scale);    
    hbkgalt2->Scale(scale);    
  }
 
  //if (dosig) hsig->Scale(datasum/hsig->GetSumOfWeights());
  if (dosig) hsig->Scale(100.0);
 
  UInt_t numbins = hbkg->GetXaxis()->GetNbins();
  TGraphAsymmErrors *bkgalt = new TGraphAsymmErrors;
  TGraphAsymmErrors *bkgaltratio = new TGraphAsymmErrors;
  //bkgalt->SetFillColor(k
  for (int i=1; i<(numbins+1); ++i) {
    double nominal = hbkg->GetBinContent(i); 
    double alt1 = hbkgalt1->GetBinContent(i);
    double alt2 = hbkgalt2->GetBinContent(i);
    
    double xerr = 0.5*hbkg->GetXaxis()->GetBinWidth(i);
    double uperr = TMath::Max(0.0,TMath::Max(alt1,alt2)-nominal );
    double downerr = TMath::Min(0.0,TMath::Min(alt1,alt2)-nominal );
  
    double staterr = hbkg->GetBinError(i);
    
//     uperr = sqrt(uperr*uperr+staterr*staterr);
//     downerr = -sqrt(downerr*downerr+staterr*staterr);
    
    //printf("nominal=%5f, uperr=%5f, downerr=%5f, alt1=%5f, alt2=%5f\n",nominal,uperr,downerr,alt1,alt2);
    
    bkgalt->SetPoint(i-1,hbkg->GetXaxis()->GetBinCenter(i),nominal);
    bkgaltratio->SetPoint(i-1,hbkg->GetXaxis()->GetBinCenter(i),1.0);    
    
    bkgalt->SetPointError(i-1,xerr,xerr,TMath::Abs(downerr),TMath::Abs(uperr));
    
    double up = TMath::Max(nominal,TMath::Max(alt1,alt2));
    double down = TMath::Min(nominal,TMath::Min(alt1,alt2));
    
    bkgaltratio->SetPointError(i-1,xerr,xerr,1.0-down/nominal,up/nominal-1.0);
    
    //bkgaltratio->SetPointError(i-1,0,0,TMath::Abs(downerr)/nominal,TMath::Abs(uperr)/nominal);

//     bkgalt->SetPointError(i-1,0,0,5000.0,5000.0);
    
    bkgalt->SetFillStyle(3244);
    bkgalt->SetFillColor(kRed);
    bkgaltratio->SetFillStyle(3244);
    bkgaltratio->SetFillColor(kRed);    

    
  }
  
  
  hpromptprompt->SetFillColor(8);
  hpromptfake->SetFillColor(kOrange-2);
  hfakefake->SetFillColor(kRed);
  hdy->SetFillColor(kBlue-7);
  
  hpromptprompt->SetFillStyle(1001);
  hpromptfake->SetFillStyle(1001);
  hfakefake->SetFillStyle(1001);
  hdy->SetFillStyle(1001);
  

  hstack->Add(hpromptprompt);
  hstack->Add(hpromptfake);
  hstack->Add(hfakefake);
  hstack->Add(hdy);
  
  if (hdata->GetMaximum()>hstack->GetMaximum()) {
    hstack->SetMaximum(1.1*hdata->GetMaximum());
  }
     
  
  TH1D *hratio = new TH1D( (*hdata)/(*hbkg) );
  hratio->GetYaxis()->SetTitle("Data/MC");
  
  TCanvas *c1 = new TCanvas;

  bool doratio=true;
  
  TPad*c1npad1 = 0;
  TPad*c1npad2 = 0;
  if (doratio) {
    //TCanvas* c1 = new TCanvas("1","Plots1",5,30,900,600); // create new canvas
    c1->Clear();
    c1npad1 = new TPad("c1npad1","",0,0.2,1,1); // create pad for the plots
    c1npad2 = new TPad("c1npad2","",0,0,1,0.2); // create pad for the ratio
    c1->cd(); // set working area to main area of canvas
    c1npad1->Draw(); // draw larger pad in main area
    c1npad1->cd(); // change working area to inside main pad
  }
  
  hstack->Draw("HIST");
  hstack->GetXaxis()->SetTitle(xlabel);
  hstack->GetYaxis()->SetTitle(ylabel);    
  if (doalt) bkgalt->Draw("2 SAME");
  hdata->Draw("PESAME");
  
  if (draw=="mva") {
    TLine *line1 = new TLine(0.05,hstack->GetMinimum(),0.05,hstack->GetMaximum());
    TLine *line2 = new TLine(0.55,hstack->GetMinimum(),0.55,hstack->GetMaximum());
    TLine *line3 = new TLine(0.72,hstack->GetMinimum(),0.72,hstack->GetMaximum());
    TLine *line4 = new TLine(0.89,hstack->GetMinimum(),0.89,hstack->GetMaximum());
    
    line1->SetLineWidth(2); line1->SetLineStyle(9);
    line2->SetLineWidth(2); line2->SetLineStyle(9);
    line3->SetLineWidth(2); line3->SetLineStyle(9);
    line4->SetLineWidth(2); line4->SetLineStyle(9);
    
    line1->Draw();
    line2->Draw();
    line3->Draw();
    line4->Draw(); 
  }  
  
  hsig->SetLineColor(kBlue);
  if (dosig) hsig->Draw("HISTSAME");
   
  
  TLegend *leg = 0;
  if (legleft) leg = new TLegend(0.20,0.7,0.47,0.9);
  else leg = new TLegend(0.65,0.7,0.92,0.9);
  
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);  
  
  
  leg->AddEntry(hdata,"Data","LPE");
  leg->AddEntry(hdy,"Drell-Yan","LF");
  if (splitsingle) {
    leg->AddEntry(hfakefake,"Fake","LF");
    leg->AddEntry(hpromptprompt,"Prompt","LF");    
  }
  else {
    leg->AddEntry(hfakefake,"Fake-Fake","LF");
    leg->AddEntry(hpromptfake,"Prompt-Fake","LF");
    leg->AddEntry(hpromptprompt,"Prompt-Prompt","LF");
  }
  
  if (doalt) leg->AddEntry(bkgalt,altlabel,"F");
  if (dosig) leg->AddEntry(hsig,"100x m_h=125 GeV","L");
  
  leg->Draw();
    
  if (doratio) {
    c1->cd(); // go back to working area of whole canvas
    c1npad2->Draw(); // draw smaller pad in main canvas
    c1npad2->cd(); // set working area to smaller pad
    hratio->GetYaxis()->SetTitle("Ratio");
    hratio->GetYaxis()->SetTitleSize(0.12);
    hratio->GetYaxis()->SetTitleOffset(0.3);
    hratio->GetYaxis()->SetLabelSize(0.1);
    hratio->GetXaxis()->SetLabelSize(0.1);
    //hratio->GetYaxis()->SetRangeUser(0.8,1.2);
    hratio->GetYaxis()->SetNdivisions(4,kTRUE);
    hratio->SetTickLength(0.01,"Y");
    //hratio->SetMarkerStyle(7);  
    if (hratio->GetMaximum()>1.5) hratio->SetMaximum(1.5);
    if (hratio->GetMinimum()<0.8) hratio->SetMinimum(0.8);
    hratio->Draw("PE");
    if (doalt) {
      bkgaltratio->Draw("2 Same");
      hratio->Draw("PESame");
    }
    TLine *line = new TLine(xmin,1.,xmax,1.);
    line->Draw();
    
    if (draw=="mva") {
      hratio->SetMaximum(1.15);
      hratio->SetMinimum(0.85);
      
      TLine *line1 = new TLine(0.05,hratio->GetMinimum(),0.05,hratio->GetMaximum());
      TLine *line2 = new TLine(0.55,hratio->GetMinimum(),0.55,hratio->GetMaximum());
      TLine *line3 = new TLine(0.72,hratio->GetMinimum(),0.72,hratio->GetMaximum());
      TLine *line4 = new TLine(0.89,hratio->GetMinimum(),0.89,hratio->GetMaximum());
      
      line1->SetLineWidth(2); line1->SetLineStyle(9);
      line2->SetLineWidth(2); line2->SetLineStyle(9);
      line3->SetLineWidth(2); line3->SetLineStyle(9);
      line4->SetLineWidth(2); line4->SetLineStyle(9);
      
      line1->Draw();
      line2->Draw();
      line3->Draw();
      line4->Draw(); 
    }    
    
    c1->cd(); // go back to main canvas
  }
  
//   else {
//     hstack->Draw();
//     hstack->GetXaxis()->SetTitle(xlabel);
//     hstack->GetYaxis()->SetTitle(ylabel);    
//     if (doalt) bkgaltratio->Draw("4 SAME");
//     hdata->Draw("PESAME");
//   }
  

  
  TLatex lat(1.5,0.9*hstack->GetMaximum(),"#scale[0.7]{#splitline{CMS preliminary}{#sqrt{s} = 7 TeV L = 4.76 fb^{-1}}}");
  lat.Draw();  
  
  printf("data entries = %5f\n",hdata->GetSumOfWeights());
  
  if (dolog) {
    if (doratio) {
      c1npad1->SetLogy();
    }
    else {
      c1->SetLogy();
    }
    hstack->SetMinimum(0.5);
  }
  
  c1->SaveAs(name+".eps");
  c1->SaveAs(name+".pdf");
  
  //hpromptprompt->Draw();
  
}

void stackplotsmva() {
  
  ROOT::Cintex::Cintex::Enable();   

  gROOT->Macro("MitStyle.C");
  gStyle->SetErrorX(0); 
  gStyle->SetOptStat(0);
  gROOT->ForceStyle();    
  
  gSystem->cd("./scaleplotJun17");
  
  //-------------set up the weights-------------------
  //pu target
  //TFile *filepuest = new TFile("/scratch/bendavid/root/puweightsNov13/2011_0100_73500.pileup.root","READ");
  TFile *filepuest = new TFile("/scratch/bendavid/root/puweightsJun14/rereco_prompt_2012.json.69000.observed.pileup.root","READ");
  TH1D *hpuest = (TH1D*) filepuest->Get("pileup");

  
  TChain *tree = new TChain("RunLumiSelectionMod/MCProcessSelectionMod/HLTModP/GoodPVFilterMod/PhotonMvaMod/JetPub/JetCorrectionMod/PhotonPairSelectorPresel/PhotonTreeWriterPresel/hPhotonTree");
  //TChain("RunLumiSelectionMod/MCProcessSelectionMod/HLTModP/GoodPVFilterMod/PhotonMvaMod/JetPub/JetCorrectionMod/PhotonPairSelectorPreselInvertEleVeto/PhotonTreeWriterPreselInvertEleVeto/hPhotonTree");
  tree->Add("/scratch/bendavid/cms/hist/hgg-v0/merged/hgg-v0_s12-h120gg-gf-v9_noskim.root");
  tree->Add("/scratch/bendavid/cms/hist/hgg-v0/merged/hgg-v0_s12-h120gg-vbf-v9_noskim.root");
  tree->Add("/scratch/bendavid/cms/hist/hgg-v0/merged/hgg-v0_s12-h120gg-vh-v9_noskim.root");
  tree->Add("/scratch/bendavid/cms/hist/hgg-v0/merged/hgg-v0_s12-h120gg-tt-v9_noskim.root");
  tree->Add("/scratch/bendavid/cms/hist/hgg-v0/merged/hgg-v0_s12-2pibx10_25-v9_noskim.root");
  tree->Add("/scratch/bendavid/cms/hist/hgg-v0/merged/hgg-v0_s12-2pibx25_250-v9_noskim.root");
  tree->Add("/scratch/bendavid/cms/hist/hgg-v0/merged/hgg-v0_s12-2pibx250-v9_noskim.root");
  tree->Add("/scratch/bendavid/cms/hist/hgg-v0/merged/hgg-v0_s12-diphoj-v9_noskim.root");
  tree->Add("/scratch/bendavid/cms/hist/hgg-v0/merged/hgg-v0_s12-qcd-2em3040-v9_noskim.root");
  tree->Add("/scratch/bendavid/cms/hist/hgg-v0/merged/hgg-v0_s12-qcd-2em40-v9_noskim.root");
  tree->Add("/scratch/bendavid/cms/hist/hgg-v0/merged/hgg-v0_s12-pj20_40-2em-v9_noskim.root");
  tree->Add("/scratch/bendavid/cms/hist/hgg-v0/merged/hgg-v0_s12-pj40-2em-v9_noskim.root");
  tree->Add("/scratch/bendavid/cms/hist/hgg-v0/merged/hgg-v0_s12-zllm50-2-v9_noskim.root");
  tree->Add("/scratch/bendavid/cms/hist/hgg-v0/MergedPhoton2012.root");
  
  tree->SetCacheSize(64*1024*1024);
  
  
  
  float xsecs[50];
  xsecs[0] = 21.20*2.25e-03;
  xsecs[1] = 1.632*2.25e-03;
  xsecs[2] = (0.7966+0.4483)*2.25e-03;
  xsecs[3] = 0.147*2.25e-03;
  xsecs[4] = 424.8;
  xsecs[5] = 15.54;
  xsecs[6] = 0.0011805;
  xsecs[7] = 81.16;
  xsecs[8] = 0.000235*5.195e+07;
  xsecs[9] = 0.002175*2.365e+07;
  xsecs[10] = 0.001835*81930.0;
  xsecs[11] = 0.05387*8884.0;
  xsecs[12] = 2950.0;

  
  const double lumi = 3058.;  
  initweights(tree,hpuest,xsecs,lumi);
  puweights[13] = 0;
  xsecweights[13] = 1.0;
   

  TFile *fdata = TFile::Open("/tmp/bendavid/treesJun15/MergedPhoton2012.root");
  TDirectory *ddata = (TDirectory*)fdata->FindObjectAny("PhotonTreeWriterPreselInvertEleVeto");
  TTree *hdata = (TTree*)ddata->Get("hPhotonTree");

  TFile *fmc = TFile::Open("/tmp/bendavid/treesJun15/hgg-v0_s12-zllm50-2-v9_noskim.root");
  TDirectory *dmc = (TDirectory*)fmc->FindObjectAny("PhotonTreeWriterPreselInvertEleVeto");
  TTree *hmc = (TTree*)dmc->Get("hPhotonTree");  
  puweights[15] = getpuweights(fmc,hpuest);
  xsecweights[15] = getweight(fmc,lumi*2950.);
  
  TFile *freg = TFile::Open("/home/bendavid/cms/cmssw/026/CMSSW_5_2_3_patch4/src/MitPhysics/data/gbrv3ph_52x.root");
  GBRForest *forestvareb = (GBRForest*)freg->Get("EBUncertainty");
  GBRForest *forestvaree = (GBRForest*)freg->Get("EEUncertainty");
  std::vector<std::string> *ebvars = (std::vector<std::string>*)freg->Get("varlisteb");
  std::vector<std::string> *eevars = (std::vector<std::string>*)freg->Get("varlistee");
  
  std::vector<std::string> ebvarsnew(*ebvars);
  std::vector<std::string> eevarsnew(*eevars);
  
  for (int i=0; i<ebvarsnew.size(); ++i) {
    printf("ebvarsnew[%i] = \"%s\";\n",i,ebvarsnew.at(i).c_str());
  }

  ebvarsnew[4] = "TMath::Min(1.0,(1.0+ismc*0.0016)*ph.e5x5/ph.scrawe)";  
  //ebvarsnew[13] = "(1.0+ismc*0.0016)*ph.eseed/ph.scrawe";  
  //ebvarsnew[14] = "(1.0+ismc*0045)*ph.e3x3seed/ph.eseed + ismc*0.00";  
  ebvarsnew[14] = "ph.r9*ph.scrawe/ph.eseed";  
  ebvarsnew[15] = "(TMath::Min(1.0,1.0+ismc*0.0016)*ph.e5x5seed/ph.eseed)";   
  ebvarsnew[16] = "(1.0 - ismc*0.108168)*ph.sigietaietaseed + ismc*0.0009133";
  ebvarsnew[17] = "(1.0 - ismc*0.007)*ph.sigiphiphiseed";
  ebvarsnew[19] = "(1.0 + ismc*0.012)*ph.emaxseed/ph.eseed";
  ebvarsnew[20] = "(1.0 + ismc*0.0)*ph.e2ndseed/ph.eseed";
  ebvarsnew[21] = "(1.0-ismc*0.06)*ph.etopseed/ph.eseed";
  ebvarsnew[22] = "(1.0-ismc*0.06)*ph.ebottomseed/ph.eseed";
  ebvarsnew[23] = "(1.0-ismc*0.06)*ph.eleftseed/ph.eseed";
  ebvarsnew[24] = "(1.0-ismc*0.06)*ph.erightseed/ph.eseed";
  ebvarsnew[25] = "(1.0+ismc*0.006)*ph.e2x5maxseed/ph.eseed";
  ebvarsnew[26] = "(1.0+ismc*0.09)*ph.e2x5topseed/ph.eseed";
  ebvarsnew[27] = "(1.0+ismc*0.09)*ph.e2x5bottomseed/ph.eseed";
  ebvarsnew[28] = "(1.0+ismc*0.09)*ph.e2x5leftseed/ph.eseed";
  ebvarsnew[29] = "(1.0+ismc*0.09)*ph.e2x5rightseed/ph.eseed";
 
  eevarsnew[4] = "TMath::Min(1.0,(1.0+ismc*0.0022)*ph.e5x5/ph.scrawe)";  
  eevarsnew[13] = "TMath::Min(1.0,(1.0+ismc*0.0022)*ph.eseed/ph.scrawe)"; 
  eevarsnew[14] = "ph.r9*ph.scrawe/ph.eseed";
  //eevarsnew[14] = "(1.0+ismc*0.0022)*ph.e5x5seed/ph.eseed";  
  //eevarsnew[14] = "(1.0 + ismc*0.007)*ph.e3x3seed/ph.eseed";  
  //eevarsnew[15] = "(1.0+ismc*0.0022)*ph.e5x5seed/ph.eseed";   
  eevarsnew[16] = "(1.0 - ismc*0.0053)*ph.sigietaietaseed + ismc*0.00003";
  eevarsnew[17] = "ph.sigiphiphiseed";
  eevarsnew[19] = "(1.0+ismc*0.005)*ph.emaxseed/ph.eseed";
  eevarsnew[20] = "(1.0+ismc*0.02)*ph.e2ndseed/ph.eseed";
  eevarsnew[21] = "(1.0-ismc*0.04)*ph.etopseed/ph.eseed";
  eevarsnew[22] = "(1.0-ismc*0.04)*ph.ebottomseed/ph.eseed";
  eevarsnew[23] = "(1.0-ismc*0.04)*ph.eleftseed/ph.eseed";
  eevarsnew[24] = "(1.0-ismc*0.04)*ph.erightseed/ph.eseed";
  eevarsnew[25] = "(1.0+ismc*0.0075)*ph.e2x5maxseed/ph.eseed";
  eevarsnew[26] = "(1.0+ismc*0.13)*ph.e2x5topseed/ph.eseed";
  eevarsnew[27] = "(1.0+ismc*0.13)*ph.e2x5bottomseed/ph.eseed";
  eevarsnew[28] = "(1.0+ismc*0.13)*ph.e2x5leftseed/ph.eseed";
  eevarsnew[29] = "(1.0+ismc*0.13)*ph.e2x5rightseed/ph.eseed";  
  
  std::vector<std::string> ebvarsnew1;
  std::vector<std::string> ebvarsnew2;
  
  std::vector<std::string> eevarsnew1;
  std::vector<std::string> eevarsnew2;
  
  for (int i=0; i<ebvarsnew.size(); ++i) {
    TString str1 = ebvarsnew[i];
    TString str2 = ebvarsnew[i];
    str1.ReplaceAll("ph.","ph1.");
    str2.ReplaceAll("ph.","ph2.");
    
    ebvarsnew1.push_back(std::string(str1));
    ebvarsnew2.push_back(std::string(str2));
  }
  
  for (int i=0; i<eevarsnew.size(); ++i) {
    TString str1 = eevarsnew[i];
    TString str2 = eevarsnew[i];
    str1.ReplaceAll("ph.","ph1.");
    str2.ReplaceAll("ph.","ph2.");
    
    eevarsnew1.push_back(std::string(str1));
    eevarsnew2.push_back(std::string(str2));
  }
  

  
  GBRApply<GBRForest> apply;
//   apply.ApplyAsFriend(hmc,forestvareb,ebvarsnew1,"ph1.eerrmodebraw");
//   apply.ApplyAsFriend(hmc,forestvareb,ebvarsnew2,"ph2.eerrmodebraw");
//   apply.ApplyAsFriend(hmc,forestvaree,eevarsnew1,"ph1.eerrmodeeraw");
//   apply.ApplyAsFriend(hmc,forestvaree,eevarsnew2,"ph2.eerrmodeeraw");
  apply.ApplyAsFriend(tree,forestvareb,ebvarsnew1,"ph1.eerrmodebraw");
  apply.ApplyAsFriend(tree,forestvareb,ebvarsnew2,"ph2.eerrmodebraw");
  apply.ApplyAsFriend(tree,forestvaree,eevarsnew1,"ph1.eerrmodeeraw");
  apply.ApplyAsFriend(tree,forestvaree,eevarsnew2,"ph2.eerrmodeeraw");  
  
  tree->SetAlias("ph1.eerrmod","ismc*(ph1.isbarrel*ph1.scrawe*ph1.eerrmodebraw + (!ph1.isbarrel)*(ph1.scrawe+ph1.scpse)*ph1.eerrmodeeraw) + (!ismc)*ph1.eerr");
  tree->SetAlias("ph2.eerrmod","ismc*(ph2.isbarrel*ph2.scrawe*ph2.eerrmodebraw + (!ph2.isbarrel)*(ph2.scrawe+ph2.scpse)*ph2.eerrmodeeraw) + (!ismc)*ph2.eerr");  
  
  hmc->SetAlias("ph1.eerrmod","ismc*(ph1.isbarrel*ph1.scrawe*ph1.eerrmodebraw + (!ph1.isbarrel)*(ph1.scrawe+ph1.scpse)*ph1.eerrmodeeraw) + (!ismc)*ph1.eerr");
  hmc->SetAlias("ph2.eerrmod","ismc*(ph2.isbarrel*ph2.scrawe*ph2.eerrmodebraw + (!ph2.isbarrel)*(ph2.scrawe+ph2.scpse)*ph2.eerrmodeeraw) + (!ismc)*ph2.eerr");
  
  hmc->SetAlias("massmvacorerrmod","0.5*mass*sqrt( ph1.eerrmod*ph1.eerrmod/ph1.e/ph1.e + ph1.esmearing*ph1.esmearing/ph1.e/ph1.e + ph2.eerrmod*ph2.eerrmod/ph2.e/ph2.e + ph2.esmearing*ph2.esmearing/ph2.e/ph2.e )");
  
  hdata->SetAlias("ph1.eerrmod","ph1.eerr");
  hdata->SetAlias("ph2.eerrmod","ph2.eerr");
  hdata->SetAlias("massmvacorerrmod","masserrsmeared");
  
//  return;
  
  TCut kcut = "(  (procidx>3 && procidx<12 && procidx>3 && ph1.ispromptgen && ph2.ispromptgen)*1.15 + (procidx>3 && procidx<12 && !ph1.ispromptgen && !ph2.ispromptgen)*1.0 + ( procidx>3 && procidx<12 && (ph1.ispromptgen || ph2.ispromptgen) && !(ph1.ispromptgen && ph2.ispromptgen) )*1.3  + (procidx==12)*1.15 + (procidx==13)*1.0 )";
  
  TCut diphojCut = "!( (procidx==7) && ((evt>=270000 && evt<=390000) || (evt>=480000 && evt<=510000) || (evt>=1170000 && evt<=1200000) || (evt>=1230000 && evt<=1260000)  || (evt>=1890000 && evt<=1920000)  || (evt>=2670000 && evt<=2730000) || (evt>=3090000 && evt<=3210000) || (evt>=3270000 && evt<=3300000) || (evt>=3510000 && evt<=3600000) || (evt>=4890000 && evt<=5010000) || (evt>=5340000 && evt<=5370000) || (evt>=5430000 && evt<=5500000)))";  
  
  TCut cut("xsecweight(procidx)*puweight(numPU,procidx)*(ph1.pt>(mass/3.0) && ph2.pt>(mass/4.0) && ph1.idmva>0.05 && ph2.idmva>-0.3)");
  //TCut cut("xsecweight(15)*puweight(numPU,15)*(ph1.pt>(mass/3.0) && ph2.pt>(mass/4.0) && ph1.idmva>-0.3 && ph2.idmva>-0.3)");
  
  
  //TCut cut("xsecweight(0)*puweight(numPU,0)*(ph1.pt>(mass/3.0) && ph2.pt>(mass/4.0) && ph1.idmva>-0.3 && ph2.idmva>-0.3)");
  //TCut cut("1.0*puweight(numPU,0)*(ph1.pt>(mass/3.0) && ph2.pt>(mass/4.0) && ph1.idmva>-0.3 && ph2.idmva>-0.3)");
  //TCut cut("xsecweight(procidx)*puweight(numPU,procidx)*(ph1.pt>(mass/3.0) && ph2.pt>(mass/4.0) && ph1.idmva>0.05 && ph2.idmva>-0.3)");
  //TCut cut("xsecweight(procidx)*puweight(numPU,procidx)*(ph1.idmva>0.05 && ph2.idmva>-0.3)");
  //TCut cut("xsecweight(procidx)*(ph1.pt>(mass/3.0) && ph2.pt>(mass/4.0) && ph1.idmva>-0.3 && ph2.idmva>-0.3)");
  //TCut cut("xsecweight(procidx)*(ph1.pt>(mass/3.0) && ph2.pt>(mass/4.0) && ph2.idmva>-0.3");
  //TCut cut("2.0*xsecweight(procidx)*(ph1.pt>(mass/3.0) && ph2.pt>(mass/4.0) && ph1.idmva>-0.3 && ph2.idmva>-0.3 && mass>100. && mass<180. && evt%2!=0)");
  //TCut cut("2.0*xsecweight(procidx)*(ph1.pt>(mass/3.0) && ph2.pt>(mass/4.0) && ph1.idmva>-0.3 && ph2.idmva>-0.3 && mass>100. && mass<180. && evt%2!=0)");

  //TCut masscut= kcut*diphojCut*cut*TCut("(mass>150)");
  
  TCut masscut= kcut*diphojCut*cut*TCut("(mass>160)");
  TCut masscutlow= kcut*diphojCut*cut*TCut("(mass>100 && mass<160)");

//   TCut masscut= cut*TCut("(mass>65 && mass<115 && ph1.pt>80)");
//   TCut masscutlow= cut*TCut("(mass>65 && mass<115)");
  
  //TCut masscutlow= kcut*diphojCut*cut*TCut("(mass>70 && mass<115)");
  TCut ebcut = masscut*TCut("ph1.isbarrel");
  TCut eecut = masscut*TCut("!ph1.isbarrel");
  
  TCut ebcutlow = masscutlow*TCut("ph1.isbarrel");
  TCut eecutlow = masscutlow*TCut("!ph1.isbarrel");  
  
  TCut ebebcutlow = masscutlow*TCut("(ph1.isbarrel&&ph2.isbarrel)");
  TCut nebcutlow = masscutlow*TCut("!(ph1.isbarrel&&ph2.isbarrel)");
  
//   plotstack("pt1barrel",tree,"ph1.pt",ebcut,"lead p_{T} (GeV)","# of events/GeV",0.0,250.0,50,true);
//   plotstack("pt1endcap",tree,"ph1.pt",eecut,"lead p_{T} (GeV)","# of events/GeV",0.0,250.0,50,true);

//   plotstack("sige1barrelalt",tree,"((1.0+ismc*0.06)*ph1.eerr - ismc*0.005)/ph1.e",ebcut,"#sigma_{E}/E","# of events",0.0,0.02,100,true,true);
//   plotstack("sige1barrellowalt",tree,"((1.0+ismc*0.06)*ph1.eerr - ismc*0.005)/ph1.e",ebcutlow,"#sigma_{E}/E","# of events",0.0,0.02,100,true,true);

  
//   plotstack("e55",tree,"ph1.eerr/ph1.e",ebcutlow,"e2x5left/eseed","# of events",0.0,0.02,100,true,true);  
//   plotstack("e55",tree,"ph1.eerr/ph1.e",eecutlow,"e2x5left/eseed","# of events",0.0,0.04,100,true,true);  
//   
//   plotstack("e55",tree,"ph1.eerr/ph1.e",ebcut,"e2x5left/eseed","# of events",0.0,0.02,100,true,true);  
//   plotstack("e55",tree,"ph1.eerr/ph1.e",eecut,"e2x5left/eseed","# of events",0.0,0.04,100,true,true);  
  
//   plotstack("pt1barrel",hdata,hmc,"massmvacorerrmod",ebebcutlow,"lead p_{T} (GeV)","# of events/GeV",0.0,5.0,100,true);
//   plotstack("pt1barrel",hdata,hmc,"massmvacorerrmod",nebcutlow,"lead p_{T} (GeV)","# of events/GeV",0.0,5.0,100,true);
// 
//   //return;
//   cbfits("massreseb", ebebcutlow, hmc, hdata, "massmvacorerrmod", "massmvacorerrmod", 20, 0., 3.5, "x", "#sigma_{m} (GeV)",1.0);
//   cbfits("massresneb", nebcutlow, hmc, hdata, "massmvacorerrmod", "massmvacorerrmod", 20, 1.2, 5.0, "x", "#sigma_{m} (GeV)",3.0);
//   
//   cbfits("massresebnom", ebebcutlow, hmc, hdata, "masserrsmeared", "masserrsmeared", 20, 0., 3.5, "x", "#sigma_{m} (GeV)",1.0);
//   cbfits("massresnebnom", nebcutlow, hmc, hdata, "masserrsmeared", "masserrsmeared", 20, 1.2, 5.0, "x", "#sigma_{m} (GeV)",3.0);
// 
//   
//   return;
// 
//   
//   plotstack("e55",tree,"ph1.eerrmod/ph1.e",ebcutlow,"e2x5left/eseed","# of events",0.0,0.02,100,true,true);  
//   plotstack("e55",tree,"ph1.eerrmod/ph1.e",eecutlow,"e2x5left/eseed","# of events",0.0,0.04,100,true,true);  
//   
//   plotstack("e55",tree,"ph1.eerrmod/ph1.e",ebcut,"e2x5left/eseed","# of events",0.0,0.02,100,true,true);  
//   plotstack("e55",tree,"ph1.eerrmod/ph1.e",eecut,"e2x5left/eseed","# of events",0.0,0.04,100,true,true);  
//   
//   return;
  

/*  plotstack("e55",hdata,hmc,"ph1.eerr/ph1.e",ebcutlow,"#sigma_{E}/E","# of events",0.0,0.02,100,true);  
  plotstack("e55",hdata,hmc,"ph1.eerr/ph1.e",eecutlow,"#sigma_{E}/E","# of events",0.0,0.04,100,true);    
  
  plotstack("e55",hdata,hmc,"ph1.eerr/ph1.e",ebcut,"#sigma_{E}/E","# of events",0.0,0.02,100,true);  
  plotstack("e55",hdata,hmc,"ph1.eerr/ph1.e",eecut,"#sigma_{E}/E","# of events",0.0,0.04,100,true);*/   
  
//   plotstack("sigmaEbarrelEle",hdata,hmc,"ph1.eerrmod/ph1.e",ebcutlow,"#sigma_{E}/E","# of events",0.0,0.02,100,true);  
//   plotstack("sigmaEendcapEle",hdata,hmc,"ph1.eerrmod/ph1.e",eecutlow,"#sigma_{E}/E","# of events",0.0,0.04,100,true);  
//   
//   plotstack("sigmaEbarrelEleHighPt",hdata,hmc,"ph1.eerrmod/ph1.e",ebcut,"#sigma_{E}/E","# of events",0.0,0.02,100,true);  
//   plotstack("sigmaEendcapEleHighPt",hdata,hmc,"ph1.eerrmod/ph1.e",eecut,"#sigma_{E}/E","# of events",0.0,0.04,100,true);  
  
  plotstack("sigmaEbarrelPho100_160",tree,"ph1.eerrmod/ph1.e",ebcutlow,"#sigma_{E}/E","# of events",0.0,0.02,100,true,true);  
  plotstack("sigmaEendcapPho100_160",tree,"ph1.eerrmod/ph1.e",eecutlow,"#sigma_{E}/E","# of events",0.0,0.04,100,true,true);  
  
  plotstack("sigmaEbarrelPho160_inf",tree,"ph1.eerrmod/ph1.e",ebcut,"#sigma_{E}/E","# of events",0.0,0.02,100,true,true);  
  plotstack("sigmaEendcapPho160_inf",tree,"ph1.eerrmod/ph1.e",eecut,"#sigma_{E}/E","# of events",0.0,0.04,100,true,true);   
  
  return;
  
     plotstack("e55",hdata,hmc,"(1.0+ismc*0.0016)*ph1.eseed/ph1.scrawe",ebcutlow,"e55seed/eseed","# of events",0.92,1.1,100,true);
     plotstack("e55",hdata,hmc,"(1.0+ismc*0.00)*ph1.eseed/ph1.scrawe",ebcutlow,"e55seed/eseed","# of events",0.92,1.1,100,true);  

  
    plotstack("e55",hdata,hmc,"ph1.e3x3seed/ph1.scrawe",ebcutlow,"e55seed/eseed","# of events",0.92,1.1,100,true);
    plotstack("e55",hdata,hmc,"ph1.e3x3seed/ph1.scrawe",eecutlow,"e55seed/eseed","# of events",0.92,1.1,100,true);

  
   plotstack("e55",hdata,hmc,"ph1.e5x5seed/ph1.scrawe",ebcutlow,"e55seed/eseed","# of events",0.92,1.1,100,true);
  
   plotstack("e55",hdata,hmc,"ph1.e5x5seed/ph1.scrawe",eecutlow,"e55seed/eseed","# of events",0.92,1.1,100,true);
   plotstack("e55",hdata,hmc,"ph1.eseed/ph1.scrawe",eecutlow,"e55seed/eseed","# of events",0.92,1.1,100,true);

   return;
  
  
   plotstack("e55",hdata,hmc,"(1.0+ismc*0.00)*ph1.e3x3seed/ph1.eseed + 0.00*ismc",ebcutlow,"e55seed/eseed","# of events",0.8,1.0,100,true);
  
   plotstack("e55",hdata,hmc,"(1.0+ismc*0.0045)*ph1.e3x3seed/ph1.eseed + 0.00*ismc",ebcutlow,"e55seed/eseed","# of events",0.8,1.0,100,true);
   
   return;
   //plotstack("e55",hdata,hmc,"(1.0+ismc*0.00)*ph1.e3x3seed/ph1.eseed - 0.00*ismc",eecutlow,"e55seed/eseed","# of events",0.8,1.0,100,true);  

   //plotstack("e55",hdata,hmc,"(1.0+ismc*0.0086)*ph1.e3x3seed/ph1.eseed - 0.000*ismc",eecutlow,"e55seed/eseed","# of events",0.8,1.0,100,true);
   plotstack("e55",hdata,hmc,"(1.0 + ismc*0.007)*ph1.e3x3seed/ph1.eseed + 0.00*ismc",eecutlow,"e55seed/eseed","# of events",0.8,1.0,100,true);  
  
   plotstack("e55",hdata,hmc,"ph1.e5x5seed/ph1.eseed",ebcutlow,"e55seed/eseed","# of events",0.92,1.1,100,true);

   
   plotstack("e55",hdata,hmc,"TMath::Min(1.0,(1.0+ismc*0.0016)*ph1.e5x5seed/ph1.eseed)",ebcutlow,"e55seed/eseed","# of events",0.92,1.1,100,true);
   plotstack("e55",hdata,hmc,"(1.0+ismc*0.00)*ph1.e5x5seed/ph1.eseed",eecutlow,"e55seed/eseed","# of events",0.92,1.1,100,true);

   plotstack("e55",hdata,hmc,"TMath::Min(1.0,(1.0+ismc*0.00)*ph1.eseed/ph1.scrawe)",ebcutlow,"e55seed/eseed","# of events",0.92,1.1,100,true);
   plotstack("e55",hdata,hmc,"TMath::Min(1.0,(1.0+ismc*0.0022)*ph1.eseed/ph1.scrawe)",eecutlow,"e55seed/eseed","# of events",0.92,1.1,100,true);   
   
return;



  
  
  //good scaling numbers
//   plotstack("e55",hdata,hmc,"(1.0+ismc*0.0016)*ph1.e5x5/ph1.scrawe",ebcutlow,"e55/escrawe","# of events",0.92,1.0,100,true);
//   plotstack("e55",hdata,hmc,"(1.0+ismc*0.0022)*ph1.e5x5/ph1.scrawe",eecutlow,"e55/escrawe","# of events",0.92,1.0,100,true);
// 
//   plotstack("e55",hdata,hmc,"(1.0+ismc*0.0021)*ph1.e5x5/ph1.scrawe - ismc*0.001",ebcutlow,"e55/escrawe","# of events",0.92,1.0,100,true);
//   plotstack("e55",hdata,hmc,"(1.0+ismc*0.004)*ph1.e5x5/ph1.scrawe - ismc*0.0015",eecutlow,"e55/escrawe","# of events",0.92,1.0,100,true);  
//   
//   plotstack("e55",hdata,hmc,"(1.0+ismc*0.012)*ph1.emaxseed/ph1.eseed",ebcutlow,"emax/eseed","# of events",0.0,1.0,100,true);  
//   plotstack("e55",hdata,hmc,"(1.0+ismc*0.005)*ph1.emaxseed/ph1.eseed",eecutlow,"emax/eseed","# of events",0.0,1.0,100,true);  
//  
//   plotstack("e55",hdata,hmc,"(1.0-ismc*0.06)*ph1.eleftseed/ph1.eseed",ebcutlow,"emax/eseed","# of events",0.0,0.15,100,true);  
//   plotstack("e55",hdata,hmc,"(1.0-ismc*0.04)*ph1.eleftseed/ph1.eseed",eecutlow,"emax/eseed","# of events",0.0,0.15,100,true);    
//   
//   plotstack("e55",hdata,hmc,"(1.0+ismc*0.006)*ph1.e2x5maxseed/ph1.eseed",ebcutlow,"e2x5max/eseed","# of events",0.8,1.0,100,true);  
//   plotstack("e55",hdata,hmc,"(1.0+ismc*0.0075)*ph1.e2x5maxseed/ph1.eseed",eecutlow,"e2x5max/eseed","# of events",0.8,1.0,100,true);   
//  
//   plotstack("e55",hdata,hmc,"(1.0-ismc*0.09)*ph1.e2x5leftseed/ph1.eseed",ebcutlow,"e2x5left/eseed","# of events",0.0,0.3,100,true);  
//  plotstack("e55",hdata,hmc,"(1.0-ismc*0.13)*ph1.e2x5leftseed/ph1.eseed",eecutlow,"e2x5left/eseed","# of events",0.0,0.3,100,true);     

/*  plotstack("e55",hdata,hmc,"(1.0-ismc*0.007)*ph1.sigiphiphiseed",ebcutlow,"e2x5left/eseed","# of events",0.0,0.025,100,true);  
  plotstack("e55",hdata,hmc,"(1.0-ismc*0.00)*ph1.sigiphiphiseed",eecutlow,"e2x5left/eseed","# of events",0.0,0.08,100,true);   */  
  
//   plotstack("e55",hdata,hmc,"(1.0-ismc*0.00)*ph1.e2ndseed/ph1.eseed",ebcutlow,"e2x5left/eseed","# of events",0.0,0.45,100,true);  
//   plotstack("e55",hdata,hmc,"(1.0+ismc*0.02)*ph1.e2ndseed/ph1.eseed",eecutlow,"e2x5left/eseed","# of events",0.0,0.45,100,true);  
  
  return;
  

  
  return;
  
  plotstack("massbarrel",tree,"mass",ebcut,"#sigma_{E}/E","# of events",70.0,125,100,false,true);
  plotstack("massbarrellow",tree,"mass",ebcutlow,"#sigma_{E}/E","# of events",70.0,125,100,false,true);
  plotstack("massendcap",tree,"mass",eecut,"#sigma_{E}/E","# of events",70.0,125,100,false,true);
  plotstack("massendcaplow",tree,"mass",eecutlow,"#sigma_{E}/E","# of events",70.0,125,100,false,true);
  
    return;
  
  //plotstack("sige1barrelalt",tree,"((1.0+ismc*0.05)*ph1.eerr - ismc*0.005)/ph1.e",ebcut,"#sigma_{E}/E","# of events",0.0,0.02,100,true,true);
  
  plotstack("sige1barrel",tree,"((1.0+ismc*0.02693)*ph1.eerr - ismc*0.0042793)/ph1.e",ebcut,"#sigma_{E}/E","# of events",0.0,0.02,100,true,true,false,false,"1.1*((1.0+ismc*0.02693)*ph1.eerr - ismc*0.0042793)/ph1.e","0.9*((1.0+ismc*0.02693)*ph1.eerr - ismc*0.0042793)/ph1.e");
  
  plotstack("sige1barrelalt",tree,"((1.0+ismc*0.05)*ph1.eerr - ismc*0.005)/ph1.e",ebcut,"#sigma_{E}/E","# of events",0.0,0.02,100,true,true,false,false,"1.1*((1.0+ismc*0.15)*ph1.eerr - ismc*0.005)/ph1.e","0.9*((1.0+ismc*0.05)*ph1.eerr - ismc*0.005)/ph1.e");  
  
  plotstack("sige1endcap",tree,"((1.0+ismc*0.01372)*ph1.eerr + ismc*0.000156943)/ph1.e",eecut,"#sigma_{E}/E","# of events",0.0,0.04,100,true,true,false,false,"1.1*((1.0+ismc*0.01372)*ph1.eerr + ismc*0.000156943)/ph1.e","0.9*((1.0+ismc*0.01372)*ph1.eerr + ismc*0.000156943)/ph1.e");  

  
  plotstack("sige1barrellow",tree,"((1.0+ismc*0.02693)*ph1.eerr - ismc*0.0042793)/ph1.e",ebcutlow,"#sigma_{E}/E","# of events",0.0,0.02,100,true,true,false,false,"1.1*((1.0+ismc*0.02693)*ph1.eerr - ismc*0.0042793)/ph1.e","0.9*((1.0+ismc*0.02693)*ph1.eerr - ismc*0.0042793)/ph1.e");
  
  plotstack("sige1barrellowalt",tree,"((1.0+ismc*0.05)*ph1.eerr - ismc*0.005)/ph1.e",ebcutlow,"#sigma_{E}/E","# of events",0.0,0.02,100,true,true,false,false,"1.1*((1.0+ismc*0.15)*ph1.eerr - ismc*0.005)/ph1.e","0.9*((1.0+ismc*0.05)*ph1.eerr - ismc*0.005)/ph1.e");  
  
  plotstack("sige1endcaplow",tree,"((1.0+ismc*0.01372)*ph1.eerr + ismc*0.000156943)/ph1.e",eecutlow,"#sigma_{E}/E","# of events",0.0,0.04,100,true,true,false,false,"1.1*((1.0+ismc*0.01372)*ph1.eerr + ismc*0.000156943)/ph1.e","0.9*((1.0+ismc*0.01372)*ph1.eerr + ismc*0.000156943)/ph1.e");  
  
 // plotstack("sige1barrellowalt",tree,"((1.0+ismc*0.05)*ph1.eerr - ismc*0.005)/ph1.e",ebcutlow,"#sigma_{E}/E","# of events",0.0,0.02,100,true,true,false,false,"((1.0+ismc*0.15)*ph1.eerr - ismc*0.005)/ph1.e","((1.0-ismc*0.05)*ph1.eerr - ismc*0.005)/ph1.e");
 

  
 // plotstack("sige1barrelalt2",tree,"((1.0+ismc*0.10)*ph1.eerr - ismc*0.006)/ph1.e",ebcut,"#sigma_{E}/E","# of events",0.0,0.02,100,true,true,false,false,"((1.0+ismc*0.2)*ph1.eerr - ismc*0.005)/ph1.e","((1.0-ismc*0)*ph1.eerr - ismc*0.005)/ph1.e");  
  
 // plotstack("sige1barrellowalt2",tree,"((1.0+ismc*0.10)*ph1.eerr - ismc*0.006)/ph1.e",ebcutlow,"#sigma_{E}/E","# of events",0.0,0.02,100,true,true,false,false,"((1.0+ismc*0.2)*ph1.eerr - ismc*0.005)/ph1.e","((1.0-ismc*0)*ph1.eerr - ismc*0.005)/ph1.e");  
  
//   plotstack("sige1barrellowalt",tree,"((1.0+ismc*0.05)*ph1.eerr - ismc*0.005)/ph1.e",ebcutlow,"#sigma_{E}/E","# of events",0.0,0.02,100,true,true,false,false,"((1.0+ismc*0.15)*ph1.eerr - ismc*0.005)/ph1.e","((1.0-ismc*0.05)*ph1.eerr - ismc*0.005)/ph1.e");
  
  

  return;
  
  plotstack("sige1barrel",tree,"((1.0+ismc*0.02693)*ph1.eerr - ismc*0.0042793)/ph1.e",ebcut,"#sigma_{E}/E","# of events",0.0,0.02,100,true,true);
  plotstack("sige1endcap",tree,"((1.0+ismc*0.01372)*ph1.eerr + ismc*0.000156943)/ph1.e",eecut,"#sigma_{E}/E","# of events",0.0,0.04,100,true,true);

  plotstack("sige1barrellow",tree,"((1.0+ismc*0.02693)*ph1.eerr - ismc*0.0042793)/ph1.e",ebcutlow,"#sigma_{E}/E","# of events",0.0,0.02,100,true,true);
  plotstack("sige1endcaplow",tree,"((1.0+ismc*0.01372)*ph1.eerr + ismc*0.000156943)/ph1.e",eecutlow,"#sigma_{E}/E","# of events",0.0,0.04,100,true,true);  
  
//   plotstack("idmvabarrel",tree,"ph1.idmva",ebcut,"#sigma_{E}/E","# of events",-0.5,0.5,100,true,true);
//   plotstack("idmvaendcap",tree,"ph1.idmva",eecut,"#sigma_{E}/E","# of events",-0.5,0.5,100,true,true);
 
  
  return;
}
