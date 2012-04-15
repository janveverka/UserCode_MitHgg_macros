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

void fitbkgdataglobep3(bool dobands=false, bool dosig=false) { 
  
  
  gROOT->Macro("MitStyle.C");
  gStyle->SetErrorX(0); 
  gStyle->SetOptStat(0);
  gROOT->ForceStyle();  


  gSystem->cd("./bambuOutputFeb16smmvavbf/databkgoverlay");
  //gSystem->cd("./globeOutputDec8/databkgoverlay");
  //gSystem->cd("./bambuOutputNov20smetar9/databkgcheb");
//  gSystem->cd("./test");

    //TFile *fdata = new TFile("/scratch/bendavid/root/globeOutputNov15/CMS-HGG_4686pb_r9scale.root","READ");
    
   
    //TFile *fdata = new TFile("/scratch/bendavid/root/globeOutputDec8/CMS-HGG_4763pb.root","READ");
    
    
    TFile *fdata = new TFile("/scratch/bendavid/root/bambuOutputFeb16smmvavbf/CMS-HGG-data.root","READ");
    
    RooWorkspace* win = (RooWorkspace*)fdata->Get("cms_hgg_workspace_data");
    RooRealVar *IntLumi = win->var("IntLumi");

    
    RooWorkspace *wcomb=0;
    
    if (dosig) {
      TFile *fmc = new TFile("/scratch/bendavid/root/bambuOutputFeb16smmvavbf/signalplots/combsignal.root","READ");
      wcomb = (RooWorkspace*)fmc->Get("wcomb");
    }
    
  //  win->Print();
  //  return;
    
    const double massmax = 180.0;
  
    RooRealVar *hmass = win->var("CMS_hgg_mass");
    //RooRealVar *hmass = win->var("mass");
    hmass->setRange(100,massmax);
    hmass->setBins(4.0*(massmax-100.0));
    hmass->SetTitle("m_{#gamma#gamma}");
    hmass->setUnit("GeV/c^{2}");
    hmass->setRange("fitrange",100,massmax);
    
    std::vector<TString> catdesc;
    //Feb2/16 smmvavbf
    catdesc.push_back("#scale[0.5]{BDT >= 0.89}");
    catdesc.push_back("#scale[0.5]{0.72 <= BDT < 0.89}");
    catdesc.push_back("#scale[0.5]{0.55 <= BDT < 0.72}");
    catdesc.push_back("#scale[0.5]{0.05 <= BDT < 0.55}");
    catdesc.push_back("#scale[0.5]{BDT >= 0.05 VBF Tag}");
    catdesc.push_back("#scale[0.5]{All Categories Combined}");  
    
    std::vector<TString> catnames;  
    catnames.push_back("cat0");
    catnames.push_back("cat1");
    catnames.push_back("cat2");
    catnames.push_back("cat3");
    catnames.push_back("cat4");

     
    std::vector<RooAbsData*> datav;
    std::vector<RooAbsPdf*> pdfv;
    std::vector<RooAbsPdf*> pdfuv;

    RooArgList updfs;
    RooArgList aconsts;
    
    
    std::vector<RooRealVar*> coeffv;
    std::vector<RooRealVar*> normv;
    std::vector<RooAbsReal*> normuv;

    
    RooWorkspace *w = new RooWorkspace("wbkg","wbkg") ;
  
    RooRealVar *p1first=0;
    RooRealVar *p2first=0;
    RooRealVar *nfirst=0;
    
    RooCategory finalcat("finalcat","finalcat") ;  
    RooSimultaneous fullbkgpdf("fullbkgpdf","fullbkgpdf",finalcat);
    RooDataSet datacomb("datacomb","",RooArgList(*hmass,finalcat)) ;
    RooDataSet *datacombcat = new RooDataSet("data_combcat","",RooArgList(*hmass)) ;
    RooAbsPdf *combcatpdf = 0;
    catnames.push_back("combcat");    
    RooFitResult *combresult = 0;
    for (UInt_t icat=0; icat<catnames.size(); ++icat) {
      TString catname = catnames.at(icat);
      finalcat.defineType(catname);
      
      RooDataSet *indata = (RooDataSet*)win->data(TString("data_mass_")+catname);
      RooDataSet *data = 0;
      if (indata) {
        data = new RooDataSet(TString("data_")+catname,"",*hmass,Import(*indata));
        RooDataSet *datacat = new RooDataSet(TString("datacat")+catname,"",*hmass,Index(finalcat),Import(catname,*data)) ;
        datacomb.append(*datacat);
        datacombcat->append(*data);
      }
      else {
        data = datacombcat;
      }
      
      RooRealVar *p1 = new RooRealVar(TString::Format("CMS_hgg_%s_p1",catname.Data()),"",0.1,-10.0,10.0);
      RooRealVar *p2 = new RooRealVar(TString::Format("CMS_hgg_%s_p2",catname.Data()),"",0.1,-10.0,10.0);
      RooRealVar *p3 = new RooRealVar(TString::Format("CMS_hgg_%s_p3",catname.Data()),"",0.1,-10.0,10.0);
      RooRealVar *p4 = new RooRealVar(TString::Format("CMS_hgg_%s_p4",catname.Data()),"",0.1,-10.0,10.0);      
      RooRealVar *p5 = new RooRealVar(TString::Format("CMS_hgg_%s_p5",catname.Data()),"",0.1,-10.0,10.0);      
      

      
      coeffv.push_back(p1);
      coeffv.push_back(p2);
      coeffv.push_back(p3);      
      coeffv.push_back(p4);      
      coeffv.push_back(p5);      
      

      
      RooFormulaVar *p1mod = new RooFormulaVar(TString("p1mod"+catname),"","@0*@0",*p1);
      RooFormulaVar *p2mod = new RooFormulaVar(TString("p2mod"+catname),"","@0*@0",*p2);
      RooFormulaVar *p3mod = new RooFormulaVar(TString("p3mod"+catname),"","@0*@0",*p3);
      RooFormulaVar *p4mod = new RooFormulaVar(TString("p4mod"+catname),"","@0*@0",*p4);
      RooFormulaVar *p5mod = new RooFormulaVar(TString("p5mod"+catname),"","@0*@0",*p5);

          
      RooRealVar *nbkg = new RooRealVar(TString::Format("CMS_hgg_%s_bkgshape_norm",catname.Data()),"",800.0,0.0,25e3);
      normv.push_back(nbkg);

      RooRealVar *cbkg = new RooRealVar(TString::Format("cbkg%s",catname.Data()),"",0.0,0.0,1e3);
      //cbkg->setVal(data->sumEntries());
      
      
      
      
      RooAbsPdf *bkgshape = 0;
      
      if (icat==4) {
        bkgshape = new RooBernstein(TString::Format("CMS_hgg_%s_bkgshape",catname.Data()),"",*hmass,RooArgList(RooConst(1.0),*p1mod,*p2mod,*p3mod));
      }      
      else if (icat==0) {
        bkgshape = new RooBernstein(TString::Format("CMS_hgg_%s_bkgshape",catname.Data()),"",*hmass,RooArgList(RooConst(1.0),*p1mod,*p2mod,*p3mod,*p4mod));
      }      
      else {
        bkgshape = new RooBernstein(TString::Format("CMS_hgg_%s_bkgshape",catname.Data()),"",*hmass,RooArgList(RooConst(1.0),*p1mod,*p2mod,*p3mod,*p4mod,*p5mod));
      }
      
      
      
      RooAbsPdf *bkgpdf = 0;
    
      if (data==datacombcat) {
        nbkg->setVal(8000.0);
        nbkg->SetName(TString::Format("CMS_hgg_%s_bkgshapesingle_norm",catname.Data()));
        bkgshape->SetName(TString::Format("CMS_hgg_%s_bkgshapesingle",catname.Data()));
        bkgpdf = new RooExtendPdf(TString("bkgpdfsingle")+catname,"",*bkgshape,*nbkg);
        bkgpdf->fitTo(*data,Strategy(1),Minos(kFALSE),Save(kTRUE));
        w->import(*bkgpdf);
        
        
        RooArgList subpdfs;
        for (int ipdf=0; ipdf<pdfv.size(); ++ipdf) {
          subpdfs.add(*pdfv.at(ipdf));
        }
        bkgpdf = new RooAddPdf(TString("bkgpdf")+catname,"",subpdfs);
        bkgshape = bkgpdf;
        combcatpdf = bkgpdf;
        
      }
      else {
        bkgpdf = new RooExtendPdf(TString("bkgpdf")+catname,"",*bkgshape,*nbkg);
        fullbkgpdf.addPdf(*bkgpdf,catname);
        updfs.add(*bkgshape);
      }
      
      pdfuv.push_back(bkgshape);
      
      if (icat<(catnames.size()-2)) {
        normuv.push_back(cbkg);
        aconsts.add(*cbkg);
      }
      
      
      
      datav.push_back(data);
      pdfv.push_back(bkgpdf);
      
      RooDataHist *databinned = new RooDataHist(TString("databinned_")+catname,"",*hmass,*data);
      
      w->import(*data);
      w->import(*databinned);
      
      if (icat==0) {
        p1first=p1;
        p2first=p2;
        nfirst=nbkg;
      }
      
      
    }
    
  gStyle->SetOptStat(1110);
  
  fullbkgpdf.fitTo(datacomb,Strategy(1),Minos(kFALSE),Save(kTRUE));
  RooFitResult *fullbkgfitres = fullbkgpdf.fitTo(datacomb,Strategy(2),Minos(kFALSE),Save(kTRUE));

  //return;
  
  const UInt_t nbins = massmax-100.0;
  
  std::vector<RooCurve*> fitcurves;
  
    
  for (UInt_t icat=0; icat<catnames.size(); ++icat) {
    TString catname = catnames.at(icat);
    
    RooFitResult *plotresult=0;
    //if (catname=="combcat") plotresult = combresult;
    //else plotresult = fullbkgfitres;
    plotresult = fullbkgfitres;
    

    
    
    
    
    //testpdf.fitTo(*datav.at(icat));
    
    TCanvas *cbkg = new TCanvas;
    RooPlot *plot = hmass->frame(Bins(nbins),Range("fitrange"));
    datav.at(icat)->plotOn(plot,LineColor(kWhite),MarkerColor(kWhite),Invisible());
//    pdfv.at(icat)->plotOn(plot,FillColor(kYellow),Range("fitrange"),NormRange("fitrange"),VisualizeError(*plotresult,2.0,kTRUE),Normalization(1.0));    
//    pdfv.at(icat)->plotOn(plot,FillColor(kGreen),Range("fitrange"),NormRange("fitrange"),VisualizeError(*plotresult,1.0,kTRUE),Normalization(1.0));
//     pdfv.at(icat)->plotOn(plot,FillColor(kYellow),Range("fitrange"),NormRange("fitrange"),VisualizeError(*fullbkgfitres,RooArgSet(*p1first,*p2first),2.0,kFALSE));    
//     pdfv.at(icat)->plotOn(plot,FillColor(kGreen),Range("fitrange"),NormRange("fitrange"),VisualizeError(*fullbkgfitres,RooArgSet(*p1first,*p2first),1.0,kFALSE));
    pdfv.at(icat)->plotOn(plot,RooFit::LineColor(kRed),Range("fitrange"),NormRange("fitrange"));
    
    datav.at(icat)->plotOn(plot);    

    
    plot->SetTitle("");      
    plot->SetMinimum(0.0);
    plot->SetMaximum(1.40*plot->GetMaximum());
    plot->GetXaxis()->SetTitle("m_{#gamma#gamma} (GeV/c^{2})");
    plot->Draw();       
    
    
    
    //bool dobands = true;
    
    TGraphAsymmErrors *onesigma = new TGraphAsymmErrors();
    TGraphAsymmErrors *twosigma = new TGraphAsymmErrors();
    
    if (dobands) {
      
      
      
      
      
      RooAbsPdf *cpdf = pdfuv.at(icat); 
      RooRealVar *nlim = new RooRealVar(TString::Format("nlim%s",catnames.at(icat).Data()),"",0.0,0.0,10.0);
      nlim->removeRange();
      
      RooAddition sumcatsnm1("sumcatsnm1","",aconsts);
      RooFormulaVar *nlast = new RooFormulaVar("nlast","","TMath::Max(0.1,@0-@1)",RooArgList(*nlim,sumcatsnm1));
      
     
      
      RooCurve *nomcurve = dynamic_cast<RooCurve*>(plot->getObject(1));
      fitcurves.push_back(nomcurve);
      
      bool iscombcat = true;
      RooAbsData *datanorm = datav.at(icat);
      if (catname == "combcat") {
        datanorm = &datacomb;
        normuv.push_back(nlast);
        iscombcat = true;
        //normv.at(0)->setConstant();
      }
      
      for (int i=1; i<(plot->GetXaxis()->GetNbins()+1); ++i) {
      //for (int i=1; i<3; ++i) {
        double lowedge = plot->GetXaxis()->GetBinLowEdge(i);
        double upedge = plot->GetXaxis()->GetBinUpEdge(i);
        double center = plot->GetXaxis()->GetBinCenter(i);
        
        double nombkg = nomcurve->interpolate(center);
        nlim->setVal(nombkg);
        hmass->setRange("errRange",lowedge,upedge);
        RooAbsPdf *epdf = 0;
        if (catname == "combcat") {
          epdf = new RooSimultaneous("epdf","",finalcat);
          for (int jcat=0; jcat<(catnames.size()-1); ++jcat) {
            RooRealVar *rvar = dynamic_cast<RooRealVar*>(normuv.at(jcat));
            if (rvar) rvar->setVal(fitcurves.at(jcat)->interpolate(center));
            RooExtendPdf *ecpdf = new RooExtendPdf(TString::Format("ecpdf%s",catnames.at(jcat).Data()),"",*pdfuv.at(jcat),*normuv.at(jcat),"errRange");
            static_cast<RooSimultaneous*>(epdf)->addPdf(*ecpdf,catnames.at(jcat));
          }
        }
        else {
          epdf = new RooExtendPdf("epdf","",*cpdf,*nlim,"errRange");
        }
        
        RooAbsReal *nll = epdf->createNLL(*datanorm,Extended(),NumCPU(16));
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
        
        //epdf->fitTo(*datav.at(icat),Strategy(2),Minos(*nlim));
        
        printf("errlo = %5f, errhi = %5f\n",nlim->getErrorLo(),nlim->getErrorHi());
        
        onesigma->SetPoint(i-1,center,nombkg);
        onesigma->SetPointError(i-1,0.,0.,-nlim->getErrorLo(),nlim->getErrorHi());
        
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
        
        
        delete nll;
        delete epdf;
      }
      hmass->setRange("errRange",100,massmax);
      //normv.at(0)->setConstant(kFALSE);
      //delete nlim;
      
      onesigma->Print("V");
      

    
    
      
        
      //testpdf.plotOn(plot,LineColor(kBlue),LineStyle(kDashed));
      //plot->SetTitle(catname);
  

      twosigma->SetLineColor(kGreen);
      twosigma->SetFillColor(kGreen);
      twosigma->SetMarkerColor(kGreen);
      twosigma->Draw("L3 SAME");     
      
      onesigma->SetLineColor(kYellow);
      onesigma->SetFillColor(kYellow);
      onesigma->SetMarkerColor(kYellow);
      onesigma->Draw("L3 SAME");
        
  //     pdfv.at(icat)->plotOn(plot,RooFit::LineColor(kRed),Range("fitrange"),NormRange("fitrange"));
  //     
  //     
      plot->Draw("SAME");
    }
    
  
    TLegend *legmc = new TLegend(0.68,0.70,0.97,0.90);  
    legmc->AddEntry(plot->getObject(2),"Data","LPE");
    legmc->AddEntry(plot->getObject(1),"Bkg Model","L");
    if (dobands) {
      legmc->AddEntry(onesigma,"#pm1 #sigma","F");  
      legmc->AddEntry(twosigma,"#pm2 #sigma","F");  
    }
//     legmc->AddEntry(plot->getObject(2),"#pm1 #sigma","F");  
//     legmc->AddEntry(plot->getObject(1),"#pm2 #sigma","F");  
    //legmc->AddEntry(plot->getObject(5),"Chebychev","L");  
    legmc->SetBorderSize(0);
    legmc->SetFillStyle(0);
    legmc->Draw();   
    
    
    
    TLatex lat(103.0,0.9*plot->GetMaximum(),"#scale[0.7]{#splitline{CMS preliminary}{#sqrt{s} = 7 TeV L = 4.76 fb^{-1}}}");
    //TLatex lat(103.0,0.9*plot->GetMaximum(),"#scale[0.7]{#splitline{CMS preliminary}{#sqrt{s} = 7 TeV L = 4.72 fb^{-1}}}");
    
    //TLatex lat(103.0,0.9*plot->GetMaximum(),"#scale[0.7]{#splitline{CMS preliminary}{#sqrt{s} = 7 TeV L = 0.51 fb^{-1}}}");
    lat.Draw();
    
    TLatex lat2(103.0,0.75*plot->GetMaximum(),catdesc.at(icat));
    lat2.Draw();
    
    cbkg->SaveAs(TString("databkg") + catname + TString(".pdf"));
    cbkg->SaveAs(TString("databkg") + catname + TString(".eps"));
    cbkg->SaveAs(TString("databkg") + catname + TString(".root"));
    
    
    if (dosig) {
      RooAbsPdf *sigpdf = wcomb->pdf(TString::Format("fullpdf_%s",catnames.at(icat).Data()));
      
      TCanvas *cbkgsig = new TCanvas;
      RooPlot *plotsig = hmass->frame(Bins(nbins),Range("fitrange"));
      datav.at(icat)->plotOn(plotsig,LineColor(kWhite),MarkerColor(kWhite),Invisible());
      
      double lumi = 4763.;
      //double lumi = 4717.;
      double norm = 1.0*lumi*sigpdf->expectedEvents(*hmass);
      sigpdf->plotOn(plotsig,Normalization(norm,RooAbsPdf::NumEvent),DrawOption("F"),LineColor(kBlue),FillStyle(1001),FillColor(19));
      sigpdf->plotOn(plotsig,Normalization(norm,RooAbsPdf::NumEvent),LineColor(kBlue));
      
      pdfv.at(icat)->plotOn(plotsig,RooFit::LineColor(kRed),Range("fitrange"),NormRange("fitrange"));
      
      datav.at(icat)->plotOn(plotsig);    

      
      plotsig->SetTitle("");      
      plotsig->SetMinimum(0.0);
      plotsig->SetMaximum(1.40*plotsig->GetMaximum());
      plotsig->GetXaxis()->SetTitle("m_{#gamma#gamma} (GeV/c^{2})");
      plotsig->Draw();       
      
      if (dobands) {
        twosigma->Draw("L3 SAME");     
        onesigma->Draw("L3 SAME");     
        
        plotsig->Draw("SAME");
      }
      
      TLegend *legsig = new TLegend(0.68,0.70,0.97,0.90);  
      legsig->AddEntry(plotsig->getObject(4),"Data","LPE");
      legsig->AddEntry(plotsig->getObject(3),"Bkg Model","L");
      if (dobands) {
        legsig->AddEntry(onesigma,"#pm1 #sigma","F");  
        legsig->AddEntry(twosigma,"#pm2 #sigma","F");  
      }
      legsig->AddEntry(plotsig->getObject(1),"1xSM m_{H}=120 GeV");            
      legsig->SetBorderSize(0);
      legsig->SetFillStyle(0);
      legsig->Draw();   

      lat.Draw();
      lat2.Draw();

      cbkgsig->SaveAs(TString("databkgoverlay") + catname + TString(".pdf"));
      cbkgsig->SaveAs(TString("databkgoverlay") + catname + TString(".eps"));
      cbkgsig->SaveAs(TString("databkgoverlay") + catname + TString(".root"));
          
      
    
    }
    
    
    //cbkg->SaveAs(TString("databkg") + catname + TString(".png"));
    //return;
    
    normv.at(icat)->setRange(0.0,normv.at(icat)->getVal()+8.0*sqrt(normv.at(icat)->getVal()));
    //normv.at(icat)->setConstant();
  }
  
//   for (UInt_t i=0; i<coeffv.size(); ++i) {
//     coeffv.at(i)->setRange(0.0,20.0);
//   }
    
  
    
  w->import(datacomb);  
  w->import(fullbkgpdf);
  w->import(*fullbkgfitres);
  w->import(*combcatpdf,RecycleConflictNodes());  
  w->Print();
  w->writeToFile("bkgdatawithfit.root") ;  
  
  printf("IntLumi = %5f\n",IntLumi->getVal());
  printf("ndata:\n");
  for (UInt_t icat=0; icat<catnames.size(); ++icat) {    
    //printf("%u: ndata = %i\n",icat,datav.at(icat)->numEntries());
    printf("%i ",datav.at(icat)->numEntries());
    
  }   
  printf("\n");
  
}


