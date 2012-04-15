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
#include "TLine.h"
#include "RooPoisson.h"
#include "RooNLLVar.h"
#include "RooPolynomial.h"
#include "RooChebychev.h"
#include "TText.h"
#include "TPaveText.h"
#include "RooNDKeysPdf.h"
#include "RooArgusBG.h"
#include "RooLandau.h"
#include "RooLognormal.h"
#include "RooGamma.h"
#include "RooBernstein.h"
#include "TGraphAsymmErrors.h"
#include "RooRandom.h"
#include <exception>

using namespace RooFit;



void makedatacard() { 
  
    FILE *file = fopen ("testcard.txt","w");  

    TString wslocation = "/scratch/bendavid/root/bambuOutputFeb16smmvavbf";
    
    const double lumi = 4763.;    
    
    const double lumiuncert = 1.045;
    const double triguncert = 1.01;
    const double vtxuncert = 0.005;
    
    std::vector<TString> catnames;    
    catnames.push_back("cat0");
    catnames.push_back("cat1");
    catnames.push_back("cat2");
    catnames.push_back("cat3");
    catnames.push_back("cat4");    
 
    std::vector<TString> procnames;
    procnames.push_back("ggH");
    procnames.push_back("qqH");
    procnames.push_back("VH");
    procnames.push_back("ttH");
    procnames.push_back("bkg_mass");
    
    const int ncats = catnames.size();
    
    printf("intro\n");
    
    //introductory stuff
    fprintf(file,"CMS-HGG DataCard for Unbinned Limit Setting\n");
    fprintf(file,"Run with: combine\n");
    fprintf(file,"---------------------------------------------\n");
    fprintf(file, "imax %i\n", ncats);
    fprintf(file, "jmax *\n");
    fprintf(file, "kmax *\n");
    fprintf(file,"---------------------------------------------\n");    
    
    fprintf(file,"\n\n\n\n");

    printf("done intro\n");
    
    //data
    fprintf(file, "shapes data_obs * %s/databkg/bkgdatawithfit.root wbkg:data_$CHANNEL\n", wslocation.Data());
    
    //bkg pdfs
    fprintf(file, "shapes bkg_mass * %s/databkg/bkgdatawithfit.root wbkg:CMS_hgg_$CHANNEL_bkgshape\n", wslocation.Data());
    
    //signal pdfs
    fprintf(file, "shapes ggH * %s/ggh/ubersignalmodel.root wsig:hggpdfsmrel_$CHANNEL_ggh\n", wslocation.Data());
    fprintf(file, "shapes qqH * %s/vbf/ubersignalmodel.root wsig:hggpdfsmrel_$CHANNEL_vbf\n", wslocation.Data());
    fprintf(file, "shapes VH  * %s/wzh/ubersignalmodel.root wsig:hggpdfsmrel_$CHANNEL_wzh\n", wslocation.Data());
    fprintf(file, "shapes ttH * %s/tth/ubersignalmodel.root wsig:hggpdfsmrel_$CHANNEL_tth", wslocation.Data());

    fprintf(file,"\n\n\n\n");
    
    
    //bin and observation lines
    fprintf(file, "bin         ");
    for (int icat=0; icat<ncats; ++icat) {
      fprintf(file, "%s ",catnames.at(icat).Data());
    }
    fprintf(file,"\n");
    
    fprintf(file, "observation ");
    for (int icat=0; icat<ncats; ++icat) {
      fprintf(file, "-1 ");
    }    
    fprintf(file,"\n");
    
    //rates
    fprintf(file,"bin     ");
    for (int icat=0; icat<ncats; ++icat) {
      for (int iproc=0; iproc<procnames.size(); ++iproc) {
        fprintf(file,"%s ", catnames.at(icat).Data());
      }
    }
    fprintf(file,"\n");

    fprintf(file,"process ");
    for (int icat=0; icat<ncats; ++icat) {
      for (int iproc=0; iproc<procnames.size(); ++iproc) {
        fprintf(file,"%s ", procnames.at(iproc).Data());
      }
    }
    fprintf(file,"\n");    
 
    fprintf(file,"process ");
    for (int icat=0; icat<ncats; ++icat) {
      for (int iproc=0; iproc<procnames.size(); ++iproc) {
	int idx;
	if (iproc==(procnames.size()-1)) idx = 1;
	else idx = -iproc;
        fprintf(file, "%i ", idx);
      }
    }
    fprintf(file,"\n");    
    
    fprintf(file,"rate ");
    for (int icat=0; icat<ncats; ++icat) {
      for (int iproc=0; iproc<procnames.size(); ++iproc) {
	float rate;
	if (iproc==(procnames.size()-1)) rate = 1.0;
	else rate = lumi;
        fprintf(file, "%3f ", rate);
      }
    }
    fprintf(file,"\n");      
 
    fprintf(file,"\n\n\n\n");
    

    printf("done rates\n");
    
    //theory systematics (overall normalization)
    fprintf(file, "QCDscale_ggH         lnN  ");
    for (int icat=0; icat<ncats; ++icat) {
      fprintf(file, "0.918/1.125 - - - - ");
    }
    fprintf(file,"\n");      
    
    fprintf(file, "PDF_gg               lnN  ");
    for (int icat=0; icat<ncats; ++icat) {
      fprintf(file, "0.923/1.079 - - 0.915/1.085 - ");
    }
    fprintf(file,"\n");      
    
    fprintf(file, "QCDscale_qqH         lnN  ");
    for (int icat=0; icat<ncats; ++icat) {
      fprintf(file, "- 0.997/1.005 - - - ");
    }
    fprintf(file,"\n");      
    
    fprintf(file, "PDF_qqbar            lnN  ");
    for (int icat=0; icat<ncats; ++icat) {
      fprintf(file, "- 0.979/1.027 0.958/1.042 - - ");
    }
    fprintf(file,"\n");      

    fprintf(file, "QCDscale_VH          lnN  ");
    for (int icat=0; icat<ncats; ++icat) {
      fprintf(file, "- - 0.982/1.018 - - ");
    }
    fprintf(file,"\n");      
    
    fprintf(file, "QCDscale_ttH         lnN  ");
    for (int icat=0; icat<ncats; ++icat) {
      fprintf(file, "- - - 0.905/1.036 - ");
    }
    fprintf(file,"\n");      
    
    fprintf(file,"\n\n\n\n");

    //lumi uncertainty
    fprintf(file,"lumi                  lnN  ");
    for (int icat=0; icat<ncats; ++icat) {
      for (int iproc=0; iproc<procnames.size(); ++iproc) {
	float rate;
	if (iproc==(procnames.size()-1)) fprintf(file, "%3f ", lumiuncert);
	else fprintf(file, "- ");      
      }
    }
    fprintf(file,"\n");     
    
    //trigger efficiency uncertainty
    fprintf(file,"CMS_hgg_eff_trig lnN  ");
    for (int icat=0; icat<ncats; ++icat) {
      for (int iproc=0; iproc<procnames.size(); ++iproc) {
	float rate;
	if (iproc==(procnames.size()-1)) fprintf(file, "%3f ", triguncert);
	else fprintf(file, "- ");      
      }
    }
    fprintf(file,"\n");   
    
    //vertex selection fraction uncertainty
    fprintf(file, "CMS_hgg_nuissancedeltafracright param 1.0 %3f\n", vtxuncert);
    
    
    //load workspace
    TFile *fws = new TFile("/scratch/bendavid/root/bambuOutputFeb16smmvavbf/CMS-HGG-mclow.root","READ");
    RooWorkspace *ws = (RooWorkspace*)fws->FindObjectAny("cms_hgg_workspace_mclow");    
    
    std::vector<RooAbsData*> dsets;
    for (int icat=0; icat<ncats; ++icat) {
      RooDataSet *dcat = (RooDataSet*)ws->data(TString::Format("sig_mass_m120_%s",catnames.at(icat).Data()));
      dsets.push_back(dcat);
    }
        
    //load mva application file to do photon id and per-event resolution estimate systematics
    TFile *mvafile = new TFile("/scratch/bendavid/root/hggmvaFeb16/appOutput_SM_Feb16.root","READ");
    TCut kinsel = "pt1>(mass/3.0) && pt2>(mass/4.0) && pt1>(100.0/3.0) && pt2>(100.0/4.0) && mass > 100. && mass < 180.";
    TCut idsel = "idmva1>-0.3 && idmva2>-0.3";
    TCut idselup = "idmva1>-0.325 && idmva2>-0.325";
    TCut idseldown = "idmva1>-0.275 && idmva2>-0.275";    
    TCut msel = kinsel && idsel;
    TTree* tup = (TTree*) mvafile->FindObjectAny("MITMVAtuple");

    std::vector<TCut> catcuts;
    std::vector<TCut> catcutsup;
    std::vector<TCut> catcutsdown;
    std::vector<TCut> catcutsidup;
    std::vector<TCut> catcutsiddown;
    
    TCut vbfcut = "(jet1pt>30.0 && jet2pt>20.0 && abs(jet1eta-jet2eta)>3.5 && dijetmass>350.0 && zeppenfeld<2.5 && abs(dphidijetgg)>2.6 && pt1>(55.0*mass/120.0))";
    
    //feb 16 smmvavbf
    TCut cat0 = msel && !vbfcut && "mva>= 0.89";
    TCut cat1 = msel && !vbfcut && "mva>= 0.72 && mva<0.89";
    TCut cat2 = msel && !vbfcut && "mva>= 0.55 && mva<0.72";
    TCut cat3 = msel && !vbfcut && "mva>= 0.05 && mva<0.55";    
    TCut cat4 = msel &&  vbfcut && "mva>= 0.05";
    catcuts.push_back(cat0);
    catcuts.push_back(cat1);
    catcuts.push_back(cat2);
    catcuts.push_back(cat3);
    catcuts.push_back(cat4);    
    
    TCut cat0up = msel && !vbfcut && "mvaup>= 0.89";
    TCut cat1up = msel && !vbfcut && "mvaup>= 0.72 && mvaup<0.89";
    TCut cat2up = msel && !vbfcut && "mvaup>= 0.55 && mvaup<0.72";
    TCut cat3up = msel && !vbfcut && "mvaup>= 0.05 && mvaup<0.55";    
    TCut cat4up = msel &&  vbfcut && "mvaup>= 0.05";
    catcutsup.push_back(cat0up);
    catcutsup.push_back(cat1up);
    catcutsup.push_back(cat2up);
    catcutsup.push_back(cat3up);
    catcutsup.push_back(cat4up);
    

    TCut cat0down = msel && !vbfcut && "mvadown>= 0.89";
    TCut cat1down = msel && !vbfcut && "mvadown>= 0.72 && mvadown<0.89";
    TCut cat2down = msel && !vbfcut && "mvadown>= 0.55 && mvadown<0.72";
    TCut cat3down = msel && !vbfcut && "mvadown>= 0.05 && mvadown<0.55";    
    TCut cat4down = msel &&  vbfcut && "mvadown>= 0.05";
    catcutsdown.push_back(cat0down);
    catcutsdown.push_back(cat1down);
    catcutsdown.push_back(cat2down);
    catcutsdown.push_back(cat3down);
    catcutsdown.push_back(cat4down);
    
    TCut cat0idup = kinsel && idselup && !vbfcut && "mvaidup>= 0.89";
    TCut cat1idup = kinsel && idselup && !vbfcut && "mvaidup>= 0.72 && mvaidup<0.89";
    TCut cat2idup = kinsel && idselup && !vbfcut && "mvaidup>= 0.55 && mvaidup<0.72";
    TCut cat3idup = kinsel && idselup && !vbfcut && "mvaidup>= 0.05 && mvaidup<0.55";    
    TCut cat4idup = kinsel && idselup &&  vbfcut && "mvaidup>= 0.05";    
    catcutsidup.push_back(cat0idup);
    catcutsidup.push_back(cat1idup);
    catcutsidup.push_back(cat2idup);
    catcutsidup.push_back(cat3idup);
    catcutsidup.push_back(cat4idup);

    TCut cat0iddown = kinsel && idseldown && !vbfcut && "mvaiddown>= 0.89";
    TCut cat1iddown = kinsel && idseldown && !vbfcut && "mvaiddown>= 0.72 && mvaiddown<0.89";
    TCut cat2iddown = kinsel && idseldown && !vbfcut && "mvaiddown>= 0.55 && mvaiddown<0.72";
    TCut cat3iddown = kinsel && idseldown && !vbfcut && "mvaiddown>= 0.05 && mvaiddown<0.55";    
    TCut cat4iddown = kinsel && idseldown &&  vbfcut && "mvaiddown>= 0.05";
    catcutsiddown.push_back(cat0iddown);
    catcutsiddown.push_back(cat1iddown);
    catcutsiddown.push_back(cat2iddown);
    catcutsiddown.push_back(cat3iddown);
    catcutsiddown.push_back(cat4iddown);

    std::vector<TCut> procuts;
    procuts.push_back("proc==0");
    procuts.push_back("proc==1");
    procuts.push_back("proc==2");
    procuts.push_back("proc==3");

    TH1D *hmmva = new TH1D("hmmva","",100,-1.0,1.0);
    TH1D *hmmvaup = new TH1D("hmmvaup","",100,-1.0,1.0);
    TH1D *hmmvadown = new TH1D("hmmvadown","",100,-1.0,1.0);
    TCut mvasel = msel && "proc>0 && proc<3";
    tup->Draw("mva>>hmmva",mvasel);
    tup->Draw("mvaup>>hmmvaup",mvasel);
    tup->Draw("mvadown>>hmmvadown",mvasel);
    
    hmmvaup->SetLineColor(kBlue);
    hmmvadown->SetLineColor(kRed);
    new TCanvas;
    hmmva->Draw("HIST");
    hmmvaup->Draw("HISTSAME");
    hmmvadown->Draw("HISTSAME");
    //return;
    
    new TCanvas;
    
    TH1D *hmasscount = new TH1D("hmasscount","",160,100.,180.);
    
    printf("sigmae loop\n");
    
    fprintf(file, "CMS_hgg_n_sigmae lnN ");
    for (int i=0; i<catcuts.size(); ++i) {
      printf("cat %i\n",i);
      for (int j=0; j<procuts.size(); ++j) {
	printf("proc %i\n",j);
        TCut catcut = catcuts[i] && procuts[j];
        TCut catcutup = catcutsup[i] && procuts[j];
        TCut catcutdown = catcutsdown[i] && procuts[j];        
        
        hmasscount->Reset();
        tup->Draw("mass>>hmasscount",catcut);
        double nom = hmasscount->GetSumOfWeights();
        
        hmasscount->Reset();
        tup->Draw("mass>>hmasscount",catcutup);
        double up = hmasscount->GetSumOfWeights();
        
        hmasscount->Reset();
        tup->Draw("mass>>hmasscount",catcutdown);
        double down = hmasscount->GetSumOfWeights();        
        
        fprintf(file, "%3f/%3f ",up/nom,down/nom);
        
      }
      fprintf(file,"- ");
    }
    fprintf(file, "\n");
    
    
    printf("id mva loop\n");
    fprintf(file, "CMS_hgg_n_id lnN ");
    for (int i=0; i<catcuts.size(); ++i) {
      printf("cat %i\n",i);
      for (int j=0; j<procuts.size(); ++j) {
        printf("proc %i\n",j);
        TCut catcut = catcuts[i] && procuts[j];
        TCut catcutidup = catcutsidup[i] && procuts[j];
        TCut catcutiddown = catcutsiddown[i] && procuts[j];        
        
        hmasscount->Reset();
        tup->Draw("mass>>hmasscount",catcut);
        double nom = hmasscount->GetSumOfWeights();
        
        hmasscount->Reset();
        tup->Draw("mass>>hmasscount",catcutidup);
        double idup = hmasscount->GetSumOfWeights();
        
        hmasscount->Reset();
        tup->Draw("mass>>hmasscount",catcutiddown);
        double iddown = hmasscount->GetSumOfWeights();        
        
        fprintf(file, "%3f/%3f ",idup/nom,iddown/nom);
        
      }
      fprintf(file, "- ");
    }
    fprintf(file, "\n");    
 
    fprintf(file,"\n\n\n\n");
    
    
    printf("done shape uncertainties\n");
    
    //energy scale/smearing single photon categories
    TCut ph1cut1 = "ph1.phcat==1 && abs(ph1.sceta)<=1.0";
    TCut ph1cut2 = "ph1.phcat==2 && abs(ph1.sceta)<=1.0";
    TCut ph1cut3 = "ph1.phcat==1 && abs(ph1.sceta)>1.0";
    TCut ph1cut4 = "ph1.phcat==2 && abs(ph1.sceta)>1.0";
    TCut ph1cut5 = "ph1.phcat==3 && abs(ph1.sceta)<=2.0";
    TCut ph1cut6 = "ph1.phcat==4 && abs(ph1.sceta)<=2.0";    
    TCut ph1cut7 = "ph1.phcat==3 && abs(ph1.sceta)>2.0";
    TCut ph1cut8 = "ph1.phcat==4 && abs(ph1.sceta)>2.0";

    TCut ph2cut1 = "ph2.phcat==1 && abs(ph2.sceta)<=1.0";
    TCut ph2cut2 = "ph2.phcat==2 && abs(ph2.sceta)<=1.0";    
    TCut ph2cut3 = "ph2.phcat==1 && abs(ph2.sceta)>1.0";
    TCut ph2cut4 = "ph2.phcat==2 && abs(ph2.sceta)>1.0";
    TCut ph2cut5 = "ph2.phcat==3 && abs(ph2.sceta)<=2.0";
    TCut ph2cut6 = "ph2.phcat==4 && abs(ph2.sceta)<=2.0";    
    TCut ph2cut7 = "ph2.phcat==3 && abs(ph2.sceta)>2.0";
    TCut ph2cut8 = "ph2.phcat==4 && abs(ph2.sceta)>2.0";
    
    std::vector<TCut> ph1cuts;
    std::vector<TCut> ph2cuts;
    
    ph1cuts.push_back(ph1cut1);
    ph1cuts.push_back(ph1cut2);
    ph1cuts.push_back(ph1cut3);
    ph1cuts.push_back(ph1cut4);
    ph1cuts.push_back(ph1cut5);
    ph1cuts.push_back(ph1cut6);
    ph1cuts.push_back(ph1cut7);
    ph1cuts.push_back(ph1cut8);

    ph2cuts.push_back(ph2cut1);
    ph2cuts.push_back(ph2cut2);
    ph2cuts.push_back(ph2cut3);
    ph2cuts.push_back(ph2cut4);      
    ph2cuts.push_back(ph2cut5);
    ph2cuts.push_back(ph2cut6);
    ph2cuts.push_back(ph2cut7);
    ph2cuts.push_back(ph2cut8);      

    //single photon smearing numbers
    std::vector<double> smearsingles;   
    smearsingles.push_back(0.0045);
    smearsingles.push_back(0.0109);    
    smearsingles.push_back(0.0156);
    smearsingles.push_back(0.0203);   
    smearsingles.push_back(0.0303);
    smearsingles.push_back(0.0326);    
    smearsingles.push_back(0.0318);
    smearsingles.push_back(0.0331);  
    
    //single photon smearing uncertainty
    std::vector<double> dsmearsingles;
    dsmearsingles.push_back(0.22e-2);
    dsmearsingles.push_back(0.24e-2);    
    dsmearsingles.push_back(0.60e-2);
    dsmearsingles.push_back(0.59e-2);
    dsmearsingles.push_back(0.90e-2);
    dsmearsingles.push_back(0.30e-2);    
    dsmearsingles.push_back(0.34e-2);
    dsmearsingles.push_back(0.52e-2);
    
    //single photon energy scale (approximated as 1.0 for calculating systematics)
    std::vector<double> scalesingles;
    scalesingles.push_back(1.0);
    scalesingles.push_back(1.0);
    scalesingles.push_back(1.0);
    scalesingles.push_back(1.0);
    scalesingles.push_back(1.0);
    scalesingles.push_back(1.0);
    scalesingles.push_back(1.0);
    scalesingles.push_back(1.0);
    
    //single photon energy scale uncertainties
    std::vector<double> dscalesingles;
    dscalesingles.push_back(0.19e-2);
    dscalesingles.push_back(0.13e-2);    
    dscalesingles.push_back(0.71e-2);
    dscalesingles.push_back(0.51e-2);    
    dscalesingles.push_back(0.88e-2);
    dscalesingles.push_back(0.18e-2);    
    dscalesingles.push_back(0.19e-2);
    dscalesingles.push_back(0.28e-2);        

    std::vector<double> smeardoubles;
    std::vector<double> dsmeardoubles;
    std::vector<double> dscaledoubles;
    
    for (int icat=0; icat<dsets.size(); ++icat) {
      RooAbsData *dset = dsets.at(icat);
      double ntot = dset->sumEntries();
      double smear = 0.0;
      double dsmearsq = 0.0;
      double scale = 0.0;
      double dscalesq = 0.0;
      
      for (int i=0; i<ph1cuts.size(); ++i) {
        TCut ph1cuti = ph1cuts.at(i);
        TCut ph2cuti = ph2cuts.at(i);

        double smeari = smearsingles.at(i);
        double dsmeari = dsmearsingles.at(i);
        double scalei = scalesingles.at(i);
        double dscalei = dscalesingles.at(i);
        
        double dsmearpartial = 0.0;
        double dscalepartial = 0.0;
        for (int j=0; j<ph1cuts.size(); ++j) {
          TCut ph1cutj = ph1cuts.at(j);
          TCut ph2cutj = ph2cuts.at(j);          
          
          double smearj = smearsingles.at(j);
          double scalej = scalesingles.at(j);
          
          TCut subcut = (ph1cuti && ph2cutj) || (ph1cutj && ph2cuti);
          double frac = dset->sumEntries(subcut)/ntot;
          
          //printf("icat %i, i %i j %i frac %5f\n",icat,i,j,frac);
          
          if (j>=i) {
            smear += 0.5*frac*sqrt(smeari*smeari + smearj*smearj);
            scale += frac*sqrt(scalei*scalej);
          }
          
          if (j==i) {
            //smear += 0.5*sqrt(2.0)*frac*smeari;
            dsmearpartial += 0.5*sqrt(2.0)*frac;
            //scale += frac*scalei;
            dscalepartial += frac;
          }
          else {
            //smear += 0.5*frac*sqrt(smeari*smeari+smearj*smearj);
            dsmearpartial += 0.5*frac*smeari/sqrt(smeari*smeari+smearj*smearj);
            //scale += frac*sqrt(scalei*scalej);
            dscalepartial += frac*scalej/sqrt(scalei*scalej)/2.0;
          }
          
        }
        dsmearsq += dsmearpartial*dsmearpartial*dsmeari*dsmeari;
        dscalesq += dscalepartial*dscalepartial*dscalei*dscalei;
        
      }
      
      //printf("smear = %5f +- %5f\n",smear,sqrt(dsmearsq));
      //printf("scale = %5f +- %5f\n",scale,sqrt(dscalesq));
      
      smeardoubles.push_back(smear);
      dsmeardoubles.push_back(sqrt(dsmearsq));
      dscaledoubles.push_back(sqrt(dscalesq));
    }
    
    for (int icat=0; icat<ncats; ++icat) {
      fprintf(file, "CMS_hgg_nuissancedeltasmear%s param 0.0 %5f\n" , catnames.at(icat).Data(), dsmeardoubles.at(icat));
      printf("%s: smear = %5f +- %5f\n",catnames.at(icat).Data(),smeardoubles.at(icat),dsmeardoubles.at(icat));
    }
    
    fprintf(file,"\n\n\n\n");
    
    
    for (int icat=0; icat<ncats; ++icat) {
      fprintf(file, "CMS_hgg_nuissancedeltam%s param 0.0 %5f\n" , catnames.at(icat).Data(), dscaledoubles.at(icat));
      printf("%s: scale = 1.0 +- %5f\n",catnames.at(icat).Data(),dscaledoubles.at(icat));
    }    

    fprintf(file,"\n\n\n\n");
    
    
    //photon efficiency uncertainties
    std::vector<double> deffsinglesr;
    deffsinglesr.push_back(1.010);
    deffsinglesr.push_back(1.026);
    
    std::vector<TCut> singlecutsr;
    singlecutsr.push_back("(abs(ph1.sceta)<1.5 || abs(ph2.sceta)<1.5) && !(abs(ph1.sceta)<1.5 && abs(ph2.sceta)<1.5)");
    singlecutsr.push_back("(abs(ph1.sceta)>1.5 || abs(ph2.sceta)>1.5) && !(abs(ph1.sceta)>1.5 && abs(ph2.sceta)>1.5)");

    std::vector<TCut> doublecutsr;
    doublecutsr.push_back("(abs(ph1.sceta)<1.5 && abs(ph2.sceta)<1.5)");
    doublecutsr.push_back("(abs(ph1.sceta)>1.5 && abs(ph2.sceta)>1.5)");

    std::vector<TCut> procutsr;
    procutsr.push_back("procidx==0");
    procutsr.push_back("procidx==1");
    procutsr.push_back("procidx==2");
    procutsr.push_back("procidx==3");    
    
    printf("starting 2eff loop\n");
    
    for (int j=0; j<deffsinglesr.size(); ++j) {
      std::vector<double> catuncerts;
      double uncert = deffsinglesr.at(j);
      fprintf(file, "CMS_eff_g_%i      lnN  ",j);
      for (int icat=0; icat<dsets.size(); ++icat) {
        RooAbsData *dset = dsets.at(icat);
	for (int k=0; k<procutsr.size(); ++k) {
	  double ntot = dset->sumEntries(procutsr.at(k));
	  double singlefrac = dset->sumEntries(singlecutsr.at(j)&&procutsr.at(k))/ntot;
	  double doublefrac = dset->sumEntries(doublecutsr.at(j)&&procutsr.at(k))/ntot;
	  double catuncert = 1.0 + singlefrac*(uncert-1.0) + doublefrac*(uncert*uncert-1.0);
	  fprintf(file, "%3f ",catuncert);
	}
        //printf("uncert = %5f,  ntot = %5f, singlefrac = %5f, doublefrac = %5f\n",uncert,ntot,singlefrac,doublefrac);
        //printf("%3f %3f %3f %3f - ",catuncert,catuncert,catuncert,catuncert);
        fprintf(file, "- ");
      }      
      fprintf(file, "\n");
      
    }
 
   fprintf(file,"\n\n\n\n");

 
    printf("ending 2eff loop\n");
 
 
    //vbf uncertainties
    std::vector<double> vbfgguncerts;
    std::vector<double> vbfuncerts;
    
    
    vbfgguncerts.push_back(1.20);
    vbfgguncerts.push_back(1.15);
    vbfgguncerts.push_back(1.70);
    
    vbfuncerts.push_back(1.06);
    vbfuncerts.push_back(1.04);
    vbfuncerts.push_back(1.08);    
   
    std::vector<TCut> mvacuts;
    mvacuts.push_back("bdt>=0.89");
    mvacuts.push_back("bdt>=0.72 && bdt<0.89");
    mvacuts.push_back("bdt>=0.55 && bdt<0.72");
    mvacuts.push_back("bdt>=0.05 && bdt<0.55");
      
    RooAbsData *tagdset = dsets.back();
    RooAbsData *tagdset2 = 0;
    //RooAbsData *tagdset2 = dsets.at(dsets.size()-2);
    for (int j=0; j<vbfuncerts.size(); ++j) {
      double vbfuncert = vbfuncerts.at(j);
      double gguncert = vbfgguncerts.at(j);      
      fprintf(file, "CMS_hgg_vbf_%i      lnN  ",j);
      for (int icat=0; (icat<dsets.size()-1); ++icat) {
        RooAbsData *dset = dsets.at(icat);
	for (int k=0; k<procutsr.size(); ++k) {
	  double uncert;
	  if (k==1) uncert = vbfuncert;
	  else uncert = gguncert;
	  double ntot = dset->sumEntries(procutsr.at(k));
	  double ntagged = tagdset->sumEntries(procutsr.at(k)&&mvacuts.at(icat));
          if (tagdset2) ntagged += tagdset2->sumEntries(procutsr.at(k)&&mvacuts.at(icat));
	  double catuncert = (ntot - (uncert-1.0)*ntagged)/ntot;
	  fprintf(file, "%3f ",catuncert);
	}
        //printf("uncert = %5f,  ntot = %5f, singlefrac = %5f, doublefrac = %5f\n",uncert,ntot,singlefrac,doublefrac);
        //printf("%3f %3f %3f %3f - ",catuncert,catuncert,catuncert,catuncert);
        fprintf(file, "- ");
      }
      fprintf(file, "%3f %3f %3f %3f - ",gguncert,vbfuncert,gguncert,gguncert);
      //printf("%3f %3f - ",vbfuncert,gguncert);
      fprintf(file, "\n");
      
    }
    return;

  
}


