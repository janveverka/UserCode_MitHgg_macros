#include <sstream>
#include <iostream>
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

using namespace RooFit;

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


void fithmasssgglobesplit5(int filterproc = -1) { 
  
  
  gROOT->Macro("MitStyle.C");
  gStyle->SetErrorX(0); 
  gStyle->SetOptStat(0);
  gROOT->ForceStyle();  
  
  

  //load pt-weights
  TFile *fileptweight = new TFile("/home/bendavid/cms/root/pudists/Kfactors.root","READ");*/
  std::vector<TString> procnames;
  procnames.push_back("ggh");
  procnames.push_back("vbf");
  procnames.push_back("wzh");
  procnames.push_back("tth");

  
  gSystem->cd(TString::Format("./bambuOutputMar19smmet/%s",procnames.at(filterproc).Data())); 
  gStyle->SetOptStat(1110);
  
 
  
  
  TFile *fdata = new TFile("/scratch/bendavid/root/bambuOutputMar19smmet/CMS-HGG-mclow.root","READ");
  RooWorkspace* win1 = (RooWorkspace*)fdata->Get("cms_hgg_workspace_mclow");

  TFile *fdata2 = new TFile("/scratch/bendavid/root/bambuOutputMar19smmet/CMS-HGG-mchigh.root","READ");
  RooWorkspace* win2 = (RooWorkspace*)fdata2->Get("cms_hgg_workspace_mchigh");  
  
  RooWorkspace *win = win1;
  
  RooWorkspace *w = new RooWorkspace("wsig","") ;

  const double massmax = 180.0;
  RooRealVar *hmass = win->var("CMS_hgg_mass");
  //RooRealVar *hmass = win->var("mass");
  hmass->setRange(100,massmax);
  hmass->setBins(4.0*(massmax-100.0));
  hmass->SetTitle("m_{#gamma#gamma}");
  hmass->setUnit("GeV");  
  

  
//   win->Print("V");
//   return;
    

  std::vector<double> smearingv;
  //jan27/feb2/feb16 smmvavbf
  smearingv.push_back(0.005432);
  smearingv.push_back(0.005196);
  smearingv.push_back(0.007464);
  smearingv.push_back(0.012978);
  smearingv.push_back(0.008722);      

  
    
  std::vector<TString> catnamesbase;    
  catnamesbase.push_back("cat0");
  catnamesbase.push_back("cat1");
  catnamesbase.push_back("cat2");
  catnamesbase.push_back("cat3");
  catnamesbase.push_back("cat4");
//   catnamesbase.push_back("cat5");
//   catnamesbase.push_back("cat6");  
  
  //define categories
  std::vector<TString> catnames;
  std::vector<TCut> catcuts;
  

  
  catnames.push_back(catnamesbase.at(0)+TString("_")+procnames.at(filterproc));
  catnames.push_back(catnamesbase.at(1)+TString("_")+procnames.at(filterproc));
  catnames.push_back(catnamesbase.at(2)+TString("_")+procnames.at(filterproc));
  catnames.push_back(catnamesbase.at(3)+TString("_")+procnames.at(filterproc));
  catnames.push_back(catnamesbase.at(4)+TString("_")+procnames.at(filterproc));
//   catnames.push_back(catnamesbase.at(5)+TString("_")+procnames.at(filterproc));
//   catnames.push_back(catnamesbase.at(6)+TString("_")+procnames.at(filterproc));

  
  std::set<TString> complexset;
//  complexset.insert(catnamesbase.at(0)+TString("_")+procnames.at(filterproc));
  complexset.insert(catnamesbase.at(1)+TString("_")+procnames.at(filterproc));
  complexset.insert(catnamesbase.at(2)+TString("_")+procnames.at(filterproc));
  complexset.insert(catnamesbase.at(3)+TString("_")+procnames.at(filterproc));
//   complexset.insert(catnamesbase.at(4)+TString("_")+procnames.at(filterproc));
  //complexset.insert(catnamesbase.at(5)+TString("_")+procnames.at(filterproc));
  //complexset.insert(catnamesbase.at(6)+TString("_")+procnames.at(filterproc));


  
  std::set<TString> simpleset;
  //simpleset.insert("cat6");
  //simpleset.insert("cat7");

  
  std::vector<RooRealVar*> fitparms;
  std::vector<double> fitparmsinit;

  
  //RooRealVar mnom("mnom","m_{h}",110.0,100.0,200.0,"GeV");
  RooRealVar mnom("MH","m_{h}",110.0,100.0,200.0,"GeV");
  mnom.setConstant();
 
  //mnom.setRange("plotrange",110,140);
  
  //define signal pdf variables
  
  const double dm1init = 0.2;
  RooRealVar dm1("dm1","",dm1init,-1.0,1.0);
  //dm1.removeRange();
  
  const double dm2init = -1.0;
  RooRealVar dm2("dm2","",dm2init,-5.0,0.0);
  dm2.removeRange();
  
  const double dm3init = -2.0;
  RooRealVar dm3("dm3","",dm3init,-9.0,0.0);
  dm3.removeRange();  
  
  RooFormulaVar mean1("mean1","","@0+@1",RooArgList(mnom,dm1));
  RooFormulaVar mean2("mean2","","@0+@1",RooArgList(mnom,dm2));
  RooFormulaVar mean3("mean3","","@0+@1",RooArgList(mnom,dm3));
  
  const double sigma1init = 0.5;
  RooRealVar sigma1("sigma1","",sigma1init,0.4,5.0);
  sigma1.removeRange();
  
  const double sigma2init = 0.8;
  RooRealVar sigma2("sigma2","",sigma2init,0.8,7.0);
  sigma2.removeRange();
  
  const double sigma3init = 2.0;
  RooRealVar sigma3("sigma3","",sigma3init,1.0,10.0);
  sigma3.removeRange();  
  
  const double f1init = 0.60;
  RooRealVar f1("f1","",f1init,0.52,1.0);
  
  const double f2init = 0.90;
  RooRealVar f2("f2","",f2init,0.0,1.0);
  
  RooGaussian g1("g1","",*hmass,mean1,sigma1);
  RooGaussian g2("g2","",*hmass,mean2,sigma2);
  RooGaussian g3("g3","",*hmass,mean3,sigma3);
  
  RooAddPdf combh("combh","",RooArgList(g1,g2,g3),RooArgList(f1,f2),kTRUE);
  RooAddPdf combhmin("combhmin","",RooArgList(g1,g2),RooArgList(f1),kTRUE);
 // RooGaussian &combh = g1;
  
  const double wdm1init = 0.5;
  RooRealVar wdm1("wdm1","",wdm1init,-12.0,5.0);
  wdm1.removeRange();
  
  const double wdm2init = -0.7;
  RooRealVar wdm2("wdm2","",wdm2init,-12.0,5.0);
  wdm2.removeRange();
  
  
  RooFormulaVar wmean1("wmean1","","@0+@1",RooArgList(mnom,wdm1));
  RooFormulaVar wmean2("wmean2","","@0+@1",RooArgList(mnom,wdm2));
  
  const double wsigma1init = 4.0;
  RooRealVar wsigma1("wsigma1","",wsigma1init,0.0,10.0); //2.0
  wsigma1.removeRange();

  const double wsigma2init = 2.0;  
  RooRealVar wsigma2("wsigma2","",wsigma2init,0.0,10.0); //3.0
  wsigma2.removeRange();
  
  const double wf1init = 0.6;
  RooRealVar wf1("wf1","",wf1init,0.52,1.0);
  
  RooGaussian wg1("wg1","",*hmass,wmean1,wsigma1);
  RooGaussian wg2("wg2","",*hmass,wmean2,wsigma2);
  
  RooAddPdf combhwrong("combhwrong","",RooArgList(wg1,wg2),RooArgList(wf1),kTRUE);
  //RooGaussian &combhwrong = wg1;

  
     
  const double fracrightinit = 0.9;
  RooRealVar fracright("fracright","fracright",fracrightinit,0.0,1.0);

  const double effaccinit = 0.9;
  RooRealVar effacc("effacc","effacc",effaccinit,0.0,1.0);  
  
  //signal pdf
  
  RooAddPdf combhvtx("combhvtx","combhvtx",RooArgList(combh,combhwrong),RooArgList(fracright));
  RooAddPdf combhvtxmin("combhvtxmin","combhvtxmin",RooArgList(combhmin,wg1),RooArgList(fracright));
  RooAddPdf combhvtxsimple("combhvtxsimple","combhvtxsimple",RooArgList(g1,wg1),RooArgList(fracright));
  
  fitparms.push_back(&dm1);
  fitparms.push_back(&dm2);  
  fitparms.push_back(&dm3);
  fitparms.push_back(&sigma1);
  fitparms.push_back(&sigma2);
  fitparms.push_back(&sigma3);
  fitparms.push_back(&f1);
  fitparms.push_back(&f2);
  fitparms.push_back(&wdm1);
  fitparms.push_back(&wdm2);  
  fitparms.push_back(&wsigma1);
  fitparms.push_back(&wsigma2);
  fitparms.push_back(&wf1);
  fitparms.push_back(&fracright);
  fitparms.push_back(&effacc);

  fitparmsinit.push_back(dm1init);
  fitparmsinit.push_back(dm2init);  
  fitparmsinit.push_back(dm3init);
  fitparmsinit.push_back(sigma1init);
  fitparmsinit.push_back(sigma2init);
  fitparmsinit.push_back(sigma3init);
  fitparmsinit.push_back(f1init);
  fitparmsinit.push_back(f2init);
  fitparmsinit.push_back(wdm1init);
  fitparmsinit.push_back(wdm2init);  
  fitparmsinit.push_back(wsigma1init);
  fitparmsinit.push_back(wsigma2init);
  fitparmsinit.push_back(wf1init);
  fitparmsinit.push_back(fracrightinit);
  fitparmsinit.push_back(effaccinit);
    
  //define mass points and cross sections
  std::vector<int> mhs;
  std::vector<double> gfxsecs;
  std::vector<double> vbfxsecs;
  std::vector<double> vtthxsecs;
  std::vector<double> vhxsecs;
  std::vector<double> ffxsecs;


/*  mhs.push_back(90);
  gfxsecs.push_back(0.03625);
  vbfxsecs.push_back(0.00210);
  vtthxsecs.push_back(0.00334); 
  vhxsecs.push_back(0.00307);

  
  mhs.push_back(100);
  gfxsecs.push_back(0.03819);
  vbfxsecs.push_back(0.00246);
  vtthxsecs.push_back(0.00315);
  vhxsecs.push_back(0.00289);
*/
  

  mhs.push_back(110);
  gfxsecs.push_back(0.03908);
  vbfxsecs.push_back(0.00275);
  vtthxsecs.push_back(0.00290);  
  vhxsecs.push_back(0.00265);
  ffxsecs.push_back(0.08323692+0.08023015);


  mhs.push_back(115);
  gfxsecs.push_back(0.03862);
  vbfxsecs.push_back(0.00284);
  vtthxsecs.push_back(0.00272);  
  vhxsecs.push_back(0.00248);
  ffxsecs.push_back(0.0481518+0.042125595);
   
  
  mhs.push_back(120);
  gfxsecs.push_back(0.03742);
  vbfxsecs.push_back(0.00286);
  vtthxsecs.push_back(0.00251);
  vhxsecs.push_back(0.00229);
  ffxsecs.push_back(0.02927583+0.023436813);
  
  //mhs.push_back(123);
  //mhs.push_back(125);
  
  
  mhs.push_back(130);
  gfxsecs.push_back(0.03191);
  vbfxsecs.push_back(0.00261);
  vtthxsecs.push_back(0.00193);  
  vhxsecs.push_back(0.00176);
  ffxsecs.push_back(0.01224394+0.008260946);
  
  //mhs.push_back(135);

  
  mhs.push_back(140);
  gfxsecs.push_back(0.02353);
  vbfxsecs.push_back(0.00204);
  vtthxsecs.push_back(0.00129);    
  vhxsecs.push_back(0.00117);
  ffxsecs.push_back(0.005656604+0.0032417933);
  
  //mhs.push_back(145);
  mhs.push_back(150);

  
  RooRealVar *weight = new RooRealVar("weight","",1.0);
  weight->removeRange();
  
  RooRealVar *procidx = win->var("procidx");  
  
  
  RooRealVar *IntLumi = win->var("IntLumi");
  //RooConstVar *IntLumi = new RooConstVar("IntLumi","",951.0);

  const int numsmpoints = 81;
  double smmasses[numsmpoints] = {110.0,110.5,111.0,111.5,112.0,112.5,113.0,113.5,114.0,114.5,115.0,115.5,116.0,116.5,117.0,117.5,118.0,118.5,119.0,119.5,120.0,120.5,121.0,121.5,122.0,122.5,123.0,123.5,124.0,124.5,125.0,125.5,126.0,126.5,127.0,127.5,128.0,128.5,129.0,129.5,130.0,130.5,131.0,131.5,132.0,132.5,133.0,133.5,134.0,134.5,135.0,135.5,136.0,136.5,137.0,137.5,138.0,138.5,139.0,139.5,140.0,140.5,141.0,141.5,142.0,142.5,143.0,143.5,144.0,144.5,145.0,145.5,146.0,146.5,147.0,147.5,148.0,148.5,149.0,149.5,150.0};  
  double smxsbr[numsmpoints]  = {0.04474106,0.04476087,0.04458980,0.04461816,0.04465764,0.04445241,0.04450231,0.04452202,0.04432575,0.04435697,0.04417173,0.04419175,0.04402058,0.04382749,0.04365297,0.04367220,0.04349037,0.04330598,0.04314091,0.04297301,0.04277804,0.04260475,0.04226444,0.04208721,0.04190975,0.04154534,0.04122657,0.04106506,0.04072321,0.04038161,0.04006593,0.03972509,0.03940765,0.03909069,0.03877650,0.03846030,0.03797801,0.03766622,0.03721147,0.03689953,0.03645195,0.03600491,0.03556311,0.03514415,0.03470779,0.03427222,0.03370990,0.03328123,0.03290089,0.03235105,0.03180601,0.03141397,0.03087743,0.03034568,0.02981916,0.02933985,0.02882227,0.02830964,0.02782357,0.02733986,0.02672205,0.02624931,0.02578068,0.02525736,0.02473893,0.02414999,0.02356720,0.02307409,0.02258560,0.02202963,0.02147946,0.02101546,0.02055579,0.02003015,0.01950998,0.01893346,0.01836331,0.01786838,0.01737859,0.01683392,0.01629523};
  double ffxsbr[numsmpoints]   = {0.16346707,0.153639388,0.144474184,0.135978954,0.128011936,0.120580658,0.113698278,0.1072188,0.101204298,0.095574901,0.090277395,0.085366199,0.080723132,0.076390615,0.072348276,0.068536512,0.064981344,0.061626793,0.058491355,0.055518848,0.052712643,0.050108436,0.047635539,0.045289044,0.043111872,0.041038726,0.0390816,0.03724344,0.035498398,0.033842835,0.032293386,0.03081212,0.029395632,0.028062999,0.026825547,0.025612992,0.024479542,0.023414238,0.02239302,0.021406322,0.020504886,0.01962496,0.018787808,0.0179886462,0.0172365732,0.0165092732,0.0158259735,0.0151657605,0.014543112,0.0139496598,0.01338225,0.0128386206,0.012320084,0.0118226735,0.0113473074,0.010898992,0.0104629713,0.010045404,0.0096510688,0.0092662975,0.0088983973,0.0085471461,0.0082114608,0.007885017,0.007570577,0.0072704658,0.0069803288,0.0067006628,0.0064322028,0.00617446755,0.005927448,0.005686585,0.00545616,0.005232885,0.0050181952,0.00480948195,0.004609108,0.0044133128,0.0042256148,0.0040421231,0.0038650707};
  
  double smbr[numsmpoints] = {0.00197,0.00199,0.002,0.00202,0.00204,0.00205,0.00207,0.00209,0.0021,0.00212,0.00213,0.00215,0.00216,0.00217,0.00218,0.0022,0.00221,0.00222,0.00223,0.00224,0.00225,0.00226,0.00226,0.00227,0.00228,0.00228,0.00228,0.00229,0.00229,0.00229,0.00229,0.00229,0.00229,0.00229,0.00229,0.00229,0.00228,0.00228,0.00227,0.00227,0.00226,0.00225,0.00224,0.00223,0.00222,0.00221,0.00219,0.00218,0.00217,0.00215,0.00213,0.00212,0.0021,0.00208,0.00206,0.00204,0.00202,0.002,0.00198,0.00196,0.00193,0.00191,0.00189,0.001865,0.00184,0.00181,0.00178,0.001755,0.00173,0.0017,0.00167,0.001645,0.00162,0.00159,0.00156,0.001525,0.00149,0.00146,0.00143,0.001395,0.00136};
  double ffbr[numsmpoints] = {0.059540,0.056510,0.053660,0.050980,0.048460,0.046090,0.043860,0.041760,0.039780,0.037910,0.036150,0.034490,0.032920,0.031430,0.030030,0.028710,0.027460,0.026270,0.025150,0.024080,0.023070,0.022120,0.021210,0.020340,0.019520,0.018740,0.018000,0.017300,0.016630,0.015990,0.015380,0.014800,0.014240,0.013710,0.013210,0.012720,0.012260,0.011820,0.011400,0.010990,0.010610,0.010240,0.009880,0.009539,0.009212,0.008897,0.008595,0.008305,0.008026,0.007758,0.007500,0.007251,0.007012,0.006781,0.006558,0.006344,0.006137,0.005937,0.005744,0.005557,0.005377,0.005202,0.005034,0.004870,0.004711,0.004558,0.004409,0.004264,0.004124,0.003987,0.003855,0.003725,0.003600,0.003477,0.003358,0.003241,0.003128,0.003016,0.002908,0.002801,0.002697};
  //double sm4br[numsmpoints] = {1.72E-004,0.0001747,0.0001774,0.0001801,0.0001828,0.0001855,0.0001882,0.0001909,0.0001936,0.0001963,0.000199,0.0002017,0.0002044,0.0002071,0.0002098,0.0002125,0.0002152,0.0002179,0.0002206,0.0002233,2.26E-004,0.00022945,0.0002329,0.00023635,0.0002398,0.00024325,0.0002467,0.00025015,0.0002536,0.00025705,0.0002605,0.00026395,0.0002674,0.00027085,0.0002743,0.00027775,0.0002812,0.00028465,0.0002881,0.00029155,2.95E-004,0.0002993,0.0003036,0.0003079,0.0003122,0.0003165,0.0003208,0.0003251,0.0003294,0.0003337,0.000338,0.0003423,0.0003466,0.0003509,0.0003552,0.0003595,0.0003638,0.0003681,0.0003724,0.0003767,3.81E-004,0.0003858,0.0003906,0.0003954,0.0004002,0.000405,0.0004098,0.0004146,0.0004194,0.0004242,4.29E-004,0.0004335,0.000438,0.0004425,0.000447,0.0004515,0.000456,0.0004605,0.000465,0.0004695,4.74E-004};
  double sm4br[numsmpoints] = {4.40E-005,0.000043685,0.00004337,0.000043055,0.00004274,0.000042425,0.00004211,0.000041795,0.00004148,0.000041165,0.00004085,0.000040535,0.00004022,0.000039905,0.00003959,0.000039275,0.00003896,0.000038645,0.00003833,0.000038015,3.77E-005,0.00003717,0.00003664,0.00003611,0.00003558,0.00003505,0.00003452,0.00003399,0.00003346,0.00003293,0.0000324,0.00003187,0.00003134,0.00003081,0.00003028,0.00002975,0.00002922,0.00002869,0.00002816,0.00002763,2.71E-005,0.000026395,0.00002569,0.000024985,0.00002428,0.000023575,0.00002287,0.000022165,0.00002146,0.000020755,0.00002005,0.000019345,0.00001864,0.000017935,0.00001723,0.000016525,0.00001582,0.000015115,0.00001441,0.000013705,1.30E-005,0.000012421,0.000011842,0.000011263,0.000010684,0.000010105,0.000009526,0.000008947,0.000008368,0.000007789,0.00000721,0.000006631,0.000006052,0.000005473,0.000004894,0.000004315,0.000003736,0.000003157,0.000002578,0.000001999,1.42E-006};

  
  double gghxsec[numsmpoints] = {19.840,19.650,19.480,19.300,19.130,18.950,18.790,18.620,18.450,18.290,18.130,17.970,17.820,17.660,17.510,17.360,17.210,17.060,16.920,16.780,16.630,16.490,16.360,16.220,16.080,15.940,15.820,15.690,15.560,15.430,15.310,15.180,15.060,14.940,14.820,14.700,14.580,14.460,14.350,14.230,14.120,14.010,13.900,13.800,13.690,13.580,13.480,13.370,13.280,13.180,13.080,12.980,12.880,12.780,12.680,12.600,12.500,12.400,12.310,12.220,12.130,12.0,11.950,11.9,11.780,11.7,11.600,11.5,11.440,11.4,11.270,11.2,11.120,11.0,10.960,10.9,10.800,10.7,10.650,10.6,10.500};
  double vbfxsec[numsmpoints] = {1.3980,1.3910,1.3840,1.3780,1.3710,1.3640,1.3580,1.3510,1.3450,1.3390,1.3320,1.3260,1.3190,1.3130,1.3070,1.3000,1.2940,1.2880,1.2820,1.2760,1.2690,1.2630,1.2570,1.2510,1.2460,1.2400,1.2340,1.2280,1.2220,1.2160,1.2110,1.2050,1.1990,1.1930,1.1880,1.1820,1.1760,1.1710,1.1650,1.1590,1.1540,1.1480,1.1430,1.1370,1.1320,1.1260,1.1210,1.1150,1.1100,1.1050,1.1000,1.0950,1.0900,1.0850,1.0800,1.0760,1.0710,1.0660,1.0620,1.0570,1.0520,1.0475,1.0430,1.0380,1.0330,1.0280,1.0230,1.0180,1.0130,1.0085,1.0040,0.9996,0.9951,0.9909,0.9866,0.9824,0.9782,0.9741,0.9699,0.9658,0.9617};
  double whxsec[numsmpoints] = {0.8754,0.8623,0.8495,0.8368,0.8244,0.8122,0.8003,0.7885,0.7770,0.7657,0.7546,0.7439,0.7333,0.7230,0.7129,0.7030,0.6933,0.6837,0.6744,0.6651,0.6561,0.6472,0.6384,0.6297,0.6212,0.6129,0.6046,0.5965,0.5885,0.5806,0.5729,0.5652,0.5576,0.5501,0.5428,0.5355,0.5284,0.5213,0.5144,0.5075,0.5008,0.4942,0.4877,0.4813,0.4749,0.4687,0.4626,0.4566,0.4506,0.4448,0.4390,0.4333,0.4277,0.4221,0.4167,0.4113,0.4060,0.4008,0.3957,0.3907,0.3857,0.3809,0.3761,0.3715,0.3669,0.3624,0.3579,0.3535,0.3491,0.3449,0.3406,0.3364,0.3321,0.3280,0.3238,0.3198,0.3157,0.3118,0.3078,0.3040,0.3001};
  double zhxsec[numsmpoints] = {0.4721,0.4655,0.4589,0.4525,0.4462,0.4400,0.4340,0.4280,0.4221,0.4164,0.4107,0.4052,0.3998,0.3945,0.3893,0.3842,0.3791,0.3742,0.3693,0.3645,0.3598,0.3551,0.3505,0.3459,0.3414,0.3370,0.3326,0.3283,0.3241,0.3199,0.3158,0.3117,0.3077,0.3038,0.2999,0.2961,0.2923,0.2886,0.2849,0.2813,0.2778,0.2743,0.2709,0.2675,0.2642,0.2609,0.2577,0.2545,0.2514,0.2483,0.2453,0.2423,0.2393,0.2364,0.2336,0.2307,0.2279,0.2252,0.2225,0.2198,0.2172,0.2147,0.2121,0.2096,0.2071,0.2047,0.2023,0.2000,0.1976,0.1953,0.1930,0.1907,0.1884,0.1862,0.1840,0.1818,0.1796,0.1775,0.1754,0.1734,0.1713};
  double tthxsec[numsmpoints] = {0.12570,0.12410,0.12250,0.12090,0.11940,0.11790,0.11640,0.11490,0.11340,0.11200,0.11060,0.10920,0.10780,0.10650,0.10510,0.10380,0.10250,0.10130,0.10000,0.09878,0.09756,0.09636,0.09518,0.09402,0.09287,0.09174,0.09063,0.08954,0.08846,0.08739,0.08634,0.08530,0.08428,0.08327,0.08227,0.08129,0.08032,0.07937,0.07842,0.07750,0.07658,0.07568,0.07479,0.07391,0.07304,0.07219,0.07135,0.07052,0.06970,0.06890,0.06810,0.06731,0.06654,0.06577,0.06502,0.06428,0.06355,0.06282,0.06211,0.06141,0.06072,0.06005,0.05937,0.05872,0.05807,0.05744,0.05680,0.05618,0.05556,0.05496,0.05435,0.05376,0.05316,0.05258,0.05200,0.05144,0.05087,0.05032,0.04976,0.04923,0.04869};

  double sm4gghxsec[numsmpoints] = {199,197.3,195.6,193.9,192.2,190.5,188.8,187.1,185.4,183.7,182,180.3,178.6,176.9,175.2,173.5,171.8,170.1,168.4,166.7,165,163.65,162.3,160.95,159.6,158.25,156.9,155.55,154.2,152.85,151.5,150.15,148.8,147.45,146.1,144.75,143.4,142.05,140.7,139.35,138,136.95,135.9,134.85,133.8,132.75,131.7,130.65,129.6,128.55,127.5,126.45,125.4,124.35,123.3,122.25,121.2,120.15,119.1,118.05,117,116,115,114,113,112,111,110,109,108,107,106.22,105.44,104.66,103.88,103.1,102.32,101.54,100.76,99.98,99.2};
  
  TH1D *hsmbrs = new TH1D("hsmbr","",81,109.75,150.25);
  TH1D *hffbrs = new TH1D("hffbrs","",81,109.75,150.25);
  TH1D *hsm4brs = new TH1D("hsm4brs","",81,109.75,150.25);
 
  
  TH1D *hgghxsecs = new TH1D("hgghxsecs","",81,109.75,150.25);
  TH1D *hvbfxsecs = new TH1D("hvbfxsecs","",81,109.75,150.25);
  TH1D *hwzhxsecs = new TH1D("hwzhxsecs","",81,109.75,150.25);
  TH1D *htthxsecs = new TH1D("htthxsecs","",81,109.75,150.25);
  
  TH1D *hsm4gghxsecs = new TH1D("hsm4gghxsecs","",81,109.75,150.25);
  
  
  for (int ipoint=0; ipoint<numsmpoints; ++ipoint) {
    hsmbrs->Fill(smmasses[ipoint],smbr[ipoint]);
    hffbrs->Fill(smmasses[ipoint],ffbr[ipoint]);
    hsm4brs->Fill(smmasses[ipoint],sm4br[ipoint]);    
    
    hgghxsecs->Fill(smmasses[ipoint],gghxsec[ipoint]);
    hvbfxsecs->Fill(smmasses[ipoint],vbfxsec[ipoint]);
    hwzhxsecs->Fill(smmasses[ipoint],whxsec[ipoint]+zhxsec[ipoint]);
    htthxsecs->Fill(smmasses[ipoint],tthxsec[ipoint]);
    
    hsm4gghxsecs->Fill(smmasses[ipoint],sm4gghxsec[ipoint]);    
  }
  RooDataHist *dsmbrs = new RooDataHist("dsmbrs","",RooArgList(mnom),hsmbrs);
  RooHistFunc *fsmbrs = new RooHistFunc("fsmbrs","",RooArgList(mnom),*dsmbrs,1);

  RooDataHist *dffbrs = new RooDataHist("dffbrs","",RooArgList(mnom),hffbrs);
  RooHistFunc *fffbrs = new RooHistFunc("fffbrs","",RooArgList(mnom),*dffbrs,1);  
  
  RooDataHist *dsm4brs = new RooDataHist("dsm4brs","",RooArgList(mnom),hsm4brs);
  RooHistFunc *fsm4brs = new RooHistFunc("fsm4brs","",RooArgList(mnom),*dsm4brs,1);  
  
  RooDataHist *dgghxsecs = new RooDataHist("dgghxsecs","",RooArgList(mnom),hgghxsecs);
  RooHistFunc *fgghxsecs = new RooHistFunc("fgghxsecs","",RooArgList(mnom),*dgghxsecs,1);  
 
  RooDataHist *dvbfxsecs = new RooDataHist("dvbfxsecs","",RooArgList(mnom),hvbfxsecs);
  RooHistFunc *fvbfxsecs = new RooHistFunc("fvbfxsecs","",RooArgList(mnom),*dvbfxsecs,1);    
 
  RooDataHist *dwzhxsecs = new RooDataHist("dwzhxsecs","",RooArgList(mnom),hwzhxsecs);
  RooHistFunc *fwzhxsecs = new RooHistFunc("fwzhxsecs","",RooArgList(mnom),*dwzhxsecs,1);    
  
  RooDataHist *dtthxsecs = new RooDataHist("dtthxsecs","",RooArgList(mnom),htthxsecs);
  RooHistFunc *ftthxsecs = new RooHistFunc("ftthxsecs","",RooArgList(mnom),*dtthxsecs,1);    
  
  RooAddition *fsmxsecs = new RooAddition("fsmxsecs","",RooArgList(*fgghxsecs,*fvbfxsecs,*fwzhxsecs,*ftthxsecs));
  RooAddition *fffxsecs = new RooAddition("fffxsecs","",RooArgList(*fvbfxsecs,*fwzhxsecs));
  
  RooDataHist *dsm4gghxsecs = new RooDataHist("dsm4gghxsecs","",RooArgList(mnom),hsm4gghxsecs);
  RooHistFunc *fsm4gghxsecs = new RooHistFunc("fsm4gghxsecs","",RooArgList(mnom),*dsm4gghxsecs,1);  
  
  RooAddition *fsm4xsecs = new RooAddition("fsm4xsecs","",RooArgList(*fsm4gghxsecs,*fvbfxsecs,*fwzhxsecs,*ftthxsecs));

  
  std::vector<RooAbsReal*> procxseclist;
  procxseclist.push_back(fgghxsecs);
  procxseclist.push_back(fvbfxsecs);
  procxseclist.push_back(fwzhxsecs);
  procxseclist.push_back(ftthxsecs);
  
  RooAbsReal *procxsecs = procxseclist.at(filterproc);

  std::vector<RooAbsReal*> procxsecsm4list;
  procxsecsm4list.push_back(fsm4gghxsecs);
  procxsecsm4list.push_back(fvbfxsecs);
  procxsecsm4list.push_back(fwzhxsecs);
  procxsecsm4list.push_back(ftthxsecs);
  
  RooAbsReal *procxsecssm4 = procxsecsm4list.at(filterproc);  
  
//   TCanvas *csmxsecs = new TCanvas;
//   RooPlot *plotsmxsecs = mnom.frame();
//   dsmxsecs->plotOn(plotsmxsecs);
//   fsmxsecs->plotOn(plotsmxsecs);
//   plotsmxsecs->Draw();
//   csmxsecs->SaveAs("xsecs.eps");
//   
//   RooAbsReal *xsecnorm = fsmxsecs;
    
    //loop over categories for fits
  RooWorkspace *wextra = new RooWorkspace("wextra","") ;
    
  TString procname = procnames.at(filterproc);
  
  for (UInt_t i=0; i<mhs.size(); ++i) {
    std::stringstream numstringstr;
    numstringstr<<mhs.at(i);
    TString numstring(numstringstr.str());

    
    if (mhs.at(i)<=120) win = win1;
    else win = win2;
    
    RooDataSet *rvdata = new RooDataSet(TString::Format("sig_mass_rv_m%i_combcat_%s",mhs.at(i),procname.Data()),"",RooArgList(*hmass,*weight),weight->GetName());
    RooDataSet *wvdata = new RooDataSet(TString::Format("sig_mass_wv_m%i_combcat_%s",mhs.at(i),procname.Data()),"",RooArgList(*hmass,*weight),weight->GetName());
    RooDataSet *alldata = new RooDataSet(TString::Format("sig_mass_m%i_combcat_%s",mhs.at(i),procname.Data()),"",RooArgList(*hmass,*weight),weight->GetName());
    for (UInt_t icat=0; icat<catnames.size(); ++icat) {
      RooDataSet *mcsigdata = (RooDataSet*)win->data(TString::Format("sig_%s_mass_m",procname.Data())+numstring+TString("_rv_")+catnamesbase.at(icat));
      RooDataSet *mcsigwrongdata = (RooDataSet*)win->data(TString::Format("sig_%s_mass_m",procname.Data())+numstring+TString("_wv_")+catnamesbase.at(icat));
      RooDataSet *mcsigalldata = (RooDataSet*)win->data(TString::Format("sig_%s_mass_m",procname.Data())+numstring+TString("_")+catnamesbase.at(icat));

      printf("mcsigdata = %p, mcsigwrongdata = %p, mcsigalldata = %p\n",(void*)mcsigdata,(void*)mcsigwrongdata,(void*)mcsigalldata);
      
      //rvdata->append(*mcsigdata);
      //wvdata->append(*mcsigwrongdata);
      //alldata->append(*mcsigalldata);
      appendcwd(rvdata,mcsigdata,hmass,weight,1.0);
      appendcwd(wvdata,mcsigwrongdata,hmass,weight,1.0);
      appendcwd(alldata,mcsigalldata,hmass,weight,1.0);


    }
    
    win->import(*rvdata);
    win->import(*wvdata);
    win->import(*alldata);
    
    wextra->import(*rvdata);
    wextra->import(*wvdata);
    wextra->import(*alldata);    

    
  }
  //catnames.push_back("singlecat");  
  
  
  const uint ncats = catnames.size();  
  const uint nparms = fitparms.size();
    
  //define histograms to keep track of fit parameters for each mass point
  TH1F **fitparmhists = new TH1F*[ncats*nparms];
  for (UInt_t icat=0; icat<catnames.size(); ++icat) {
    for (UInt_t iparm=0; iparm<fitparms.size(); ++iparm) {
      TString histname = TString("hist") + TString(fitparms.at(iparm)->GetName()) + catnames.at(icat);
      fitparmhists[icat*nparms + iparm] = new TH1F(histname,histname,(massmax-110.0)/10.0,105,massmax-5.0);
    }
  }
  

  TH2F *histfitstatus = new TH2F("histfitstatus","histfitstatus",(massmax-110.0)/10.0,105,massmax-5.0,catnames.size(),-0.5,catnames.size()-0.5);
  TH2F *histfitstatuswrong = new TH2F("histfitstatuswrong","histfitstatuswrong",(massmax-110.0)/10.0,105,massmax-5.0,catnames.size(),-0.5,catnames.size()-0.5);
  
  std::vector<RooDataSet*> testdsets;
  std::vector<RooDataSet*> testdsetswrong;
  std::vector<RooDataSet*> testdsetsall;
  
  

  
  
  
  //loop over categories for fits
  for (UInt_t icat=0; icat<catnames.size(); ++icat) {
  
    
    //printf("intlumi = %5f\n",IntLumi->getVal());
    //return;
    //TCut rightcut = catcuts.at(icat) && rightvtx;
    //TCut wrongcut = catcuts.at(icat) && wrongvtx;
    //TCut allcut = catcuts.at(icat);
   
    //loop over mass points for fits
    for (UInt_t i=0; i<mhs.size(); ++i) {

      if (mhs[i]<=120) win = win1;
      else win = win2;
      
      std::stringstream numstringstr;
      numstringstr<<mhs.at(i);
      TString numstring(numstringstr.str());
      
      TString vbfstring = numstring;
      
      //remap names of mislabeled MC samples
//       if (mhs.at(i)==90) vbfstring = "130";
//       else if (mhs.at(i)==100) vbfstring = "90";
//       else if (mhs.at(i)==110) vbfstring = "100";
//       else if (mhs.at(i)==115) vbfstring = "105";    
//       else if (mhs.at(i)==120) vbfstring = "110";
//       else if (mhs.at(i)==130) vbfstring = "115";
//       else if (mhs.at(i)==140) vbfstring = "120";
//       
//       TString vtthstring = numstring;
//       if (mhs.at(i)==90) vtthstring = "95";
//       if (mhs.at(i)==95) vtthstring = "90";
      
      //load correct pt weights
      //ptweights= (TH1D*) fileptweight->Get(TString("kfactors_") + numstring);

    
      //TFile *fgf = new TFile(TString(" /home/bendavid/cms/hist/hgg-v0-May19/local/filefi/merged/hgg-v0_p11-h") + numstring + TString("gg-gf-v1g1-pu_noskim.root"),"READ");
      //TFile *fvbf = new TFile(TString(" /home/bendavid/cms/hist/hgg-v0-May19/local/filefi/merged/hgg-v0_p11-h") + vbfstring + TString("gg-vbf-v1g1-pu_noskim.root"),"READ");
      //TFile *fvtth = new TFile(TString(" /home/bendavid/cms/hist/hgg-v0-May19/local/filefi/merged/hgg-v0_p11-h") + vtthstring + TString("gg-vtth-v1g1-pu_noskim.root"),"READ");

      
      //define aggregate weight, so far using xsec, pileup, pt-reweighting and efficiency scale factors
      //RooFormulaVar totweight("totweight","totweight","xsecweight*puweight(ngenvtx-1)*ptweight(genhpt,isgf)*effweight(iseb1,iseb2)",RooArgList(xsecweight,ngenvtx,genhpt,isgf,iseb1,iseb2));     
      //const double weightscale = 3e3;
      
      //right vertex
      RooDataSet *mcsigdata = (RooDataSet*)win->data(TString::Format("sig_%s_mass_m",procname.Data())+numstring+TString("_rv_")+catnamesbase.at(icat));
      const double weightscale = static_cast<double>(mcsigdata->numEntries())/mcsigdata->sumEntries();
      RooDataSet *mcsigwdata = cwdset(mcsigdata,hmass,weight,procidx,TString("mcsigwdata") + numstring+catnames.at(icat),weightscale,filterproc);
       
      
      //wrong vertex
      RooDataSet *mcsigwrongdata = (RooDataSet*)win->data(TString::Format("sig_%s_mass_m",procname.Data())+numstring+TString("_wv_")+catnamesbase.at(icat));
      RooDataSet *mcsigwrongwdata = cwdset(mcsigwrongdata,hmass,weight,procidx,TString("mcsigwrongwdata") + numstring+catnames.at(icat),weightscale,filterproc);
            
      //combined
      RooDataSet *mcsigalldata = (RooDataSet*)win->data(TString::Format("sig_%s_mass_m",procname.Data())+numstring+TString("_")+catnamesbase.at(icat));
      RooDataSet *mcsigallwdata = cwdset(mcsigalldata,hmass,weight,procidx,TString("mcsigallwdata") + numstring+catnames.at(icat),weightscale,filterproc);

            
      printf("sigdata = %i, sigwrong = %i, sigall = %i\n",mcsigdata!=0, mcsigwrongdata!=0,mcsigalldata!=0);
      
//       printf("right = %5f, wrong = %5f, all = %5f\n",mcsigdata->sumEntries(),mcsigwrongdata->sumEntries(), mcsigalldata->sumEntries());
//       printf("right = %5f, wrong = %5f, all = %5f\n",mcsigwdata->sumEntries(),mcsigwrongwdata->sumEntries(), mcsigallwdata->sumEntries());
//       return;
//       
      //track some test datasets for later tests
     if (mhs.at(i)==115) {        
        testdsets.push_back(mcsigwdata);
        testdsetswrong.push_back(mcsigwrongwdata);
        testdsetsall.push_back(mcsigallwdata);
        
        //continue;
      }
      
      //if (mhs.at(i)%10 != 0) continue;
      
      //set higgs mass and fit range
      mnom.setVal(mhs.at(i));
      hmass->setRange("higgsrange",100.0,massmax);
      hmass->setRange("plotrange",TMath::Max(100.0,mhs.at(i)-30.0),TMath::Min(massmax,mhs.at(i)+20.0));
      
      //reset parameters to initial values
      for (UInt_t iparm=0; iparm<fitparms.size(); ++iparm) {
        fitparms.at(iparm)->setVal(fitparmsinit.at(iparm));       
      }
      RooAbsPdf *rightpdf = 0;
      RooAbsPdf *wrongpdf = 0;
      RooAbsPdf *allpdf = 0;
      //usecomplexmodelw = catnames.at(icat)=="cat0" || catnames.at(icat)=="cat1"  || catnames.at(icat)=="cat4" || catnames.at(icat)=="cat5" || catnames.at(icat)=="singlecat"; 
      //Bool_t usecomplexmodel = catnames.at(icat)=="cat0" || catnames.at(icat)=="cat1"  || catnames.at(icat)=="cat4" || catnames.at(icat)=="cat5" || catnames.at(icat)=="singlecat";
      //Bool_t usecomplexmodel = catnames.at(icat)=="cat0" || catnames.at(icat)=="cat4" || catnames.at(icat)=="singlecat"; 
      //Bool_t usecomplexmodel = catnames.at(icat)=="cat0" || catnames.at(icat)=="cat1"  || catnames.at(icat)=="cat4" || catnames.at(icat)=="cat5";
      //Bool_t usesimplemodel = catnames.at(icat)=="cat2" || catnames.at(icat)=="cat3"  || catnames.at(icat)=="cat6" || catnames.at(icat)=="cat7";
      //Bool_t usesimplemodel = false;
      Bool_t usecomplexmodel = complexset.count(catnames.at(icat));
      Bool_t usesimplemodel = simpleset.count(catnames.at(icat));
      if (usecomplexmodel) {
        rightpdf = &combh;
        wrongpdf = &combhwrong;
        allpdf = &combhvtx;
        f1.setRange(0.52,1.0);
      }
      else if (usesimplemodel) {
        rightpdf = &g1;
        wrongpdf = &wg1;
        allpdf = &g1;
        sigma1.setVal(5.0);
      }
      else {
        rightpdf = &combhmin;
        wrongpdf = &wg1;
        allpdf = &combhvtxmin;
        sigma1.setVal(2.0);
        sigma2.setVal(3.0);
        f1.setRange(0.0,1.0);
      }
      

      
      if (mhs.at(i)==110 && icat>0) {
        f1.setVal(0.8);
        f1.setRange(0.6,1.0);
//         sigma1.setVal(3.0);
//         dm2.setVal(-3.0);
        dm2.setRange(-6.0,0.0);
      }      
      
      if (mhs.at(i)==110 && !usecomplexmodel) {
        f1.setVal(0.95);
        f1.setRange(0.9,1.0);
        sigma1.setVal(3.0);
        dm2.setVal(-3.0);
        dm2.setRange(-8.0,0.0);
      }
      
      if (mhs.at(i)==110 && icat==1) {
        sigma3.setRange(3.0,12.0);
      }
      
      if (icat==0) {
        sigma1.removeRange();
        sigma1.setVal(0.7);
        sigma1.removeRange();
        
//         dm1.removeRange();
//         dm1.setVal(0.05); 
      }
      
/*      if (icat==0 && mhs[i]==150) {
        sigma1.removeRange();
        sigma1.setVal(0.9);
        sigma1.removeRange();
        
//         dm1.removeRange();
//         dm1.setVal(0.1); 
      }   */   
      
/*      if (mhs.at(i)==110 && icat==0) {
        dm1.setRange(-1.0,1.0);
        dm1.setVal(0.3);
      }   */   
      
//     if (mhs.at(i)==110 &&) {
//       dm2.setRange(-6.0,0.0);
//     }
      
//       if (mhs.at(i)==150) {
//         dm1.setVal(0.0);
//         dm1.setRange(-2.0,1.0);        
//       }
//       if (mhs.at(i)==110) {
//         dm1.setRange(-3.0,2.0);
//         dm2.setRange(-8.0,5.0);
//       }      
//       if (mhs.at(i)==110) {
//         dm1.setVal(0.0);
//       }      
//       if (mhs.at(i)==150 || mhs.at(i)==130) {
//         dm1.setVal(0.0);
//         dm1.setVal(-0.5);
//         dm2.setVal(-1.5);
//       }
//       if (mhs.at(i)==150) {
//         dm2.setRange(-8.0,5.0);
//       }

      if (icat==0) {
        sigma1.setRange(0.55,1.2);
      }

      if (icat==0 && !usesimplemodel) {
        sigma1.setVal(1.0);
        sigma2.setVal(1.5);
        sigma3.setVal(2.5);    
        if (filterproc==0) sigma3.setVal(4.5);
      }

      if (icat>0 && !usesimplemodel) {
        sigma1.setVal(1.0);
        sigma2.setVal(1.5);
        sigma3.setVal(2.5);
      }
      
      //if (icat==3 && mhs.at(i)==130 && usecomplexmodel) {
//       if (icat>=0) {
//         sigma1.setVal(3.0);
//         sigma2.setVal(3.5);
//         sigma3.setVal(4.0);
//       }
      
      
      
//       if (icat==0 && mhs.at(i)==110) {
//         dm1.setVal(0.0);
//         dm1.setRange(-1.0,0.0);
//         dm2.setRange(-5.0,0.0);
//         dm3.setRange(-8.0,0.0);        
//       }
      
//       if ((icat==1||icat==3) && mhs.at(i)==110) {
//         f1.setVal(0.90);
//         f1.setRange(0.65,1.0);
//         sigma1.setVal(3.0);
//         dm2.setVal(-3.0);
//         dm2.setRange(-8.0,0.0); 
//         dm1.setVal(0.0);
//         dm1.setRange(-1.0,1.0);
//         dm2.setRange(-5.0,0.0);
//         dm3.setRange(-8.0,0.0);         
//       }
//       else {
//         dm1.removeRange();
//       }
      
      
// 
//       
//       if (icat==1 && mhs.at(i)==140) {
//         sigma1.setVal(1.2);
//         sigma2.setVal(2.0);
//         sigma3.setVal(3.0);
//       }   
//       
//       if (icat==2 && mhs.at(i)==130) {
//         sigma1.setVal(1.2);
//         sigma2.setVal(2.0);
//         sigma3.setVal(3.0);
//       }     
//       
//       
//       if (icat==0) {
//         sigma1.setRange(0.5,3.0);
//         sigma2.setRange(1.0,5.0);
//         dm1.setVal(0.0);
//       }
      if (icat>=2) {
        sigma1.setRange(1.0,4.0);
        sigma2.setRange(1.5,7.0);        
      }
      
      
      
      RooFitResult *fitres = 0;
      RooFitResult *fitreswrong = 0;
      
      if (!usesimplemodel) {
        fitres = rightpdf->fitTo(*mcsigwdata,Strategy(1),Minimizer("Minuit2",""),Minos(kFALSE),SumW2Error(kFALSE), Save(kTRUE),NumCPU(8)); 
        if (!fitres->status()) fitres = rightpdf->fitTo(*mcsigwdata,Strategy(0),Minimizer("Minuit2",""),Minos(kFALSE),SumW2Error(kFALSE), Save(kTRUE),NumCPU(8)); 
        if (!fitres->status()) fitres = rightpdf->fitTo(*mcsigwdata,Strategy(2),Minimizer("Minuit2",""),Minos(kFALSE),SumW2Error(kFALSE), Save(kTRUE),NumCPU(8)); 
        fitreswrong = wrongpdf->fitTo(*mcsigwrongwdata,Strategy(1),Minimizer("Minuit2",""),Minos(kFALSE),SumW2Error(kFALSE), Save(kTRUE),NumCPU(8)); 
        if (!fitreswrong->status()) fitreswrong = wrongpdf->fitTo(*mcsigwrongwdata,Strategy(0),Minimizer("Minuit2",""),Minos(kFALSE),SumW2Error(kFALSE), Save(kTRUE),NumCPU(8)); 
      }
      
      else {
        fitres = rightpdf->fitTo(*mcsigallwdata,Strategy(0),Minimizer("Minuit2",""),Minos(kFALSE),SumW2Error(kFALSE), Save(kTRUE),NumCPU(8)); 
        if (!fitres->status()) fitres = rightpdf->fitTo(*mcsigallwdata,Strategy(1),Minimizer("Minuit2",""),Minos(kFALSE),SumW2Error(kFALSE), Save(kTRUE),NumCPU(8)); 
        if (!fitres->status()) fitres = rightpdf->fitTo(*mcsigallwdata,Strategy(2),Minimizer("Minuit2",""),Minos(kFALSE),SumW2Error(kFALSE), Save(kTRUE),NumCPU(8));         
        fitreswrong = fitres;
      }
      
      //dm1.removeRange();
      dm1.setRange(-0.5,1.5);
      dm2.removeRange();
      
      if (usecomplexmodel) {
        if (TMath::Abs(sigma2.getVal())>TMath::Abs(sigma3.getVal()) || f2.getVal()<0.5) {
          double sigma2tmp = sigma2.getVal();
          double dm2tmp = dm2.getVal();
          
          sigma2.setVal(sigma3.getVal());
          dm2.setVal(dm3.getVal());
          sigma3.setVal(sigma2tmp);
          dm3.setVal(dm2tmp);
          f2.setVal(1.0-f2.getVal());
          
        }
        
        
      }
      if (!usesimplemodel) {
        //printf("testing sigmas, sigma1 = %5f, sigma2 = %5f\n",sigma1.getVal(),sigma2.getVal());
        //if (TMath::Abs(sigma1.getVal())>TMath::Abs(sigma2.getVal()) || f1.getVal()<0.5 || dm1.getVal()<dm2.getVal()) {
        if (f1.getVal()<0.5 || (usecomplexmodel && (TMath::Abs(sigma1.getVal())>TMath::Abs(sigma2.getVal()) || f1.getVal()<0.5 || dm1.getVal()<dm2.getVal()))   ) {
          //printf("swapping gaussians\n");
          double sigma1tmp = sigma1.getVal();
          double dm1tmp = dm1.getVal();
          
          sigma1.setVal(sigma2.getVal());
          dm1.setVal(dm2.getVal());
          sigma2.setVal(sigma1tmp);
          dm2.setVal(dm1tmp);
          f1.setVal(1.0-f1.getVal());
        }
        
      }
      
      RooRealVar *mtotalxsec = win->var(TString::Format("XSBR_%s_%i",procname.Data(),mhs.at(i)));
      //RooRealVar *mtotalxsec = win->var(TString("ff_XSBR_")+numstring);
      //RooRealVar *mtotalxsec = win->var(TString("ggH_XSBR_")+numstring);
      
      printf("lumi = %5f, xsec = %5f\n",IntLumi->getVal(),mtotalxsec->getVal());
      //return;
      //RooConstVar *mtotalxsec = new RooConstVar(TString("XSBR_")+numstring+catnames.at(icat),"",ffxsecs.at(i));
      
/*      if(icat==1) { 
        printf("xsec = %5f\n",mtotalxsec->getVal());
        return;
      }*/
      
      //compute acceptance*efficiency and right vertex fraction
      double eaccnum = mcsigwdata->sumEntries()+mcsigwrongwdata->sumEntries();
      double eaccden = IntLumi->getVal()*mtotalxsec->getVal()*weightscale;
      //double eaccden = (gfxsecs.at(i)+vbfxsecs.at(i)+vhxsecs.at(i))*1e6;
      //double eaccden = (vbfxsecs.at(i)+vhxsecs.at(i))*1e6;

      double eacc = eaccnum/eaccden;
      printf("eacc = %5f, eaccnum = %5f, eaccden = %5f\n",eacc,eaccnum,eaccden);
      double eaccerrlo = TEfficiency::ClopperPearson(Int_t(eaccden), Int_t(eaccnum), 0.683, kFALSE) - eacc;
      double eaccerrhi = TEfficiency::ClopperPearson(Int_t(eaccden), Int_t(eaccnum), 0.683, kTRUE) - eacc;
      printf("eff done\n");
      
      
      //double fright = mcsigwdata->sumEntries()/(mcsigwdata->sumEntries()+mcsigwrongwdata->sumEntries());
      double frightnum = mcsigwdata->sumEntries();
      double frightden = mcsigwdata->sumEntries()+mcsigwrongwdata->sumEntries();
      double fright = frightnum/frightden;
      double frighterrlo = TEfficiency::ClopperPearson(Int_t(frightden), Int_t(frightnum), 0.683, kFALSE) - fright;
      double frighterrhi = TEfficiency::ClopperPearson(Int_t(frightden), Int_t(frightnum), 0.683, kTRUE) - fright;
      
      
      fracright.setVal(fright);
      fracright.setAsymError(frighterrlo,frighterrhi);
      
      effacc.setVal(eacc);
      effacc.setAsymError(eaccerrlo,eaccerrhi);
      
      //correct negative resoltion terms which screw up interpolation
      fitparms.at(3)->setVal(TMath::Abs(fitparms.at(3)->getVal()));
      fitparms.at(4)->setVal(TMath::Abs(fitparms.at(4)->getVal()));
      fitparms.at(5)->setVal(TMath::Abs(fitparms.at(5)->getVal()));
      fitparms.at(10)->setVal(TMath::Abs(fitparms.at(10)->getVal()));
      fitparms.at(11)->setVal(TMath::Abs(fitparms.at(11)->getVal()));
      
      //fill histograms with fit parameters in 5 GeV steps
      if ( (mhs.at(i)%10) == 0) {
      //if (1) {
      
        for (UInt_t iparm=0; iparm<fitparms.size(); ++iparm) {
          fitparmhists[icat*nparms + iparm]->Fill(mhs.at(i),fitparms.at(iparm)->getVal()); 
          fitparmhists[icat*nparms + iparm]->SetBinError(fitparmhists[icat*nparms + iparm]->FindFixBin(mhs.at(i)),fitparms.at(iparm)->getError()); 
          if (fitparms.at(iparm)->hasAsymError()) {
            double avgerror = (fitparms.at(iparm)->getErrorLo() + fitparms.at(iparm)->getErrorHi())/2.0;
            fitparmhists[icat*nparms + iparm]->SetBinError(fitparmhists[icat*nparms + iparm]->FindFixBin(mhs.at(i)),avgerror); 
          }
          printf("filling histogram named: %s, variable named %s, val = %5f, err = %5f\n",fitparmhists[icat*nparms + iparm]->GetName(),fitparms.at(iparm)->GetName(),fitparms.at(iparm)->getVal(),fitparms.at(iparm)->getError());
        }
                
        histfitstatus->Fill(mhs.at(i),icat,fitres->status());
        histfitstatuswrong->Fill(mhs.at(i),icat,fitreswrong->status());
        
      }
      
      //plot fit results for this mass point
      TCanvas *chfit = new TCanvas;
      TString plotname = TString("rightvtx") + numstring + catnames.at(icat) + TString(".eps");      
      RooPlot *hplot = hmass->frame(Bins(100),Range("plotrange"));
      mcsigwdata->plotOn(hplot);
      rightpdf->plotOn(hplot,Components(g1),LineColor(kOrange),Range("higgsrange"),NormRange("higgsrange"));
      rightpdf->plotOn(hplot,Components(g2),LineColor(kMagenta),Range("higgsrange"),NormRange("higgsrange"));
      rightpdf->plotOn(hplot,Components(g3),LineColor(kRed),Range("higgsrange"),NormRange("higgsrange"));
 
      rightpdf->plotOn(hplot,RooFit::LineColor(kBlue),Range("higgsrange"),NormRange("higgsrange"));  
      hplot->SetTitle("");
      hplot->Draw();  
      chfit->SaveAs(plotname);

      TCanvas *chfitwrong = new TCanvas;
      TString plotnamewrong = TString("wrongvtx") + numstring + catnames.at(icat) + TString(".eps");            
      RooPlot *hplotwrong = hmass->frame(Bins(40),Range("plotrange"));
      mcsigwrongwdata->plotOn(hplotwrong);
      wrongpdf->plotOn(hplotwrong,Components(wg1),LineColor(kOrange),Range("higgsrange"),NormRange("higgsrange"));
      wrongpdf->plotOn(hplotwrong,Components(wg2),LineColor(kMagenta),Range("higgsrange"),NormRange("higgsrange"));      
      wrongpdf->plotOn(hplotwrong,RooFit::LineColor(kBlue),Range("higgsrange"),NormRange("higgsrange"));  
      hplotwrong->SetTitle("");
      hplotwrong->Draw();        
      chfitwrong->SaveAs(plotnamewrong);
      
      TCanvas *chfitall = new TCanvas;
      TString plotnameall = TString("allvtx") + numstring + catnames.at(icat) + TString(".eps");                  
      RooPlot *hplotall = hmass->frame(Bins(100),Range("plotrange"));
      mcsigallwdata->plotOn(hplotall);
      allpdf->plotOn(hplotall,RooFit::LineColor(kBlue),Range("higgsrange"),NormRange("higgsrange"));  
      hplotall->SetTitle("");
      hplotall->Draw();   
      chfitall->SaveAs(plotnameall);
      
      printf ("right = %5f, wrong = %5f, all = %5f, right+wrong = %5f\n",mcsigwdata->sumEntries(),mcsigwrongwdata->sumEntries(), mcsigallwdata->sumEntries(),mcsigwdata->sumEntries()+mcsigwrongwdata->sumEntries());    
      printf("data weights = %5e\n",mcsigwdata->sumEntries());
      printf("mass = %i, status = %i, statuswrong = %i\n",mhs.at(i),fitres->status(),fitreswrong->status());
      
      //return;
      
    }
  }
  
  //construct interpolation construction for each fit parameter in each category
  RooDataHist **fitparmdatas = new RooDataHist*[ncats*nparms];
  RooHistFunc **fitparmfuncs = new RooHistFunc*[ncats*nparms];
  for (UInt_t icat=0; icat<catnames.size(); ++icat) {
    for (UInt_t iparm=0; iparm<fitparms.size(); ++iparm) {
      TString dataname = TString("data") + TString(fitparms.at(iparm)->GetName()) + catnames.at(icat);
      TString funcname = TString("func") + TString(fitparms.at(iparm)->GetName()) + catnames.at(icat);
      
      fitparmdatas[icat*nparms + iparm] = new RooDataHist(dataname,dataname,RooArgList(mnom),fitparmhists[icat*nparms + iparm]);
      fitparmfuncs[icat*nparms + iparm] = new RooHistFunc(funcname,funcname,RooArgList(mnom),*fitparmdatas[icat*nparms + iparm],1);      
    }
  }

  //plot evolution of each fit parameter in each category as function of higgs mass
  for (UInt_t icat=0; icat<catnames.size(); ++icat) {
    for (UInt_t iparm=0; iparm<fitparms.size(); ++iparm) {
      TString plotname = TString("func") + TString(fitparms.at(iparm)->GetName()) + catnames.at(icat) + TString(".eps");
      
      TCanvas *cfunctest = new TCanvas;
    //  fitparmhists[icat*nparms + iparm]->Draw();
      RooPlot *hploteffacc = mnom.frame(Bins(100),Range(105,155));
      fitparmdatas[icat*nparms + iparm]->plotOn(hploteffacc);  
      fitparmfuncs[icat*nparms + iparm]->plotOn(hploteffacc,RooFit::LineColor(kBlue));  
      hploteffacc->SetTitle("");
      hploteffacc->GetYaxis()->SetTitle(fitparms.at(iparm)->GetName());
      hploteffacc->Draw(); 
      cfunctest->SaveAs(plotname);
      delete cfunctest;
    } 
  }
  
  TCanvas *cfitstatus = new TCanvas;
  histfitstatus->Draw("COL");
  cfitstatus->SaveAs("fitstatus.eps");

  TCanvas *cfitstatuswrong = new TCanvas;
  histfitstatuswrong->Draw("COL");
  cfitstatuswrong->SaveAs("fitstatuswrong.eps");
  
  
 
 
  

  RooRealVar nuissancedeltafracright("CMS_hgg_nuissancedeltafracright","",1.0,0.1,10.0);
  nuissancedeltafracright.setConstant();
  
//   RooRealVar *nuissancedeltaescalebarrel = new RooRealVar("nuissancedeltaescalebarrel","",0.0,-0.1,0.1);
//   nuissanceescalebarrel->setConstant();
//   RooRealVar *nuissancedeltaescaleendcap = new RooRealVar("nuissancedeltaescaleendcap","",0.0,-0.1,0.1);
//   nuissanceescaleendcap->setConstant();
//   
//   RooFormulaVar *deltambarrel = new RooRealVar("deltambarrel","","@0*@1",RooArgList(mnom,*nuissancedeltaescalebarrel));
//   RooFormulaVar *deltammixedhighr9 = new RooFormulaVar("deltammixedhighr9","","@0*(0.7*(sqrt((@0-1.0)*(@1-1.0))-1.0) + 0.3*@2)",RooArgList(mnom,*nuissancedeltaescalebarrel,*nuissancedeltaescaleendcap));
//   RooFormulaVar *deltammixedmixedr9 = new RooFormulaVar("deltammixedmixedr9","","@0*(0.8*(sqrt((@0-1.0)*(@1-1.0))-1.0) + 0.2*@2)",RooArgList(mnom,*nuissancedeltaescalebarrel,*nuissancedeltaescaleendcap));
//   RooRealVar **nuissancedeltamcors = new RooRealVar*[8];
//   nuissancedeltamcors[0] = deltambarrel;
//   nuissancedeltamcors[1] = deltambarrel;
//   nuissancedeltamcors[2] = deltammixedhighr9;
//   nuissancedeltamcors[3] = deltammixedmixedr9;
//   nuissancedeltamcors[4] = deltambarrel;
//   nuissancedeltamcors[5] = deltambarrel;
//   nuissancedeltamcors[6] = deltammixedhighr9;
//   nuissancedeltamcors[7] = deltammixedmixedr9;
  
  
  //define final pdfs in each category
  RooConstVar **smears = new RooConstVar*[catnames.size()];
  RooRealVar **nuissancedeltasmears = new RooRealVar*[catnames.size()];
  RooRealVar **nuissancedeltams = new RooRealVar*[catnames.size()];
  RooFormulaVar **smearmods = new RooFormulaVar*[catnames.size()];
  RooFormulaVar **mean1slides = new RooFormulaVar*[catnames.size()];
  RooFormulaVar **mean2slides = new RooFormulaVar*[catnames.size()];
  RooFormulaVar **mean3slides = new RooFormulaVar*[catnames.size()];
  RooFormulaVar **sigma1slides = new RooFormulaVar*[catnames.size()];
  RooFormulaVar **sigma2slides = new RooFormulaVar*[catnames.size()];
  RooFormulaVar **sigma3slides = new RooFormulaVar*[catnames.size()];
  RooFormulaVar **wmean1slides = new RooFormulaVar*[catnames.size()];
  RooFormulaVar **wmean2slides = new RooFormulaVar*[catnames.size()];
  RooFormulaVar **wsigma1slides = new RooFormulaVar*[catnames.size()];
  RooFormulaVar **wsigma2slides = new RooFormulaVar*[catnames.size()];
  RooGaussian **g1slides = new RooGaussian*[catnames.size()];
  RooGaussian **g2slides = new RooGaussian*[catnames.size()];
  RooGaussian **g3slides = new RooGaussian*[catnames.size()];
  RooAddPdf **combhslides = new RooAddPdf*[catnames.size()];
  RooAddPdf **combhminslides = new RooAddPdf*[catnames.size()];
  RooGaussian **wg1slides = new RooGaussian*[catnames.size()];
  RooGaussian **wg2slides = new RooGaussian*[catnames.size()];
  RooAddPdf **combhwrongslides = new RooAddPdf*[catnames.size()];
  RooFormulaVar **fracrightmodslides = new RooFormulaVar*[catnames.size()];  
  RooAddPdf **combhvtxslides = new RooAddPdf*[catnames.size()];
  RooAddPdf **combhvtxminslides = new RooAddPdf*[catnames.size()];
  RooAbsPdf **combhvtxsimpleslides = new RooAbsPdf*[catnames.size()];
  RooAbsPdf **finalpdfslides = new RooAbsPdf*[catnames.size()];
  RooAbsReal **finalnormslides = new RooAbsReal*[catnames.size()]; 
  for (UInt_t icat=0; icat<catnames.size(); ++icat) {
    smears[icat] = new RooConstVar(TString("smear")+catnames.at(icat),"",smearingv.at(icat));

    nuissancedeltasmears[icat] = new RooRealVar(TString("CMS_hgg_nuissancedeltasmear")+catnamesbase.at(icat),"",0.0, -smearingv.at(icat),smearingv.at(icat));
    nuissancedeltasmears[icat]->setConstant();
    nuissancedeltams[icat] = new RooRealVar(TString("CMS_hgg_nuissancedeltam")+catnamesbase.at(icat),"",0.0,-5.0,5.0);
    nuissancedeltams[icat]->setConstant();
      
    smearmods[icat] = new RooFormulaVar(TString("smearmod")+catnames.at(icat),"","@0*(@1 + @2)",RooArgList(mnom,*smears[icat],*nuissancedeltasmears[icat]));
    
    mean1slides[icat] = new RooFormulaVar(TString("mean1slide")+catnames.at(icat),"","@0 + @1 + @0*@2",RooArgList(mnom,*fitparmfuncs[icat*nparms+0],*nuissancedeltams[icat]));
    mean2slides[icat] = new RooFormulaVar(TString("mean2slide")+catnames.at(icat),"","@0 + @1 + @0*@2",RooArgList(mnom,*fitparmfuncs[icat*nparms+1],*nuissancedeltams[icat]));
    mean3slides[icat] = new RooFormulaVar(TString("mean3slide")+catnames.at(icat),"","@0 + @1 + @0*@2",RooArgList(mnom,*fitparmfuncs[icat*nparms+2],*nuissancedeltams[icat]));
    sigma1slides[icat] = new RooFormulaVar(TString("sigma1slide")+catnames.at(icat),"","TMath::Max(0.01,sqrt(@0*@0-@3*@3*@2*@2 +@1*@1))",RooArgList(*fitparmfuncs[icat*nparms+3],*smearmods[icat],*smears[icat],mnom));
    sigma2slides[icat] = new RooFormulaVar(TString("sigma2slide")+catnames.at(icat),"","TMath::Max(0.01,sqrt(@0*@0-@3*@3*@2*@2 +@1*@1))",RooArgList(*fitparmfuncs[icat*nparms+4],*smearmods[icat],*smears[icat],mnom));
    sigma3slides[icat] = new RooFormulaVar(TString("sigma3slide")+catnames.at(icat),"","TMath::Max(0.01,sqrt(@0*@0-@3*@3*@2*@2 +@1*@1))",RooArgList(*fitparmfuncs[icat*nparms+5],*smearmods[icat],*smears[icat],mnom));
    wmean1slides[icat] = new RooFormulaVar(TString("wmean1slide")+catnames.at(icat),"","@0 + @1 + @0*@2",RooArgList(mnom,*fitparmfuncs[icat*nparms+8],*nuissancedeltams[icat]));
    wmean2slides[icat] = new RooFormulaVar(TString("wmean2slide")+catnames.at(icat),"","@0 + @1 + @0*@2",RooArgList(mnom,*fitparmfuncs[icat*nparms+9],*nuissancedeltams[icat]));
    wsigma1slides[icat] = new RooFormulaVar(TString("wsigma1slide")+catnames.at(icat),"","TMath::Max(0.01,sqrt(@0*@0-@3*@3*@2*@2 +@1*@1))",RooArgList(*fitparmfuncs[icat*nparms+10],*smearmods[icat],*smears[icat],mnom));
    wsigma2slides[icat] = new RooFormulaVar(TString("wsigma2slide")+catnames.at(icat),"","TMath::Max(0.01,sqrt(@0*@0-@3*@3*@2*@2 +@1*@1))",RooArgList(*fitparmfuncs[icat*nparms+11],*smearmods[icat],*smears[icat],mnom));
    
    g1slides[icat] = new RooGaussian(TString("g1slide")+catnames.at(icat),"",*hmass,*mean1slides[icat],*sigma1slides[icat]);
    g2slides[icat] = new RooGaussian(TString("g2slide")+catnames.at(icat),"",*hmass,*mean2slides[icat],*sigma2slides[icat]);
    g3slides[icat] = new RooGaussian(TString("g3slide")+catnames.at(icat),"",*hmass,*mean3slides[icat],*sigma3slides[icat]);
    combhslides[icat] = new RooAddPdf(TString("combhslide")+catnames.at(icat),"",RooArgList(*g1slides[icat],*g2slides[icat],*g3slides[icat]),RooArgList(*fitparmfuncs[icat*nparms+6],*fitparmfuncs[icat*nparms+7]),kTRUE);  
    combhminslides[icat] = new RooAddPdf(TString("combhminslide")+catnames.at(icat),"",RooArgList(*g1slides[icat],*g2slides[icat]),RooArgList(*fitparmfuncs[icat*nparms+6]),kTRUE);  
    
    wg1slides[icat] = new RooGaussian(TString("wg1slide")+catnames.at(icat),"",*hmass,*wmean1slides[icat],*wsigma1slides[icat]);
    wg2slides[icat] = new RooGaussian(TString("wg2slide")+catnames.at(icat),"",*hmass,*wmean2slides[icat],*wsigma2slides[icat]);
    combhwrongslides[icat] = new RooAddPdf(TString("combhwrongslide")+catnames.at(icat),"",RooArgList(*wg1slides[icat],*wg2slides[icat]),RooArgList(*fitparmfuncs[icat*nparms+12]));      
    
    fracrightmodslides[icat] = new RooFormulaVar(TString("fracrightmodslide")+catnames.at(icat),"","TMath::Min(@0*@1,1.0)",RooArgList(nuissancedeltafracright,*fitparmfuncs[icat*nparms+13]));
    combhvtxslides[icat] = new RooAddPdf(TString("combhvtxslide")+catnames.at(icat),"",RooArgList(*combhslides[icat],*combhwrongslides[icat]),RooArgList(*fracrightmodslides[icat]),kTRUE);
    combhvtxminslides[icat] = new RooAddPdf(TString("combhvtxminslide")+catnames.at(icat),"",RooArgList(*combhminslides[icat],*wg1slides[icat]),RooArgList(*fracrightmodslides[icat]),kTRUE);
    //combhvtxsimpleslides[icat] = new RooAddPdf(TString("combhvtxsimpleslide")+catnames.at(icat),"",RooArgList(*g1slides[icat],*wg1slides[icat]),RooArgList(*fracrightmodslides[icat]),kTRUE);
    combhvtxsimpleslides[icat] = g1slides[icat];
    
    
    finalnormslides[icat] = fitparmfuncs[icat*nparms+14];
    //Bool_t usecomplexmodel = catnames.at(icat)=="cat0" || catnames.at(icat)=="cat1"  || catnames.at(icat)=="cat4" || catnames.at(icat)=="cat5" || catnames.at(icat)=="singlecat"; 
    //  Bool_t usecomplexmodel = catnames.at(icat)=="cat0" || catnames.at(icat)=="cat1"  || catnames.at(icat)=="cat4" || catnames.at(icat)=="cat5";
    //Bool_t usecomplexmodel = catnames.at(icat)=="cat0" || catnames.at(icat)=="cat4" || catnames.at(icat)=="singlecat"; 
    //Bool_t usesimplemodel = catnames.at(icat)=="cat2" || catnames.at(icat)=="cat3"  || catnames.at(icat)=="cat6" || catnames.at(icat)=="cat7";
    //Bool_t usesimplemodel = false;

    
    Bool_t usecomplexmodel = complexset.count(catnames.at(icat));
    Bool_t usesimplemodel = simpleset.count(catnames.at(icat));    
    
    if (usecomplexmodel) {
      finalpdfslides[icat] = combhvtxslides[icat];
    }
    else if (usesimplemodel) {
      finalpdfslides[icat] = combhvtxsimpleslides[icat];
    }
    else {
      finalpdfslides[icat] = combhvtxminslides[icat];
    }

  }
  
  //comparison tests of interpolated shape vs actual fit
  for (UInt_t icat=0; icat<catnames.size(); ++icat) {

    //nuissancedeltasmears[icat]->setVal(-smearingv[icat]);
    
    for (UInt_t iparm=0; iparm<fitparms.size(); ++iparm) {
      fitparms.at(iparm)->setVal(fitparmsinit.at(iparm));       
    }

    RooAbsPdf *rightpdf = 0;
    RooAbsPdf *wrongpdf = 0;
    RooAbsPdf *allpdf = 0;
    RooAbsPdf *interpdf = 0;
    //Bool_t usecomplexmodel = catnames.at(icat)=="cat0" || catnames.at(icat)=="cat1"  || catnames.at(icat)=="cat4" || catnames.at(icat)=="cat5" || catnames.at(icat)=="singlecat"; 
   // Bool_t usecomplexmodel = catnames.at(icat)=="cat0" || catnames.at(icat)=="cat4" || catnames.at(icat)=="singlecat";
    //Bool_t usesimplemodel = catnames.at(icat)=="cat2" || catnames.at(icat)=="cat3"  || catnames.at(icat)=="cat6" || catnames.at(icat)=="cat7";
    //Bool_t usesimplemodel = false;
    //Bool_t usecomplexmodel = catnames.at(icat)=="cat0" || catnames.at(icat)=="cat1"  || catnames.at(icat)=="cat4" || catnames.at(icat)=="cat5";
    Bool_t usecomplexmodel = complexset.count(catnames.at(icat));
    Bool_t usesimplemodel = simpleset.count(catnames.at(icat));    
    if (usecomplexmodel) {
      rightpdf = &combh;
      wrongpdf = &combhwrong;
      allpdf = &combhvtx;
      interpdf = combhvtxslides[icat];
      f1.setRange(0.52,1.0);
    }
    else if (usesimplemodel) {
      rightpdf = &g1;
      wrongpdf = &wg1;
      allpdf = &g1;
      interpdf = combhvtxsimpleslides[icat];
    }    
    else {
      rightpdf = &combhmin;
      wrongpdf = &wg1;
      allpdf = &combhvtxmin;
      interpdf = combhvtxminslides[icat];
      sigma1.setVal(2.0);
      sigma2.setVal(3.0);
      f1.setRange(0.0,1.0);
    }
 
 
    mnom.setVal(115);  
    hmass->setRange("higgsrange",100.0,massmax);
    hmass->setRange("plotrange",100.0,115+20.0);
    RooFitResult *fitres = rightpdf->fitTo(*testdsets.at(icat),Strategy(1),Minimizer("Minuit2",""),Minos(kFALSE),Range("higgsrange"),SumW2Error(kTRUE), Save(kTRUE),NumCPU(8));    
    RooFitResult *fitreswrong = wrongpdf->fitTo(*testdsetswrong.at(icat),Strategy(1),Minimizer("Minuit2",""),Minos(kFALSE),Range("higgsrange"),SumW2Error(kTRUE), Save(kTRUE),NumCPU(8));    
    
    double fright = testdsets.at(icat)->sumEntries()/(testdsets.at(icat)->sumEntries()+testdsetswrong.at(icat)->sumEntries());
    fracright.setVal(fright);    
    
    TCanvas *ccompint = new TCanvas;
    TString plotname = TString("inttest") + catnames.at(icat) + TString(".eps");
    RooPlot *hplotcompint = hmass->frame(Bins(100),Range("plotrange"));
    testdsetsall.at(icat)->plotOn(hplotcompint);
    allpdf->plotOn(hplotcompint,RooFit::LineColor(kBlue),NormRange("higgsrange"),Range("higgsrange"));      
    interpdf->plotOn(hplotcompint,RooFit::LineColor(kRed),RooFit::LineStyle(kDashed),NormRange("higgsrange"),Range("higgsrange"));    
    hplotcompint->SetTitle("");
    hplotcompint->Draw();    
    
    TLegend *legmc = new TLegend(0.62,0.75,0.92,0.9);  
    legmc->AddEntry(hplotcompint->getObject(0),"MC","LPE");
    legmc->AddEntry(hplotcompint->getObject(1),"Fit","L");
    legmc->AddEntry(hplotcompint->getObject(2),"Interpolated Fit","L");
    legmc->SetBorderSize(0);
    legmc->SetFillStyle(0);
    legmc->Draw();          
    
    ccompint->SaveAs(plotname);

    
    TCanvas *csmear = new TCanvas;
    TString smearplotname = TString("smeartest") + catnames.at(icat) + TString(".eps");
    RooPlot *hplotsmear = hmass->frame(Bins(100),Range("plotrange"));
    testdsetsall.at(icat)->plotOn(hplotsmear);
    interpdf->plotOn(hplotsmear,RooFit::LineColor(kBlue),NormRange("higgsrange"),Range("higgsrange"));    
    nuissancedeltasmears[icat]->setVal(smearingv[icat]);
    //nuissancedeltasmears[icat]->setVal(0.0);
    interpdf->plotOn(hplotsmear,RooFit::LineColor(kRed),NormRange("higgsrange"),Range("higgsrange"));    
    hplotsmear->Draw();
    
    TLegend *legsmear = new TLegend(0.62,0.75,0.92,0.9);  
    legsmear->AddEntry(hplotsmear->getObject(0),"MC","LPE");
    legsmear->AddEntry(hplotsmear->getObject(1),"Fit","L");
    legsmear->AddEntry(hplotsmear->getObject(2),"Smeared Fit","L");
    legsmear->SetBorderSize(0);
    legsmear->SetFillStyle(0);
    legsmear->Draw();  
    csmear->SaveAs(smearplotname);
    
    nuissancedeltasmears[icat]->setVal(0.0);
    
  }
  
  //r9 migration effect on effacc in categories
  RooRealVar nuissancedeltar9fracbarrel("CMS_hgg_nuissancedeltar9fracbarrel","",1.0,0.1,10.0);
  nuissancedeltar9fracbarrel.setConstant();
  RooRealVar nuissancedeltar9fracmixed("CMS_hgg_nuissancedeltar9fracmixed","",1.0,0.1,10.0);
  nuissancedeltar9fracmixed.setConstant();  
  
  //RooFormulaVar *r9fracbarrel = new RooFormulaVar("r9fracbarrel","","@0*@1/(@1+@2)",RooArgList(nuissancedeltar9fracbarrel,*finalnormslides[0],*finalnormslides[1]));
  //RooFormulaVar *r9fracmixed =  new RooFormulaVar("r9fracmixed","", "@0*@1/(@1+@2)",RooArgList(nuissancedeltar9fracmixed,*finalnormslides[2],*finalnormslides[3]));  
  
  std::vector<RooAbsReal*> nsigcats;
  RooFormulaVar nsigcat0(TString::Format("nsig%s",catnames.at(0).Data()),"","@0*@1",RooArgList(nuissancedeltar9fracbarrel,*finalnormslides[0]));
  RooFormulaVar nsigcat1(TString::Format("nsig%s",catnames.at(1).Data()),"","(1.0-@0)*@1 + @2",RooArgList(nuissancedeltar9fracbarrel,*finalnormslides[0],*finalnormslides[1]));
  RooFormulaVar nsigcat2(TString::Format("nsig%s",catnames.at(2).Data()),"","@0*@1",RooArgList(nuissancedeltar9fracmixed,*finalnormslides[2]));
  RooFormulaVar nsigcat3(TString::Format("nsig%s",catnames.at(3).Data()),"","(1.0-@0)*@1 + @2",RooArgList(nuissancedeltar9fracmixed,*finalnormslides[2],*finalnormslides[3]));

  RooFormulaVar nsigcat4(TString::Format("nsig%s",catnames.at(4).Data()),"","@0",RooArgList(*finalnormslides[4]));
//   RooFormulaVar nsigcat5(TString::Format("nsig%s",catnames.at(5).Data()),"","@0",RooArgList(*finalnormslides[5]));
//   RooFormulaVar nsigcat6(TString::Format("nsig%s",catnames.at(6).Data()),"","@0",RooArgList(*finalnormslides[6]));
  
  
  nsigcats.push_back(&nsigcat0);
  nsigcats.push_back(&nsigcat1);  
  nsigcats.push_back(&nsigcat2);
  nsigcats.push_back(&nsigcat3);
  nsigcats.push_back(&nsigcat4);  
//   nsigcats.push_back(&nsigcat5);  
//   nsigcats.push_back(&nsigcat6);  
  
  
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
    RooAbsPdf *hggpdfsmrel = (RooAbsPdf*)finalpdfslides[icat]->Clone(TString::Format("hggpdfsmrel_%s",catnames.at(icat).Data()));
    RooAbsPdf *hggpdfffabs = (RooAbsPdf*)finalpdfslides[icat]->Clone(TString::Format("hggpdfffabs_%s",catnames.at(icat).Data()));
    RooAbsPdf *hggpdfffrel = (RooAbsPdf*)finalpdfslides[icat]->Clone(TString::Format("hggpdfffrel_%s",catnames.at(icat).Data()));
    RooAbsPdf *hggpdfsm4abs = (RooAbsPdf*)finalpdfslides[icat]->Clone(TString::Format("hggpdfsm4abs_%s",catnames.at(icat).Data()));    
    RooAbsPdf *hggpdfsm4rel = (RooAbsPdf*)finalpdfslides[icat]->Clone(TString::Format("hggpdfsm4rel_%s",catnames.at(icat).Data()));
   
    
    RooFormulaVar *nsigsmabs = new RooFormulaVar(TString::Format("hggpdfsmabs_%s_norm",catnames.at(icat).Data()),"","@0*@1/@2",RooArgList(*procxsecs,*nsigcats[icat],*fsmxsecs));
    RooFormulaVar *nsigsmrel = new RooFormulaVar(TString::Format("hggpdfsmrel_%s_norm",catnames.at(icat).Data()),"","@0*@1*@2",RooArgList(*fsmbrs,*procxsecs,*nsigcats[icat]));
    RooFormulaVar *nsigffabs = new RooFormulaVar(TString::Format("hggpdfffabs_%s_norm",catnames.at(icat).Data()),"","@0*@1/@2",RooArgList(*procxsecs,*nsigcats[icat],*fffxsecs));
    RooFormulaVar *nsigffrel = new RooFormulaVar(TString::Format("hggpdfffrel_%s_norm",catnames.at(icat).Data()),"","@0*@1*@2",RooArgList(*fffbrs,*procxsecs,*nsigcats[icat]));    
    RooFormulaVar *nsigsm4abs = new RooFormulaVar(TString::Format("hggpdfsm4abs_%s_norm",catnames.at(icat).Data()),"","@0*@1/@2",RooArgList(*procxsecssm4,*nsigcats[icat],*fsm4xsecs));
    RooFormulaVar *nsigsm4rel = new RooFormulaVar(TString::Format("hggpdfsm4rel_%s_norm",catnames.at(icat).Data()),"","@0*@1*@2",RooArgList(*fsm4brs,*procxsecssm4,*nsigcats[icat]));
    
    RooExtendPdf *sigpdfsmabs = new RooExtendPdf(TString::Format("sigpdfsmabs%s",catnames.at(icat).Data()),"",*hggpdfsmabs,*nsigsmabs);
    RooExtendPdf *sigpdfsmrel = new RooExtendPdf(TString::Format("sigpdfsmrel%s",catnames.at(icat).Data()),"",*hggpdfsmrel,*nsigsmrel);
    RooExtendPdf *sigpdfffabs = new RooExtendPdf(TString::Format("sigpdfffabs%s",catnames.at(icat).Data()),"",*hggpdfffabs,*nsigffabs);
    RooExtendPdf *sigpdfffrel = new RooExtendPdf(TString::Format("sigpdfffrel%s",catnames.at(icat).Data()),"",*hggpdfffrel,*nsigffrel);
    RooExtendPdf *sigpdfsm4abs = new RooExtendPdf(TString::Format("sigpdfsm4abs%s",catnames.at(icat).Data()),"",*hggpdfsm4abs,*nsigsm4abs);
    RooExtendPdf *sigpdfsm4rel = new RooExtendPdf(TString::Format("sigpdfsm4rel%s",catnames.at(icat).Data()),"",*hggpdfsm4rel,*nsigsm4rel);
    
    w->import(*sigpdfsmabs,RecycleConflictNodes());
    w->import(*sigpdfsmrel,RecycleConflictNodes());
    w->import(*sigpdfffabs,RecycleConflictNodes());
    w->import(*sigpdfffrel,RecycleConflictNodes());    
    w->import(*sigpdfsm4abs,RecycleConflictNodes());
    w->import(*sigpdfsm4rel,RecycleConflictNodes());    
  
    addnorm.add(*nsigsmrel);
    combnorms.push_back(nsigsmrel);
    
    addnormff.add(*nsigffrel);
    combnormsff.push_back(nsigffrel);    
    
  }
  
    
  mnom.setVal(130);
  
  //save everything to file with RooWorkspace

  w->Print();
  w->writeToFile("ubersignalmodel.root") ;
  
  wextra->writeToFile("extra.root");

     
  return;
 
  
  
}

