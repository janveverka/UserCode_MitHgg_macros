//--------------------------------------------------------------------------------------------------
// Perform a plot task using a specified set of samples. Nice and clean.
//
// Authors: C.Paus                                                                        (Aug 2010)
//--------------------------------------------------------------------------------------------------
#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
#include <iostream>
#include "TH1.h"
#include "TFile.h"
#include "TLatex.h"
#include "TPostScript.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "MitPlots/Style/interface/MitStyle.h"
#include "MitPlots/Input/interface/TaskSamples.h"
#include "MitPlots/Plot/interface/PlotTask.h"
#endif

using namespace std;
using namespace mithep;

void plot(const char *name, const char* title, int logy,
          double xmin, double xmax, double ymin, double ymax,
          int nRebin, double lumi, TString draw="", TString cut="", int nbins=100);

//==================================================================================================
void plotsHgg(double lumi = 3000.)
{
  // setup graphics stuff before starting
  MitStyle::Init();
  gROOT->Macro("$CMSSW_BASE/src/MitHgg/macros/plot.C+");

  TString draw = "hmass";
 // TString cut = "pt1>40.0 && pt2>30.0 && (!isconverted1 || (eoverp1<2.0 && eoverp1>0.8 && abs(detaconv1)<0.01 && abs(dphiconv1)<0.01) ) && (!isconverted2 || (eoverp2<2.0 && eoverp2>0.8 && abs(detaconv2)<0.01 && abs(dphiconv2)<0.01) )";
  //TString cut = "pt1>40.0 && pt2>30.0 && isconverted1 && isconverted2";
  //TString cut = "pt1>40.0 && pt2>30.0 && (isconverted1 && (eoverp1<2.0 && eoverp1>0.8 && abs(detaconv1)<0.01 && abs(dphiconv1)<0.01) ) || (isconverted2 && (eoverp2<2.0 && eoverp2>0.8 && abs(detaconv2)<0.01 && abs(dphiconv2)<0.01) )";
// TString cut = "pt1>40.0 && pt2>30.0 && ( (!isconverted1 && !haspixelseed1) || (isconverted1 && eoverp1>0.8 && eoverp1<2.0 && abs(detaconv1)<0.01 && abs(dphiconv1)<0.01) ) &&  ( (!isconverted2 && !haspixelseed2) || (isconverted2 && eoverp2>0.8 && eoverp2<2.0 && abs(detaconv2)<0.01 && abs(dphiconv2)<0.01) )";
 //TString cut = "pt1>40.0 && pt2>30.0 && ( (!isconverted1 && !haspixelseed1) || (isconverted1 && eoverp1>0.8 && eoverp1<2.0 && abs(detaconv1)<0.01 && abs(dphiconv1)<0.01) ) &&  ( (!isconverted2 && !haspixelseed2) || (isconverted2 && eoverp2>0.8 && eoverp2<2.0 && abs(detaconv2)<0.01 && abs(dphiconv2)<0.01) ) && (isconverted1+isconverted2)>0 && !haspixelseed1 && !haspixelseed2";
//TString cut = "pt1>40.0 && pt2>30.0 && (isconverted1 && eoverp1<2.0 && eoverp1>0.8 && abs(detaconv1)<0.01 && abs(dphiconv1)<0.01 && abs(eta1)<1.5) && !haspixelseed2 && abs(eta2)<1.5";
  TString cut = "pt1>40.0 && pt2>30.0 && isconverted2 && !haspixelseed1";



 plot("hHggNtuple","di-photon mass [GeV/c^{2}]",  0, 0.,  10., 0., -1.,1,lumi,"eoverp2",cut,100);
  //plot("hHggNtuple","di-photon mass [GeV/c^{2}]",  0, -0.5,  2.5, 0., -1.,1,lumi,"isconverted1+isconverted2",cut,3);
  //plot("hHggNtuple","di-photon mass [GeV/c^{2}]",  0, 0.,  5., 0., -1.,1,lumi,"abs(eta1-eta2)",cut,20);
  //plot("hHggNtuple","di-photon mass [GeV/c^{2}]",  0, 0.,  100., 0., -1.,1,lumi,"hpt",cut,100);

 // plot("hHggNtuple","di-photon mass [GeV/c^{2}]",  0, 100.,  200., 0., -1.,1,lumi,"hmass",cut,100);


//  plot("hHggNtuple","di-photon mass [GeV/c^{2}]",  0, 100.,  200., 0., -1.,1,lumi,"hmass",cut,100);
  //plot("hHggNtuple","#Delta#eta",  0, -0.5,  2.5., 0., -1.,1,lumi,"isconverted1+isconverted2",cut,3);
  //plot("hHggNtuple","#Delta#eta",  0, -1.0,  1.0., 0., -1.,1,lumi,"cosdbetaboost",cut,100);



  return;
  plot("h2TrigPhotonMass","di-photon mass [GeV/c^{2}]",  0, 0.,  0., 0., -1.,8,lumi);
  return;

  plot("hPhotonEta1",     "photon #eta_{1}",             0, 0.,  0., 0., 50.,2,lumi);
  plot("hPhotonEta2",     "photon #eta_{2}",             0, 0.,  0., 0., 50.,2,lumi);
  plot("hPhotonPhi1",     "photon #phi_{1}",             0, 0.,  0., 0., 75.,4,lumi);
  plot("hPhotonPhi2",     "photon #phi_{2}",             0, 0.,  0., 0., 75.,4,lumi);
  plot("hPhotonDelR",     "di-elec #Delta R",            0, 0.,  0., 0., -1.,2,lumi);
  plot("hPhotonEt1",      "E_{T,1} [GeV]",               0, 0.,200., 0., -1.,1,lumi);
  plot("hPhotonEt2",      "E_{T,2} [GeV]",               0, 0.,200., 0., -1.,1,lumi);
  plot("hPhotonR91",      "R9_{1}",                      0, 0.,  0., 0.,-1.0,2,lumi);
  plot("hPhotonR92",      "R9_{2}",                      0, 0.,  0., 0.,-1.0,2,lumi);
  plot("h2R9PhotonMass",  "di-photon mass [GeV/c^{2}]",  0, 0.,  0., 0.,-1.0,8,lumi);
  plot("h2PhotonPt",      "di-photon p_{T} [GeV/c]",     0, 0.,  0., 0., -1.,4,lumi);
  plot("h2PhotonMass",    "di-photon mass [GeV/c^{2}]",  0, 0.,  0., 0., -1.,8,lumi);
  plot("h2TrigPhotonMass","di-photon mass [GeV/c^{2}]",  0, 0.,  0., 0., -1.,8,lumi);
  //plot("h2SelePhotonMass","di-photon mass [GeV/c^{2}]",  0, 0.,  0., 0., -1.,8,lumi);


  plot("hPhotonEt1",      "E_{T,1} [GeV]",               1, 0.,200., 0.005, -1.,1,lumi);
  plot("hPhotonEt2",      "E_{T,2} [GeV]",               1, 0.,200., 0.005, -1.,1,lumi);
  plot("h2PhotonPt",      "di-photon p_{T} [GeV/c]",     1, 0.,200., 0.005, -1.,2,lumi);

  plot("h2PhotonMass",    "di-photon mass [GeV/c^{2}]",  1, 0.,  0., 0.005, -1.,8,lumi);

  return;
}
