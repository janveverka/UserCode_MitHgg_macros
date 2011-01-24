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

void plot(const char *prod, const char *name, const char* title, int logy,
	  double xmin, double xmax, double ymin, double ymax,
	  int nRebin, double lumi);


//==================================================================================================
void plotsZmm(const char *prod = "hgg-v0", double lumi = 36.1)
{
  // setup graphics stuff before starting
  MitStyle::Init();
  gROOT->Macro("$CMSSW_BASE/src/MitHgg/macros/plot.C+");
  
  plot(prod,"h2MuonMass",    "di-muon mass [GeV/c^{2}]",1,60.,120.,0.,-1.,2,lumi);
  return;

  plot(prod,"hMuonEta1",     "muon #eta_{1}",           0, 0.,  0., 0., 50.,2,lumi);
  plot(prod,"hMuonEta2",     "muon #eta_{2}",           0, 0.,  0., 0., 50.,2,lumi);
  plot(prod,"hMuonPhi1",     "muon #phi_{1}",           0, 0.,  0., 0., 75.,4,lumi);
  plot(prod,"hMuonPhi2",     "muon #phi_{2}",           0, 0.,  0., 0., 75.,4,lumi);
  plot(prod,"hMuonDelR",     "di-muon #Delta R",        0, 0.,  0., 0., -1.,2,lumi);
  plot(prod,"hMuonEt1",      "E_{T,1} [GeV]",           0, 0.,200., 0., -1.,1,lumi);
  plot(prod,"hMuonEt2",      "E_{T,2} [GeV]",           0, 0.,200., 0., -1.,1,lumi);
  plot(prod,"h2MuonPt",      "di-muon p_{T} [GeV/c]",   0, 0.,200., 0., -1.,2,lumi);
  plot(prod,"h2MuonMass",    "di-muon mass [GeV/c^{2}]",0,60.,120., 0., -1.,2,lumi);
  plot(prod,"h2TrigMuonMass","di-muon mass [GeV/c^{2}]",0,60.,120., 0., -1.,2,lumi);
  plot(prod,"h2SeleMuonMass","di-muon mass [GeV/c^{2}]",0,60.,120., 0., -1.,2,lumi);

  plot(prod,"hMuonEt1",      "E_{T,1} [GeV]",           1, 0.,200., 0.005, -1.,1,lumi);
  plot(prod,"hMuonEt2",      "E_{T,2} [GeV]",           1, 0.,200., 0.005, -1.,1,lumi);
  plot(prod,"h2MuonPt",      "di-muon p_{T} [GeV/c]",   1, 0.,200., 0.005, -1.,2,lumi);

  plot(prod,"h2MuonMass",    "di-muon mass [GeV/c^{2}]",1, 0.,  0., 0.005, -1.,2,lumi);

  return;
}
