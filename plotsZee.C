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
	  int nRebin, double lumi);

//==================================================================================================
void plotsZee(const char *prod = "hgg-v0", double lumi = 36.1)
{
  // setup graphics stuff before starting
  MitStyle::Init();
  gROOT->Macro("$CMSSW_BASE/src/MitHgg/macros/plot.C+");

  plot("h2ElecMass",    "di-electron mass [GeV/c^{2}]",0,60.,120.,0.,-1.,2,lumi);
  return;

  plot("hElecEta1",     "electron #eta_{1}",           0, 0.,  0., 0., 50.,2,lumi);
  plot("hElecEta2",     "electron #eta_{2}",           0, 0.,  0., 0., 50.,2,lumi);
  plot("hElecPhi1",     "electron #phi_{1}",           0, 0.,  0., 0., 75.,4,lumi);
  plot("hElecPhi2",     "electron #phi_{2}",           0, 0.,  0., 0., 75.,4,lumi);
  plot("hElecDelR",     "di-elec #Delta R",            0, 0.,  0., 0., -1.,2,lumi);
  plot("hElecEt1",      "E_{T,1} [GeV]",               0, 0.,200., 0., -1.,1,lumi);
  plot("hElecEt2",      "E_{T,2} [GeV]",               0, 0.,200., 0., -1.,1,lumi);
  //plot("hElecR91",      "R9_{1}",                      0, 0.,  0., 0.,-1.0,2,lumi);
  //plot("hElecR92",      "R9_{2}",                      0, 0.,  0., 0.,-1.0,2,lumi);
  //plot("h2R9ElecMass",  "di-electron mass [GeV/c^{2}]",0, 0.,  0., 0.,-1.0,8,lumi);
  plot("h2ElecPt",      "di-electron p_{T} [GeV/c]",   0, 0.,200., 0., -1.,2,lumi);
  plot("h2ElecMass",    "di-electron mass [GeV/c^{2}]",0,60.,120., 0., -1.,2,lumi);
  plot("h2TrigElecMass","di-electron mass [GeV/c^{2}]",0,60.,120., 0., -1.,2,lumi);
  plot("h2SeleElecMass","di-electron mass [GeV/c^{2}]",0,60.,120., 0., -1.,2,lumi);


  plot("hElecEt1",      "E_{T,1} [GeV]",               1, 0.,200., 0.005, -1.,1,lumi);
  plot("hElecEt2",      "E_{T,2} [GeV]",               1, 0.,200., 0.005, -1.,1,lumi);
  plot("h2ElecPt",      "di-electron p_{T} [GeV/c]",   1, 0.,200., 0.005, -1.,2,lumi);

  plot("h2ElecMass",    "di-electron mass [GeV/c^{2}]",1, 0.,  0., 0.005, -1.,2,lumi);

  return;
}
