//--------------------------------------------------------------------------------------------------
// Perform a plot task using a specified set of samples. Nice and clean.
//
// Authors: Ch.Paus                                                                       (Aug 2010)
//--------------------------------------------------------------------------------------------------
#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
#include <iostream>
#include "TH1.h"
#include "TFile.h"
#include "TSystem.h"
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

TString getEnv(const char* name);

//==================================================================================================
void plot()
{
  printf("\n Hey Higgs to Gamma Gamma folks: WELCOME TO THE PLOT PACKAGE ! \n\n");
  return;
}

TCanvas *plot(const char *name, const char *filename, const char* title, int logy,
	      double xmin, double xmax, double ymin, double ymax,
	      int nRebin, double lumi, TString draw="", TString cut="", int nbins=100,
	      const TH1D *putarget=0, bool stack = kTRUE)
{

  // use logarithmic scale?
  //TCanvas *canvas = MitStyle::MakeCanvas("c","c");
  TCanvas *canvas = new TCanvas;
  canvas->SetLogy(logy);

  // read all environment variables
  TString home   = getEnv("HOME");
  TString mitHgg = getEnv("MIT_HGG_DIR");
  TString hstDir = getEnv("MIT_ANA_HIST");
  TString anaCfg = getEnv("MIT_ANA_CFG");
  TString prdCfg = getEnv("MIT_PROD_CFG");

  // define sample
  TaskSamples* samples = new TaskSamples(prdCfg.Data(),hstDir.Data());
  samples->SetNameTxt(anaCfg.Data());
  samples->ReadFile((mitHgg + TString("/config")).Data());

  // plot what we want
  PlotTask   *plotTask = new PlotTask(samples,lumi);
  plotTask->SetPuTarget(putarget);
  // adjust ranges if needed
  if (xmin!=0)
    plotTask->SetHistXMinimum(xmin);
  if (xmax!=0)
    plotTask->SetHistXMaximum(xmax);
  if (ymin!=0)
    plotTask->SetHistMinimum(ymin);
  if (ymax>0)
    plotTask->SetHistMaximum(ymax);
  // rebinning
  plotTask->SetNRebin     (nRebin);

  plotTask->SetNBins(nbins);
  plotTask->SetDrawExp(draw,cut);

  // set the titles
  plotTask->SetAxisTitles(title,"Number of Events");

  if (stack)
    plotTask->PlotStack("",name);
  else
    plotTask->PlotContributions("",name);

  canvas->SaveAs((TString(filename)+TString(".eps")).Data());

  return canvas;
}

TString getEnv(const char* name)
{
  if (! gSystem->Getenv(name)) {
    printf(" Environment variable: %s  not defined. EXIT!\n",name);
    return TString("");
  } 
  return TString(gSystem->Getenv(name));  
}
