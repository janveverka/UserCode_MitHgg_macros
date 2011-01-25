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

//==================================================================================================
void plot()
{
  printf("\n Hey Higgs to Gamma Gamma folks.\n\n   Welcome to the plot package.\n\n");
}

void plot(const char *prod, const char *name, const char* title, int logy,
	  double xmin, double xmax, double ymin, double ymax,
	  int nRebin, double lumi)
{

  // use logarithmic scale?
  TCanvas *canvas = MitStyle::MakeCanvas("c","c");
  canvas->SetLogy(logy);

  // define the task with its relevant samples
  char homeDir[1024];
  if (gSystem->Getenv("HOME"))
    sprintf(homeDir,"%s",gSystem->Getenv("HOME"));
  else {
    printf(" HOME not defined. EXIT! (export HOME=/home/$USER ?)\n");
    return;
  } 
  TString home = TString(homeDir);

  // now define sample
  TaskSamples* samples = new TaskSamples(prod,dir.Data());
  samples->SetNameTxt("hgg");
  samples->ReadFile((mitHgg + TString("/config")).Data());
  // plot what we want
  PlotTask   *plotTask = new PlotTask(samples,lumi);
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
  // set the titles
  plotTask->SetAxisTitles(title,"Number of Events");
  plotTask->PlotStack("",name);

  if (logy == 1)
    canvas->SaveAs((TString("png/")+TString(name)+TString("_log.png")).Data());
  else
    canvas->SaveAs((TString("png/")+TString(name)+TString("_lin.png")).Data());

  return;
}
  //  // adjust the default plot styles to our liking
  //  HistStyles *styles   = new HistStyles();
  //  HistStyle  *s        = 0;
  //  // start with a clean slate
  //  styles->Clear();
  //  // add the default Monte Carlo styles
  //  s = styles->AddStyle();
  //  s->SetColor      (kBlack);
  //  s->SetFillStyle  (0);
  //  s = styles->AddStyle();
  //  s->SetColor      (kMagenta);
  //  s->SetFillStyle  (3001);
  //  s = styles->AddStyle();
  //  s->SetColor      (kRed);
  //  s->SetFillStyle  (3004);
  //  s = styles->AddStyle();
  //  s->SetColor      (kYellow);
  //  s->SetFillStyle  (3007);
  //  s = styles->AddStyle();
  //  s->SetColor      (kCyan);
  //  s->SetFillStyle  (3010);
  //  // add the default data style
  //  s = styles->SetDataStyle();
  //  s->SetColor      (kBlue);
  //  s->SetFillStyle  (0);
  //  s->SetMarkerStyle(20);
  //  s->SetMarkerSize (1.0);
  //  // set the new styles as the current styles
  //  plotTask->SetHistStyles(styles);
