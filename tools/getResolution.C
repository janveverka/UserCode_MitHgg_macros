using namespace RooFit;

TH1D* puweights = NULL;

static Float_t puweight(float npu) {

  //std::cout<<"#PU = "<<npu<<std::endl;

  if(npu<0) return 1.;
  return puweights->GetBinContent(puweights->FindBin(npu));  
}


void getResolution(int phCat=1, bool MIT=true,  TString fileDir = "/scratch/fabstoec/cms/hist/hgg-PHRES/t2mit/filefi/merged") {

  RooMsgService::instance().getStream(1).removeTopic(RooFit::Caching);
  RooMsgService::instance().getStream(0).removeTopic(RooFit::Minimization);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::Minimization);
  RooMsgService::instance().getStream(0).removeTopic(RooFit::Plotting);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::Plotting);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::Fitting);
  RooMsgService::instance().getStream(0).removeTopic(RooFit::Eval);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::Eval);

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetLabelColor(1, "XYZ");
  gStyle->SetLabelFont(42, "XYZ");
  gStyle->SetLabelOffset(0.007, "XYZ");
  gStyle->SetLabelSize(0.035, "XYZ");
  
  gStyle->SetTitleColor(1, "XYZ");
  gStyle->SetTitleFont(42, "XYZ");
  gStyle->SetTitleOffset(1., "XYZ");
  gStyle->SetTitleSize(0.04, "XYZ");
  
  gStyle->SetPalette(1);
  
  TString labelRight;
  TString labelLeft;
  TString labelRightPass;
  TString labelRightFail;

  TString plotName;

  // -----------------------------------------------------------
  TString base = " ( mass > 80. && mass < 100. ";
  switch(phCat) {
  case 1:
    base+= TString (" && ph1Cat==1 && ph2Cat==1"); 
    labelRight="Barrel Photons, R9 > 0.94";
    break;
  case 2:
    base+= TString (" && ph1Cat==2 && ph2Cat==2"); 
    labelRight="Barrel Photons, R9 < 0.94";
    break;
  case 3:
    base+= TString (" && ph1Cat==3 && ph2Cat==3"); 
    labelRight="Endcap Photons, R9 > 0.94";
    break;
  case 4:
    base+= TString (" && ph1Cat==4 && ph2Cat==4"); 
    labelRight="Endcap Photons, R9 < 0.94";
    break;
  default:
    return;
  }

  labelRightPass=labelRight+" Data";
  labelRightFail=labelRight+" MC";

  TString probe = " && mass > 25. ";
  TString pass = " && mass > 0 )";
  TString fail = " && mass > 0 )";

  RooRealVar mass("mass","m_{e-#gamma}  [GeV]",80.,100.);
  mass.setBins(80);
  mass.setBins(10000,"cache");
  RooRealVar ph1Cat("ph1Cat","",-10000.,10000.);
  RooRealVar ph2Cat("ph2Cat","",-10000.,10000.);
  RooRealVar numPU("numPU","",-10000.,10000.);
  
  RooArgSet treeSet(RooArgList(mass,ph1Cat,ph2Cat,numPU));
  
  TString fileNameMC;
  TString fileNameDATA;

  labelLeft="CMS Z #rightarrow ee";
  plotName="res_";


  fileNameMC = "/hgg-PHRES_s11-zjetsm50-mg-v11-pu_noskim.root";
  fileNameDATA = "/hgg-PHRES_r11a-del-m10-v1_noskim.root";


  TFile* dataFileMC   = new TFile( (fileDir+fileNameMC).Data() );  
  TFile* dataFileDATA = new TFile( (fileDir+fileNameDATA).Data() );  
  TNtuple* tupleMC;
  TNtuple* tupleDATA;
  if(MIT) {
    tupleMC   = (TNtuple*) dataFileMC->FindObjectAny("MITTuple");
    tupleDATA = (TNtuple*) dataFileDATA->FindObjectAny("MITTuple");
  } else {
    tupleMC   = (TNtuple*) dataFileMC->FindObjectAny("CiCTuple");
    tupleDATA = (TNtuple*) dataFileDATA->FindObjectAny("CiCTuple");
  }

  // --------------------------------------------------------------
  // get the stupid PU weights...
  TFile* filepu = new TFile("./tools/PUreweighting/2011_May10thReReco.root","READ");
  puweights=(TH1D*) filepu->Get("pileup");
  // normalize...
  puweights->Scale(1./puweights->Integral());
  
  // get the denominator
  TH1D* puDenom = (TH1D*) dataFileMC->FindObjectAny("hNPU");
  puDenom->Scale(1./puDenom->Integral());
  
  for(Int_t i=1; i<=puweights->GetNbinsX(); ++i) 
    puweights->SetBinContent(i, puweights->GetBinContent(i)/puDenom->GetBinContent(i));
  // --------------------------------------------------------------
  
  TH1D* dummyMuMuPass = new TH1D("dummyMuMuPass","",80,80.,100.);
  dummyMuMuPass->GetYaxis()->SetTitle("# Events");
  dummyMuMuPass->GetXaxis()->SetTitle("m_{e-#gamma}   [GeV]");
  
  tupleDATA->Draw("mass >> dummyMuMuPass", (base+probe+pass).Data());
  dummyMuMuPass = (TH1D*) gPad->GetPrimitive("dummyMuMuPass");
  RooDataSet* dataPass = new RooDataSet("dataPass","",tupleDATA,treeSet,(base+probe+pass).Data());
  
  RooPlot* framePass = mass.frame();

  TH1D* dummyMuMuFail = new TH1D("dummyMuMuFail","",80,80.,100.);
  dummyMuMuFail->GetYaxis()->SetTitle("# Events");
  dummyMuMuFail->GetXaxis()->SetTitle("m_{e-#gamma}   [GeV]");

  tupleMC->Draw("mass >> dummyMuMuFail", (base+probe+fail).Data());
  dummyMuMuFail = (TH1D*) gPad->GetPrimitive("dummyMuMuFail");
  RooDataSet* dataFail = new RooDataSet("dataFail","",tupleMC,treeSet,(base+probe+fail).Data());  

  RooPlot* frameFail = mass.frame();

  // try doing the fits binned...
  // efficinecy
  //RooRealVar nSIGtot("nSIGtot","",(dataFail->sumEntries()+dataPass->sumEntries()),0.,10000000.);
  //RooRealVar nSIGtot("nSIGtot","",(dummyMuMuFail->Integral()+dummyMuMuPass->Integral()));
  //nSIGtot.removeRange();

  RooRealVar delM("delM","",0.,-1,1.);
  delM.removeRange();

  RooRealVar nSIGP("nSIGP","", dataPass->sumEntries(),0,10000000.);
  RooRealVar nSIGF("nSIGF","", dataFail->sumEntries(),0,10000000.);

  nSIGP.removeRange();
  nSIGF.removeRange();

  // ---------------------------------------------------------------------
  // FAIL MODEL
  
  // SIGNAL MODEL
  RooRealVar SIGF_CBmean("SIGF_CBmean","",-0.0112,-10.,10.);
  SIGF_CBmean.removeRange();

  RooRealVar SIGF_CBwidth("SIGF_CBwidth","",1.58, -10., 10.);
  SIGF_CBwidth.removeRange();

  RooRealVar nF("nF","",1.15,0.0,1000.);
  nF.removeRange();

  RooRealVar alphaF("alphaF","",0.762, 0.,10.);
  alphaF.removeRange();

  RooCBShape SIGF_CB("SIGF_CB","",mass,SIGF_CBmean,SIGF_CBwidth,alphaF,nF);
  
  // ... the second part of the signal we use a Gaussian.
  //     as mean we again use the mass, as sigma the width from MC.
  RooRealVar SIGF_BWwidth("SIGF_BWwidth","",2.495);
  RooRealVar SIGF_BWmean("SIGF_BWmean","",91.188);
  RooBreitWigner SIGF_BW("SIGF_BW","",mass,SIGF_BWmean,SIGF_BWwidth);
  
  RooFFTConvPdf SIGFmod("SIGFmod","",mass,SIGF_BW,SIGF_CB);
  SIGFmod.setBufferFraction(2.5);  

  RooRealVar c1F("c1F","",0.0051,-100.,100.);
  c1F.removeRange();
  RooRealVar nBGF("nBGF","",690.,0.,100000.);
  //nBGF.removeRange();
  RooPolynomial BGF("BGF","",mass,RooArgList(c1F));
  //RooChebychev BGF("BGF","",mass,RooArgList(c1F));

  // --------------- PASS VARIABLES ------------------
  // SIGNAL MODEL
  RooRealVar SIGP_CBmean("SIGP_CBmean","",0.,-10.,10.);
  SIGP_CBmean.removeRange();
  //RooRealVar SIGP_CBmean("SIGP_CBmean","",0.);

  RooRealVar SIGP_CBwidth("SIGP_CBwidth","",1.64,-10.,10.);
  SIGP_CBwidth.removeRange();
  
  //RooRealVar nP("nP","",2.11);
  RooRealVar nP("nP","",2.87,0.0,1000.);
  nP.removeRange();

  //RooRealVar alphaP("alphaP","",2.);
  RooRealVar alphaP("alphaP","",1.053,0.,10.);
  alphaP.removeRange();

  RooCBShape SIGP_CB("SIGP_CB","",mass,SIGP_CBmean,SIGP_CBwidth,alphaP,nP);
  //RooGaussian SIGP_CB("SIGP_CB","",mass,SIGP_CBmean,SIGP_CBwidth);
  
  // ... the second part of the signal we use a Gaussian.

  //     as mean we again use the mass, as sigma the width from MC.
  RooRealVar SIGP_BWwidth("SIGP_BWwidth","",2.495);
  //RooRealVar SIGP_BWmean("SIGP_BWmean","",91.4, 85.,95.);
  RooRealVar SIGP_BWmean("SIGP_BWmean","",91.188);
  RooBreitWigner SIGP_BW("SIGP_BW","",mass,SIGP_BWmean,SIGP_BWwidth);
  
  RooFFTConvPdf SIGPmod("SIGPmod","",mass,SIGP_BW,SIGP_CB);
  SIGPmod.setBufferFraction(2.5);
  
  RooRealVar c1P("c1P","",0.,-100.,100.);
  c1P.removeRange();
  RooRealVar nBGP("nBGP","",1000.,0.,100000.);
  nBGP.removeRange();
  RooPolynomial BGP("BGP","",mass,RooArgList(c1P));
  //RooChebychev BGP("BGP","",mass,RooArgList(c1P));

  //=================================================================
  // Generic Erf PDF
  RooRealVar erfP0("erfP0","",0.14,0.001,1.);
  RooRealVar erfP1("erfP1","",-59.,-64.,-10.);
  erfP0.removeRange();
  erfP1.removeRange();
  RooGenericPdf erfPPdf("erfPPdf","","TMath::Erf(erfP0*(mass+erfP1))",RooArgList(erfP0,erfP1,mass));

  RooRealVar erfF0("erfF0","",0.14,0.001,1.);
  RooRealVar erfF1("erfF1","",-59.,-64.,-10.);
  erfF0.removeRange();
  erfF1.removeRange();
  RooGenericPdf erfFPdf("erfFPdf","","TMath::Erf(erfF0*(mass+erfF1))",RooArgList(erfF0,erfF1,mass));

  
  RooProdPdf finPass("finPass","",erfPPdf,SIGPmod);
  RooProdPdf finFail("finFail","",erfFPdf,SIGFmod);

  RooProdPdf finPass_BG("finPass_BG","",erfPPdf,BGP);
  RooProdPdf finFail_BG("finFail_BG","",erfFPdf,BGF);

  //RooAddPdf totPMod("totPMod","",RooArgList(finPass),RooArgList(nSIGP));
  //RooAddPdf totPMod("totPMod","",RooArgList(finPass,finPass_BG),RooArgList(nSIGP,nBGP));
  //RooAddPdf totFMod("totFMod","",RooArgList(finFail,finFail_BG),RooArgList(nSIGF,nBGF));

  //RooAddPdf totPMod("totPMod","",RooArgList(finPass),RooArgList(nSIGP));
  //RooAddPdf totFMod("totFMod","",RooArgList(finFail),RooArgList(nSIGF));

  RooAddPdf totPMod("totPMod","",RooArgList(SIGPmod),RooArgList(nSIGP));
  RooAddPdf totFMod("totFMod","",RooArgList(SIGFmod),RooArgList(nSIGF));
  
   // combine Data... or die trying
  RooCategory sample("sample","") ;
  sample.defineType("Pass", 1) ;
  sample.defineType("Fail", 2) ;

  RooFormulaVar puweightf("puweightv","puweightv","puweight(numPU)",RooArgList(numPU));
  RooRealVar *puweightvP = (RooRealVar*) dataPass->addColumn(puweightf);
  RooRealVar *puweightvF = (RooRealVar*) dataFail->addColumn(puweightf);

  RooDataSet dataPassW("dataPassW","dataPassW",*(dataPass->get()),RooFit::Import(*dataPass),RooFit::WeightVar(*puweightvP));
  RooDataSet dataFailW("dataFailW","dataFailW",*(dataFail->get()),RooFit::Import(*dataFail),RooFit::WeightVar(*puweightvF));

  dataPassW.plotOn(framePass);
  dataFailW.plotOn(frameFail);

  RooDataSet *dataCombined = new RooDataSet("dataCombined","dataCombined", RooArgList(mass,numPU), RooFit::Index(sample), 
  					    RooFit::Import("Pass",*dataPass),
 					    RooFit::Import("Fail",*dataFail));  


  RooRealVar *puweightvA = (RooRealVar*) dataCombined->addColumn(puweightf);
  RooDataSet dataCombinedW("dataCombinedW","dataCombinedW",*(dataCombined->get()),RooFit::Import(*dataCombined),RooFit::WeightVar(*puweightvA));

  // PDF for simultaneous fit
  RooSimultaneous totalPdf("totalPdf","totalPdf", sample);
  totalPdf.addPdf(totPMod,"Pass");
  totalPdf.addPdf(totFMod,"Fail");


  // definr the formulas for our results
  RooFormulaVar delE("delE","(SIGP_CBmean-SIGF_CBmean)/(91.188+SIGF_CBmean)", RooArgList(SIGP_CBmean,SIGF_CBmean));
  RooFormulaVar delS2("delS2","2*(SIGP_CBwidth*SIGP_CBwidth-SIGF_CBwidth*SIGF_CBwidth)/(91.188+SIGF_CBmean)/(91.188+SIGF_CBmean)", RooArgList(SIGP_CBwidth,SIGF_CBwidth,SIGF_CBmean));

  totPMod.fitTo(dataPassW, RooFit::Save(true),RooFit::Extended(true), RooFit::PrintLevel(0),RooFit::Strategy(2),RooFit::SumW2Error(false));
  totFMod.fitTo(dataFailW, RooFit::Save(true),RooFit::Extended(true), RooFit::PrintLevel(0),RooFit::Strategy(2),RooFit::SumW2Error(false));

  RooFitResult* fitR = totalPdf.fitTo(dataCombinedW, RooFit::Save(true),RooFit::Extended(true), RooFit::PrintLevel(0),RooFit::Strategy(2),RooFit::SumW2Error(false));

  totPMod.plotOn(framePass,LineStyle(kDashed));
  totFMod.plotOn(frameFail,LineStyle(kDashed),LineColor(kRed));

  // TEXT -----------------------------------------------------------------

  TPavesText* result = new TPavesText(0.135, 0.65, 0.4, 0.85, 1,"NDC");
  result->SetTextFont(42);
  result->SetTextSize(0.035);// dflt=28
  result->SetFillColor(0);
  result->SetLineColor(0);
  result->SetTextAlign(12);

  std::stringstream catSt;
  catSt<<"Photon Category "<<phCat;

  std::stringstream plotSt;
  plotSt<<plotName<<"Cat"<<phCat;
  plotName=TString(plotSt.str().c_str());


  std::stringstream catEff;
  char catEffChar[15];
  sprintf(catEffChar,"%.2f #pm %.2f",delE.getVal()*100.,delE.getPropagatedError(*fitR)*100.);
  catEff<<"#DeltaE/E_{MC} = ( "<<catEffChar<<" ) %";

  std::stringstream catEff2;
  char catEffChar2[15];

  catEff2<<"#Delta#sigma/E_{MC} = ( ";

  double sM = delS2.getVal();
  double sME = TMath::Sqrt(TMath::Abs(delS2.getPropagatedError(*fitR)))*100.;
  if(sM < 0.)
    catEff2<<"-";
  sM=TMath::Abs(sM);
  sM=TMath::Sqrt(sM)*100.;
  sprintf(catEffChar2,"%.2f #pm %.2f",sM,sME);

  catEff2<<catEffChar2<<" ) %";

  result->AddText(catSt.str().c_str());
  result->AddText(catEff.str().c_str());
  result->AddText(catEff2.str().c_str());

  TLatex* text1P=new TLatex(3.5,23.5,labelRightPass.Data());
  text1P->SetNDC();
  text1P->SetTextAlign(33);
  text1P->SetX(0.9);//(0.940);
  text1P->SetY(0.95);
  text1P->SetTextFont(42);
  text1P->SetTextSize(0.035);// dflt=28

  TLatex* text1F=new TLatex(3.5,23.5,labelRightFail.Data());
  text1F->SetNDC();
  text1F->SetTextAlign(33);
  text1F->SetX(0.9);//(0.940);
  text1F->SetY(0.95);
  text1F->SetTextFont(42);
  text1F->SetTextSize(0.035);// dflt=28
  
  TLatex* text2=new TLatex(3.5,23.5,labelLeft.Data());
  text2->SetNDC();
  text2->SetTextAlign(13);
  text2->SetX(0.1);//(0.940);
  text2->SetY(0.95);
  text2->SetTextFont(42);
  text2->SetTextSize(0.035);// dflt=28


  TCanvas* canPass = new TCanvas();
  canPass->cd();
  dummyMuMuPass->Draw();

  //framePass->Draw("SAME");
  framePass->Draw();
  text2->Draw();
  text1P->Draw();

  canPass->SaveAs((TString("plots/")+plotName+TString("_Data.eps")).Data());

  TCanvas* canFail = new TCanvas();
  canFail->cd();
  dummyMuMuFail->Draw();
  //frameFail->Draw("SAME");
  frameFail->Draw();
  text2->Draw();
  text1F->Draw();
  result->Draw();


  canFail->SaveAs((TString("plots/")+plotName+TString("_MC.eps")).Data());

  return;

}


