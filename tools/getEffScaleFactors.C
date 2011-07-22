using namespace RooFit;

TH1D* puweights = NULL;

static Float_t puweight(float npu) {

  //std::cout<<"#PU = "<<npu<<std::endl;

  if(npu<0) return 1.;
  return puweights->GetBinContent(puweights->FindBin(npu));  
}


void getEffScaleFactors(bool MC=false, bool highR9=true, bool EB=true, bool MIT=true,  TString fileDir = "/scratch/fabstoec/cms/hist/hgg-MCSF/t2mit/filefi/merged") {

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
  //TString base = "puweight(numPU) * ( !(TMath::Abs(phEta)>1.4442 && TMath::Abs(phEta)<1.556) && invMass > 60. && invMass < 120.";
  TString base = " ( !(TMath::Abs(phEta)>1.4442 && TMath::Abs(phEta)<1.556) && invMass > 60. && invMass < 120. ";
  Int_t phCat=1;
  if(EB) {
    base+= TString (" && ( TMath::Abs(phEta) < 1.5)"); 
    labelRight="Barrel Photons E_{T}>20 GeV";
  } else {
    phCat+=2;
    base+= TString (" && ( TMath::Abs(phEta) > 1.5)"); 
    labelRight="Endcap Photons E_{T}>20 GeV";
  }
  if(highR9) {
    base+= TString (" && ( phR9 > 0.94)"); 
    labelRight+=", R9 > 0.94";
  } else {
    phCat++;
    base+= TString (" && ( phR9 < 0.94)"); 
    labelRight+=", R9 < 0.94";
  }
  labelRightPass=labelRight+" Pass";
  labelRightFail=labelRight+" Fail";

  TString probe = " && phEt > 25. ";
  TString pass = " && phTrig > 0 )";
  TString fail = " && phTrig < 0 )";

  RooRealVar invMass("invMass","m_{e-#gamma}  [GeV]",60.,120.);
  invMass.setBins(120);
  invMass.setBins(10000,"cache");
  RooRealVar phEta("phEta","",-10000.,10000.);
  RooRealVar phR9("phR9","",-10000.,10000.);
  RooRealVar elEta("elEta","",-10000.,10000.);
  RooRealVar phEt("phEt","",-10000.,10000.);
  RooRealVar phTrig("phTrig","",-10000.,10000.);
  RooRealVar numPV("numPV","",-10000.,10000.);
  RooRealVar numPU("numPU","",-10000.,10000.);
  
  RooArgSet treeSet(RooArgList(invMass,phEta,elEta,phR9,phEt,phTrig,numPV,numPU));
  
  TString fileName;
  if(MC) {
    fileName = "/hgg-MCSF_s11-zjetsm50-mg-v11-pu_noskim.root";
    labelLeft="CMS Z #rightarrow ee MonteCarlo";
    plotName="eff_MC_";
  }  else {
    fileName = "/hgg-MCSF_r11a-del-m10-v1_noskim.root";
    labelLeft="CMS preliminary 209 pb^{-1}";
    plotName="eff_DATA_";
  }

  TFile* dataFile = new TFile( (fileDir+fileName).Data() );  
  TNtuple* tuple;
  if(MIT)
    tuple = (TNtuple*) dataFile->FindObjectAny("MITtuple");
  else
    tuple = (TNtuple*) dataFile->FindObjectAny("CiCtuple");

  // --------------------------------------------------------------
  // get the stupid PU weights...
  TFile* filepu = new TFile("./tools/PUreweighting/2011_May10thReReco.root","READ");
  puweights=(TH1D*) filepu->Get("pileup");
  // normalize...
  puweights->Scale(1./puweights->Integral());
  
  // get the denominator
  TH1D* puDenom = (TH1D*) dataFile->FindObjectAny("hNPU");
  puDenom->Scale(1./puDenom->Integral());
  
  for(Int_t i=1; i<=puweights->GetNbinsX(); ++i) 
    puweights->SetBinContent(i, puweights->GetBinContent(i)/puDenom->GetBinContent(i));
  // --------------------------------------------------------------
  
  TH1D* dummyMuMuPass = new TH1D("dummyMuMuPass","",120,60.,120.);
  dummyMuMuPass->GetYaxis()->SetTitle("# Events");
  dummyMuMuPass->GetXaxis()->SetTitle("m_{e-#gamma}   [GeV]");
  
  tuple->Draw("invMass >> dummyMuMuPass", (base+probe+pass).Data());
  dummyMuMuPass = (TH1D*) gPad->GetPrimitive("dummyMuMuPass");
  //RooDataHist* dataPass = new RooDataHist("dataPass","",invMass,dummyMuMuPass);
  RooDataSet* dataPass = new RooDataSet("dataPass","",tuple,treeSet,(base+probe+pass).Data());
  
  RooPlot* framePass = invMass.frame();
  //dataPass->plotOn(framePass);
  //framePass->Draw("SAME");

  TH1D* dummyMuMuFail = new TH1D("dummyMuMuFail","",120,60.,120.);
  dummyMuMuFail->GetYaxis()->SetTitle("# Events");
  dummyMuMuFail->GetXaxis()->SetTitle("m_{e-#gamma}   [GeV]");

  tuple->Draw("invMass >> dummyMuMuFail", (base+probe+fail).Data());
  dummyMuMuFail = (TH1D*) gPad->GetPrimitive("dummyMuMuFail");
  //RooDataHist* dataFail = new RooDataHist("dataFail","",invMass,dummyMuMuFail);
  RooDataSet* dataFail = new RooDataSet("dataFail","",tuple,treeSet,(base+probe+fail).Data());  

  RooPlot* frameFail = invMass.frame();
  //dataFail->plotOn(frameFail);

  // try doing the fits binned...
  // efficinecy
  RooRealVar nSIGtot("nSIGtot","",(dataFail->sumEntries()+dataPass->sumEntries()),0.,10000000.);
  //RooRealVar nSIGtot("nSIGtot","",(dummyMuMuFail->Integral()+dummyMuMuPass->Integral()));
  nSIGtot.removeRange();
  RooRealVar theEff("theEff","",0.9,0.,1.);
  theEff.removeRange();
  //RooRealVar theEff("theEff","",0.7);
  RooFormulaVar nSIGP("nSIGP","theEff*nSIGtot", RooArgList(theEff,nSIGtot));
  RooFormulaVar nSIGF("nSIGF","(1.0-theEff)*nSIGtot", RooArgList(theEff,nSIGtot));
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

  RooCBShape SIGP_CB("SIGP_CB","",invMass,SIGP_CBmean,SIGP_CBwidth,alphaP,nP);
  //RooGaussian SIGP_CB("SIGP_CB","",invMass,SIGP_CBmean,SIGP_CBwidth);
  
  // ... the second part of the signal we use a Gaussian.

  //     as mean we again use the mass, as sigma the width from MC.
  RooRealVar SIGP_BWwidth("SIGP_BWwidth","",2.495);
  //RooRealVar SIGP_BWmean("SIGP_BWmean","",91.4, 85.,95.);
  RooRealVar SIGP_BWmean("SIGP_BWmean","",91.188);
  RooBreitWigner SIGP_BW("SIGP_BW","",invMass,SIGP_BWmean,SIGP_BWwidth);
  
  RooFFTConvPdf SIGPmod("SIGPmod","",invMass,SIGP_BW,SIGP_CB);
  SIGPmod.setBufferFraction(2.5);
  
  RooRealVar c1P("c1P","",0.,-100.,100.);
  c1P.removeRange();
  RooRealVar nBGP("nBGP","",1000.,0.,100000.);
  nBGP.removeRange();
  RooPolynomial BGP("BGP","",invMass,RooArgList(c1P));
  //RooChebychev BGP("BGP","",invMass,RooArgList(c1P));


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

  RooCBShape SIGF_CB("SIGF_CB","",invMass,SIGF_CBmean,SIGF_CBwidth,alphaF,nF);
  
  // ... the second part of the signal we use a Gaussian.
  //     as mean we again use the mass, as sigma the width from MC.
  RooRealVar SIGF_BWwidth("SIGF_BWwidth","",2.495);
  RooRealVar SIGF_BWmean("SIGF_BWmean","",91.188);
  RooBreitWigner SIGF_BW("SIGF_BW","",invMass,SIGF_BWmean,SIGF_BWwidth);
  
  RooFFTConvPdf SIGFmod("SIGFmod","",invMass,SIGF_BW,SIGF_CB);
  SIGFmod.setBufferFraction(2.5);  

  RooRealVar c1F("c1F","",0.0051,-100.,100.);
  c1F.removeRange();
  RooRealVar nBGF("nBGF","",690.,0.,100000.);
  //nBGF.removeRange();
  RooPolynomial BGF("BGF","",invMass,RooArgList(c1F));
  //RooChebychev BGF("BGF","",invMass,RooArgList(c1F));
  //=================================================================
  // Generic Erf PDF
  RooRealVar erfP0("erfP0","",0.14,0.001,1.);
  RooRealVar erfP1("erfP1","",-59.,-64.,-10.);
  erfP0.removeRange();
  erfP1.removeRange();
  RooGenericPdf erfPPdf("erfPPdf","","TMath::Erf(erfP0*(invMass+erfP1))",RooArgList(erfP0,erfP1,invMass));

  RooRealVar erfF0("erfF0","",0.14,0.001,1.);
  RooRealVar erfF1("erfF1","",-59.,-64.,-10.);
  erfF0.removeRange();
  erfF1.removeRange();
  RooGenericPdf erfFPdf("erfFPdf","","TMath::Erf(erfF0*(invMass+erfF1))",RooArgList(erfF0,erfF1,invMass));

  
  RooProdPdf finPass("finPass","",erfPPdf,SIGPmod);
  RooProdPdf finFail("finFail","",erfFPdf,SIGFmod);

  RooProdPdf finPass_BG("finPass_BG","",erfPPdf,BGP);
  RooProdPdf finFail_BG("finFail_BG","",erfFPdf,BGF);

  //RooAddPdf totPMod("totPMod","",RooArgList(finPass),RooArgList(nSIGP));
  RooAddPdf totPMod("totPMod","",RooArgList(finPass,finPass_BG),RooArgList(nSIGP,nBGP));

  RooAddPdf totFMod("totFMod","",RooArgList(finFail,finFail_BG),RooArgList(nSIGF,nBGF));
  
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

  std::cout<<dataPassW.isWeighted()<<std::endl;
  std::cout<<dataPassW.numEntries()<<std::endl;
  std::cout<<dataPassW.sumEntries()<<std::endl;

  std::cout<<dataFailW.isWeighted()<<std::endl;
  std::cout<<dataFailW.numEntries()<<std::endl;
  std::cout<<dataFailW.sumEntries()<<std::endl;


//   RooDataSet *dataCombined = new RooDataSet("dataCombined","dataCombined", RooArgList(invMass), RooFit::Index(sample), 
//   					    RooFit::Import("Pass",dataPassW),
//   					    RooFit::Import("Fail",dataFailW),RooFit::WeightVar(*puweightvP),RooFit::WeightVar(*puweightvF));

  RooDataSet *dataCombined = new RooDataSet("dataCombined","dataCombined", RooArgList(invMass,numPU), RooFit::Index(sample), 
  					    RooFit::Import("Pass",*dataPass),
 					    RooFit::Import("Fail",*dataFail));  


  RooRealVar *puweightvA = (RooRealVar*) dataCombined->addColumn(puweightf);
  RooDataSet dataCombinedW("dataCombinedW","dataCombinedW",*(dataCombined->get()),RooFit::Import(*dataCombined),RooFit::WeightVar(*puweightvA));

  std::cout<<dataCombined->isWeighted()<<std::endl;
  std::cout<<dataCombined->numEntries()<<std::endl;
  std::cout<<dataCombined->sumEntries()<<std::endl;

  std::cout<<dataCombinedW.isWeighted()<<std::endl;
  std::cout<<dataCombinedW.numEntries()<<std::endl;
  std::cout<<dataCombinedW.sumEntries()<<std::endl;

  //return;


//   for(int i=0; i<dataCombinedW.numEntries(); ++i){
//     std::cout<<i<<std::endl;
//     std::cout<<dataCombinedW.get(i)->Print()<<std::endl;
//   }


  // PDF for simultaneous fit
  RooSimultaneous totalPdf("totalPdf","totalPdf", sample);
  totalPdf.addPdf(totPMod,"Pass");
  totalPdf.addPdf(totFMod,"Fail");

  totPMod.fitTo(*dataPass);
  totFMod.fitTo(*dataFail);
  totalPdf.fitTo(dataCombinedW, RooFit::Save(true),RooFit::Extended(true), RooFit::PrintLevel(0),RooFit::Strategy(2),RooFit::SumW2Error(false));

  if(!MC) {
    totPMod.plotOn(framePass,LineStyle(kDashed));
    totFMod.plotOn(frameFail,LineStyle(kDashed));
  } else {
    totPMod.plotOn(framePass,LineStyle(kDashed),LineColor(kRed));
    totFMod.plotOn(frameFail,LineStyle(kDashed),LineColor(kRed));
  }

  // TEXT -----------------------------------------------------------------

  TPavesText* result = new TPavesText(0.15, 0.65, 0.4, 0.85, 1,"NDC");
  result->SetTextFont(42);
  result->SetTextSize(0.037);// dflt=28
  result->SetFillColor(0);
  result->SetLineColor(0);
  result->SetTextAlign(12);

  std::stringstream catSt;
  catSt<<"Photon Category "<<phCat;

  std::stringstream plotSt;
  plotSt<<plotName<<"Cat"<<phCat;
  plotName=TString(plotSt.str().c_str());

  std::cout<<"======================================="<<std::endl;
  std::cout<<theEff.getVal()<<std::endl;
  std::cout<<theEff.getError()<<std::endl;
  std::cout<<theEff.hasError()<<std::endl;
  std::cout<<theEff.hasAsymError()<<std::endl;

  std::stringstream catEff;
  char catEffChar[15];
  sprintf(catEffChar,"%.2e #pm %.2e",theEff.getVal()*100.,theEff.getError()*100.);
  catEff<<"#epsilon = ( "<<catEffChar<<" ) %";

  if(MC)
    result->AddText("Monte Carlo");
  else 
    result->AddText("CMS Data");

  result->AddText(catSt.str().c_str());
  result->AddText(catEff.str().c_str());

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

  framePass->Draw("SAME");
  framePass->Draw();
  text2->Draw();
  text1P->Draw();
  result->Draw();

  canPass->SaveAs((TString("plots/")+plotName+TString("_pass.eps")).Data());

  TCanvas* canFail = new TCanvas();
  canFail->cd();
  dummyMuMuFail->Draw();
  frameFail->Draw("SAME");
  text2->Draw();
  text1F->Draw();

  canFail->SaveAs((TString("plots/")+plotName+TString("_fail.eps")).Data());

  return;

}


