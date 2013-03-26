#include "doFit_Chi2.C"

TString doFit_AllChi2(int fitMode, int maxOrder, TString whichCat, int& truthOrder) {
  
  double lastNll = 0.;
  double chi2    = 0.;
  double pval    = 0.;
  //for(int iOrder = 1; iOrder <= maxOrder; ++iOrder) {
  int iOrder = 1;

  while( maxOrder < 0 || iOrder <= maxOrder ) {
    double thisNll = doFit_Chi2(fitMode, iOrder, whichCat);
    pval = TMath::Prob(TMath::Abs(2*(lastNll-thisNll)), ( fitMode < 3 ? 2 : 1));
    chi2 = TMath::Abs(2*(lastNll-thisNll));
    std::cout<<" p(val) = "<<pval<<std::endl;
    std::cout<<" chi2   = "<<chi2<<std::endl;
    if( pval > 0.05 ) break;
    lastNll = thisNll;
    iOrder++;
  }

  // for the last chosen truth model, redo the fit and plot the bands
  double dummyNll = doFit_Chi2(fitMode, iOrder - 1, whichCat,true);
  truthOrder = (iOrder - 1);

  int numDF = 0;  
  TString modelName = "Exp";
  
  fstHgg_getNameAndNDF(fitMode, iOrder, modelName, numDF);
    
  TString returnString = TString::Format("\t%d%s\t&\t%d\t&\t%.2f\t&\t%.2f\t&\t%.2f\t\\\\\n", (iOrder-1), modelName.Data(), numDF, lastNll, chi2, pval);  

  return returnString;
}
