void dumptxt
{
 TString pathAndFile    =TString("/scratch/bendavid/root/hggmvaMar19/appOutput_SMmet_Feb16.root");
  
  TFile* file = TFile::Open(pathAndFile.Data());
  TTree* tup = (TTree*) file->FindObjectAny("MITMVAtuple");

  TCut scancut = "proc==11 && pt1>(mass/3.0) && pt2>(mass/4.0) && pt1>(100.0/3.0) && pt2>(100.0/4.0) && mass > 100. && idmva1>-0.3 && idmva2>-0.3 && mass>100 && mass<180";
  
  tup->Scan("run:lumi:evt:mass:ptgg:pt1/mass:pt2/mass:eta1:eta2:cos(phi1-phi2):vtxprob:res:reswrong:idmva1:idmva2:e1:e2:relresraw1:relresraw2:etcorecaliso1:etcorecaliso2:etcorhcaliso1:etcorhcaliso2:etcortrkiso1:etcortrkiso2:abstrkisocic1:abstrkisocic2:mva",scancut,"colsize=10");
 return;

}
