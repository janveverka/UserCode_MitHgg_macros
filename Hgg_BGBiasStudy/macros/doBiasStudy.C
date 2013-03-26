#include "doFit_AllChi2.C"
#include "doFit_Workspace.C"

TString caseToString(int iModel) {

  TString modelName;

  switch(iModel) {
  case 1:
    modelName="Exp";
    break;
  case 2:
    modelName="Pow";
    break;
  case 3:
    modelName="Pol";
    break;
  case 4:
    modelName="Lau";
    break;
  case 5:
    modelName="Ber";
    break;
  default:
    std::cout<<" Mode not implemenyed."<<std::endl;
    modelName="NaN";
  }

  return modelName;
}


void doBiasStudy(TString catName = "hgg_8TeV_2013moriond_bdt0", TString catLabel = "BDT Class 0") {
  
  // ------------------------------------------------------------
  // Truth models to be tested in F-Test
  std::vector<int> truthModels;

  truthModels.push_back(1); /// Exp
  truthModels.push_back(2); /// Pow
  truthModels.push_back(5); /// Ber
  truthModels.push_back(4); /// Lau

  // in this vector we store the order to be used as truth model for each function family
  std::vector<int> truthOrders;
  truthOrders.resize(truthModels.size());

  // ------------------------------------------------------------
  // Fit models to test, hard coded here... (sorry)
  std::map<int, std::vector<int>* > fitModels;

  std::vector<int>* expOrders = new std::vector<int>();
  expOrders->push_back(1); // 1Exp
  expOrders->push_back(2); // 2Exp
  fitModels.insert( std::pair<int, std::vector<int>* > ( 1, expOrders) );

  std::vector<int>* powOrders = new std::vector<int>();
  powOrders->push_back(1); // 1Pow
  fitModels.insert( std::pair<int, std::vector<int>* > ( 2, powOrders) );

  std::vector<int>* berOrders = new std::vector<int>();
  berOrders->push_back(2); // 2Ber
  berOrders->push_back(3); // 3Ber
  berOrders->push_back(4); // 4Ber
  berOrders->push_back(5); // 5Ber
  berOrders->push_back(6); // 6Ber
  fitModels.insert( std::pair<int, std::vector<int>* > ( 5, berOrders) );

  std::vector<int>* lauOrders = new std::vector<int>();
  lauOrders->push_back(1); // 1Lau
  lauOrders->push_back(2); // 2Lau
  lauOrders->push_back(3); // 3Lau
  fitModels.insert( std::pair<int, std::vector<int>* > ( 4, lauOrders) );


  // ------------------------------------------------------------

  TString texTable = TString::Format("\n\\begin{table}\n\\caption{ Results for the bias in %s.}\n\\centering\n\\begin{tabular}{|c|c|c|c|c|}\n\\hline\n\\multicolumn{5}{|c|}{ {\\bf %s } }\\\\\n\\hline\n\\multicolumn{5}{|c|}{ Test for suitable Truth Models } \\\\\n\\hline\n Truth Model  & df & $NLL_N$ & $\\chi^2(\\Delta NLL_{N+1})$ & $p(\\chi^2 > \chi^2(\\Delta NLL_{N+1}))$ \\\\\n\\hline\n",catLabel.Data(), catLabel.Data());
  
  // this round we loop over the models and find the order for each according to the F-Test
  for(unsigned int iModel = 0; iModel < truthModels.size(); ++iModel)
    texTable += doFit_AllChi2( truthModels[iModel], -1, catName, truthOrders[iModel]);

  texTable += TString("\\hline\n\\multicolumn{5}{|c|}{Test for Maximum Bias} \\\\\n\\hline\n Fit Model");
  for(unsigned int iModel = 0; iModel < truthModels.size(); ++iModel)
    texTable += TString::Format("\t&\t%d%s",truthOrders[iModel],caseToString(truthModels[iModel]).Data());
  texTable += TString("\t\\\\\n\\hline\n");
  
  // now we loop over all the fit models (hard coded, many of them probably useless and stupid)
  for( std::map<int, std::vector<int>* >::iterator iFamily = fitModels.begin(); iFamily != fitModels.end(); ++iFamily ) {
    // for each family loop over the orders hard coded above
    std::vector<int>* orders = iFamily->second;
    for(unsigned int iOrder = 0; iOrder < orders->size(); ++iOrder) {
      int theModel = iFamily->first;
      int theOrder = (*orders)[iOrder];

      texTable += TString::Format("\t%d%s",theOrder,caseToString(theModel).Data());
      
      // now loop over all truth models...
      for(unsigned int iModel = 0; iModel < truthModels.size(); ++iModel) {
	int truthOrder = truthOrders[iModel];

	double maxBias = doFit_Workspace( truthModels[iModel],   // mode fot the GEN funtion
					  truthOrder,
					  theModel,
					  theOrder,
					  catName,
					  true,     // plot all of them...
					  false,    // don't be verbose....
					  5,       // number of mass points to test for bias
					  "/home/fabstoec/Hgg_BGBiasStudy/MoriondHggWorkspaces" ); // directory for the WS

	if( TMath::Abs(maxBias) < 0.2 )
	  texTable += TString::Format("\t&\t%.2f",maxBias);
	else
	  texTable += TString::Format("\t&\t{\\bf%.2f}",maxBias);
	
      }
      texTable += TString("\\\\\n");
    }
  }
  
  texTable += TString("\\hline\n\\end{tabular}\n\\label{tab:XXXFIXMEXXX}\n\\end{table}\n");
  
  std::cout<<" -------------- BEGIN TEX TABLE --------------- "<<std::endl;
  
  std::cout<<std::endl<<texTable<<std::endl;
  
  std::cout<<" -------------- END   TEX TABLE --------------- "<<std::endl;
  
  return;
}
