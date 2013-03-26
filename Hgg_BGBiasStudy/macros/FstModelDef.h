#ifndef FST_MODEL_DEF
#define FST_MODEL_DEF

bool fstHgg_getNameAndNDF(int mode, int iOrder, TString& label, int& ndf) {
  
  switch(mode) {
  case 1:
    label="Exp";
    ndf = 2*(iOrder-1);
    break;
  case 2:
    label="Pow";
    ndf = 2*(iOrder-1);
    break;
  case 3:
    label="Pol";
    ndf = iOrder;
    break;
  case 4:
    label="Lau";
    ndf = iOrder-1;
    break;
  case 5:
    label="Ber";
    ndf = iOrder;
    break;
  default:
    std::cout<<" Mode not implemenyed."<<std::endl;
    return false;
  }

  return true;
}

fstModel* fstHgg_getModel(RooRealVar* mass, int iModel, int order, int iMod, TString name) {
  
  fstModel* model = NULL;
  
  switch(iModel) {
  case 1:
    model = new fstExpModel(mass,order,iMod,name.Data());    
    break;     						       
  case 2:      						       
    model = new fstPowModel(mass,order,iMod,name.Data());      
    break;     						       
  case 3:      						       
    model = new fstPolModel(mass,order,iMod,name.Data());      
    break;     						       
  case 4:      						       
    model = new fstLauModel(mass,order,iMod,name.Data());      
    break;     						       
  case 5:      						       
    model = new fstBernModel(mass,order,iMod,name.Data());     
    break;     	
  default:
    std::cout<<"  Fitting Mode = "<<iModel<<" not implemented."<<std::endl;
    return NULL;
  }  

  return model;
}

#endif
