#ifndef FSTMODELS_H
#define FSTMODELS_H

// ===========================================================================================

class fstModel {
  
public:
  fstModel() {};
  
  virtual void resetValues() {};
  virtual void storeValues() {};
  virtual RooAbsPdf* getPdf() = 0;
  
private:
  
};

// =============== BERNSTEIN PDF ============================================================================

class fstBernModel : public fstModel {
public:
  fstBernModel() {};
  fstBernModel(RooRealVar* mass, int order, int number, TString base, TString base2 = "hgg");
  
  virtual void resetValues();
  virtual void storeValues();
  virtual RooAbsPdf* getPdf();
  
private:

  RooRealVar*      _theMass;

  int _enum;
  int _order;
  
  RooBernstein*       _thePdf;
  
  RooConstVar*     _bernconstvar;
  RooRealVar**     _bernval;
  RooFormulaVar**  _bernvalsq;
  RooArgList*      _bernval_list;
  
  double*          _bernval_start;
};

fstBernModel::fstBernModel(RooRealVar* mass, int order, int number, TString base, TString base2):
  _theMass(mass), _enum(number), _order(order) {
  
  // second order polinomial
  _bernval       = new RooRealVar    *[_order];
  _bernvalsq     = new RooFormulaVar *[_order];
  _bernval_start = new double         [_order];

  stringstream pSS_bernName;  
  pSS_bernName.str("");
  pSS_bernName << "fstbern_bernval_list_" << base.Data() << "_" << _enum;
  _bernval_list = new RooArgList(pSS_bernName.str().c_str());  

  pSS_bernName.str("");
  pSS_bernName << "fstbern_bernval_constvar_" << base.Data() << "_" << _enum;
  _bernconstvar = new RooConstVar(pSS_bernName.str().c_str(),"",1.0);

  _bernval_list->add(*_bernconstvar);


  for(int iOrder = 0; iOrder < _order; ++iOrder) {
    pSS_bernName.str("");
    pSS_bernName << "CMS_"<<base2.Data()<<"_" << base.Data() << "_p" << (iOrder+1);
    _bernval[iOrder] = new RooRealVar(pSS_bernName.str().c_str(),"",-10.,10);
    _bernval[iOrder]->setVal((iOrder+1)*0.1);
    //_bernval[iOrder] -> removeRange();
    _bernval_start[iOrder] = _bernval[iOrder]->getVal();

    TString formula = TString(pSS_bernName.str().c_str())+TString("*")+TString(pSS_bernName.str().c_str());

    pSS_bernName.str("");
    pSS_bernName << "fstbern_bernvalsq_" << base.Data() << "_" << _enum << "_" << iOrder;    
    _bernvalsq[iOrder] = new RooFormulaVar(pSS_bernName.str().c_str(),"@0*@0",RooArgList(*_bernval[iOrder]));    
    //_bernvalsq[iOrder] = new RooFormulaVar(pSS_bernName.str().c_str(),formula.Data(),RooArgList(*_bernval[iOrder]));    
    _bernval_list -> add(*_bernvalsq[iOrder]);
  }
  
  pSS_bernName.str("");

  pSS_bernName << "CMS_"<<base2.Data()<<"_" << base.Data() << "_bkgshape";
 _thePdf  = new RooBernstein(pSS_bernName.str().c_str(),"",*_theMass, *_bernval_list);
  
  return;
}

void fstBernModel::resetValues() {

  for(int iOrder = 0; iOrder < _order; ++iOrder) {
    _bernval[iOrder]->setVal(_bernval_start[iOrder]);
  }

  return;
}

void fstBernModel::storeValues() {
  for(int iOrder = 0; iOrder < _order; ++iOrder) {
    _bernval_start[iOrder] = _bernval[iOrder]->getVal();
  }
  return;
}

RooAbsPdf* fstBernModel::getPdf() {
  return (RooAbsPdf*) _thePdf;
}

#endif
