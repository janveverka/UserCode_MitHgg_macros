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

// =============== EXPONENTIAL PDF ============================================================================

class fstErfExpModel : public fstModel {
public:
  fstErfExpModel() {};
  fstErfExpModel(RooRealVar* mass, int order, int number, TString base);
  
  virtual void resetValues();
  virtual void storeValues();
  virtual RooAbsPdf* getPdf();

  void setErfMean(double);
  
private:
  
  RooRealVar*      _theMass;
  
  int _enum;
  int _order;
  
  RooGenericPdf*       _thePdf;
  
  RooRealVar*     _expoval1;
  RooRealVar*     _expoval2;
  
  RooRealVar*     _exporatio;

  RooRealVar*     _erfMean;
  RooRealVar*     _erfWidth;
  
  double*          _expoval_start;

  double           _exporatio_start;

  double          _erfMean_start;
  double          _erfWidth_start;
  

};

void fstErfExpModel::setErfMean(double mean) {
  _erfMean->setVal(mean);
  _erfMean->setConstant();
  return;
}

fstErfExpModel::fstErfExpModel(RooRealVar* mass, int order, int number, TString base):
  _theMass(mass), _enum(number), _order(order) {
  
  _expoval_start   = new double[2];

  _exporatio_start = 1.0;

  _erfMean_start   = 115.;
  _erfWidth_start  = 1.;
  
  _expoval1 = new RooRealVar(TString::Format("fstexp_erfexpoval_1_%s_%d",base.Data(),_enum),"",-100.,100);
  _expoval2 = new RooRealVar(TString::Format("fstexp_erfexpoval_2_%s_%d",base.Data(),_enum),"",-100.,100);

  _exporatio = new RooRealVar(TString::Format("fstexp_erfexporatio_2_%s_%d",base.Data(),_enum),"",0.,1);
  
  _expoval1->setVal(-0.027);
  _expoval2->setVal(-0.005);

  //_expoval2->setConstant();

  _exporatio->setVal(0.5);
  //_exporatio->setConstant();

  _erfMean  = new RooRealVar(TString::Format("fstexp_erfmean_%s_%d",base.Data(),_enum),"",110.,150.);
  _erfWidth  = new RooRealVar(TString::Format("fstexp_erfwidth_%s_%d",base.Data(),_enum),"",-100.,100);
  
  _erfMean->setVal(_erfMean_start);
  _erfWidth->setVal(_erfWidth_start);
  
  _erfMean->setConstant();
  //_erfWidth->setConstant();

  _thePdf  = new RooGenericPdf(TString::Format("fsterfexp_%s_%d",base.Data(),_enum),"",
  			       "@5*TMath::Erf((@0-@1)/TMath::Abs(@2))*TMath::Exp(@3*@0)+(1.-@5)*(1.-TMath::Erf((@0-@1)/TMath::Abs(@2)))*TMath::Exp(@4*@0)",
			       RooArgList(*_theMass,*_erfMean,*_erfWidth,*_expoval1, *_expoval2, *_exporatio));
  
  
  //   _thePdf  = new RooGenericPdf(TString::Format("fsterfexp_%s_%d",base.Data(),_enum),"",
// 			       "TMath::Exp(@1*@0)",
// 			       RooArgList(*_theMass,*_expoval1));

  return;
}

void fstErfExpModel::resetValues() {
  return;
}

void fstErfExpModel::storeValues() {
  return;
}

RooAbsPdf* fstErfExpModel::getPdf() {
  return (RooAbsPdf*) _thePdf;
}


class fstExpModel : public fstModel 
{
public:
  fstExpModel() {};
  fstExpModel(RooRealVar* mass, int order, int number, TString base);
  
  virtual void resetValues();
  virtual void storeValues();
  virtual RooAbsPdf* getPdf();
  
private:

  RooRealVar*      _theMass;

  int _enum;
  int _order;
  
  RooAddPdf*       _thePdf;
  
  RooExponential** _expo;
  RooRealVar**     _expoval;
  RooRealVar**     _exporatio;
  
  RooArgList*      _expo_list;
  RooArgList*      _exporatio_list;
  
  double*          _expoval_start;
  double*          _exporatio_start;
};

fstExpModel::fstExpModel(RooRealVar* mass, int order, int number, TString base):
  _theMass(mass), _enum(number), _order(order) {
  
  _expo      = new RooExponential*[_order];
  _expoval   = new RooRealVar    *[_order];
  _exporatio = new RooRealVar    *[_order - 1];
  
  _expoval_start   = new double[_order];
  _exporatio_start = new double[_order - 1];
  
  stringstream pSS_expName;  
  pSS_expName.str("");
  pSS_expName << "fstexp_expo_list_" << base.Data() << "_" << _enum;
  _expo_list = new RooArgList(pSS_expName.str().c_str());
  pSS_expName.str("");
  pSS_expName << "fstexp_exporatio_list_" << base.Data() << "_" << _enum;
  _exporatio_list = new RooArgList(pSS_expName.str().c_str());

  for(int iOrder = 0; iOrder < _order; ++iOrder) {
    pSS_expName.str("");
    pSS_expName << "fstexp_expoval_" << base.Data() << "_" << _enum << "_" << iOrder;
    _expoval[iOrder] = new RooRealVar(pSS_expName.str().c_str(),"",-100.,100);
    _expoval[iOrder] -> setVal( -(double)(iOrder+1)/10.);
    _expoval[iOrder] -> removeRange();
    _expoval_start[iOrder] = _expoval[iOrder]->getVal();

    if( iOrder < _order - 1 ) {
      pSS_expName.str("");
      pSS_expName << "fstexp_exporatio_" << base.Data() << "_" << _enum << "_" << iOrder;
      _exporatio[iOrder] = new RooRealVar(pSS_expName.str().c_str(),"",-100.,100.);
      _exporatio[iOrder] -> setVal ( 1./(double) order );
      _exporatio[iOrder] -> removeRange();
      _exporatio_start[iOrder] = _exporatio[iOrder]->getVal();
      _exporatio_list -> add(*_exporatio[iOrder]);
    }
    
    pSS_expName.str("");
    pSS_expName << "fstexp_expo_" << base.Data() << "_" << _enum << "_" << iOrder;
    _expo[iOrder] = new RooExponential(pSS_expName.str().c_str(),"",*_theMass,*_expoval[iOrder]);
    _expo_list    -> add(*_expo[iOrder]);
  }
  
  pSS_expName.str("");
  pSS_expName << "fstexp_theNexp_" << base.Data() << "_" << _enum;
  _thePdf  = new RooAddPdf(pSS_expName.str().c_str(),"",*_expo_list, *_exporatio_list);

  return;
}

void fstExpModel::resetValues() {

  for(int iOrder = 0; iOrder < _order; ++iOrder) {
    _expoval[iOrder]->setVal(_expoval_start[iOrder]);
    if( iOrder < _order - 1 ) {      
      _exporatio[iOrder]->setVal(_exporatio_start[iOrder]);
    }
  }

  return;
}

void fstExpModel::storeValues() {
  for(int iOrder = 0; iOrder < _order; ++iOrder) {
    _expoval_start[iOrder] = _expoval[iOrder]->getVal();
    if( iOrder < _order - 1 ) {
      _exporatio_start[iOrder] = _exporatio[iOrder]->getVal();
    }
  }
  return;
}

RooAbsPdf* fstExpModel::getPdf() {
  return (RooAbsPdf*) _thePdf;
}

// =============== CHEBYCHEV PDF ============================================================================

class fstPolModel : public fstModel {
public:
  fstPolModel() {};
  fstPolModel(RooRealVar* mass, int order, int number, TString base);
  
  virtual void resetValues();
  virtual void storeValues();
  virtual RooAbsPdf* getPdf();
  
  void setAndFixValues(std::vector<double>& values);
  
private:
  
  RooRealVar*      _theMass;
  
  int _enum;
  int _order;
  
  //RooAbsPdf*       _thePdf;
  RooChebychev*        _thePdf;
  
  RooRealVar**     _polval;
  RooArgList*      _polval_list;
  
  double*          _polval_start;
};

fstPolModel::fstPolModel(RooRealVar* mass, int order, int number, TString base):
  _theMass(mass), _enum(number), _order(order) {
  
  // second order polinomial
  _polval       = new RooRealVar    *[_order];
  _polval_start = new double         [_order];

  stringstream pSS_polName;  
  pSS_polName.str("");
  pSS_polName << "fstpol_polval_list_" << base.Data() << "_" << _enum;
  _polval_list = new RooArgList(pSS_polName.str().c_str());
  
  for(int iOrder = 0; iOrder < _order; ++iOrder) {
    pSS_polName.str("");
    pSS_polName << "fstpol_polval_" << base.Data() << "_" << _enum << "_" << iOrder;
    _polval[iOrder] = new RooRealVar(pSS_polName.str().c_str(),"",-100.,100);
    _polval[iOrder] -> removeRange();
    _polval_start[iOrder] = _polval[iOrder]->getVal();
    _polval_list -> add(*_polval[iOrder]);
  }
  
  pSS_polName.str("");
  pSS_polName << "fstpol_theNpol_" << base.Data() << "_" << _enum;
  _thePdf  = new RooChebychev(pSS_polName.str().c_str(),"",*_theMass, *_polval_list);
  //_thePdf  = new RooPolynomial(pSS_polName.str().c_str(),"",*_theMass, *_polval_list);
  
  return;
}

void fstPolModel::resetValues() {

  for(int iOrder = 0; iOrder < _order; ++iOrder) {
    _polval[iOrder]->setVal(_polval_start[iOrder]);
  }

  return;
}

void fstPolModel::storeValues() {
  for(int iOrder = 0; iOrder < _order; ++iOrder) {
    _polval_start[iOrder] = _polval[iOrder]->getVal();
  }
  return;
}

void fstPolModel::setAndFixValues(std::vector<double>& values) {
  if( (int) values.size() != _order ) return;
  else {
    for(unsigned int iO = 0; iO < values.size(); ++iO){
      _polval_start[iO] = values[iO];
      _polval[iO]       ->setVal(values[iO]);
      _polval[iO]       ->setConstant();
    }
  }
  return;
}

RooAbsPdf* fstPolModel::getPdf() {
  return (RooAbsPdf*) _thePdf;
}

// =============== LAURENT PDF ============================================================================

class fstLauModel : public fstModel {
public:
  fstLauModel() {};
  fstLauModel(RooRealVar* mass, int order, int number, TString base);
  
  virtual void resetValues();
  virtual void storeValues();
  virtual RooAbsPdf* getPdf();

private:

  RooRealVar*      _theMass;

  int _enum;
  int _order;
  
  RooGenericPdf*       _thePdf;
  
  RooRealVar**     _lauratio;
  RooArgList*      _lau_list;
  
  double*          _lauratio_start;
};

fstLauModel::fstLauModel(RooRealVar* mass, int order, int number, TString base):
  _theMass(mass), _enum(number), _order(order) {
  
  // second order polinomial
  _lauratio       = new RooRealVar    *[_order - 1];
  _lauratio_start = new double         [_order - 1];

  stringstream pSS_lauName;  
  pSS_lauName.str("");
  pSS_lauName << "fstlau_lau_list_" << base.Data() << "_" << _enum;
  _lau_list = new RooArgList(pSS_lauName.str().c_str());
  
  int massIndex = ( _order - 1 );

  stringstream lauFunc;

  double lauPower = -4.;
  
  for(int iOrder = 0; iOrder < _order; ++iOrder) {
  
    if( (iOrder%2) )
      lauPower -= (double) iOrder;
    else
      lauPower += (double) iOrder;

    if( iOrder > 0 ){
      pSS_lauName.str("");
      pSS_lauName << "fstlau_lauratio_" << base.Data() << "_" << _enum << "_" << iOrder;
      _lauratio[iOrder-1] = new RooRealVar(pSS_lauName.str().c_str(),"",-100.,100.);
      _lauratio[iOrder-1] -> setVal ( (double) (iOrder) / 10. );
      _lauratio[iOrder-1] -> removeRange();
      _lauratio_start[iOrder-1] = _lauratio[iOrder-1]->getVal();
      _lau_list -> add(*_lauratio[iOrder-1]);
      lauFunc << "+@" << (iOrder - 1) << "*pow(@"<<massIndex<<","<<lauPower<<")";
    } else
      lauFunc << "pow(@"<<massIndex<<","<<lauPower<<")";    
  }
  
  _lau_list -> add(*_theMass);
  
  pSS_lauName.str("");
  pSS_lauName << "fstlau_theNlau_" << base.Data() << "_" << _enum;
  //std::cout<<lauFunc.str().c_str()<<std::endl;
  _thePdf  = new RooGenericPdf(pSS_lauName.str().c_str(),"",lauFunc.str().c_str(),*_lau_list);
  
  return;
}

void fstLauModel::resetValues() {

  for(int iOrder = 0; iOrder < _order - 1; ++iOrder) {
    _lauratio[iOrder]->setVal(_lauratio_start[iOrder]);
  }

  return;
}

void fstLauModel::storeValues() {
  for(int iOrder = 0; iOrder < _order - 1; ++iOrder) {
      _lauratio_start[iOrder] = _lauratio[iOrder]->getVal();
  }
  return;
}

RooAbsPdf* fstLauModel::getPdf() {
  return (RooAbsPdf*) _thePdf;
}

// =============== POWER PDF ============================================================================

class fstPowModel : public fstModel {
public:
  fstPowModel() {};
  fstPowModel(RooRealVar* mass, int order, int number, TString base);
  
  virtual void resetValues();
  virtual void storeValues();
  virtual RooAbsPdf* getPdf();
  
private:

  RooRealVar*      _theMass;

  int _enum;
  int _order;
  
  RooGenericPdf*       _thePdf;
  
  RooRealVar**     _powval;
  RooRealVar**     _powratio;
  RooArgList*      _pow_list;
  
  double*          _powval_start;
  double*          _powratio_start;
};

fstPowModel::fstPowModel(RooRealVar* mass, int order, int number, TString base):
  _theMass(mass), _enum(number), _order(order) {
  
  // second order polinomial
  _powval         = new RooRealVar    *[_order];
  _powratio       = new RooRealVar    *[_order - 1];

  _powval_start   = new double         [_order];
  _powratio_start = new double         [_order - 1];

  stringstream pSS_powName;  
  pSS_powName.str("");
  pSS_powName << "fstpow_pow_list_" << base.Data() << "_" << _enum;
  _pow_list = new RooArgList(pSS_powName.str().c_str());
  
  int massIndex = (_order - 1)*2 + 1;

  stringstream powFunc;
  
  for(int iOrder = 0; iOrder < _order; ++iOrder) {
    
    if( iOrder > 0 ){
      pSS_powName.str("");
      pSS_powName << "fstpow_powratio_" << base.Data() << "_" << _enum << "_" << iOrder;
      _powratio[iOrder-1] = new RooRealVar(pSS_powName.str().c_str(),"",-100.,100.);
      _powratio[iOrder-1] ->setVal( 1./ (double) order);
      _powratio[iOrder-1] -> removeRange();
      _powratio_start[iOrder-1] = _powval[iOrder-1]->getVal();
      _pow_list -> add(*_powratio[iOrder-1]);
      powFunc << "+@" << (iOrder - 1)*2+1<< "*TMath::Power(@"<<massIndex<<",@"<<iOrder *2<<")";
    } else {
      powFunc << "TMath::Power(@"<<massIndex<<",@"<<iOrder<<")";
      //powFunc << "TMath::Exp(@"<<massIndex<<"*@"<<iOrder<<")";
    }
    
    pSS_powName.str("");
    pSS_powName << "fstpow_powval_" << base.Data() << "_" << _enum << "_" << iOrder;
    _powval[iOrder] = new RooRealVar(pSS_powName.str().c_str(),"",-100.,0.);
    _powval[iOrder] ->setVal( -(double)(iOrder+1) );
    _powval[iOrder] -> removeRange();
    _powval_start[iOrder] = _powval[iOrder]->getVal();
    _pow_list -> add(*_powval[iOrder]);    

  }
  
  _pow_list -> add(*_theMass);

  pSS_powName.str("");
  pSS_powName << "fstpow_theNpow_" << base.Data() << "_" << _enum;
  _thePdf  = new RooGenericPdf(pSS_powName.str().c_str(),"",powFunc.str().c_str(),*_pow_list);

  //std::cout<<powFunc.str().c_str()<<std::endl;
  
  return;
}

void fstPowModel::resetValues() {

  for(int iOrder = 0; iOrder < _order; ++iOrder) {
    _powval[iOrder]->setVal(_powval_start[iOrder]);
    if( iOrder < _order - 1) {
      _powratio[iOrder]->setVal(_powratio_start[iOrder]);
    }
  }

  return;
}

void fstPowModel::storeValues() {
  for(int iOrder = 0; iOrder < _order; ++iOrder) {
    _powval_start[iOrder] = _powval[iOrder]->getVal();
    if( iOrder < _order - 1) {
      _powratio_start[iOrder] = _powratio[iOrder]->getVal();
    }
  }
  return;
}

RooAbsPdf* fstPowModel::getPdf() {
  return (RooAbsPdf*) _thePdf;
}


// =============== BERNSTEIN PDF ============================================================================


class fstBernModel : public fstModel {
public:
  fstBernModel() {};
  fstBernModel(RooRealVar* mass, int order, int number, TString base);
  
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

fstBernModel::fstBernModel(RooRealVar* mass, int order, int number, TString base):
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
    pSS_bernName << "fstbern_bernval_" << base.Data() << "_" << _enum << "_" << iOrder;
    _bernval[iOrder] = new RooRealVar(pSS_bernName.str().c_str(),"",-100.,100);
    _bernval[iOrder] -> removeRange();
    _bernval_start[iOrder] = _bernval[iOrder]->getVal();

    TString formula = TString(pSS_bernName.str().c_str())+TString("*")+TString(pSS_bernName.str().c_str());

    pSS_bernName.str("");
    pSS_bernName << "fstbern_bernvalsq_" << base.Data() << "_" << _enum << "_" << iOrder;    
    _bernvalsq[iOrder] = new RooFormulaVar(pSS_bernName.str().c_str(),formula.Data(),RooArgList(*_bernval[iOrder]));    
    _bernval_list -> add(*_bernvalsq[iOrder]);
  }
  
  pSS_bernName.str("");
  pSS_bernName << "fstbern_theNbern_" << base.Data() << "_" << _enum;
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


// =============== BERNSTEIN PDF, FFT DAMPED  ============================================================================


class fstBdamModel : public fstModel {
public:
  fstBdamModel() {};
  fstBdamModel(RooRealVar* mass, int order, int number, TString base);
  
  virtual void resetValues();
  virtual void storeValues();
  virtual RooAbsPdf* getPdf();
  
private:

  RooRealVar*      _theMass;

  int _enum;
  int _order;
  
  RooBernstein*    _thePdf_pol;
  RooGaussian*     _thePdf_gau;

  RooFFTConvPdf*   _thePdf;

  RooConstVar*     _bdamconstvar;
  RooRealVar**     _bdamval;
  RooFormulaVar**  _bdamvalsq;
  RooArgList*      _bdamval_list;

  RooConstVar*     _bdam_gaumean;
  RooConstVar*     _bdam_gauwidth;
  
  double*          _bdamval_start;
};

fstBdamModel::fstBdamModel(RooRealVar* mass, int order, int number, TString base):
  _theMass(mass), _enum(number), _order(order) {
  
  // second order polinomial
  _bdamval       = new RooRealVar    *[_order];
  _bdamvalsq     = new RooFormulaVar *[_order];
  _bdamval_start = new double         [_order];

  stringstream pSS_bdamName;  
  pSS_bdamName.str("");
  pSS_bdamName << "fstbdam_bdamval_list_" << base.Data() << "_" << _enum;
  _bdamval_list = new RooArgList(pSS_bdamName.str().c_str());  
  pSS_bdamName.str("");
  pSS_bdamName << "fstbdam_bdamval_constvar_" << base.Data() << "_" << _enum;
  _bdamconstvar = new RooConstVar(pSS_bdamName.str().c_str(),"",1.0);
  _bdamval_list->add(*_bdamconstvar);

  for(int iOrder = 0; iOrder < _order; ++iOrder) {
    pSS_bdamName.str("");
    pSS_bdamName << "fstbdam_bdamval_" << base.Data() << "_" << _enum << "_" << iOrder;
    _bdamval[iOrder] = new RooRealVar(pSS_bdamName.str().c_str(),"",-100.,100);
    _bdamval[iOrder] -> removeRange();
    _bdamval_start[iOrder] = _bdamval[iOrder]->getVal();

    TString formula = TString(pSS_bdamName.str().c_str())+TString("*")+TString(pSS_bdamName.str().c_str());

    pSS_bdamName.str("");
    pSS_bdamName << "fstbdam_bdamvalsq_" << base.Data() << "_" << _enum << "_" << iOrder;    
    _bdamvalsq[iOrder] = new RooFormulaVar(pSS_bdamName.str().c_str(),formula.Data(),RooArgList(*_bdamval[iOrder]));    
    _bdamval_list -> add(*_bdamvalsq[iOrder]);
  }
  
  pSS_bdamName.str("");
  pSS_bdamName << "fstbdam_theNbdam_pol_" << base.Data() << "_" << _enum;
  _thePdf_pol  = new RooBernstein(pSS_bdamName.str().c_str(),"",*_theMass, *_bdamval_list);


  // generate the damping gaussian
  pSS_bdamName.str("");
  pSS_bdamName << "fstbdam_bdam_gaumean_" << base.Data() << "_" << _enum;
  _bdam_gaumean = new RooConstVar(pSS_bdamName.str().c_str(),"",0.0);
  pSS_bdamName.str("");
  pSS_bdamName << "fstbdam_bdam_gauwidth_" << base.Data() << "_" << _enum;
  _bdam_gauwidth = new RooConstVar(pSS_bdamName.str().c_str(),"",20.0);
  
  pSS_bdamName.str("");
  pSS_bdamName << "fstbdam_theNbdam_gau_" << base.Data() << "_" << _enum;
  _thePdf_gau  = new RooGaussian(pSS_bdamName.str().c_str(),"",*_theMass,*_bdam_gaumean,*_bdam_gauwidth);
  
  pSS_bdamName.str("");
  pSS_bdamName << "fstbdam_theNbdam_" << base.Data() << "_" << _enum;
  
  _thePdf = new RooFFTConvPdf(pSS_bdamName.str().c_str(),"",*_theMass,*_thePdf_pol,*_thePdf_gau);

  return;
}

void fstBdamModel::resetValues() {

  for(int iOrder = 0; iOrder < _order; ++iOrder) {
    _bdamval[iOrder]->setVal(_bdamval_start[iOrder]);
  }

  return;
}

void fstBdamModel::storeValues() {
  for(int iOrder = 0; iOrder < _order; ++iOrder) {
    _bdamval_start[iOrder] = _bdamval[iOrder]->getVal();
  }
  return;
}

RooAbsPdf* fstBdamModel::getPdf() {
  return (RooAbsPdf*) _thePdf;
}

// =============== EXP(-CHEBYCHEV) PDF OF ORDER N ============================================================================

class fstNeyModel : public fstModel {
public:
  fstNeyModel() {};
  fstNeyModel(RooRealVar* mass, int order, int number, TString base);
  
  virtual void resetValues();
  virtual void storeValues();
  virtual RooAbsPdf* getPdf();
  
private:

  RooRealVar*      _theMass;

  int _enum;
  int _order;
  
  RooGenericPdf*       _thePdf;
  
  RooRealVar**     _neyval;
  RooArgList*      _ney_list;
  
  double*          _neyval_start;
};

fstNeyModel::fstNeyModel(RooRealVar* mass, int order, int number, TString base):
  _theMass(mass), _enum(number), _order(order) {
  
  TString chebyPol[6];
  chebyPol[0] = "@0";
  chebyPol[1] = "(2.*pow(@0,2)-1.)";
  chebyPol[2] = "(4.*pow(@0,3)-3.*@0)";
  chebyPol[3] = "(8.*pow(@0,4)-8.*pow(@0,2)+1.)";
  chebyPol[4] = "(16.*pow(@0,5)-20.*pow(@0,3)+5.*@0)";
  chebyPol[5] = "(32.*pow(@0,6)-48.*pow(@0,4)+18.*pow(@0,2)-1.)";

  // second order polinomial
  _neyval         = new RooRealVar    *[_order];
  _neyval_start   = new double         [_order];

  stringstream pSS_neyName;  
  pSS_neyName.str("");
  pSS_neyName << "fstney_ney_list_" << base.Data() << "_" << _enum;
  _ney_list = new RooArgList(pSS_neyName.str().c_str());
  
  _ney_list-> add(*_theMass);
  
  stringstream neyFunc;
  
  neyFunc << "TMath::Exp(-(1." ;
  
  for(int iOrder = 0; iOrder < _order; ++iOrder) {
    
    neyFunc << "+@"<<iOrder+1<<"*"<<chebyPol[iOrder];
    
    pSS_neyName.str("");
    pSS_neyName << "fstney_neyval_" << base.Data() << "_" << _enum << "_" << iOrder;
    _neyval[iOrder] = new RooRealVar(pSS_neyName.str().c_str(),"",-100.,100.);
    if ( !(iOrder % 2) )
      _neyval[iOrder] -> setVal( (double) (iOrder+1)/TMath::Power(10.,iOrder+1) );
    else
      _neyval[iOrder] -> setVal( - (double) (iOrder+1)/TMath::Power(10.,iOrder+1) );
    _neyval[iOrder] -> removeRange();
    _neyval_start[iOrder] = _neyval[iOrder]->getVal();
    _ney_list -> add(*_neyval[iOrder]);    
  }

  neyFunc << "))";

  pSS_neyName.str("");
  pSS_neyName << "fstney_theNney_" << base.Data() << "_" << _enum;
  //std::cout<<neyFunc.str().c_str()<<std::endl;
  _thePdf  = new RooGenericPdf(pSS_neyName.str().c_str(),"",neyFunc.str().c_str(),*_ney_list);
  
  return;
}

void fstNeyModel::resetValues() {

  for(int iOrder = 0; iOrder < _order; ++iOrder) {
    _neyval[iOrder]->setVal(_neyval_start[iOrder]);
  }

  return;
}

void fstNeyModel::storeValues() {
  for(int iOrder = 0; iOrder < _order; ++iOrder) {
    _neyval_start[iOrder] = _neyval[iOrder]->getVal();
  }
  return;
}

RooAbsPdf* fstNeyModel::getPdf() {
  return (RooAbsPdf*) _thePdf;
}


// =============== POW*EXP(-CHEBYCHEV) PDF OF ORDER N ============================================================================

class fstPeyModel : public fstModel {
public:
  fstPeyModel() {};
  fstPeyModel(RooRealVar* mass, int order, int number, TString base);
  
  virtual void resetValues();
  virtual void storeValues();
  virtual RooAbsPdf* getPdf();
  
private:

  RooRealVar*      _theMass;

  int _enum;
  int _order;
  
  RooGenericPdf*       _thePdf;
  
  RooRealVar**     _peyval;
  RooArgList*      _pey_list;
  
  double*          _peyval_start;
};

fstPeyModel::fstPeyModel(RooRealVar* mass, int order, int number, TString base):
  _theMass(mass), _enum(number), _order(order) {
  
  TString chebyPol[6];
  chebyPol[0] = "@0";
  chebyPol[1] = "(2.*pow(@0,2)-1.)";
  chebyPol[2] = "(4.*pow(@0,3)-3.*@0)";
  chebyPol[3] = "(8.*pow(@0,4)-8.*pow(@0,2)+1.)";
  chebyPol[4] = "(16.*pow(@0,5)-20.*pow(@0,3)+5.*@0)";
  chebyPol[5] = "(32.*pow(@0,6)-48.*pow(@0,4)+18.*pow(@0,2)-1.)";

  // second order polinomial
  _peyval         = new RooRealVar    *[_order+1];
  _peyval_start   = new double         [_order+1];

  stringstream pSS_peyName;  
  pSS_peyName.str("");
  pSS_peyName << "fstpey_pey_list_" << base.Data() << "_" << _enum;
  _pey_list = new RooArgList(pSS_peyName.str().c_str());
  
  _pey_list-> add(*_theMass);
  
  stringstream peyFunc;
  
  peyFunc << "pow(@0,@1)*TMath::Exp(-(1." ;
  
  pSS_peyName.str("");
  pSS_peyName << "fstpey_peypowval_" << base.Data() << "_" << _enum;
  _peyval[0] = new RooRealVar(pSS_peyName.str().c_str(),"",-100.,0.);
  _peyval_start[0] = _peyval[0]->getVal();
  _pey_list -> add(*_peyval[0]);
  
  for(int iOrder = 0; iOrder < _order; ++iOrder) {
    
    peyFunc << "+@"<<iOrder+2<<"*"<<chebyPol[iOrder];
    
    pSS_peyName.str("");
    pSS_peyName << "fstpey_peyval_" << base.Data() << "_" << _enum << "_" << iOrder;
    _peyval[iOrder+1] = new RooRealVar(pSS_peyName.str().c_str(),"",-100.,100.);
    _peyval[iOrder+1] -> removeRange();
    _peyval_start[iOrder+1] = _peyval[iOrder+1]->getVal();
    _pey_list -> add(*_peyval[iOrder+1]);    
  }

  peyFunc << "))";
  
  pSS_peyName.str("");
  pSS_peyName << "fstpey_theNpey_" << base.Data() << "_" << _enum;
  //std::cout<<peyFunc.str().c_str()<<std::endl;
  _thePdf  = new RooGenericPdf(pSS_peyName.str().c_str(),"",peyFunc.str().c_str(),*_pey_list);
  
  return;
}

void fstPeyModel::resetValues() {

  for(int iOrder = 0; iOrder <= _order; ++iOrder) {
    _peyval[iOrder]->setVal(_peyval_start[iOrder]);
  }

  return;
}

void fstPeyModel::storeValues() {
  for(int iOrder = 0; iOrder <= _order; ++iOrder) {
    _peyval_start[iOrder] = _peyval[iOrder]->getVal();
  }
  return;
}

RooAbsPdf* fstPeyModel::getPdf() {
  return (RooAbsPdf*) _thePdf;
}


// ============ GREGs PDFs

class fstGregModel : public fstModel {
public:
  fstGregModel() {};
  fstGregModel(RooRealVar* mass, int order, int number, TString base);
  
  virtual void resetValues();
  virtual void storeValues();
  virtual RooAbsPdf* getPdf();
  
private:

  RooRealVar*      _theMass;

  int _enum;
  int _order;
  
  RooGenericPdf*       _thePdf;
  
  RooRealVar**     _gregval;
  RooArgList*      _greg_list;
  
  double*          _gregval_start;
};

fstGregModel::fstGregModel(RooRealVar* mass, int order, int number, TString base):
  _theMass(mass), _enum(number) {
  
  int rOrder = 3;

  switch(order) {
  case 1:
    _gregval         = new RooRealVar    *[3];
    _gregval_start   = new double         [3];
    break;
  case 2:
    _gregval         = new RooRealVar    *[3];
    _gregval_start   = new double         [3];
    break;
  default :
    _gregval         = new RooRealVar    *[2];
    _gregval_start   = new double         [2];
    rOrder = 2;
    break;
  }
  
  stringstream pSS_gregName;  
  pSS_gregName.str("");
  pSS_gregName << "fstgreg_greg_list_" << base.Data() << "_" << _enum;
  _greg_list = new RooArgList(pSS_gregName.str().c_str());
  
  _greg_list-> add(*_theMass);
  
  stringstream gregFunc;
  
  switch(order) {
  case 1:
    gregFunc << "(pow((1+@0),@1))/pow(@0,(@2+@3*TMath::Log(@0)))";
    break;
  case 2:
    gregFunc << "1./pow((@1+@2*pow(@0,2)),@3)";
    break;
  default:
    gregFunc << "1./pow((@1+@0),@2)";
    break;
  }

    
  for(int iOrder = 0; iOrder <rOrder; ++iOrder) {
    
    pSS_gregName.str("");
    pSS_gregName << "fstgreg_gregval_" << base.Data() << "_" << _enum << "_" << iOrder;
    _gregval[iOrder] = new RooRealVar(pSS_gregName.str().c_str(),"",-100.,100.);
    _gregval[iOrder] -> removeRange();
    _gregval_start[iOrder] = _gregval[iOrder]->getVal();
    _greg_list -> add(*_gregval[iOrder]);    
  }
  
  _order = rOrder;

  pSS_gregName.str("");
  pSS_gregName << "fstgreg_theNgreg_" << base.Data() << "_" << _enum;
  //std::cout<<gregFunc.str().c_str()<<std::endl;
  _thePdf  = new RooGenericPdf(pSS_gregName.str().c_str(),"",gregFunc.str().c_str(),*_greg_list);
  
  return;
}

void fstGregModel::resetValues() {

  for(int iOrder = 0; iOrder < _order; ++iOrder) {
    _gregval[iOrder]->setVal(_gregval_start[iOrder]);
  }

  return;
}

void fstGregModel::storeValues() {
  for(int iOrder = 0; iOrder < _order; ++iOrder) {
    _gregval_start[iOrder] = _gregval[iOrder]->getVal();
  }
  return;
}

RooAbsPdf* fstGregModel::getPdf() {
  return (RooAbsPdf*) _thePdf;
}


#endif
