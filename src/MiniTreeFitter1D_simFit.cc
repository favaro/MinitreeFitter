#include "interface/MiniTreeFitter1D.hh"

#include <TString.h>
#include <TTree.h>
#include <TROOT.h>
#include <TFile.h>
#include <TCanvas.h>

#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooSimultaneous.h"
#include "RooFitResult.h"
#include "RooAddPdf.h"
#include "RooGaussian.h"
#include "RooPlot.h"
using namespace std;
using namespace RooFit;

vector<RooRealVar*> _fsigFit[2];

void MiniTreeFitter1D::scaleLumiExp( float scaleSig, float scaleBkg ) {
  _scaleSig = scaleSig;
  if( scaleBkg < 0 ) _scaleBkg = scaleSig; /// usually scale both sig and bkg
  else               _scaleBkg = scaleBkg; /// but can also increase sig and bkg differently
}


void MiniTreeFitter1D::resetBackgroundPdf( void ) {  
  int sHyp = 0;
  if( _bkgStudyPdf.size() == 0 ) {
    cout << " ---- creating signal + background pdfs function " << endl;
    for( unsigned is = 0; is <  _nSigExp.size(); is++ ) delete _nSigExp[is];
    _nSigExp.erase( _nSigExp.begin(),  _nSigExp.end() );

    vector<RooAbsPdf*> bkgPoly;
    vector<RooAbsPdf*> bkgExpo;
    vector<RooAbsPdf*> bkgPower;
    
    for( unsigned c = 0; c < _categories.size(); c++ ) {
      TString sigNameCat = _pdfSigName[sHyp] + TString::Format("_cat%d",c);
      TString bkgExtNameCat    = _pdfBkgName + TString::Format("_cat%d",c);
      TString bkgExtNameCatPol = _pdfBkgName + TString::Format("_cat%d",c) + "Poly";
      TString bkgExtNameCatExp = _pdfBkgName + TString::Format("_cat%d",c) + "Expo";
   
      RooAbsPdf *pdfSig      = _hlf->GetWs()->pdf( sigNameCat );      
      RooAbsPdf *pdfBkg      = _hlf->GetWs()->pdf( bkgExtNameCat    );
      RooAbsPdf *pdfBkgPol   = _hlf->GetWs()->pdf( bkgExtNameCatPol );
      RooAbsPdf *pdfBkgExp   = _hlf->GetWs()->pdf( bkgExtNameCatExp );
      
      RooRealVar *nBkgCat    = _hlf->GetWs()->var("nBkg" + bkgExtNameCat    );
      RooRealVar *nBkgPolCat = _hlf->GetWs()->var("nBkg" + bkgExtNameCatPol );
      RooRealVar *nBkgExpCat = _hlf->GetWs()->var("nBkg" + bkgExtNameCatExp );
      RooRealVar *nSigCat    = new RooRealVar("nSig_" + sigNameCat,"Nsig",10,-500,5000);
      
      TString allExtNameCat  = TString::Format("CMS_hgg_allExt_%s",sigNameCat.Data());
      RooAddPdf *pdfAll      = new RooAddPdf(allExtNameCat,allExtNameCat,RooArgList(*pdfSig,*pdfBkg)   ,RooArgList(*nSigCat,*nBkgCat));
      RooAddPdf *pdfAllPol   = new RooAddPdf(allExtNameCat + "Pol",allExtNameCat,RooArgList(*pdfSig,*pdfBkgPol),RooArgList(*nSigCat,*nBkgPolCat));
      RooAddPdf *pdfAllExp   = new RooAddPdf(allExtNameCat + "Exp",allExtNameCat,RooArgList(*pdfSig,*pdfBkgExp),RooArgList(*nSigCat,*nBkgExpCat));
      bkgPower.push_back( pdfAll );    
      bkgPoly .push_back( pdfAllPol ); 
      bkgExpo .push_back( pdfAllExp );    
      
      _nSigExp.push_back( nSigCat );    
    }
    _bkgStudyPdf.push_back( bkgPoly  );  
    _bkgStudyPdf.push_back( bkgPower ); 
    _bkgStudyPdf.push_back( bkgExpo  ); 
  }

  cout << " ----- reset background power law --------- " << endl;
  if( _slopeBkg.size() == _categories.size() )
  for( unsigned c = 0 ; c < _categories.size(); c++ ) {
    TString bkgExtNameCat = _pdfBkgName + TString::Format("_cat%d",c);      
    RooRealVar *slopeBkg = _hlf->GetWs()->var("slope"  + bkgExtNameCat );
    if(slopeBkg) slopeBkg->setVal( _slopeBkg[c]);
    RooRealVar *nBkgCat  = _hlf->GetWs()->var("nBkg"  + bkgExtNameCat );
    if(nBkgCat) nBkgCat->setVal(  _nBkg[c] );
  }

  cout << " ----- reset background exponential --------- " << endl;
  if( _texpBkg.size() == _categories.size() )
  for( unsigned c = 0 ; c < _categories.size(); c++ ) {
    TString bkgExtNameCat = _pdfBkgName + TString::Format("_cat%d",c) + "Expo";      
    RooRealVar *texpBkg = _hlf->GetWs()->var("texp"  + bkgExtNameCat );
    if(texpBkg) texpBkg->setVal( _texpBkg[c]);
    RooRealVar *nBkgCat  = _hlf->GetWs()->var("nBkg"  + bkgExtNameCat );
    if(nBkgCat) nBkgCat->setVal(  _nBkgExp[c] );
  }

  cout << " ----- reset background polynomial --------- " << endl;
  cout << "   -->  # of polynomials: " << _polBkg.size()  << endl;
  if( _polBkg.size() == _categories.size() )
  for( unsigned c = 0 ; c < _categories.size(); c++ ) {
    TString bkgExtNameCat = _pdfBkgName + TString::Format("_cat%d",c) + "Poly";    
    cout << "  -- set polynomial to: " << endl;
    for( unsigned ip = 0 ; ip < _polBkg[c].size(); ip++ ) {
      RooRealVar *polN = _hlf->GetWs()->var( TString::Format("hgg_p%d_",ip) + bkgExtNameCat );
      if( polN ) polN->setVal( _polBkg[c][ip] );
      cout << "    - pol" << ip << "[ cat: " << c << "] = " << _polBkg[c][ip] << endl;
    }
    RooRealVar *nBkgCat = _hlf->GetWs()->var("nBkg"  + bkgExtNameCat );
    if( nBkgCat ) nBkgCat->setVal( _nBkgPol[c] );
  }

  cout << " ----- reset signal number --------- " << endl;
  for( unsigned c = 0; c < _categories.size(); c++ ) {
    _nSigExp[c]->setVal(_scaleSig*_datasetSignalWiXS[sHyp][c]->sumEntries());
  }
  
  cout << " =====> all backgrounds reset to data fitted value <====== " << endl;

}



double MiniTreeFitter1D::simFitPreparePdf(int sHyp) {
  float nSigTot = 0;

  for( unsigned c = 0 ; c < _categories.size(); c++ ) nSigTot += _datasetSignalWiXS[sHyp][c]->sumEntries()*_scaleSig;

  if( sHyp < int(_simPdf.size()) ) {
    // /// pdf exists, just reset
    // for( unsigned c = 0 ; c < _categories.size(); c++ ) {
    //   TString bkgExtNameCat = _pdfBkgName + TString::Format("_cat%d",c);      
    //   RooRealVar *slopeBkg = _hlf->GetWs()->var("slope"  + bkgExtNameCat );
    //   slopeBkg->setVal( _slopeBkg[c]);
    //   RooRealVar *nBkgCat  = _hlf->GetWs()->var("nBkg"  + bkgExtNameCat );
    //   nBkgCat->setVal( _dataset[c]->sumEntries()*_scaleBkg );
    // }

    cout << " ----- reset background power law --------- " << endl;
    for( unsigned c = 0 ; c < _categories.size(); c++ ) {
      TString bkgExtNameCat = _pdfBkgName + TString::Format("_cat%d",c);      
      RooRealVar *slopeBkg = _hlf->GetWs()->var("slope"  + bkgExtNameCat );
      if(slopeBkg) slopeBkg->setVal( _slopeBkg[c]);
      RooRealVar *nBkgCat  = _hlf->GetWs()->var("nBkg"  + bkgExtNameCat );
       if(nBkgCat) nBkgCat->setVal(  _nBkg[c]*_scaleBkg  );
    }

    cout << " ----- reset background exponential --------- " << endl;
    for( unsigned c = 0 ; c < _categories.size(); c++ ) {
      TString bkgExtNameCat = _pdfBkgName + TString::Format("_cat%d",c) + "Expo";      
      RooRealVar *texpBkg = _hlf->GetWs()->var("texp"  + bkgExtNameCat );
      if(texpBkg) texpBkg->setVal( _texpBkg[c]);
      RooRealVar *nBkgCat  = _hlf->GetWs()->var("nBkg"  + bkgExtNameCat );
      if(nBkgCat) nBkgCat->setVal(  _nBkgExp[c]*_scaleBkg );
    }
    
    cout << " ----- reset background polynomial --------- " << endl;
    cout << "   -->  # of polynomials: " << _polBkg.size()  << endl;
    if( _polBkg.size() == _categories.size() )
    for( unsigned c = 0 ; c < _categories.size(); c++ ) {
      TString bkgExtNameCat = _pdfBkgName + TString::Format("_cat%d",c) + "Poly";    
      for( unsigned ip = 0 ; ip < _polBkg[c].size(); ip++ ) {
	RooRealVar *polN = _hlf->GetWs()->var( TString::Format("hgg_p%d_",ip) + bkgExtNameCat );
	if( polN ) polN->setVal( _polBkg[c][ip] );
	//	cout << "    - pol" << ip << "[ cat: " << c << "] = " << _polBkg[c][ip] << endl;
      }
      RooRealVar *nBkgCat = _hlf->GetWs()->var("nBkg"  + bkgExtNameCat );
      if( nBkgCat ) nBkgCat->setVal( _nBkgPol[c]*_scaleBkg );
    }
    

    for( unsigned c = 0 ; c < _categories.size(); c++ ) {
      // TString sigNameCat = _pdfSigName[sHyp] + TString::Format("_cat%d",c);
      // RooRealVar *fsigFit = _hlf->GetWs()->var( "fsigFit_"  + sigNameCat );
      //      if( fsigFit ) {
      //      cout << " fsigFit[" << c << "] = " << _fsigFit[sHyp][c]->getVal() << "  -- reset to 1" << endl;
      _fsigFit[sHyp][c]->setVal(1);
      //      } 
    }
    _nSigExp[sHyp]->setVal(nSigTot);
    return nSigTot;
  }    

  if( _category == 0 ) {
    // Define category to distinguish physics and control samples events
    _category = new RooCategory("category","catgeory") ;
    for( unsigned c = 0 ; c < _categories.size(); c++ )
      _category->defineType(TString::Format("cat%d",c));
  }

  if( _simPdf.size() == 0 && _nSigExp.size() != 0 ) _nSigExp.erase(_nSigExp.begin(),_nSigExp.end());
  
  // build combine pdf 
  _simPdf .push_back( new RooSimultaneous(("CMS_hgg_simPdf_"+_pdfSigName[sHyp]).c_str(),"simultaneous pdf",*_category) );
  _nSigExp.push_back( new RooRealVar(("CMS_nSigExp_" + _pdfSigName[sHyp]).c_str(),"# expected signal",nSigTot,-10*nSigTot,100*nSigTot) );

  for( unsigned c = 0 ; c < _categories.size(); c++ ) {
    TString sigNameCat = _pdfSigName[sHyp] + TString::Format("_cat%d",c);
    float fsigErr = 0.05;

    RooRealVar  *fsigMean  = new RooRealVar(  "fsigMean_" + sigNameCat,"fsig mean",1);
    RooRealVar  *fsigSigma = new RooRealVar(  "fsigErr_"  + sigNameCat,"fsig sigma",fsigErr);
    RooRealVar  *fsigFit   = new RooRealVar(  "fsigFit_"  + sigNameCat,"fsig sigma",1,0,5);
    RooGaussian *pdfFSIG   = new RooGaussian( "pdfFSIG_"  + sigNameCat,"fsig fit",*fsigFit,*fsigMean,*fsigSigma ); 
    
    if( !_addSyst ) fsigFit->setConstant(true);

    RooAbsPdf *pdfSig = _hlf->GetWs()->pdf( sigNameCat );
    _fsigFit[sHyp].push_back(fsigFit);
    // _hlf->GetWs()->import( *fsigFit );
    // _hlf->GetWs()->import( *pdfFSIG );

    TString bkgExtNameCat = _pdfBkgName + "???";
    if     ( _bkgModel == 0 ) bkgExtNameCat = _pdfBkgName + TString::Format("_cat%d",c);
    else if( _bkgModel == 1 ) bkgExtNameCat = _pdfBkgName + TString::Format("_cat%d",c) + "Poly";
    else if( _bkgModel == 2 ) bkgExtNameCat = _pdfBkgName + TString::Format("_cat%d",c) + "Expo";

    //TString bkgExtNameCat = _pdfBkgName + TString::Format("_cat%d",c) ;
    RooAbsPdf *pdfBkg    = _hlf->GetWs()->pdf( bkgExtNameCat + "_pdf" );

    RooRealVar    *fsig    = new RooRealVar( "fsig_"    + sigNameCat,"fsig",
					    _datasetSignalWiXS[sHyp][c]->sumEntries()/nSigTot*_scaleSig );

    RooFormulaVar *nSigCat = new RooFormulaVar("CMS_hgg_nSigCat_" + sigNameCat,"fsig","@0*@1*@2",RooArgList(*fsig,*_nSigExp[sHyp],*fsigFit));
    RooRealVar    *nBkgCat = _hlf->GetWs()->var("nBkg"  + bkgExtNameCat );
    nBkgCat->setVal( _dataset[c]->sumEntries()*_scaleBkg );
    TString allExtNameCat = TString::Format("CMS_hgg_allExt_%s",sigNameCat.Data());    
    RooAddPdf  *pdfAllFSIG = new RooAddPdf(allExtNameCat+"_toto",allExtNameCat,RooArgList(*pdfSig,*pdfBkg),RooArgList(*nSigCat,*nBkgCat));
    RooProdPdf *pdfAllExt  = new RooProdPdf(allExtNameCat,allExtNameCat,*pdfAllFSIG,*pdfFSIG);
    
    _simPdf[sHyp]->addPdf( *pdfAllExt, TString::Format("cat%d",c));
  }
  // _hlf->GetWs()->import(*_nSigExp[sHyp]);
  // _hlf->GetWs()->import(*_simPdf[sHyp]);
  return _nSigExp[sHyp]->getVal();
}


void MiniTreeFitter1D::preparePseudoData(int sHyp, vector<RooDataSet*>  &pseudoData_c) {
  RooRealVar* mass = (RooRealVar*) &_variables[massVarName().c_str()];  
  mass->setMin( _fitMassMin );
  mass->setMax( _fitMassMax );

  simFitPreparePdf(sHyp);

  for( unsigned c = 0 ; c < _categories.size(); c++ ) {    
    _category->setIndex(c);
    pseudoData_c.push_back(_simPdf[sHyp]->getPdf(_category->getLabel())->generate(*mass));
  }
}

void MiniTreeFitter1D::preparePseudoDataBinned(int sHyp, vector<RooDataHist*>  &pseudoData_c) {
  RooRealVar* mass = (RooRealVar*) &_variables[massVarName().c_str()];  
  mass->setMin( _fitMassMin );
  mass->setMax( _fitMassMax );
  mass->setBins(160);
  simFitPreparePdf(sHyp);

  for( unsigned c = 0 ; c < _categories.size(); c++ ) {    
    _category->setIndex(c);
    pseudoData_c.push_back(_simPdf[sHyp]->getPdf(_category->getLabel())->generateBinned(*mass));
  }
}



RooDataSet * MiniTreeFitter1D::prepareCombinedData( vector<RooDataSet*> &data ) {
  RooRealVar* mass = (RooRealVar*) &_variables[massVarName().c_str()];  
  mass->setMin( _fitMassMin );
  mass->setMax( _fitMassMax );
  mass->setBins(160);
  RooDataSet *combData = 0;
 
  combData = new RooDataSet("combData","combined data",*mass,Index(*_category),
			     Import("cat0",*data[0]));
  for( unsigned icat = 1; icat < _categories.size(); icat++ ) {
    RooDataSet tmp(TString::Format("combData_%d",icat),"combined data",*mass,Index(*_category),
		   Import(TString::Format("cat%d",icat),*data[icat])) ;
    combData->append( tmp );
  }
  

  return combData;
}

RooDataHist * MiniTreeFitter1D::prepareCombinedData( vector<RooDataHist*> &data ) {
  RooRealVar* mass = (RooRealVar*) &_variables[massVarName().c_str()];  
  mass->setMin( _fitMassMin );
  mass->setMax( _fitMassMax );
  mass->setBins(160);
  RooDataHist *combData = 0;
 
  combData = new RooDataHist("combData","combined data",*mass,Index(*_category),
			     Import("cat0",*data[0]));
  for( unsigned icat = 1; icat < _categories.size(); icat++ ) {
    RooDataHist tmp(TString::Format("combData_%d",icat),"combined data",*mass,Index(*_category),
		   Import(TString::Format("cat%d",icat),*data[icat])) ;
    combData->add( tmp );
  }
  

  return combData;
}



void MiniTreeFitter1D::simultaneousFit(int sHyp) {
  string rootfileout = _plotDirectory + "SimFitOut_" + _pdfSigName[sHyp] + ".root";
  TFile *f = TFile::Open(rootfileout.c_str(),"recreate");
  TTree t("SimFitResult","simultaneous fit result");
  float nsig[2], e_nsig[2], nsigExp[2],nll[2];
  int isdata(-1),status[2];
  t.Branch("nll"     ,nll     ,"nll[2]/F");
  t.Branch("nSig"    ,nsig    ,"nSig[2]/F");
  t.Branch("nSigExp" ,nsigExp ,"nSigExp[2]/F");
  t.Branch("err_nSig",e_nsig  ,"err_nSig[2]/F");
  t.Branch("status"  ,status  ,"status[2]/I");
  t.Branch("data"    ,&isdata ,"data/I");

  RooRealVar* mass = (RooRealVar*) &_variables[massVarName().c_str()];  
  
  for( unsigned stest = 0 ; stest < 2; stest++ ) {
    simFitPreparePdf(stest);
    nsigExp[stest] = _nSigExp[stest]->getVal();
  }

  TCanvas* canvas = new TCanvas( "canvasSimFit2","Simulatneous Fit after fit",1000,800);
  if     ( _categories.size() <= 3 )  canvas->Divide(_categories.size(),1);
  else if( _categories.size() <= 6 )  canvas->Divide(int(TMath::Ceil(_categories.size()/2.)),2);
  else if( _categories.size() <= 9 )  canvas->Divide(int(TMath::Ceil(_categories.size()/3.)),3);
  else                                canvas->Divide(int(TMath::Ceil(_categories.size()/4.)),4);

  int color[2] = {kOrange,kAzure};
  /// itoy == 0 = real data;
  unsigned ntoy = 1000;// ntoy = 2;
  for( unsigned itoy = 0 ; itoy <= ntoy ; itoy++ ) { 

    _simPdf.erase(  _simPdf.begin() , _simPdf.end()  );
    _nSigExp.erase( _nSigExp.begin(), _nSigExp.end() );
    for( unsigned stest = 0 ; stest < 2; stest++ ) simFitPreparePdf(stest);


    cout << "  --- prepare pseudo-data toy: " << itoy << " -----" << endl;    
    RooDataSet *combData = 0;
    //RooDataHist *combData = 0;
    if( itoy == 0 ) { isdata = 1; combData = prepareCombinedData( _dataset );}
    else            { 
      isdata = 0; 
      simFitPreparePdf(sHyp);
      vector<RooDataSet*> pseudoData_c;  preparePseudoData(sHyp, pseudoData_c );
      //vector<RooDataHist*> pseudoData_c; preparePseudoDataBinned(sHyp, pseudoData_c );
      combData = prepareCombinedData( pseudoData_c ); 
      for( unsigned c= 0; c<pseudoData_c.size(); c++ ) delete pseudoData_c[c];
    }


    for( int stest = 0 ; stest < 2; stest++ ) {
      // simFitPreparePdf(stest);
      // _nSigExp[stest]->setVal(_nSigExp[stest]->getVal()*2);
      cout << " before fit: " << endl;
      RooFitResult *r =  _simPdf[stest]->fitTo(*combData,NumCPU(10),Save(),
					       Minos(false),Hesse(false),
					       PrintEvalErrors(0),PrintLevel(-1));
      cout << " after fit: " << endl;
      status[stest] = r->status();
      if( status[stest] == 0 ) {
	nll[stest]    = r->minNll();
	nsig[stest]   = _nSigExp[stest]->getVal();
	e_nsig[stest] = _nSigExp[stest]->getError();
      } else {
	nll[stest]    = 999999;
	nsig[stest]   = -9999;
	e_nsig[stest] = +9999;
      }
      delete r;
    }
   
    t.Fill();
    for( int stest = 0 ; stest < 2; stest++ )
      cout << " nsig[is = "<<stest<<"] = " << nsig[stest] << endl;
    
    if( ! gROOT->IsBatch() && (itoy-1)%1 == 0 && itoy > 0 ) {
      for( unsigned c = 0 ; c < _categories.size(); c++ ) {
	  TString bkgExtNameCat = _pdfBkgName + TString::Format("_cat%d",c);
	  canvas->cd(c+1);
	  int nbins = int(_fitMassMax-_fitMassMin) / 1 ; /// Entries / GeV
	  RooPlot* plotMgg = mass->frame(nbins);
	  combData->plotOn(plotMgg,Cut("category==category::"+TString::Format("cat%d",c)));
	  for( int stest = 0 ; stest < 2; stest++ ) {
	    
	    _simPdf[stest]->plotOn(plotMgg,LineColor(color[stest]),Slice(*_category,TString::Format("cat%d",c)),ProjWData(*_category,*combData)) ;
	    _simPdf[stest]->plotOn(plotMgg,LineColor(color[stest]),LineStyle(kDashed),
	     			   Slice(*_category,TString::Format("cat%d",c)),ProjWData(*_category,*combData),
	     			   Components(bkgExtNameCat+"_pdf")) ;
	  }
	  plotMgg->Draw();
      }
      //      return;
      canvas->Update();
    }
    delete combData;
  }
  f->cd();
  t.Write();
  //  f->Close();
  cout << " ---- simultaneous fit results saved in: " << f->GetName() << endl;
}




//-------------------------------------------------------------------------------------------------------//
//----------------------------------------  fit single signal -------------------------------------------//
//-------------------------------------------------------------------------------------------------------//
void MiniTreeFitter1D::simultaneousFitOnlyOneSig(int sHyp) {
  string rootfileout = _plotDirectory + "SimFitOut_" + _pdfSigName[sHyp] + ".root";
  TFile *f = TFile::Open(rootfileout.c_str(),"recreate");
  TTree t("SimFitResult","simultaneous fit result");
  float nsig, e_nsig, nsigExp,nll;
  int isdata(-1);
  t.Branch("nll"     ,&nll     ,"nll/F");
  t.Branch("nSig"    ,&nsig    ,"nSig/F");
  t.Branch("nSigExp" ,&nsigExp ,"nSigExp/F");
  t.Branch("err_nSig",&e_nsig  ,"err_nSig/F");
  t.Branch("data"    ,&isdata ,"data/I");

  RooRealVar* mass = (RooRealVar*) &_variables[massVarName().c_str()];  
  // mass->setMin( _fitMassMin );
  // mass->setMax( _fitMassMax );
  // mass->setRange("SB",_fitMassMin,_fitMassMax);

  for( unsigned stest = 0 ; stest < _pdfSigName.size(); stest++ )
    simFitPreparePdf(stest);

  simFitPreparePdf(sHyp);
  nsigExp = _nSigExp[sHyp]->getVal();  

  TCanvas* canvas = new TCanvas( "canvasSimFit2","Simulatneous Fit after fit",1000,500);
  if     ( _categories.size() <= 3 )  canvas->Divide(_categories.size(),1);
  else if( _categories.size() <= 6 )  canvas->Divide(int(TMath::Ceil(_categories.size()/2.)),2);
  else if( _categories.size() <= 9 )  canvas->Divide(int(TMath::Ceil(_categories.size()/3.)),3);
  else                                canvas->Divide(int(TMath::Ceil(_categories.size()/4.)),4);

  int color = kAzure;
  /// itoy == 0 = real data;
  unsigned ntoy = 2000; ntoy = 1;
  float nSigData = -1;
  float nSigDataErr = -1;
  for( unsigned itoy = 0 ; itoy < ntoy ; itoy++ ) { 
    cout << "  --- prepare pseudo-data toy: " << itoy << " -----" << endl;    
    simFitPreparePdf(sHyp);
    RooDataSet *combData = 0;
    if( itoy == 0 ) { isdata = 1; combData = prepareCombinedData( _dataset );}
    else            { 
      isdata = 0; 
      vector<RooDataSet*> pseudoData_c;
      preparePseudoData(sHyp, pseudoData_c );
      combData = prepareCombinedData( pseudoData_c ); 
      for( unsigned c= 0; c<pseudoData_c.size(); c++ ) delete pseudoData_c[c];
    }

    _nSigExp[sHyp]->setVal(_nSigExp[sHyp]->getVal()*5); 
    RooFitResult *r =  _simPdf[sHyp]->fitTo(*combData,NumCPU(10),Save());
    nll    = r->minNll();
    nsig   = _nSigExp[sHyp]->getVal();
    e_nsig = _nSigExp[sHyp]->getError();
    cout << "  -- nsig = " << nsig << endl;
    delete r;
    t.Fill();
    if( isdata == 1 ) {
      nSigData = nsig / nsigExp;
      nSigDataErr = e_nsig /nsigExp;
    }
    if(  ! gROOT->IsBatch() && itoy%1 == 0 ) {
      for( unsigned c = 0 ; c < _categories.size(); c++ ) {
	TString bkgExtNameCat = _pdfBkgName + TString::Format("_cat%d",c);
	canvas->cd(c+1);
	int nbins = int(_fitMassMax-_fitMassMin) / 1 ; /// Entries / GeV
	RooPlot* plotMgg = mass->frame(nbins);
	combData->plotOn(plotMgg,Cut("category==category::"+TString::Format("cat%d",c)));
	_simPdf[sHyp]->plotOn(plotMgg,LineColor(color),Slice(*_category,TString::Format("cat%d",c)),ProjWData(*_category,*combData)) ;
	_simPdf[sHyp]->plotOn(plotMgg,LineColor(color),LineStyle(kDashed),
			      Slice(*_category,TString::Format("cat%d",c)),ProjWData(*_category,*combData),
			      Components(bkgExtNameCat+"_pdf")) ;	
	plotMgg->Draw();
      }
      canvas->Update();
    }
    delete combData;
  }
  cout << " ------ Data FIT RESULT -------- " << endl;
  cout << " muData = " << nSigData << " +/- " << nSigDataErr << endl;
  f->cd();
  t.Write();
  // f->Close();
  cout << " ---- simultaneous fit results saved in: " << f->GetName() << endl;
}




//-------------------------------------------------------------------------------------------------------//
//------------------------------------  Background fit study  -------------------------------------------//
//-------------------------------------------------------------------------------------------------------//
void MiniTreeFitter1D::backgroundFitResut(int icat,  unsigned iBkgModelRef ) {
  string rootfileout = _plotDirectory + "BkgFitOut" + _pdfSigName[0] + TString::Format("_cat%d",icat).Data() + string(".root");
  TFile *f = TFile::Open(rootfileout.c_str(),"recreate");
  TTree t("BkgFitResult","simultaneous fit result");
  float nll[5], mu125[5], nSigFit[5], nSigErr[5];
  float nSigExp;
  int isdata(-1), status[5];
  t.Branch("nll"     , nll     ,"nll[5]/F");
  t.Branch("mu125"   , mu125   ,"mu125[5]/F");
  t.Branch("nSigFit" , nSigFit ,"nSigFit[5]/F");
  t.Branch("nSigErr" , nSigErr ,"nSigErr[5]/F");
  t.Branch("status"  , status  ,"status[5]/I");

  t.Branch("data"    ,&isdata  ,"data/I");
  t.Branch("nSigExp" ,&nSigExp ,"nSigExp/F");

  int color[] = { kRed, kAzure, kGray+2 };
  int style[] = { 1, 2, 2 };

  RooRealVar* mass = (RooRealVar*) &_variables[massVarName().c_str()];  
  mass->setMin( _fitMassMin );
  mass->setMax( _fitMassMax );
  mass->setRange("SB" ,_fitMassMin,_fitMassMax);

  vector<string> fitRange;
  fitRange.push_back("SBpol");
  fitRange.push_back("SBpow");
  vector<float> xmin, xmax;
  xmin.push_back( _fitMassMin ); xmin.push_back(_fitMassMin); xmin.push_back(_fitMassMin);
  xmax.push_back( _fitMassMax ); xmax.push_back(_fitMassMax); xmax.push_back(_fitMassMax);

  TCanvas* canvas = new TCanvas( "canvasSimFit2","Simulatneous Fit after fit",1000,500);
  if     ( _categories.size() <= 3 )  canvas->Divide(_categories.size(),1);
  else if( _categories.size() <= 6 )  canvas->Divide(int(TMath::Ceil(_categories.size()/2.)),2);
  else if( _categories.size() <= 9 )  canvas->Divide(int(TMath::Ceil(_categories.size()/3.)),3);
  else                                canvas->Divide(int(TMath::Ceil(_categories.size()/4.)),4);


  /// itoy == 0 = real data;
  unsigned ntoy = 100;

  resetBackgroundPdf();
  for( unsigned c = 0 ; c < _categories.size(); c++ ) {
    if( icat >= 0 && int(c) != icat ) continue;
    TString bkgExtNameCat = _pdfBkgName + TString::Format("_cat%d",c);
 


    nSigExp = _scaleSig*_datasetSignalWiXS[0][c]->sumEntries();
    for( unsigned itoy = 1; itoy < ntoy ; itoy++ ) { 

      bool doPlot =  !gROOT->IsBatch() ;
      cout << "  --- prepare pseudo-data toy: " << itoy << " -----" << endl;    
    
      _bkgStudyPdf.erase( _bkgStudyPdf.begin(),  _bkgStudyPdf.end() );
      resetBackgroundPdf();
      vector<RooAbsPdf*> pdfs;
      vector<float> locmin, locmax;
      /// put the truth background model in front
      pdfs.push_back(  _bkgStudyPdf[iBkgModelRef][c] );
      locmin.push_back( xmin[iBkgModelRef] );
      locmax.push_back( xmax[iBkgModelRef] );
      
      for( unsigned ipdf = 0 ; ipdf < _bkgStudyPdf.size(); ipdf++ ) 
	if( ipdf != iBkgModelRef )  {
	  pdfs.push_back(  _bkgStudyPdf[ipdf][c] );
	  locmin.push_back( xmin[ipdf] );
	  locmax.push_back( xmax[ipdf] );
	}

      RooDataSet *pseudoData = 0;
      if( itoy == 0 )  {
	isdata = 1;
	pseudoData = _dataset[c];
      } else {
	isdata = 0;
	pseudoData = pdfs[0]->generate(*mass);
      }

      RooPlot* plotMgg = 0;
      if( doPlot ) {
	canvas->cd(c+1);
	int nbins = int(_fitMassMax-_fitMassMin) / 1 ; /// Entries / GeV
	//plotMgg = mass->frame(115,135,20);
	plotMgg = mass->frame(nbins);
	pseudoData->plotOn(plotMgg);
      }

      for( unsigned ipdf = 0 ; ipdf < pdfs.size(); ipdf++ ) {
	//	RooFitResult *r =  pdfs[ipdf]->fitTo(*pseudoData, NumCPU(10), Save(), Range(locmin[ipdf],locmax[ipdf]) );
	RooFitResult *r =  pdfs[ipdf]->fitTo(*pseudoData, NumCPU(10), Save() );
	nll[ipdf]    = r->minNll();
	status[ipdf] = r->status();
	delete r;
	
	mass->setVal(125);
	mu125[ipdf] = pdfs[ipdf]->getVal(RooArgSet(*mass));

	nSigFit[ipdf] = _nSigExp[c]->getVal();
	nSigErr[ipdf] = _nSigExp[c]->getError();
	cout << " mu125[" << ipdf << "] = " << mu125[ipdf] << endl;
	cout << "     --> nsig   = " << nSigFit[ipdf] << " +/- " << nSigErr[ipdf] << endl;
	cout << "     --> status = " << status[ipdf] << endl;
	cout << "     --> nll    = " << nll[ipdf] << endl; 

	//	cout << "      ---> nll check: " << pdfs[ipdf]->createNLL(*pseudoData, NumCPU(5), Save(), Range(locmin[ipdf],locmax[ipdf])  )->getVal() << endl;
	cout << "     range: [" << locmin[ipdf] << " ; " << locmax[ipdf] << " ] " << endl;
	
        if( ipdf == 0 && isdata ) {
	  for( unsigned ip = 0 ; ip < _polBkg[c].size(); ip++ ) {
	    RooRealVar *polN = _hlf->GetWs()->var( TString::Format("hgg_p%d_",ip) + bkgExtNameCat );
	    _polBkg[c][ip]  = polN->getVal();
	  }
	  RooRealVar *slopeBkg = _hlf->GetWs()->var("slope"  + bkgExtNameCat );
	  _slopeBkg[c] = slopeBkg->getVal();	  
	  RooRealVar *texpBkg = _hlf->GetWs()->var("texp"  + bkgExtNameCat );
	  _texpBkg[c] = texpBkg->getVal();	  
	}
      
	if( doPlot )  pdfs[ipdf]->plotOn(plotMgg,LineColor(color[ipdf]  ), LineStyle(style[ipdf]), Normalization(1), Range(_fitMassMin,_fitMassMax) );
      }
      for( unsigned ipdf = 0 ; ipdf < pdfs.size(); ipdf++ )
	cout << " Nsig["  << ipdf << "] = " << nSigFit[ipdf] <<" +/- " << nSigErr[ipdf] << " for nll = " << nll[ipdf] << endl;     
      if( doPlot ) {
	plotMgg->Draw();
	canvas->Update();
	//	getchar();
      }
      t.Fill();
      
      if( itoy != 0 ) delete pseudoData;
    }
  }
  f->cd();
  t.Write();
  // f->Close();
  cout << " ---- simultaneous fit results saved in: " << f->GetName() << endl;
}
