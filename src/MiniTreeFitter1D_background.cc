#include "interface/MiniTreeFitter1D.hh"

#include "interface/RooPower.hh"
#include "RooFormulaVar.h"
#include "RooPlot.h"
#include "RooAddPdf.h"
#include "RooNumIntConfig.h"
#include "RooExponential.h"
#include "RooBernstein.h"
#include "RooExtendPdf.h"

#include <TFile.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TMath.h>
#include <TF1.h>
#include <TH1.h>

#include <iomanip>
#include <iostream>
#include <fstream>
using namespace std;
using namespace RooFit;



////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// ----------- add data set ----------- ////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
void MiniTreeFitter1D::addData( string dataset ) {
  if( !_massVarSet ) addMassVar();

  TFile dataFile(dataset.c_str(),"read'");  
  TTree* dataTree     = (TTree*) dataFile.Get(_minitreeName.c_str());

  RooDataSet Data("Data","dataset",dataTree,_variables,_mainCut); 
  for (unsigned c = 0; c < _categories.size(); ++c ) {
    _dataset.push_back( (RooDataSet*) Data.reduce(RooArgList(_variables[massVarName().c_str()]), _mainCut && _categories[c]) ); 
    _hlf->GetWs()->import(*_dataset[c],Rename(TString::Format("Data_cat%d",c)));
  }

  cout << "----------- Data yields --------" << endl;
  for (unsigned c = 0; c < _categories.size(); ++c ) {
    cout << " cat" << c << ": " << _dataset[c]->sumEntries() << endl;
  }
  cout << "--------------------------------" << endl;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// ----------- fit background ----------- //////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
vector<double> MiniTreeFitter1D::modelBackground( float testMass, string suffixName ) {
  _bkgModel = 0;
  TString filePlotXcheck = _plotDirectory.c_str();
  filePlotXcheck += TString::Format("backgroundPlotXchecks_%s_m%3.1f.pdf",suffixName.c_str(),testMass);

  if( testMass == 0 ) testMass = 125;
   
  TCanvas* canvas = new TCanvas( "canvasBkgPower","Background Inclusive",1000,500);
  if     ( _categories.size() <= 3 )  canvas->Divide(_categories.size(),1);
  else if( _categories.size() <= 6 )  canvas->Divide(int(TMath::Ceil(_categories.size()/2.)),2);
  else if( _categories.size() <= 9 )  canvas->Divide(int(TMath::Ceil(_categories.size()/3.)),3);
  else  canvas->Divide(int(TMath::Ceil(_categories.size()/4.)),4);

  RooRealVar* mass = (RooRealVar*) &_variables[massVarName().c_str()];
  mass->SetTitle("M_{#gamma#gamma}");
  mass->setUnit("GeV");
  mass->setMin( _fitMassMin );
  mass->setMax( _fitMassMax );
  mass->setRange("SB1",_fitMassMin,120);
  mass->setRange("SB2",130,_fitMassMax);
  mass->setRange("SB" ,_fitMassMin,_fitMassMax);

  mass->printTree(cout);
  mass->Print();
  string rangeFit  = "SB";
  string rangePlot = "SB1,SB2";
  if( _unblind ) rangePlot = "SB";
  _pdfBkgName = "MggBkg" + suffixName;
  double signiTot = 0;
  vector<TH1F*> hMassPlot;
  vector<double> signis;
  vector<double> nSig;
  vector<double> nBkg;
  for( unsigned c = 0; c < _categories.size(); ++c) {
    unsigned ncx = 1+c;
    TString pdfBkgNameCat = _pdfBkgName + TString::Format("_cat%d",c);
    RooDataSet *data_c = _dataset[c];
    TLatex tex; tex.SetNDC();
    int nbins = int(_fitMassMax-_fitMassMin) / 1.6 ; /// Entries / GeV
    RooPlot* plotMgg = mass->frame(nbins);
    hMassPlot.push_back( (TH1F*)data_c->createHistogram(massVarName().c_str(),nbins) );


    vector<RooRealVar*> varLocalFit;
    varLocalFit.push_back( new RooRealVar("nBkg"   + pdfBkgNameCat,"# bkg" ,5000,0,1e6));
    varLocalFit.push_back( new RooRealVar("slope"  + pdfBkgNameCat,"slope" , 0.,-20.,2.));
    varLocalFit[1]->setVal(-4); 
    /// fit Mgg in MVA slices
    RooPower     MggBkg(pdfBkgNameCat+"_pdf", "",*mass, *varLocalFit[1] );
    RooExtendPdf MggBkgExt(pdfBkgNameCat,"",MggBkg,*varLocalFit[0]);

    MggBkgExt.fitTo( *data_c, Range(rangeFit.c_str()) );
    //    RooNumIntConfig* intConfig = MggBkg.getIntegratorConfig() ;
    //    intConfig->printMultiline(cout,10000,true);


    _slopeBkg.push_back( varLocalFit[1]->getVal() );
    _nBkg.push_back(  varLocalFit[0]->getVal() );
    data_c->plotOn(plotMgg, CutRange(rangePlot.c_str()) );
    MggBkgExt.plotOn(plotMgg, LineColor(kBlue), Normalization(1),Range("SB")); 

    cout << " # sHyp = " << _sigFitDone.size() << endl;

    for( unsigned sHyp = 0; sHyp < _sigFitDone.size(); sHyp++ )
    if ( _sigFitDone[sHyp] ) {
      // double norm =  5.0*_datasetSignalWiXS[sHyp][c]->sumEntries();
      double norm =  _datasetSignalWiXS[sHyp][c]->sumEntries();
      if( sHyp > 0 ) norm *= _altSigAccCorr;
      RooAbsPdf *sig = (RooAbsPdf*) _hlf->GetWs()->pdf(TString::Format("%s_cat%d",_pdfSigName[sHyp].c_str(),c));
      if( sHyp == 0 ) sig->plotOn(plotMgg,Normalization(norm,RooAbsPdf::NumEvent), LineColor(kRed) );
      if( sHyp == 1 ) sig->plotOn(plotMgg,Normalization(norm,RooAbsPdf::NumEvent), LineColor(kGreen), LineStyle(kDashed) );
      double nS = _datasetSignalWiXS[sHyp][c]->sumEntries();
      double nB = varLocalFit[0]->getVal();

      double alpha = varLocalFit[1]->getVal();// alpha = 0;
      double alphaP1 = alpha+1;
      
      double mmin  = 115; double mmax = 135;
      double Abkg = alphaP1/(TMath::Power(mmax/100.,alphaP1)-TMath::Power(mmin/100.,alphaP1));
      double Atot = alphaP1/(TMath::Power(_fitMassMax/100.,alphaP1)-TMath::Power(_fitMassMin/100.,alphaP1));
      nB *= Atot / Abkg;
      cout << " NB[tot] from fit = " <<  varLocalFit[0]->getVal() << " extrapolated in signal window --> " << nB << endl;

      TF1 f(TString::Format("fIntegral_%d_%d",c,sHyp),"TMath::Power(x/100.0,-[0])*TMath::Gaus(x,[1],[2],1)",mmin,mmax);
      f.SetParameters(alpha,125,_signalSigmaEff[sHyp][c]/sqrt(2.));
      // double SoverSqrtB = nS / sqrt(nB) *
      // 	sqrt(f.Integral(mmin,mmax) / TMath::TwoPi()) /_signalSigmaEff[sHyp][c] * sqrt(100./Abkg);
      double SoverSqrtB = nS / sqrt(nB) *
	sqrt(f.Integral(mmin,mmax) / ( 2* sqrt(TMath::Pi()) * _signalSigmaEff[sHyp][c] ) ) * sqrt(100./Abkg);
      
      signiTot += SoverSqrtB*SoverSqrtB;
      signis.push_back(SoverSqrtB);
      nSig.push_back( nS );
      nBkg.push_back( nB* 3.4*_signalSigmaEff[sHyp][c] / (_fitMassMax - _fitMassMin) );
      cout << " mmin = " << mmin << " - max = " << mmax << endl;
      cout << " Integral(of fitting function)[ cat: " << c << "] = " << f.Integral(mmin,mmax)  << 
	" vs " << sqrt(TMath::TwoPi())/4.*_signalSigmaEff[sHyp][c] << endl;
      cout << " Abkg[ cat: " << c << "] = " << Abkg << endl;
      cout << " expected signal yield (from xsec)[" << c << "] for sHyp[" << sHyp << "] = " << nS << endl;
      cout << " Ns/sqrt(Nb)[ cat: " << c << "] = " << nS/sqrt(nB)  << endl;
      cout << " significance[ cat: " << c << "] = " << SoverSqrtB << endl;
    }
    canvas->cd(ncx); plotMgg->Draw();  canvas->Update();

    /// *** export background pdf
    RooArgSet vForWS;
    for( unsigned iv = 0 ; iv < varLocalFit.size(); iv++ ) {
      _hlf->GetWs()->import( *varLocalFit[iv] );
      vForWS.add( *varLocalFit[iv] );
    }
    _hlf->GetWs()->defineSet(TString::Format("BkgParam_cat%d",c), vForWS );
    _hlf->GetWs()->import(MggBkgExt);
  }  

  signiTot = sqrt(signiTot);
  cout << " significance[ TOTAL ] = " << signiTot  << endl;
  canvas->Print(filePlotXcheck);
  
  // TCanvas *canWeigthed = new TCanvas("cWeigthed","weigthed mass plot",1000,600);
  // canWeigthed->Divide(2,1);
  // TH1F *hSoSB    = (TH1F*) hMassPlot[0]->Clone(); hSoSB->Reset();    hSoSB->Sumw2();
  // TH1F *hLogSoSB = (TH1F*) hMassPlot[0]->Clone(); h LogSoSB->Reset(); hLogSoSB->Sumw2();
  // TF1 *f = new TF1("fFit","TMath::Power(x,[0])*[1]+[2]*TMath::Gaus(x,[3],[4])",_fitMassMin,_fitMassMax);
  // f->SetParameters(-4,100,1,126,2);  
  // f->SetParLimits(0,-10,1);
  // f->FixParameter(3,126.5);
  // f->FixParameter(4,1.8);
  // double sumSoSB    = 0;
  // double sumLogSoSB = 0;
  // double sumSig     = 0;
  // for( unsigned ic = 0 ; ic < hMassPlot.size(); ic++ ) {
  //   cout << " cat[ " << ic << " ]: nSig = " << nSig[ic] << "  nBkg = " << nBkg[ic] << endl;
  //   sumSig     += nSig[ic];
  //   sumSoSB    += nSig[ic] / (nSig[ic]+nBkg[ic]);
  //   sumLogSoSB += log( 1 + nSig[ic] / nBkg[ic] );
  // }
  // for( unsigned ic = 0 ; ic < hMassPlot.size(); ic++ ) {
  //   hSoSB   ->Add( hMassPlot[ic], nSig[ic] / (nSig[ic]+nBkg[ic]) / sumSoSB );
  //   hLogSoSB->Add( hMassPlot[ic], log( 1 + nSig[ic] / nBkg[ic] ) / sumLogSoSB );
  // }
  // hSoSB->SetMinimum(0);
  // hLogSoSB->SetMinimum(0);
  // canWeigthed->cd(1); hSoSB->DrawCopy("e");
  // canWeigthed->cd(2); hLogSoSB->Fit(f,"MHE R"); hLogSoSB->DrawCopy("e");
  // canWeigthed->Print(filePlotXcheck);
  return signis;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// ----------- fit background ----------- //////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
double MiniTreeFitter1D::modelBackgroundExp( float testMass, string suffixName ) {
  _bkgModel = 2;
  TString filePlotXcheck = _plotDirectory.c_str();
  filePlotXcheck += TString::Format("backgroundExpPlotXchecks_%s_m%3.1f.pdf",suffixName.c_str(),testMass);

  if( testMass == 0 ) testMass = 125;
   
  TCanvas* canvas = new TCanvas( "canvasBkgExpo","Background Inclusive",1000,500);
  if     ( _categories.size() <= 3 )  canvas->Divide(_categories.size(),1);
  else if( _categories.size() <= 6 )  canvas->Divide(int(TMath::Ceil(_categories.size()/2.)),2);
  else if( _categories.size() <= 9 )  canvas->Divide(int(TMath::Ceil(_categories.size()/3.)),3);
  else  canvas->Divide(int(TMath::Ceil(_categories.size()/4.)),4);

  RooRealVar* mass = (RooRealVar*) &_variables[massVarName().c_str()];
  mass->SetTitle("M_{#gamma#gamma}");
  mass->setUnit("GeV");
  mass->setMin( _fitMassMin );
  mass->setMax( _fitMassMax );
  mass->setRange("SB1",_fitMassMin,120);
  mass->setRange("SB2",130,_fitMassMax);
  mass->setRange("SB" ,_fitMassMin,_fitMassMax);

  mass->printTree(cout);
  mass->Print();
  string rangeFit  = "SB";
  string rangePlot = "SB1,SB2";
  if( _unblind ) rangePlot = "SB";
  _pdfBkgName = "MggBkg" + suffixName;

  for( unsigned c = 0; c < _categories.size(); ++c) {
    unsigned ncx = 1+c;
    TString pdfBkgNameCat = _pdfBkgName + TString::Format("_cat%d",c) + "Expo";
    RooDataSet *data_c = _dataset[c];
    vector<RooRealVar*> varLocalFit;
    varLocalFit.push_back( new RooRealVar("nBkg"  + pdfBkgNameCat,"# bkg", 100,0,1e6));
    varLocalFit.push_back( new RooRealVar("texp1" + pdfBkgNameCat,"texp1" ,-0.1,-20.,0.));
    varLocalFit.push_back( new RooRealVar("texp2" + pdfBkgNameCat,"texp2" ,-0.1,-20.,0.));
    varLocalFit.push_back( new RooRealVar("a1"    + pdfBkgNameCat,"a1"    ,0.5,0.,1.));


    TLatex tex; tex.SetNDC();
    int nbins = int(_fitMassMax-_fitMassMin) / 1 ; /// Entries / GeV
    RooPlot* plotMgg = mass->frame(nbins);
    varLocalFit[1]->setVal(-0.1); 
    /// fit Mgg in MVA slices
    RooExponential MggExp1( pdfBkgNameCat+"_pdf1", "",*mass, *varLocalFit[1]);
    RooExponential MggExp2( pdfBkgNameCat+"_pdf2", "",*mass, *varLocalFit[2]);
    RooAddPdf      MggBkg(  pdfBkgNameCat+"_pdf", "",MggExp1,MggExp2,*varLocalFit[3]);
    RooExtendPdf MggBkgExt(pdfBkgNameCat,"",MggBkg,*varLocalFit[0]);

    MggBkgExt.fitTo( *data_c, Range(rangeFit.c_str()) );

    _texpBkg.push_back( varLocalFit[1]->getVal() );
    _nBkgExp.push_back( varLocalFit[0]->getVal() );
    data_c->plotOn(plotMgg, CutRange(rangePlot.c_str()) );
    MggBkgExt.plotOn(plotMgg, LineColor(kBlue), Normalization(1),Range("SB")); 

    for( unsigned sHyp = 0; sHyp < _sigFitDone.size(); sHyp++ )
    if ( _sigFitDone[sHyp] ) {
      double norm =  5.0*_datasetSignalWiXS[sHyp][c]->sumEntries();
      if( sHyp > 0 ) norm *= _altSigAccCorr;
      RooAbsPdf *sig = (RooAbsPdf*) _hlf->GetWs()->pdf(TString::Format("%s_cat%d",_pdfSigName[sHyp].c_str(),c));
      if( sHyp == 0 ) sig->plotOn(plotMgg,Normalization(norm,RooAbsPdf::NumEvent), LineColor(kRed) );
      if( sHyp == 1 ) sig->plotOn(plotMgg,Normalization(norm,RooAbsPdf::NumEvent), LineColor(kGreen), LineStyle(kDashed) );
    }
    canvas->cd(ncx); plotMgg->Draw();  canvas->Update();

    /// *** export background pdf
    RooArgSet vForWS;
    for( unsigned iv = 0 ; iv < varLocalFit.size(); iv++ ) {
      _hlf->GetWs()->import( *varLocalFit[iv] );
      vForWS.add( *varLocalFit[iv] );
    }
    _hlf->GetWs()->defineSet(TString::Format("BkgParam_cat%d",c), vForWS );
    _hlf->GetWs()->import(MggBkgExt);
  }  
  canvas->Print(filePlotXcheck);
  return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////// ----------- fit background pol ----------- /////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
double MiniTreeFitter1D::modelBackgroundPol( float testMass, string suffixName ) {
  _bkgModel = 1;
  TString filePlotXcheck = _plotDirectory.c_str();
  filePlotXcheck += TString::Format("backgroundPolPlotXchecks_%s_m%3.1f.eps",suffixName.c_str(),testMass);

  if( testMass == 0 ) testMass = 125;
   
  TCanvas* canvas = new TCanvas( "canvasBkgPoly","Background Inclusive",900,500);
  if     ( _categories.size() <= 3 )  canvas->Divide(_categories.size(),1);
  else if( _categories.size() <= 6 )  canvas->Divide(int(TMath::Ceil(_categories.size()/2.)),2);
  else if( _categories.size() <= 9 )  canvas->Divide(int(TMath::Ceil(_categories.size()/3.)),3);
  else  canvas->Divide(4,int(TMath::Ceil(_categories.size()/4.)));

  RooRealVar* mass = (RooRealVar*) &_variables[massVarName().c_str()];
  mass->SetTitle("M_{#gamma#gamma}");
  mass->setUnit("GeV");
  mass->setMin( _fitMassMin );
  mass->setMax( _fitMassMax );
  mass->setRange("SB1",_fitMassMin,110);
  mass->setRange("SB2",150,_fitMassMax);
  mass->setRange("SB" ,_fitMassMin,_fitMassMax);

  mass->printTree(cout);
  mass->Print();
  string rangeFit  = "SB";
  string rangePlot = "SB1,SB2";
  if( _unblind ) rangePlot = "SB";
  _pdfBkgName = "MggBkg" + suffixName;
  double signiTot = 0;
  assert( _polynomialOrder.size() == _categories.size() );
  for( unsigned c = 0; c < _categories.size(); ++c) {
    unsigned ncx = 1+c;
    TString pdfBkgNameCat = _pdfBkgName + TString::Format("_cat%d",c) + "Poly";
    RooDataSet *data_c = _dataset[c];
    
    //   RooAbsPdf *MggBkg = _hlf->GetWs()->pdf( "MggBkg_example" );
    vector<RooRealVar*> varLocalFit;
    varLocalFit.push_back( new RooRealVar("nBkg"    + pdfBkgNameCat,"# bkg" ,data_c->sumEntries(),
					  data_c->sumEntries()/2.,data_c->sumEntries()*2) );

    varLocalFit.push_back( new RooRealVar("hgg_p0_" + pdfBkgNameCat,"p0",1.));
    varLocalFit.push_back( new RooRealVar("hgg_p1_" + pdfBkgNameCat,"p1",0.70,0,2));
    if( _polynomialOrder[c] >= 2 ) varLocalFit.push_back( new RooRealVar("hgg_p2_" + pdfBkgNameCat,"p2",0.20,0,2));
    if( _polynomialOrder[c] >= 3 ) varLocalFit.push_back( new RooRealVar("hgg_p3_" + pdfBkgNameCat,"p3",0.0,0,2));
    if( _polynomialOrder[c] >= 4 ) varLocalFit.push_back( new RooRealVar("hgg_p4_" + pdfBkgNameCat,"p4",0.0,0,2));
    if( _polynomialOrder[c] >= 5 ) varLocalFit.push_back( new RooRealVar("hgg_p5_" + pdfBkgNameCat,"p5",0.0,0,2));
    if( _polynomialOrder[c] >= 6 ) varLocalFit.push_back( new RooRealVar("hgg_p6_" + pdfBkgNameCat,"p6",0.0,0,2));
    if( _polynomialOrder[c] >= 7 ) varLocalFit.push_back( new RooRealVar("hgg_p7_" + pdfBkgNameCat,"p7",0.0,0,2));
    if( _polynomialOrder[c] >= 8 ) varLocalFit.push_back( new RooRealVar("hgg_p8_" + pdfBkgNameCat,"p8",0.0,0,2));

    TLatex tex; tex.SetNDC();
    int nbins = int(_fitMassMax-_fitMassMin) / 1 ; /// Entries / GeV
    RooPlot* plotMgg = mass->frame(nbins);
    for( unsigned ip = 2; ip < varLocalFit.size(); ip++ ) varLocalFit[0]->setVal(0); 

    /// fit Mgg in MVA slices
    RooArgList poln;
    for( unsigned ip = 1; ip < varLocalFit.size(); ip++ ) poln.add( *varLocalFit[ip] );

    RooBernstein MggBkg(    pdfBkgNameCat + "_pdf", "",*mass, poln);
    RooExtendPdf MggBkgExt( pdfBkgNameCat,"",MggBkg,*varLocalFit[0]);

    RooFitResult *fitresult = MggBkgExt.fitTo( *data_c, Range(rangeFit.c_str()), Save() );
    vector<float> polyVal;
    for( unsigned ip = 1; ip < varLocalFit.size(); ip++ ) polyVal.push_back( varLocalFit[ip]->getVal() );
    _polBkg.push_back( polyVal );
    _nBkgPol.push_back(  varLocalFit[0]->getVal() );
    data_c->plotOn(plotMgg, CutRange(rangePlot.c_str()) );
    MggBkgExt.plotOn(plotMgg, FillColor(kGreen+2) ,Range("SB"), VisualizeError(*fitresult, 2.0,kFALSE));
    MggBkgExt.plotOn(plotMgg, FillColor(kYellow+1),Range("SB"), VisualizeError(*fitresult, 1.0,kFALSE));
    MggBkgExt.plotOn(plotMgg, LineColor(kRed), LineStyle(kDashed), Normalization(1),Range("SB")); 
    data_c->plotOn(plotMgg, CutRange(rangePlot.c_str()) );

    for( unsigned sHyp = 0; sHyp < _sigFitDone.size(); sHyp++ )
    if ( _sigFitDone[sHyp] ) {
      double norm =  5.0*_datasetSignalWiXS[sHyp][c]->sumEntries();
      if( sHyp > 0 ) norm *= _altSigAccCorr;
      RooAbsPdf *sig = (RooAbsPdf*) _hlf->GetWs()->pdf(TString::Format("%s_cat%d",_pdfSigName[sHyp].c_str(),c));
      if( sHyp == 0 ) sig->plotOn(plotMgg,Normalization(norm,RooAbsPdf::NumEvent), LineColor(kBlue) );
      if( sHyp == 1 ) sig->plotOn(plotMgg,Normalization(norm,RooAbsPdf::NumEvent), LineColor(kGreen+1), LineStyle(kDashed) );
    }

    canvas->cd(ncx); 
    plotMgg->Draw();  
    tex.DrawLatex(0.4 ,0.80,"CMS preliminary");
    tex.DrawLatex(0.4, 0.65,"#sqrt{s} = 8 TeV, L = 19.5 fb^{-1}");
    tex.DrawLatex(0.6, 0.50, TString::Format("cat%d",c));
    canvas->Update();

    /// *** export background pdf
    RooArgSet vForWS;
    for( unsigned iv = 0 ; iv < varLocalFit.size(); iv++ ) {
      _hlf->GetWs()->import( *varLocalFit[iv] );
      vForWS.add( *varLocalFit[iv] );
    }
    _hlf->GetWs()->defineSet(TString::Format("BkgParamPol_cat%d",c), vForWS );
    _hlf->GetWs()->import(MggBkgExt);
  }  
  signiTot = sqrt(signiTot);
  canvas->Print(filePlotXcheck);
  return signiTot;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////// --------- prepare bkg workspace --------- ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
void MiniTreeFitter1D::makeBackgroundWorkspace( string bkg ) {
    
  ///---- create a new workspace containing the signal pdfs and the shape syst.
  RooWorkspace *wAll = new RooWorkspace("w_all","w_all");
  for( unsigned c = 0; c < _categories.size(); ++c) {
    TString bkgNameCat = _pdfBkgName + TString::Format("_cat%d",c); 

    // (0) import dataset
    wAll->import(*_dataset[c],Rename(TString::Format("data_obs_cat%d",c)));;

    if( bkg == "Pow" ) {
      // (1) import signal P.D.F.s
      wAll->import( *_hlf->GetWs()->pdf(bkgNameCat) );
      
      // // (3) Systematics on resolution: create new sigmasa     
      wAll->factory("CMS_hgg_slope" + bkgNameCat + TString::Format("[%f,-20,0]", _hlf->GetWs()->var( "slope" + bkgNameCat )->getVal())  );
      wAll->factory("CMS_hgg_nBkg"  + bkgNameCat + TString::Format("[%f,0,1e6]", _hlf->GetWs()->var( "nBkg"  + bkgNameCat )->getVal())  ); 
   
    // // (4) do reparametrization of signal
    wAll->factory( "EDIT::CMS_hgg_"  + bkgNameCat + "("      + bkgNameCat + "," +
     		   "slope"           + bkgNameCat + "= CMS_hgg_slope" + bkgNameCat + ","
     		   "nBkg"            + bkgNameCat + "= CMS_hgg_nBkg"  + bkgNameCat +
		   ")" );
    } else if (bkg == "Pol" ) {
      bkgNameCat += "Poly";
      wAll->import( *_hlf->GetWs()->pdf(bkgNameCat) );
      
      // // (2) Systematics on resolution: create new sigmasa
      wAll->factory( "CMS_hgg_p0_"+ bkgNameCat + "[1]" );
      for( int pOrder = 1;  pOrder <= _polynomialOrder[c]; pOrder++ ) 
	wAll->factory( TString::Format("CMS_hgg_p%d_",pOrder ) + bkgNameCat + 
		       TString::Format("[%f,0,+2]", _hlf->GetWs()->var( TString::Format("hgg_p%d_",pOrder) + bkgNameCat )->getVal())  );
      wAll->factory("CMS_hgg_nBkg"  + bkgNameCat + TString::Format("[%f,0,1e6]", _hlf->GetWs()->var( "nBkg"  + bkgNameCat )->getVal())  ); 
      
      // // (4) do reparametrization of signal
      TString factoBuild =  "EDIT::CMS_hgg_"  + bkgNameCat + "(" + bkgNameCat + "," ;
      for( int pOrder = 0;  pOrder <= _polynomialOrder[c]; pOrder++ ) 
	factoBuild += TString::Format("hgg_p%d_",pOrder) +  bkgNameCat + 
	  TString::Format("=CMS_hgg_p%d_",pOrder) + bkgNameCat + ",";
      factoBuild += "nBkg" + bkgNameCat + "= CMS_hgg_nBkg"  + bkgNameCat + ")";
      cout << " Factory building: " << factoBuild << endl;
      wAll->factory(factoBuild );
    }
  }  
  //---- write the workspace
  TString filename( (_plotDirectory + _pdfBkgName +".inputbkg.root").c_str());
  wAll->writeToFile(filename);
  wAll->Print();
  delete wAll;
  cout << "Write background workspace in: " << filename << " file" << endl;
}




////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// ----------- create data card ----------- /////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
void MiniTreeFitter1D::createDataCard( string cardname, bool useAltSig) {

  int iSig = 0;
  //  if( useAltSig ) iSig = 1;

  //  useAltSig = false;

  ///--------- quick summary output  
  cout << "===================== Data Events Yields =======================" << endl;  
  cout << "......... Lumi = " << _hlf->GetWs()->var("lumi")->getVal() << " /fb ............................" << endl;  
  for( unsigned c = 0; c < _categories.size(); ++c ) 
    cout << TString::Format("#Events data cat%d:   ",c) << _dataset[c]->numEntries()  << endl;
  
  cout << "================= Expected Signal Yields =======================" << endl;  
  cout << ".........Expected Signal for L = " << _hlf->GetWs()->var("lumi")->getVal() << " /fb ............................" << endl; 
  float nSigTot = 0;
  for( unsigned c = 0; c < _categories.size(); ++c ) {
    cout   << TString::Format("#Events Def. Signal cat%d: ",c) << _datasetSignalWiXS[iSig][c]->sumEntries() << endl;
    nSigTot +=  _datasetSignalWiXS[iSig][c]->sumEntries();
    // if( useAltSig ) 
    //   cout << TString::Format("#Events Alt. Signal cat%d: ",c) << _datasetSignalWiXS[1][c]->sumEntries() << endl;
  }
  cout << "   ..... total expected: " << nSigTot << endl;
  cout << "======================= +++++++++++++ ==========================" << endl;  


  
  ///--------- actually write the datacard
  string carddir = _plotDirectory;
  ofstream outFile( (carddir + cardname ).c_str() );

  outFile << "#CMS-HGG DataCard for Unbinned Limit Setting " << endl;
  outFile << "# Lumi =  " << _hlf->GetWs()->var("lumi")->getVal() << " fb-1" << endl;
  outFile << "imax " << _categories.size() << endl;
  outFile << "jmax " << 1+int(useAltSig) << endl;
  //outFile << "jmax " << 1 << endl;
  outFile << "kmax *" << endl;
  outFile << "---------------" << endl;

  string bkgName    = "HggBkg"; // _pdfBkgName
  string sigName    = "HggSig"; // _pdfSigName[0]
  string sigAltName = "HggSig_ALT"; // _pdfSigName[1]
  if( iSig == 1 ) sigName = sigAltName; 
  

  // outFile   << "shapes data_obs  * " <<_pdfBkgName    + ".inputbkg.root" << " w_all:data_obs_$CHANNEL" << endl;
  // outFile   << "shapes " << _pdfBkgName    << " * " << _pdfBkgName    + ".inputbkg.root" 
  // 	    << " w_all:CMS_hgg_"<< _pdfBkgName    <<"_$CHANNEL" << endl;
  // outFile   << "shapes " << _pdfSigName[0] << " * " << _pdfSigName[0] + ".inputsig.root" 
  // 	    << " w_all:CMS_hgg_"<< _pdfSigName[0] <<"_$CHANNEL" << endl;
  // if( useAltSig ) 
  //   outFile << "shapes " << _pdfSigName[1] << " * " <<  _pdfSigName[1] +".inputsig.root" 
  // 	    << " w_all:CMS_hgg_"<< _pdfSigName[1] <<"_$CHANNEL" << endl;

  string bkgModel = "Poly";
  if     ( _bkgModel == 0 ) bkgModel = "";
  else if( _bkgModel == 1 ) bkgModel = "Poly";
  else if( _bkgModel == 2 ) bkgModel = "Expo";

  outFile   << "shapes data_obs  * " <<_pdfBkgName    + ".inputbkg.root" << " w_all:data_obs_$CHANNEL" << endl;
  outFile   << "shapes " << bkgName    << " * " << _pdfBkgName    + ".inputbkg.root" 
	    << " w_all:CMS_hgg_"<< _pdfBkgName  <<"_$CHANNEL" << bkgModel << endl;
  outFile   << "shapes " << sigName << " * " << _pdfSigName[iSig] + ".inputsig.root" 
	    << " w_all:CMS_hgg_"<< _pdfSigName[iSig] <<"_$CHANNEL" << endl;
  if( useAltSig ) 
    outFile << "shapes " << sigAltName << " * " <<  _pdfSigName[1] +".inputsig.root" 
  	    << " w_all:CMS_hgg_"<< _pdfSigName[1] <<"_$CHANNEL" << endl;



  outFile << "---------------" << endl;
  outFile << "bin           ";
  for( unsigned c = 0; c < _categories.size(); ++c ) outFile << "cat" << c << "   ";
  outFile << endl;

  outFile << setprecision(5) << "observation   "; 
  for( unsigned c = 0; c < _categories.size(); ++c ) outFile << _dataset[c]->numEntries() << "  ";
  outFile << endl;
						       
  outFile << "------------------------------" << endl;
  outFile << "bin               ";
  for( unsigned c = 0; c < _categories.size(); ++c ) {
    outFile << "cat" << c << "        " << "cat" << c << "        " ;
    if( useAltSig ) outFile << "cat" << c << "        " ;
  }
  outFile << endl;
  outFile << "process           ";
  for( unsigned c = 0; c < _categories.size(); ++c ) {
    outFile << sigName <<  "      " ;
    outFile << bkgName <<  "      " ;
    if( useAltSig ) outFile << sigAltName << "  " ;
  }
  outFile << endl;

  outFile << "process           ";
  for( unsigned c = 0; c < _categories.size(); ++c ) {
    outFile << "  0  " << "       " ;
    outFile << "  1  " << "       " ;
    if( useAltSig ) outFile << "   -1      " << "  " ;
  }
  outFile << endl;
  outFile << setprecision(4) << "rate           ";
  for( unsigned c = 0; c < _categories.size(); ++c ) {
    outFile << _datasetSignalWiXS[iSig][c]->sumEntries() << "  ";
    outFile << _dataset[c]->numEntries() << "  ";
    if( useAltSig ) outFile << _datasetSignalWiXS[1][c]->sumEntries()*_altSigAccCorr << "  ";
  }
  outFile << endl;
  outFile << endl;
  outFile << "--------------------------------" << endl;
  outFile << "lumi            lnN  ";
  for( unsigned c = 0; c < _categories.size(); ++c ) {
    outFile << " 0.978/1.022 " ;
    outFile << "      -      ";
    if( useAltSig ) outFile <<  " 0.978/1.022 " ;
  } 
  outFile << endl;

  // outFile << "QCDscale_qqH    lnN  ";
  // for( unsigned c = 0; c < _categories.size(); ++c ) {
  //   outFile << " 0.939/1.062 " ;
  //   outFile << "     -       " ;
  //   if( useAltSig ) outFile <<  " 0.939/1.062 " ;
  // } 
  // outFile << endl;

  // outFile << "pdf_qqbar       lnN  ";
  // for( unsigned c = 0; c < _categories.size(); ++c ) {
  //   outFile << " 0.968/1.035 " ;
  //   outFile << "      -      ";
  //   if( useAltSig ) outFile <<  " 0.968/1.035 " ;
  // } 
  // outFile << endl;

  /// for SpinAnalysis only
  //  if( iSig == 0 ) outFile << "CMS_hgg_MuSig_spin0 flatParam #uncertainty on signal norm" << endl;
  //  if( iSig == 1 ) outFile << "CMS_hgg_MuSig_spin2 flatParam #uncertainty on signal norm" << endl;


  //// =========================== Do I really need the following syst ???? =========================== ////
  // outFile << "CMS_hgg_bkgNorm   lnN  ";
  // for( unsigned c = 0; c < _categories.size(); ++c ) {
  //   outFile << " - " ;
  //   outFile << "  1.05   ";
  //   if( useAltSig ) outFile <<  " - " ;
  // } 
  // outFile << endl;

  for( unsigned c = 0; c < _categories.size(); ++c) {
    TString bkgNameCat = _pdfBkgName + TString::Format("_cat%d",c); 
    if( bkgModel == "Poly" ) {
      bkgNameCat += "Poly";
      outFile << "CMS_hgg_nBkg"  + bkgNameCat  << " flatParam  # uncertainty on background norm" << endl;
      for( int pOrder = 1;  pOrder <= _polynomialOrder[c]; pOrder++ ) 
	outFile << TString::Format("CMS_hgg_p%d_",pOrder ) + bkgNameCat <<" flatParam  # polynomial coefs" << endl;
    } else {
      outFile << "CMS_hgg_nBkg"  + bkgNameCat  << " flatParam  # uncertainty on background norm" << endl;
      outFile << "CMS_hgg_slope" + bkgNameCat  << " flatParam  # uncertainty on background norm" << endl;
    }
  }

  // outFile << "CMS_hgg_eff_g         lnN  1.020      -        1.020       -         1.032    -       1.032   -   # Photon efficiency" << endl;
  // outFile << "CMS_hgg_ntrigbarrel   lnN  1.01       -        1.01        -         -        -       -       -   # Trigger efficiency in barrel" << endl;
  // outFile << "CMS_hgg_ntrigmixed    lnN  -          -        -           -         1.01     -       1.01    -   # Trigger efficiency in endcap " << endl;
  // outFile << "CMS_hgg_migr9         lnN  1.040      -        0.967       -         1.065    -       0.94   -   # Migration R9 accounts N(A) * (kappa(A) - 1) + N(B) * (kap

  outFile.close();
  cout << " - Data card written in: " << carddir + cardname << endl;
}




////////////////////////// dumping function
void MiniTreeFitter1D::dumpBkgFitParam() {
  cout << "-------------- power law fit: " << endl;
  for( unsigned ic = 0; ic < _nBkg.size(); ic++ ) {
    cout << " Nbkg[cat" << ic << "] = " << _nBkg[ic] << endl;
    cout << "    - power[cat" << ic << "] = " << _slopeBkg[ic] << endl;
  }

  cout << "-------------- exponential fit: " << endl;
  for( unsigned ic = 0; ic < _nBkgExp.size(); ic++ ) {
    cout << " Nbkg[cat" << ic << "] = " << _nBkgExp[ic] << endl;
    cout << "    - tau[cat" << ic << "]*125 = " << _texpBkg[ic]*125 << endl;
  }

  cout << "-------------- polynomial fit: " << endl;
  for( unsigned ic = 0; ic < _nBkgPol.size(); ic++ ) {
    cout << " Nbkg[cat" << ic << "] = " << _nBkgPol[ic] << endl;
    for( unsigned ip = 0; ip < _polBkg[ic].size(); ip++ ) {
      cout << "    - p" << ip << "[cat" << ic << "] = " << _polBkg[ic][ip] << endl;
    }
  }
}

