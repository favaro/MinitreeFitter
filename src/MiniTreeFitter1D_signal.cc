#include "interface/MiniTreeFitter1D.hh"

#include "RooFormulaVar.h"
#include "RooPlot.h"
#include "RooGaussian.h"
#include "RooAddPdf.h"
#include "RooExtendPdf.h"
#include "RooHist.h"

#include <TFile.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TMath.h>
#include <TH1.h>

#include <iostream>
using namespace std;
using namespace RooFit;


double sigmaEffective( const TH1F & );


////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////// ----------- add signal samples ----------- /////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
void MiniTreeFitter1D::addSigSamples(vector<string> samples    , vector<float> xsec, 
				     vector<string> signalNames, float br, float luminosity) {

  if( !_massVarSet ) addMassVar();

  RooRealVar lumi("lumi","lumi",luminosity);
  _hlf->GetWs()->import(lumi); 
  
  /// test signal samples and compatibility
  if( samples.size() != xsec.size() ||samples.size() != signalNames .size() ||
      signalNames .size() != xsec.size() ) { 
    cerr << "MiniTreeFitter1D[ERROR]: numbers of signal samples and correspondind xsec do not match " 
	 << "   # signal samples = " << samples.size() << endl
	 << "   # signal xsec    = " << xsec.size() << endl
	 << "   # signal names   = " << signalNames.size() << endl;
    cerr << "...... I will crash for sure ....." << endl;
    return;
  }

  float totalScale = 0;
  vector<float> scaleMC;
  for( unsigned is = 0; is < samples.size(); is++ ) {
    scaleMC.push_back( xsec[is] * luminosity * br / 100000);
    cout << " %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% " << endl;
    cout << is << " xsec = " << xsec[is] << "  scalefactor = " << scaleMC[is] << endl;
    totalScale +=  xsec[is] * luminosity * br / 100000;
  }

  

  vector<RooDataSet> datasetSignalWoXS,datasetSignalWiXS; 
  unsigned sHyp = _datasetSignalWoXS.size();

  for( unsigned isig = 0 ; isig < samples.size(); isig++ ) {
    RooRealVar *scaleName = new RooRealVar( TString::Format("scaleXS_%d_%d",sHyp,isig),"",0,1000);
    cout << " --- getting signal minitree " << _minitreeName << " from file: " << samples[isig] << endl;
    TFile file( samples[isig].c_str(),"read");
    TTree* sigTree = (TTree*) file.Get(_minitreeName.c_str());
    if( sigTree == 0 ) { cerr << " can not find minitree: " <<_minitreeName << " in file: " << samples[isig] << endl; }
    RooFormulaVar scaleXS( scaleName->GetName() , "@0*@1",
			   RooArgList( _variables["wei"], _variables["wei_xs"] ) );
    _variables.setRealValue("wei_xs", scaleMC[isig] );  
    RooDataSet tmpWiXS("tmp","tmp",sigTree,_variables,_mainCut); tmpWiXS.addColumn(scaleXS);
    _variables.setRealValue("wei_xs", scaleMC[isig] / totalScale );
    RooDataSet tmpWoXS("tmp","tmp",sigTree,_variables,_mainCut); tmpWoXS.addColumn(scaleXS);
    _variables.add( *scaleName);

    RooDataSet dWiXS( TString::Format("sig%d_ScaledWiXs_%s",sHyp, signalNames[isig].c_str()),"", 
		      _variables, Import(tmpWiXS),WeightVar(scaleXS.GetName()) );
    RooDataSet dWoXS( TString::Format("sig%d_ScaledWoXs_%s",sHyp, signalNames[isig].c_str()),"",
		      _variables, Import(tmpWoXS),WeightVar(scaleXS.GetName()) );

    datasetSignalWoXS.push_back( dWoXS );
    datasetSignalWiXS.push_back( dWiXS );
    cout << " NSig[" << sHyp << ", " << signalNames[isig] << "] = " << dWiXS.sumEntries() << endl;
  }

  /// sum up signal
  RooDataSet sigSumWiXs(datasetSignalWiXS[0]), sigSumWoXs(datasetSignalWoXS[0]);
  for( unsigned isig = 1 ; isig < samples.size(); isig++ ) {
    sigSumWiXs.append( datasetSignalWiXS[isig] );
    sigSumWoXs.append( datasetSignalWoXS[isig] );
  }

  /// split in different categories
  vector<RooDataSet*> dCatSignalWoXS, dCatSignalWiXS;
  for( unsigned c = 0 ; c < _categories.size(); c++ )
    cout << " category: " << _categories[c] << endl;
  for( unsigned c = 0 ; c < _categories.size(); c++ ) {
    dCatSignalWiXS.push_back( (RooDataSet*) sigSumWiXs.reduce(RooArgList(_variables[massVarName().c_str()]), _mainCut && _categories[c] ) );
    dCatSignalWoXS.push_back( (RooDataSet*) sigSumWoXs.reduce(RooArgList(_variables[massVarName().c_str()]), _mainCut && _categories[c] ) );
    _hlf->GetWs()->import(*dCatSignalWiXS[c],Rename( TString::Format("sig%d_ScaledWiXs_cat%d",sHyp,c) ));
    _hlf->GetWs()->import(*dCatSignalWoXS[c],Rename( TString::Format("sig%d_ScaledWoXs_cat%d",sHyp,c) ));    
  }

  /// can handle several signal hypothesis
  _datasetSignalWoXS.push_back( dCatSignalWoXS );
  _datasetSignalWiXS.push_back( dCatSignalWiXS );

  for( unsigned isig = 0 ; isig < samples.size(); isig++ ) {
    cout << " ----------- Signal[" << signalNames[isig] << "] contribution --------------------" << endl;
    for (unsigned c = 0; c < _categories.size(); ++c) 
      cout << "    - nCAT[" << c << "] = "<< datasetSignalWiXS[isig].reduce( RooArgList(_variables[massVarName().c_str()]),
									     _mainCut && _categories[c])->sumEntries() << endl;   
  }
  cout << "****************************************************************" << endl
       << "----------------------- eff x Acc numbers ----------------------" << endl
       << "****************************************************************" << endl;
  for( unsigned isig = 0 ; isig < samples.size(); isig++ ) {
    cout << " ----------- Signal[" << signalNames[isig] << "] contribution --------------------" << endl;
    for (unsigned c = 0; c < _categories.size(); ++c) 
      cout << "    - effxAcc[" << c << "] = "<< datasetSignalWiXS[isig].reduce( RooArgList(_variables[massVarName().c_str()]),
										_mainCut && _categories[c])->sumEntries()
	/( xsec[isig] * luminosity * br)*100 << endl;   
  }

  cout << "****************************************************************" << endl
       << "------------------- eff x Acc numbers overall ------------------" << endl
       << "****************************************************************" << endl;

  for( unsigned icat = 0 ; icat < _categories.size(); icat++ ) {
    double Nexp = 0;
    double Xsec = 0;
    for( unsigned isig = 0 ; isig < samples.size(); isig++ ) {
      Nexp += datasetSignalWiXS[isig].reduce( RooArgList(_variables[massVarName().c_str()]),
					       _mainCut && _categories[icat])->sumEntries();
      Xsec += xsec[isig] * luminosity * br;
    }
    if( icat%4 == 0 ) cout << "\\hline" << endl; 
    printf( "$%1.1f\\leq|\\cos\\theta^*|<%1.1f$, cat%d & $%2.2f \\%%$ & $%2.1f$ \\\\ \n", icat/4*0.2, (icat/4+1)*0.2, icat%4, Nexp / Xsec * 100,  Nexp );
  }


  /// for each signal specify that the fit is not yet done
  _sigFitDone.push_back(false);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////// --------- fitting signal samples --------- /////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
void MiniTreeFitter1D::modelSignal( float testMass, string suffixName, int sHyp ) {
  if( sHyp < 0 ) sHyp = _datasetSignalWoXS.size() -1;
  if( sHyp < 0 ) {
    cerr << "MiniTreeFitter1D::fitSigSamples[ERROR]: I do not have any signal sample to fit. bailing out..." << endl;
    return;
  }
  
  if( suffixName == "" ) suffixName = TString::Format("SignalHypo%d",sHyp).Data();
  _pdfSigName.push_back(string(TString::Format("MggSig%s_m%3.1f",suffixName.c_str(),testMass)));
  
  TString filePlotXcheck =  (_plotDirectory + "signalPlotXchecks_" + suffixName).c_str() + TString::Format("_m%.1f.eps",testMass);
  if( testMass == 0 ) testMass = 125;
  
  float minMassFit(testMass-15),maxMassFit(testMass+15);
 
  RooRealVar* mass = (RooRealVar*) &_variables[massVarName().c_str()];
  mass->SetTitle("M_{#gamma#gamma}");
  mass->setUnit("GeV");
  mass->setMin( _fitMassMin );
  mass->setMax( _fitMassMax );

  TCanvas* canvas = new TCanvas("canvasSig_AllCat","Signal Inclusive",1000,500);
  if     ( _categories.size() <= 3 )  canvas->Divide(_categories.size(),1);
  else if( _categories.size() <= 6 )  canvas->Divide(int(TMath::Ceil(_categories.size()/2.)),2);
  else if( _categories.size() <= 9 )  canvas->Divide(int(TMath::Ceil(_categories.size()/3.)),3);
  else  canvas->Divide(int(TMath::Ceil(_categories.size()/4.)),4);

  float resolutionError[] = { 0.21, 0.21, 0.40, 0.29 };     // ?
  vector<double> sigEffCat;
  for (unsigned c = 0; c < _categories.size(); ++c) {
    TString sigNameCat = _pdfSigName[sHyp] + TString::Format("_cat%d",c);

    /// take signa dataset without xs weight to get meaningful error bars.
    RooDataSet *data_c = _datasetSignalWoXS[sHyp][c];      // ?

    unsigned ncx = c+1;
    canvas->cd(ncx); 
    vector<RooRealVar*>  varLocalFit;
    varLocalFit.push_back( new RooRealVar("mu" + sigNameCat,"#mu_{core} [GeV]",testMass,100,170));
    varLocalFit.push_back( new RooRealVar("dmu"+ sigNameCat,"#mu_{core} - #mu_{broad} [GeV]",-0.6,-5,5));
    varLocalFit.push_back( new RooRealVar("sigmaCore" + sigNameCat,"sigma_{core} [GeV]",1,0.7,4));
    varLocalFit.push_back( new RooRealVar("sigmaRatio" + sigNameCat,"#sigma_{broad}/#sigma_{core}",1.5,1,6));
    varLocalFit.push_back( new RooRealVar("fCore" + sigNameCat,"f_{core}",0.5,0.1,0.99));
    
    minMassFit = TMath::Max(float(testMass-15),float(mass->getMin()));
    maxMassFit = TMath::Min(float(testMass+15),float(mass->getMax()));

    RooPlot *plotCat = mass->frame(testMass-15,testMass+15,60);    
    data_c->plotOn( plotCat );
    varLocalFit[0]->setVal(testMass); 
    varLocalFit[1]->setVal(-0.6); 
    varLocalFit[2]->setVal(0.01*testMass);
    varLocalFit[3]->setVal(1.5); 
    varLocalFit[4]->setVal(0.5); 
    

    RooFormulaVar muWide(   "muWide"    + sigNameCat,"@0+@1",RooArgList(*varLocalFit[0],*varLocalFit[1]));    // sum of means
    RooFormulaVar sigmaWide("sigmaWide" + sigNameCat,"@0*@1",RooArgList(*varLocalFit[3],*varLocalFit[2]));    // product of widths 
    RooGaussian g1("gaus1" + sigNameCat,"",*mass,*varLocalFit[0],*varLocalFit[2]);
    RooGaussian g2("gaus2" + sigNameCat,"",*mass,muWide,sigmaWide);

    RooAddPdf   MggSig( sigNameCat,"",RooArgList(g1,g2),RooArgList(*varLocalFit[4]));     // final signal pdf

    RooRealVar *nSig = new RooRealVar("nSig"  + sigNameCat,"nSig"    ,_datasetSignalWiXS[sHyp][c]->sumEntries());
    RooExtendPdf MggSigExt( sigNameCat + "Ext","",MggSig,*nSig);           // extend the fit

    MggSig.fitTo(*data_c, SumW2Error(kTRUE), Range( minMassFit, maxMassFit) );      // FIT!
    MggSig.plotOn(plotCat, LineColor(kBlue)  ); 
    MggSig.plotOn(plotCat, Components("gaus2" + sigNameCat),LineColor(kOrange+1),LineStyle(kDashed) ); 
    plotCat->SetMinimum(0);
    plotCat->SetMaximum(plotCat->GetMaximum()*1.20);
    plotCat->Draw(); 
    
    // TH1F *hist = new TH1F("hist", "hist", 580, _fitMassMin, _fitMassMax );
    // data_c->Draw("mass>>hist");
    TH1F *hist = (TH1F*)data_c->createHistogram(massVarName().c_str(),580);
    double sigmaEff = sigmaEffective(*hist);
    sigEffCat.push_back( sigmaEff );

    TLatex tex; tex.SetNDC();    
    tex.SetTextSize(0.08);
    //    tex.DrawLatex(0.20,0.8,TString::Format("#sigma_{eff} = %2.2f GeV",sigmaEff));
    tex.DrawLatex(0.16,0.825,"CMS preliminary: #sqrt{s} = 8 TeV, L = 19.5 fb^{-1}");
    tex.DrawLatex(0.16,0.70, TString::Format("cat%d: #sigma_{eff} = %2.2f GeV",c,sigmaEff));

    canvas->Update();
    
    /// *** export signal pdf
    RooArgSet vForWS;
    for( unsigned iv = 0 ; iv < varLocalFit.size(); iv++ ) {
      //  if( iv == 0 )  varLocalFit[iv]->setVal(125.5);
      //  if( iv == 2 )  varLocalFit[iv]->setVal(varLocalFit[iv]->getVal()*0.8);
      varLocalFit[iv]->setConstant(true);
      _hlf->GetWs()->import( *varLocalFit[iv] );
      vForWS.add( *varLocalFit[iv] );
    }
    _hlf->GetWs()->defineSet("SigParam" + sigNameCat, vForWS );
    _hlf->GetWs()->import(MggSigExt);  
  }
  _signalSigmaEff.push_back( sigEffCat );

  canvas->Print(filePlotXcheck);

  _sigFitDone[sHyp] = true;

}



////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////// --------- prepare signal workspace --------- ////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
void MiniTreeFitter1D::makeSignalWorkspace(void) {
  TString  muName[] = { "CMS_hgg_MuSig_spin0", "CMS_hgg_MuSig_spin2", "CMS_unknown1", "CMS_unknown2" };

  ///---- create a new workspace containing the signal pdfs and the shape syst.
  for( unsigned sHyp = 0; sHyp < _pdfSigName.size(); sHyp++ ) {
    RooWorkspace *wAll = new RooWorkspace("w_all","w_all");
    

    for( unsigned c = 0; c < _categories.size(); ++c) {
      
      TString sigNameCat = _pdfSigName[sHyp] + TString::Format("_cat%d",c);
      // (1) import signal P.D.F.s
      //      wAll->import( *_hlf->GetWs()->pdf(sigNameCat + "Ext") );
      wAll->import( *_hlf->GetWs()->pdf(sigNameCat) );
      
      // (2) Systematics on energy scale and resolution
      TString systName_EnScale = "CMS_hgg_EnScale" + sigNameCat   + "[1]";
      TString systName_EnResol = "CMS_hgg_EnResol" + sigNameCat   + "[1]";
      TString systName_MuSig   = muName[sHyp] + "[1,-100,100.0]";


      // (3) Systematics on resolution: create new sigmasa      
      wAll->factory("prod::CMS_hgg_mu"        + sigNameCat + "(mu"        + sigNameCat + "," + systName_EnScale + ")");
      wAll->factory("prod::CMS_hgg_sigmaCore" + sigNameCat + "(sigmaCore" + sigNameCat + "," + systName_EnResol + ")");
      //      wAll->factory("prod::CMS_hgg_nSig"      + sigNameCat + "(nSig"      + sigNameCat + "," + systName_MuSig   + ")");

      // (4) do reparametrization of signal
      wAll->factory( "EDIT::CMS_hgg_"  + sigNameCat + "(" + sigNameCat + 
		     //		     TString("Ext") + "," +"nSig"      + sigNameCat + "= CMS_hgg_nSig"      + sigNameCat + 
		     "," +
		     "mu"        + sigNameCat + "= CMS_hgg_mu"        + sigNameCat + "," +
		     "sigmaCore" + sigNameCat + "= CMS_hgg_sigmaCore" + sigNameCat + ")" );      
    }  
    //---- write the workspace
    TString filename( (_plotDirectory + _pdfSigName[sHyp] +".inputsig.root").c_str());
    wAll->writeToFile(filename);
    wAll->Print();
    delete wAll;
    cout << "Write signal workspace in: " << filename << " file" << endl;
  }
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// --------- get sigma effective --------- ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
double sigmaEffective( const TH1F & hist) {
  TAxis *xaxis = hist.GetXaxis();
  Int_t nb = xaxis->GetNbins();
  if(nb < 10) {
    std::cout << "effsigma: Not a valid histo. nbins = " << nb << std::endl;
    return 0.;
  }

  Double_t bwid = xaxis->GetBinWidth(1);
  if(bwid == 0) {
    std::cout << "effsigma: Not a valid histo. bwid = " << bwid << std::endl;
    return 0.;
  }

  Double_t xmin = xaxis->GetXmin();
  Double_t ave = hist.GetMean();
  Double_t rms = hist.GetRMS();

  Double_t total=0.;
  for(Int_t i=0; i<nb+2; i++) {
    total+=hist.GetBinContent(i);
  }

  if(total < 10.) {     // CF changes this to get an estimate of sigmaEff with few events
    std::cout << "effsigma: Too few entries " << total << std::endl;
    return 0.;
  }
  Int_t ierr=0;
  Int_t ismin=999;

  Double_t rlim=0.683*total;
  Int_t nrms=rms/(bwid);    // Set scan size to +/- rms
  if(nrms > nb/10) nrms=nb/10; // Could be tuned...

  Double_t widmin=9999999.;
  for(Int_t iscan=-nrms;iscan<nrms+1;iscan++) { // Scan window centre
    Int_t ibm=(ave-xmin)/bwid+1+iscan;
    Double_t x=(ibm-0.5)*bwid+xmin;
    Double_t xj=x;
    Double_t xk=x;
    Int_t jbm=ibm;
    Int_t kbm=ibm;
    Double_t bin=hist.GetBinContent(ibm);
    total=bin;
    for(Int_t j=1;j<nb;j++){
      if(jbm < nb) {
        jbm++;
        xj+=bwid;
        bin=hist.GetBinContent(jbm);
        total+=bin;
        if(total > rlim) break;
      }
      else ierr=1;
      if(kbm > 0) {
        kbm--;
        xk-=bwid;
        bin=hist.GetBinContent(kbm);
        total+=bin;
        if(total > rlim) break;
      }
      else ierr=1;
    }
    Double_t dxf=(total-rlim)*bwid/bin;
    Double_t wid=(xj-xk+bwid-dxf)*0.5;
    if(wid < widmin) {
      widmin=wid;
      ismin=iscan;
    }
  }
  if(ismin == nrms || ismin == -nrms) ierr=3;
  if(ierr != 0) std::cout << "effsigma: Error of type " << ierr << std::endl;

  return widmin;
}

