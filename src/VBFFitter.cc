#include "HiggsCrossSectionReader.hh"
#include "MiniTreeFitter1D.hh"


#include <TFile.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TRint.h>
#include <TH1.h>

#include <iomanip>
#include <boost/program_options.hpp>
using namespace std;


float AccTot(string rootfile, const vector<TCut>& cuts, string treename = "HToGG" );

int main( int nargc, char **argv ) {
    
  namespace po = boost::program_options;
  string config_file;
  string dirAFS;
  float mMin, mMax, mh;
  int  bkgModel;
  int  cat;
  bool addSig, addBkg, simFit;
  bool doFits = true;
  bool rootInteractive(true);
  po::options_description config("Configuration");
  config.add_options()
    ("help,h"   ,"help")
    ("mh"  , po::value<float>(&mh)       ->default_value(125),"Higgs Mass" )
    ("bkg" , po::value<int>(&bkgModel)   ->default_value(0)  ,"bkg model")
    ("mMin",po::value<float>(&mMin)      ->default_value(100.),"minimum mass cut")
    ("mMax",po::value<float>(&mMax)      ->default_value(180.),"maximum mass cut")
    ("inDir,d",po::value<string>(&dirAFS)->default_value("../diphoton2012_mvaSel_cms53x_v10/"),"directory with input")
    ("addSig", po::value<bool>(&addSig)  ->default_value(true ),"add signal" )
    ("addBkg", po::value<bool>(&addBkg)  ->default_value(true ),"add background" )
    ("simFit", po::value<bool>(&simFit)  ->default_value(false),"simultaneous fit" )
    ("cat"   , po::value<int>(&cat)      ->default_value(-1),"category" )
    ("doFits", po::value<bool>(&doFits)      ->default_value(1),"do fits" )
    ("interactive,i", po::value<bool>(&rootInteractive)->default_value(1),"root interactive" )

    ;
  
  po::variables_map vm;
  po::store(po::command_line_parser(nargc, argv).
	    options(config).run(), vm);
  po::notify(vm);
  
  if( vm.count("help") ) {
    cout << config << endl;
    cout << "Usage: ./bin/MiniTreeFitter [options]" << endl;
    return 1;
  }
  
  TRint *interactive = 0;
  if( rootInteractive ) interactive = new TRint("rootInteractive",(int*)0,0);

  string categorisation = "notDefined";
  if( dirAFS.find( "mva" ) != string::npos ||  dirAFS.find( "MVA" ) != string::npos ) categorisation = "mva";
  if( dirAFS.find( "cic" ) != string::npos ||  dirAFS.find( "CiC" ) != string::npos ||
      dirAFS.find( "CIC" ) != string::npos )  categorisation = "cicpf";
  
  if( categorisation == "notDefined" ) {
    cout << "  categorisation not defined: directory name should contain mva or cic or cicpf" << endl;
    return 1;
  }



  // ----  Define the categories ---- // *************************************************************

  vector<TCut> smCategories;
  vector<TCut> melaCategories;
  vector<int>  polOrder;
 
  // temporary MELA categories, based on likelihood scatter plots
  melaCategories.push_back( "tagCat == 8 && abs(dEtaJJ) > 3.0 && mela_VBFvsgg > 0.6 && mela_SMvsPS_VBF < 0.2" );   // dominated by VBF 0m 
  melaCategories.push_back( "tagCat == 8 && abs(dEtaJJ) > 3.0 && ( (mela_VBFvsgg > 0.8 && mela_SMvsPS_VBF > 0.4 && mela_SMvsPS_VBF < 0.8) || (mela_VBFvsgg < 0.4 && mela_SMvsPS_VBF < 0.4) )" );   // dominated by ggH
  melaCategories.push_back( "tagCat == 8 && abs(dEtaJJ) > 3.0 && mela_VBFvsgg > 0.8 && mela_SMvsPS_VBF > 0.8" );   // dominated by VBF SM+ggH
 
   
  if( categorisation == "mva" ) {
    smCategories.push_back( "catMva == 0" ); polOrder.push_back( 5 );
    smCategories.push_back( "catMva == 1" ); polOrder.push_back( 5 );
    smCategories.push_back( "catMva == 2" ); polOrder.push_back( 5 );
    smCategories.push_back( "catMva == 3" ); polOrder.push_back( 5 );   
  } else if( categorisation == "cicpf" ) {
    for( int imelacat = 0; imelacat < melaCategories.size(); imelacat++ ) {
      smCategories.push_back( "maxSCEta < 1.49" && melaCategories[imelacat] );                  polOrder.push_back( 3 );        
      smCategories.push_back( "maxSCEta > 1.49 && minR9 > 0.94" && melaCategories[imelacat] );  polOrder.push_back( 3 );
    }
  }

  vector<int> polOrderCuts;
  vector<TCut> smCatCuts;
  if( cat >=  0 && cat < int(smCategories.size()) ) {    
    smCatCuts.push_back( smCategories[cat] );
    polOrderCuts.push_back(polOrder[cat]);
  }
  else  {
    smCatCuts = smCategories;
    polOrderCuts = polOrder;
  }

  // ----  prepare the fitter ---- // *************************************************************

  string hlFactoryCard = "etc/workspaceConfig/mva_sm_model_cond.rs";
  MiniTreeFitter1D fitter(hlFactoryCard);
  fitter.setPlotDirectory("workspace_VBF/" + categorisation + "/");
  // fitter.setMainCut( basicCut );

  string mName = "mass";
  //  RooRealVar *catMva   = new RooRealVar( "catMva"  ,"",-999,10.0); fitter.addVariable( catMva   );
  //  RooRealVar *diphoMva = new RooRealVar( "diphoMva","",-3  ,10.0); fitter.addVariable( diphoMva );
  RooRealVar *catBase   = new RooRealVar( "tagCat" ,"",-999,10.0); fitter.addVariable( catBase  );
  RooRealVar *dEtaJJ    = new RooRealVar( "dEtaJJ" ,"",-999,10.0); fitter.addVariable( dEtaJJ  );
  RooRealVar *minR9     = new RooRealVar( "minR9", "", -1.,10.);   fitter.addVariable( minR9 );
  RooRealVar *maxSCEta  = new RooRealVar( "maxSCEta", "", -1.,10.);  fitter.addVariable( maxSCEta );
  RooRealVar *mela_VBFvsgg   = new RooRealVar( "mela_VBFvsgg", "", -1.,10.);   fitter.addVariable( mela_VBFvsgg );
  RooRealVar *mela_SMvsPS_VBF   = new RooRealVar( "mela_SMvsPS_VBF", "", -1.,10.);   fitter.addVariable( mela_SMvsPS_VBF );
  RooRealVar *mass      = new RooRealVar( mName.c_str(),"",mMin,mMax); fitter.addVariable( mass     );
  fitter.setMassVarName( mName );
  fitter.setMassVarSet(true);

  fitter.setCategories( smCatCuts );
  fitter.setMassMin(mMin);
  fitter.setMassMax(mMax);

  string escaleSyst = "";  
  string categorySuffix = "_cat" + itostr(cat);
  string dirData = dirAFS + "/data/";
  string dirMC   = dirAFS + "/mc/" + escaleSyst;

  if( addSig ) {
    SMHiggsCrossSection HiggsXS; HiggsXS.is8TeV();
    float xsec_ggh = 1000*HiggsXS.HiggsSMxsec_ggh(mh); //fb
    float xsec_vbf = 1000*HiggsXS.HiggsSMxsec_vbf(mh); //fb
    float xsec_vh  = 1000*(HiggsXS.HiggsSMxsec_wh(mh) + HiggsXS.HiggsSMxsec_zh(mh)); //fb
    float xsec_tth = 1000*HiggsXS.HiggsSMxsec_tth(mh); //fb
    float br = HiggsXS.HiggsBR(mh);
    float lumi = 19.5; //fb-1
    cout << " br = " << br << " - xsec_ggh = " << xsec_ggh << endl;

    vector<string> sFiles1, sNames1;
    vector<float>  sXsec1;
    vector<string> sFiles2, sNames2;
    vector<float>  sXsec2;
    sFiles1.push_back( dirMC + "minitree_jhu_8TeV_SM0p_VBF_125p6_v3.root" );  sNames1.push_back( "VBF0p"); sXsec1.push_back(xsec_vbf);
    sFiles2.push_back( dirMC + "minitree_jhu_8TeV_0m_VBF_125p6_v3.root"   );  sNames2.push_back( "VBF0m"); sXsec2.push_back(xsec_vbf);

    float accSM = 0; float accPS = 0;
    for( size_t ifile = 0; ifile < sFiles1.size(); ifile++ ) {
      accSM += AccTot( sFiles1[ifile], smCatCuts  )*sXsec1[ifile];
    }
    for( size_t ifile = 0; ifile < sFiles2.size(); ifile++ ) {
      accPS += AccTot( sFiles2[ifile], smCatCuts )*sXsec2[ifile];
    } 
    for( size_t ifile = 0; ifile < sFiles2.size(); ifile++) {
      sXsec2[ifile] *= accSM/accPS;
    } 
    fitter.addSigSamples(sFiles1,sXsec1,sNames1,br,lumi);
    fitter.addSigSamples(sFiles2,sXsec2,sNames2,br,lumi);
    if( doFits) {    
      fitter.modelSignal(125,categorySuffix + sNames1[0], 0);
      fitter.modelSignal(125,categorySuffix + sNames2[0], 1);
      fitter.makeSignalWorkspace();
    }
  }
  
  float signi = -1;
  vector<double> signis;
  if( addBkg ) {
    string dataset = dirData + "data_8TeV_skimMVA_runABCD_withMELA.root";

    fitter.unblind();
    fitter.addData(dataset);
    if( doFits ) {
      if( bkgModel == 0 ){
	signis = fitter.modelBackground(125.,categorySuffix); 
	fitter.makeBackgroundWorkspace();
      }
      else if (bkgModel == 1) { 
	fitter.setPolynomialOrder(  polOrderCuts );
	fitter.modelBackground(   125.,categorySuffix );
	fitter.modelBackgroundPol(125.,categorySuffix );
	fitter.modelBackgroundExp(125.,categorySuffix );
	// fitter.makeBackgroundWorkspace("Pol");
	  // fitter.makeBackgroundWorkspace();
	fitter.dumpBkgFitParam();
	if( simFit ) {
	  fitter.backgroundFitResut(cat,0);
	}
      }
    }
  }
  //  return;
  if( addSig && addBkg && simFit ) fitter.simultaneousFitOnlyOneSig(0);   // CF: WHY ONLY ONE SIGNAL HYPOTHESIS?!  what does this do exactly?

  if( addSig && addBkg && doFits ) {
    string datacard = "MVAFact_SM_mh125_cat" + itostr(cat) + ".datacard.txt";
    fitter.createDataCard(datacard);    
  }
  signi = 0;
  for( unsigned ic = 0 ;ic < signis.size(); ic++ ) {
    cout << " signi cat[" << ic << "] = " << signis[ic] << endl;
    signi += signis[ic]*signis[ic];
  }
  cout << " SigniTot: " << sqrt(signi) << endl;

  if( interactive ) {
    interactive->Run();
    delete interactive;
  }
  
  return 0;
}



float AccTot(string rootfile, const vector<TCut>& cuts, string treename ) {
  TTree *t =0; t = (TTree*) TFile::Open(rootfile.c_str())->Get(treename.c_str());
  float ntot = 0;
  if( t ) {
    for( unsigned ic = 0 ; ic < cuts.size()     ; ++ic ) {
      TH1F h("htmpMass","htmpMass",80,100,180);
      t->Draw("mass>>htmpMass","wei"*cuts[ic],"goff");
      ntot += h.Integral();
      cout << "AccCat: " << h.Integral()/(1e5) << "  <-->  " << cuts[ic] << endl;
    }
  }
  return ntot/float(1e5);
}

