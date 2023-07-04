//g++ -O3 -o tmva-test-allYears tmva-test-allYears.cc  `root-config --cflags --libs` -lRooFit -lGenVector -lMathCore -lTMVA -l TMVAGui 
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sys/stat.h>
#include <dirent.h>
#include "Riostream.h"
#include <map>
#include <string>
#include "RooGlobalFunc.h"
#include <vector>
#include <math.h>
#include <TGenericClassInfo.h> 
#include <TH1.h>
#include <TH2.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TROOT.h>
#include <TEnv.h>
#include <TSystem.h>
#include <TTree.h>
#include "TBranch.h"
#include <TApplication.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TStyle.h> 
#include <TF1.h>  
#include <TF2.h> 
#include <TLorentzVector.h>
#include <TVector3.h>
#include "TDSet.h"
#include "TChain.h"
#include <time.h> 
#include <TSystemDirectory.h>
#include <TAttLine.h>
#include "TRatioPlot.h"
#include "TMVA/Config.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Factory.h"
#include "TMVA/BDT.h"
#include "TMVA/BDT_Reg.h"
#include "TMVA/MethodBase.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Reader.h"
#include "TMVA/PDF.h"
#include <sys/time.h>
#include <sys/times.h>
#include <iostream>
// System stuff
#include <fstream> 
#include <sys/time.h>
#include <sys/times.h>
timeval startTime, stopTime, totalTime;
clock_t startCPU, stopCPU; 
tms startProc, stopProc; 
void tmva_test_dnn();
void tmva_evaluate_dnn(bool save=false);
void replaceAll(std::string& str, const std::string& from, const std::string& to) ;
std::map<std::string, std::string> splitFeatures(std::string FEATURES);

double mc_sigma = 0.040;
double mc_mass  = 5.27783;
double JPsiMass = 3.096916;
double nSigma_psiRej = 3.;
double tagged_mass=0;
double nsig_sw=0;
double mumuMass=-99;
double mumuMassE=-99;
double kstTrk1Pt;
double kstTrk2Pt;
double kstTrk1Eta;
double kstTrk2Eta;
float  weight=-99;
double DRweight=-99;
double ScoreDNN=-1.;
int    xcut=-99;
double dxcut=-99;
int    passB0Psi_jpsi=-99;

//
//  Long64_t CutEve =800000000;
Long64_t CutEve =1000000;
double CutEffDATA = 0.036;
double CutEffMC   = 0.093/2;
//double CutEffDATA = 0.00060324;
//double CutEffMC   = 0.00035352;
//Long64_t CutEve =500000000000000000; 
//Long64_t CutEve =5000; 
//
   
int BDTNumTree = 500;
double kstTrk1PtD=-99;
double kstTrk2PtD=-99;
double kstTrk1EtaD=-99;
double kstTrk2EtaD=-99;
Double_t truthMatchMum=-99 ;
Double_t truthMatchMup=-99 ; 
Double_t truthMatchTrkm =-99;
Double_t truthMatchTrkp=-99 ;
Double_t tagB0=-99 ;
Double_t genSignal=-99 ;
Double_t trig=-99 ;
Long64_t eventN=-99;
int pass_preselection=-99;
int num_threads=30; 
float dnncut=0.97;   
//
std::string BookTXT="";
//
std::string isFeature="";
std::string year=""; 
std::string year_default="2016";
std::string NameFileMCp0 = \
"/gwteray/users/fiorendi/p5prime/data2018/flat_ntuples/pre_selection_forPaolo/2018_MC_LMNR_add_vars_fixDR_passPreselection.root";
std::string NameFileMCp1 = \
"/gwteray/users/fiorendi/p5prime/data2018/flat_ntuples/pre_selection_forPaolo/2018_MC_LMNR_add_vars_fixDR_passPreselection.root";
std::string NameFileDatap0 = \
"/gwteray/users/fiorendi/p5prime/data2018/flat_ntuples/pre_selection_forPaolo/data_2018_LMNR_plus_Charmonium_add_vars.root";
std::string NameFileDatap1 = \
"/gwteray/users/fiorendi/p5prime/data2018/flat_ntuples/pre_selection_forPaolo/data_2018_LMNR_plus_Charmonium_add_vars.root";
// std::string NameFileModel =  "/gwpool/users/dini/p5prime/DRweights/2016MC_JPSI_PDGinputs.root";
std::string datasetYear;
std::string dirfilexml;
std::string datasetname = "dataset";
std::string filexml= "/weights/TMVAClassification_DNN_CPU.weights.xml";
std::string OutputFileName="TMVA_ClassificationOutput.root";
/* //d::vector<std::string> features = {"bVtxCL", "bLBS", "bLBSE" , "bDCABS","bDCABSE","kstTrk1Pt", "kstTrk2Pt",\
//                                   "kstTrk1Eta", "kstTrk2Eta","kstTrk1DCABS","kstTrk1DCABSE",\
//				     "kstTrk2DCABS", "kstTrk2DCABSE","mu1Pt","mu2Pt","mu1Eta","mu2Eta","sum_isopt_04"};
 */ 
/*    std::vector<std::string> features = {"bVtxCL", "bLBS", "bLBSE" ,"bCosAlphaBS", "bDCABS","bDCABSE","kstTrk1Pt", "kstTrk2Pt",\
                                        "kstTrk1Eta", "kstTrk2Eta","kstTrk1DCABS","kstTrk1DCABSE",\
   				     "kstTrk2DCABS", "kstTrk2DCABSE","mu1Pt","mu2Pt","mu1Eta","mu2Eta","sum_isopt_04"};
 */   
std::vector<std::string> spectatorF = {"mumuMass","xcut","trig","truthMatchMum","truthMatchMup","truthMatchTrkm","truthMatchTrkp"};
//std::vector<std::string> spectatorF = {"mumuMass","tagged_mass","xcut","trig","truthMatchMum","truthMatchMup","truthMatchTrkm","truthMatchTrkp"};
std::vector<std::string> spectatorI = {"eventN","trig","truthMatchMum","truthMatchMup","truthMatchTrkm","truthMatchTrkp"};
//std::vector<std::string> spectatorI = {"eventN","pass_preselection","passB0Psi_jpsi","trig","truthMatchMum","truthMatchMup","truthMatchTrkm","truthMatchTrkp"};
std::vector<std::string> variables = {};
std::vector<std::string> vartested = {"tagged_mass"};
//std::vector<std::string> vartested = {"cos_theta_l","cos_theta_k","phi_kst_mumu"};
   

  std::vector<std::string> features = {"kstarmass","bVtxCL" , "bCosAlphaBS" ,\
  				   "bLBS/bLBSE","bDCABS/bDCABSE", "kstTrkmDCABS/kstTrkmDCABSE",
  				   "kstTrkpDCABS/kstTrkpDCABSE","sum_isopt_04"};
/*
//std::vector<std::string> features = {"kstarmass","bVtxCL" , "bCosAlphaBS" ,\
//				       "bLBS","bDCABS", "kstTrkmDCABS",
//				       "kstTrkpDCABS","sum_isopt_04"};
//
//   				     "mu1Pt","mu2Pt", "mu1Eta","mu2Eta","kstTrk2Eta","kstTrk1Eta",
// std::vector<std::string> vartested = {"cos_theta_l","cos_theta_k","phi_kst_mumu","bCosAlphaBS"};
//    
// 
// std::vector<std::string> features = {"kstTrk1Pt" , "kstTrk2Pt" ,\
//                                      "kstTrk1Eta", "kstTrk2Eta",\
//    				     "mu1Pt","mu2Pt","mu1Eta","mu2Eta",
// 				     "bCosAlphaBS","kstTrk1DCABS","kstTrk2DCABS","sum_isopt_04"};
 */


/* std::vector<std::string> features = {"bVtxCL","bLBS",  "bCosAlphaBS", "bDCABS","kstTrk1Pt", "kstTrk2Pt",\
                                        "kstTrk1Eta", "kstTrk2Eta","kstTrk1DCABS","kstTrk1DCABSE",\
   				     "kstTrk2DCABS", "kstTrk2DCABSE","mu1Pt","mu2Pt","mu1Eta","mu2Eta","sum_isopt_04"};
 */
 //  std::vector<std::string> features = {"kstTrk1Pt", "kstTrk2Pt","kstTrk1Eta", "kstTrk2Eta"};
//  std::vector<std::string> features = {"kstTrk1Pt", "kstTrk2Pt","kstTrk1Eta", "kstTrk2Eta","bCosAlphaBS"};
std::vector<TCut> RemoveNone;
std::vector<TH1D*> HistData;
std::vector<TH1D*> HistMC;
std::vector<TH1D*> HistData_Cut;
std::vector<TH1D*> HistMC_Cut;
std::vector<TH1D*> HistMCW;
std::vector<TCanvas*> cstudies;			     
std::vector<TRatioPlot*> RatiosDataMC;			     
std::vector<TRatioPlot*> RatiosDataMCW;

TCut mycuts;
TCut mycutb;
TCut mycuts0;
TCut mycutb0;
TCut Cut_sig_mass;    
TCut Cut_bkg_mass;    
TCut Cut_ct	 ;    
TCut Cut_truth_match; 

std::map<unsigned int, std::map<std::string,std::string>> mapFeat;

std::vector<int> iFeatDiv(20,-1);

			     
int main (int argc, char** argv) {


  if( argc<=1 ){
    std::cout<<Form("Usage: %s [year=2016,2017,2018] {dnn}",argv[0])<<std::endl;
    std::cout<<Form("Please, set the year (at least)")<<std::endl;
    std::cout<<Form("example: %s 2016 bdt; if you want to train a new dnn & plots\n",argv[0])<<std::endl;
    std::cout<<Form("         %s 2016	 ; if you want just produce the  plots\n",argv[0])<<std::endl;
    exit(0);
  }else{
   if ((strcmp(argv[1],"2016") == 0 || \
        strcmp(argv[1],"2017") == 0 || \
	strcmp(argv[1],"2018") == 0)){
      year=argv[1];
//      replaceAll( NameFileModel ,  year_default, year);
      replaceAll( NameFileMCp0  ,  year_default, year);
      replaceAll( NameFileMCp1  ,  year_default, year);
      replaceAll( NameFileDatap0,  year_default, year);
      replaceAll( NameFileDatap1,  year_default, year);
      datasetYear=datasetname+year;
      dirfilexml=datasetYear+filexml;
      ROOT::EnableImplicitMT(num_threads);
      TMVA::Tools::Instance();
      std::cout<<Form("Setting year=%s",year.c_str())<<std::endl;
      std::cout<<Form("Setting dataset year=%s",datasetYear.c_str())<<std::endl;
       if (strcmp(argv[1],"2016") != 0){
	BDTNumTree=1000;
       }
    }else{
     std::cout<<Form("not recognize year=%s",year.c_str())<<std::endl;
     exit(0);
    }
	
  }  
   
  if( argc>2 ){
   if ((strcmp(argv[2],"dnn") == 0)){
    std::cout<<Form("Start training the DNN \n")<<std::endl;
    tmva_test_dnn();
    tmva_evaluate_dnn(false);
   }
   if ((strcmp(argv[2],"save") == 0)){
    std::cout<<Form("Save weights \n")<<std::endl;
    tmva_evaluate_dnn(true);
   }
  }else{
    tmva_evaluate_dnn(false);
  }
}
 void tmva_test_dnn(){
  gROOT ->Reset();
  gROOT->SetStyle("Plain");
  ROOT::EnableImplicitMT(num_threads); 
  TMVA::Tools::Instance();
  std::cout<<Form("Warning, analyze only %lli events!!!!\n",CutEve)<<std::endl;
  std::cout<<Form("Warning, analyze only %lli events!!!!\n",CutEve)<<std::endl;
  std::cout<<Form("Warning, analyze only %lli events!!!!\n",CutEve)<<std::endl;
  std::cout<<Form("Warning, analyze only %lli events!!!!\n",CutEve)<<std::endl;
  std::cout<<Form("Warning, analyze only %lli events!!!!\n",CutEve)<<std::endl;




//   TFile *fModel = new TFile(NameFileModel.c_str(),"READ");
//   TTree *TreeModel     = (TTree*)fModel->Get("ntuple_DRw");
//   std::cout<<Form("Open Model File :%s",NameFileModel.c_str())<<std::endl;
//   int nentriesModel = (int)TreeModel->GetEntries();
//   std::cout<<Form("Found Model entries= = %d\n",nentriesModel)<<std::endl;


//  TH1D* HxMassWn     = new TH1D( "HxMassWn"          , "B^{0} Mass sPlot <S0",100, 5., 5.6);
  
  TFile *fData = new TFile(NameFileDatap0.c_str(),"READ");
  TTree *TreeData     = (TTree*)fData->Get("ntuple");
  std::cout<<Form("Open Data File :%s",NameFileDatap0.c_str())<<std::endl;
  int nentriesData = (int)TreeData->GetEntries();
  std::cout<<Form("Found Data entries= = %d\n",nentriesData)<<std::endl;
  TFile *fMC   = new TFile(NameFileMCp0.c_str(),"READ");
  std::cout<<Form("Open MC File :%s",NameFileMCp0.c_str())<<std::endl;
  TTree *TreeMC     = (TTree*)fMC->Get("ntuple");
  int nentriesMC = (int)TreeMC->GetEntries();
//  TreeMC->AddFriend(TreeModel);
//  TreeMC->AddFriend("ntuple_DRw",NameFileModel.c_str());
  std::cout<<Form("Found MC entries= = %d\n",nentriesMC)<<std::endl;
  DIR* dir = opendir(datasetYear.c_str());
  
  if(!dir) {
   gSystem->Exec(Form("mkdir -p %s",datasetYear.c_str()));
   std::cout<<Form("Create Dir  %s",datasetYear.c_str())<<std::endl;
  } 
  auto outputFile = TFile::Open(Form("%s/%s",datasetYear.c_str(),OutputFileName.c_str()), "RECREATE");

//  TMVA::Factory factory("TMVAClassification", outputFile,"!V:ROC:!Silent:Color:!DrawProgressBar:AnalysisType=Classification" ); 
//  TMVA::Factory factory =  TMVA::Factory( "TMVAClassification", outputFile,"!V:!Silent:Color:DrawProgressBar:Transformations=None:AnalysisType=multiclass" );
//TMVA::Factory factory =  TMVA::Factory( "TMVAClassification", outputFile,"!V:!Silent:Color:DrawProgressBar:Transformations=I;P;G:AnalysisType=multiclass" );

//  TMVA::Factory factory =  TMVA::Factory( "TMVAClassification", outputFile,"!V:!Silent:Color:DrawProgressBar:Transformations=None:AnalysisType=Auto" );
  TMVA::Factory *factory =  new TMVA::Factory( "TMVAClassification", outputFile,"V:!Silent:Color:DrawProgressBar:Transformations=None:AnalysisType=Auto");
//  TMVA::Factory factory =  TMVA::Factory( "TMVAClassification", outputFile,"!V:!Silent:Color:DrawProgressBar:Transformations=None:AnalysisType=Auto");
//  TMVA::Factory factory =  TMVA::Factory( "TMVAClassification", outputFile,"!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Auto");

//  TMVA::Factory factory =  TMVA::Factory( "TMVAClassification", outputFile,"!V:!Silent:Color:DrawProgressBar:Transformations=None:AnalysisType=multiclass" );
//    TMVA::Factory factory =  TMVA::Factory( "TMVAClassification", outputFile,"!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=multiclass" );
//  TMVA::Factory factory =  TMVA::Factory( "TMVAMulticlass", outputFile,"!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=multiclass" );
//  TMVA::Factory factory("TMVAClassification", outputFile,"!V:ROC:!Silent:Color:!DrawProgressBar:AnalysisType=Classification" ); 
  TMVA::DataLoader * loader = new TMVA::DataLoader(datasetYear.c_str());
//  TMVA::DataLoader * loader = new TMVA::DataLoader(Form("dataset");
//global event weights per tree (see below for setting event-wise weights)
  Double_t signalWeight     = 1.0;
  Double_t backgroundWeight = 1.0;
  
double mc_sigma = 0.0400;
double mc_mass  = 5.27783;   
// You can add an arbitrary number of signal or background trees
Long64_t CutEveMC=(Long64_t)CutEve/CutEffMC;
Long64_t CutEveDATA=(Long64_t)CutEve/CutEffDATA;
std::cout<<Form("mycuteveMC   eventN<%lld",CutEveMC) <<std::endl; 
std::cout<<Form("mycuteveDATA eventN<%lld",CutEveDATA) <<std::endl; 
TCut mycuteveMC   = Form("eventN<%lld",CutEveMC); 
TCut mycuteveDATA = Form("eventN<%lld",CutEveDATA);

Cut_sig_mass	 = Form("(tagged_mass > %f - 2.5*%f) && (tagged_mass < %f + 2.5*%f)",mc_mass,mc_sigma,mc_mass,mc_sigma);

Cut_ct  	 = "( ( (tagB0==1) && (genSignal==1)) || ( (tagB0==0) && (genSignal==2) ) )";

Cut_bkg_mass	 = Form("( ( (tagged_mass > %f-7.*%f ) && (tagged_mass < %f-3.*%f ) ) || ( (tagged_mass > %f+3.*%f ) && (tagged_mass < %f+7*%f ) ) )",mc_mass,mc_sigma,mc_mass,mc_sigma,mc_mass,mc_sigma,mc_mass,mc_sigma );

Cut_truth_match  = "( (truthMatchMum==1) && (truthMatchMup ==1) && (truthMatchTrkm==1) && (truthMatchTrkp==1) )";

mycuts0 	 = Cut_sig_mass  && Cut_truth_match &&  Cut_ct && "trig == 0 && (mumuMass < 2.702)";

mycutb0 	 = Cut_bkg_mass &&  "(mumuMass < 2.702)";

//   TCut Jpsi_Cut="(mumuMass*mumuMass<8.68  || mumuMass*mumuMass > 10.09)";
//   TCut Psi2SCut="(mumuMass*mumuMass<12.86 || mumuMass*mumuMass > 14.18)";
// 
//   TCut PreCut = "(tagged_mass > 5.0 && tagged_mass < 5.6)"&&Jpsi_Cut&&Psi2SCut; 
//   TCut SBCut  = "(tagged_mass > 5.0 && tagged_mass < 5.1)||(tagged_mass > 5.44 && tagged_mass < 5.6)"; 
// //   TCut mycuts0 = mycuteveMC  &&PreCut&&"int(eventN%2)==0&&(trig==1) && (truthMatchMum==1) && (truthMatchMup==1) && (truthMatchTrkm==1) && (truthMatchTrkp==1)";
// //   TCut mycutb0 = mycuteveDATA&&PreCut&&"int(eventN)%2==0"&&SBCut; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
// //  TCut mycuts0 = PreCut&&"int(eventN)%2==0";
//   TCut mycuts0 = mycuteveMC  &&PreCut&&"(int(trig)==1) && (int(truthMatchMum)==1) && (int(truthMatchMup)==1) && (int(truthMatchTrkm)==1) && (int(truthMatchTrkp)==1)";
//   TCut mycutb0 = mycuteveDATA&&PreCut&&SBCut; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
//   TCut mycuts0 = mycuteveMC  &&PreCut&&"int(eventN)%2==0&&(int(trig)==1) && (int(truthMatchMum)==1) && (int(truthMatchMup)==1) && (int(truthMatchTrkm)==1) && (int(truthMatchTrkp)==1)";
//   TCut mycutb0 = mycuteveDATA&&PreCut&&"int(eventN)%2==0"&&SBCut; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
  TCut RMNone="";
  for (unsigned int i=0;i<=features.size()-1;i++) {
   RemoveNone.push_back(TCut(Form("!TMath::IsNaN(%s)",features[i].c_str())));
   if(i==0) RMNone=RemoveNone[i];
   if(i>0)  RMNone=RMNone&&RemoveNone[i];
  }
  TCut mycutNEve="int(eventN)<2000000";
//   TCut mycuts=mycuts0&&RMNone&&mycutNEve;
//   TCut mycutb=mycutb0&&RMNone&&mycutNEve;
  mycuts=mycuts0&&RMNone;
  mycutb=mycutb0&&RMNone;
//  TCut mycutb=mycutb0&&RMNone&&"!TMath::IsNaN(DRweight)";
  TTree *TreeData_Cut = TreeData->CopyTree(mycutb);
  TTree *TreeMC_Cut   = TreeMC  ->CopyTree(mycuts);
  int nentriesData_Cut= TreeData_Cut->GetEntries();
  int nentriesMC_Cut  = TreeMC_Cut->GetEntries();
  std::cout<<Form("Found MC   entries after Cuts = %d",nentriesMC_Cut)<<std::endl;
  std::cout<<Form("Found data entries after Cuts = %d",nentriesData_Cut)<<std::endl;
  std::cout<<Form("Cuts [MC]   = %s \n", mycuts.GetTitle ())<<std::endl;
  std::cout<<Form("Cuts [Data] = %s \n", mycutb.GetTitle ())<<std::endl;
//   loader->AddTree( TreeMC    ,"Signal"	  , signalWeight    , mycuts );
//   loader->AddTree( TreeData  ,"Background",backgroundWeight , mycutb );
  loader->AddTree( TreeMC_Cut    ,"Signal"	  , signalWeight    );
  loader->AddTree( TreeData_Cut  ,"Background",backgroundWeight );
//  loader->AddTree( TreeData,"Signal"	, signalWeight , mycuts   );
//  loader->AddTree( TreeMC, "Background",backgroundWeight ,mycutb );
  for (unsigned int i=0;i<=features.size()-1;i++) {
   loader->AddVariable( features[i].c_str(), 'F' );
  }
//  TreeMC_Cut->Print(mycuts);
//  exit(1);
//   for (unsigned int i=0;i<=spectatorF.size()-1;i++) {
//    loader->AddSpectator( spectatorF[i].c_str(), 'F' );
//   }
//   for (unsigned int i=0;i<=spectatorI.size()-1;i++) {
//    loader->AddSpectator( spectatorI[i].c_str(), 'I' );
//  }
//   
//   loader->AddVariable( "kstTrk1Pt", 'D' );
//   loader->AddVariable( "kstTrk2Pt", 'D' );
//   loader->AddVariable( "kstTrk1Eta", 'D' );
//   loader->AddVariable( "kstTrk2Eta", 'D' );
  
//   loader->AddVariable( "kstTrk1Pt", 'F' );
//   loader->AddVariable( "kstTrk2Pt", 'F' );
//   loader->AddVariable( "kstTrk1Eta", 'F' );
//   loader->AddVariable( "kstTrk2Eta", 'F' );
//   loader->AddSpectator( "truthMatchMum",  'F' );
//   loader->AddSpectator( "truthMatchMup",  'F' );
////////////  loader->SetSignalWeightExpression("weight");
////////////   loader->SetBackgroundWeightExpression(backgroundWeight);
//  loader->AddVariable( "nsig_sw", 'D' );
// Apply additional cuts on the signal and background samples (can be different)
// Tell the factory how to use the training and testing events
//
// If no numbers of events are given, half of the events in the tree are used 
// for training, and the other half for testing:
//    loader->PrepareTrainingAndTestTree( mycut, "SplitMode=random:!V" );
// To also specify the number of testing events, use:
//    loader->PrepareTrainingAndTestTree( mycut,
//                                         "NSigTrain=3000:NBkgTrain=3000:NSigTest=3000:NBkgTest=3000:SplitMode=Random:!V" );
//   loader->PrepareTrainingAndTestTree( "", "",
//   loader->PrepareTrainingAndTestTree( mycuts, mycutb,
//                                  "nTrain_Signal=5000:nTrain_Background=5000:SplitMode=Random:NormMode=NumEvents:!V" );
//                                    "nTrain_Signal=500000:nTrain_Background=500000:SplitMode=Random:NormMode=NumEvents:!V" );
//                                      "SplitMode=Random:NormMode=NumEvents:!V" );
   std::cout<<"==== PrepareTrainingAndTestTree ===="<<std::endl;
//   loader->PrepareTrainingAndTestTree( "","",
   loader->PrepareTrainingAndTestTree( "", "",
//   loader->PrepareTrainingAndTestTree( mycuts, mycutb,
//                                  "nTrain_Signal=5000:nTrain_Background=5000:SplitMode=Random:NormMode=NumEvents:!V" );
//                                    "nTrain_Signal=5000:nTrain_Background=5000:SplitMode=Random:NormMode=NumEvents:!V" );
                                      "SplitMode=Random:NormMode=NumEvents:V" );
////                                      "SplitMode=Random:NormMode=EqualNumEvents:V" );
//Boosted Decision Trees
   BookTXT =\
/* //Form("H:V:NTrees=%d:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.50:nCuts=20:MaxDepth=2:DoBoostMonitor:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=100:NsmoothMVAPdf=10",BDTNumTree);
// Form("V:ErrorStrategy=CROSSENTROPY:VarTransform=N:WeightInitialization=XAVIERUNIFORM:Layout=TANH|128,TANH|128,TANH|128,LINEAR:TrainingStrategy=LearningRate=1e-2,Momentum=0.9,\
//                                         ConvergenceSteps=20,BatchSize=100,TestRepetitions=1,\
//                                         WeightDecay=1e-4,Regularization=None,\
//                                         DropConfig=0.0+0.5+0.5+0.5:nTrain_Signal=10000:nTrain_Background=10000:nTest_Signal=6000:nTest_Background=6000:Architecture=CPU:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=100:NsmoothMVAPdf=10");
 */Form("V:ErrorStrategy=CROSSENTROPY:VarTransform=N:WeightInitialization=XAVIERUNIFORM:Layout=TANH|128,TANH|128,TANH|128,LINEAR:TrainingStrategy=LearningRate=1e-2,Momentum=0.9,\
                                        ConvergenceSteps=20,BatchSize=128,TestRepetitions=1,\
                                        WeightDecay=1e-4,Regularization=None,\
                                        DropConfig=0.0+0.0+0.0+0.0:Architecture=CPU");

      TString inputLayoutString = "InputLayout=0|0|0"; 
      TString batchLayoutString= "BatchLayout=0|0|0";
      TString layoutString ("Layout=DENSE|64|TANH,BNORM|0.99|0.0001,DENSE|64|TANH,BNORM|0.99|0.0001,DENSE|64|TANH,BNORM|0.99|0.0001,DENSE|64|TANH,DENSE,BNORM|0.99|0.0001,DENSE|1|LINEAR");
      // Define Training strategies 
      // one can catenate several training strategies 
      TString training1("LearningRate=1e-6,Momentum=0.9,Repetitions=1,"
                        "ConvergenceSteps=30,BatchSize=128,TestRepetitions=1,"
                        "MaxEpochs=100,WeightDecay=1e-4,Regularization=None,"
                        "Optimizer=ADAM,DropConfig=0.0+0.0+0.0+0.:ValidationSize=0.50");
      //     TString training2("LearningRate=1e-3,Momentum=0.9,Repetitions=1,"
      //                       "ConvergenceSteps=10,BatchSize=128,TestRepetitions=1,"
      //                       "MaxEpochs=20,WeightDecay=1e-4,Regularization=None,"
      //                       "Optimizer=SGD,DropConfig=0.0+0.0+0.0+0.");
  
      TString trainingStrategyString ("TrainingStrategy=");
      trainingStrategyString += training1; // + "|" + training2;

      // General Options.

      TString dnnOptions ("H:V:ErrorStrategy=CROSSENTROPY:VarTransform=G:"
                          "WeightInitialization=XAVIER");
      dnnOptions.Append (":"); dnnOptions.Append (inputLayoutString);
      dnnOptions.Append (":"); dnnOptions.Append (batchLayoutString);
      dnnOptions.Append (":"); dnnOptions.Append (layoutString);
      dnnOptions.Append (":"); dnnOptions.Append (trainingStrategyString);
      dnnOptions += ":Architecture=CPU:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10";
//      dnnOptions += ":Architecture=GPU:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10";

   std::cout<<dnnOptions<<std::endl;
   std::cout<<"==== BookMethod ===="<<std::endl;
   factory->BookMethod( loader,  TMVA::Types::kDL, "DNN_CPU",BookTXT.c_str() );
//   factory->BookMethod( loader,  TMVA::Types::kDNN, "DNN_CPU",BookTXT.c_str() );
//   factory->BookMethod( loader,  TMVA::Types::kDNN, "DNN",dnnOptions);
//   factory->BookMethod( loader,  TMVA::Types::kDNN, "DNN_CPU",BookTXT.c_str() );



///////////////   factory->OptimizeAllMethods();
   std::cout<<"==== Start Training ===="<<std::endl;
   factory->TrainAllMethods();

   std::cout<<"==== Start Testing ===="<<std::endl;
   factory->TestAllMethods();  

   std::cout<<"==== Start Evalueting ===="<<std::endl;
   factory->EvaluateAllMethods();  
   TCanvas *TRoc = factory->GetROCCurve(loader);
//    c1->Draw();
//    c1->Print("c1-tmva-test-allYears-2016.pdf");
   std::cout<<"==== Close Output file ===="<<std::endl;
   outputFile->Close();   
//   TCanvas *TRoc = factory->GetROCCurve (loader);
   TRoc->Print(Form("%s/roc-tmva-test-allYears-%s.pdf",datasetYear.c_str(),year.c_str()));
 } 
//
//=================================================================================================================================================================
//   
 void tmva_evaluate_dnn(bool save){
 
   TFile *foutData = 0;
   TFile *foutMC   = 0;
   TFile *pout	   = 0;
   TTree *toutMC   = 0;
   TTree *toutData = 0;
//   TBranch* bDataScore=0;
//   TBranch* bMCScore=0;
   if(save) {
    std::cout<<"Save the weights!!!"<<std::endl;
    std::cout<<"The parity is not set!!!"<<std::endl;
    foutMC = new TFile(Form("%sMC_LMNR-ScoreDNN.root",year.c_str()),"RECREATE");
    toutMC = new TTree("sDNN", "score DNN");
    toutMC->Branch("eventN", &eventN,"eventN/L");
    toutMC->Branch("ScoreDNN", &ScoreDNN,"ScoreDNN/D");
    foutData = new TFile(Form("%sData_LMNR-ScoreDNN.root",year.c_str()),"RECREATE");
    toutData = new TTree("sDNN", "score DNN");
    toutData->Branch("eventN", &eventN,"eventN/L");
    toutData->Branch("ScoreDNN", &ScoreDNN,"ScoreDNN/D");
   }else{
    std::cout<<"Plots with the opposite parity!!!"<<std::endl;
    pout = new TFile(Form("%sDATA_LMNR-Plots.root",year.c_str()),"RECREATE");
   }
 
   gROOT ->Reset();
   gROOT->SetStyle("Plain");

   TH1D* HxDNNData   	 = new TH1D( "HxDNNData"	, Form("DNN Output Data %s",year.c_str()) ,100, 0., 1.0);
   TH1D* HxDNNMC   	 = new TH1D( "HxDNNMC"		, Form("DNN Output MC %s",year.c_str()) ,100, 0., 1.0);
//   TH1D* HxDNNMCW   	 = new TH1D( "HxDNNMCW"		, Form("Bdt Output MC %s reweighted",year.c_str()) ,100, 0., 1.0);
   std::vector<double> VarD;			     
   std::vector<float> VarF;

//    ROOT::EnableImplicitMT(num_threads);
//    TMVA::Tools::Instance();
//

   for (unsigned int i=0;i<=features.size()-1;i++) {
    VarD.push_back(0);
    VarF.push_back(0);
     std::cout <<features[i]<< std::endl;
     std::map<std::string,std::string> temap= splitFeatures(features[i]);
     mapFeat.insert(make_pair(i,temap));
     
   }
   TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );
   
   int count=0;
   int idivi=0;
   for (std::map<unsigned int, std::map<std::string,std::string>>::iterator imap = mapFeat.begin();
   			      imap != mapFeat.end();
   			      ++imap)
   {
    for ( std::map<std::string,std::string>::iterator imap2 = (*imap).second.begin();
   			      imap2 != (*imap).second.end();
   			      ++imap2)
    {
     variables.push_back((*imap2).first);
     std::cout<<"Insert in variables (1): "<<(*imap2).first<<std::endl;
     iFeatDiv[idivi]=idivi;
     if ((*imap2).second!="") {
      VarD.push_back(0);
      iFeatDiv[idivi]=idivi+1;
      variables.push_back((*imap2).second);
      std::cout<<"Insert in variables (2): "<<(*imap2).second<<std::endl;
      idivi++;
     } 
     cstudies.push_back(new TCanvas(Form("c_%s",features[count].c_str()),Form("MC %s DNN studies %s (DNN feature)",features[count].c_str(),year.c_str()),200,10,900,780));
     count++;
     idivi++;
    } 
//    std::cout <<"mappa->"<< (*imap).first<<" = "<<(*imap).second << std::endl;
   }
  
//    for (unsigned int i=0;i<=features.size()-1;i++) { 
//     VarD.push_back(0);
//     VarF.push_back(0);
//     variables.push_back(features[i]);
//     cstudies.push_back(new TCanvas(Form("c_%s",features[i].c_str()),Form("MC %s DNN studies %s (DNN feature)",features[i].c_str(),year.c_str()),200,10,900,780));
//    }
   for (unsigned int i=0;i<=vartested.size()-1;i++) { 
    VarD.push_back(0);
    variables.push_back(vartested[i]);
    cstudies.push_back(new TCanvas(Form("c_%s",vartested[i].c_str()),Form("MC %s DNN studies %s (variable check)",vartested[i].c_str(),year.c_str()),200,10,900,780));
   }
   
   for (unsigned int i=0;i<=features.size()-1;i++) { 
    reader->AddVariable( features[i].c_str(), &VarF[i]);
//    cstudies.push_back(new TCanvas(Form("c_%s",features[i].c_str()),Form("MC %s reweighting studies %s",features[i].c_str(),year.c_str()),200,10,900,780));
   }

//   TString weightfile = dir + prefix  + TString(".weights.xml");
   FILE *file = fopen(dirfilexml.c_str(), "r");
   if (file){
    std::cout<<"Reading DNN file "<<dirfilexml.c_str()<<std::endl;
    reader->BookMVA( "DNN",dirfilexml.c_str() );
   }else{
    std::cout<<dirfilexml.c_str()<<" not found; exit"<<std::endl;
    exit(1);
   } 
  std::cout<<Form("Size VarD %d sould be equal to %d",int(VarD.size()),int(variables.size()))<<std::endl;
  std::cout<<Form("Size VarF %d should be equal to Size Features %d ",int(VarF.size()),int(features.size()))<<std::endl;

//=============== Model TTree ================== 

//   TFile *fModel = new TFile(NameFileModel.c_str(),"READ");
//   TTree *TreeModel     = (TTree*)fModel->Get("ntuple_DRw");
//   std::cout<<Form("Open Model File :%s \n",NameFileModel.c_str())<<std::endl;
//   int nentriesModel = (int)TreeModel->GetEntries();
//   std::cout<<Form("Found Model entries= = %d",nentriesModel)<<std::endl;
//   TreeModel->SetBranchAddress("DRweight"     ,&DRweight );

//=============== MC TTree ================== 
   TFile *fMC   = new TFile(NameFileMCp1.c_str(),"READ");
   std::cout<<Form("Opening MC File :%s \n",NameFileMCp1.c_str())<<std::endl;
   TTree *TreeMC     = (TTree*)fMC->Get("ntuple");
//   if (save) bMCScore = TreeMC->Branch("ScoreDNN", &ScoreDNN,"ScoreDNN/D");
   
   

   for (unsigned int i=0;i<=variables.size()-1;i++) { 
//   for (unsigned int i=0;i<=features.size()-1;i++) { 
     TreeMC->SetBranchAddress(variables[i].c_str()     ,&VarD[i] );
     float hMin = TreeMC->GetMinimum(variables[i].c_str());
     float hMax = TreeMC->GetMaximum(variables[i].c_str());
     
     
     if(variables[i]=="tagged_mass") hMin=5.0;
     if(variables[i]=="tagged_mass") hMax=5.6;
     if(variables[i]=="kstTrk1Pt") hMin=0;
     if(variables[i]=="kstTrk1Pt") hMax=30;
     if(variables[i]=="kstTrk2Pt") hMin=0;
     if(variables[i]=="kstTrk2Pt") hMax=20;
     if(variables[i]=="mu1Pt"|| variables[i]=="mu2Pt") hMin=0;
     if(variables[i]=="mu1Pt"|| variables[i]=="mu2Pt") hMax=30;
     if(variables[i]=="kstTrk1Eta"|| variables[i]=="kstTrk2Eta") hMin=-2.8;
     if(variables[i]=="kstTrk1Eta"|| variables[i]=="kstTrk2Eta") hMax= 2.8;
     if(variables[i]=="sum_isopt_04") hMin=0;
     if(variables[i]=="sum_isopt_04") hMax=10;
     if(variables[i]=="bLBS")  hMin=0.;
     if(variables[i]=="bLBS")  hMax=2.;
     if(variables[i]=="bLBSE") hMin=0.;
     if(variables[i]=="bLBSE") hMax=0.02;
     if(variables[i]=="bDCABS") hMin=-0.01;
     if(variables[i]=="bDCABS") hMax=0.01;
     if(variables[i]=="bDCABSE") hMin=0.;
     if(variables[i]=="bDCABSE") hMax=0.004;
     if(variables[i]=="bCosAlphaBS") hMin=0.999;
     if(variables[i]=="bCosAlphaBS") hMax=1.;
     if(variables[i]=="bVtxCL") hMin=0.;
     if(variables[i]=="bVtxCL") hMax=1.;
     if(variables[i]=="kstTrk1DCABS") hMin=-0.4;
     if(variables[i]=="kstTrk1DCABS") hMax= 0.4;
     if(variables[i]=="kstTrk2DCABS") hMin=-0.4;
     if(variables[i]=="kstTrk2DCABS") hMax= 0.4;
     if(variables[i]=="kstTrk1DCABSE") hMin= 0.;
     if(variables[i]=="kstTrk1DCABSE") hMax= 0.02;
     if(variables[i]=="kstTrk2DCABSE") hMin= 0.;
     if(variables[i]=="kstTrk2DCABSE") hMax= 0.02;
     
     if(variables[i]==features[i]){
        std::cout<<Form("%s is a DNN features",variables[i].c_str())<<std::endl;
        isFeature ="(DNN feature)";
     }	
     
     if(iFeatDiv[i]>int(i)) {
       std::cout<<Form("%s/%s is a DNN features",variables[i].c_str(),variables[iFeatDiv[i]].c_str())<<std::endl;
       isFeature ="(DNN feature)";
     }  
     std::cout<<Form("iFeatDiv[%i]=%i",i,iFeatDiv[i])<<std::endl;
     isFeature="";
     std::cout<<Form("Booking Histograms for %s",variables[i].c_str())<<std::endl;
     
//   float hMin = 0;
//   float hMax = 35;
     HistData_Cut.push_back(new TH1D(  Form("Hx%s_Data_Cut" ,variables[i].c_str()), Form("%s Data %s %s DNN CUT"	 ,variables[i].c_str(),year.c_str(),isFeature.c_str()),100,hMin,hMax));
     HistData.push_back(new TH1D(  Form("Hx%s_Data" ,variables[i].c_str()), Form("%s Data %s %s"	 ,variables[i].c_str(),year.c_str(),isFeature.c_str()),100,hMin,hMax));
     HistMC.push_back(  new TH1D(  Form("Hx%s_MC"   ,variables[i].c_str()), Form("%s MC %s %s" 	         ,variables[i].c_str(),year.c_str(),isFeature.c_str()),100,hMin,hMax));
     HistMC_Cut.push_back(  new TH1D(  Form("Hx%s_MC_Cut"   ,variables[i].c_str()), Form("%s MC %s %s" 	         ,variables[i].c_str(),year.c_str(),isFeature.c_str()),100,hMin,hMax));
//     HistMCW.push_back( new TH1D(  Form("Hx%s_MCW"  ,variables[i].c_str()), Form("%s MC %s Reweighted %s",variables[i].c_str(),year.c_str(),isFeature.c_str()),100,hMin,hMax));
   }
//   fMC->cd();
   TreeMC->SetBranchAddress("eventN"	            ,&eventN );
//   TreeMC->SetBranchAddress("pass_preselection"     ,&pass_preselection );
//   TreeMC->SetBranchAddress("tagged_mass"	    ,&tagged_mass );
   TreeMC->SetBranchAddress("tagB0"  	            ,&tagB0 );
   TreeMC->SetBranchAddress("genSignal"  	    ,&genSignal );
   TreeMC->SetBranchAddress("mumuMass"  	    ,&mumuMass );
   TreeMC->SetBranchAddress("mumuMassE" 	    ,&mumuMassE );
   TreeMC->SetBranchAddress("truthMatchMum"	    ,&truthMatchMum );
   TreeMC->SetBranchAddress("truthMatchMup"	    ,&truthMatchMup );
   TreeMC->SetBranchAddress("truthMatchTrkm"	    ,&truthMatchTrkm );
   TreeMC->SetBranchAddress("truthMatchTrkp"	    ,&truthMatchTrkp );
   TreeMC->SetBranchAddress("trig"	            ,&trig );
   TreeMC->SetBranchAddress("weight"	            ,&weight );
   if(save) TreeMC->AddFriend(toutMC);

//   TreeMC->SetBranchAddress("xcut"	            ,&xcut );
//   TreeMC->SetBranchAddress("passB0Psi_jpsi"	    ,&passB0Psi_jpsi );
// add friend
//   TreeMC->AddFriend(TreeModel);
// add friend
//   float val_mva_B=0;
//   int shift = variables.size()-features.size();

   Long64_t nentriesMC = (Long64_t)TreeMC->GetEntries(); 
   std::cout<<Form("Found MC entries= = %lld",nentriesMC)<<std::endl;
   Long64_t MaxEveMC = (Long64_t)std::min(CutEve,(Long64_t)nentriesMC);
   std::cout<<Form("Reading num. MC entries= = %lld",MaxEveMC)<<std::endl;
//   for (Long64_t i=0;i<10000;i++) {
   for (Long64_t i=0;i<nentriesMC;i++) {
//   for (Long64_t i=0;i<nentriesMC;i++) {
//    for (Long64_t i=0;i<10000;i++) {
//     for (Long64_t i=0;i<MaxEveMC;i++) {
//     for (Long64_t i=0;i<nentriesMC;i++) {
       ScoreDNN = -1.;
       TreeMC->GetEntry(i);
       if ( i%100000==0 ) {
          std::cout<<Form("Event %lld",i)<<std::endl;
       }
       if( !save ){
       
        if(!(Cut_truth_match &&  Cut_ct && trig == 0 && (mumuMass < 2.702))) continue;
       }	
        bool skip = false;
	int Indx =0;
 	for (unsigned int j=0;j<=variables.size()-1;j++) {
	
	 if(Indx>int(VarF.size())) {
	  std::cout<<Form("This should never occur!! Indx=%i>VarF.size=%i",int(VarF.size()),Indx)<<std::endl;
	 }
//	   if(fabs(float(VarD[j]))<TMath::Exp(-60)) std::cout<<Form("warning in MC  %s=%f<10^-60",features[j].c_str(),VarD[j])<<std::endl;
 	 if(TMath::IsNaN(VarD[j])) {
//	  std::cout<<Form("warning in MC  %s is NaN",features[j].c_str())<<std::endl;
 	  skip = true;
 	  break;
 	 }
	 int iTest=iFeatDiv[j];
 	 if(iTest!=-1){
 	  if(iTest==int(j)){
 	   VarF[Indx]=float(VarD[j]);
	   Indx++;
	  }else{
	   if(TMath::IsNaN(VarD[iTest])){
   	    skip = true;
 	    break;
	   } 
 	   VarF[Indx]=float(VarD[j]/VarD[iTest]);
	   Indx++;
  	  }
 	 } 
 	}
	if(!skip){
 	 ScoreDNN = reader->EvaluateMVA( "DNN"    );
    	 HxDNNMC->Fill(ScoreDNN);
  	 for (unsigned int j=0;j<=variables.size()-1;j++) {
	  HistMC[j] ->Fill(VarD[j],weight);
	  if(ScoreDNN>dnncut) HistMC_Cut[j] ->Fill(VarD[j],weight);
  	 }
 	} 
       if(save) toutMC->Fill();
//       if(save) bMCScore->Fill();
     } 
     if(save){
      std::cout<<Form("Saving the MC   result in %s",foutMC->GetName())<<std::endl;
      foutMC->cd();
      toutMC->Write();
//      TreeMC->Write("", TObject::kOverwrite);
      foutMC->Close();
      delete foutMC;
     } 

 // //=============== Data TTree ================== 
    TFile *fData   = new TFile(NameFileDatap1.c_str(),"READ");
    std::cout<<Form("Opening Data File :%s \n",NameFileDatap1.c_str())<<std::endl;
    fData->cd();
    TTree *TreeData     = (TTree*)fData->Get("ntuple");  
//    if (save)  bDataScore = TreeData->Branch("ScoreDNN", &ScoreDNN,"ScoreDNN/D");

    for (unsigned int i=0;i<=variables.size()-1;i++) { 
     TreeData->SetBranchAddress(variables[i].c_str()     ,&VarD[i] );
    }
//    TreeData->SetBranchAddress("pass_preselection" ,&pass_preselection );
//    TreeData->SetBranchAddress("tagged_mass"	   ,&tagged_mass );
//    TreeData->SetBranchAddress("nsig_sw"	   ,&nsig_sw );
    TreeData->SetBranchAddress("mumuMass"	   ,&mumuMass );
    TreeData->SetBranchAddress("mumuMassE"	   ,&mumuMassE );
    TreeData->SetBranchAddress("eventN"	           ,&eventN );   
    if(save) TreeData->AddFriend(toutData);

//    TreeData->SetBranchAddress("xcut"	           ,&xcut );
//    TreeData->SetBranchAddress("passB0Psi_jpsi"	   ,&passB0Psi_jpsi );
    Long64_t nentriesData = (Long64_t)TreeData->GetEntries();
    std::cout<<Form("Found Data entries= = %d",(int)nentriesData)<<std::endl;
    Long64_t MaxEveData = (Long64_t)std::min(CutEve,nentriesData);
    std::cout<<Form("Reading num. Data entries= = %lld",MaxEveData)<<std::endl;
    for (Long64_t i=0;i<nentriesData;i++) {
//    for (Long64_t i=0;i<MaxEveData;i++) {
//    for (Long64_t i=0;i<10000;i++) {
           ScoreDNN = -1.;
    	   TreeData->GetEntry(i);
 	   if ( i%100000==0 ) {
 	      std::cout<<Form("Event %lld",i)<<std::endl;
 	    }
            if(!save){
	     if(mumuMass > 2.702) continue;
//	     if( eventN%2==0 ) continue; // parity!!!!!!!!!!!!!!!!!!!!!!111
//   	     if(tagged_mass < 5.0 || tagged_mass > 5.6) continue;
//   	     if(mumuMass*mumuMass>8.68  && mumuMass*mumuMass < 10.09) continue;
//   	     if(mumuMass*mumuMass>12.86 && mumuMass*mumuMass < 14.18) continue;
	    }
	    bool skip = false;
	    int Indx =0;
 	    for (unsigned int j=0;j<=variables.size()-1;j++) {
	
	     if(Indx>int(VarF.size())) {
	      std::cout<<Form("This should never occur!! Indx=%i>VarF.size=%i",int(VarF.size()),Indx)<<std::endl;
	      exit(1);
	     }
//	       if(fabs(float(VarD[j]))<TMath::Exp(-60)) std::cout<<Form("warning in MC  %s=%f<10^-60",features[j].c_str(),VarD[j])<<std::endl;
 	     if(TMath::IsNaN(VarD[j])) {
//	      std::cout<<Form("warning in MC  %s is NaN",features[j].c_str())<<std::endl;
 	      skip = true;
 	      break;
 	     }
	     int iTest=iFeatDiv[j];
 	     if(iTest!=-1){
 	      if(iTest==int(j)){
 	       VarF[Indx]=float(VarD[j]);
	       Indx++;
	      }else{
	       if(TMath::IsNaN(VarD[iTest])){
 	    	skip = true;
 	    	break;
	       }
 	       VarF[Indx]=float(VarD[j]/VarD[iTest]);
	       Indx++;
 	      }
 	     }
 	    }
	    if(!skip){
 	     ScoreDNN = reader->EvaluateMVA( "DNN"    );
    	     HxDNNData->Fill(ScoreDNN);
 	     for (unsigned int j=0;j<=variables.size()-1;j++) {
 	      HistData[j]->Fill(VarD[j]);
	      if(ScoreDNN>dnncut) HistData_Cut[j] ->Fill(VarD[j]);
 	     }
     	    }
        if(save) toutData->Fill();
//        if(save) bDataScore->Fill();
    }
    if(save){ 
     std::cout<<Form("Saving the Data result in %s",foutData->GetName())<<std::endl;
     foutData->cd();
     toutData->Write();
//     TreeData->Write("", TObject::kOverwrite);
     foutData->Close();
    }else{
     std::cout<<Form("Saving the plots in %s",pout->GetName())<<std::endl;
     pout->cd();
     HxDNNData->Write();
     HxDNNMC->Write();
     for (unsigned int i=0;i<=variables.size()-1;i++) {
       HistData[i]->Write();
       HistMC[i]->Write();
       HistData_Cut[i]->Write();
       HistMC_Cut[i]->Write();
     }  
     pout->Close();
    } 

}
//===============================================================================================================
void replaceAll(std::string& str, const std::string& from, const std::string& to) {
    if(from.empty())
        return;
    size_t start_pos = 0;
    while((start_pos = str.find(from, start_pos)) != std::string::npos) {
        str.replace(start_pos, from.length(), to);
        start_pos += to.length(); // In case 'to' contains 'from', like replacing 'x' with 'yx'
    }
}
//===============================================================================================================

std::map<std::string, std::string> splitFeatures(std::string FEATURES){
   std::vector<std::string> split( char *str, char c = ' ');
   std::stringstream indata;
   std::map<std::string, std::string> mappa;
   std::string line;
   std::vector<std::string>vstring ={{"",""}} ;
   
//
   indata<<FEATURES;
   if(!indata) { // file couldn't be opened
   	std::cout <<" "<<FEATURES<< " Error:  can not be opened" << std::endl;
   	exit(1);
   }
   while(std::getline(indata, line)) {
	 line.erase(std::remove(line.begin(), line.end(), '\t'), line.end());
	 line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
	 line.erase(std::remove(line.begin(), line.end(), ' ' ), line.end());

 	 char *cstr = new char [line.size()+1];


 	 strcpy (cstr, line.c_str());
	 std::cout <<"stringa->"<< cstr << std::endl;
	 vstring = split(cstr,'/');
	 mappa.insert( std::pair<std::string,std::string>(vstring[0],vstring[1]) );
// 	 std::cout <<"vstring[0]"<< vstring[0] << std::endl;
// 	 std::cout <<"vstring[1]"<< vstring[1] << std::endl;
	 
    }
    std::cout<<"//////////////////////////////////////////////////////////////////////////////////////////////////////////////"<<std::endl;
    for (std::map<std::string,std::string>::iterator imap = mappa.begin();
    			       imap != mappa.end();
    			       ++imap)
    {
   	std::cout <<"mappa->"<< (*imap).first<<" = "<<(*imap).second << std::endl;
    }
    std::cout<<"/////////////////////////////////////////////////////////////////////////////////////////////////////////////"<<std::endl;
  return mappa ;
}   
//===============================================================================================================
//===============================================================================================================
std::vector<std::string> split( char *str, char c = ' ')
{
    std::vector<std::string> result;

    while(1)
    {
         char *begin = str;

        while(*str != c && *str)
                str++;

        result.push_back(std::string(begin, str));

        if(0 == *str++)
                break;
    }

//	 std::cout <<"result[0]"<< result[0] << std::endl;
	 if (result.size()<2) result.push_back(""); 
//	 std::cout <<"result[1]"<< result[1] << std::endl;
    return result;
}
//===============================================================================================================

