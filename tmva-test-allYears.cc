//g++ -O3 -std=c++17 complete.cc -o complete     `root-config --cflags --libs`     -lRooFitCore -lRooFit -lRooStats     -lTMVA -lTMVAGui     -lstdc++fs
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
#include <filesystem>
//TestPlot.cc
#include "TGraphErrors.h"
#include <TText.h>
#include <RooRealVar.h>
#include <RooConstVar.h>
#include <RooPlot.h>
#include <RooDataSet.h>
#include <RooDataHist.h>
#include <RooHistPdf.h>
#include <RooCrystalBall.h>
#include <RooExponential.h>
#include <RooExtendPdf.h>
#include <RooFitResult.h>

//


using namespace RooFit;

timeval startTime, stopTime, totalTime;
clock_t startCPU, stopCPU; 
tms startProc, stopProc; 
void tmva_test_dnn();
void tmva_evaluate_dnn(bool save=false);
void test_evaluate_dnn_steps(int ncut, float xxstart);
void testPlot(std::string year, int ncut, float xxstart);
void replaceAll(std::string& str, const std::string& from, const std::string& to) ;
std::map<std::string, std::string> splitFeatures(std::string FEATURES);

double mc_sigma = 0.033;
double mc_mass  = 5.280;
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
int ncut = 0;
float xxstart = 0.90 ;
//
//Long64_t CutEve =14000;
//double CutEffDATA = 0.014;
double CutEffMC   = 0.266;
//
   
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
std::string year_default="2024";
std::string NameFileMCp0 = \
"/gwdata/users/dini/BPH/Run3/flatNtuple/ntuple_flat_2024MC_MC_alladd_vars.root";//modello
std::string NameFileMCp1 = \
"/gwdata/users/dini/BPH/Run3/flatNtuple/ntuple_flat_2024MC_MC_alladd_vars.root";
std::string NameFileDatap0 = \
"/gwdata/users/dini/BPH/Run3/flatNtuple/ntuple_flat_2024C012356_2024E0346_v2_data_alladd_vars.root";//modello
std::string NameFileDatap1 = \
"/gwdata/users/dini/BPH/Run3/flatNtuple/ntuple_flat_2024C012356_2024E0346_v2_data_alladd_vars.root"; //più statistica
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
std::vector<double> cutDNN;	

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
    std::cout<<Form("Usage: %s [year=2016,2017,2018, 2024] {dnn/save/steps/plot} {ncut} {xxstart}",argv[0])<<std::endl;
    std::cout<<Form("Please, set the year (at least)")<<std::endl;
    std::cout<<Form("example: %s 2024 dnn; if you want to train a new dnn\n",argv[0])<<std::endl;
    std::cout<<Form("         %s 2024 save ; if you want just produce the  score\n",argv[0])<<std::endl;
    std::cout<<Form("         %s 2024 steps ncut xxstart ; if you want to produce the file to plot\n",argv[0])<<std::endl;
    std::cout<<Form("         %s 2024 plot ncut xxstart ; if you want to produce the plot\n",argv[0])<<std::endl;
    exit(0);
  }else{
   if ( strcmp(argv[1],"2016") == 0 || \
        strcmp(argv[1],"2017") == 0 || \
	strcmp(argv[1],"2018") == 0 || \
	strcmp(argv[1],"2022") == 0 || \
	strcmp(argv[1],"2023") == 0 || \
	strcmp(argv[1],"2024") == 0 || \
	strcmp(argv[1],"2025") == 0      ){
      year=argv[1];
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
      //std::cout<<Form("FIRST CHECK")<<std::endl;
    }else{
     std::cout<<Form("not recognize year=%s",year.c_str())<<std::endl;
     exit(0);
    }
	
  }  
   
  if( argc > 2 ){
    
      if ((strcmp(argv[2],"dnn") == 0)){
      std::cout<<Form("Start training the DNN \n")<<std::endl;
      tmva_test_dnn();
     //tmva_evaluate_dnn(false);
      //std::cout<<Form("SECOND CHECK")<<std::endl;
      }else if ((strcmp(argv[2],"save") == 0)){
      std::cout<<Form("Save weights \n")<<std::endl;
      //std::cout<<Form("THIRD CHECK")<<std::endl;
      tmva_evaluate_dnn(true);
      }else if((strcmp(argv[2], "steps") == 0) && (argc>3)){
		   ncut = atoi(argv[3]);
		   std::cout<<Form("Set number of cut: %d",ncut)<<std::endl;
         //std::cout<<Form("FOURTH CHECK")<<std::endl;	
			if (argc > 4){
			xxstart = atof(argv[4]);
         //std::cout<<Form("FIFTH CHECK")<<std::endl;
			std::cout<<Form("YOU-Set xxstart: %f", xxstart)<<std::endl;
			std::cout<<Form("Making histo with cut")<<std::endl;
                        test_evaluate_dnn_steps(ncut, xxstart);
			std::cout<<Form("%s", year.c_str())<<std::endl;
			}else{
         //std::cout<<Form("SIXTH CHECK")<<std::endl;
			std::cout<<Form("AUTO-set xxstart: %f", xxstart)<<std::endl; 
			std::cout<<Form("Making histo with cut")<<std::endl;
			test_evaluate_dnn_steps(ncut, xxstart);
         }
		}else if((strcmp(argv[2], "plot") == 0) && (argc>3)){
                        ncut = atoi(argv[3]);
                        std::cout<<Form("Set number of cut: %d",ncut)<<std::endl;
         //std::cout<<Form("FOURTH CHECK")<<std::endl;	
		      if (argc > 4){
		       xxstart = atof(argv[4]);
         //std::cout<<Form("FIFTH CHECK")<<std::endl;
		       std::cout<<Form("YOU-Set xxstart: %f", xxstart)<<std::endl;
		       std::cout<<Form("%s", year.c_str())<<std::endl;
		       std::cout<<"STARTING PLOT"<<std::endl;
 		       testPlot(year, ncut, xxstart);
	              }else{
          //std::cout<<Form("SIXTH CHECK")<<std::endl;
		       std::cout<<Form("AUTO-set xxstart: %f", xxstart)<<std::endl;
		       std::cout<<"STARTING PLOT"<<std::endl;
		       testPlot(year, ncut, xxstart);
         }
      } 
   }  

   else{
    std::cout<<Form("Option not found! Exit... \n")<<std::endl;
    exit(0);
   }

}
// Here the DNN stuff 
void tmva_test_dnn(){
  gROOT ->Reset();
  gROOT->SetStyle("Plain");
  ROOT::EnableImplicitMT(num_threads); 
  TMVA::Tools::Instance();


  TFile *fData = new TFile(NameFileDatap0.c_str(),"READ");
  TTree *TreeData     = (TTree*)fData->Get("ntuple");
  std::cout<<Form("Open Data File :%s",NameFileDatap0.c_str())<<std::endl;
  int nentriesData = (int)TreeData->GetEntries();
  std::cout<<Form("Found Data entries= = %d\n",nentriesData)<<std::endl;
  TFile *fMC   = new TFile(NameFileMCp0.c_str(),"READ");
  std::cout<<Form("Open MC File :%s",NameFileMCp0.c_str())<<std::endl;
  TTree *TreeMC     = (TTree*)fMC->Get("ntuple");
  int nentriesMC = (int)TreeMC->GetEntries();
  std::cout<<Form("Found MC entries= = %d\n",nentriesMC)<<std::endl;
  DIR* dir = opendir(datasetYear.c_str());
  
  if(!dir) {
   gSystem->Exec(Form("mkdir -p %s",datasetYear.c_str()));
   std::cout<<Form("Create Dir  %s",datasetYear.c_str())<<std::endl;
  } 
  auto outputFile = TFile::Open(Form("%s/%s",datasetYear.c_str(),OutputFileName.c_str()), "RECREATE");

  TMVA::Factory *factory    =  new TMVA::Factory( "TMVAClassification", outputFile,"V:!Silent:Color:DrawProgressBar:Transformations=None:AnalysisType=Auto");
  TMVA::DataLoader * loader = new TMVA::DataLoader(datasetYear.c_str());
  Double_t signalWeight     = 1.0;
  Double_t backgroundWeight = 1.0;
  
  double mc_sigma = 0.04;
  double mc_mass  = 5.280;
  // You can add an arbitrary number of signal or background trees
  //Long64_t CutEveMC=(Long64_t)CutEve/CutEffMC;
  //Long64_t CutEveDATA=(Long64_t)CutEve/CutEffDATA;
  //std::cout<<Form("mycuteveMC	eventN<%lld",CutEveMC) <<std::endl;
  //std::cout<<Form("mycuteveDATA eventN<%lld",CutEveDATA) <<std::endl;
  //TCut mycuteveMC   = Form("eventN<%lld",CutEveMC);
  //TCut mycuteveDATA = Form("eventN<%lld",CutEveDATA);

  Cut_sig_mass     = Form("(tagged_mass > %f - 2.5*%f) && (tagged_mass < %f + 2.5*%f)",mc_mass,mc_sigma,mc_mass,mc_sigma); //perchè 2.5 e non 3 ? 

  Cut_ct	   = "( ( (tagB0==1) && (genSignal==1)) || ( (tagB0==0) && (genSignal==2) ) )";

  Cut_bkg_mass     = Form("( ( (tagged_mass > %f-7.*%f ) && (tagged_mass < %f-3.*%f ) ) || ( (tagged_mass > %f+3.*%f ) && (tagged_mass < %f+7*%f ) ) )",mc_mass,mc_sigma,mc_mass,mc_sigma,mc_mass,mc_sigma,mc_mass,mc_sigma );

  Cut_truth_match  = "( (truthMatchMum==1) && (truthMatchMup ==1) && (truthMatchTrkm==1) && (truthMatchTrkp==1) )";

  //double mc_mass_limit_i = mc_mass - 2.5 * mc_sigma;
  //double mc_mass_limit_u = mc_mass + 2.5 * mc_sigma;

  mycuts0	   = Cut_sig_mass  && "( (truthMatchMum==1) && (truthMatchMup ==1) && (truthMatchTrkm==1) && (truthMatchTrkp==1) ) &&  ( ( (tagB0==1) && (genSignal==1)) || ( (tagB0==0) && (genSignal==2) ) ) && trig == 0 && ((mumuMass < 2.702) || (mumuMass > 4.00) )"; //per MC 

  mycutb0	   = Cut_bkg_mass &&  "((mumuMass < 2.702) || (mumuMass > 4.00) )";

  TCut RMNone="";
  for (unsigned int i=0;i<=features.size()-1;i++) {
   RemoveNone.push_back(TCut(Form("!TMath::IsNaN(%s)",features[i].c_str())));
   if(i==0) RMNone=RemoveNone[i];
   if(i>0)  RMNone=RMNone&&RemoveNone[i];
  }
  TCut mycutNEve="int(eventN)<2000000";
  mycuts=mycuts0&&RMNone;
  mycutb=mycutb0&&RMNone;
  TTree *TreeData_Cut = TreeData->CopyTree(mycutb);
  int nentriesData_Cut= TreeData_Cut->GetEntries();

  std::cout<<Form("Warning, analyze only %i events!!!!\n",nentriesData_Cut)<<std::endl;
  std::cout<<Form("Warning, analyze only %i events!!!!\n",nentriesData_Cut)<<std::endl;
  std::cout<<Form("Warning, analyze only %i events!!!!\n",nentriesData_Cut)<<std::endl;
  std::cout<<Form("Warning, analyze only %i events!!!!\n",nentriesData_Cut)<<std::endl;
  std::cout<<Form("Warning, analyze only %i events!!!!\n",nentriesData_Cut)<<std::endl;

  Long64_t CutEveMC=(Long64_t)nentriesData_Cut/CutEffMC;
  TTree *TreeMC_Cut   = TreeMC  ->CopyTree(mycuts,  "" , CutEveMC, 0);
  int nentriesMC_Cut  = TreeMC_Cut->GetEntries();
  std::cout<<Form("Found MC   entries after Cuts = %d",nentriesMC_Cut)<<std::endl;
  std::cout<<Form("Found data entries after Cuts = %d",nentriesData_Cut)<<std::endl;
  std::cout<<Form("Cuts [MC]   = %s \n", mycuts.GetTitle ())<<std::endl;
  std::cout<<Form("Cuts [Data] = %s \n", mycutb.GetTitle ())<<std::endl;
  loader->AddTree( TreeMC_Cut    ,"Signal"	  , signalWeight    );
  loader->AddTree( TreeData_Cut  ,"Background",backgroundWeight );
  for (unsigned int i=0;i<=features.size()-1;i++) {
   loader->AddVariable( features[i].c_str(), 'F' );
  }
   std::cout<<"==== PrepareTrainingAndTestTree ===="<<std::endl;
   loader->PrepareTrainingAndTestTree( "", "",
                                      "SplitMode=Random:NormMode=NumEvents:V" );
////                                      "SplitMode=Random:NormMode=EqualNumEvents:V" );
//Boosted Decision Trees
   BookTXT =\
    Form("V:ErrorStrategy=CROSSENTROPY:VarTransform=N:WeightInitialization=XAVIERUNIFORM:Layout=TANH|128,TANH|128,TANH|128,LINEAR:TrainingStrategy=LearningRate=1e-5,Momentum=0.9,\
                                        ConvergenceSteps=100,BatchSize=128,TestRepetitions=7,\
                                        WeightDecay=1e-4,Regularization=None,\
                                        DropConfig=0.0+0.0+0.0+0.0:ValidationSize=0.20:Architecture=CPU");

//       TString inputLayoutString = "InputLayout=0|0|0"; 
//       TString batchLayoutString= "BatchLayout=0|0|0";
//       TString layoutString ("Layout=DENSE|64|TANH,BNORM|0.99|0.0001,DENSE|64|TANH,BNORM|0.99|0.0001,DENSE|64|TANH,BNORM|0.99|0.0001,DENSE|64|TANH,DENSE,BNORM|0.99|0.0001,DENSE|1|LINEAR");
//       // Define Training strategies 
//       // one can catenate several training strategies 
//       TString training1("LearningRate=1e-6,Momentum=0.9,Repetitions=1,"
//                         "ConvergenceSteps=30,BatchSize=128,TestRepetitions=1,"
//                         "MaxEpochs=100,WeightDecay=1e-4,Regularization=None,"
//                         "Optimizer=ADAM,DropConfig=0.0+0.0+0.0+0.:ValidationSize=0.50");
//       //     TString training2("LearningRate=1e-3,Momentum=0.9,Repetitions=1,"
//       //                       "ConvergenceSteps=10,BatchSize=128,TestRepetitions=1,"
//       //                       "MaxEpochs=20,WeightDecay=1e-4,Regularization=None,"
//       //                       "Optimizer=SGD,DropConfig=0.0+0.0+0.0+0.");
//   
//       TString trainingStrategyString ("TrainingStrategy=");
//       trainingStrategyString += training1; // + "|" + training2;
// 
//       // General Options.
// 
//       TString dnnOptions ("H:V:ErrorStrategy=CROSSENTROPY:VarTransform=G:"
//                           "WeightInitialization=XAVIER");
//       dnnOptions.Append (":"); dnnOptions.Append (inputLayoutString);
//       dnnOptions.Append (":"); dnnOptions.Append (batchLayoutString);
//       dnnOptions.Append (":"); dnnOptions.Append (layoutString);
//       dnnOptions.Append (":"); dnnOptions.Append (trainingStrategyString);
//       dnnOptions += ":Architecture=CPU:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10";
// //      dnnOptions += ":Architecture=GPU:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10";

//   std::cout<<dnnOptions<<std::endl;
   std::cout<<"==== BookMethod ===="<<std::endl;
   factory->BookMethod( loader,  TMVA::Types::kDL, "DNN_CPU",BookTXT.c_str() );



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
   TRoc->Print(Form("%s/roc-tmva-test-allYears-%s.pdf",datasetYear.c_str(),year.c_str()));
} 
//
//=================================================================================================================================================================
//   
// Here Evaluate DNN score x event
void tmva_evaluate_dnn(bool save){   //produco i file score, che vengono usati dallla funzione test_steps
 

   ROOT::EnableImplicitMT(num_threads);
   TMVA::Tools::Instance();

   TFile *foutData = 0;
   TFile *foutMC   = 0;
   TFile *pout	   = 0;
   TTree *toutMC   = 0;
   TTree *toutData = 0;
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
   TH1D* HxDNNMC   	 = new TH1D( "HxDNNMC"	, Form("DNN Output MC %s",year.c_str()) ,100, 0., 1.0);
   std::vector<double> VarD;			     
   std::vector<float> VarF;

//

   for (unsigned int i=0;i<=features.size()-1;i++) {
    VarD.push_back(0);
    VarF.push_back(0);
     std::cout <<features[i]<< std::endl;
     std::map<std::string,std::string> temap= splitFeatures(features[i]);
     mapFeat.insert(make_pair(i,temap));
     
   }

// DNN reader init   
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
  
   for (unsigned int i=0;i<=vartested.size()-1;i++) { 
    VarD.push_back(0);
    variables.push_back(vartested[i]);
    cstudies.push_back(new TCanvas(Form("c_%s",vartested[i].c_str()),Form("MC %s DNN studies %s (variable check)",vartested[i].c_str(),year.c_str()),200,10,900,780));
   }
   
   for (unsigned int i=0;i<=features.size()-1;i++) { 
    reader->AddVariable( features[i].c_str(), &VarF[i]);
//    cstudies.push_back(new TCanvas(Form("c_%s",features[i].c_str()),Form("MC %s reweighting studies %s",features[i].c_str(),year.c_str()),200,10,900,780));
   }

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

//
//=============== Read  MC TTree ================== 
//
   TFile *fMC   = new TFile(NameFileMCp1.c_str(),"READ");
   std::cout<<Form("Opening MC File :%s \n",NameFileMCp1.c_str())<<std::endl;
   TTree *TreeMC     = (TTree*)fMC->Get("ntuple");
   

   for (unsigned int i=0;i<=variables.size()-1;i++) { 
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
   }
//   fMC->cd();
   TreeMC->SetBranchAddress("eventN"	            ,&eventN );
   TreeMC->SetBranchAddress("pass_preselection"     ,&pass_preselection );
//   TreeMC->SetBranchAddress("tagged_mass"	  ,&tagged_mass );
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


   Long64_t nentriesMC = (Long64_t)TreeMC->GetEntries(); 
   std::cout<<Form("Found MC entries= = %lld",nentriesMC)<<std::endl;
   //Long64_t MaxEveMC = (Long64_t)std::min(CutEve,(Long64_t)nentriesMC);
   //std::cout<<Form("Reading num. MC entries= = %lld",MaxEveMC)<<std::endl;
   for (Long64_t i=0;i<nentriesMC;i++) {
       ScoreDNN = -1.;
       TreeMC->GetEntry(i);
       if ( i%100000==0 ) {
          std::cout<<Form("Event %lld",i)<<std::endl;
       }
      //  if( !save ){
      //  
      //   if(!(Cut_truth_match &&  Cut_ct && trig == 0)) continue;
      //   if(!((mumuMass < 2.702) || (mumuMass > 4.))) continue;
      //  }	
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
     } 
     if(save){
      std::cout<<Form("Saving the MC   result in %s",foutMC->GetName())<<std::endl;
      foutMC->cd();
      toutMC->Write();
      foutMC->Close();
      delete foutMC;
     } 
//
//=============== Read Data TTree ================== 
//
    TFile *fData   = new TFile(NameFileDatap1.c_str(),"READ");
    std::cout<<Form("Opening Data File :%s \n",NameFileDatap1.c_str())<<std::endl;
    fData->cd();
    TTree *TreeData     = (TTree*)fData->Get("ntuple");  
//    if (save)  bDataScore = TreeData->Branch("ScoreDNN", &ScoreDNN,"ScoreDNN/D");

    for (unsigned int i=0;i<=variables.size()-1;i++) { 
     TreeData->SetBranchAddress(variables[i].c_str()     ,&VarD[i] );
    }
    TreeData->SetBranchAddress("pass_preselection" ,&pass_preselection );
//    TreeData->SetBranchAddress("tagged_mass"	   ,&tagged_mass );
    TreeData->SetBranchAddress("mumuMass"	   ,&mumuMass );
    TreeData->SetBranchAddress("mumuMassE"	   ,&mumuMassE );
    TreeData->SetBranchAddress("eventN"	           ,&eventN );   
    if(save) TreeData->AddFriend(toutData);

    Long64_t nentriesData = (Long64_t)TreeData->GetEntries();
    std::cout<<Form("Found Data entries= = %d",(int)nentriesData)<<std::endl;
    //Long64_t MaxEveData = (Long64_t)std::min(CutEve,nentriesData);
    //std::cout<<Form("Reading num. Data entries= = %lld",MaxEveData)<<std::endl;
    for (Long64_t i=0;i<nentriesData;i++) {
      {
           ScoreDNN = -1.;
    	   TreeData->GetEntry(i);
 	   if ( i%100000==0 ) {
 	      std::cout<<Form("Event %lld",i)<<std::endl;
 	    }
      // if(!save){
	   //   if(!((mumuMass < 2.702) || (mumuMass > 4.))) continue;
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


void test_evaluate_dnn_steps(int ncut, float xxstart){
 
   ROOT::EnableImplicitMT(num_threads); 
   TMVA::Tools::Instance();

   float xxcut=0.;
   for (int  i=0;i<ncut;i++) {
    xxcut=xxstart+i*(1.-xxstart)/float(ncut);
    cutDNN.push_back( xxcut);
    HistMC_Cut.push_back  (  new TH1D(  Form("Hxmass_MC_Cut%i"   ,i), Form("Mass MC   ScoreDNN>%f",xxcut),70,5.0,5.6));
    HistData_Cut.push_back(  new TH1D(  Form("Hxmass_Data_Cut%i" ,i), Form("Mass Data ScoreDNN>%f",xxcut),70,5.0,5.6));
   }
   TFile *foutData = 0;
   TFile *foutMC   = 0;
   TFile *pout	   = 0;
   TTree *toutMC   = 0;
   TTree *toutData = 0;
   foutMC = new TFile(Form("%sMC_LMNR-ScoreDNN.root",year.c_str()),"READ");
   std::cout<<Form("Opening MC Score DNN File :%s \n",foutMC->GetName())<<std::endl;
   toutMC     = (TTree*)foutMC->Get("sDNN");
   if(toutMC==0) {
     std::cout<<Form("Ntupla sDNN MC not found!")<<std::endl;
    exit(0);
   } 
   toutMC->SetBranchAddress("eventN", &eventN);
   toutMC->SetBranchAddress("ScoreDNN", &ScoreDNN);
   foutData = new TFile(Form("%sData_LMNR-ScoreDNN.root",year.c_str()),"READ");
   std::cout<<Form("Opening Data Score DNN File :%s \n",foutData->GetName())<<std::endl;
   toutData     = (TTree*)foutData->Get("sDNN");
   if(toutData==0) {
     std::cout<<Form("Ntupla sDNN Data not found!")<<std::endl;
    exit(0);
   } 
   toutData->SetBranchAddress("eventN", &eventN);
   toutData->SetBranchAddress("ScoreDNN", &ScoreDNN);
   pout = new TFile(Form("test-%sDATA_LMNR-Plots-%dCUT.root",year.c_str(), ncut),"RECREATE"); //genero i file che poi devono essere plottati dall'alto programma testplot.cc
 
   gROOT ->Reset();
   gROOT->SetStyle("Plain");

   TH1D* HxDNNData   	 = new TH1D( "HxDNNData"	, Form("DNN Output Data %s",year.c_str()) ,100, 0., 1.0);
   TH1D* HxDNNMC   	 = new TH1D( "HxDNNMC"		, Form("DNN Output MC %s",year.c_str()) ,100, 0., 1.0);

   

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
   
   
   TreeMC->SetBranchAddress("eventN"	            ,&eventN );
   TreeMC->SetBranchAddress("pass_preselection"     ,&pass_preselection );
   TreeMC->SetBranchAddress("tagged_mass"	    ,&tagged_mass );
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
   TreeMC->AddFriend(toutMC);


   Long64_t nentriesMC = (Long64_t)TreeMC->GetEntries(); 
   std::cout<<Form("Found MC entries= = %lld",nentriesMC)<<std::endl;
//   for (Long64_t i=0;i<10000;i++) {
   for (Long64_t i=0;i<nentriesMC;i++) {
       TreeMC->GetEntry(i);
       if ( i%100000==0 ) std::cout<<Form("Event %lld",i)<<std::endl;
       if(tagged_mass<5.0||tagged_mass>5.6) continue;
       //if (!( ( (tagB0==1) && (genSignal==1)) || ( (tagB0==0) && (genSignal==2) )) ) continue; //tCut esplicitato
       if (!(trig == 0 && pass_preselection==1 )) continue;
       if (!((mumuMass < 2.702) || (mumuMass > 4.))) continue;
//       if(!(Cut_truth_match &&  Cut_ct && trig == 0 && (mumuMass < 2.702) && pass_preselection==1 )) continue;
       
	 HxDNNMC->Fill(ScoreDNN);
  	 for (unsigned int j=0;j<cutDNN.size()-1;j++) {
	  if(ScoreDNN>cutDNN[j]) HistMC_Cut[j] ->Fill(tagged_mass,weight);
  	 }
     } 
 // //=============== Data TTree ================== 
    TFile *fData   = new TFile(NameFileDatap1.c_str(),"READ");
    std::cout<<Form("Opening Data File :%s \n",NameFileDatap1.c_str())<<std::endl;
    fData->cd();
    TTree *TreeData     = (TTree*)fData->Get("ntuple");  

    TreeData->SetBranchAddress("pass_preselection" ,&pass_preselection );
    TreeData->SetBranchAddress("tagged_mass"	   ,&tagged_mass );
    TreeData->SetBranchAddress("mumuMass"	   ,&mumuMass );
    TreeData->SetBranchAddress("mumuMassE"	   ,&mumuMassE );
    TreeData->SetBranchAddress("eventN"	           ,&eventN );   
    TreeData->AddFriend(toutData);

    Long64_t nentriesData = (Long64_t)TreeData->GetEntries();
    std::cout<<Form("Found Data entries= = %d",(int)nentriesData)<<std::endl;
    //Long64_t MaxEveData = (Long64_t)std::min(CutEve,nentriesData);
    //std::cout<<Form("Reading num. Data entries= = %lld",MaxEveData)<<std::endl;
    for (Long64_t i=0;i<nentriesData;i++) {
           ScoreDNN = -1.;
    	   TreeData->GetEntry(i);
 	   if ( i%100000==0 ) std::cout<<Form("Event %lld",i)<<std::endl;
	   if(tagged_mass<5.0||tagged_mass>5.6) continue;
           if(pass_preselection!=1) continue;
           if(!((mumuMass < 2.702) || (mumuMass > 4.))) continue;
    	   HxDNNData->Fill(ScoreDNN);
  	   for (unsigned int j=0;j<=cutDNN.size()-1;j++) {
	    if(ScoreDNN>cutDNN[j]) HistData_Cut[j] ->Fill(tagged_mass);
	   } 
    } 
    std::cout<<Form("Saving the plots in %s",pout->GetName())<<std::endl;
    pout->cd();
    HxDNNData->Write();
    HxDNNMC->Write();
    for (unsigned int i=0;i<cutDNN.size()-1;i++) {
      HistData_Cut[i]->Write();
      HistMC_Cut[i]->Write();
    }
    pout->Close();

}

//=================================================================================================


void testPlot(const std::string year, int ncut, float xxstart) {


  std::vector<RooDataHist*> MC_Hist;
  std::vector<RooDataHist*> Data_Hist;
  std::vector<TCanvas*> cMCfitted;
  std::vector<TCanvas*> cDatafitted;

  // mean mass and sigma values taken from previous analysis work
  const double mc_sigma = 0.033;
  const double mc_mass  = 5.280;

  const double B0Mass_ = 5.27958; // nominal mass for bkg section

  RooRealVar mass ("mass", "mass", 5., 5.6); // declaring var
  mass.setRange("full", 5.0,5.6);
  // MC: declaring pars
  // RooRealVar mean ("mean", "mean of CryBall", 5.0, 5.6);
  // RooRealVar sigma("sigma", "CryBall_sigma", 0.01, 0.2);
  // RooRealVar aL("L_tail", "L_gauss_tail", 1., 2.5);
  // RooRealVar nL("L_norm", "L_norm", 0.9, 10.);
  // RooRealVar aR("R_tail", "R_gauss_tail", 1., 2.5);
  // RooRealVar nR("R_norm", "R_norm", 0.9, 10.);
  // RooRealVar mean ("mean", "mean of CryBall", B0Mass_, 5.15, 5.4);
  // RooRealVar sigma("sigma", "CryBall_sigma", 0.03, 0.01, 0.2);
  // RooRealVar aL("L_tail", "L_gauss_tail", 1.2, 1., 1.5);
  // RooRealVar nL("L_norm", "L_norm", 0.9, 10.);
  // RooRealVar aR("R_tail", "R_gauss_tail", 1.4, 1., 1.8);
  // RooRealVar nR("R_norm", "R_norm", 0.9, 15.);
 
  // RooRealVar mean ("mean", "mean of CryBall", B0Mass_,  5.15, 5.4);
  // RooRealVar sigma("sigma", "CryBall_sigma", 0.03, 0.01, 0.2);
  // RooRealVar aL("L_tail", "L_gauss_tail", 1.2, 1., 1.5);
  // RooRealVar nL("L_norm", "L_norm", 0.9, 10.);
  // RooRealVar aR("R_tail", "R_gauss_tail", 1.3, 1., 1.6);
  // RooRealVar nR("R_norm", "R_norm", 1., 15.);

  RooRealVar mean ("mean", "mean of CryBall", B0Mass_,  5.2, 5.35);
  RooRealVar sigma("sigma", "CryBall_sigma", 0.03, 0.001, 0.3);
  RooRealVar aL("L_tail", "L_gauss_tail", 1.2, 1., 2.5);
  RooRealVar nL("L_norm", "L_norm", 0.9, 10.);
  RooRealVar aR("R_tail", "R_gauss_tail", 1.3, 1., 2.5);
  RooRealVar nR("R_norm", "R_norm", 0.8, 20.);

  // Data: declaring pars
  RooRealVar alpha("alpha", "exp_alpha", -10., 0.);

  // setting exponential fit range ( x < mc_mass - 3mc_sigma || x> mc_mass + 3mc_sigma )
  mass.setRange("R1", 5.00, mc_mass - 3.*mc_sigma);
  mass.setRange("R2", mc_mass + 3.*mc_sigma, 5.6);



  //significance vars
  double n_signal;
  double n_background;
  double n_bkg;
  double fact;
  std::vector<double> signifDNN;
  std::vector<double> errSignifDNN;


  double index = 0.;
  double max_signi = 0.;
  double max_cut = 0.;
  double max_B = 0;
  double max_S = 0;

  // scaling vars

  //double corrFactorDataonMCJpsi = 1.3;
  double corrFactorDataonMCJpsi = 0.51; //proviene dal calcolo della Jpsi
  double lumiMC = 684.1*98.86/9.4/corrFactorDataonMCJpsi;

  double fraz_C = 6./8.;
  double fraz_E = 4./16.;
  double lumiDatasetC =  7.24;// tot lumi (1/fb) 2024C
  //double lumiDatasetD =  7.96 ;// tot lumi (1/fb) 2024D
  double lumiDatasetE =  11.32 ;// tot lumi (1/fb) 2024E
  //double lumiDatasetF =  27.76 ;// tot lumi (1/fb) 2024F
  double lumiData = fraz_C * lumiDatasetC + fraz_E * lumiDatasetE;

  double scale24 = lumiData / lumiMC ;
  //double scale24 = 1.;
  //double scale24 = 0.000580292;
  std::cout << "Scale factor 2024: " << scale24 << std::endl;


   std::string OutDir = Form("TestPlot_results_%d", ncut);
   
   if(!gSystem->AccessPathName((const char *)OutDir.c_str())){
    std::cout<<Form("Exist output Dir %s, ...save and make it again",OutDir.c_str() )<<std::endl;
    gSystem->Exec(Form("mv %s %s.`/usr/bin/date -Iminutes -r %s`",OutDir.c_str(),OutDir.c_str(),OutDir.c_str()));
    gSystem->Exec(Form("mkdir %s",OutDir.c_str()));
   }else{
    std::cout<<Form("output Dir %s doesn't exist, let's make it",OutDir.c_str() )<<std::endl;
    gSystem->Exec(Form("mkdir %s",OutDir.c_str()));
   }


  std::ofstream yield_evaluate;
  std::ofstream ratio_evaluate; //file with s/sqrt(s+b) from various cuts
  ratio_evaluate.open(Form("%s/ratio_results_%d.txt", OutDir.c_str(), ncut));
  yield_evaluate.open(Form("%s/yield_results_%d.txt", OutDir.c_str(), ncut));


// file will be closed at the end of loop

// loop on cuts 

//RooRealVar nsig("nsig"         , "signal frac"   ,   40000,     0,   1000000);
//RooRealVar nbkg("nbkg", "bkg n"  ,  1000.,       0.,   550000.);

   RooRealVar nsig("nsig", "signal yield MC"  ,    3500000,     0,   5000000);
//RooRealVar nbkg("nbkg", "bckg yield Data"  ,  8000.,     0.,   500000.);

   
   std::string NameFileInput = Form("test-%sDATA_LMNR-Plots-%dCUT.root", year.c_str(), ncut);
   TFile *FileInput =  TFile::Open(NameFileInput.c_str() ,"READ");
   if( !FileInput){
    std::cout<<Form("File %s doesn't exist",NameFileInput.c_str() )<<std::endl;
    std::cout<<Form("you should execute this command before: ./tmva-test-allYears %s steps %d %d ",year.c_str(),ncut,int(xxstart) )<<std::endl;
    exit(0);
   }
   std::vector<TH1D*> HistData_Cut;
   std::vector<TH1D*> HistMC_Cut;
   std::vector<double> cutDNN;		     
   std::vector<TCanvas*> cstudies;	
// summary canvas 
   TCanvas c1("c_allCuts","summary of the cut plots",200,10,1800,780);
   c1.Divide(2,1);


// create saving directory
// std::string folder_name = Form("TestPlot_results_%d", ncut);;
// if(std::filesystem::exists(folder_name)==1) {
//     
//     if(std::filesystem::exists(folder_name + ".save")==1){
//        std::cout<<"REMOVING PREVIOUS SAVING"<<std::endl;
//       std::filesystem::remove_all(folder_name + ".save");
//     }
// 
//     std::cout<<"FOUND FOULDER; MOVING IN .SAVE"<<std::endl;
//     
//     std::string new_name = folder_name + ".save";
//     std::filesystem::copy(folder_name, new_name);
//     std::filesystem::remove_all(folder_name);
//     std::filesystem::create_directory(folder_name);
// }
// else {
//     std::cout<<"NOT FOUND FOULDER; CREATE NEW ONE"<<std::endl;
//     std::filesystem::create_directory(folder_name);
// }

   gStyle->SetOptStat(000000);
// loop on cuts
   float xxcut=0.;
   for (int  i=0;i<ncut-1;i++) {
    xxcut=xxstart+i*(1.-xxstart)/float(ncut);
    cutDNN.push_back( xxcut);
    std::cout<<Form("Loop Histo cut scoreDNN>%f",xxcut)<<std::endl;
//  cut canvas 
    cstudies.push_back(new TCanvas(Form("c_cut_%d",i),Form("ScoreDNN>%f",xxcut),200,10,1800,780));
    cstudies[i]->Divide (2,1);
//  histograms retrive 
    HistData_Cut.push_back((TH1D*)FileInput->Get(Form("Hxmass_Data_Cut%i"   ,i)));
    HistMC_Cut.push_back((TH1D*)FileInput->Get(Form("Hxmass_MC_Cut%i"    ,i)));
//  plot cuts' histos
    cstudies[i]->cd(1);
    HistMC_Cut[i]->Draw();
    cstudies[i]->cd(2);
    HistData_Cut[i]->Draw();
//  save cuts' histos
    cstudies[i]->Print(Form("%s/testPlot-cut%d.pdf", OutDir.c_str(), i));
// summary histos
    c1.cd(1);
    HistMC_Cut[i]->SetLineColor(2+i);
    HistMC_Cut[i]->Draw("same");
    c1.cd(2);
    HistData_Cut[i]->SetLineColor(1+i);
    HistData_Cut[i]->Draw("same");
   }
//  save summary histos
   c1.Print(Form("%s/testPlot-allCuts.pdf", OutDir.c_str()));
   std::cout<<"End of the loop"<<std::endl;
   
// Here fitting the two pdfs, using Roofit 

// MC   --> CrystalBall pdf with symmetrical gaussian core & asymmetrical tails 
//    Data --> Exponential pdf   

for (int i = 0; i < ncut-1; i++) {
	xxcut=xxstart+i*(1.-xxstart)/float(ncut);
	
	// MC area 
	MC_Hist.push_back(new RooDataHist(Form("MC_hist_%d", i), Form("ScoreDNN>%f",xxcut),mass, HistMC_Cut[i])); // loading the dataset
	cMCfitted.push_back(new TCanvas(Form("cMCfit_cut_%d",i),Form("ScoreDNN>%f",xxcut),200,10,1000,780));
	MC_Hist[i]->Draw();
	RooCrystalBall cryball("cryball", "crystalball pdf", mass, mean, sigma, aL, nL, aR, nR); // creating the Crystal Ball pdf
        RooExtendPdf signal("signal","signal"    ,  cryball, nsig, "full");
	
	std::cout<<Form("Here fit results from MC for %d cut\n", i)<<std::endl;
	//RooFitResult* mc_r = signal.fitTo(*MC_Hist[i],Extended(true),SumW2Error(true), Save()); //fitting and saving MC results
        RooFitResult* mc_r = signal.fitTo(*MC_Hist[i],Extended(true),Strategy(1), MaxCalls(50000), Minimizer("Minuit"),Save());
	mc_r->Print(); //printing MC fit results
	std::cout<<Form("End of fit results from MC for %d cut\n", i)<<std::endl;
	
   //integral of the signal
//        RooAbsReal* SignalRangeMassS = cryball.createIntegral(mass, mass, "signRangeMC" );
//        RooAbsReal* SignalFullMassS = cryball.createIntegral(mass, mass, "full" );

	//drawing the MC plots
	RooPlot* frame_mc = mass.frame(); 
	MC_Hist[i]->plotOn(frame_mc);

	frame_mc->GetYaxis()->SetTitleOffset(1.55);
        frame_mc->GetXaxis()->SetTitle("#mu^{+}#mu^{-} K^{#pm} #pi{ #mp} mass (GeV/c^{2})");
        frame_mc->GetYaxis()->SetTitle(Form("Events/%4.4f (GeV/c^{2})",HistMC_Cut[i]->GetBinWidth(1)));
	signal.plotOn(frame_mc);
	signal.paramOn(frame_mc, Layout(0.55, 0.99));
	//MC_Hist.statOn(frame_mc, Layout(0.01, 0.40));
	// printing chisquare
//	TLatex Tl;
//	Tl.SetTextSize(0.04);
	double chiSq = frame_mc->chiSquare(7);
	
        n_signal = nsig.getVal()*scale24;
        double err_signal = nsig.getError() * scale24;

	frame_mc->Draw();
	frame_mc->SetTitle(Form("MassMC_ScoreDNN>%f chi^2 = %f ", xxcut, chiSq));
	cMCfitted[i]->Update();
//	Tl.DrawLatexNDC(0.2, 0.86, Form("#chi^{2}/N_{bins} = %.2f", chiSq));
	cMCfitted[i]->Print(Form("%s/MC_fitted_result_%d.pdf",OutDir.c_str(), i), "pdf");
	

	// Data area 
	Data_Hist.push_back(new RooDataHist(Form("Data_hist_%d", i), Form("ScoreDNN>%f",xxcut),mass, HistData_Cut[i]));  	
	cDatafitted.push_back(new TCanvas(Form("cDatafit_cut_%d",i),Form("ScoreDNN>%f",xxcut),200,10,1000,780));
    	Data_Hist[i]->Draw();
	
	mass.setRange("signRangeMC", B0Mass_ - 3.*sigma.getVal(), B0Mass_ + 3.*sigma.getVal() );
     //   RooRealVar nbkg("nbkg", "bkg n"  ,  1000.,       0.,   550000.);
        RooRealVar nbkg("nbkg", "bkg n"  ,  5000.,       0.,   500000.);
	RooExponential expon("exponential", "exponential pdf", mass, alpha); // creating the exponential pdf
        expon.setNormRange("R1,R2");
        RooExtendPdf bkg("bkg","ebkg"  ,  expon, nbkg, "R1,R2");
        
	std::cout<<Form("Here fit results from data for %d cut\n", i)<<std::endl;
	//RooFitResult* d_r = bkg.fitTo(*Data_Hist[i],Extended(true), Range("R1,R2"), NormRange("full"), Save(), PrintLevel(-1)); //fitting and saving data results 
         RooFitResult* d_r = bkg.fitTo(*Data_Hist[i],Extended(true),Strategy(2), MaxCalls(50000), Range("R1,R2"), Save(), PrintLevel(-1));
//	RooFitResult* d_r = expon.fitTo(*Data_Hist[i], Range("R1,R2"), NormRange("R1,R2"), Save()); //fitting and saving data results 
	d_r->Print(); //printing data fit results
	std::cout<<Form("End of fit results from data for %d cut\n", i)<<std::endl;


	mass.setRange("signRangeMC", B0Mass_ - 3.*sigma.getVal(), B0Mass_ + 3.*sigma.getVal() ); //redef mass in dep from signal itself
	RooExponential expon1("expon1", "exponential pdf", mass, alpha); // creating the exponential pdf
        expon1.setNormRange("full");
	RooAbsReal* BckgRangeMassS = expon1.createIntegral(mass, mass, "signRangeMC"); // integral of expon under the signal
	RooAbsReal* BckgFullMassS = expon1.createIntegral(mass, mass, "full");
	RooAbsReal* BckgRight= expon1.createIntegral(mass, mass, "R1");
	RooAbsReal* BckgLeft = expon1.createIntegral(mass, mass, "R2");
	

//	Data_Hist[i]->GetXAxis()->SetTitle("#mu^{+}#mu^{-}K#pi mass (GeV/c^{2})");
//	Data_Hist[i]->GetYAxis()->SetTitle(Form("Events/%4.4f (GeV/c^{2})",Data_Hist[i]->GetBinWidth));
	


        fact = BckgFullMassS->getVal() / (BckgRight->getVal() + BckgLeft->getVal()) ; //rescaling del background dalle side band al full range (5.0-5.6)
        n_bkg = nbkg.getVal() * fact;

	//drawing the data plots
	RooPlot* frame_d = mass.frame(); 
	Data_Hist[i]->plotOn(frame_d);
        frame_d->GetXaxis()->SetTitle("#mu^{+}#mu^{-}K#pi mass (GeV/c^{2})");
        frame_d->GetYaxis()->SetTitle(Form("Events/%4.4f (GeV/c^{2})",HistData_Cut[i]->GetBinWidth(1)));
	//expon.plotOn(frame_d,Range("R1,R2"),NormRange("R1,R2"));
   //expon.plotOn(frame_d,Range("R1,R2"),Normalization(nbkg.getVal(), RooAbsReal::NumEvent));
        bkg.plotOn(frame_d,Range("R1,R2"),NormRange("R1,R2"));
        chiSq = frame_d->chiSquare(1);
   //expon.plotOn(frame_d,Range("full"),Normalization(n_bkg, RooAbsReal::NumEvent), LineStyle(7));
	bkg.plotOn(frame_d,Range("full"), LineStyle(7));
   //expon.paramOn(frame_d, Layout(0.55, 0.99));
        bkg.paramOn(frame_d, Layout(0.55, 0.99));
	//data_Hist.statOn(frame_mc, Layout(0.01, 0.40));
	
	
   
	frame_d->Draw();
        frame_d->SetTitle(Form("MassData_ScoreDNN>%f chiSq = %f", xxcut, chiSq));
   //	Tl.DrawLatexNDC(0.2, 0.86, Form("#chi^{2}/N_{bins} = %.2f", chiSq));
	cMCfitted[i]->Print(Form("%s/MC_fitted_result_%d.pdf",  OutDir.c_str(), i), "pdf");
	cDatafitted[i]->Update();
	cDatafitted[i]->Print(Form("%s/Data_fitted_result_%d.pdf", OutDir.c_str(), i), "pdf");
 	
	
	// evaluating the significance and printing on file
   //n_signal = nsig.getVal()*scale24; //scaling factor added to compensate differences in Luminosity
   //n_signal = nsig.getVal();
   //double err_signal = nsig.getError() * scale24;
   //n_signal = nsig.getVal()*scale24;

	//std::cout << Form("nsig cut %d = %lf\n", i, n_signal)<<std::endl;
	//n_background = ( BckgFullMassS->getVal()*nbkg.getVal() ) / (BckgRight->getVal() + BckgLeft->getVal());

   
   
   
      n_background = BckgRangeMassS->getVal() * n_bkg; 
//    n_background = BckgFullMassS->getVal() * n_bkg;
      double err_background = sqrt( pow( nbkg.getVal() * BckgFullMassS->getPropagatedError(*d_r) ,2.) + pow( nbkg.getError() * BckgFullMassS->getVal() ,2.));
	//n_background = nbkg.getVal();
	//std::cout << Form("nbkg cut %d = %lf\n", i, n_background)<<std::endl;
      double dS=((sqrt(n_signal + n_background)- n_signal / (2. * sqrt(n_signal + n_background)))/(n_signal + n_background));
      double dB=(n_signal/(2*(pow((n_signal + n_background),(3./2.)))));
      double err = sqrt((pow((dS*err_signal),2.))+(pow((dB*err_background),2.)));

      signifDNN.push_back(n_signal / (sqrt(n_signal + n_background))); 
   //signifDNN.push_back(n_signal / (sqrt(n_signal + n_background))); 
      errSignifDNN.push_back(err);
      ratio_evaluate << xxcut << " "<< signifDNN[i] <<std::endl;

      if (signifDNN[i]>max_signi){

         index = i;
         max_signi = signifDNN[i];
         max_cut = xxcut;
         max_S = n_signal;
         max_B = n_background;

      }


//      ratio_evaluate << xxcut << " " << signifDNN[i]<<" 0 "<< errSignifDNN[i] << std::endl; 
      yield_evaluate << xxcut << " " << n_signal << " " << n_background << std::endl;


}

std::cout<<"MAX info: " << index << " " << max_cut << " " << max_signi << " " << max_S << " " << max_B << std::endl;

ratio_evaluate.close();
yield_evaluate.close();


}


