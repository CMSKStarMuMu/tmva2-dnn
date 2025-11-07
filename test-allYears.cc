//g++ -O3 -o test-allYears test-allYears.cc  `root-config --cflags --libs` -lRooFit -lGenVector -lMathCore -lTMVA -l TMVAGui 
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
void test_evaluate_dnn(int ncut, float xxstart);
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
std::vector<TCut> RemoveNone;
std::vector<TH1D*> HistData;
std::vector<TH1D*> HistMC;
std::vector<TH1D*> HistData_Cut;
std::vector<TH1D*> HistMC_Cut;
std::vector<TH1D*> HistMCW;
std::vector<TCanvas*> cstudies;	
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

int ncut=10;
float xxstart=0.90;
			     
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
      std::cout<<Form("Setting year=%s",year.c_str())<<std::endl;
      std::cout<<Form("Setting dataset year=%s",datasetYear.c_str())<<std::endl;
    }else{
     std::cout<<Form("not recognize year=%s",year.c_str())<<std::endl;
     exit(0);
    }
	
  }  
   
  if( argc>2 ){
   if ((strcmp(argv[2],"dnn") == 0)){
    std::cout<<Form("Start training the DNN \n")<<std::endl;
    tmva_test_dnn();
    test_evaluate_dnn(ncut, xxstart);
   }
   if ((strcmp(argv[2],"save") == 0)){
    std::cout<<Form("Save weights \n")<<std::endl;
    test_evaluate_dnn(ncut, xxstart);
   }
  }else{
    test_evaluate_dnn(ncut, xxstart);
  }
}
 void tmva_test_dnn(){
  gROOT ->Reset();
  gROOT->SetStyle("Plain");
  ROOT::EnableImplicitMT(num_threads); 




  double mc_sigma = 0.0400;
  double mc_mass  = 5.27783;
  // You can add an arbitrary number of signal or background trees
  Long64_t CutEveMC=(Long64_t)CutEve/CutEffMC;
  Long64_t CutEveDATA=(Long64_t)CutEve/CutEffDATA;
  std::cout<<Form("mycuteveMC	eventN<%lld",CutEveMC) <<std::endl;
  std::cout<<Form("mycuteveDATA eventN<%lld",CutEveDATA) <<std::endl;
  TCut mycuteveMC   = Form("eventN<%lld",CutEveMC);
  TCut mycuteveDATA = Form("eventN<%lld",CutEveDATA);

  Cut_sig_mass     = Form("(tagged_mass > %f - 2.5*%f) && (tagged_mass < %f + 2.5*%f)",mc_mass,mc_sigma,mc_mass,mc_sigma);

  Cut_ct	   = "( ( (tagB0==1) && (genSignal==1)) || ( (tagB0==0) && (genSignal==2) ) )";

  Cut_bkg_mass     = Form("( ( (tagged_mass > %f-7.*%f ) && (tagged_mass < %f-3.*%f ) ) || ( (tagged_mass > %f+3.*%f ) && (tagged_mass < %f+7*%f ) ) )",mc_mass,mc_sigma,mc_mass,mc_sigma,mc_mass,mc_sigma,mc_mass,mc_sigma );

  Cut_truth_match  = "( (truthMatchMum==1) && (truthMatchMup ==1) && (truthMatchTrkm==1) && (truthMatchTrkp==1) )";

  mycuts0	   = Cut_sig_mass  && Cut_truth_match &&  Cut_ct && "trig == 0 && (mumuMass < 2.702)";

  mycutb0	   = Cut_bkg_mass &&  "(mumuMass < 2.702)";

  TCut RMNone="";
  mycuts=mycuts0&&RMNone;
  mycutb=mycutb0&&RMNone;
  std::cout<<Form("Cuts [MC]   = %s \n", mycuts.GetTitle ())<<std::endl;
  std::cout<<Form("Cuts [Data] = %s \n", mycutb.GetTitle ())<<std::endl;
 } 
//
//=================================================================================================================================================================
//   
 void test_evaluate_dnn(int ncu, float xxstart){
 
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
   foutMC = new TFile(Form("../%sMC_LMNR-ScoreDNN.root",year.c_str()),"READ");
   std::cout<<Form("Opening MC Score DNN File :%s \n",foutMC->GetName())<<std::endl;
   toutMC     = (TTree*)foutMC->Get("sDNN");
   if(toutMC==0) {
     std::cout<<Form("Ntupla sDNN MC not found!")<<std::endl;
    exit(0);
   } 
   toutMC->SetBranchAddress("eventN", &eventN);
   toutMC->SetBranchAddress("ScoreDNN", &ScoreDNN);
   foutData = new TFile(Form("../%sData_LMNR-ScoreDNN.root",year.c_str()),"READ");
   std::cout<<Form("Opening Data Score DNN File :%s \n",foutData->GetName())<<std::endl;
   toutData     = (TTree*)foutData->Get("sDNN");
   if(toutData==0) {
     std::cout<<Form("Ntupla sDNN Data not found!")<<std::endl;
    exit(0);
   } 
   toutData->SetBranchAddress("eventN", &eventN);
   toutData->SetBranchAddress("ScoreDNN", &ScoreDNN);
   pout = new TFile(Form("test-%sDATA_LMNR-Plots.root",year.c_str()),"RECREATE");
 
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
       if(!(Cut_ct && trig == 0 && (mumuMass < 2.702) && pass_preselection==1 )) continue;
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
    Long64_t MaxEveData = (Long64_t)std::min(CutEve,nentriesData);
    std::cout<<Form("Reading num. Data entries= = %lld",MaxEveData)<<std::endl;
    for (Long64_t i=0;i<nentriesData;i++) {
           ScoreDNN = -1.;
    	   TreeData->GetEntry(i);
 	   if ( i%100000==0 ) std::cout<<Form("Event %lld",i)<<std::endl;
	   if(tagged_mass<5.0||tagged_mass>5.6) continue;
           if((mumuMass > 2.702 || pass_preselection!=1)) continue;
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

