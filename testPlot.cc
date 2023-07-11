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
void testPlot(int ncut=10) {
   float xxstart=0.9;
   TFile *FileInput = new TFile("~fgiovenzana/tmva2-dnn/test-2018DATA_LMNR-Plots.root","READ");
   std::vector<TH1D*> HistData_Cut;
   std::vector<TH1D*> HistMC_Cut;
   std::vector<double> cutDNN;		     
   std::vector<TCanvas*> cstudies;	
// summary canvas 
   TCanvas c1("c_allCuts","summary of the cut plots",200,10,1800,780);
   c1.Divide(2,1);
//
   float xxcut=0.;
   for (int  i=0;i<ncut-1;i++) {
    xxcut=xxstart+i*(1.-xxstart)/float(ncut);
    cutDNN.push_back( xxcut);
    std:cout<<Form("Loop Histo cut scoreDNN>%f",xxcut)<<std::endl;
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
    cstudies[i]->Print(Form("testPlot-cut%d.pdf",i));
    c1.cd(1);
    HistMC_Cut[i]->SetLineColor(20+i);
    HistMC_Cut[i]->Draw("same");
    c1.cd(2);
    HistData_Cut[i]->SetLineColor(10+i);
    HistData_Cut[i]->Draw("same");
   }
//  save summary histos
   c1.Print(Form("testPlot-allCuts.pdf"));
   std::cout<<"End of the loop"<<std::endl;

}
