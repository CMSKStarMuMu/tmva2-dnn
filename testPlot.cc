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
#include "TGraphErrors.h"
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
#include <TText.h>
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
#include <RooRealVar.h>
#include <RooConstVar.h>
#include <RooPlot.h>
#include <RooDataSet.h>
#include <RooDataHist.h>
#include <RooHistPdf.h>
#include <RooCrystalBall.h>
#include <RooExponential.h>
#include <RooExtendPdf.h>
void testPlot(int ncut=20);
using namespace RooFit;
int main (int argc, char** argv) {
 testPlot();
}

void testPlot(int ncut) {
   float xxstart=0.96;
   TFile *FileInput = new TFile("~fgiovenzana/Tesi/new_bins_tmva2-dnn/test-2018DATA_LMNR-Plots.root","READ");
   std::vector<TH1D*> HistData_Cut;
   std::vector<TH1D*> HistMC_Cut;
   std::vector<double> cutDNN;		     
   std::vector<TCanvas*> cstudies;	
// summary canvas 
   TCanvas c1("c_allCuts","summary of the cut plots",200,10,1800,780);
   c1.Divide(2,1);
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
    cstudies[i]->Print(Form("TestPlot_results/testPlot-cut%d.pdf",i));
// summary histos
    c1.cd(1);
    HistMC_Cut[i]->SetLineColor(2+i);
    HistMC_Cut[i]->Draw("same");
    c1.cd(2);
    HistData_Cut[i]->SetLineColor(1+i);
    HistData_Cut[i]->Draw("same");
   }
//  save summary histos
   c1.Print(Form("TestPlot_results/testPlot-allCuts.pdf"));
   std::cout<<"End of the loop"<<std::endl;
   
/* Here fitting the two pdfs, using Roofit */

/*  MC   --> CrystalBall pdf with symmetrical gaussian core & asymmetrical tails 
    Data --> Exponential pdf */  

std::vector<RooDataHist*> MC_Hist;
std::vector<RooDataHist*> Data_Hist;
std::vector<TCanvas*> cMCfitted;
std::vector<TCanvas*> cDatafitted;

// mean mass and sigma values taken from previous analysis work
const double mc_sigma = 0.0400;
const double mc_mass  = 5.27783;

const double B0Mass_ = 5.27958; // nominal mass for bkg section

RooRealVar mass ("mass", "mass", 5., 5.6); // declaring var
mass.setRange("full", 5.0,5.6);
// MC: declaring pars
RooRealVar mean ("mean", "mean of CryBall", 5.0, 5.6);
RooRealVar sigma("sigma", "CryBall_sigma", 0.01, 0.2);
RooRealVar aL("L_tail", "L_gauss_tail", 1., 2.5);
RooRealVar nL("L_norm", "L_norm", 0.9, 10.);
RooRealVar aR("R_tail", "R_gauss_tail", 1., 2.5);
RooRealVar nR("R_norm", "R_norm", 0.9, 10.);
// Data: declaring pars
RooRealVar alpha("alpha", "exp_alpha", -10., 0.);

// setting exponential fit range ( x < mc_mass - 3mc_sigma || x> mc_mass + 3mc_sigma )
mass.setRange("R1", 5.00, mc_mass - 3.*mc_sigma);
mass.setRange("R2", mc_mass + 3.*mc_sigma, 5.6);



//significance vars
double n_signal;
double n_background;
std::vector<double> signifDNN;

// scaling vars
double scale18 = (14.741 + 7.149 + 6.899 + 32.347)/12321/1.30; 
double rescale = 1./sqrt(11); //1/sqrt(11) is due to the n° of subsamples

std::ofstream ratio_rslts; //file with s/sqrt(s+b) from various cuts
ratio_rslts.open("TestPlot_results/ratio_results.txt");
// file will be closed at the end of loop

// loop on cuts 
for (int i = 0; i < ncut-1; i++) {
	xxcut=xxstart+i*(1.-xxstart)/float(ncut);
	
        RooRealVar nsig("Yield"         , "signal frac"   ,   40000,     0,   1000000);
	
	/* MC area */
	MC_Hist.push_back(new RooDataHist(Form("MC_hist_%d", i), Form("ScoreDNN>%f",xxcut),mass, HistMC_Cut[i])); // loading the dataset
	cMCfitted.push_back(new TCanvas(Form("cMCfit_cut_%d",i),Form("ScoreDNN>%f",xxcut),200,10,1000,780));
	MC_Hist[i]->Draw();
	RooCrystalBall cryball("cryball", "crystalball pdf", mass, mean, sigma, aL, nL, aR, nR); // creating the Crystal Ball pdf
        RooExtendPdf signal("signal","signal"    ,  cryball, nsig, "full");
	
	std::cout<<Form("Here fit results from MC for %d cut\n", i)<<std::endl;
	RooFitResult* mc_r = signal.fitTo(*MC_Hist[i],Extended(true),SumW2Error(true), Save()); //fitting and saving MC results
	//mc_r->Print(); //printing MC fit results
	std::cout<<Form("End of fit results from MC for %d cut\n", i)<<std::endl;
	
	//drawing the MC plots
	RooPlot* frame_mc = mass.frame(); 
	MC_Hist[i]->plotOn(frame_mc);
	frame_mc->GetYaxis()->SetTitleOffset(1.55);
        frame_mc->GetXaxis()->SetTitle("#mu^{+}#mu^{-}K#pi mass (GeV/c^{2})");
        frame_mc->GetYaxis()->SetTitle(Form("Events/%4.4f (GeV/c^{2})",HistMC_Cut[i]->GetBinWidth(1)));
	cryball.plotOn(frame_mc);
	cryball.paramOn(frame_mc, Layout(0.55, 0.99));
	//MC_Hist.statOn(frame_mc, Layout(0.01, 0.40));
	// printing chisquare
//	TLatex Tl;
//	Tl.SetTextSize(0.04);
	double chiSq = frame_mc->chiSquare();
	
	
	
	
	frame_mc->Draw();
	frame_mc->SetTitle(Form("MassMC_ScoreDNN>%f", xxcut));
	cMCfitted[i]->Update();
//	Tl.DrawLatexNDC(0.2, 0.86, Form("#chi^{2}/N_{bins} = %.2f", chiSq));
	cMCfitted[i]->Print(Form("TestPlot_results/MC_fitted_result_%d.pdf", i), "pdf");
	
	
	/* Data area */
	Data_Hist.push_back(new RooDataHist(Form("Data_hist_%d", i), Form("ScoreDNN>%f",xxcut),mass, HistData_Cut[i]));  	
	cDatafitted.push_back(new TCanvas(Form("cDatafit_cut_%d",i),Form("ScoreDNN>%f",xxcut),200,10,1000,780));
    	Data_Hist[i]->Draw();
	
	mass.setRange("signRangeMC", B0Mass_ - 3.*sigma.getVal(), B0Mass_ + 3.*sigma.getVal() );
        RooRealVar nbkg("nbkg", "bkg n"  ,  1000.,       0.,   550000.);
	RooExponential expon("exponential", "exponential pdf", mass, alpha); // creating the exponential pdf
        RooExtendPdf bkg("bkg","ebkg"    ,  expon, nbkg, "full");
	std::cout<<Form("Here fit results from data for %d cut\n", i)<<std::endl;
	RooFitResult* d_r = bkg.fitTo(*Data_Hist[i],Extended(true), Range("R1,R2"), NormRange("full"), Save()); //fitting and saving data results 
//	RooFitResult* d_r = expon.fitTo(*Data_Hist[i], Range("R1,R2"), NormRange("R1,R2"), Save()); //fitting and saving data results 
	//d_r->Print(); //printing data fit results
	std::cout<<Form("End of fit results from data for %d cut\n", i)<<std::endl;


	mass.setRange("signRangeMC", B0Mass_ - 3.*sigma.getVal(), B0Mass_ + 3.*sigma.getVal() ); //redef mass in dep from signal itself
	RooAbsReal* BckgFullMassS = expon.createIntegral(mass, mass, "signRangeMC"); // integral of expon under the signal
	
	BckgFullMassS->getVal()*nbkg.getVal();
//	Data_Hist[i]->GetXAxis()->SetTitle("#mu^{+}#mu^{-}K#pi mass (GeV/c^{2})");
//	Data_Hist[i]->GetYAxis()->SetTitle(Form("Events/%4.4f (GeV/c^{2})",Data_Hist[i]->GetBinWidth));
	
	//drawing the data plots
	RooPlot* frame_d = mass.frame(); 
	Data_Hist[i]->plotOn(frame_d);
        frame_d->GetXaxis()->SetTitle("#mu^{+}#mu^{-}K#pi mass (GeV/c^{2})");
        frame_d->GetYaxis()->SetTitle(Form("Events/%4.4f (GeV/c^{2})",HistData_Cut[i]->GetBinWidth(1)));
	expon.plotOn(frame_d,Range("R1,R2"),NormRange("R1,R2"));
	expon.paramOn(frame_d, Layout(0.55, 0.99));
	//data_Hist.statOn(frame_mc, Layout(0.01, 0.40));
	
	chiSq = frame_d->chiSquare();

	frame_d->Draw();
	frame_d->SetTitle(Form("MassData_ScoreDNN>%f", xxcut));
//	Tl.DrawLatexNDC(0.2, 0.86, Form("#chi^{2}/N_{bins} = %.2f", chiSq));
	cMCfitted[i]->Print(Form("TestPlot_results/MC_fitted_result_%d.pdf", i), "pdf");
	cDatafitted[i]->Update();
	cDatafitted[i]->Print(Form("TestPlot_results/Data_fitted_result_%d.pdf", i), "pdf");
	
	
	// evaluating the significance and printing on file
	n_signal = nsig.getVal()*scale18; //scaling factor added to compensate differences in Luminosity
	//std::cout << Form("nsig cut %d = %lf\n", i, n_signal)<<std::endl;
	n_background = BckgFullMassS->getVal()*nbkg.getVal();
	//n_background = nbkg.getVal();
	//std::cout << Form("nbkg cut %d = %lf\n", i, n_background)<<std::endl;
	signifDNN.push_back(rescale*n_signal / (sqrt(n_signal + n_background))); 
	ratio_rslts << xxcut << " "<< signifDNN[i] <<std::endl;
}

// closing significance file
ratio_rslts.close();



}


