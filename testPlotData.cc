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
#include <filesystem>
#include <RooFit.h>
#include <RooAddPdf.h>
#include "RooFitResult.h"

#include "TLegend.h"

using namespace RooFit;


void testPlot(std::string year, int ncut, float xxstart);



int main (int argc, char** argv) {


if(argc<4){
std::cout<<"use: ./programm year ncut xxstart"<<std::endl;
exit(0);
}else{
std::string year = argv[1];
int ncut = atoi(argv[2]);
float xxstart = atof(argv[3]);

testPlot(year, ncut, xxstart);
	}
}

void testPlot(const std::string year, int ncut, float xxstart) {
   
   TFile *FileInput = new TFile(Form("test-%sDATA_LMNR-Plots-%dCUT.root", year.c_str(), ncut) ,"READ"); 
   std::vector<TH1D*> HistData_Cut;
   std::vector<TH1D*> HistMC_Cut;
   std::vector<double> cutDNN;		     
   std::vector<TCanvas*> cstudies;	
// summary canvas 
   TCanvas c1("c_allCuts","summary of the cut plots",200,10,1800,780);
   c1.Divide(2,1);


double corrFactorDataonMCJpsi = 1.;
double lumiMC = 684.1*98.86/9.4/corrFactorDataonMCJpsi;
double fraz_C = 6./8.;
double fraz_E = 4./16.;
double lumiDatasetC =  7.24;// tot lumi (1/fb) 2024C
double lumiDatasetE =  11.32 ;// tot lumi (1/fb) 2024E
double lumiData = fraz_C * lumiDatasetC + fraz_E * lumiDatasetE;
double scale24_calc = lumiData / lumiMC ; 


   std::string OutDir = Form("DATA_TestPlot_results_%d", ncut);
   
   if(!gSystem->AccessPathName((const char *)OutDir.c_str())){
    std::cout<<Form("Exist output Dir %s, ...save and make it again",OutDir.c_str() )<<std::endl;
    gSystem->Exec(Form("mv %s %s.`/usr/bin/date -Iminutes -r %s`",OutDir.c_str(),OutDir.c_str(),OutDir.c_str()));
    gSystem->Exec(Form("mkdir %s",OutDir.c_str()));
   }else{
    std::cout<<Form("output Dir %s doesn't exist, let's make it",OutDir.c_str() )<<std::endl;
    gSystem->Exec(Form("mkdir %s",OutDir.c_str()));
   }


// create saving directory
// std::string folder_name = Form("DATA_TestPlot_results_%d", ncut);;
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
// Data --> Exponential pdf   

std::vector<RooDataHist*> MC_Hist;
std::vector<RooDataHist*> Data_Hist;
std::vector<TCanvas*> cMCfitted;
std::vector<TCanvas*> cDatafitted;

// mean mass and sigma values taken from previous analysis work
const double mc_sigma = 0.033;
const double mc_mass  = 5.280;

const double B0Mass_ = 5.27958; // nominal mass for bkg section

RooRealVar mass ("mass", "mass", 5.0, 5.6); // declaring var
mass.setRange("full", 5.0,5.6);

RooRealVar mean ("mean", "mean of CryBall", B0Mass_,  5.2, 5.35);
RooRealVar sigma("sigma", "CryBall_sigma", 0.03, 0.001, 0.3);
RooRealVar aL("L_tail", "L_gauss_tail", 1.2, 1., 2.5);
RooRealVar nL("L_norm", "L_norm", 0.9, 10.);
RooRealVar aR("R_tail", "R_gauss_tail", 1.3, 1., 2.5);
RooRealVar nR("R_norm", "R_norm", 0.8, 20.);

// Data: declaring pars
RooRealVar alpha("alpha", "exp_alpha", -9, 0.);


//RooRealVar mean_d("mean_d", "mean_d", B0Mass_, 5.0, 5.6);

// setting exponential fit range ( x < mc_mass - 3mc_sigma || x> mc_mass + 3mc_sigma )
mass.setRange("R1", 5.00, mc_mass - 3.*mc_sigma);
mass.setRange("R2", mc_mass + 3.*mc_sigma, 5.6);
mass.setRange("signRangeMC", B0Mass_ - 3.*sigma.getVal(), B0Mass_ + 3.*sigma.getVal() );

//significance vars
double S;
double B;
std::vector<double> signifDNN;
std::vector<double> errSignifDNN;
// scaling vars
//double scale18 = (14.741 + 7.149 + 6.899 + 32.347)/12321/1.30; 
//double rescale = 1./sqrt(11); //1/sqrt(11) is due to the n° of subsamples

std::ofstream ratio_rslts; //file with s/sqrt(s+b) from various cuts
std::ofstream yield;
std::ofstream scale;
ratio_rslts.open(Form("DATA_TestPlot_results_%d/ratio_results_%d.txt", ncut, ncut));
yield.open(Form("DATA_TestPlot_results_%d/yeald_%d.txt", ncut, ncut));
scale.open(Form("DATA_TestPlot_results_%d/scale.txt", ncut));
// file will be closed at the end of loop

// loop on cuts 

RooRealVar nsig_MC("nsig_MC" , "signal frac"   ,  3500000.,     0.,   5000000.);
RooRealVar nsig("nsig", "signal frac"   ,    2000.,     0,   100000.);
//RooRealVar nbkg("nbkg", "bckg fraction"  ,  8000.,       0.,   500000.);

double scale24 = 0; 
double scale_lum = 0;

double index = 0.;
double max_signi = 0;
double max_cut = 0;
double max_S = 0;
double max_B = 0;
 

for (int i = 0; i < ncut-1; i++) {
	xxcut=xxstart+i*(1.-xxstart)/float(ncut);
	
	sigma.setConstant(kFALSE);
	//mean.setConstant(kFALSE);
	aL.setConstant(kFALSE);
	aR.setConstant(kFALSE);
	nL.setConstant(kFALSE);
	nR.setConstant(kFALSE);
	mean.setConstant(kFALSE);
	//mean.setVal(mc_mass);

	// MC area 
	MC_Hist.push_back(new RooDataHist(Form("MC_hist_%d", i), Form("ScoreDNN>%f",xxcut),mass, HistMC_Cut[i])); // loading the dataset
	cMCfitted.push_back(new TCanvas(Form("cMCfit_cut_%d",i),Form("ScoreDNN>%f",xxcut),200,10,1000,780));
	//MC_Hist[i]->Draw();
	RooCrystalBall cryball_MC("cryball", "crystalball pdf", mass, mean, sigma, aL, nL, aR, nR); // creating the Crystal Ball pdf
	RooExtendPdf signal_MC("signal_MC","signal"    ,  cryball_MC, nsig_MC, "full");
	
	std::cout<<Form("Here fit results from MC for %d cut\n", i)<<std::endl;
	RooFitResult* mc_r = signal_MC.fitTo(*MC_Hist[i], Extended(true),Strategy(1),MaxCalls(50000), Minimizer("Minuit"), Save());	
	mc_r -> Print();

	
	//update crystal ball parameter for the Data Fit.
	sigma.setConstant(kTRUE);
	//mean.setConstant(kTRUE);
	aL.setConstant(kTRUE);
 	aR.setConstant(kTRUE);
 	nL.setConstant(kTRUE);
 	nR.setConstant(kTRUE);


	// Data area 
	Data_Hist.push_back(new RooDataHist(Form("Data_hist_%d", i), Form("ScoreDNN>%f",xxcut),mass, HistData_Cut[i]));  	
	cDatafitted.push_back(new TCanvas(Form("cDatafit_cut_%d",i),Form("ScoreDNN>%f",xxcut),200,10,1000,780));
        Data_Hist[i]->Draw();


	RooRealVar nbkg("nbkg", "bckg fraction"  ,  8000.,       0.,   500000.);

	RooExponential expon("exponential", "exponential pdf", mass, alpha); // creating the exponential pdf
	RooCrystalBall cryball("cryball", "crystalball pdf", mass, mean, sigma, aL, nL, aR, nR); // creating the Crystal Ball pdf for the data. 

	RooExtendPdf signal("signal","signal"    ,  cryball, nsig, "full");
	RooExtendPdf background("background","background", expon, nbkg, "full");
	
	RooAddPdf model("model", "model", RooArgList(signal, background), RooArgList(nsig, nbkg));

	std::cout<<Form("Here fit results from data for %d cut\n", i)<<std::endl;
	RooFitResult* d_r = model.fitTo(*Data_Hist[i],Extended(true),Strategy(2),MaxCalls(50000),Save());

	d_r -> Print();


	//int status = d_r ->status();

	double err_nsig = nsig.getError();
	double err_nbkg = nbkg.getError();

	std::cout<<Form("End of fit results from data for %d cut\n", i)<<std::endl;

	mass.setRange("signRangeMC", B0Mass_ - 3.*sigma.getVal(), B0Mass_ + 3.*sigma.getVal() ); //redef mass in dep from signal itself
	//RooAbsReal* BckgFullMassS = expon.createIntegral(mass, mass, "signRangeMC"); // integral of expon under the signal
	RooAbsReal* BckgRangeMassS = expon.createIntegral(mass, mass, "signRangeMC" );
//	RooAbsReal* SignalRangeMassS = cryball.createIntegral(mass, mass, "signRangeMC" );

	double err_BckgRangeMassS = BckgRangeMassS->getPropagatedError(*d_r);
//	double err_SignalRangeMassS = SignalRangeMassS->getPropagatedError(*d_r);

//	RooAbsReal* BckgFullMassS = background.createIntegral(mass, mass, "full" );
	RooAbsReal* SignalFullMassS = signal.createIntegral(mass, mass, "full" );

	double err_SignalFullMassS = SignalFullMassS->getPropagatedError(*d_r);

	mass.setRange("R1", 5.00, mc_mass - 3.*mc_sigma);
	mass.setRange("R2", mc_mass + 3.*mc_sigma, 5.6);
	
	
	//drawing the data plots
	RooPlot* frame_d = mass.frame(); 
	Data_Hist[i]->plotOn(frame_d);
        frame_d->GetXaxis()->SetTitle("#mu^{+}#mu^{-}K#pi mass (GeV/c^{2})");
        frame_d->GetYaxis()->SetTitle(Form("Events/%4.4f (GeV/c^{2})",HistData_Cut[i]->GetBinWidth(1)));
	//expon.plotOn(frame_d,Range("R1,R2"),NormRange("R1,R2"));
	//expon.paramOn(frame_d, Layout(0.1, 0.99));
	//cryball.paramOn(frame_d, Layout(0.55, 0.99));
	//model.plotOn(frame_d);
	//data_Hist.statOn(frame_mc, Layout(0.01, 0.40));

	// Componenti riempite
	model.plotOn(frame_d, Components(signal),
             FillColor(kGreen), FillStyle(3001), DrawOption("F"), Name("sig_fill"));

	model.plotOn(frame_d, Components(background),
             FillColor(kRed), FillStyle(3002), DrawOption("F"), Name("back_fill"));

// Curva totale sopra
	background.plotOn(frame_d, LineColor(kRed), LineStyle(9), Normalization(nbkg.getVal(), RooAbsReal::NumEvent), Name("back_line"));
	model.plotOn(frame_d, LineColor(kBlue), Name("curve_model"));


	
	
	double chiSq = frame_d->chiSquare(3); //se inserisco il numero di parametri del fit (ndf) ottengo quello ridotto
	
	
	B = BckgRangeMassS->getVal()*nbkg.getVal();
	S = SignalFullMassS->getVal()*nsig.getVal();

	//std::cout << "Signal full mass s " << SignalFullMassS->getVal() << std::endl;

	double err_B = sqrt(pow(nbkg.getVal()*err_BckgRangeMassS,2.)+pow(BckgRangeMassS->getVal()*err_nbkg,2.));
	double err_S = sqrt(pow(nsig.getVal()*err_SignalFullMassS,2.)+pow(SignalFullMassS->getVal()*err_nsig,2.));

	

	//error evaluate (nbkg && nsig with no error)
	double dS=((sqrt(S + B)- S / (2. * sqrt(S + B )))/(S + B));
	double dB=(S/(2*(pow((S + B),(3./2.)))));
	double err = sqrt((pow((dS*err_S),2.))+(pow((dB*err_B),2.)));

	signifDNN.push_back( S / (sqrt(S + B)));
	errSignifDNN.push_back(err);
	
	ratio_rslts << xxcut << " " << signifDNN[i]<<" 0 "<< errSignifDNN[i] << std::endl; 
	yield << xxcut<< " " << S << " " << B << std::endl; 
	
	scale << xxcut << " " << S/nsig_MC.getVal() << std::endl; 
	scale24 = scale24 +  S/nsig_MC.getVal();
	scale_lum = scale_lum +  nsig.getVal() / (scale24_calc * nsig_MC.getVal()) ;
	 

	 if (signifDNN[i]>max_signi){ 

	  index = i;	
 	  max_signi = signifDNN[i];
 	  max_cut = xxcut;
 	  max_S = S;
 	  max_B = B;
        }

	frame_d->Draw();
	frame_d->SetTitle(Form("DNN Score>%4.2f S/(sqrt(S + B)) = %4.2f  #chi^{2}/ndf = %4.2f", xxcut, signifDNN[i], chiSq));

	TLegend *leg = new TLegend(0.65, 0.65, 0.88, 0.88);

	leg->AddEntry(frame_d->findObject("sig_fill"),
              "Signal (area)", "f");

	leg->AddEntry(frame_d->findObject("back_fill"),
              "Background (area)", "f");

	leg->AddEntry(frame_d->findObject("back_line"),
              "Background fit", "l");

	leg->AddEntry(frame_d->findObject("curve_model"),
              "Total fit", "l");

	leg->Draw();

	cDatafitted[i]->Update();
	cDatafitted[i]->Print(Form("%s/Data_fitted_result_%d.pdf", OutDir.c_str(), i), "pdf");
}

double scale_mean = scale24 / ncut;
double mean_scale_lum = scale_lum / ncut ;

std::cout << "mean scale factor: "<<scale_mean << " exstimated scale factor (luminosity) " << mean_scale_lum <<  std::endl;
std::cout<<"MAX info: " << index << " " <<max_cut << " " << max_signi << " " << max_S << " " << max_B << std::endl;

// closing significance file
ratio_rslts.close();
yield.close();
scale.close();



}


