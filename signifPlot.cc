// This simple macro takes signif .txt from testPlot.cc, in order to combine different runs in a single graph
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sys/stat.h>
#include <dirent.h>
#include "Riostream.h"
#include <map>
#include <string>
//#include "RooGlobalFunc.h"
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
// #include <RooRealVar.h>
// #include <RooConstVar.h>
// #include <RooPlot.h>
// #include <RooDataSet.h>
// #include <RooDataHist.h>
// #include <RooHistPdf.h>
// #include <RooCrystalBall.h>
// #include <RooExponential.h>
// #include <RooExtendPdf.h>

void signifPlot();

int main () {
signifPlot();
}

void signifPlot() {

// reading file 1: 20 pts from 0.9
std::ifstream in_file_1 ("TestPlot_ddf90_results/ratio_results.txt");
std::vector<double> v_col1;
std::vector<double> v_col2;
	
double input1_col1 = 0.;
double input1_col2 = 0.;
	
while (true) {
	if (in_file_1.eof() == true) break;
	//std::cout << "file1" << std::endl;
	in_file_1 >> input1_col1;
	in_file_1 >> input1_col2;
	v_col1.push_back (input1_col1);
	v_col2.push_back (input1_col2);
}

in_file_1.close();


// reading file 2: 20 pts from 0.96
std::ifstream in_file_2 ("TestPlot_ddf96_results/ratio_results.txt");
	
double input2_col1 = 0.;
double input2_col2 = 0.;
	
while (true) {
	if (in_file_2.eof() == true) break;
	//std::cout << "file2" << std::endl;
	in_file_2 >> input2_col1;
	in_file_2 >> input2_col2;
	v_col1.push_back (input2_col1);
	v_col2.push_back (input2_col2);
}

in_file_2.close();

// plotting the significance

TGraphErrors cuts_signif (v_col1.size(), &v_col1[0], &v_col2[0]);

cuts_signif.SetMarkerStyle(20);
cuts_signif.SetMarkerColor(kBlue);
cuts_signif.GetHistogram() ->GetXaxis()->SetTitle("DNN Score");
cuts_signif.GetHistogram() ->GetYaxis()->SetTitle("#frac{s}{#sqrt{s+b}}");
cuts_signif.SetTitle("Figure of merit DNN");

TCanvas cSign;
cuts_signif.Draw("AP");
cSign.SetGrid(1, 10);
cSign.Update();
cSign.Print("WorkPointStudy_CombCuts.pdf");

}
