//Fitting histlist.root output form Edep excutable
using namespace std;

#include "langaufun.C"

#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TStyle.h>
#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TLegend.h>
#include <TColor.h>
#include <TProfile.h>
#include <TMath.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include <string>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <sys/stat.h>
#include "GPlottingTools.hh"

#ifdef USE_BOOST_PROGRAM_OPTIONS
#include "GOptionParser.hh"
#include "GFileIO.hh"
#endif

int main(int argc, char *argv[]){

GOptionParser* parser = GOptionParser::GetInstance();
parser->AddProgramDescription("Minimal Reproducable Example for Extracing Data from Reco Data");
parser->AddCommandLineOption<string>("in_path", "path to instrument data files", "./*", "i");
parser->AddCommandLineOption<string>("out_file", "name of output root file", "out.root", "o");
parser->ParseCommandLine(argc, argv);
parser->Parse();

string reco_path = parser->GetOption<string>("in_path");

char FilenameRoot[400];
sprintf(FilenameRoot,"%s/histlist.root",reco_path.c_str());

cout << FilenameRoot << endl;

TFile *f = TFile::Open(FilenameRoot);

double xlow = 0.3; //Low range for histogram MeV
double xhigh = 4; //High range for histogram MeV
double fitlow = 0.5;
double fithigh = 3.5;
const Int_t NBins = 50;

TH1F * hpstrip = new TH1F ("h0_l0r0m0s0","h0_l0r0m0s0", NBins, xlow,xhigh);
hpstrip = (TH1F*)f->Get("h0_l0r0m0s0");


double min = 0, max = 10; //Change this to your preferences
int numParameter = 4;

TF1*  f_landau_mod = new TF1("f_landau_gauss",langaufun,min,max,numParameter);

f_landau_mod->SetParameter(0, 1);
f_landau_mod->SetParameter(1, 1);
f_landau_mod->SetParameter(2, 1);
f_landau_mod->SetParameter(3, 1);

hpstrip->Fit(f_landau_mod,"R");

//hpstrip->Draw();


TCanvas * c1 = new TCanvas("c1", "c1", 200, 10, 900, 900);
c1->SetLeftMargin(0.12);
c1->SetRightMargin(0.1);
c1->SetTopMargin(0.1);
c1->SetBottomMargin(0.1);

hpstrip->SetLineColor(1);

hpstrip->GetXaxis()->SetRangeUser(xlow, xhigh);
hpstrip->GetXaxis()->SetTitle("Energy Deposition of Hit #times  Cos(#theta) [MeV]");
hpstrip->GetYaxis()->SetTitle("Number of Hits on Track");


gStyle->SetOptStat(0);
gPad->SetLogy(1);
hpstrip->Draw();
hpstrip->SetTitle("");
hpstrip->SaveAs("p4hLGFit.root");
c1->SaveAs("Test.png");

}
