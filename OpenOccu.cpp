//New executable for opening the root files after the fit is done
//Example how to use, navigate into directory of .root file:
// /home/kelsey/GAPS_TrackerEnergyDep/build/ArbOpen2Droot -r hsig_pulsed -i . -m 1 -x 10


//Fitting histlist.root output form Edep excutable
using namespace std;

#include "langaufun.C"
#include "KYtools.C"

#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TStyle.h>
#include <TLine.h>
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
parser->AddCommandLineOption<string>("in_path", "path to instrument data files", "", "i");
parser->AddCommandLineOption<string>("out_path", "name of output root file", "", "o");
parser->AddCommandLineOption<string>("root_title", "name of root histogram", "HTofUmbOccu", "r");
parser->AddCommandLineOption<int>("zmin", "Minimum Z value 2D Histogram", 0, "m");
parser->AddCommandLineOption<int>("zmax", "Maximum Z value 2D Histogram", 1000, "x");
parser->ParseCommandLine(argc, argv);
parser->Parse();

string reco_path = parser->GetOption<string>("in_path");
string out_path = parser->GetOption<string>("out_path");
string root_name = parser->GetOption<string>("root_title");
int Zmin = (parser->GetOption<int>("zmin"));
int Zmax = (parser->GetOption<int>("zmax"));

cout << Zmin << endl;
cout << Zmax << endl;

//char FilenameRoot[400];
//sprintf(FilenameRoot,"",reco_path.c_str(),root_name.c_str());

//cout << FilenameRoot << endl;

if(reco_path != "" && reco_path[reco_path.length()-1] != '/' ){ cout <<  "reco path no slash!" << endl; reco_path = reco_path + '/'; }
if(out_path != "" && out_path[out_path.length()-1] != '/' ){ cout <<  "out path no slash!" << endl; out_path = out_path + '/'; }


TFile *f = new TFile((reco_path + root_name + ".root").c_str());

TH2D *h = (TH2D*)f->Get(root_name.c_str());

//TH2D *h = (TH2D*)f->Get(root_name.c_str());
h->SetMinimum(Zmin);
h->SetMaximum(Zmax);
histplot2d("c1", h, root_name.c_str(),"X Location Hit (mm)","Y Location Hit (mm)","NEntries", out_path + "OccuTest");
//TFile *f = TFile::Open(FilenameRoot);
//TH2D* h = new TH2D(root_name.c_str(), 168, -2000, 2000, 168, -2000, 2000, "rec. hit position x [mm]", "rec. hit position y [mm]", "events", Zmin, Zmax));

//h = (TH2F*)f->Get((root_name + ";1").c_str());

/*
TCanvas * c1 = new TCanvas("c1", "c1", 200, 10, 900, 1100);
c1->SetLeftMargin(0.1);
c1->SetRightMargin(0.1);
c1->SetTopMargin(0.05);
c1->SetBottomMargin(0.1);
h->SetBit(TH1::kNoStats);
h->GetXaxis()->SetTitle("row(0-5)*32 + strp (0-31)");
h->GetYaxis()->SetRangeUser(0,42);
h->GetYaxis()->SetTitle("layer(0-6)*6 + mod(0-5)");
//h->GetZaxis()->SetTitle("Energy Deposition MPV * Cos(#theta) (MeV)");
h->GetZaxis()->SetTitleOffset(1.8);
h->Draw("COLZ");
h->SetMaximum(Zmax);
h->SetMinimum(Zmin);

char histname[400];
sprintf(histname, "%s.png",root_name.c_str());
c1->SaveAs(histname);
*/

}
