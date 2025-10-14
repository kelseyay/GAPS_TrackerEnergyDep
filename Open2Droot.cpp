//New executable for opening the root files after the fit is done

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
sprintf(FilenameRoot,"%s/hcol21.root",reco_path.c_str());

char FilenameNEntries[400];
sprintf(FilenameNEntries,"%s/hnentries.root",reco_path.c_str());
cout << FilenameRoot << endl;
cout << FilenameNEntries << endl;

double mpvmin = 0.4;
double mpvmax = 0.8;

//Number of layers and strips
const int nstrips = 32;
const int nmods = 6;
const int nrows = 6;
const int nlayers = 6;
int dt[4] = {3,4,1,2};

TFile *f = TFile::Open(FilenameRoot);
TFile *fnentries = TFile::Open(FilenameNEntries);

auto hcol = new TH2F("hcol","MPV Full Tracker",nrows*nstrips,0,nrows*nstrips,nlayers*nmods,0,nlayers*nmods);
auto hn = new TH2F("hn","Full Tracker Strip-Level NHits",nrows*nstrips,0,nrows*nstrips,nlayers*nmods,0,nlayers*nmods);

hcol = (TH2F*)f->Get("hcol21");
hn = (TH2F*)fnentries->Get("hnentries");

TCanvas * c1 = new TCanvas("c1", "c1", 200, 10, 900, 900);
c1->SetLeftMargin(0.1);
c1->SetRightMargin(0.23);
c1->SetTopMargin(0.1);
c1->SetBottomMargin(0.1);
hcol->SetBit(TH1::kNoStats);
hcol->GetXaxis()->SetTitle("row(0-5)*32 + strp (0-31)");
hcol->GetYaxis()->SetRangeUser(0,36);
hcol->GetYaxis()->SetTitle("layer(0-5)*6 + mod(0-5)");
hcol->GetZaxis()->SetTitle("Energy Deposition MPV * Cos(#theta) (MeV)");
hcol->GetZaxis()->SetTitleOffset(1.8);
hcol->Draw("COLZ");
hcol->SetMaximum(mpvmax);
hcol->SetMinimum(mpvmin);

string title = "MPVFullTracker";
char histname[400];
sprintf(histname, "%s.png",title.c_str());
c1->SaveAs(histname);

//Histogram for NEntries at a strip level
TCanvas * c2 = new TCanvas("c2", "c2", 200, 10, 900, 900);
c2->SetLeftMargin(0.1);
c2->SetRightMargin(0.23);
c2->SetTopMargin(0.1);
c2->SetBottomMargin(0.1);
hn->SetBit(TH1::kNoStats);
hn->GetXaxis()->SetTitle("row(0-5)*32 + strp (0-31)");
hn->GetYaxis()->SetTitle("layer(0-5)*6 + mod(0-5)");
hn->GetYaxis()->SetRangeUser(0,36);
hn->GetZaxis()->SetTitle("Number of Hits");
hn->GetZaxis()->SetTitleOffset(1.8);
hn->Draw("COLZ");

title = "NEntries";
sprintf(histname, "%s.png",title.c_str());
c2->SaveAs(histname);


}
