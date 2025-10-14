//New executable for opening the root files after the fit is done
//Example how to use, navigate into directory of .root file:
// /home/kelsey/GAPS_TrackerEnergyDep/build/ArbOpen2Droot -r hsig_pulsed -i . -m 1 -x 10


//Fitting histlist.root output form Edep excutable
using namespace std;

#include "langaufun.C"

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
parser->AddCommandLineOption<string>("in_path", "path to instrument data files", "./*", "i");
parser->AddCommandLineOption<string>("out_file", "name of output root file", "out.root", "o");
parser->AddCommandLineOption<string>("root_title", "name of root histogram", "h.root", "r");
parser->AddCommandLineOption<string>("zmin", "Minimum Z value 2D Histogram", "", "m");
parser->AddCommandLineOption<string>("zmax", "Maximum Z value 2D Histogram", "", "x");
parser->ParseCommandLine(argc, argv);
parser->Parse();

string reco_path = parser->GetOption<string>("in_path");
string root_name = parser->GetOption<string>("root_title");
int Zmin = stoi(parser->GetOption<string>("zmin"));
int Zmax = stoi(parser->GetOption<string>("zmax"));

char FilenameRoot[400];
sprintf(FilenameRoot,"%s/%s.root",reco_path.c_str(),root_name.c_str());

cout << FilenameRoot << endl;

//Number of layers and strips
const int nstrips = 32;
const int nmods = 6;
const int nrows = 6;
const int nlayers = 7;
int dt[4] = {3,4,1,2};

vector<TLine*> lines;
for (uint i=0; i<nrows; i++)  lines.push_back( new TLine((i+1)*nstrips,            0,(i+1)*nstrips,nlayers*nmods));
for (uint i=0; i<nlayers; i++)lines.push_back( new TLine(             0,(i+1)*nmods, nrows*nstrips, (i+1)*nmods));

TFile *f = TFile::Open(FilenameRoot);

auto h = new TH2F(root_name.c_str(),"Histogram",nrows*nstrips,0,nrows*nstrips,nlayers*nmods,0,nlayers*nmods);

h = (TH2F*)f->Get((root_name + ";1").c_str());

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
for( uint i=0; i<lines.size(); i++)lines.at(i)->Draw("same");
h->SetMaximum(Zmax);
h->SetMinimum(Zmin);

char histname[400];
sprintf(histname, "%s.png",root_name.c_str());
c1->SaveAs(histname);

}
