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
int numParameter = 4;

//Number of layers and strips
const int nstrips = 32;
const int nmods = 6;
const int nrows = 6;
const int nlayers = 6;
int dt[4] = {3,4,1,2};

//Full tracker histogram range
double mpvmin = 0.55;
double mpvmax = 0.75;

int lyr[nlayers];
int rw[nrows];
int md[nmods];
int strps[nstrips];

for(int l = 0; l < nlayers; l++){lyr[l] = l;}
for(int r = 0; r < nrows; r++){rw[r] = r;}
for(int m = 0; m < nmods; m++){md[m] = m;}
for(int s = 0; s < nstrips; s++){strps[s] = s;}

//Need a histogram and fitting function for every strip
TH1F * h[nlayers][nrows][nmods][nstrips];
TF1 * g1[nlayers][nrows][nmods][nstrips];

auto hcol21 = new TH2F("hcol21","MPV Full Tracker",nrows*nstrips,0,nrows*nstrips,nlayers*nmods,0,nlayers*nmods);
auto hnentries = new TH2F("hnentries","Full Tracker Strip-Level NHits",nrows*nstrips,0,nrows*nstrips,nlayers*nmods,0,nlayers*nmods);

for(int l = 0; l<nlayers;l++){
        for(int r = 0; r<nrows;r++){
                for(int m = 0; m < nmods; m++){
                        for(int j = 0; j < ceil(nstrips/8); j++){
                                cout << "j is " << j << endl;
                                TCanvas * EdepCompare = new TCanvas("EdepCompare", "EdepCompare", 200, 10, 1800, 900);
                                EdepCompare->SetLeftMargin(0.11);
                                EdepCompare->SetRightMargin(0.04);
                                EdepCompare->SetTopMargin(0.04);
                                TLegend* LegEdepCompare = new TLegend(0.5, 0.75, 0.95, 0.95);
                                LegEdepCompare->SetFillColor(0);

                                EdepCompare->Divide(4,2);

                                for(int s = 0 + 8*j; s < 8 + 8*j;s++){
                                        h[l][r][m][s] = new TH1F (TString::Format("hfit0_l%ir%im%is%i",l,r,m,s), ("Edep l" + to_string(lyr[l]) + "r" + to_string(rw[r]) + "m" + to_string(md[m]) + "s" + to_string(strps[s])).c_str(), NBins, xlow,xhigh);
                                        h[l][r][m][s] = (TH1F*)f->Get(TString::Format("h0_l%ir%im%is%i",l,r,m,s));
                                        h[l][r][m][s]->SetLineColor(1);
                                        h[l][r][m][s]->GetXaxis()->SetTitle("Energy Deposition of Hit * Cos(#theta) (MeV)");
                                        h[l][r][m][s]->GetYaxis()->SetTitle("Number of Events");

                                        g1[l][r][m][s] = new TF1("f_landau_gauss",langaufun,fitlow,fithigh,numParameter);
                                        g1[l][r][m][s]->SetParameter(0, 0.1);
                                        g1[l][r][m][s]->SetParameter(1, 0);
                                        g1[l][r][m][s]->SetParameter(2, 1000);
                                        g1[l][r][m][s]->SetParameter(3, 1);

                                        EdepCompare->cd((s % 8)+1);
                                        h[l][r][m][s]->Fit(g1[l][r][m][s],"R");
                                        gPad->SetGridx(1);
                                        gPad->SetGridy(1);
                                        gPad->SetLogy(1);
                                        gStyle->SetTitleW(0.9);
                                        gStyle->SetOptFit();
                                        h[l][r][m][s]->Draw();
                                        //h[l][r][m][s]->SaveAs(TString::Format("hstrip/h_l%ir%im%is%i.root",l,r,m,s));

                                        hcol21->Fill(r*32 + s,l*6 + m,g1[l][r][m][s]->GetParameter(1));
                                        cout << "NEntries should be " << h[l][r][m][s]->GetEntries() << endl;
                                        cout << "X bin should be " << r*32+s << endl;
                                        cout << "Y bin should be " << l*6+m << endl;
                                        hnentries->Fill(r*32+s, l*6+m, h[l][r][m][s]->GetEntries());

                                        //mpv[r*32+s][l*6+m] = g1[l][r][m][s]->GetParameter(1); //Save the calculated MPV, it will be used for the histogram
                                        //myfile << (TString::Format(    "%i \t %i \t %i \t %i \t %f \t %f \t %i \n",l,r,m,s,g1[l][r][m][s]->GetParameter(1),g1[l][r][m][s]->GetParameter(2), static_cast<int>(h[l][r][m][s]->GetEntries())   ));

                                }

                                string title = "FullEdepl" + to_string(lyr[l]) + "r" + to_string(rw[r]) + "m" + to_string(md[m]) + "d" + to_string(dt[j]);
                                char name[400];
                                //sprintf(name, "%s.root",title.c_str());
                                //EdepCompare->SaveAs(name);
                                sprintf(name, "%s.png",title.c_str());
                                EdepCompare->SaveAs(name);
                                delete EdepCompare;
                        }
                }
        }
}

TCanvas * c1 = new TCanvas("c1", "c1", 200, 10, 900, 900);
c1->SetLeftMargin(0.1);
c1->SetRightMargin(0.23);
c1->SetTopMargin(0.1);
c1->SetBottomMargin(0.1);
hcol21->SetBit(TH1::kNoStats);
hcol21->GetXaxis()->SetTitle("row(0-5)*32 + strp (0-31)");
hcol21->GetYaxis()->SetTitle("layer(0-5)*6 + mod(0-5)");
hcol21->GetZaxis()->SetTitle("Energy Deposition MPV * Cos(#theta) (MeV)");
hcol21->GetZaxis()->SetTitleOffset(1.8);
hcol21->Draw("COLZ");
hcol21->SetMaximum(mpvmax);
hcol21->SetMinimum(mpvmin);
hcol21->SaveAs("hcol21.root");


string title = "HistFullTrackerMPV";
char histname[400];
//sprintf(histname, "%s.root",title.c_str());
//c1->SaveAs(histname);
sprintf(histname, "%s.png",title.c_str());
c1->SaveAs(histname);

//Histogram for NEntries at a strip level
TCanvas * c2 = new TCanvas("c2", "c2", 200, 10, 900, 900);
c2->SetLeftMargin(0.1);
c2->SetRightMargin(0.23);
c2->SetTopMargin(0.1);
c2->SetBottomMargin(0.1);
hnentries->SetBit(TH1::kNoStats);
hnentries->GetXaxis()->SetTitle("row(0-5)*32 + strp (0-31)");
hnentries->GetYaxis()->SetTitle("layer(0-5)*6 + mod(0-5)");
hnentries->GetZaxis()->SetTitle("Number of Hits");
hnentries->GetZaxis()->SetTitleOffset(1.8);
hnentries->Draw("COLZ");
hnentries->SaveAs("hnentries.root");

title = "HistFullTrackerNEntriesTest";
sprintf(histname, "%s.root",title.c_str());
c2->SaveAs(histname);
sprintf(histname, "%s.png",title.c_str());
c2->SaveAs(histname);

}
