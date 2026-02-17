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
parser->ParseCommandLine(argc, argv);
parser->Parse();

char hsigFilename[400];
string reco_path = parser->GetOption<string>("in_path");
sprintf(hsigFilename,"%s/%s.root",reco_path.c_str(),"hsig");
char hsig_nocorrFilename[400];
sprintf(hsig_nocorrFilename,"%s/%s.root",reco_path.c_str(),"hsig_nocorr");
char hsig_pulsedFilename[400];
sprintf(hsig_pulsedFilename,"%s/%s.root",reco_path.c_str(),"hsig_pulsed");

char hpulserFilename[400];
sprintf(hpulserFilename,"%s/%s.root",reco_path.c_str(),"hpulser");

cout << "hsig filename " << hsigFilename << endl;
cout << "hsig_nocorr filename " << hsig_nocorrFilename << endl;
cout << "hsig_pulsed filename " << hsig_pulsedFilename << endl;
cout << "hpulsedFilename " << hpulserFilename << endl;

//Number of layers and strips
const int nstrips = 32;
const int nmods = 6;
const int nrows = 6;
const int nlayers = 7;
int dt[4] = {3,4,1,2};

const int xthresh = 4; //ADC threshold for X-rays is 4 ADC
const int tthresh = 50; //ADC threshold for tracker is 50 ADC

vector<TLine*> lines;
for (uint i=0; i<nrows; i++)  lines.push_back( new TLine((i+1)*nstrips,            0,(i+1)*nstrips,nlayers*nmods));
for (uint i=0; i<nlayers; i++)lines.push_back( new TLine(             0,(i+1)*nmods, nrows*nstrips, (i+1)*nmods));

TFile *fsig = TFile::Open(hsigFilename);
TFile *fsig_nocorr = TFile::Open(hsig_nocorrFilename);
TFile *fsig_pulsed = TFile::Open(hsig_pulsedFilename);
TFile *f_pulser = TFile::Open(hpulserFilename);

auto hsig = new TH2F("hsig","Histogram",nrows*nstrips,0,nrows*nstrips,nlayers*nmods,0,nlayers*nmods);
auto hsig_pulsed = new TH2F("hsig_pulsed","Histogram",nrows*nstrips,0,nrows*nstrips,nlayers*nmods,0,nlayers*nmods);
auto hsig_nocorr = new TH2F("hsig_nocorr","Histogram",nrows*nstrips,0,nrows*nstrips,nlayers*nmods,0,nlayers*nmods);

auto hpulser = new TH2F("hpulser","Histogram",nrows*nstrips,0,nrows*nstrips,nlayers*nmods,0,nlayers*nmods);

//This is the histogram that categorizes strips based on performance (pulsed channel is categorized)
auto hsigcat = new TH2F("hsigcat","Histogram",nrows*nstrips,0,nrows*nstrips,nlayers*nmods,0,nlayers*nmods);
auto hsig_pulsedcat = new TH2F("hsig_pulsedcat","Histogram",nrows*nstrips,0,nrows*nstrips,nlayers*nmods,0,nlayers*nmods);
auto hsig_nocorrcat = new TH2F("hsig_nocorrcat","Histogram",nrows*nstrips,0,nrows*nstrips,nlayers*nmods,0,nlayers*nmods);
auto hpulsertrack = new TH2F("hpulsertrack","Histogram",nrows*nstrips,0,nrows*nstrips,nlayers*nmods,0,nlayers*nmods);


//This is the histogram that shows the quality of change made by the pulsed channel by changing the color of the pulsed channel
auto hdiff = new TH2F("hdiff","Histogram",nrows*nstrips,0,nrows*nstrips,nlayers*nmods,0,nlayers*nmods);

hsig = (TH2F*)fsig->Get("hsig;1");
hsig_nocorr = (TH2F*)fsig_nocorr->Get("hsig_nocorr;1");
hsig_pulsed = (TH2F*)fsig_pulsed->Get("hsig_pulsed;1");
hpulser = (TH2F*)f_pulser->Get("hpulser;1");

int ncn = 2;
bool pulsed[nlayers][nrows][nmods][ncn];

//cout << "Behold, I will attempt to show the bin content 5,5 = " << h->GetBinContent(5,5) << endl; //This did work :o Now need to convert the 2D histo back to LRMS

for(int x = 0; x < nrows*nstrips;x++){
    for(int y = 0; y < nlayers*nmods+1;y++){

        int strip = floor(x % nstrips);
        int row = floor( (x-strip) / nstrips);
        int mod = floor(y % nmods);
        int layer = floor( (y-mod) / nmods);

        if(hpulser->GetBinContent(x,y) > 2.5){ cout << "Pulsed ch L" <<layer << "R" << row << "M" << mod << "S" << strip << endl; }

    }
}


for(int x = 0; x < nrows*nstrips+1;x++){
    for(int y = 0; y < nlayers*nmods+1;y++){

        //int strip = floor(x-1 % nstrips);
        //int row = floor( (x-1-strip) / nstrips);
        //int mod = floor(y-1 % nmods);
        //int layer = floor( (y-1-mod) / nmods);

        if(hsig_nocorr->GetBinContent(x,y) != 0) hsig_nocorrcat->Fill(x-1,y-1);
        if(hsig_nocorr->GetBinContent(x,y) > xthresh) hsig_nocorrcat->Fill(x-1,y-1);
        if(hsig_nocorr->GetBinContent(x,y) > tthresh) hsig_nocorrcat->Fill(x-1,y-1);

        if(hsig_pulsed->GetBinContent(x,y) == 0 && hsig_nocorrcat->GetBinContent(x,y) != 0){
            //hsig_pulsedcat->SetBinContent(x,y,10);
            hsig_pulsedcat->SetBinContent(x,y,hsig_nocorrcat->GetBinContent(x,y));
        }
        if(hsig_pulsed->GetBinContent(x,y) != 0) hsig_pulsedcat->Fill(x-1,y-1);
        if(hsig_pulsed->GetBinContent(x,y) > xthresh) hsig_pulsedcat->Fill(x-1,y-1);
        if(hsig_pulsed->GetBinContent(x,y) > tthresh) hsig_pulsedcat->Fill(x-1,y-1);

        //cout << "x = " << x << " y = " << y << endl;
        //cout << "l" << layer << "r" << row << "m" << mod << "s" << strip << " sigma " << hsig->GetBinContent(x,y) << endl;

    }
}

TCanvas * c1 = new TCanvas("c1", "c1", 200, 10, 900, 1100);
c1->SetLeftMargin(0.1);
c1->SetRightMargin(0.1);
c1->SetTopMargin(0.05);
c1->SetBottomMargin(0.1);
hsig_pulsedcat->SetBit(TH1::kNoStats);
hsig_pulsedcat->GetXaxis()->SetTitle("row(0-5)*32 + strp (0-31)");
hsig_pulsedcat->GetYaxis()->SetRangeUser(0,42);
hsig_pulsedcat->GetYaxis()->SetTitle("layer(0-6)*6 + mod(0-5)");
hsig_pulsedcat->GetZaxis()->SetTitleOffset(1.8);
gStyle->SetPalette(kRainBow);
hsig_pulsedcat->SetMinimum(1);
hsig_pulsedcat->SetMaximum(10);
hsig_pulsedcat->Draw("COLZ");
for( uint i=0; i<lines.size(); i++)lines.at(i)->Draw("same");

char histname[400];
sprintf(histname, "hsig_pulsedcat.png");
c1->SaveAs(histname);

TCanvas * c2 = new TCanvas("c2", "c2", 200, 10, 900, 1100);
c2->SetLeftMargin(0.1);
c2->SetRightMargin(0.1);
c2->SetTopMargin(0.05);
c2->SetBottomMargin(0.1);
hsig_nocorrcat->SetBit(TH1::kNoStats);
hsig_nocorrcat->GetXaxis()->SetTitle("row(0-5)*32 + strp (0-31)");
hsig_nocorrcat->GetYaxis()->SetRangeUser(0,42);
hsig_nocorrcat->GetYaxis()->SetTitle("layer(0-6)*6 + mod(0-5)");
hsig_nocorrcat->GetZaxis()->SetTitleOffset(1.8);
hsig_nocorrcat->SetMinimum(1);
hsig_nocorrcat->SetMaximum(10);
hsig_nocorrcat->Draw("COLZ");
for( uint i=0; i<lines.size(); i++)lines.at(i)->Draw("same");

sprintf(histname, "hsic_nocorrcat.png");
c2->SaveAs(histname);

}
