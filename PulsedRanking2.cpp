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
int ncn = 2;

vector<TLine*> lines;
//for (uint i=0; i<nrows; i++)  lines.push_back( new TLine((i+1)*nstrips,            0,(i+1)*nstrips,nlayers*nmods));
for (uint i=0; i<nrows; i++)  lines.push_back( new TLine((i+1)*ncn,            0,(i+1)*ncn,nlayers*nmods));
for (uint i=0; i<nlayers; i++)lines.push_back( new TLine(             0,(i+1)*nmods, nrows*nstrips, (i+1)*nmods));

TFile *fsig = TFile::Open(hsigFilename);
TFile *fsig_nocorr = TFile::Open(hsig_nocorrFilename);
TFile *fsig_pulsed = TFile::Open(hsig_pulsedFilename);
TFile *f_pulser = TFile::Open(hpulserFilename);

auto hsig = new TH2F("hsig","Histogram",nrows*nstrips,0,nrows*nstrips,nlayers*nmods,0,nlayers*nmods);
auto hsig_pulsed = new TH2F("hsig_pulsed","Histogram",nrows*nstrips,0,nrows*nstrips,nlayers*nmods,0,nlayers*nmods);
auto hsig_nocorr = new TH2F("hsig_nocorr","Histogram",nrows*nstrips,0,nrows*nstrips,nlayers*nmods,0,nlayers*nmods);

auto hdiff = new TH2F("hdiff","",ncn*nrows,0,ncn*nrows,nlayers*nmods,0,nlayers*nmods);

//auto hpulser = new TH2F("hpulser","Histogram",nrows*nstrips,0,nrows*nstrips,nlayers*nmods,0,nlayers*nmods);
//auto hpulsertrack = new TH2F("hpulsertrack","Histogram",nrows*nstrips,0,nrows*nstrips,nlayers*nmods,0,nlayers*nmods);

//This is the histogram that shows the quality of change made by the pulsed channel by changing the color of the pulsed channel

hsig = (TH2F*)fsig->Get("hsig;1");
hsig_nocorr = (TH2F*)fsig_nocorr->Get("hsig_nocorr;1");
hsig_pulsed = (TH2F*)fsig_pulsed->Get("hsig_pulsed;1");
//hpulser = (TH2F*)f_pulser->Get("hpulser;1");



for(int x = 0; x < nrows*nstrips;x++){
    for(int y = 0; y < nlayers*nmods;y++){

        int strip = floor(x % nstrips);
        int row = floor( (x-strip) / nstrips);
        int mod = floor(y % nmods);
        int layer = floor( (y-mod) / nmods);
        int ncn_val = floor(strip/(nstrips/ncn));

        if(hsig_nocorr->GetBinContent(x+1,y+1)!= 0){hdiff->Fill(row*ncn+ncn_val,layer*nmods+mod);}

    }
}

for(int x = 0; x < nrows*nstrips;x++){
    for(int y = 0; y < nlayers*nmods;y++){

        int strip = floor(x % nstrips);
        int row = floor( (x-strip) / nstrips);
        int mod = floor(y % nmods);
        int layer = floor( (y-mod) / nmods);
        int ncn_val = floor(strip/(nstrips/ncn));
        //cout << "L" <<layer << "R" << row << "M" << mod << "C" << strip << "ncn" << ncn_val << endl;

        if(hsig_nocorr->GetBinContent(x+1,y+1) > xthresh && hsig_pulsed->GetBinContent(x+1,y+1) < xthresh && hsig_pulsed->GetBinContent(x+1,y+1) > 0){
            hdiff->Fill(row*ncn+ncn_val,layer*nmods+mod);
            //It's a point!!
        }

    }
}


TCanvas * c1 = new TCanvas("c1", "c1", 200, 10, 900, 1100);
c1->SetLeftMargin(0.1);
c1->SetRightMargin(0.1);
c1->SetTopMargin(0.05);
c1->SetBottomMargin(0.1);
hdiff->SetBit(TH1::kNoStats);
hdiff->GetXaxis()->SetTitle("row(0-5)*32 + strp (0-31)");
hdiff->GetYaxis()->SetRangeUser(0,42);
hdiff->GetYaxis()->SetTitle("layer(0-6)*6 + mod(0-5)");
hdiff->GetZaxis()->SetTitleOffset(1.8);
gStyle->SetPalette(kRainBow);
hdiff->SetMinimum(1);
hdiff->SetMaximum(10+(nstrips/ncn));
hdiff->Draw("COLZ");
for( uint i=0; i<lines.size(); i++)lines.at(i)->Draw("same");

char histname[400];
sprintf(histname, "hdiff.png");
c1->SaveAs(histname);


}
