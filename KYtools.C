using namespace std;

#include <TColor.h>
#include <TProfile.h>
#include <TMath.h>
#include <TChain.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TLatex.h>
#include <TMinuit.h>
#include <vector>
#include <map>
#include <string>
#include <stdlib.h>
#include "TH1.h"
#include "TCanvas.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>


//FIXME: does this work on mac?
#include <sys/stat.h>

//#include "CRawTrk.hh"

#include "CEventMc.hh"
#include "CAnalysisManager.hh"
#include "GAnalysisIdentification.hh"
#include "GBasicTrigger.hh"
#include "GSimulationParameter.hh"
#include "GPreselection.hh"
#include "CraneConstants.hh"
#include "CraneLogging.hh"
#include "GPlottingTools.hh"
#include "CNet.hh"
#include "CBackpropagation.hh"

#include "GGeometry.hh"

#ifdef USE_BOOST_PROGRAM_OPTIONS
#include "GOptionParser.hh"
#include "GFileIO.hh"
#endif

using namespace Crane::Analysis;
namespace ca = Crane::Analysis;
namespace cl = Crane::Common;
//using Crane::Calibration;

//Function for rounding a double or float to some number of decimal points before making it a string.
string roundstr_d(double value, int precision){
    ostringstream oss;
    oss << std::fixed << std::setprecision(precision) << value;
    return oss.str();
}

//This function will take in a VolumeID number and then two other numbers that specify which part of the VolumeID you want to interrogate
//so ideally volspec(VolumeID,0,3) will give you the first three numbers of VolumeID which can tell you which CBE part it is.
//How to use, volspec(12345,1,4) outputs 234 as an integer for your comparing needs
int volspec(int volnum,int a, int b){
        stringstream ss;
        ss << volnum;
        return atoi(ss.str().substr(a, b).c_str());
}

//The two functions below aim to take the SimpleDet mod from 0-35 and translate that into a row and module
int getmod(int sdlayer, int sdmod){
        if(sdlayer % 2 == 0){
                return 5 - (sdmod % 6);
        }else{
                return 5 - (floor(sdmod/6));
        }
}

//This function takes the SimpleDet layer adn mod and outputs the row.
int getrow(int sdlayer, int sdmod){
        if(sdlayer % 2 == 0){
                return floor(sdmod/6);
        }else{
                return sdmod % 6;
        }
}

//This function aims to take layer, sddet, and sdstrip and translate that into a channel 0-31
int getch(int sdlayer, int sddet, int sdstrp){
	if(sdlayer % 2 == 0){
		if(sddet == 0){return sdstrp + 24;}
		if(sddet == 1){return 23-sdstrp;}
		if(sddet == 2){return sdstrp;}
		if(sddet == 3){return 15-sdstrp;}
		else{return -1;}
	}
	if(sdlayer % 2 != 0){
                if(sddet == 0){return 7 - sdstrp;}
                if(sddet == 1){return 31 - sdstrp;}
                if(sddet == 2){return sdstrp + 8;}
                if(sddet == 3){return sdstrp + 16;}
		else{return -1;}
	}
	else{return -1;}
}

//The function below plots TH1D's without all the blocks and blocks of text. More variables can be added if margins, etc... want to change.
void histplot1d(string ctitle, TH1D* h1, string title, string xtitle, string ytitle, string savename){
    TCanvas * c1 = new TCanvas(ctitle.c_str(), ctitle.c_str(), 200, 10, 900, 900);
    c1->SetLeftMargin(0.1);
    c1->SetRightMargin(0.1);
    c1->SetTopMargin(0.1);
    c1->SetBottomMargin(0.1);

    h1->SetTitle(title.c_str());
    h1->GetXaxis()->SetTitle(xtitle.c_str());
    h1->GetYaxis()->SetTitle(ytitle.c_str());
    h1->Draw();

    gPad->SetGridx(1);
    gPad->SetGridy(1);
    gPad->SetLogy(1);

    char histname[400];
    sprintf(histname, "%s.png",savename.c_str());
    c1->SaveAs(histname);
}

//The function below plots TH1D's without all the blocks and blocks of text. More variables can be added if margins, etc... want to change.
void histplot1f(string ctitle, TH1F* h1, string title, string xtitle, string ytitle, string savename){
    TCanvas * c1 = new TCanvas(ctitle.c_str(), ctitle.c_str(), 200, 10, 900, 900);
    c1->SetLeftMargin(0.1);
    c1->SetRightMargin(0.1);
    c1->SetTopMargin(0.1);
    c1->SetBottomMargin(0.1);

    h1->SetTitle(title.c_str());
    h1->GetXaxis()->SetTitle(xtitle.c_str());
    h1->GetYaxis()->SetTitle(ytitle.c_str());
    h1->Draw();

    gPad->SetGridx(1);
    gPad->SetGridy(1);
    gPad->SetLogy(1);

    char histname[400];
    sprintf(histname, "%s.png",savename.c_str());
    c1->SaveAs(histname);
}

void histplot2d(string ctitle, TH2D* h1, string title, string xtitle, string ytitle, string ztitle, string savename){
    TCanvas * c1 = new TCanvas(ctitle.c_str(), ctitle.c_str(), 200, 10, 900, 900);
    c1->SetLeftMargin(0.14);
    c1->SetRightMargin(0.16);
    c1->SetTopMargin(0.1);
    c1->SetBottomMargin(0.1);

    h1->SetTitle(title.c_str());
    h1->SetBit(TH1::kNoStats);
    h1->GetXaxis()->SetTitle(xtitle.c_str());
    h1->GetYaxis()->SetTitle(ytitle.c_str());
    h1->GetZaxis()->SetTitle(ztitle.c_str());

    h1->Draw("COLZ");
    gPad->SetLogz();

    char histname[400];
    sprintf(histname, "%s.png",savename.c_str());
    c1->SaveAs(histname);
}

void histplot2f(string ctitle, TH2F* h1, string title, string xtitle, string ytitle, string ztitle, string savename){
    TCanvas * c1 = new TCanvas(ctitle.c_str(), ctitle.c_str(), 200, 10, 900, 900);
    c1->SetLeftMargin(0.14);
    c1->SetRightMargin(0.16);
    c1->SetTopMargin(0.1);
    c1->SetBottomMargin(0.1);

    h1->SetTitle(title.c_str());
    h1->SetBit(TH1::kNoStats);
    h1->GetXaxis()->SetTitle(xtitle.c_str());
    h1->GetYaxis()->SetTitle(ytitle.c_str());
    h1->GetZaxis()->SetTitle(ztitle.c_str());

    h1->Draw("COLZ");
    gPad->SetLogz();

    char histname[400];
    sprintf(histname, "%s.png",savename.c_str());
    c1->SaveAs(histname);
}

void histplot2f_pts(string ctitle, TH2F* h1, TGraph * pts, string title, string xtitle, string ytitle, string ztitle, string savename){
    TCanvas * c1 = new TCanvas(ctitle.c_str(), ctitle.c_str(), 200, 10, 900, 900);
    c1->SetLeftMargin(0.14);
    c1->SetRightMargin(0.16);
    c1->SetTopMargin(0.1);
    c1->SetBottomMargin(0.1);

    h1->SetTitle(title.c_str());
    h1->SetBit(TH1::kNoStats);
    h1->GetXaxis()->SetTitle(xtitle.c_str());
    h1->GetYaxis()->SetTitle(ytitle.c_str());
    h1->GetZaxis()->SetTitle(ztitle.c_str());

    h1->Draw("COLZ");
    pts->Draw("P same");
    gPad->SetLogz();

    char histname[400];
    sprintf(histname, "%s.png",savename.c_str());
    c1->SaveAs(histname);
}


void histplot2f_pts_errs(string ctitle, TH2F* h1, TGraph * pts, TGraph * errs, string title, string xtitle, string ytitle, string ztitle, string savename){
    TCanvas * c1 = new TCanvas(ctitle.c_str(), ctitle.c_str(), 200, 10, 900, 900);
    c1->SetLeftMargin(0.14);
    c1->SetRightMargin(0.16);
    c1->SetTopMargin(0.1);
    c1->SetBottomMargin(0.1);

    h1->SetTitle(title.c_str());
    h1->SetBit(TH1::kNoStats);
    h1->GetXaxis()->SetTitle(xtitle.c_str());
    h1->GetYaxis()->SetTitle(ytitle.c_str());
    h1->GetZaxis()->SetTitle(ztitle.c_str());

    h1->Draw("COLZ");
    pts->Draw("P same");
    errs->Draw("P same");
    gPad->SetLogz();

    char histname[400];
    sprintf(histname, "%s.png",savename.c_str());
    c1->SaveAs(histname);
}

//Function that takes a 2D histogram, pts, and errors. Finds the max histogram bin in each column and puts it in pts
//Finds the FWHM and puts it in errs
//Let's see if this idea even works.

void label_2Dhisto(TH2F * h2, TGraph *pts, TGraph *errs) {
    int cols = h2->GetNbinsX();
    int rows = h2->GetNbinsY();

    for(int i = 1; i <= cols; i++){
         double maxVal = 0; //For each column, there will be a maximum value
         int FWHMlow = 0; //FWHM low bin
         int FWHMhi = 0; //FWHM high bin
         int maxBin = 0;

        for(int j = 0; j < rows; j++){ //Iterate over all of the bins and find the maximum bin for TOF and TKR
             double content = h2->GetBinContent(i, j);
             if (content > maxVal && j != 0){maxVal = content; maxBin = j;}
        }
        //cout << "Maxval " << maxVal << " at max bin " << maxBin << endl;

        for(int j = 0; j < maxBin; j++){ //Find low bin by starting at max bin of TKR and step down
            double content = h2->GetBinContent(i, maxBin-j);
            if (content < maxVal/2){ FWHMlow = maxBin-j; break;}
        }

        for(int j = maxBin; j < rows - maxBin; j++){ //Find high bin by starting at max bin of TKR and step down
            double content = h2->GetBinContent(i, j);
            if (content < maxVal/2){ FWHMhi = j; break;}
        }

        pts->SetPoint(i, h2->GetXaxis()->GetBinCenter(i), h2->GetYaxis()->GetBinCenter(maxBin));
        errs->SetPoint(i, h2->GetXaxis()->GetBinCenter(i), h2->GetYaxis()->GetBinCenter(FWHMhi));
        errs->SetPoint(i+rows, h2->GetXaxis()->GetBinCenter(i), h2->GetYaxis()->GetBinCenter(FWHMlow));

    }

    // Set pts appearance (do a different set for the FWHMs)
    pts->SetMarkerStyle(20); // Solid circle
    pts->SetMarkerColor(kBlack);
    pts->SetMarkerSize(1.5);

    errs->SetMarkerStyle(3); // Stars
    errs->SetMarkerColor(kBlack);
    errs->SetMarkerSize(1.5);
}
