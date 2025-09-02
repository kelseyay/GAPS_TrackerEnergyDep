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
//To run, do ./Occupancy -i /home/kelsey/simulations/simdat/simnew/mu-_gaps_triggerlevel1_FTFP_BERT_1744342800_rec.root
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

int main(int argc, char *argv[]){

GOptionParser* parser = GOptionParser::GetInstance();
parser->AddProgramDescription("Minimal Reproducable Example for Extracing Data from Reco Data");
parser->AddCommandLineOption<string>("in_path", "path to instrument data files", "./*", "i");
parser->ParseCommandLine(argc, argv);
parser->Parse();

string reco_path = parser->GetOption<string>("in_path");

cout << reco_path << endl;

char FilenameRoot[400];
sprintf(FilenameRoot,"%s*.root",reco_path.c_str());
cout << FilenameRoot << endl;

CEventRec* Event = new CEventRec(); //New reconstructed event
TChain * TreeRec = new TChain("TreeRec"); //New TreeRec Tchain object (this is new to me)
TreeRec->SetBranchAddress("Rec", &Event); //Set the branch address using Event (defined above)
TreeRec->Add(FilenameRoot);

int MainLoopScaleFactor = 1; //Set this number to scale the step size. Larger means runs faster and fewer events
double TrackerCut = 0.4; //Threshold for an energy deposition to be considered a hit
double TofCutLow = 0;
double xlow = 0.3; //Hist range
double xhigh = 6; //Hist range
double coshigh = 0; //0.92; //0.995; //0.92 //0.54 is the highest angle that can hit UMB, CBEtop, CBEbot
double coslow = 1; //0.62 //0.8
double eventbetacut = 0.2; //Cut that is applied to all events
double betacut = 0.9; //Separation of quandrants in the beta plot
const Int_t NBins = 50;

//Prepare cuts:
map<int, unsigned int> TofIndexVolumeIdMap;

//Not a proper tell me anymore lol
TH1F * h;
h = new TH1F ("Edep l", "Edep", NBins, xlow,xhigh);

//For Plotting purposes
ca::GPlottingTools Plotting;
char text[400]; //This variable is used later to name the plots

//All of the plots are declared here
//("Title",Number of bins,xmin,xmax,"xlabel","ylabel",ymin,ymax)

//Now we can go over the loop
TreeRec->GetEntry(0);

cout << "Total Number of events / Mainscale Factor = " << TreeRec->GetEntries()/MainLoopScaleFactor << endl;

//Using i to loop over every event in the tree
for(unsigned int i = 0; i < TreeRec->GetEntries(); i+=MainLoopScaleFactor){
        TreeRec->GetEntry(i);

        //Cuts are implemented in this chunk:
        if(Event->GetNTracks() == 1){  //First select the single track event

                //cout << "Single Track Event " << endl;
                bool Umbflag = 0;
                bool CBEtopflag = 0;
                bool CBEbotflag = 0;

                CTrackRec* pt = Event->GetPrimaryTrack();
                uint pt_index = 0;
                for( ; pt_index < Event->GetNTracks(); pt_index++) if( Event->GetTrack(pt_index)->IsPrimary() ) break;
                //Apply event-level cuts
                if(pt != nullptr && -fabs(Event->GetPrimaryMomentumDirection().CosTheta()) > -coslow && -fabs(Event->GetPrimaryMomentumDirection().CosTheta()) < -coshigh && Event->GetPrimaryBeta()*Event->GetPrimaryMomentumDirection()[2] < 0 && fabs(Event->GetPrimaryBeta()) > eventbetacut){

                        cout << endl << "Event is " << i << endl;

                        //First loop over the event for your needed flags and variables
                        for(unsigned int isig = 0; isig < Event->GetTrack(0)->GetEnergyDeposition().size(); isig++){
                                unsigned int VolumeId  = Event->GetTrack(0)->GetVolumeId(isig); //Event->GetVolumeId().at(isig); //Check the VolumeId of the event
                                //cout << "Volume ID is " << VolumeId << " Energy Dep is " << Event->GetTrack(0)->GetEnergyDeposition(isig) << " MeV "<< endl;
                                if(volspec(VolumeId,0,3) == 100){ Umbflag = 1;} // cout << "UMB hit!" << endl;
                                if(volspec(VolumeId,0,3) == 110) {CBEtopflag = 1;}// cout << "CBE top hit!" << endl;
                                if(volspec(VolumeId,0,3) == 111) {CBEbotflag = 1;}// cout << "CBE bot hit!" << endl;
                                if(GGeometryObject::IsTrackerVolume(VolumeId)){h->Fill(Event->GetTrack(0)->GetEnergyDeposition(isig));}
                                if(GGeometryObject::IsTrackerVolume(VolumeId)){

                              		int layer = GGeometryObject::GetTrackerLayer(VolumeId);
                                    int sdmod = GGeometryObject::GetLayerModule(VolumeId);
						            int det = GGeometryObject::GetModuleDetector(VolumeId);
						            int sdstrip = GGeometryObject::GetDetectorStrip(VolumeId);

						            int row = getrow(layer,sdmod);
						            int mod = getmod(layer,sdmod);
						            int strip = getch(layer, det, sdstrip);
                                    cout << "Volume ID is " << VolumeId << " lrms is : " << layer << row << mod << " " << strip << " Energy Dep is " << Event->GetTrack(0)->GetEnergyDeposition(isig) << " MeV "<< endl;
                                }
                        }

                        //If the desired flags are checked, proceed to fill the relevant histograms.
                        if(Umbflag && CBEtopflag && CBEbotflag /*&& (pt->GetChi2()/pt->GetNdof()) < 3.2*/){
                                //cout << "Event is " << i << " passes the cuts" << endl;

                                //Loop over the events on track 0
                                for(unsigned int k = 0; k < Event->GetTrack(0)->GetEnergyDeposition().size(); k++){
                                unsigned int VolumeId = Event->GetTrack(0)->GetVolumeId(k); //For each hit, check the volume ID
                                        //Check if the volume is a TOF volume and make sure you pass the low energy criteria


                                        }
                                }

                }
        }
}


TCanvas * c1 = new TCanvas("c1", "c1", 200, 10, 900, 900);
c1->SetLeftMargin(0.1);
c1->SetRightMargin(0.16);
c1->SetTopMargin(0.1);
c1->SetBottomMargin(0.1);

gPad->SetGridx(1);
gPad->SetGridy(1);
gPad->SetLogy(1);

h->SetTitle("Full Tracker Energy Depositions (MeV)");
h->GetXaxis()->SetTitle("Energy Deposition (no Cos factor) MeV");
h->GetYaxis()->SetTitle("Number of Hits");
h->GetYaxis()->SetRangeUser(100, 500000);

h->SetLineColor(2);
//h->Draw("hist");

//sprintf(text, "EdepTest.png");
//c1->SaveAs(text);

cout << endl << "I am done" << endl;
return 1;

}
