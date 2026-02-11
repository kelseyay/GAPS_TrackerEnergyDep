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
//To run, do ./EdepHist -i /home/kelsey/simulations/simdat/simnew/mu-_gaps_triggerlevel1_FTFP_BERT_1744342800_rec.root
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
parser->AddCommandLineOption<string>("out_file", "name of output root file", "out.root", "o");
parser->ParseCommandLine(argc, argv);
parser->Parse();

string reco_path = parser->GetOption<string>("in_path");
string out_path = parser->GetOption<string>("out_file");
cout << reco_path << endl;
//sprintf(FilenameRoot,"%s/%s*.root",argv[1], argv[2]);
//cout << argv[1] << endl;

char FilenameRoot[400];
sprintf(FilenameRoot,"%s*.root",reco_path.c_str());
//sprintf(FilenameRoot,reco_path.c_str()); //210 simu data on my computer!
cout << FilenameRoot << endl;

CEventRec* Event = new CEventRec(); //New reconstructed event
TChain * TreeRec = new TChain("TreeRec"); //New TreeRec Tchain object (this is new to me)
TreeRec->SetBranchAddress("Rec", &Event); //Set the branch address using Event (defined above)
TreeRec->Add(FilenameRoot);

int MainLoopScaleFactor = 1; //Set this number to scale the step size. Larger means runs faster and fewer events
double TrackerCut = 0.3; //Threshold for an energy deposition to be considered a hit

double xlow = 0.3; //Low range for histogram MeV
double xhigh = 5; //High range for histogram MeV

double coshigh = 0.54; //0.995; //0.92 //0.54 is the highest angle that can hit UMB, CBEtop, CBEbot
double coslow = 1; //0.62 //0.8
double fitlow = 0.5;
double fithigh = 3.5;
const Int_t NBins = 50;
double betahigh = 0.75;
double betacut = 0.5; //Currently we're only doing a beta > 0 cutoff for real data. Beta > 0.8 recommended for sim

//Full tracker histogram range
double mpvmin = 0.66;
double mpvmax = 0.75;

TH1F * hedep;
hedep = new TH1F ("h0", ("Edep l Beta " + to_string(betacut) + " - " + to_string(betahigh) ).c_str(), NBins, xlow,xhigh);

//Prepare cuts:
map<int, unsigned int> TofIndexVolumeIdMap;

//For Plotting purposes
ca::GPlottingTools Plotting;
char text[400]; //This variable is used later to name the plots

//How many entries:
cout << "Total Number of events / Mainscale Factor = " << TreeRec->GetEntries()/MainLoopScaleFactor << endl;
cout << "Beta high " << betahigh << endl;
cout << "Beta cut " << betacut << endl;
cout << "MPV range on 2D Full Tracker will be " << mpvmin << " - " << mpvmax << endl;

//Now we can go over the loop
TreeRec->GetEntry(0);

cout << "Total Number of events / Mainscale Factor = " << TreeRec->GetEntries()/MainLoopScaleFactor << endl;

//Using i to loop over every event in the tree
//for(unsigned int i = 0; i < 100; i+=MainLoopScaleFactor){
for(unsigned int i = 0; i < TreeRec->GetEntries(); i+=MainLoopScaleFactor){
    TreeRec->GetEntry(i);

	//Cuts are implemented in this chunk:
	if(Event->GetNTracks() == 1){  //First select the single track event
		bool Umbflag = 0;
		bool CBEtopflag = 0;
		bool CBEbotflag = 0;

		if( ((int)i % (int)ceil(TreeRec->GetEntries()/10)) == 0){
		    cout << "Event number " << i << endl;
		}

		CTrackRec* pt = Event->GetPrimaryTrack();
		uint pt_index = 0;
      	        for( ; pt_index < Event->GetNTracks(); pt_index++) if( Event->GetTrack(pt_index)->IsPrimary() ) break;

		//Note downwards beta enforced by Event->GetPrimaryBeta() (should be positive) multiplied by Event->GetPrimaryMomentumDirection()[2] (z trajectory of particle)
		if(pt != nullptr && -fabs(Event->GetPrimaryMomentumDirection().CosTheta()) > -coslow && -fabs(Event->GetPrimaryMomentumDirection().CosTheta()) < -coshigh && Event->GetPrimaryBeta()*Event->GetPrimaryMomentumDirection()[2] < 0 && fabs(Event->GetPrimaryBeta()) >  betacut && fabs(Event->GetPrimaryBeta()) <  betahigh ){
			//cout << "Event is " << i << endl;

			//-----------EVENT LEVEL CUT APPLIED

			for(uint isig=0; isig<Event->GetTrack(0)->GetEnergyDeposition().size(); isig++){
                unsigned int VolumeId  = Event->GetTrack(0)->GetVolumeId(isig); //Check the VolumeId of the event
                if(volspec(VolumeId,0,3) == 100){ Umbflag = 1;} // cout << "UMB hit!" <<endl ;
                if(volspec(VolumeId,0,3) == 110) {CBEtopflag = 1;}// cout << "CBE top hit!" << endl;
                if(volspec(VolumeId,0,3) == 111) {CBEbotflag = 1;}// cout << "CBE bot hit!" << endl;
            }

			if(Umbflag && CBEtopflag /*&& CBEbotflag*/ && (pt->GetChi2()/pt->GetNdof()) < 3.2 ){
				//cout << "Event number " << i << " passes the cuts!" << endl;
				for(uint isig=0; isig<Event->GetTrack(0)->GetEnergyDeposition().size(); isig++){
					unsigned int VolumeId  = Event->GetTrack(0)->GetVolumeId(isig);
					if(GGeometryObject::IsTrackerVolume(VolumeId) && Event->GetTrack(0)->GetEnergyDeposition(isig) > TrackerCut){

							hedep->Fill(  (Event->GetTrack(0)->GetEnergyDeposition(isig)*fabs(Event->GetPrimaryMomentumDirection().CosTheta()))   );

					} //Closed bracket for Tracker volume and tracker cutoff

				} //Closed bracket for iteration over event with TOF cuts

			} //Closed bracket for if statement for cuts


			//-----------EVENT LEVEL CUTS END


		} //Closed bracket for event level cut

	} //Closed bracket for single track cut

}  //Closed bracket for iteration through tree events, move on to the next event i



//Histogram section
//--------------------------------------

//Histogram for NEntries at a strip level
TCanvas * c1 = new TCanvas("c2", "c2", 200, 10, 900, 900);
c1->SetLeftMargin(0.1);
c1->SetRightMargin(0.16);
c1->SetTopMargin(0.1);
c1->SetBottomMargin(0.1);

hedep->GetXaxis()->SetTitle("row(0-6)*32 + strp (0-31)");
hedep->GetYaxis()->SetTitle("layer(0-7)*6 + mod(0-6)");

gPad->SetGridx(1);
gPad->SetGridy(1);
gPad->SetLogy(1);

hedep->Draw();

string title = "OneHistFullTracker";
char histname[400];
sprintf(histname, "%s.root",title.c_str());
c1->SaveAs(histname);
sprintf(histname, "%s.png",title.c_str());
c1->SaveAs(histname);

//--------------------------------------

cout << endl << "I am done" << endl;

return 1;

}
