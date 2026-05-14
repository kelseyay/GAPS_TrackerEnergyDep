using namespace std;

#include "KYtools.C"

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

int main(int argc, char *argv[]){

GOptionParser* parser = GOptionParser::GetInstance();
parser->AddProgramDescription("Minimal Reproducable Example for Extracing Data from Reco Data");
parser->AddCommandLineOption<string>("in_path", "path to instrument data files", "./*", "i");
parser->AddCommandLineOption<string>("out_file", "name of output root file", "out.root", "o");
parser->AddCommandLineOption<double>("beta_low", "low Beta Cut",0.8,"l");
parser->AddCommandLineOption<double>("beta_high", "upper Beta Cut",1,"u");
parser->AddCommandLineOption<int>("MainloopScale", "Main loop scale factor",1,"m");
parser->ParseCommandLine(argc, argv);
parser->Parse();

double betacut = parser->GetOption<double>("beta_low");
if(betacut <= 0 || betacut >=1){ betacut = 0.8; cout << "Error with low beta choice. Setting Beta low to 0.8" << endl; }
cout << "beta cut = " << betacut << endl;

double betahigh = parser->GetOption<double>("beta_high");
if(betahigh <= 0 || betahigh >=2 || betahigh < betacut){ betahigh = 1; cout << "Error with high/low beta choice! Setting Beta upper to 1" << endl; }
cout << "beta high = " << betahigh << endl;

int MainLoopScaleFactor = parser->GetOption<int>("MainloopScale");

//int MainLoopScaleFactor = 1; //Set this number to scale the step size. Larger means runs faster and fewer events

string reco_path = parser->GetOption<string>("in_path");
string out_path = parser->GetOption<string>("out_file");
cout << reco_path << endl;

if(out_path != "" && out_path[out_path.length()-1] != '/' ){ cout <<  "out path no slash!" << endl; out_path = out_path + '/'; }

char FilenameRoot[400];
sprintf(FilenameRoot,"%s*.root",reco_path.c_str());
cout << FilenameRoot << endl;

//Prepare reconstructed event
CEventRec* Event = new CEventRec(); //New reconstructed event
TChain * TreeRec = new TChain("TreeRec"); //New TreeRec Tchain object (this is new to me)
TreeRec->SetBranchAddress("Rec", &Event); //Set the branch address using Event (defined above)
TreeRec->Add(FilenameRoot);

//Prepare MC event
CEventMc* MCEvent = new CEventMc(); //New reconstructed event
TChain * TreeMC = new TChain("TreeMc"); //New TreeMC Tchain object (this is new to me)
TreeMC->SetBranchAddress("Mc", &MCEvent); //Set the branch address using Event (defined above)
TreeMC->Add(FilenameRoot);

double TofCutLow = 0; //No low Tof cut right now
double TrackerCut = 0.3; //Threshold for an energy deposition to be considered a hit

double xlow = 0.3; //Low range for histogram MeV
double xhigh = 20; //High range for histogram MeV

double coshigh = 0.54; //0.995; //0.92 //0.54 is the highest angle that can hit UMB, CBEtop, CBEbot
double coslow = 1; //0.62 //0.8
double fitlow = 0.5;
double fithigh = 3.5;
const Int_t NBins = 50;

//Full tracker histogram range
double mpvmin = 0.66;
double mpvmax = 0.75;

//From fitted peaks from flight! Check over!
double TofAngleCorrectedMip = 0.82;
double TrackerAngleCorrectedMip = 0.57;

int alphactr = 0;
int pctr = 0;
int muctr = 0;

const int nlayers = 5;

//Trying simulated:
//double TofAngleCorrectedMip = 1;
//double TrackerAngleCorrectedMip = 1;

//Prepare cuts:
map<int, unsigned int> TofIndexVolumeIdMap;

//For Plotting purposes
ca::GPlottingTools Plotting;
char text[400]; //This variable is used later to name the plots

TH1F * HUMBHit = new TH1F ("h0", ("Edep UMB: Beta " + roundstr_d(betacut,2) + " - " + roundstr_d(betahigh,2) ).c_str(), NBins, xlow,xhigh);
TH1F * HCBE_topHit = new TH1F ("h1", ("Edep CBE_top: Beta " + roundstr_d(betacut,2) + " - " + roundstr_d(betahigh,2) ).c_str(), NBins, xlow,xhigh);

TH1F * HLHit[nlayers];
for(int l = 0;l<nlayers;l++){
    HLHit[l] = new TH1F (TString::Format("h0_l%i",l), ("Edep L" + to_string(l) + ": Beta " + roundstr_d(betacut,2) + " - " + roundstr_d(betahigh,2) ).c_str(), NBins, xlow,xhigh);
}

//How many entries:
cout << "Total Number of events / Mainscale Factor = " << TreeRec->GetEntries()/MainLoopScaleFactor << endl;
cout << "Beta high " << betahigh << endl;
cout << "Beta cut " << betacut << endl;
cout << "MPV range on 2D Full Tracker will be " << mpvmin << " - " << mpvmax << endl;

//Now we can go over the loop
TreeRec->GetEntry(0);
TreeMC->GetEntry(0);

cout << "Total Number of events / Mainscale Factor = " << TreeRec->GetEntries()/MainLoopScaleFactor << endl;

//Using i to loop over every event in the tree
//for(unsigned int i = 0; i < 1000; i+=MainLoopScaleFactor){
//for(unsigned int i = 0; i < TreeRec->GetEntries(); i+=MainLoopScaleFactor){
for(unsigned int i = 0; i < TreeRec->GetEntries()/MainLoopScaleFactor; i++){
    TreeRec->GetEntry(i);
    TreeMC->GetEntry(i);

	//Cuts are implemented in this chunk:
	if(Event->GetNTracks() == 1){  //First select the single track event
		bool Umbflag = 0;
		bool CBEtopflag = 0;
		bool CBEbotflag = 0;

		if( ((int)i % (int)ceil(TreeRec->GetEntries()/(MainLoopScaleFactor*10))) == 0){
		    cout << "Event number " << i << endl;
		}

		CTrackRec* pt = Event->GetPrimaryTrack();
		uint pt_index = 0;
      	        for( ; pt_index < Event->GetNTracks(); pt_index++) if( Event->GetTrack(pt_index)->IsPrimary() ) break;

		//Note downwards beta enforced by Event->GetPrimaryBeta() (should be positive) multiplied by Event->GetPrimaryMomentumDirection()[2] (z trajectory of particle)
		if(pt != nullptr && (pt->GetChi2()/pt->GetNdof()) < 3.2 && -fabs(MCEvent->GetPrimaryMomentumDirection().CosTheta()) > -coslow && -fabs(MCEvent->GetPrimaryMomentumDirection().CosTheta()) < -coshigh && MCEvent->GetPrimaryBeta()*MCEvent->GetPrimaryMomentumDirection()[2] < 0 && fabs(MCEvent->GetPrimaryBeta()) >  betacut && fabs(MCEvent->GetPrimaryBeta()) <  betahigh ){
		//if(pt != nullptr && -fabs(Event->GetPrimaryMomentumDirection().CosTheta()) > -coslow && -fabs(Event->GetPrimaryMomentumDirection().CosTheta()) < -coshigh && Event->GetPrimaryBeta()*MCEvent->GetPrimaryMomentumDirection()[2] < 0 && fabs(Event->GetPrimaryBeta()) >  betacut && fabs(Event->GetPrimaryBeta()) <  betahigh ){
		//if(pt != nullptr && -fabs(Event->GetPrimaryMomentumDirectionGenerated().CosTheta()) > -coslow && -fabs(Event->GetPrimaryMomentumDirectionGenerated().CosTheta()) < -coshigh && Event->GetPrimaryBetaGenerated()*Event->GetPrimaryMomentumDirectionGenerated()[2] < 0 && fabs(Event->GetPrimaryBetaGenerated()) >  betacut && fabs(Event->GetPrimaryBetaGenerated()) <  betahigh ){
		    //cout << "Event number " << i << endl;
			//-----------EVENT LEVEL CUT APPLIED

			//First iteration over event for flags and EnergyDepositionMip vector filling

			for(uint isig=0; isig<Event->GetTrack(0)->GetEnergyDeposition().size(); isig++){
                unsigned int VolumeId  = Event->GetTrack(0)->GetVolumeId(isig); //Check the VolumeId of the event

                if(volspec(VolumeId,0,3) == 100){
                    HUMBHit->Fill(Event->GetTrack(0)->GetEnergyDeposition(isig)*fabs(Event->GetPrimaryMomentumDirection().CosTheta()));
                    Umbflag = 1;
                    //cout << "UMB HIT! " <<endl;
                }

                if(volspec(VolumeId,0,3) == 110) {
                    HCBE_topHit->Fill(Event->GetTrack(0)->GetEnergyDeposition(isig)*fabs(Event->GetPrimaryMomentumDirection().CosTheta()));
                    CBEtopflag = 1;
                    //cout << "CBETOP HIT! " << endl;
                }
                //if(volspec(VolumeId,0,2) == 10 && volspec(VolumeId,2,1) != 0 ){ } //cout << "COR hit! " << endl;
                //if(volspec(VolumeId,0,3) == 111) { CBEbotflag = 1; }// cout << "CBE bot hit!" << endl;

			}

			//Do we want to only run this on certain tracks? Yeah probably. Can remove the TOF flags. Just run on whatever lol.
			if(Umbflag && CBEtopflag){
			    //cout << "UMB AND CBETOP HIT " << endl;
				//cout << "Event number " << i << " passes the cuts!" << endl;
				//cout << "Particle species " << MCEvent->GetTrack(0)->GetPdg() << endl;
				if(MCEvent->GetTrack(0)->GetPdg() == 1000020040){
				    alphactr++;
				}
				if(MCEvent->GetTrack(0)->GetPdg() == 2212){
				    pctr++;
				}
				if(MCEvent->GetTrack(0)->GetPdg() == 13){
				    muctr++;
				}

				for(uint isig=0; isig<Event->GetTrack(0)->GetEnergyDeposition().size(); isig++){
					unsigned int VolumeId  = Event->GetTrack(0)->GetVolumeId(isig);
					if(GGeometryObject::IsTrackerVolume(VolumeId) && Event->GetTrack(0)->GetEnergyDeposition(isig) > TrackerCut){
					    //cout << "TRACKER HIT! Layer " <<  GGeometryObject::GetTrackerLayer(VolumeId) << endl;
						int layer = GGeometryObject::GetTrackerLayer(VolumeId);
						//Fill histos! Woot woot
						if(layer < nlayers){ HLHit[layer]->Fill((Event->GetTrack(0)->GetEnergyDeposition(isig)) * fabs(Event->GetPrimaryMomentumDirection().CosTheta()) ); }

					} //Closed bracket for Tracker volume and tracker cutoff

				} //Closed bracket for iteration over event with TOF cuts

			} //Closed bracket for if statement for cuts


			//-----------EVENT LEVEL CUTS END


		} //Closed bracket for event level cut

	} //Closed bracket for single track cut

}  //Closed bracket for iteration through tree events, move on to the next event i



//Histogram section
//--------------------------------------
//histplot1f(string ctitle, TH1F* h1, string title, string xtitle, string ytitle, string savename)
histplot1f("c1",HUMBHit,"Edep UMB: Beta " +  roundstr_d(betacut,2)  + " - " + roundstr_d(betahigh,2), "Energy Deposition x Cos(theta)","NEntries",out_path + "EdepUMB" + roundstr_d(betacut,2)  + "to" + roundstr_d(betahigh,2) );
histplot1f("c2",HCBE_topHit,"Edep CBE_top: Beta " +  roundstr_d(betacut,2)  + " - " + roundstr_d(betahigh,2), "Energy Deposition x Cos(theta)","NEntries",out_path + "EdepCBE_top"+ roundstr_d(betacut,2)  + "to" + roundstr_d(betahigh,2) );

for(int j = 0;j < nlayers;j++){
    histplot1f("c"+to_string(j+3),HLHit[j], "Edep L" + to_string(j) + ": Beta " + roundstr_d(betacut,2)  + " - " + roundstr_d(betahigh,2) , "Energy Deposition x Cos(theta)","NEntries", out_path + "EdepL"+ to_string(j) + "B" + roundstr_d(betacut,2) + "to" + roundstr_d(betahigh,2)  );
}

cout << endl << "I am done" << endl;

return 1;

}
