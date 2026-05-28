using namespace std;

#include "KYtools.C"

#include "GDataEvent.hh"
#include "GGeometryTools.hh"
#include "GDataPoint.hh"
#include "GDataTrack.hh"
#include "GDataVertex.hh"

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


int main(int argc, char *argv[]){

GOptionParser* parser = GOptionParser::GetInstance();
parser->AddProgramDescription("Minimal Reproducable Example for Extracing Data from Reco Data");
parser->AddCommandLineOption<string>("in_path", "path to instrument data files", "./*", "i");
parser->AddCommandLineOption<int>("event", "specific event", 0, "e");
parser->ParseCommandLine(argc, argv);
parser->Parse();

string reco_path = parser->GetOption<string>("in_path");
int sp_event = parser->GetOption<int>("event");


cout << reco_path << endl;

char FilenameRoot[400];
sprintf(FilenameRoot,"%s*.root",reco_path.c_str());
cout << FilenameRoot << endl;

CEventRec* Event = new CEventRec(); //New reconstructed event
TChain * TreeRec = new TChain("TreeRec"); //New TreeRec Tchain object (this is new to me)
TreeRec->SetBranchAddress("Rec", &Event); //Set the branch address using Event (defined above)
TreeRec->Add(FilenameRoot);

//Prepare Reconstruction variable:
Crane::Reconstruction::TrackFit::GDataEvent * reco_data_event_ = new Crane::Reconstruction::TrackFit::GDataEvent();
TChain * TreeGReco = new TChain("TreeGReco");
TreeGReco->SetBranchAddress("FindPrimaryStarIterative", &reco_data_event_); //Set the branch address using Event (defined above)
TreeGReco->Add(FilenameRoot);

double tracker_ZBaricenter = 734; //mm
double zTolerance = 500; //mm
double yTolerance = 600; //mm

/*
//Prepare MC event
CEventMc* MCEvent = new CEventMc(); //New reconstructed event
TChain * TreeMC = new TChain("TreeMc"); //New TreeMC Tchain object (this is new to me)
TreeMC->SetBranchAddress("Mc", &MCEvent); //Set the branch address using Event (defined above)
TreeMC->Add(FilenameRoot);
*/

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
//for(unsigned int i = 0; i < TreeRec->GetEntries(); i+=MainLoopScaleFactor){
for(unsigned int i = sp_event; i < sp_event+1; i+=MainLoopScaleFactor){
    TreeRec->GetEntry(i);
    TreeGReco->GetEntry(i);
    //TreeMC->GetEntry(i);

    //cout << endl << "Event is " << i << endl;
    //cout << "Number of tracks " << Event->GetNTracks() << endl;

    //Energy deposition information for all of the tracks!

    for(uint t = 0; t < Event->GetNTracks(); t++){
        cout << "Track is " << t << endl;
        if(Event->GetTrack(t)->IsPrimary()) cout << "I am the Primary track! " << endl;
        for(uint isig=0; isig<Event->GetTrack(t)->GetEnergyDeposition().size(); isig++){
            cout << "Edep " << isig << " is " << Event->GetTrack(t)->GetEnergyDeposition(isig) << " at " << Event->GetTrack(t)->GetVolumeId(isig) << endl;
        }
    }

    //cout << "Event ID? " << Event->GetEventId() << endl;
    //cout << "Event Number? " << Event->GetEventNumber() << endl; //Gviz2D is for sure pulling Event number!
    //cout << "Reconstruction used: " << Event->GetActiveReconstruction() << endl;

/*
    //Searching for slowing down event in flight data.
    if(fabs(Event->GetPrimaryBeta()) > 0.9 && fabs(Event->GetPrimaryBeta()) < 1.0 && Event->GetNTracks() == 1 && -fabs(Event->GetPrimaryMomentumDirection().CosTheta()) > -0.75 && -fabs(Event->GetPrimaryMomentumDirection().CosTheta()) < -0.54){
        int UMBflag = 0;
        int CORflag = 0;
        int CBEtopflag = 0;
        int CBEbotflag = 0;
        int CBEsideflag = 0;
        int tkrflag = 0;
        double beta = Event->GetPrimaryBeta();

        //First iteration over events to check for TOF hits
        for(uint isig=0; isig<Event->GetTrack(0)->GetEnergyDeposition().size(); isig++){
        //for(uint isig=0; isig<Event->GetVolumeId().size(); isig++){
            //unsigned int VolumeId = Event->GetVolumeId().at(isig);
            unsigned int VolumeId = Event->GetTrack(0)->GetVolumeId(isig);
            if(volspec(VolumeId,0,2) == 20){
                if(GGeometryObject::GetTrackerLayer(VolumeId))tkrflag++;
                //if(GGeometryObject::GetTrackerLayer(VolumeId) < 4)tkrflag++; //Slow-stop
            }
            if(volspec(VolumeId,0,3) == 100)UMBflag++;
            if(volspec(VolumeId,0,3) == 110){ CBEtopflag++; }
            if(volspec(VolumeId,0,3) == 111){ CBEbotflag++;  }
            if(volspec(VolumeId,0,3) == 112 || volspec(VolumeId,0,3) == 113 || volspec(VolumeId,0,3) == 114 || volspec(VolumeId,0,3) == 115 || volspec(VolumeId,0,3) == 116)CBEsideflag++;
            if(volspec(VolumeId,0,3) == 102 || volspec(VolumeId,0,3) == 103 || volspec(VolumeId,0,3) == 104 || volspec(VolumeId,0,3) == 105 || volspec(VolumeId,0,3) == 106)CORflag++;
        }

        CTrackRec* pt = Event->GetPrimaryTrack();
		uint pt_index = 0;
        for( ; pt_index < Event->GetNTracks(); pt_index++) if( Event->GetTrack(pt_index)->IsPrimary() ) break;

        //Gr8 single track muons
        if(UMBflag > 0 && CBEtopflag > 0 && CBEbotflag > 0 && tkrflag > 4 && (pt->GetChi2()/pt->GetNdof()) < 2){cout << "Event " << i << " is unbelieveably rad " << endl;
            cout << "Event ID? " << Event->GetEventId() << endl;
            cout << "Event Number? " << Event->GetEventNumber() << endl;
            for(uint isig=0; isig<Event->GetTrack(0)->GetEnergyDeposition().size(); isig++){
                 unsigned int VolumeId = Event->GetTrack(0)->GetVolumeId(isig);
                 cout << "Energy deposition " << isig << " is " << Event->GetTrack(0)->GetEnergyDeposition(isig) << " at " << VolumeId << endl;
             }
        }


        //Stopping particles!
        if(UMBflag > 0 && CBEtopflag > 0 && CBEbotflag == 0 && tkrflag > 2 && tkrflag < 5 && CORflag == 0 && CBEsideflag == 0 &&  Event->GetTriggerVolumeId().size() < 3){
                cout << "Event " << i << " is a cool stopping(?) friend " << endl;
               	for(unsigned int k = 0; k < Event->GetTriggerVolumeId().size(); k++){
                    unsigned int LGVolumeId = Event->GetTriggerVolumeId().at(k);
                    cout << " LG Trigger VID " << k << ": " << LGVolumeId << endl;
                }
        }
    }*/

    //End stopping friend search

    //cout << endl << "All hits? " << endl;

    //All hits, remove on track requirement
    /*
    for(uint isig=0; isig<Event->GetVolumeId().size(); isig++){
        cout << "Hit " << isig << " VolumrId " << Event->GetVolumeId().at(isig) << " Edep " << Event->GetHitSeries().at(isig).GetTotalEnergyDeposition() << endl;
        cout << "Position " << Event->GetHitSeries().at(isig).GetPosition().X() << endl;
    }
    */

    /*
	cout << endl << "MC info " << endl;
	cout << "Number of tracks " << MCEvent->GetNTracks() << endl;
	vector<unsigned int> MCVolid;
	vector<double> MCEdep;

    for(uint t = 0; t < MCEvent->GetNTracks();t++){
        cout << "Track is " << t << endl;
        for(uint isig=0; isig<MCEvent->GetTrack(t)->GetEnergyDeposition().size(); isig++){
            unsigned int VolumeId  = MCEvent->GetTrack(t)->GetVolumeId(isig);
            //cout << "VID = " << VolumeId << endl;
            //cout << MCVolid.size() << endl;

            //int vsize = MCVolid.size();
            int index = -1;

            for (int k = 0; k < MCVolid.size(); k++) {
                if (MCVolid[k] == VolumeId) { //If volumeid is in the vector already, add the energy deposition to the vector counting edeps there
                    MCEdep[k] = MCEdep[k] + MCEvent->GetTrack(t)->GetEnergyDeposition(isig);
                    index = k;
                }
            }
            if(index == -1){ //If VolumeId is not in the vector, add it to the vector and add the energy deposition to the energy desposition vector
                //cout << "VID not found, adding" << endl;
                MCVolid.push_back(VolumeId);
                MCEdep.push_back(MCEvent->GetTrack(t)->GetEnergyDeposition(isig));
            }

            //cout << "MC: Hit " << isig << " Edep " << MCEvent->GetTrack(t)->GetEnergyDeposition(isig) << " at " << VolumeId << endl;
            //cout << "Method 2? " << Event->GetTrack(0)->GetEnergyDeposition().at(isig) << endl;
        }
    } //End loop over tracks
    */

    /*
    cout << endl;
    for (int k = 0; k < MCVolid.size(); k++) {
        cout << "VID " << MCVolid[k] << " total edep " << MCEdep[k] << endl;
    } //This seems to be a fine way of determining total edeps in the instrument from MC
    cout << endl;
    */

    /*
	if(MCEvent->GetNTracks() >= 0){ //Just output everything right now
	    //Mmk, got it! Of course the Gviz and Gviewer are code and code is readable! In the SimpleDet tools, the .cc files have the info :3
	    //If GetNTracks > 1, then there are tracks to iterate over!
		for(uint t = 0; t < MCEvent->GetNTracks();t++){
		    //What hits needed to call this track worthwhile? More than one non-zero energy deposition, right?
		    for(unsigned int isig = 0; isig < MCEvent->GetTrack(t)->GetEnergyDeposition().size(); isig++){
				if(MCEvent->GetTrack(t)->GetEnergyDeposition(isig) > 0){ inttrk++; }
			}
			//if(inttrk > 1){cool = 1;}
			//inttrk = 0;


		    if(cool){
				cout << "Track is " << t << endl;
				cout << "GetTrackId()? " << MCEvent->GetTrack(t)->GetTrackId()<< endl; //t is not the same as GetTrackId, interesting!
				    //Try to pick at some information?
				cout << "IsXray? " << MCEvent->GetTrack(t)->IsXray() << endl;
				cout << "GetParentId()? " << MCEvent->GetTrack(t)->GetParentId()<< endl;
				cout << "GetProcessType()? " << MCEvent->GetTrack(t)->GetProcessType()<< endl;
				//if(MCEvent->GetTrack(t)->GetProcessType() == "ELECTROMAGNETIC"){ cout << 'Hmm?' << endl;}
				cout << "GetPdg()? " << MCEvent->GetTrack(t)->GetPdg()<< endl;
				//cout << "IsMc()? " << MCEvent->GetTrack(t)->IsMc()<< endl;
				//cout << "IsXray? " << MCEvent->GetTrack(t)->IsXray() << endl;
				for(unsigned int isig = 0; isig < MCEvent->GetTrack(t)->GetEnergyDeposition().size(); isig++){
				    //Loop over each energy deposition on each track
				    cout << "Energy deposition " << MCEvent->GetTrack(t)->GetEnergyDeposition(isig) << " at Volume " << MCEvent->GetTrack(t)->GetVolumeId(isig) << endl;
								//The volume ID is weird? Multiple different energy deposits at the same volumeID? Maybe that's a simulations thing...?
					//So there's a ton of "tracks" but so many just have one or a bunch of Energy depositions = 0...
					//How can there be a track also with only one point too haha?
					cool = 1;
				}
			}
		}

	}*/

}

cout << endl << "I am done" << endl;
return 1;

}
