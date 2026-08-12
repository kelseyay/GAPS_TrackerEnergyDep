//Maybe write in some options so you can ask what to tell me lol

using namespace std;

#include "KYtools.C"

#include "GDataEvent.hh"
#include "GGeometryTools.hh"
#include "GDataPoint.hh"
#include "GDataTrack.hh"
#include "GDataVertex.hh"

using namespace Crane::Analysis;
namespace ca = Crane::Analysis;
namespace cl = Crane::Common;
//using Crane::Calibration;


int main(int argc, char *argv[]){

GOptionParser* parser = GOptionParser::GetInstance();
parser->AddProgramDescription("Minimal Reproducable Example for Extracing Data from Reco Data");
parser->AddCommandLineOption<string>("in_path", "path to instrument data files", "./*", "i");
parser->AddCommandLineOption<int>("event", "specific event", 0, "e");
parser->AddCommandLineOption<bool>("MC_truth", "specific event", 0, "m");
parser->AddCommandLineOption<bool>("dEdx", "dEdx information", 0, "x");
parser->AddCommandLineOption<bool>("Edep", "Just energy depositions", 0, "p");
parser->AddCommandLineOption<bool>("low_gain", "Low gain hit info", 0, "l");
parser->AddCommandLineOption<bool>("reco_info", "reconstruction info", 0, "r");
parser->AddCommandLineOption<string>("out_file", "name of output root file", "", "o");
parser->ParseCommandLine(argc, argv);
parser->Parse();

string reco_path = parser->GetOption<string>("in_path");
bool RECO = parser->GetOption<bool>("reco_info");
bool EDEP = parser->GetOption<bool>("Edep");
bool MC = parser->GetOption<bool>("MC_truth");
bool LG = parser->GetOption<bool>("low_gain");
bool DEDX = parser->GetOption<bool>("dEdx");
int sp_event = parser->GetOption<int>("event");

string out_path = parser->GetOption<string>("out_file");
cout << "out path: " << out_path << endl;
//cout << "out path last string " << out_path[out_path.length()-1] << endl;

if(out_path != "" && out_path[out_path.length()-1] != '/' ){ cout <<  "out path no slash! Adding! " << endl; out_path = out_path + '/'; }


string compare = ".root";

//str1.compare(0, 5, str2) == 0
cout << "reco_path: " << reco_path << endl;
//Compare == 0 if they are the same,
if(reco_path.compare(reco_path.length()-5,reco_path.length(),compare) == 0){ cout << ".root at the end of the reco path! Deleting!" << endl; reco_path = reco_path.substr(0,reco_path.length()-5); }

//if(reco_path.compare(reco_path.length()-5,reco_path.length(),"*.root") == 1){ cout << ".root at the end of the reco path! Deleting!" << endl; }

char FilenameRoot[400];
sprintf(FilenameRoot,"%s*.root",reco_path.c_str());
cout << FilenameRoot << endl;

CEventRec* Event = new CEventRec(); //New reconstructed event
TChain * TreeRec = new TChain("TreeRec"); //New TreeRec Tchain object (this is new to me)
TreeRec->SetBranchAddress("Rec", &Event); //Set the branch address using Event (defined above)
TreeRec->Add(FilenameRoot);

//Perfect! I wanted to output the first root file if a whole bunch were added with a star and this works!
TreeRec->LoadTree(0);
cout << "File 0 is: " << TreeRec->GetCurrentFile()->GetName() << endl;

//Prepare Reconstruction variable:
Crane::Reconstruction::TrackFit::GDataEvent * reco_data_event_ = new Crane::Reconstruction::TrackFit::GDataEvent();
TChain * TreeGReco = new TChain("TreeGReco");
TreeGReco->SetBranchAddress("FindPrimaryStarIterative", &reco_data_event_); //Set the branch address using Event (defined above)
TreeGReco->Add(FilenameRoot);

double tracker_ZBaricenter = 734; //mm
double zTolerance = 500; //mm
double yTolerance = 600; //mm

//Prepare MC event
CEventMc* MCEvent = new CEventMc(); //New reconstructed event
TChain * TreeMC = new TChain("TreeMc"); //New TreeMC Tchain object (this is new to me)
TreeMC->SetBranchAddress("Mc", &MCEvent); //Set the branch address using Event (defined above)
TreeMC->Add(FilenameRoot);

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

float Ztof = 5.574;
float Atof = 10.3;
float Ltof = 0.635;
float rtof = 1.032;
//float Gtofnew = tf; //1/0.9; //(1/0.71);
//cout << "sqrt(Gtofnew) = " << sqrt(tf) << endl;

float Ztkr = 14;
float Atkr = 28.3;
float Ltkr = 0.22; //GAPS Tracker //0.25; //MC
float rtkr = 2.33;
//float Gtkrnew = tk; //1; // (1/0.83);
//cout << "sqrt(Gtkrnew) = " << sqrt(tk) << endl;

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

int ev_prev = 0;

cout << "Total Number of events / Mainscale Factor = " << TreeRec->GetEntries()/MainLoopScaleFactor << endl;

//Using i to loop over every event in the tree
for(unsigned int i = 0; i < 100; i+=MainLoopScaleFactor){
//for(unsigned int i = 0; i < TreeRec->GetEntries(); i+=MainLoopScaleFactor){
//for(unsigned int i = sp_event; i < sp_event+1; i+=MainLoopScaleFactor){
    TreeRec->GetEntry(i);
    TreeGReco->GetEntry(i);
    TreeMC->GetEntry(i);

    if(Event->GetEventNumber() <= ev_prev){
        cout << "Event is " << i << "Event id is " << Event->GetEventNumber() << endl;
        cout << "Previous event id is " << ev_prev << endl;
    }

    ev_prev = Event->GetEventNumber();

    //cout << endl << "Event is " << i << endl;
    //cout << "Number of tracks " << Event->GetNTracks() << endl;

    //Energy deposition information for all of the tracks!


    //Tell me about the LG hits!
    if(LG){
        for(unsigned int k = 0; k < Event->GetTriggerVolumeId().size(); k++){
            unsigned int VolumeId = Event->GetTriggerVolumeId().at(k);
            cout << "LG Hit " << k << " at Volid " << VolumeId << endl;
        }
    }

    //Tell me about the energy depositions of all the tracks
    /*
    for(uint t = 0; t < Event->GetNTracks(); t++){
        cout << "Track is " << t << endl;
        if(Event->GetTrack(t)->IsPrimary()) cout << "I am the Primary track! " << endl;
        for(uint isig=0; isig<Event->GetTrack(t)->GetEnergyDeposition().size(); isig++){
            cout << "Edep " << isig << " is " << Event->GetTrack(t)->GetEnergyDeposition(isig) << " at " << Event->GetTrack(t)->GetVolumeId(isig) << endl;
        }
    }*/

    //Reco information (which reconstruction used?)
    if(RECO){
        cout << "Trigger source of event: " << (int)Event->GetTriggerSources().at(0) << endl;
        cout << "Reconstruction used: " << Event->GetActiveReconstruction() << endl;
    }

    if(DEDX || EDEP){

        //Tell me about the dE/dx of all tracks where dE/dx = Energy deposit / (rho_tkr/tof * length traveled through detector  )
        //Length traveled through flat detector (Si(Li), CBE_top,bot, UMB) = L_tof/tkr / cos(theta)
        //Length traveled through vertical detector CBE_sides, COR = L_tof/tkr / sqrt(1 - cos^2(theta))

    for(uint t = 0; t < Event->GetNTracks(); t++){
        cout << "Track is " << t << endl;
        if(Event->GetTrack(t)->IsPrimary()) cout << "I am the Primary track! " << endl;
        //cout << "Primary Track Momentum Direction is " << Event->GetPrimaryMomentumDirection()<< endl;
        //cout << "Primary: Cos(Theta) is " << Event->GetPrimaryMomentumDirection().CosTheta() << endl;
        //cout << "Track: " << t << " Momentum Direction[0] is " << Event->GetTrack(t)->GetMomentumDirection()[0] << endl;
        cout << "Track: " << t << " Momentum Direction[0][0] " << Event->GetTrack(t)->GetMomentumDirection()[0][0] << endl;
        cout << "Track: " << t << " Momentum Direction[0][1] " << Event->GetTrack(t)->GetMomentumDirection()[0][1] << endl;
        cout << "Track: " << t << " Momentum Direction[0][2] " << Event->GetTrack(t)->GetMomentumDirection()[0][2] << endl;

        float costheta = fabs(Event->GetTrack(t)->GetMomentumDirection()[0][2]);

        for(uint isig=0; isig<Event->GetTrack(t)->GetEnergyDeposition().size(); isig++){
            unsigned int VolumeId  = Event->GetTrack(t)->GetVolumeId(isig);
            //cout << "Edep " << isig << " is " << Event->GetTrack(t)->GetEnergyDeposition(isig) << " at " << VolumeId<< endl;
            if(GGeometryObject::IsTofVolume(VolumeId) && Event->GetTrack(t)->GetEnergyDeposition(isig) > TofCutLow){
                if(volspec(VolumeId,2,1) == 0 || volspec(VolumeId,2,1) == 1){
                    //cout << "Volid is " << VolumeId << " it's a flat paddle! " << endl;
                    //cout << "Step length is " << Ltof/costheta << endl;
                    if(EDEP) cout << "Edep " << isig  << " is " << Event->GetTrack(t)->GetEnergyDeposition(isig) << " at " << VolumeId << endl; // " at TOF " << volspec(VolumeId,0,3) << endl;
                    if(DEDX) cout << "dE/dx " << isig  << " is " << Event->GetTrack(t)->GetEnergyDeposition(isig)/(rtof*Ltof/costheta) << " at TOF " << volspec(VolumeId,0,3) << endl;
                }else{
                    //cout << "Volid is " << VolumeId << " it's a vertical paddle! " << endl;
                    //cout << "Step length is " << Ltof/sqrt(1 - pow(costheta,2)) << endl;
                    if(EDEP)cout << "Edep " << isig << " is " << Event->GetTrack(t)->GetEnergyDeposition(isig) << " at " << VolumeId << endl; //" at TOF " << volspec(VolumeId,0,3) << endl;
                    if(DEDX)cout << "dE/dx " << isig << " is " << Event->GetTrack(t)->GetEnergyDeposition(isig)/(rtof*Ltof/sqrt(1 - pow(costheta,2)) )  << " at TOF " << volspec(VolumeId,0,3) << endl;
                }
            }

            if(GGeometryObject::IsTrackerVolume(VolumeId) && Event->GetTrack(t)->GetEnergyDeposition(isig) > TrackerCut){
                //cout << "Volid is " << VolumeId << " it's the tracker " << endl;
                int layer = GGeometryObject::GetTrackerLayer(VolumeId);
                //cout << "Step length is " << Ltkr/costheta << endl;
                if(EDEP)cout << "Edep " << isig << " is " << Event->GetTrack(t)->GetEnergyDeposition(isig) << " at " << VolumeId << endl; /*<< " at TKR L" << layer << endl;*/
                if(DEDX)cout << "dE/dx " << isig << " is " << Event->GetTrack(t)->GetEnergyDeposition(isig)/(rtkr*Ltkr/costheta) << " at TKR L" << layer << endl;
            }

        }
    }

    } //Closed bracket for dEdx

    //cout << "Event ID? " << Event->GetEventId() << endl;
    //cout << "Event Number? " << Event->GetEventNumber() << endl; //Gviz2D is for sure pulling Event number!
    //cout << "Reconstruction used: " << Event->GetActiveReconstruction() << endl;

    /*
    //Searching for slowing down event in flight data.
    if(fabs(Event->GetPrimaryBeta()) > 0.9 && fabs(Event->GetPrimaryBeta()) < 1.0 && Event->GetNTracks() == 1 && Event->GetPrimaryBeta()*Event->GetPrimaryMomentumDirection()[2] < 0 && -fabs(Event->GetPrimaryMomentumDirection().CosTheta()) > -1 && -fabs(Event->GetPrimaryMomentumDirection().CosTheta()) < -0.54){
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
        }*/


        /*
        //Stopping particles!
        if(UMBflag > 0 && CBEtopflag > 0 && CBEbotflag == 0 && tkrflag > 2 && tkrflag < 5 && CORflag == 0 && CBEsideflag == 0 &&  Event->GetTriggerVolumeId().size() < 3){
                cout << "Event " << i << " is a cool stopping(?) friend " << endl;
               	for(unsigned int k = 0; k < Event->GetTriggerVolumeId().size(); k++){
                    unsigned int LGVolumeId = Event->GetTriggerVolumeId().at(k);
                    cout << " LG Trigger VID " << k << ": " << LGVolumeId << endl;
                }
        }*/
    //}

    //End stopping friend search

    //cout << endl << "All hits? " << endl;

    //All hits, remove on track requirement
    /*
    for(uint isig=0; isig<Event->GetVolumeId().size(); isig++){
        cout << "Hit " << isig << " VolumrId " << Event->GetVolumeId().at(isig) << " Edep " << Event->GetHitSeries().at(isig).GetTotalEnergyDeposition() << endl;
        cout << "Position " << Event->GetHitSeries().at(isig).GetPosition().X() << endl;
    }
    */

    //For MC truth data on those rad as hay events

    if(MC){

	cout << endl << "MC info " << endl;
	cout << "Number of tracks " << MCEvent->GetNTracks() << endl;
	vector<unsigned int> MCVolid;
	vector<double> MCEdep;
	vector<int> MCSpec;

    for(uint t = 0; t < MCEvent->GetNTracks();t++){
        cout << endl << "Track is " << t << endl;
        cout << "Particle is " << MCEvent->GetTrack(t)->GetPdg() << endl;
        for(uint isig=0; isig<MCEvent->GetTrack(t)->GetEnergyDeposition().size(); isig++){
            unsigned int VolumeId  = MCEvent->GetTrack(t)->GetVolumeId(isig);
            cout << "Energy deposition is " << MCEvent->GetTrack(t)->GetEnergyDeposition(isig) << " at VID " << VolumeId << endl;
            //cout << MCVolid.size() << endl;

            //int vsize = MCVolid.size();
            int index = -1;

            for (int k = 0; k < MCVolid.size(); k++) {
                if (MCVolid[k] == VolumeId) { //If volumeid is in the vector already, add the energy deposition to the vector counting edeps there
                    MCEdep[k] = MCEdep[k] + MCEvent->GetTrack(t)->GetEnergyDeposition(isig);
                    index = k;
                    //cout << "Added to existing Volid " << endl;
                }
            }
            if(index == -1){ //If VolumeId is not in the vector, add it to the vector and add the energy deposition to the energy desposition vector
                //cout << "VID not found, adding" << endl;
                MCSpec.push_back(MCEvent->GetTrack(t)->GetPdg());
                MCVolid.push_back(VolumeId);
                MCEdep.push_back(MCEvent->GetTrack(t)->GetEnergyDeposition(isig));
                //cout << "New Volid Hit! " << endl;
            }

            //cout << "MC: Hit " << isig << " Edep " << MCEvent->GetTrack(t)->GetEnergyDeposition(isig) << " at " << VolumeId << endl;
            //cout << "Method 2? " << Event->GetTrack(0)->GetEnergyDeposition().at(isig) << endl;
        }
    } //End loop over tracks


    for (int k = 0; k < MCVolid.size(); k++) { //Check to see if there's a significant hit from a primary particle (need to choose)
        if(MCEdep[k] > 0.4) cout << "High Edep! Volid " << MCVolid[k] << " Edep " << MCEdep[k] << " Main particle: " <<  MCSpec[k] << endl;
    }

    } //End MC if statement

    /*cout << endl;
    for (int k = 0; k < MCVolid.size(); k++) {
        cout << "VID " << MCVolid[k] << " total edep " << MCEdep[k] << endl;
    } //This seems to be a fine way of determining total edeps in the instrument from MC
    cout << endl;*/


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
