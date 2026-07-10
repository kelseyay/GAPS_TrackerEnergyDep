//Hello, I am a sad piece of code that exists because I needed to do something in one day and couldn't figure out how to copy everything
//Oh noooooo
//Use this on the reconstructed files. Eventually figure out how to merge me with the other code and get rid of me.
//How to use
// ./RecSeconSearch -i /home/kelsey/simulations/simdat/antip/v.3.0.0/anti_proton_gaps_triggerlevel2_FTFP_BERT_1754120716

using namespace std;

#include "KYtools.C"

#include "GDataEvent.hh"
#include "GGeometryTools.hh"
#include "GDataPoint.hh"
#include "GDataTrack.hh"
#include "GDataVertex.hh"

//FIXME: does this work on mac?
#include <sys/stat.h>

//#include "CRawTrk.hh"

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
parser->AddCommandLineOption<bool>("save", "save the special events to a root file?",0,"s");
parser->ParseCommandLine(argc, argv);
parser->Parse();

bool SAVE = parser->GetOption<bool>("save");

string reco_path = parser->GetOption<string>("in_path");

cout << reco_path << endl;

char FilenameRoot[400];
sprintf(FilenameRoot,"%s*.root",reco_path.c_str());
cout << FilenameRoot << endl;

CEventRec* Event = new CEventRec(); //New reconstructed event
TChain * TreeRec = new TChain("TreeRec"); //New TreeRec Tchain object (this is new to me)
TreeRec->SetBranchAddress("Rec", &Event); //Set the branch address using Event (defined above)
TreeRec->Add(FilenameRoot);

//Prepare FPSI Reconstruction variable for locating vertex:
Crane::Reconstruction::TrackFit::GDataEvent * reco_data_event_ = new Crane::Reconstruction::TrackFit::GDataEvent();
TChain * TreeGReco = new TChain("TreeGReco");
TreeGReco->SetBranchAddress("FindPrimaryStarIterative", &reco_data_event_); //Set the branch address using Event (defined above)
TreeGReco->Add(FilenameRoot);

//Tracker center and tolerances
double tracker_ZBaricenter = 734; //mm
double zTolerance = 500; //mm
double yTolerance = 600; //mm

//Reconstructed root file
TFile *f_template = new TFile("/home/kelsey/simulations/simdat/flight/251221/FPSI/starlink251221_0003_FindPrimaryStarIterative_rec.root", "READ");


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


//My horrible rat child: I will need to figure out how to copy the GGeometry folder in the thing!!
//Scrub the tree clean
TFile *f = new TFile("ky_root_Rec.root", "UPDATE");

// Delete the tree from memory/disk (the ;* ensures all cycles are removed)
f->Delete("TreeRec;*");
f->Delete("TreeGReco;*");

// Write changes and close file
f->Write();
f->Close();

//Next need to add another event?

TFile f2("ky_root_Rec.root", "update");
//TFile f2("ky_root_Rec.root", "update");

TTree *Copy_GRecoTree = new TTree("TreeGReco", "GReco Tree");
TTree *Copy_RecTree = new TTree("TreeRec", "Rec Tree");
Copy_GRecoTree = TreeGReco->CloneTree(0);
Copy_RecTree = TreeRec->CloneTree(0);

TreeRec->GetEntry(0);
TreeGReco->GetEntry(0);
Copy_GRecoTree->Fill();
Copy_GRecoTree->Write();
Copy_RecTree->Fill();
Copy_RecTree->Write();

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
//for(unsigned int i = 100; i < 200; i+=MainLoopScaleFactor){
    TreeRec->GetEntry(i);
    TreeGReco->GetEntry(i);
    //TreeMC->GetEntry(i);

    //cout << endl << "Event is " << i << endl;
    //cout << "Event ID? " << Event->GetEventId() << endl;
    //cout << "Event Number? " << Event->GetEventNumber() << endl; //Gviz2D is for sure pulling Event number!
    //cout << "Reconstruction used: " << Event->GetActiveReconstruction() << endl;

    /*
	for(unsigned int k = 0; k < Event->GetTriggerVolumeId().size(); k++){
		cout << "Trigger VID ? " << Event->GetTriggerVolumeId().at(k) << endl;
	}*/

	//if( TreeMC->GetEntries() % (i+1) == 0){cout << "Time at Event " << i << " = " << Event->GetEventTime() << endl;}

 	Crane::Reconstruction::TrackFit::GDataVertex dvertex;
	dvertex   = reco_data_event_->GetVertex();
	auto vertex = dvertex.GetPosition();

	//CTrackMc* pt = Event->GetPrimaryTrack();
	CTrackRec* pt = Event->GetPrimaryTrack();
	uint pt_index = 0;

	for( ; pt_index < Event->GetNTracks(); pt_index++) if( Event->GetTrack(pt_index)->IsPrimary() ) break;

	bool vertexIsOk_Reco = false;
	if(fabs(Event->GetPrimaryBeta()*Event->GetPrimaryMomentumDirection()[2] < 0 && Event->GetPrimaryBeta()) > 0 && fabs(Event->GetPrimaryBeta()) < 1.8 && -fabs(Event->GetPrimaryMomentumDirection().CosTheta()) > -1 && -fabs(Event->GetPrimaryMomentumDirection().CosTheta()) < -0.54){

	    //cout << "Event " << i << endl;
	    //Check if vertex is in the tracker
	    if( dvertex.GetNdof() > 0. ){
        	if( (vertex.Z() < tracker_ZBaricenter + zTolerance  && vertex.Z() > tracker_ZBaricenter - zTolerance) &&
                    (vertex.X() < (yTolerance)  && vertex.X() >-(yTolerance)) &&
                        (vertex.Y() < (yTolerance)  && vertex.Y() >-(yTolerance)) ) vertexIsOk_Reco = true;
        }

		int OffHitCtr = 0;

		//Remove this for now:
		int PtrackTKR = 0; //Want events where the primary has at least one L1+ TKR hit.
        for(uint isig=0; isig<Event->GetPrimaryTrack()->GetEnergyDeposition().size(); isig++){
            unsigned int VolumeId = Event->GetPrimaryTrack()->GetVolumeId().at(isig);
            if(GGeometryObject::IsTrackerVolume(VolumeId)){
                if(GGeometryObject::GetTrackerLayer(VolumeId) > 0) PtrackTKR++;   //Flag for lower layer hit
            }
        }

		int needthree = 0;
		for(uint t = 0; t < Event->GetNTracks(); t++){
		    if(Event->GetTrack(t)->GetEnergyDeposition().size() > 2) needthree++;
		}

		//Check for hits off track? (Does this only work with tracker hits??)
		for(unsigned int k = 0; k < Event->GetTotalEnergyDeposition().size(); k++){
			unsigned int VolumeId = Event->GetVolumeId().at(k);
			if(GGeometryObject::IsTrackerVolume(VolumeId)){
				if(Event->GetHitTrackIndex().at(k) == -1) OffHitCtr++; //cout << "Off track TKR hit! " << endl;
			}
			if(GGeometryObject::IsTofVolume(VolumeId)){
				if(Event->GetHitTrackIndex().at(k) == -1) OffHitCtr++; //cout << "Off track TOF hit! " << endl;
			}
		}

		//This tell me was looking for really nice reconstructed events with many secondaries.
		if(vertexIsOk_Reco && Event->GetNTracks() > 6 && Event->GetNTracks() < 9 /*&& PtrackTKR > 0*/ && OffHitCtr < 4 && needthree == Event->GetNTracks() && (pt->GetChi2()/pt->GetNdof()) < 3.2 ){
			cout << "Event " << i << " vertex in the tracker! Reasonable Secondary Number! Not so many Off track hits!" << endl;
			if(SAVE){
			    Copy_GRecoTree->Fill();
                Copy_RecTree->Fill();
			}
		}



	} //End event level cuts

}

Copy_GRecoTree->Write();
Copy_RecTree->Write();
f2.Close();
cout << endl << "I am done" << endl;
return 1;

}
