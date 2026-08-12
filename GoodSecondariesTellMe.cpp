//How to use
// ./SeconSearch -i /home/kelsey/simulations/simdat/antip/v.3.0.0/anti_proton_gaps_triggerlevel2_FTFP_BERT_1754120716

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
parser->AddCommandLineOption<bool>("save", "save the special events to a root file?",0,"s");
parser->AddCommandLineOption<string>("out_file", "name of output path", "", "o");
parser->AddCommandLineOption<string>("sts_root_name", "name of output root file", "sts_root_test.root", "n");
parser->ParseCommandLine(argc, argv);
parser->Parse();

bool SAVE = parser->GetOption<bool>("save");

string out_path = parser->GetOption<string>("out_file");
cout << "out path: " << out_path << endl;
if(out_path != "" && out_path[out_path.length()-1] != '/' ){ cout <<  "out path no slash! Adding! " << endl; out_path = out_path + '/'; }

string reco_path = parser->GetOption<string>("in_path");

string sts_file_name = parser->GetOption<string>("sts_root_name");
if(sts_file_name.compare(sts_file_name.length()-5,sts_file_name.length(),".root") != 0){ cout << "NO .root at the end of the root name! Adding!" << endl; sts_file_name = sts_file_name + ".root"; }

cout << reco_path << endl;
if(reco_path.compare(reco_path.length()-5,reco_path.length(),".root") == 0){ cout << ".root at the end of the reco path! Deleting!" << endl; reco_path = reco_path.substr(0,reco_path.length()-5); }

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

TFile *f;
TFile *f_source;
TDirectoryFile *GOptions_copy;

TTree *Copy_GRecoTree = new TTree("TreeGReco", "GReco Tree");
TTree *Copy_RecTree = new TTree("TreeRec", "Rec Tree");
Copy_GRecoTree = TreeGReco->CloneTree(0);
Copy_RecTree = TreeRec->CloneTree(0);

//Next need to add another event?

//Prepare cuts:
map<int, unsigned int> TofIndexVolumeIdMap;

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
		if(vertexIsOk_Reco && Event->GetNTracks() > 3 && Event->GetNTracks() < 6 && PtrackTKR > 0 && OffHitCtr < 3 /*&& needthree == Event->GetNTracks()*/ && (pt->GetChi2()/pt->GetNdof()) < 3.2 ){
			cout << "Event " << i << " vertex in the tracker! Reasonable Secondary Number! Not so many Off track hits!" << endl;
			if(SAVE){
			    Copy_GRecoTree->Fill();
                Copy_RecTree->Fill();
			}
		}

	} //End event level cuts

}

if(SAVE){
    //Make a directory to save the root file in
    string outdir = out_path + "Slim_Trim_Skim_Search";
    char SaveDir[600];
    sprintf(SaveDir, "mkdir %s", outdir.c_str());
    int success = system(SaveDir);
    if (success == 0){std::cout << "Directory " << SaveDir <<" created!" << std::endl;};


    string full_title = out_path + "Slim_Trim_Skim_Search/" + sts_file_name;

    char SaveRootFile[600];
    f = new TFile(full_title.c_str(), "RECREATE");

    TreeRec->GetEntry(0);
    cout << "source root file " << TreeRec->GetCurrentFile()->GetName() << endl;

    TObject* geo_tree = TreeRec->GetCurrentFile()->Get("GGeometry");
    geo_tree->Write("GGeometry"); //This does work!
    Copy_GRecoTree->Write();
    Copy_RecTree->Write();
    f->Close();
}

cout << endl << "I am done" << endl;
return 1;

}
