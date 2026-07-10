//To use: ./OneHist -i /home/kelsey/simulations/simdat/mu/v.2.1.0/mu-_gaps_triggerlevel1_FTFP_BERT_HP_1721258929_rec -o test
//Note: I have currently added some more cuts to make this plot without "split" hits

using namespace std;

#include "KYtools.C"

using namespace Crane::Analysis;
namespace ca = Crane::Analysis;
namespace cl = Crane::Common;
//using Crane::Calibration;
//To run, do ./OneHist -i /home/kelsey/simulations/simdat/mu/v.2.1.2/mu -o test

int main(int argc, char *argv[]){

GOptionParser* parser = GOptionParser::GetInstance();
parser->AddProgramDescription("Minimal Reproducable Example for Extracing Data from Reco Data");
parser->AddCommandLineOption<string>("in_path", "path to instrument data files", "./*", "i");
parser->AddCommandLineOption<string>("out_file", "name of output root file", "", "o");
parser->AddCommandLineOption<double>("beta_low", "low Beta Cut",0.8,"l");
parser->AddCommandLineOption<double>("beta_high", "upper Beta Cut",1,"u");
parser->ParseCommandLine(argc, argv);
parser->Parse();

string reco_path = parser->GetOption<string>("in_path");
string out_path = parser->GetOption<string>("out_file");
double betacut = parser->GetOption<double>("beta_low");
double betahigh = parser->GetOption<double>("beta_high");

if(betacut <= 0 || betacut >=1){ betacut = 0.8; cout << "Error with low beta choice. Setting Beta low to 0.8" << endl; }
cout << "beta cut = " << betacut << endl;

if(betahigh <= 0 || betahigh >=2 || betahigh < betacut){ betahigh = 1; cout << "Error with high/low beta choice! Setting Beta upper to 1" << endl; }
cout << "beta high = " << betahigh << endl;

cout << reco_path << endl;
//sprintf(FilenameRoot,"%s/%s*.root",argv[1], argv[2]);
//cout << argv[1] << endl;

char FilenameRoot[400];
sprintf(FilenameRoot,"%s*.root",reco_path.c_str());
//sprintf(FilenameRoot,reco_path.c_str()); //210 simu data on my computer!
cout << FilenameRoot << endl;

if(out_path != "" && out_path[out_path.length()-1] != '/' ){ cout <<  "out path no slash!" << endl; out_path = out_path + '/'; }

CEventRec* Event = new CEventRec(); //New reconstructed event
TChain * TreeRec = new TChain("TreeRec"); //New TreeRec Tchain object (this is new to me)
TreeRec->SetBranchAddress("Rec", &Event); //Set the branch address using Event (defined above)
TreeRec->Add(FilenameRoot);

int MainLoopScaleFactor = 1; //Set this number to scale the step size. Larger means runs faster and fewer events
double TrackerCut = 0.4; //Threshold for an energy deposition to be considered a hit

double xlow = 0.2; //Low range for histogram MeV
double xhigh = 3.5; //High range for histogram MeV

double coshigh = 0.54; //0.995; //0.92 //0.54 is the highest angle that can hit UMB, CBEtop, CBEbot
double coslow = 1; //0.62 //0.8
double fitlow = 0.5;
double fithigh = 3.5;
const Int_t NBins = 50;

//Full tracker histogram range
double mpvmin = 0.66;
double mpvmax = 0.75;

TH1F * hedep;
hedep = new TH1F ("h0", ("Edep l Beta " + to_string(betacut) + " - " + to_string(betahigh) ).c_str(), NBins, xlow,xhigh);

TH1F * hedep_nocuts;
hedep_nocuts = new TH1F ("h1", ("Edep l Beta " + to_string(betacut) + " - " + to_string(betahigh) ).c_str(), NBins, xlow,xhigh);

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
//for(unsigned int i = 0; i < 914006; i+=MainLoopScaleFactor){
for(unsigned int i = 0; i < TreeRec->GetEntries(); i+=MainLoopScaleFactor){
    TreeRec->GetEntry(i);

	//Cuts are implemented in this chunk:
	if(Event->GetNTracks() == 1){  //First select the single track event
		bool Umbflag = 0;
		bool CBEtopflag = 0;
		bool CBEbotflag = 0;
		int  Layer_Hits_Tracker[7] = {}; //Seven layers

		if( ((int)i % (int)ceil(TreeRec->GetEntries()/10)) == 0){
		    cout << "Event number " << i << endl;
		}

		CTrackRec* pt = Event->GetPrimaryTrack();
		uint pt_index = 0;
      	        for( ; pt_index < Event->GetNTracks(); pt_index++) if( Event->GetTrack(pt_index)->IsPrimary() ) break;

		//Note downwards beta enforced by Event->GetPrimaryBeta() (should be positive) multiplied by Event->GetPrimaryMomentumDirection()[2] (z trajectory of particle)
		//if(pt != nullptr && -fabs(Event->GetPrimaryMomentumDirection().CosTheta()) > -coslow && -fabs(Event->GetPrimaryMomentumDirection().CosTheta()) < -coshigh && Event->GetPrimaryBeta()*Event->GetPrimaryMomentumDirection()[2] < 0 && fabs(Event->GetPrimaryBeta()) >  betacut && fabs(Event->GetPrimaryBeta()) <  betahigh ){
		if(pt != nullptr && -fabs(Event->GetPrimaryMomentumDirectionGenerated().CosTheta()) > -coslow && -fabs(Event->GetPrimaryMomentumDirectionGenerated().CosTheta()) < -coshigh && Event->GetPrimaryBetaGenerated()*Event->GetPrimaryMomentumDirectionGenerated()[2] < 0 && fabs(Event->GetPrimaryBetaGenerated()) >  betacut && fabs(Event->GetPrimaryBetaGenerated()) <  betahigh ){
		    //cout << "Event is " << i << endl;

			//-----------EVENT LEVEL CUT APPLIED

			for(uint isig=0; isig<Event->GetTrack(0)->GetEnergyDeposition().size(); isig++){
                unsigned int VolumeId  = Event->GetTrack(0)->GetVolumeId(isig); //Check the VolumeId of the event
                if(volspec(VolumeId,0,3) == 100) { Umbflag = 1;} // cout << "UMB hit!" <<endl ;
                if(volspec(VolumeId,0,3) == 110) {CBEtopflag = 1;}// cout << "CBE top hit!" << endl;
                if(volspec(VolumeId,0,3) == 111) {CBEbotflag = 1;}// cout << "CBE bot hit!" << endl;
                if(volspec(VolumeId,0,2) == 20) {
                    int layer = GGeometryObject::GetTrackerLayer(VolumeId);
                    Layer_Hits_Tracker[layer]++;
                }
            }

			/*
			for(int l=0; l<7; l++){
			    if( Layer_Hits_Tracker[l] > 1){
					cout << endl << "Event is " << i << endl << " Double layer hit " << l << endl;
					for(uint isig=0; isig<Event->GetTrack(0)->GetEnergyDeposition().size(); isig++){
				        unsigned int VolumeId  = Event->GetTrack(0)->GetVolumeId(isig);
						cout << "Hit is " << isig << " Edep is " << Event->GetTrack(0)->GetEnergyDeposition(isig)<< " at " << VolumeId << endl;
					}
				}
			}*/

			if(Umbflag && CBEtopflag && /*CBEbotflag && */ (pt->GetChi2()/pt->GetNdof()) < 3.2 ){
				//cout << "Event number " << i << " passes the cuts!" << endl;
				for(uint isig=0; isig<Event->GetTrack(0)->GetEnergyDeposition().size(); isig++){
					unsigned int VolumeId  = Event->GetTrack(0)->GetVolumeId(isig);
					if(GGeometryObject::IsTrackerVolume(VolumeId) && Event->GetTrack(0)->GetEnergyDeposition(isig) > TrackerCut){
					    hedep_nocuts->Fill(  (Event->GetTrack(0)->GetEnergyDeposition(isig)*fabs(Event->GetPrimaryMomentumDirection().CosTheta()))   );
					    int layer = GGeometryObject::GetTrackerLayer(VolumeId);
					    int sdstrip = GGeometryObject::GetDetectorStrip(VolumeId);
						if(sdstrip != 0 && sdstrip != 7 && Layer_Hits_Tracker[layer] < 2){
						    //cout << "Strip hit is " << sdstrip << endl;
							hedep->Fill(  (Event->GetTrack(0)->GetEnergyDeposition(isig)*fabs(Event->GetPrimaryMomentumDirection().CosTheta()))   );
						}
					} //Closed bracket for Tracker volume and tracker cutoff

				} //Closed bracket for iteration over event with TOF cuts

			} //Closed bracket for if statement for cuts


			//-----------EVENT LEVEL CUTS END


		} //Closed bracket for event level cut

	} //Closed bracket for single track cut

}  //Closed bracket for iteration through tree events, move on to the next event i



//Histogram section
//--------------------------------------


histplot1f("c1", hedep, "Energy Deposition All Hits Whole Tracker","Energy Deposition x Cos(theta) MeV","NEvents", out_path + "Fulltkr");
histplot1f("c1", hedep_nocuts, "Energy Deposition All Hits Whole Tracker","Energy Deposition x Cos(theta) MeV","NEvents", out_path + "Fulltkr_nocuts");
//Histogram for NEntries at a strip level

//--------------------------------------

cout << endl << "I am done" << endl;

return 1;

}
