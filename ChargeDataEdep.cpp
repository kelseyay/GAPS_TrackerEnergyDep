//How to use:

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
parser->AddCommandLineOption<string>("out_file", "name of output root file", "", "o");
parser->AddCommandLineOption<double>("beta_low", "low Beta Cut",0.8,"l");
parser->AddCommandLineOption<double>("beta_high", "upper Beta Cut",1,"u");
parser->AddCommandLineOption<double>("tkr_factor", "tkr factor",1,"k");
parser->AddCommandLineOption<double>("tof_factor", "tof factor",1,"f");
parser->AddCommandLineOption<int>("MainloopScale", "Main loop scale factor",1,"m");
parser->AddCommandLineOption<int>("TRG", "Which trigger?",0,"r");
parser->AddCommandLineOption<bool>("TKR", "tkr use",1,"s");
parser->AddCommandLineOption<bool>("TF", "tof use",1,"t");
parser->ParseCommandLine(argc, argv);
parser->Parse();

double betacut = parser->GetOption<double>("beta_low");
if(betacut <= 0 || betacut >=1){ betacut = 0.8; cout << "Error with low beta choice. Setting Beta low to 0.8" << endl; }
cout << "beta cut = " << betacut << endl;

double betahigh = parser->GetOption<double>("beta_high");
if(betahigh <= 0 || betahigh >=2 || betahigh < betacut){ betahigh = 1; cout << "Error with high/low beta choice! Setting Beta upper to 1" << endl; }
cout << "beta high = " << betahigh << endl;

string reco_path = parser->GetOption<string>("in_path");
string out_path = parser->GetOption<string>("out_file");
cout << "reco_path " << reco_path << endl;
cout << "out path " << out_path << endl;
//cout << "out path last string " << out_path[out_path.length()-1] << endl;

if(out_path != "" && out_path[out_path.length()-1] != '/' ){ cout <<  "out path no slash!" << endl; out_path = out_path + '/'; }

int MainLoopScaleFactor = parser->GetOption<int>("MainloopScale");
double tk = parser->GetOption<double>("tkr_factor");
double tf = parser->GetOption<double>("tof_factor");

char FilenameRoot[400];
sprintf(FilenameRoot,"%s*.root",reco_path.c_str());
cout << FilenameRoot << endl;

//Prepare textile for saving values
std::ofstream myfile;
myfile.open(out_path + "Zedep_Charge.txt");
myfile << TString::Format( "Filename : %s", reco_path.c_str() )  << endl;
myfile << TString::Format( "Beta High : %f", betahigh) << endl;
myfile << TString::Format( "Beta Low : %f", betacut) << endl;
myfile.close();

//Prepare reconstructed event
CEventRec* Event = new CEventRec(); //New reconstructed event
TChain * TreeRec = new TChain("TreeRec"); //New TreeRec Tchain object (this is new to me)
TreeRec->SetBranchAddress("Rec", &Event); //Set the branch address using Event (defined above)
TreeRec->Add(FilenameRoot);

bool TKR = parser->GetOption<bool>("TKR");
bool TF = parser->GetOption<bool>("TF");
int TRG = parser->GetOption<int>("TRG");

//int MainLoopScaleFactor = 1; //Set this number to scale the step size. Larger means runs faster and fewer events
double TofCutLow = 0; //No low Tof cut right now
double TrackerCut = 0.3; //Threshold for an energy deposition to be considered a hit

double xlow = 0.3; //Low range for histogram MeV
double xhigh = 5; //High range for histogram MeV

double coshigh = 0.54; //0.995; //0.92 //0.54 is the highest angle that can hit UMB, CBEtop, CBEbot
double coslow = 1; //0.62 //0.8
double fitlow = 0.5;
double fithigh = 3.5;
double acutlow = 1.5;
const Int_t NBins = 50;

//Full tracker histogram range
double mpvmin = 0.66;
double mpvmax = 0.75;

//From fitted peaks from flight! Check over!
double TofAngleCorrectedMip = 0.77/tf; //0.82/tf;
double TrackerAngleCorrectedMip = 0.57/tk;

int passctr = 0;
int strkctr = 0;

int alphactr = 0;
int pctr = 0;
int muctr = 0;
int yesid_actr = 0;
int yesid_pctr = 0;


//Prepare cuts:
map<int, unsigned int> TofIndexVolumeIdMap;

//For Plotting purposes
ca::GPlottingTools Plotting;
char text[400]; //This variable is used later to name the plots

TH1D * HTruncatedMeanEnergyDepositionMip = Plotting.DefineTH1D("HTruncatedMeanEnergyDepositionMip",100, 0, 10, "sqrt(sqrt(truncated mean E))nergy deposition downgoing MIP [MIP]", "entries", 0.5, 1e4);
TH1D * HChargeMip = Plotting.DefineTH1D("HChargeMip",200, 0.1, 4, "particle charge for downgoing MIP", "entries", 0.5, 1e4);

//TH2D * HGenB_vs_GenZ = new TH2D("HGenB_vs_GenZ","Gen_Beta * Gen_Z vs Gen_Beta",50,betacut - 0.1, betahigh + 0.1, 50, 0.5*betacut -0.1 , 1.5*betahigh*2 + 0.1 );
TH2D * HRecB_vs_RecBTrunM = new TH2D("HRecB_vs_RecBTrunM","Rec_Beta * Tr_Mean vs Rec_Beta",50,betacut - 0.1, betahigh + 0.1, 50, 0.1 , 3 );
TH2D * HTrunM_vs_RecB= new TH2D("HTrunM_vs_RecB","Tr_Mean vs Rec_Beta",50, betacut - 0.1, betahigh + 0.1,50,  0.1 , 3);

//TH1D * hedep;
//hedep = new TH1D ("h0", ("Edep l Beta " + to_string(betacut) + " - " + to_string(betahigh) ).c_str(), NBins, xlow,xhigh);


//How many entries:
cout << "Total Number of events / Mainscale Factor = " << TreeRec->GetEntries()/MainLoopScaleFactor << endl;
cout << "Beta high " << betahigh << endl;
cout << "Beta cut " << betacut << endl;
cout << "MPV range on 2D Full Tracker will be " << mpvmin << " - " << mpvmax << endl;

//Now we can go over the loop
TreeRec->GetEntry(0);

cout << "Total Number of events / Mainscale Factor = " << TreeRec->GetEntries()/MainLoopScaleFactor << endl;

//Using i to loop over every event in the tree
//for(unsigned int i = 0; i < 1000; i+=MainLoopScaleFactor){
for(unsigned int i = 0; i < TreeRec->GetEntries(); i+=MainLoopScaleFactor){
//for(unsigned int i = 0; i < TreeRec->GetEntries()/MainLoopScaleFactor; i++){ //This is not the "correct" way to do this, but it's probably fine. Should be skipping M each time, but that seems to be really slow!!
    TreeRec->GetEntry(i);

    if( ((int)i % (int)ceil(TreeRec->GetEntries()/(MainLoopScaleFactor*10))) == 0){
	    cout << "Event number " << i << endl;
	}

	//Cuts are implemented in this chunk:
	if(Event->GetNTracks() == 1){  //First select the single track event
	    strkctr++;
		bool Umbflag = 0;
		bool CBEtopflag = 0;
		bool CBEbotflag = 0;

		CTrackRec* pt = Event->GetPrimaryTrack();
		uint pt_index = 0;
      	        for( ; pt_index < Event->GetNTracks(); pt_index++) if( Event->GetTrack(pt_index)->IsPrimary() ) break;

		//Note downwards beta enforced by Event->GetPrimaryBeta() (should be positive) multiplied by Event->GetPrimaryMomentumDirection()[2] (z trajectory of particle)
		if((pt != nullptr && pt->GetChi2()/pt->GetNdof()) < 3.2 && ( (TRG == 0) || ((int)Event->GetTriggerSources().at(0) == TRG) ) && -fabs(Event->GetPrimaryMomentumDirection().CosTheta()) > -coslow && -fabs(Event->GetPrimaryMomentumDirection().CosTheta()) < -coshigh && Event->GetPrimaryBeta()*Event->GetPrimaryMomentumDirection()[2] < 0 && fabs(Event->GetPrimaryBeta()) >  betacut && fabs(Event->GetPrimaryBeta()) <  betahigh ){
		//if(pt != nullptr && -fabs(Event->GetPrimaryMomentumDirection().CosTheta()) > -coslow && -fabs(Event->GetPrimaryMomentumDirection().CosTheta()) < -coshigh && Event->GetPrimaryBeta()*MCEvent->GetPrimaryMomentumDirection()[2] < 0 && fabs(Event->GetPrimaryBeta()) >  betacut && fabs(Event->GetPrimaryBeta()) <  betahigh ){
		//if(pt != nullptr && -fabs(Event->GetPrimaryMomentumDirectionGenerated().CosTheta()) > -coslow && -fabs(Event->GetPrimaryMomentumDirectionGenerated().CosTheta()) < -coshigh && Event->GetPrimaryBetaGenerated()*Event->GetPrimaryMomentumDirectionGenerated()[2] < 0 && fabs(Event->GetPrimaryBetaGenerated()) >  betacut && fabs(Event->GetPrimaryBetaGenerated()) <  betahigh ){

			//-----------EVENT LEVEL CUT APPLIED

			//cout << "Event: " << i << endl;
			//First iteration over event for flags and EnergyDepositionMip vector filling
			vector<double> EnergyDepositionMip;
			for(uint isig=0; isig<Event->GetTrack(0)->GetEnergyDeposition().size(); isig++){
                unsigned int VolumeId  = Event->GetTrack(0)->GetVolumeId(isig); //Check the VolumeId of the event

                if(TF && GGeometryObject::IsTofVolume(VolumeId) && Event->GetTrack(0)->GetEnergyDeposition(isig) > TofCutLow){
                    //A hit in the COR or CBE_sides needs to be multiplied by sin(theta) instead of cos(theta)
                    if(volspec(VolumeId,2,1) == 0 || volspec(VolumeId,2,1) == 1){
                        //cout << Event->GetTrack(0)->GetEnergyDeposition(isig)*fabs(Event->GetPrimaryMomentumDirection().CosTheta())/TofAngleCorrectedMip << endl;
                        EnergyDepositionMip.push_back(Event->GetTrack(0)->GetEnergyDeposition(isig)*fabs(Event->GetPrimaryMomentumDirection().CosTheta())/TofAngleCorrectedMip);
                        //cout << "FLAT PADDLE HIT" << endl;
                    } else{
                        //cout << Event->GetTrack(0)->GetEnergyDeposition(isig)*sqrt(1-pow(Event->GetPrimaryMomentumDirection().CosTheta(),2))/TofAngleCorrectedMip << endl;
                        EnergyDepositionMip.push_back(Event->GetTrack(0)->GetEnergyDeposition(isig)*sqrt(1-pow(Event->GetPrimaryMomentumDirection().CosTheta(),2))/TofAngleCorrectedMip);
                        //cout << "VERTICAL PADDLE HIT" << endl;
                    }
                }


                if(TKR && GGeometryObject::IsTrackerVolume(VolumeId) && Event->GetTrack(0)->GetEnergyDeposition(isig) > TrackerCut){
                    //cout << Event->GetTrack(0)->GetEnergyDeposition(isig)*fabs(Event->GetPrimaryMomentumDirection().CosTheta())/TrackerAngleCorrectedMip << endl;
                    EnergyDepositionMip.push_back(Event->GetTrack(0)->GetEnergyDeposition(isig)*fabs(Event->GetPrimaryMomentumDirection().CosTheta())/TrackerAngleCorrectedMip);
                    //cout << "TRACKER HIT" << endl;
                }

                if(volspec(VolumeId,0,3) == 100){ Umbflag = 1; } // cout << "UMB hit!" <<endl ;
                //if(volspec(VolumeId,0,2) == 10 && volspec(VolumeId,2,1) != 0 ){ } //cout << "COR hit! " << endl;
                if(volspec(VolumeId,0,3) == 110) { CBEtopflag = 1; }// cout << "CBE top hit!" << endl;
                if(volspec(VolumeId,0,3) == 111) { CBEbotflag = 1; }// cout << "CBE bot hit!" << endl;

			}

			//Do we want to only run this on certain tracks? Yeah probably. Can remove the TOF flags. Just run on whatever lol.
			if(Umbflag && CBEtopflag /*&& CBEbotflag &&  1*/){

				//Calculate the truncated mean:
				//Prepare calculation :o
				//Want to take lowest energy deposits to avoid Landau long tail!

				std::sort(EnergyDepositionMip.begin(), EnergyDepositionMip.begin()+EnergyDepositionMip.size());
				double TruncatedMeanEnergyDepositionMip = 0;
				double CtrTruncatedMeanEnergyDepositionMip = 0;

				if(EnergyDepositionMip.size() == 0){ TruncatedMeanEnergyDepositionMip = 0;
				}else{

				    for(unsigned int isig = 0; isig < double(EnergyDepositionMip.size())/2; isig++){
						TruncatedMeanEnergyDepositionMip += EnergyDepositionMip.at(isig);
						CtrTruncatedMeanEnergyDepositionMip++;
					}

					if(CtrTruncatedMeanEnergyDepositionMip == 0){ TruncatedMeanEnergyDepositionMip = 0;
					}else{ TruncatedMeanEnergyDepositionMip /= CtrTruncatedMeanEnergyDepositionMip;}

				}

				if(EnergyDepositionMip.size() == 1) TruncatedMeanEnergyDepositionMip = EnergyDepositionMip.at(0);

				//cout << "Trm calculated " << TruncatedMeanEnergyDepositionMip << endl;
				HTruncatedMeanEnergyDepositionMip->Fill(TruncatedMeanEnergyDepositionMip);
				HChargeMip->Fill(sqrt(TruncatedMeanEnergyDepositionMip));
				//cout << "sqrt(truncated mean E) = " << TruncatedMeanEnergyDepositionMip << endl;
				HRecB_vs_RecBTrunM->Fill(Event->GetPrimaryBeta(), (Event->GetPrimaryBeta()) * sqrt(TruncatedMeanEnergyDepositionMip) );
				HTrunM_vs_RecB->Fill(Event->GetPrimaryBeta(),sqrt(TruncatedMeanEnergyDepositionMip));

			} //Closed bracket for if statement for cuts


			//-----------EVENT LEVEL CUTS END


		} //Closed bracket for event level cut

	} //Closed bracket for single track cut

}  //Closed bracket for iteration through tree events, move on to the next event i


string name = "";
if(TF && TKR) name = "Both";
if(TF && !TKR) name = "TOF";
if(!TF && TKR) name = "TKR";

myfile.open(out_path + "Zedep_Charge.txt",std::ios::app);
myfile << "Total Events/Mainscale Factor " << TreeRec->GetEntries()/MainLoopScaleFactor << endl;
myfile.close();


//Histogram section
//--------------------------------------

histplot1d("c1", HChargeMip, (name + ": MIP Charge Distribution Beta " + roundstr_d(betacut,2) + " - " + roundstr_d(betahigh,2)  + " TK=" + roundstr_d(TrackerAngleCorrectedMip,2) + " TF=" + roundstr_d(TofAngleCorrectedMip,2) ).c_str() ,"Charge","NEvents", out_path + name + "_Zedep");
histplot2d("c3", HRecB_vs_RecBTrunM, "Rec_B * sqrt(Tr_Mean) versus Rec_B","Reconstructed Beta","sqrt(truncated mean E) * Reconstructed B","NEntries", out_path + name + "_Zedep_Rec2D");
histplot2d("c6", HTrunM_vs_RecB, "sqrt(Tr_Mean) versus Rec_B","Reconstructed Beta","sqrt(truncated mean E)","NEntries", out_path + name + "_Zedep_RecBTrunM");

cout << "TOF: " << TF << " TKR: " << TKR << endl;
cout << "Max bin Center is  " << HChargeMip->GetBinCenter(HChargeMip->GetMaximumBin()) << " factor should be " << pow(1/HChargeMip->GetBinCenter(HChargeMip->GetMaximumBin()),2) << endl;


cout << endl << "I am done" << endl;

return 1;

}
