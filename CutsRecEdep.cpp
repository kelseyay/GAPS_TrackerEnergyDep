//How to use: ./RecCutEff -i /home/kelsey/simulations/simdat/mu/v.2.1.2/mu-_gaps_triggerlevel1_FTFP_BERT_1744342800_rec -o test

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
parser->AddCommandLineOption<int>("MainloopScale", "Main loop scale factor",1,"m");
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

char FilenameRoot[400];
sprintf(FilenameRoot,"%s*.root",reco_path.c_str());
cout << FilenameRoot << endl;

//Prepare textile for saving values
std::ofstream myfile;
myfile.open(out_path + "MCCharge.txt");
myfile << TString::Format( "Filename : %s", reco_path.c_str() )  << endl;
myfile << TString::Format( "Beta High : %f", betahigh) << endl;
myfile << TString::Format( "Beta Low : %f", betacut) << endl;
myfile.close();

//Prepare reconstructed event
CEventRec* Event = new CEventRec(); //New reconstructed event
TChain * TreeRec = new TChain("TreeRec"); //New TreeRec Tchain object (this is new to me)
TreeRec->SetBranchAddress("Rec", &Event); //Set the branch address using Event (defined above)
TreeRec->Add(FilenameRoot);

//int MainLoopScaleFactor = 1; //Set this number to scale the step size. Larger means runs faster and fewer events
double TofCutLow = 0; //No low Tof cut right now
double TrackerCut = 0.3; //Threshold for an energy deposition to be considered a hit

double xlow = 0.3; //Low range for histogram MeV
double xhigh = 5; //High range for histogram MeV

double coshigh = 0.54; //0.995; //0.92 //0.54 is the highest angle that can hit UMB, CBEtop, CBEbot
double coslow = 1; //0.62 //0.8
//double fitlow = 0.5;
//double fithigh = 3.5;
double acutlow = 1.5;
const Int_t NBins = 50;

//Full tracker histogram range
//double mpvmin = 0.66;
//double mpvmax = 0.75;

//From fitted peaks from flight! Check over!
double TofAngleCorrectedMip = 0.82;
double TrackerAngleCorrectedMip = 0.57;

int passB = 0; //Counter for reco pass beta
int must = 0; //Counter for reco single tracks
int trktrg = 0; //Counter for reco single tracks that pass the track trigger
int nosides = 0; //Counter for reco single tracks that pass the track trigger that have no side hits
int yesumb_ctop_cbot = 0; //Counter for reco single tracks that pass the track trigger that have no side hits YES umb, cbe_top, cbe_bot hits
int yesumb_ctop_nocbot = 0; //Counter for reco single tracks that pass the track trigger that have no side hits YES umb, cbe_top, NO cbe_bot hits

//Prepare cuts:
map<int, unsigned int> TofIndexVolumeIdMap;

//For Plotting purposes
ca::GPlottingTools Plotting;
char text[400]; //This variable is used later to name the plots

//TH1D * hedep;
//hedep = new TH1D ("h0", ("Edep l Beta " + to_string(betacut) + " - " + to_string(betahigh) ).c_str(), NBins, xlow,xhigh);

//How many entries:
cout << "Total Number of events / Mainscale Factor = " << TreeRec->GetEntries()/MainLoopScaleFactor << endl;
cout << "Beta high " << betahigh << endl;
cout << "Beta cut " << betacut << endl;
//cout << "MPV range on 2D Full Tracker will be " << mpvmin << " - " << mpvmax << endl;

//Now we can go over the loop
TreeRec->GetEntry(0);

cout << "Total Number of events / Mainscale Factor = " << TreeRec->GetEntries()/MainLoopScaleFactor << endl;

//Using i to loop over every event in the tree
//for(unsigned int i = 0; i < 40; i+=MainLoopScaleFactor){
//for(unsigned int i = 60500; i < 90600; i+=MainLoopScaleFactor){
for(unsigned int i = 0; i < TreeRec->GetEntries()/MainLoopScaleFactor; i++){ //This is not the "correct" way to do this, but it's probably fine. Should be skipping M each time, but that seems to be really slow!!
    TreeRec->GetEntry(i);

	//Cuts are implemented in this chunk:

	if( ((int)i % (int)ceil(TreeRec->GetEntries()/(MainLoopScaleFactor*10))) == 0){
	    cout << "Event number " << i << endl;
	}

	//cout << endl << "Event number " << i << endl;

	CTrackRec* pt = Event->GetPrimaryTrack();
	uint pt_index = 0;
    for( ; pt_index < Event->GetNTracks(); pt_index++) if( Event->GetTrack(pt_index)->IsPrimary() ) break;

	//Note downwards beta enforced by Event->GetPrimaryBeta() (should be positive) multiplied by Event->GetPrimaryMomentumDirection()[2] (z trajectory of particle)
	if(pt != nullptr && fabs(Event->GetPrimaryBeta()) >  betacut && fabs(Event->GetPrimaryBeta()) <  betahigh ){
	//if((pt != nullptr && pt->GetChi2()/pt->GetNdof()) < 3.2 && -fabs(MCEvent->GetPrimaryMomentumDirection().CosTheta()) > -coslow && -fabs(MCEvent->GetPrimaryMomentumDirection().CosTheta()) < -coshigh && MCEvent->GetPrimaryBeta()*MCEvent->GetPrimaryMomentumDirection()[2] < 0 && fabs(MCEvent->GetPrimaryBeta()) >  betacut && fabs(MCEvent->GetPrimaryBeta()) <  betahigh ){
	//if(pt != nullptr && -fabs(Event->GetPrimaryMomentumDirection().CosTheta()) > -coslow && -fabs(Event->GetPrimaryMomentumDirection().CosTheta()) < -coshigh && Event->GetPrimaryBeta()*MCEvent->GetPrimaryMomentumDirection()[2] < 0 && fabs(Event->GetPrimaryBeta()) >  betacut && fabs(Event->GetPrimaryBeta()) <  betahigh ){
	//if(pt != nullptr && -fabs(Event->GetPrimaryMomentumDirectionGenerated().CosTheta()) > -coslow && -fabs(Event->GetPrimaryMomentumDirectionGenerated().CosTheta()) < -coshigh && Event->GetPrimaryBetaGenerated()*Event->GetPrimaryMomentumDirectionGenerated()[2] < 0 && fabs(Event->GetPrimaryBetaGenerated()) >  betacut && fabs(Event->GetPrimaryBetaGenerated()) <  betahigh ){

	passB++;

        int UMBflag = 0;
        int CORflag = 0;
        int CBEtopflag = 0;
        int CBEbotflag = 0;
        int CBEsideflag = 0;

    if(Event->GetNTracks() == 1){
        must++;

        //Reconstructed information, fortunately only one track.
        for(uint isig=0; isig<Event->GetTrack(0)->GetEnergyDeposition().size(); isig++){
            unsigned int VolumeId  = Event->GetTrack(0)->GetVolumeId(isig);
            if(Event->GetTrack(0)->GetEnergyDeposition(isig) > 0.4){
                if(volspec(VolumeId,0,3) == 100)UMBflag++;
                if(volspec(VolumeId,0,3) == 110)CBEtopflag++;
                if(volspec(VolumeId,0,3) == 111)CBEbotflag++;
                if(volspec(VolumeId,0,3) == 112 || volspec(VolumeId,0,3) == 113 || volspec(VolumeId,0,3) == 114 || volspec(VolumeId,0,3) == 115 || volspec(VolumeId,0,3) == 116)CBEsideflag++;
                if(volspec(VolumeId,0,3) == 102 || volspec(VolumeId,0,3) == 103 || volspec(VolumeId,0,3) == 104 || volspec(VolumeId,0,3) == 105 || volspec(VolumeId,0,3) == 106)CORflag++;
            }
            //cout << "Rec: Hit " << isig << " Edep " << Event->GetTrack(0)->GetEnergyDeposition(isig) << " at " << VolumeId << endl;
        }

        //cout << "UMBflag " << UMBflag << endl;
        //cout << "CORflag " << CORflag << endl;
        //cout << "CBE_topflag " << CBEtopflag << endl;
        //cout << "CBE_botflag " << CBEbotflag << endl;
        //cout << "CBE_sideflag " << CBEsideflag << endl;
        //cout << endl;

        if( (UMBflag + CORflag > 0) && (CBEtopflag + CBEbotflag + CBEsideflag > 0)){ //Track trigger
            trktrg++; //cout << "Track trigger passed! " << endl; //Counter for MC true single tracks that pass the track trigger
            if(CBEsideflag == 0 && CORflag == 0 ){ //No side TOF hits
                nosides++; //cout << "Event: " << i << " NO SIDES HIT!" << endl;
                if( (UMBflag>0) && (CBEtopflag>0) && (CBEbotflag>0) ) {
                    yesumb_ctop_cbot++;
                    //cout << "Event: " << i << " UMB TOP BOT YES!" << endl;
                }
                if( (UMBflag>0) && (CBEtopflag>0) && (CBEbotflag==0) ){
                    yesumb_ctop_nocbot++;
                    //cout << "Event: " << i << " UMB TOP NO BOT!" << endl;
                }
            } //No sides

        } //Track trigger satisfied if one hit in outer TOF (UMB or COR) and one hit in innner TOF (any CBE)

    } //Closed bracket single track reco

			//-----------EVENT LEVEL CUTS END

		} //Closed bracket for event level cut

}  //Closed bracket for iteration through tree events, move on to the next event i


myfile.open(out_path + "RecCuts.txt",std::ios::app);
myfile.close();

cout << "Total Events in Beta Range " << passB << endl;
cout << fixed << setprecision(1) << 100*(float)passB/(float)TreeRec->GetEntries() <<"%" <<  endl;
cout << "Reco Single Track Events " << must << endl;
cout << fixed << setprecision(1) << 100*(float)must/(float)passB << "%" << endl;
cout << "Reco Track Trigger Satisfied " << trktrg << endl;
cout << fixed << setprecision(1) << 100*(float)trktrg/(float)must << "%" << endl;
cout << "Reco Track Trigger Satisfied, No sides " << nosides << endl;
cout << fixed << setprecision(1) << 100*(float)nosides/(float)trktrg << "%" << endl;
cout << "Reco Track Trigger Satisfied, No sides, UMB, CT, CB " << yesumb_ctop_cbot << endl;
cout << fixed << setprecision(1) << 100*(float)yesumb_ctop_cbot/(float)nosides <<"%" <<  endl;
cout << "Reco Track Trigger Satisfied, No sides or CB, UMB, CT " << yesumb_ctop_nocbot << endl;
cout << fixed << setprecision(1) << 100*(float)yesumb_ctop_nocbot/(float)nosides <<"%" <<  endl;
cout << endl << "I am done" << endl;

return 1;

}
