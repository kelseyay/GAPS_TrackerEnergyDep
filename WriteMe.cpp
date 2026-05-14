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
//using Crane::Calibration;


int main(int argc, char *argv[]){

GOptionParser* parser = GOptionParser::GetInstance();
parser->AddProgramDescription("Write event information to text file for comparison");
parser->AddCommandLineOption<string>("in_path", "path to instrument data files", "./*", "i");
parser->AddCommandLineOption<string>("out_file", "name of output root file", "", "o");
parser->AddCommandLineOption<double>("beta_low", "low Beta Cut",0.8,"l");
parser->AddCommandLineOption<double>("beta_high", "upper Beta Cut",1,"u");
parser->AddCommandLineOption<string>("name", "name of written text file","Test","n");
parser->AddCommandLineOption<bool>("cuts", "apply cuts or not?",0,"s"); //s for slice idk
parser->AddCommandLineOption<int>("TRG", "Which trigger?",0,"r");
parser->ParseCommandLine(argc, argv);
parser->Parse();

double betacut = parser->GetOption<double>("beta_low");
if(betacut <= 0 || betacut >=1){ betacut = 0.8; cout << "Error with low beta choice. Setting Beta low to 0.8" << endl; }
cout << "beta cut = " << betacut << endl;

double betahigh = parser->GetOption<double>("beta_high");
if(betahigh <= 0 || betahigh >=2 || betahigh < betacut){ betahigh = 1; cout << "Error with high/low beta choice! Setting Beta upper to 1" << endl; }
cout << "beta high = " << betahigh << endl;

int TRG = parser->GetOption<int>("TRG");
bool CUTS = parser->GetOption<bool>("cuts");

string name = parser->GetOption<string>("name");
if(CUTS) name = name + "_cuts";
if(!CUTS) name = name + "_no_cuts";
string reco_path = parser->GetOption<string>("in_path");
string out_path = parser->GetOption<string>("out_file");
cout << "reco_path " << reco_path << endl;
cout << "out path " << out_path << endl;

char FilenameRoot[400];
sprintf(FilenameRoot,"%s*.root",reco_path.c_str());
cout << FilenameRoot << endl;

CEventRec* Event = new CEventRec(); //New reconstructed event
TChain * TreeRec = new TChain("TreeRec"); //New TreeRec Tchain object (this is new to me)
TreeRec->SetBranchAddress("Rec", &Event); //Set the branch address using Event (defined above)
TreeRec->Add(FilenameRoot);

std::ofstream myfile;
string txtname;
if(CUTS) txtname = name + ".txt";
if(!CUTS) txtname = name + ".txt";
myfile.open(out_path + txtname);

int MainLoopScaleFactor = 1; //Set this number to scale the step size. Larger means runs faster and fewer events
double TrackerCut = 0.4; //Threshold for an energy deposition to be considered a hit
double TofCutLow = 0;
double xlow = 0.3; //Hist range
double xhigh = 6; //Hist range
double coshigh = 0; //0.92; //0.995; //0.92 //0.54 is the highest angle that can hit UMB, CBEtop, CBEbot
double coslow = 1; //0.62 //0.8
const Int_t NBins = 50;

//Prepare cuts:
map<int, unsigned int> TofIndexVolumeIdMap;
//All of the plots are declared here
//("Title",Number of bins,xmin,xmax,"xlabel","ylabel",ymin,ymax)

//Now we can go over the loop
TreeRec->GetEntry(0);

int PT_exist = 0;
int Track_trg = 0;
int d_beta = 0;
int B_range = 0;

int TT_PT = 0;
int TT_PT_DB = 0;
int TT_PT_DB_BR = 0;
int TT_PT_DB_BR_ST = 0;

TH1F * hbeta = new TH1F ("hbeta", "Beta Reconstructed", NBins, betacut,betahigh);
TH1F * hedep = new TH1F ("h0", ("Edep l Beta " + to_string(betacut) + " - " + to_string(betahigh) ).c_str(), NBins, 0,2);


cout << "Total Number of events / Mainscale Factor = " << TreeRec->GetEntries()/MainLoopScaleFactor << endl;

//Using i to loop over every event in the tree
//for(unsigned int i = 0; i < TreeRec->GetEntries(); i+=MainLoopScaleFactor)
for(unsigned int i = 0; i < TreeRec->GetEntries(); i+=MainLoopScaleFactor){
    TreeRec->GetEntry(i);

    if( ((int)i % (int)ceil(TreeRec->GetEntries()/(MainLoopScaleFactor*10))) == 0){
	    cout << "Event number " << i << endl;
	}

    //myfile << i << "\t" << Event->GetEventNumber() << "\t\t" << Event->GetPrimaryBeta() << endl;
    CTrackRec* pt = Event->GetPrimaryTrack();
	uint pt_index = 0;
    for( ; pt_index < Event->GetNTracks(); pt_index++) if( Event->GetTrack(pt_index)->IsPrimary() ) break;

    if(pt != nullptr) PT_exist++;
    //Hey this below might all go to heck ina handcart without the no null pointer PT, so add that in if it seg faults lol.
    if(Event->GetTriggerSources().at(0) == 2) Track_trg++;
    if(fabs(Event->GetPrimaryBeta()) >  betacut && fabs(Event->GetPrimaryBeta()) <  betahigh) B_range++;
    if(Event->GetPrimaryBeta()*Event->GetPrimaryMomentumDirection()[2] < 0) d_beta++;

    if(pt != nullptr && Event->GetTriggerSources().at(0) == 2) TT_PT++;
    if(pt != nullptr && Event->GetTriggerSources().at(0) == 2 && Event->GetPrimaryBeta()*Event->GetPrimaryMomentumDirection()[2] < 0 ) TT_PT_DB++;
    if(pt != nullptr && Event->GetTriggerSources().at(0) == 2 && Event->GetPrimaryBeta()*Event->GetPrimaryMomentumDirection()[2] < 0 && fabs(Event->GetPrimaryBeta()) >  betacut && fabs(Event->GetPrimaryBeta()) <  betahigh) TT_PT_DB_BR++;
    if(pt != nullptr && Event->GetTriggerSources().at(0) == 2 && Event->GetPrimaryBeta()*Event->GetPrimaryMomentumDirection()[2] < 0 && fabs(Event->GetPrimaryBeta()) >  betacut && fabs(Event->GetPrimaryBeta()) <  betahigh && Event->GetNTracks() == 1) TT_PT_DB_BR_ST++;


    if(!CUTS || ( CUTS && pt != nullptr && ((TRG == 0) || ((int)Event->GetTriggerSources().at(0) == TRG))  && Event->GetPrimaryBeta()*Event->GetPrimaryMomentumDirection()[2] < 0 && fabs(Event->GetPrimaryBeta()) > betacut && fabs(Event->GetPrimaryBeta()) <  betahigh && Event->GetNTracks() == 1 )  ){
        hbeta->Fill(Event->GetPrimaryBeta());
        for(uint isig=0; isig<Event->GetTrack(0)->GetEnergyDeposition().size(); isig++){
			unsigned int VolumeId  = Event->GetTrack(0)->GetVolumeId(isig);
			if(GGeometryObject::IsTrackerVolume(VolumeId) && Event->GetTrack(0)->GetEnergyDeposition(isig) > TrackerCut){
            hedep->Fill(  Event->GetTrack(0)->GetEnergyDeposition(isig)   );
			}
        }
        myfile << i << "\t" << Event->GetEventNumber()  << "\t\t" << Event->GetPrimaryBeta() << endl;
    }

    /*
    cout << endl << "Event is " << i << endl;
    cout << "Event ID? " << Event->GetEventId() << endl;
    cout << "Event Number? " << Event->GetEventNumber() << endl; //Gviz2D is for sure pulling Event number!
    cout << "Reconstruction used: " << Event->GetActiveReconstruction() << endl;
	//if( TreeMC->GetEntries() % (i+1) == 0){cout << "Time at Event " << i << " = " << Event->GetEventTime() << endl;}
	*/

    //cout << endl << "All hits? " << endl;

}

myfile.close();

cout << "NEntries " << TreeRec->GetEntries()/MainLoopScaleFactor << endl;
cout << "PT exists " << PT_exist << endl;
cout << "Track trigger " << Track_trg << endl;
cout << "Downwards beta " << d_beta << endl;
cout << "Beta within range " << B_range << endl;

hedep->GetYaxis()->SetRangeUser(1.0, 30000.0);

//string ctitle, TH1F* h1, string title, string xtitle, string ytitle, string savename
histplot1f("cbeta",hbeta,"Reconstructed Beta " + name, "Beta", "NEntries", out_path +name+ "beta");
histplot1f("c1", hedep, "Energy Deposition All Hits Whole Tracker NO ANGLE CORRECTION","Energy Deposition MeV","NEvents", out_path +name+ "Fulltkr_NoAng");

//Alright this may be the silliest way of doing things, but rewriting a new few lines:

std::stringstream buffer;

//Store all text into a buffer

std::ifstream inFile(out_path + txtname);
if (inFile) {
    buffer << inFile.rdbuf(); // Read whole file into buffer
}

//Add the new lines to the top of the text file loool
std::ofstream outfile(out_path + txtname);
if (outfile) {
    outfile << TString::Format( "Filename : %s", reco_path.c_str() )  << endl;
    outfile << TString::Format( "Beta High : %f", betahigh) << endl;
    outfile << TString::Format( "Beta Low : %f", betacut) << endl << endl;
    outfile << "NEntries " << TreeRec->GetEntries()/MainLoopScaleFactor << endl;
    outfile << "PT exists " << PT_exist << endl;
    outfile << "Track trigger " << Track_trg << endl;
    outfile << "Downwards beta " << d_beta << endl;
    outfile << "Beta within range " << B_range << endl << endl;

    outfile << "Track Trigger + PT Exists: " << TT_PT << endl;
    outfile << "Track Trigger + PT Exists + Down Beta: " << TT_PT_DB << endl;
    outfile << "Track Trigger + PT Exists + Down Beta + Beta in Range: " << TT_PT_DB_BR << endl;
    outfile << "Track Trigger + PT Exists + Down Beta + Beta in Range + Single Track: " << TT_PT_DB_BR_ST << endl;

    outfile << "Trigger " << TRG << endl;
    outfile << "CUTS? " << CUTS << endl << endl;
    outfile << "i\tEvt #\t\tBetaRec" << endl;
    outfile << buffer.str();
}

cout << endl << "I am done" << endl;
return 1;

}
