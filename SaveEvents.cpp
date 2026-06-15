//If you want to make a dead simple histogram of something in the Rec tree, here you go!!

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
parser->AddCommandLineOption<string>("out_file", "name of output root file", "", "o");
parser->AddCommandLineOption<double>("beta_low", "low Beta Cut",0.8,"l");
parser->AddCommandLineOption<double>("beta_high", "upper Beta Cut",1,"u");
parser->AddCommandLineOption<bool>("gen", "generated beta plots",0,"g");
parser->ParseCommandLine(argc, argv);
parser->Parse();

string reco_path = parser->GetOption<string>("in_path");
string out_path = parser->GetOption<string>("out_file");
cout << reco_path << endl;
//sprintf(FilenameRoot,"%s/%s*.root",argv[1], argv[2]);
//cout << argv[1] << endl;

bool GEN = parser->GetOption<bool>("gen");

double betacut = parser->GetOption<double>("beta_low");
if(betacut <= 0 || betacut >=1){ betacut = 0.8; cout << "Error with low beta choice. Setting Beta low to 0.8" << endl; }
cout << "beta cut = " << betacut << endl;

double betahigh = parser->GetOption<double>("beta_high");
if(betahigh <= 0 || betahigh >=2 || betahigh < betacut){ betahigh = 1; cout << "Error with high/low beta choice! Setting Beta upper to 1" << endl; }
cout << "beta high = " << betahigh << endl;

double coshigh = 0.54; //0.995; //0.92 //0.54 is the highest angle that can hit UMB, CBEtop, CBEbot
double coslow = 1; //0.62 //0.8

char FilenameRoot[400];
sprintf(FilenameRoot,"%s*.root",reco_path.c_str());
//sprintf(FilenameRoot,reco_path.c_str()); //210 simu data on my computer!
cout << FilenameRoot << endl;

CEventRec* Event = new CEventRec(); //New reconstructed event
TChain * TreeRec = new TChain("TreeRec"); //New TreeRec Tchain object (this is new to me)
TreeRec->SetBranchAddress("Rec", &Event); //Set the branch address using Event (defined above)
TreeRec->Add(FilenameRoot);

Crane::Reconstruction::TrackFit::GDataEvent * reco_data_event_ = new Crane::Reconstruction::TrackFit::GDataEvent();
TChain * TreeGReco = new TChain("TreeGReco");
TreeGReco->SetBranchAddress("FindPrimaryStarIterative", &reco_data_event_); //Set the branch address using Event (defined above)
TreeGReco->Add(FilenameRoot);

//Scrub the tree clean
TFile *f = new TFile("ky_root.root", "UPDATE");

// Delete the tree from memory/disk (the ;* ensures all cycles are removed)
f->Delete("TreeMc;*");
f->Delete("TreeRec;*");
f->Delete("TreeGReco;*");
f->Delete("SimulationParameterTree;*");

// Write changes and close file
f->Write();
f->Close();

//Next need to add another event?

TFile f2("ky_root.root", "update");
TTree *Copy_GRecoTree = new TTree("TreeGReco", "GReco Tree");
TTree *Copy_RecTree = new TTree("TreeRec", "Rec Tree");
Copy_GRecoTree = TreeGReco->CloneTree(0);
Copy_RecTree = TreeRec->CloneTree(0);
/*
TreeRec->GetEntry(740);
TreeGReco->GetEntry(740);
Copy_GRecoTree->Fill();
Copy_GRecoTree->Write();
Copy_RecTree->Fill();
Copy_RecTree->Write();

TreeRec->GetEntry(1000);
TreeGReco->GetEntry(1000);
Copy_GRecoTree->Fill();
Copy_GRecoTree->Write();
Copy_RecTree->Fill();
Copy_RecTree->Write();
*/


//f2.Close();
//SEEMS GOOD!!!!

int MainLoopScaleFactor = 1; //Set this number to scale the step size. Larger means runs faster and fewer events
double TrackerCut = 0.3; //Threshold for an energy deposition to be considered a hit

//Prepare cuts:
map<int, unsigned int> TofIndexVolumeIdMap;

//For Plotting purposes
ca::GPlottingTools Plotting;
char text[400]; //This variable is used later to name the plots

//How many entries:
cout << "Total Number of events / Mainscale Factor = " << TreeRec->GetEntries()/MainLoopScaleFactor << endl;

//Now we can go over the loop
TreeRec->GetEntry(0);



//Using i to loop over every event in the tree
for(unsigned int i = 0; i < TreeRec->GetEntries(); i+=MainLoopScaleFactor){
    TreeRec->GetEntry(i);
    TreeGReco->GetEntry(i);

    if(i == 740 || i == 1100){ //I thiiiink we are in business :)))
        Copy_GRecoTree->Fill();
        Copy_GRecoTree->Write();
        Copy_RecTree->Fill();
        Copy_RecTree->Write();
    }

    if( ((int)i % (int)ceil(TreeRec->GetEntries()/10)) == 0){
		    cout << "Event number " << i << endl;
	}

}  //Closed bracket for iteration through tree events, move on to the next event i

f2.Close();
cout << endl << "I am done" << endl;

return 1;

}
