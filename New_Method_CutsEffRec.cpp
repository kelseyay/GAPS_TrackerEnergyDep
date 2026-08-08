//I am branching out and trying to do a thing in a way that is smart!
//This code is starting by using an intelligent framework for counting numbers of events that pass certain cuts!

//How to use:
// ./BPRecCutEff -i /home/kelsey/simulations/simdat/ground/251204/25.10/ethernet251204_082 -o test

using namespace std;

#include "KYtools.C"

using namespace Crane::Analysis;
namespace ca = Crane::Analysis;
namespace cl = Crane::Common;

double EdepCutLow = 0.4; //No low Tof cut right now
double betacut = 0.2;
double betahigh = 1.2;

//Use TOF flags gathered for
void tof_efficiencies(int cuts[], int& k, bool& pass, int UMBflag, int CORflag, int CBEtopflag, int CBEbotflag, int CBEsideflag, int TKRflag){

    if(UMBflag > 0){
        cuts[k]++;k++; //No sides, yes UMB
        if(CBEtopflag > 0){
            cuts[k]++;k++; //No sides, yes UMB, yes CBEtop
            if((CBEbotflag>0) || (TKRflag>0) || (CBEsideflag>0)){
                cuts[k]++;k++; //No sides, yes UMB, yes CBEtop, yes CBEbot
                pass = 1;
                //cout << "Event Passes ALL CUTS!! " << endl;
            }
        }
    }
    k = 7;

}

//TOF flags are marked
//There is a slightly (better?) way to do this maybe, see MC_Beta_Timning_Res
void tof_flags(const CEventRec* Event, int& UMBflag, int& CORflag, int& CBEtopflag, int& CBEbotflag, int& CBEsideflag, int& TKRflag ){
    for(uint isig=0; isig<Event->GetTrack(0)->GetEnergyDeposition().size(); isig++){
        unsigned int VolumeId  = Event->GetTrack(0)->GetVolumeId(isig);
        if(Event->GetTrack(0)->GetEnergyDeposition(isig) > EdepCutLow){
            if(volspec(VolumeId,0,3) == 100)UMBflag++;
            if(volspec(VolumeId,0,3) == 110)CBEtopflag++;
            if(volspec(VolumeId,0,3) == 111)CBEbotflag++;
            if(volspec(VolumeId,0,3) == 112 || volspec(VolumeId,0,3) == 113 || volspec(VolumeId,0,3) == 114 || volspec(VolumeId,0,3) == 115 || volspec(VolumeId,0,3) == 116)CBEsideflag++;
            if(volspec(VolumeId,0,3) == 102 || volspec(VolumeId,0,3) == 103 || volspec(VolumeId,0,3) == 104 || volspec(VolumeId,0,3) == 105 || volspec(VolumeId,0,3) == 106)CORflag++;
            if(volspec(VolumeId,0,2) == 20){ TKRflag++; }
        }
        //cout << "Rec: Hit " << isig << " Edep " << Event->GetTrack(0)->GetEnergyDeposition(isig) << " at " << VolumeId << endl;
    }
}

//Standard event selection cuts, turn "pass" variable into pass if success
void std_cuts(const CEventRec* Event, bool& pass, int& k, int cuts[]){
    //cout << "pass right now is " << pass << endl;
   	CTrackRec* pt = Event->GetPrimaryTrack();
	uint pt_index = 0;
    for( ; pt_index < Event->GetNTracks(); pt_index++) if( Event->GetTrack(pt_index)->IsPrimary() ) break;

    if(pt != nullptr){
       cuts[k]++; k++; //PT found
       if(fabs(Event->GetPrimaryBeta()) >  betacut && fabs(Event->GetPrimaryBeta()) <  betahigh){
           cuts[k]++; k++;
           if(Event->GetNTracks() == 1){
               cuts[k]++;k++; //ST cut
               pass = true;
           }
       }
    }
    k = 4; //You do need to manually set k if you want to do it this way.

}


int main(int argc, char *argv[]){

GOptionParser* parser = GOptionParser::GetInstance();
parser->AddProgramDescription("Minimal Reproducable Example for Extracing Data from Reco Data");
parser->AddCommandLineOption<string>("in_path", "path to instrument data files", "./*", "i");
parser->AddCommandLineOption<string>("out_file", "name of output root file", "", "o");
parser->AddCommandLineOption<int>("TRG", "Which trigger?",0,"r");
parser->AddCommandLineOption<int>("MainloopScale", "Main loop scale factor",1,"m");
parser->ParseCommandLine(argc, argv);
parser->Parse();

/*
double betacut = parser->GetOption<double>("beta_low");
if(betacut < 0 || betacut >=1){ betacut = 0.8; cout << "Error with low beta choice. Setting Beta low to 0.8" << endl; }
cout << "beta cut = " << betacut << endl;

double betahigh = parser->GetOption<double>("beta_high");
if(betahigh < 0 || betahigh >=5 || betahigh < betacut){ betahigh = 1; cout << "Error with high/low beta choice! Setting Beta upper to 1" << endl; }
cout << "beta high = " << betahigh << endl;
*/

string reco_path = parser->GetOption<string>("in_path");
string out_path = parser->GetOption<string>("out_file");
cout << "reco_path " << reco_path << endl;
cout << "out path " << out_path << endl;
//cout << "out path last string " << out_path[out_path.length()-1] << endl;

if(out_path != "" && out_path[out_path.length()-1] != '/' ){ cout <<  "out path no slash!" << endl; out_path = out_path + '/'; }

string title = "RecCuts.txt";

int MainLoopScaleFactor = parser->GetOption<int>("MainloopScale");

char FilenameRoot[400];
sprintf(FilenameRoot,"%s*.root",reco_path.c_str());
cout << FilenameRoot << endl;

//Prepare reconstructed event
CEventRec* Event = new CEventRec(); //New reconstructed event
TChain * TreeRec = new TChain("TreeRec"); //New TreeRec Tchain object (this is new to me)
TreeRec->SetBranchAddress("Rec", &Event); //Set the branch address using Event (defined above)
TreeRec->Add(FilenameRoot);

//I'm going to try out this new method! Try to make it easy to swap cuts in and out.
//Prepare textile for saving values
std::ofstream myfile;
myfile.open(out_path + title);
myfile << TString::Format( "Filename : %s", reco_path.c_str() )  << endl;
myfile << TString::Format( "Beta High : %f", betahigh) << endl;
myfile << TString::Format( "Beta Low : %f", betacut) << endl;
myfile.close();

int TRG = parser->GetOption<int>("TRG");

//This was the correct way to do this.
//List out your cuts here
string cutnames[] = {
    "Total NEvents ",
    "Reco Total Events with Primary Track ",
    "Reco Total Events with PT found, in Beta Range ",
    "Reco PT yes, Beta Range, Single Track Events ",
    "Reco all above YES UMB ",
    "YES UMB, CT ",
    "YES UMB, CT, (CB, TKR, or CS) ",
};
const int NCuts = size(cutnames) ;
int cuts[NCuts] = {}; //An array that counts the number of cuts
int k = 0;
cuts[k] = TreeRec->GetEntries(); //Number of events (not a cut)
k++;

//Prepare cuts:
map<int, unsigned int> TofIndexVolumeIdMap;

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

	if( ((int)i % (int)ceil(TreeRec->GetEntries()/(MainLoopScaleFactor*10))) == 0){
	    cout << "Event number " << i << endl;
		cout << Event->GetActiveReconstruction() << endl;
	}

	bool pass = 0;
	k = 1;

	std_cuts(Event,pass,k,cuts);
	if(pass){

        int UMBflag = 0;
        int CORflag = 0;
        int CBEtopflag = 0;
        int CBEbotflag = 0;
        int CBEsideflag = 0;
        int TKRflag = 0;

        bool eff_pass= 0;

        tof_flags(Event, UMBflag, CORflag, CBEtopflag, CBEbotflag, CBEsideflag, TKRflag);
        //cout << "event is " << i << endl;
        tof_efficiencies(cuts,k,eff_pass,UMBflag, CORflag, CBEtopflag, CBEbotflag, CBEsideflag, TKRflag);
        //cout << "Did is pass the pass? " << eff_pass << endl;

			//-----------EVENT LEVEL CUTS END
	} //Closed bracket for standard event selection cuts.

}  //Closed bracket for iteration through tree events, move on to the next event i

myfile.open(out_path + title,std::ios::app); //Open the file for modifying
for(int k = 0; k < NCuts; k++){
    myfile << cutnames[k] << endl; //<< cuts[j] << endl;
    if(k > 0) myfile << fixed << setprecision(2) << 100*(float)cuts[k]/(float)cuts[k-1] <<"%" <<  endl; //NEntries doesn't need to divide by anything
}
myfile.close();

cout << endl << "I am done" << endl;

return 1;

}
