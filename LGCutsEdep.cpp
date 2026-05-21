//I feel the need to look at these studies with LG hits!
//Let's start with just reconstructed ground data.
//DOES NOT seem like MC has trigger VID information, which makes this a little tricky, but 2025 data does!
//How to use: ./LGCutEff -i /home/kelsey/simulations/simdat/mu/v.2.1.2/mu-_gaps_triggerlevel1_FTFP_BERT_1744342800_rec -o test

using namespace std;

#include "KYtools.C"

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
if(betacut < 0 || betacut >=1){ betacut = 0.8; cout << "Error with low beta choice. Setting Beta low to 0.8" << endl; }
cout << "beta cut = " << betacut << endl;

double betahigh = parser->GetOption<double>("beta_high");
if(betahigh < 0 || betahigh >=5 || betahigh < betacut){ betahigh = 1; cout << "Error with high/low beta choice! Setting Beta upper to 1" << endl; }
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

//Prepare reconstructed event
CEventRec* Event = new CEventRec(); //New reconstructed event
TChain * TreeRec = new TChain("TreeRec"); //New TreeRec Tchain object (this is new to me)
TreeRec->SetBranchAddress("Rec", &Event); //Set the branch address using Event (defined above)
TreeRec->Add(FilenameRoot);

//int MainLoopScaleFactor = 1; //Set this number to scale the step size. Larger means runs faster and fewer events
double TofCutLow = 0.4; //No low Tof cut right now
//double TrackerCut = 0.3; //Threshold for an energy deposition to be considered a hit

string title = "LGRecCuts.txt";

//I'm going to try out this new method! Try to make it easy to swap cuts in and out.
//Prepare textile for saving values
std::ofstream myfile;
myfile.open(out_path + title);
myfile << TString::Format( "Filename : %s", reco_path.c_str() )  << endl;
myfile << TString::Format( "Beta High : %f", betahigh) << endl;
myfile << TString::Format( "Beta Low : %f", betacut) << endl;
myfile.close();

const int NCuts = 6; //Start simple!
string cutnames[NCuts+1] = {
    "Total NEvents ",
    "Non-Empty Trigger VID Vector ",
    "Reco triggered LG No Sides ",
    "Reco LG No Sides, YES UMB ",
    "Reco LG No Sides, YES UMB, CT ",
    "Reco LG No Sides, YES UMB, CT, YES CB ",
    "Reco LG No Sides, YES UMB, CT, NO CB ",
};

int cuts[NCuts+1] = {};
cuts[0] = TreeRec->GetEntries(); //Number of events (not a cut)

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
//for(unsigned int i = 0; i < 4; i+=MainLoopScaleFactor){
//for(unsigned int i = 60500; i < 90600; i+=MainLoopScaleFactor){
for(unsigned int i = 0; i < TreeRec->GetEntries()/MainLoopScaleFactor; i++){ //This is not the "correct" way to do this, but it's probably fine. Should be skipping M each time, but that seems to be really slow!!
    TreeRec->GetEntry(i);

	//Cuts are implemented in this chunk:

	int UMBflag = 0;
    int CORflag = 0;
    int CBEtopflag = 0;
    int CBEbotflag = 0;
    int CBEsideflag = 0;

    //cout << "Event is " << i << endl;

	if( ((int)i % (int)ceil(TreeRec->GetEntries()/(MainLoopScaleFactor*10))) == 0){
	    cout << "Event number " << i << endl;
		cout << Event->GetActiveReconstruction() << endl;
	}

	if(Event->GetTriggerVolumeId().size() != 0) cuts[1]++ ; //Non zero TVId value, but maybe could be > 2 or something.

	for(unsigned int k = 0; k < Event->GetTriggerVolumeId().size(); k++){
	    unsigned int VolumeId = Event->GetTriggerVolumeId().at(k);
	    if(volspec(VolumeId,0,3) == 100)UMBflag++;
        if(volspec(VolumeId,0,3) == 110)CBEtopflag++;
        if(volspec(VolumeId,0,3) == 111)CBEbotflag++;
        if(volspec(VolumeId,0,3) == 112 || volspec(VolumeId,0,3) == 113 || volspec(VolumeId,0,3) == 114 || volspec(VolumeId,0,3) == 115 || volspec(VolumeId,0,3) == 116)CBEsideflag++;
        if(volspec(VolumeId,0,3) == 102 || volspec(VolumeId,0,3) == 103 || volspec(VolumeId,0,3) == 104 || volspec(VolumeId,0,3) == 105 || volspec(VolumeId,0,3) == 106)CORflag++;
	}

	if(CBEsideflag == 0 && CORflag == 0 ){ //No side TOF hits
        cuts[2]++; //No sides
        //cout << "No sides" << endl;
            if(UMBflag > 0){
                cuts[3]++; //No sides, yes UMB
                //cout << "yes UMB " << endl;
                if(CBEtopflag > 0){
                    cuts[4]++; //No sides, yes UMB, yes CBEtop
                        //cout << "yes CBE_top " << endl;
                    if((CBEbotflag>0)){
                        cuts[5]++; //No sides, yes UMB, yes CBEtop, yes CBEbot
                        cout << "Event is " << i << " yes CBot " << endl;
                    }else{
                        cuts[6]++; //No sides, yes UMB, yes CBEtop, NO CBEbot
                        cout << "Event is " << i << " no CBot" << endl;
                }
            }
        }
    } //No sides

}  //Closed bracket for iteration through tree events, move on to the next event i

myfile.open(out_path + title,std::ios::app); //Open the file for modifying
for(int k = 0; k < NCuts+1; k++){
    myfile << cutnames[k] << cuts[k] << endl;
    if(k > 0 && k!= NCuts) myfile << fixed << setprecision(2) << 100*(float)cuts[k]/(float)cuts[k-1] <<"%" <<  endl; //NEntries doesn't need to divide by anything
    if(k == NCuts) myfile << fixed << setprecision(2) << 100*(float)cuts[k]/(float)cuts[k-2] <<"%" <<  endl; //Final cut is a percentage of two above
}
myfile.close();

cout << endl << "I am done" << endl;

return 1;

}
