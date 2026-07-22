//Pretty simple Beta histograms (work for MC and ground/flight data both.)
//This should be combined with BetaNoCuts, and no cuts should just be a parameter in the executable, like -c 0 or -c 1 for cuts vs no cuts.
//Yeah let me go ahead and fix that next iteration.
//This will eventually join BetaHisto

//How to use:
// /home/kelsey/GAPS_TrackerEnergyDep/build/PrimaryBeta -i /home/kelsey/simulations/simdat/flight/251221/v26.03/merged251221_0 -l 0.2 -u 1.8 -r 4

using namespace std;

#include "KYtools.C"

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
parser->AddCommandLineOption<int>("TRG", "Which trigger?",0,"r");
parser->ParseCommandLine(argc, argv);
parser->Parse();

string reco_path = parser->GetOption<string>("in_path");
string out_path = parser->GetOption<string>("out_file");
cout << reco_path << endl;
//sprintf(FilenameRoot,"%s/%s*.root",argv[1], argv[2]);
//cout << argv[1] << endl;

bool GEN = parser->GetOption<bool>("gen");
int TRG = parser->GetOption<int>("TRG");

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

int MainLoopScaleFactor = 1; //Set this number to scale the step size. Larger means runs faster and fewer events
double TrackerCut = 0.3; //Threshold for an energy deposition to be considered a hit

//Prepare cuts:
map<int, unsigned int> TofIndexVolumeIdMap;

//For Plotting purposes
ca::GPlottingTools Plotting;
char text[400]; //This variable is used later to name the plots

//What histogram would you like to plot here!!
TH1D * HPrimaryBeta = Plotting.DefineTH1D("HPrimaryBeta",20, betacut, betahigh, "Reconstructed Primary Beta", "entries", 10, 100000);
TH1D * HPrimaryBeta_NoTOFCuts = Plotting.DefineTH1D("HPrimaryBeta_NoTOFCuts",20, betacut, betahigh, "Reconstructed Primary Beta", "entries", 10, 100000);
//Reconstructed beta vs generated beta for MC!
TH1D * HGenPrimaryB = Plotting.DefineTH1D("HGenPrimaryB",20, betacut, 1, "Generated Primary Beta", "entries", 10, 100000);
TH2D * HRecPrimaryB_vs_GenB= new TH2D("HRecPrimaryB_vs_GenB","Primary_Rec_B vs Gen_Beta",50, betacut, 1,50, 0.1,1.5);

int bcounts = 0;
int bcounts_noTOF = 0;

//How many entries:
cout << "Total Number of events / Mainscale Factor = " << TreeRec->GetEntries()/MainLoopScaleFactor << endl;

//Now we can go over the loop
TreeRec->GetEntry(0);

//Using i to loop over every event in the tree
for(unsigned int i = 0; i < TreeRec->GetEntries(); i+=MainLoopScaleFactor){
    TreeRec->GetEntry(i);

    if( ((int)i % (int)ceil(TreeRec->GetEntries()/10)) == 0){
		    cout << "Event number " << i << endl;
	}

CTrackRec* pt = Event->GetPrimaryTrack();
	uint pt_index = 0;
    for( ; pt_index < Event->GetNTracks(); pt_index++) if( Event->GetTrack(pt_index)->IsPrimary() ) break;

	//Cuts are implemented in this chunk:
	if(pt != nullptr && (pt->GetChi2()/pt->GetNdof()) < 3.2 && ( (TRG == 0) || ((int)Event->GetTriggerSources().at(0) == TRG) ) && -fabs(Event->GetPrimaryMomentumDirection().CosTheta()) > -coslow && -fabs(Event->GetPrimaryMomentumDirection().CosTheta()) < -coshigh && Event->GetPrimaryBeta()*Event->GetPrimaryMomentumDirection()[2] < 0 && fabs(Event->GetPrimaryBeta()) >  betacut && fabs(Event->GetPrimaryBeta()) <  betahigh ){

	    int UMBflag = 0;
        int CORflag = 0;
        int CBEtopflag = 0;
        int CBEbotflag = 0;
        int CBEsideflag = 0;

		//Note downwards beta enforced by Event->GetPrimaryBeta() (should be positive) multiplied by Event->GetPrimaryMomentumDirection()[2] (z trajectory of particle)
		//if(pt != nullptr && fabs(Event->GetPrimaryBetaGenerated()) >  betacut && fabs(Event->GetPrimaryBetaGenerated()) < betahigh){
		//if(pt != nullptr && fabs(Event->GetPrimaryBeta()) >  betacut && fabs(Event->GetPrimaryBeta()) < betahigh){

		bcounts_noTOF++;

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

        if(UMBflag > 0 && CBEtopflag > 0 /*&& CBEbotflag > 0*/ ){
			HPrimaryBeta->Fill(Event->GetPrimaryBeta());
			bcounts++;
			if(GEN) HGenPrimaryB->Fill(Event->GetPrimaryBetaGenerated());
			if(GEN) HRecPrimaryB_vs_GenB->Fill(Event->GetPrimaryBetaGenerated(),Event->GetPrimaryBeta());
        }

		} //Closed bracket for event level cut

}  //Closed bracket for iteration through tree events, move on to the next event i

cout << bcounts << endl;
HPrimaryBeta->SetMaximum(bcounts);
//HPrimaryBeta_NoTOFCuts->SetMaximum(bcounts_noTOF);

//HGenB->Scale(1.0 / HGenB->Integral(), "width");

//Histogram section
//-------------------------------------
if(GEN) histplot2d("c1",HRecPrimaryB_vs_GenB,"Rec_B versus Gen_B","Generated Beta", "Reconstructed Beta","NEntries", out_path + "Cuts_BothgenBRecPrimaryB_TRG" + to_string(TRG) );
if(GEN) histplot1d("c2",HGenPrimaryB,"Gen_B","Generated Beta","NEntries", out_path + "Cuts_GenPrimaryB_TRG" + to_string(TRG) );
histplot1d("c3",HPrimaryBeta,"Reconstructed Primary B","Reconstructed Primary Beta","NEntries", out_path + "Cuts_Rec_Pri_B_TRG" + to_string(TRG) );
//histplot1d("c4",HPrimaryBeta_NoTOFCuts,"Reconstructed Primary B No TOF Cuts","Reconstructed Primary Beta","NEntries", out_path + "NoTofCuts_Rec_Pri_B_TRG" + to_string(TRG) );

cout << endl << "I am done" << endl;

return 1;

}
