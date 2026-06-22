//Hopefully this script will provide a primary beta from energy depositions in Not That Much Time.
//This will start simple with beta from at least three tracker hits with energy deposition > 0.4

// How to use:
// ./MCEdepBeta -i /home/kelsey/simulations/simdat/mu/v.3.0.0/triggerlevel1/mu-_gaps_triggerlevel1_FTFP_BERT_1757850432_rec -l 0.4 -u 0.8


using namespace std;

#include "KYtools.C"
#include <boost/math/special_functions/lambert_w.hpp>
#include <Math/RootFinderAlgorithms.h>
#include <TF1.h>
#include <Math/RootFinder.h>
#include <Math/Functor.h>
using namespace ROOT::Math;
using namespace std;
//Headers included in KYtools are included here, so this is nice.

//FIXME: does this work on mac?
#include <sys/stat.h>

/*float Ztof = 5.574;
float Atof = 10.3;
float Ltof = 0.635;
float rtof = 1.032;
float Gtofnew = tf; //1/0.9; //(1/0.71);
//cout << "sqrt(Gtofnew) = " << sqrt(tf) << endl;
*/

float z = 1;
float Zeff = 14;
float Aeff = 28;
float rho = 2.33; //Density g/cm^2
float L = 0.25; //Mmk so for the MC, 0.95 beta (limit of solver) corresponds to 0.787 MeV

float ion = (0.000016 * pow(Zeff,0.9));
float C_1 = 0.3071/2; //MeV/ g/cm^2 #2*pi*constants not 4*pi*constants for MPV
float me = 0.511; //mass of electron * c^2

double ZOne_Tkr(double x)
{
    return C_1*pow(z,2)*(rho*L*Zeff/Aeff)*(1/pow(x,2))* ( log( (1.022 * pow(x,2)/(1 - pow(x,2)) )/ion ) + log(C_1*pow(z,2)*(rho*L*Zeff/Aeff)*(1/pow(x,2))/ion)  + 0.2 - pow(x,2) );
}

float solvermin = ZOne_Tkr(0.95); //0.787;

//Try figuring out why this worked: https://root.cern/doc/v620/classTF1.html
TF1 *Z1_tkr_Solve = new TF1("Z1_tkr_Solve", [](double *x, double *p){ return ZOne_Tkr(x[0]); }, 0.01, 0.95, 0);

//Example:
//TGraph * g = new TGraph(npointx, xvec, yvec);
//TF1 * f = new TF1("f",[&](double*x, double *p){ return p[0]*g->Eval(x[0]); }, xmin, xmax, 1);

int main(int argc, char *argv[]){

GOptionParser* parser = GOptionParser::GetInstance();
parser->AddProgramDescription("Minimal Reproducable Example for Extracing Data from Reco Data");
parser->AddCommandLineOption<string>("in_path", "path to instrument data files", "./*", "i");
parser->AddCommandLineOption<string>("out_file", "name of output root file", "", "o");
parser->AddCommandLineOption<double>("beta_low", "low Beta Cut",0.8,"l");
parser->AddCommandLineOption<double>("beta_high", "upper Beta Cut",1,"u");
parser->AddCommandLineOption<double>("tkr_factor", "tkr factor",1,"k");
parser->AddCommandLineOption<double>("tof_factor", "tof factor",1,"f");
parser->AddCommandLineOption<bool>("TKR", "tkr use",1,"s");
parser->AddCommandLineOption<bool>("TF", "tof use",1,"t");
parser->AddCommandLineOption<int>("TRG", "Which trigger?",0,"r");
parser->AddCommandLineOption<bool>("print", "print out info?",0,"n");
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

string txtname = "MCEdep_Beta.txt";

//Prepare textile for saving values
std::ofstream myfile;
myfile.open(out_path + txtname);
myfile << TString::Format( "Filename : %s", reco_path.c_str() )  << endl;
myfile << TString::Format( "Beta High : %f", betahigh) << endl;
myfile << TString::Format( "Beta Low : %f", betacut) << endl;
myfile << TString::Format( "Charge : %f", z) << endl;
myfile.close();

//Prepare reconstructed event
CEventRec* Event = new CEventRec(); //New reconstructed event
TChain * TreeRec = new TChain("TreeRec"); //New TreeRec Tchain object (this is new to me)
TreeRec->SetBranchAddress("Rec", &Event); //Set the branch address using Event (defined above)
TreeRec->Add(FilenameRoot);

//Prepare MC event
CEventMc* MCEvent = new CEventMc(); //New reconstructed event
TChain * TreeMC = new TChain("TreeMc"); //New TreeMC Tchain object (this is new to me)
TreeMC->SetBranchAddress("Mc", &MCEvent); //Set the branch address using Event (defined above)
TreeMC->Add(FilenameRoot);

//auto tkr_edep_beta_z1 = new TF1("tkr_edep_beta_z1","OneZ_TKR(x)",0.01,0.95);

//int MainLoopScaleFactor = 1; //Set this number to scale the step size. Larger means runs faster and fewer events
double TofCutLow = 0.1; //No low Tof cut right now
double TrackerCut = 0.4; //Threshold for an energy deposition to be considered a hit
int NHitsmin = 3;

double coshigh = 0.54; //0.995; //0.92 //0.54 is the highest angle that can hit UMB, CBEtop, CBEbot
double coslow = 1; //0.62 //0.8

bool TKR = parser->GetOption<bool>("TKR");
bool TF = parser->GetOption<bool>("TF");
int TRG = parser->GetOption<int>("TRG");
int print = parser->GetOption<bool>("print");

double tk = parser->GetOption<double>("tkr_factor");
double tf = parser->GetOption<double>("tof_factor");

TH2D * HRecB_vs_GenB = new TH2D("HRecB_vs_GenB","Rec_Beta * Tr_Mean vs Rec_Beta",50,betacut - 0.1, betahigh + 0.1, 50, 0.1 , 1);
TH2D * HProxB_vs_GenB = new TH2D("HProxB_vs_GenB","Prox_B vs Gen_Beta",50, betacut - 0.1, betahigh + 0.1, 50, 0.1, 1);


//Prepare cuts:
map<int, unsigned int> TofIndexVolumeIdMap;

//For Plotting purposes
ca::GPlottingTools Plotting;
char text[400]; //This variable is used later to name the plots

//How many entries:
cout << "Total Number of events / Mainscale Factor = " << TreeRec->GetEntries()/MainLoopScaleFactor << endl;
cout << "Beta high " << betahigh << endl;
cout << "Beta cut " << betacut << endl;

//Now we can go over the loop
TreeRec->GetEntry(0);

cout << "Total Number of events / Mainscale Factor = " << TreeRec->GetEntries()/MainLoopScaleFactor << endl;

//Using i to loop over every event in the tree
//for(unsigned int i = 0; i < 1000; i+=MainLoopScaleFactor){
for(unsigned int i = 0; i < TreeRec->GetEntries(); i+=MainLoopScaleFactor){
    TreeRec->GetEntry(i);
    TreeMC->GetEntry(i);

    if( ((int)i % (int)ceil(TreeRec->GetEntries()/(MainLoopScaleFactor*10))) == 0){
		    cout << "Event number " << i << endl;
	}

	//Cuts are implemented in this chunk:
	if(Event->GetNTracks() == 1){  //First select the single track event
		bool Umbflag = 0;
		bool CBEtopflag = 0;
		bool CBEbotflag = 0;
		int Layer_Hits_Tracker[7] = {}; //Seven layers

		CTrackRec* pt = Event->GetPrimaryTrack();
		uint pt_index = 0;
      	        for( ; pt_index < Event->GetNTracks(); pt_index++) if( Event->GetTrack(pt_index)->IsPrimary() ) break;

		//Note downwards beta enforced by Event->GetPrimaryBeta() (should be positive) multiplied by Event->GetPrimaryMomentumDirection()[2] (z trajectory of particle)
		if( pt != nullptr && ( (TRG == 0) || ((int)Event->GetTriggerSources().at(0) == TRG) ) && pt->GetChi2()/pt->GetNdof() < 3.2 && -fabs(Event->GetPrimaryMomentumDirection().CosTheta()) > -coslow && -fabs(Event->GetPrimaryMomentumDirection().CosTheta()) < -coshigh && Event->GetPrimaryBeta()*Event->GetPrimaryMomentumDirection()[2] < 0 && fabs(MCEvent->GetPrimaryBeta()) > betacut && fabs(MCEvent->GetPrimaryBeta()) <  betahigh ){
			//-----------EVENT LEVEL CUT APPLIED

			//cout << "Event: " << i << endl;
			//First iteration over event for flags and EnergyDepositionMip vector filling
			double Bgen = MCEvent->GetPrimaryBeta();
			//double gamma = sqrt( 1/(1 - pow(Bgen,2)) );
			//cout << "Bgen " << Bgen << endl;

			for(uint isig=0; isig<Event->GetTrack(0)->GetEnergyDeposition().size(); isig++){
                unsigned int VolumeId  = Event->GetTrack(0)->GetVolumeId(isig); //Check the VolumeId of the event
                if(volspec(VolumeId,0,2) == 20 && Event->GetTrack(0)->GetEnergyDeposition(isig) > TrackerCut){
                    int layer = GGeometryObject::GetTrackerLayer(VolumeId);
                    Layer_Hits_Tracker[layer]++;
                }
                if(volspec(VolumeId,0,3) == 100){ Umbflag = 1; } // cout << "UMB hit!" <<endl ;
                if(volspec(VolumeId,0,3) == 110) { CBEtopflag = 1; }// cout << "CBE top hit!" << endl;
                if(volspec(VolumeId,0,3) == 111) { CBEbotflag = 1; }// cout << "CBE bot hit!" << endl;
			}

			//Do we want to only run this on certain tracks? Yeah probably. Can remove the TOF flags. Just run on whatever lol.
			if(Umbflag && CBEtopflag /*&& CBEbotflag && 1*/){
				vector<double> Edep_TKR;

				for(uint isig=0; isig<Event->GetTrack(0)->GetEnergyDeposition().size(); isig++){
				    unsigned int VolumeId  = Event->GetTrack(0)->GetVolumeId(isig); //Check the VolumeId of the event
					if(volspec(VolumeId,0,2) == 20 && Event->GetTrack(0)->GetEnergyDeposition(isig) > TrackerCut){
						int layer = GGeometryObject::GetTrackerLayer(VolumeId);
						int sdstrip = GGeometryObject::GetDetectorStrip(VolumeId);
				        if(sdstrip != 0 && sdstrip != 7 && Layer_Hits_Tracker[layer] < 2){
			                Edep_TKR.push_back(Event->GetTrack(0)->GetEnergyDeposition(isig)*fabs(Event->GetPrimaryMomentumDirection().CosTheta()));
						}
				    }
				}

				if(print){
				for(uint isig=0; isig<Event->GetTrack(0)->GetEnergyDeposition().size(); isig++){
                    unsigned int VolumeId  = Event->GetTrack(0)->GetVolumeId(isig); //Check the VolumeId of the event
                    cout << "Energy deposition is " << Event->GetTrack(0)->GetEnergyDeposition(isig)*fabs(Event->GetPrimaryMomentumDirection().CosTheta()) << " at " << VolumeId << endl;
                }
				}

                if(Edep_TKR.size() > 2){
            	    if(print)cout << endl << "Event is " << i << endl;
                    if(print)cout << "Primary particle is " << MCEvent->GetTrack(0)->GetPdg() << endl;
                    if(print)cout << "Gen Beta is " << fabs(MCEvent->GetPrimaryBeta()) << endl;
                    if(print)cout << "Reco beta is " << fabs(Event->GetPrimaryBeta()) << endl;
                    std::sort(Edep_TKR.begin(), Edep_TKR.begin()+Edep_TKR.size());

                    if(print){
                    for(unsigned int isig = 0; isig < double(Edep_TKR.size()); isig++){
                        cout << "Sorted angle-corrected tracker hit " << isig << " is " << Edep_TKR.at(isig) << endl;
                    }
                    }

                    double TrEdep = 0;
                    double CtrTrEdep = 0;
                    for(unsigned int isig = 0; isig < double(Edep_TKR.size())/2; isig++){
                    //for(unsigned int isig = 0; isig < Edep_TKR.size(); isig++){
                        TrEdep += Edep_TKR.at(isig);
                        CtrTrEdep++;
                    }
                    if(CtrTrEdep == 0){ TrEdep = 0;
                    }else{ TrEdep /= CtrTrEdep;}
                    if(print)cout << "Truncated Edep = " << TrEdep << endl;

                    if(TrEdep > solvermin){
                        double root = Z1_tkr_Solve->GetX(TrEdep);
                        if(print)cout << "Gen Beta is " << fabs(MCEvent->GetPrimaryBeta()) << endl;
                        if(print)cout << "Reco beta is " << fabs(Event->GetPrimaryBeta()) << endl;
                        if(print)cout << "Calculated Beta: " << root << endl;
                        HRecB_vs_GenB->Fill(Event->GetPrimaryBetaGenerated(),Event->GetPrimaryBeta());
                        HProxB_vs_GenB->Fill(Event->GetPrimaryBetaGenerated(),root);
                    }
                } //Closed bracked for Truncated energy deposition calculation, requiring 3 or more tracker hits.

			} //Closed bracket for if statement for TOF cuts

			//-----------EVENT LEVEL CUTS END


		} //Closed bracket for event level cut (beta, cos, pt exists)

	} //Closed bracket for single track cut

}  //Closed bracket for iteration through tree events, move on to the next event i

histplot2d("c1",HRecB_vs_GenB,"Rec_B versus Gen_B","Generated Beta", "Reconstructed Beta","NEntries", out_path + "GenBRecB" + roundstr_d(betacut,2) + "-" + roundstr_d(betahigh,2) );
histplot2d("c2",HProxB_vs_GenB,"Prox_B versus Gen_B","Generated Beta", "Proxy Beta","NEntries", out_path + "GenBProxB" + roundstr_d(betacut,2) + "-" + roundstr_d(betahigh,2) );

myfile.open(out_path + txtname,std::ios::app);
myfile << "Total Events/Mainscale Factor " << TreeRec->GetEntries()/MainLoopScaleFactor << endl;
myfile.close();


//Histogram section
//--------------------------------------

cout << endl << "I am done" << endl;

return 1;

}
