//I am branching out and trying to do a thing in a way that is smart!
//This code is starting by using an intelligent framework for counting numbers of events that pass certain cuts!

//How to use:
// ./BPRecCutEff -i /home/kelsey/simulations/simdat/ground/251204/25.10/ethernet251204_082 -o test

using namespace std;

#include "KYtools.C"

using namespace Crane::Analysis;
namespace ca = Crane::Analysis;
namespace cl = Crane::Common;

double EdepCutLow = 0.4; //For the efficiency cuts, just do 0.4 MeV to keep things consistent
double betacut = 0.2;
double betahigh = 1.2;
double coslow = 1;
double coshigh = 0.54;

float Ztof = 5.574;
float Atof = 10.3;
float Ltof = 0.635;
float rtof = 1.032;
float iontof = (0.000016 * pow(Ztof,0.9));

float z = 1;
float z2 = 2;
float Ztkr = 14;
float Atkr = 28;
float rtkr = 2.33; //Density g/cm^2
float Ltkr = 0.22; //For instrument tracker, set Ltkr = 0.22

float iontkr = (0.000016 * pow(Ztkr,0.9));
float C_1 = 0.3071/2; //MeV/ g/cm^2 #2*pi*constants not 4*pi*constants for MPV
float me = 0.511; //mass of electron * c^2

double ZOne_TKR(double x)
{
    return C_1*pow(z,2)*(rtkr*Ltkr*Ztkr/Atkr)*(1/pow(x,2))* ( log( (1.022 * pow(x,2)/(1 - pow(x,2)) )/iontkr ) + log(C_1*pow(z,2)*(rtkr*Ltkr*Ztkr/Atkr)*(1/pow(x,2))/iontkr)  + 0.2 - pow(x,2) );
}

double ZOne_TOF(double x)
{
    return C_1*pow(z,2)*(rtof*Ltof*Ztof/Atof)*(1/pow(x,2))* ( log( (1.022 * pow(x,2)/(1 - pow(x,2)) )/iontof ) + log(C_1*pow(z,2)*(rtof*Ltof*Ztof/Atof)*(1/pow(x,2))/iontof)  + 0.2 - pow(x,2) );
}

float solvermin_tkr = ZOne_TKR(0.95); //0.787;
float solvermin_tof = ZOne_TOF(0.95);

//Try figuring out why this worked: https://root.cern/doc/v620/classTF1.html
TF1 *Z1_tkr_Solve = new TF1("Z1_tkr_Solve", [](double *x, double *p){ return ZOne_TKR(x[0]); }, 0.01, 0.95, 0);
TF1 *Z1_tof_Solve = new TF1("Z1_tof_Solve", [](double *x, double *p){ return ZOne_TOF(x[0]); }, 0.01, 0.95, 0);

bool charge_cut_z1(double Proxy_Beta, double Reco_Beta){
    if(Reco_Beta > 0.75 && Proxy_Beta > 0.66){ return 1; }else if(Reco_Beta < 0.75 && Proxy_Beta >  0.83*Reco_Beta + 0.0328){return 1;}else{return 0;}
}

bool charge_cut_z2(double Proxy_Beta, double Reco_Beta){
    if(Reco_Beta > 0.75 && Proxy_Beta < 0.55){ return 1; }else if(Reco_Beta < 0.75 && Proxy_Beta <  0.78*Reco_Beta - 0.033 && Proxy_Beta >  0.067*Reco_Beta + 0.12){return 1;}else{return 0;}
}

//Use TOF flags gathered for
void tof_efficiencies(int cuts[], int& k, bool& pass, int UMBflag, int CORflag, int CBEtopflag, int CBEbotflag, int CBEsideflag, int TKRflag){

    int init = k;
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
    k = init + 3; //Number of if statements.

}

//TOF flags are marked
//There is a slightly (better?) way to do this maybe, see MC_Beta_Timning_Res
void tof_flags(const CEventRec* Event, int& UMBflag, int& CORflag, int& CBEtopflag, int& CBEbotflag, int& CBEsideflag, int& TKRflag, int Layer_Hits_Tracker[] ){
    for(uint isig=0; isig<Event->GetTrack(0)->GetEnergyDeposition().size(); isig++){
        unsigned int VolumeId  = Event->GetTrack(0)->GetVolumeId(isig);
        if(Event->GetTrack(0)->GetEnergyDeposition(isig) > EdepCutLow){
            if(volspec(VolumeId,0,3) == 100)UMBflag++;
            if(volspec(VolumeId,0,3) == 110)CBEtopflag++;
            if(volspec(VolumeId,0,3) == 111)CBEbotflag++;
            if(volspec(VolumeId,0,3) == 112 || volspec(VolumeId,0,3) == 113 || volspec(VolumeId,0,3) == 114 || volspec(VolumeId,0,3) == 115 || volspec(VolumeId,0,3) == 116)CBEsideflag++;
            if(volspec(VolumeId,0,3) == 102 || volspec(VolumeId,0,3) == 103 || volspec(VolumeId,0,3) == 104 || volspec(VolumeId,0,3) == 105 || volspec(VolumeId,0,3) == 106)CORflag++;
            if(volspec(VolumeId,0,2) == 20 && Event->GetTrack(0)->GetEnergyDeposition(isig) > EdepCutLow){
                int layer = GGeometryObject::GetTrackerLayer(VolumeId);
                Layer_Hits_Tracker[layer]++;
                TKRflag = 1;
            }
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

    int init = k;

    if( (int)Event->GetTriggerSources().at(0) == 2 ){ //Track trigger
        cuts[k]++; k++;
        if(pt != nullptr){ //PT fount
            cuts[k]++; k++;
            if(fabs(Event->GetPrimaryBeta()) >  betacut && fabs(Event->GetPrimaryBeta()) <  betahigh){ //Beta range
                cuts[k]++; k++;
                if(Event->GetNTracks() == 1){ //Single Track
                    cuts[k]++;k++;
                    if(Event->GetPrimaryBeta()*Event->GetPrimaryMomentumDirection()[2] < 0){ //Down beta
                        cuts[k]++;k++;
                        if(-fabs(Event->GetPrimaryMomentumDirection().CosTheta()) > -coslow && -fabs(Event->GetPrimaryMomentumDirection().CosTheta()) < -coshigh){ //Cos range
                            cuts[k]++;k++;
                            if( pt->GetChi2()/pt->GetNdof() < 3.2 ){ //Chi^2/Ndof cut
                                cuts[k]++;k++;
                                pass = true;
                                //cout << "Event passes so many things!! " << endl;
                            }
                        }
                    }
                }
            }
        }
    }

    k = init + 7; //Number of if statements.

}


int main(int argc, char *argv[]){

GOptionParser* parser = GOptionParser::GetInstance();
parser->AddProgramDescription("Minimal Reproducable Example for Extracing Data from Reco Data");
parser->AddCommandLineOption<string>("in_path", "path to instrument data files", "./*", "i");
parser->AddCommandLineOption<string>("out_file", "name of output root file", "", "o");
parser->AddCommandLineOption<int>("TRG", "Which trigger?",0,"r");
parser->AddCommandLineOption<int>("MainloopScale", "Main loop scale factor",1,"m");
parser->AddCommandLineOption<double>("tkr_factor", "tkr factor",1,"k");
parser->AddCommandLineOption<double>("tof_factor", "tof factor",1,"f");
parser->AddCommandLineOption<bool>("print", "print out info?",0,"n");
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

double tk = parser->GetOption<double>("tkr_factor");
double tf = parser->GetOption<double>("tof_factor");
int print = parser->GetOption<bool>("print");

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
    "Track trigger ",
    "Reco Total Events with Primary Track ",
    "Reco Total Events with PT found, in Beta Range ",
    "Reco PT yes, Beta Range, Single Track Events ",
    "Reco Down beta ",
    "Reco Cos(theta) in range ",
    "Chi^2/Ndof < 3.2 ",
    "Reco all above YES UMB ",
    "YES UMB, CT ",
    "YES UMB, CT, (CB, TKR, or CS) ",
};
const int NCuts = size(cutnames) ;
int cuts[NCuts] = {}; //An array that counts the number of cuts
int k = 0;
cuts[k] = TreeRec->GetEntries(); //Number of events (not a cut)
k++;

//For this plot, only want Beta Proxy and Beta reconstruction
TH1F * HBetaProxy = new TH1F("HBetaProxy","HBetaProxy",40, betacut, betahigh);
TH1F * HBetaRec = new TH1F("HBetaRec","HBetaRec",40, betacut, betahigh);
TH2D * HRecB_vs_ProxB = new TH2D("HRecB_vs_ProxB","Rec_Beta vs Prox_B",50,betacut - 0.1, betahigh + 0.1, 50, 0.1 , 1);

int total_bp = 0;
int proton_counter = 0;
int alpha_counter = 0;

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
//for(unsigned int i = 0; i < 1000; i+=MainLoopScaleFactor){
//for(unsigned int i = 60500; i < 90600; i+=MainLoopScaleFactor){
for(unsigned int i = 0; i < TreeRec->GetEntries()/MainLoopScaleFactor; i++){ //This is not the "correct" way to do this, but it's probably fine. Should be skipping M each time, but that seems to be really slow!!
    TreeRec->GetEntry(i);

	if( ((int)i % (int)ceil(TreeRec->GetEntries()/(MainLoopScaleFactor*10))) == 0){
	    cout << "Event number " << i << endl;
		cout << Event->GetActiveReconstruction() << endl;
	}

	bool pass = 0;
	k = 1;
	//cout << "Event is " << i << endl;

	std_cuts(Event,pass,k,cuts);
	if(pass){

        int UMBflag = 0;
        int CORflag = 0;
        int CBEtopflag = 0;
        int CBEbotflag = 0;
        int CBEsideflag = 0;
        int TKRflag = 0;
        int Layer_Hits_Tracker[7] = {}; //Seven layers

        bool eff_pass= 0; //Passing the TOF efficienciencies

        tof_flags(Event, UMBflag, CORflag, CBEtopflag, CBEbotflag, CBEsideflag, TKRflag, Layer_Hits_Tracker);
        //cout << "event is " << i << endl;
        tof_efficiencies(cuts,k,eff_pass,UMBflag, CORflag, CBEtopflag, CBEbotflag, CBEsideflag, TKRflag);
        //Passing of TOF cuts determined here, so can move to do Beta Proxy if desired.
        //cout << "Did is pass the pass? " << eff_pass << endl;

        //Beta Proxy calculation starts


        if(print) cout << endl << "Event is " << i << endl;
		if(print) cout << "Breco is " << Event->GetPrimaryBeta() << endl;
		vector<double> Beta_Proxy;

		if(eff_pass){
		for(uint isig=0; isig<Event->GetTrack(0)->GetEnergyDeposition().size(); isig++){
		    unsigned int VolumeId  = Event->GetTrack(0)->GetVolumeId(isig); //Check the VolumeId of the event
			float Edep = Event->GetTrack(0)->GetEnergyDeposition(isig);
			float Ang_Edep = 0;

			if(GGeometryObject::IsTofVolume(VolumeId) && Event->GetTrack(0)->GetEnergyDeposition(isig) > EdepCutLow){
				if(volspec(VolumeId,2,1) == 0 || volspec(VolumeId,2,1) == 1){ //It's a flat paddle.
				    Ang_Edep = tf*Edep*fabs(Event->GetPrimaryMomentumDirection().CosTheta());
					if(Ang_Edep > solvermin_tof){
					    double root = Z1_tof_Solve->GetX(Ang_Edep);
						if(print)cout << "Hit is " << isig << endl;
						if(print)cout << "Energy Deposition at Volid " << VolumeId << " is " << Edep << " angle corrected is " << Ang_Edep << endl;
                        if(print)cout << "Calculated Beta at Volid = "  << root << endl;
                         Beta_Proxy.push_back(root);
					}
                }else{ //It's a vertical paddle
                    //The calculation of the correction factor depends on the orientation of the paddle
                    float xy_path_corr = 1;
                    float x_vec = Event->GetTrack(0)->GetMomentumDirection()[0][0];
                    float y_vec = Event->GetTrack(0)->GetMomentumDirection()[0][1];

                    if(volspec(VolumeId,2,1) == 2 ||volspec(VolumeId,2,1) == 3){ //+/- X paddles
                        if(print)cout << "xvec, yvec = " << x_vec << " " << y_vec << endl;
                        if(print)cout << "cos(phi) = " << fabs(x_vec) / sqrt( pow(x_vec,2) + pow(y_vec,2) ) << endl;
                        xy_path_corr = fabs(x_vec) / sqrt( pow(x_vec,2) + pow(y_vec,2) );
                    }
                    if(volspec(VolumeId,2,1) == 4 ||volspec(VolumeId,2,1) == 5){ //+/- Y paddles
                        if(print)cout << "xvec, yvec = " << x_vec << " " << y_vec << endl;
                        if(print)cout << "cos(phi) = " << fabs(y_vec) / sqrt( pow(x_vec,2) + pow(y_vec,2) ) << endl;
                        xy_path_corr = fabs(y_vec) / sqrt( pow(x_vec,2) + pow(y_vec,2) );
                    }

                    Ang_Edep = tf*Edep*fabs( sqrt(1-pow(Event->GetPrimaryMomentumDirection().CosTheta(),2)))*xy_path_corr;
                    if(Ang_Edep > solvermin_tof){
						double root = Z1_tof_Solve->GetX(Ang_Edep);
						if(print)cout << "Hit is " << isig << endl;
						if(print)cout << "Energy Deposition at Volid " << VolumeId << " is " << Edep << " angle corrected is " << Ang_Edep << endl;
                        if(print)cout << "Calculated Beta = "  << root << endl;
                        Beta_Proxy.push_back(root);
					}
                }
			} //End TOF BP Calculation

			if(volspec(VolumeId,0,2) == 20 && Event->GetTrack(0)->GetEnergyDeposition(isig) > EdepCutLow){
				int layer = GGeometryObject::GetTrackerLayer(VolumeId);
				int sdstrip = GGeometryObject::GetDetectorStrip(VolumeId);
				if(sdstrip != 0 && sdstrip != 7 && Layer_Hits_Tracker[layer] < 2){
			        Ang_Edep = tk*Edep*fabs(Event->GetPrimaryMomentumDirection().CosTheta());
					if(Ang_Edep > solvermin_tkr){
				        double root = Z1_tkr_Solve->GetX(Ang_Edep);
                        Beta_Proxy.push_back(root);
					}
				}
			} //End TKR BP Calculation

		} //Closed bracket for iteration over isig with TOF cuts
		} //End efficiencies TOF cuts pass

		if(Beta_Proxy.size() > 2){
		    std::sort(Beta_Proxy.begin(), Beta_Proxy.end(),std::greater<>());

            double TrBP = 0;
            double CtrTrBP = 0;
            //for(unsigned int isig = 0; isig < floor(double(Beta_Proxy.size())/2); isig++){
            for(unsigned int isig = 0; isig < double(Beta_Proxy.size())/2; isig++){
                TrBP += Beta_Proxy.at(isig);
                CtrTrBP++;
            }if(CtrTrBP == 0){ TrBP = 0;
            }else{ TrBP /= CtrTrBP;}

            //To floor or not to floor!!!
            //cout << "Beta Prox floor(BP/2) = " << TrBP << endl;
            if(print)cout << "Beta Prox BP/2 not floor = " << TrBP << endl;
            total_bp++;
            if(TrBP > 0 && charge_cut_z1(TrBP,Event->GetPrimaryBeta()) ){
                proton_counter++;
                HBetaProxy->Fill(TrBP);
                HBetaRec->Fill(Event->GetPrimaryBeta());
                HRecB_vs_ProxB->Fill(Event->GetPrimaryBeta(),TrBP);
            }
            if(TrBP > 0 && charge_cut_z2(TrBP,Event->GetPrimaryBeta()) ){
                alpha_counter++;
            }

		} //Closed bracket Beta Proxy calculation



        //Beta Proxy calculation ends


			//-----------EVENT LEVEL CUTS END
	} //Closed bracket for standard event selection cuts.

}  //Closed bracket for iteration through tree events, move on to the next event i

HBetaProxy->Scale( 1./HBetaProxy->Integral(),"WIDTH");
//HBetaProxy->SetMaximum(bcounts);
histplot1f("c1",HBetaProxy,"Proxy Beta","Proxy Beta","NEntries", out_path + "Rec_BetaProxy"+ "B" + roundstr_d(betacut,2) + "-" + roundstr_d(betahigh,2)+ "TF_Factor" + roundstr_d(tf,2) + "TK_Factor" + roundstr_d(tk,2)  );
histplot1f("c2",HBetaRec,"Reconstructed Beta","Reconstructed Beta","NEntries", out_path + "Rec_BetaRec " + "B" + roundstr_d(betacut,2) + "-" + roundstr_d(betahigh,2) + "TF_Factor" + roundstr_d(tf,2) + "TK_Factor" + roundstr_d(tk,2)  );
histplot2d("c2_5",HRecB_vs_ProxB,"Prox_B versus Rec_B","Reconstructed Beta", "Proxy Beta","NEntries", out_path + "RecBProxB" + "B" +  roundstr_d(betacut,2) + "-" + roundstr_d(betahigh,2)+ "TF_Factor" + roundstr_d(tf,2) + "TK_Factor" + roundstr_d(tk,2) );


myfile.open(out_path + title,std::ios::app); //Open the file for modifying
for(int k = 0; k < NCuts; k++){
    myfile << cutnames[k] << cuts[k] << endl;
    if(k > 0) myfile << fixed << setprecision(2) << 100*(float)cuts[k]/(float)cuts[k-1] <<"%" <<  endl; //NEntries doesn't need to divide by anything
}
myfile.close();

cout << "Super rough Nprotons " << proton_counter << endl;
cout << "Super rough Nalphas " << alpha_counter << endl;
cout << "Super rough Nprotons + Nalphas " << proton_counter + alpha_counter << endl;

cout << endl << "I am done" << endl;

return 1;

}
