//This script probably not necessary.
//Hey we can't have nice things.
//So the first version of this code uses Event->GetPrimaryBetaGenerated() for the generated beta (logical cool and good.)
//I still want to look at v210, which DOES NOT have MCEvent, but DOES have Event->GetPrimaryBetaGenerated() so here we flipping go.
//How to use: ./V2_MCCutEff -i /home/kelsey/simulations/simdat/mu/v.2.1.2/mu-_gaps_triggerlevel1_FTFP_BERT_1744342800_rec -w 1

//Okay horrifyingly v210 does NOT have anything good put into Event->GetPrimaryBetaGenerated()
//The last thing I could do is weight by reconstruction, but ABSOLUTELY NOT

using namespace std;

#include "KYtools.C" //Everything is already included here lol

int main(int argc, char *argv[]){

GOptionParser* parser = GOptionParser::GetInstance();
parser->AddProgramDescription("Efficiencies from cuts at MC");
parser->AddCommandLineOption<string>("in_path", "path to instrument data files", "./*", "i");
parser->AddCommandLineOption<string>("out_file", "name of output root file", "", "o");
parser->AddCommandLineOption<double>("beta_low", "low Beta Cut",0.8,"l");
parser->AddCommandLineOption<double>("beta_high", "upper Beta Cut",1,"u");
parser->AddCommandLineOption<bool>("MC_Weighting", "Monte Carlo Weighting on or off?",0,"w");
parser->AddCommandLineOption<int>("MainloopScale", "Main loop scale factor",1,"m");
parser->ParseCommandLine(argc, argv);
parser->Parse();

bool MC_Weight = parser->GetOption<bool>("MC_Weighting");

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
myfile.open(out_path + "MCTruthCuts_V2.txt");
myfile << TString::Format( "Filename : %s", reco_path.c_str() )  << endl;
myfile << TString::Format( "Beta High : %f", betahigh) << endl;
myfile << TString::Format( "Beta Low : %f", betacut) << endl << endl;
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

//int MainLoopScaleFactor = 1; //Set this number to scale the step size. Larger means runs faster and fewer events
double TofCutLow = 0; //No low Tof cut right now
double TrackerCut = 0.3; //Threshold for an energy deposition to be considered a hit

double xlow = 0.3; //Low range for histogram MeV
double xhigh = 5; //High range for histogram MeV

double coshigh = 0.54; //0.995; //0.92 //0.54 is the highest angle that can hit UMB, CBEtop, CBEbot
double coslow = 1; //0.62 //0.8
double acutlow = 1.5;
const Int_t NBins = 50;

//From fitted peaks from flight! Check over!
double TofAngleCorrectedMip = 0.82;
double TrackerAngleCorrectedMip = 0.57;

int passctr = 0;
int strkctr = 0;

int MCpassB = 0;
int MCmust = 0; //Counter for MC true single tracks
int MCnosides = 0; //Counter for MC true single tracks that pass the track trigger that have no side hits
int MCyesumb_ctop_cbot = 0; //Counter for MC true single tracks that pass the track trigger that have no side hits YES umb, cbe_top, cbe_bot hits
int MCyesumb_ctop_nocbot = 0; //Counter for MC true single tracks that pass the track trigger that have no side hits YES umb, cbe_top, NO cbe_bot hits

int must = 0; //Counter for reco single tracks
int nosides = 0; //Counter for reco single tracks that pass the track trigger that have no side hits
int yesumb_ctop_cbot = 0; //Counter for reco single tracks that pass the track trigger that have no side hits YES umb, cbe_top, cbe_bot hits
int yesumb_ctop_nocbot = 0; //Counter for reco single tracks that pass the track trigger that have no side hits YES umb, cbe_top, NO cbe_bot hits


float Total_Weighted = 0;
float MCpassB_Weighted = 0;
float MCmust_Weighted  = 0; //Counter for MC true single tracks
float MCnosides_Weighted  = 0; //Counter for MC true single tracks that pass the track trigger that have no side hits
float MCyesumb_ctop_cbot_Weighted  = 0; //Counter for MC true single tracks that pass the track trigger that have no side hits YES umb, cbe_top, cbe_bot hits
float MCyesumb_ctop_nocbot_Weighted  = 0; //Counter for MC true single tracks that pass the track trigger that have no side hits YES umb, cbe_top, NO cbe_bot hits

float must_Weighted  = 0; //Counter for reco single tracks
float nosides_Weighted  = 0; //Counter for reco single tracks that pass the track trigger that have no side hits
float yesumb_ctop_cbot_Weighted  = 0; //Counter for reco single tracks that pass the track trigger that have no side hits YES umb, cbe_top, cbe_bot hits
float yesumb_ctop_nocbot_Weighted  = 0; //Counter for reco single tracks that pass the track trigger that have no side hits YES umb, cbe_top, NO cbe_bot hits



//MC Weighting things below:
ca::GPlottingTools Plotting;
TChain*TreeSimulationParameter = new TChain("SimulationParameterTree");
TreeSimulationParameter->Add(FilenameRoot);
GSimulationParameter * Parameter = new GSimulationParameter;
TreeSimulationParameter->SetBranchAddress("SimulationParameter", &Parameter);
TreeSimulationParameter->GetEntry(0);

double FluxScaleFactor = 71.1552/32.058750*0.25;

int BetaBins = 25;
double StartingPlaneAcceptance = 1;
TH1D* HPrimaryBeta = nullptr;
double BinWidthFactor = 1;
std::vector<double> PrimaryBetaLowHigh;

vector<pair<double, double> > CosZenithCut;
CosZenithCut.push_back(make_pair(-0.75, -1));
CosZenithCut.push_back(make_pair(-0.5, -0.75));
CosZenithCut.push_back(make_pair(-0.25, -0.5));
CosZenithCut.push_back(make_pair(0, -0.25));

vector<TGraph*> GMuonTotalFluxUnscaled;
GMuonTotalFluxUnscaled.push_back(Plotting.ConvertEnergyFluxToBetaFlux(Plotting.GetTH1D(getenv("GAPS") + string("/resources/fluxes/total_fluxes_coszenith_100_m_-13_antarctica.root"), "c6a", "p_total_altitude_zenith_energy_0.875"), 0.1057));
GMuonTotalFluxUnscaled.push_back(Plotting.ConvertEnergyFluxToBetaFlux(Plotting.GetTH1D(getenv("GAPS") + string("/resources/fluxes/total_fluxes_coszenith_100_m_-13_antarctica.root"), "c6a", "p_total_altitude_zenith_energy_0.625"), 0.1057));
GMuonTotalFluxUnscaled.push_back(Plotting.ConvertEnergyFluxToBetaFlux(Plotting.GetTH1D(getenv("GAPS") + string("/resources/fluxes/total_fluxes_coszenith_100_m_-13_antarctica.root"), "c6a", "p_total_altitude_zenith_energy_0.375"), 0.1057));
GMuonTotalFluxUnscaled.push_back(Plotting.ConvertEnergyFluxToBetaFlux(Plotting.GetTH1D(getenv("GAPS") + string("/resources/fluxes/total_fluxes_coszenith_100_m_-13_antarctica.root"), "c6a", "p_total_altitude_zenith_energy_0.125"), 0.1057));

CAnalysisManager AnalysisManagerRec;
AnalysisManagerRec.SetGSimulationParameterTChain(TreeSimulationParameter);
StartingPlaneAcceptance = AnalysisManagerRec.GetStartingPlaneAcceptance();
PrimaryBetaLowHigh = AnalysisManagerRec.GetPrimaryBetaLowHigh();
HPrimaryBeta = AnalysisManagerRec.GetHPrimaryBeta();
BinWidthFactor = (PrimaryBetaLowHigh.at(1)-PrimaryBetaLowHigh.at(0))/double(BetaBins) / HPrimaryBeta->GetBinWidth(1);

//End MC Weighting things setup!


//Prepare cuts:
map<int, unsigned int> TofIndexVolumeIdMap;

//For Plotting purposes
//ca::GPlottingTools Plotting;
//char text[400]; //This variable is used later to name the plots

//TH1D * hedep;
//hedep = new TH1D ("h0", ("Edep l Beta " + to_string(betacut) + " - " + to_string(betahigh) ).c_str(), NBins, xlow,xhigh);

//How many entries:
cout << "Total Number of events / Mainscale Factor = " << TreeRec->GetEntries()/MainLoopScaleFactor << endl;
cout << "Beta high " << betahigh << endl;
cout << "Beta cut " << betacut << endl;
//cout << "MPV range on 2D Full Tracker will be " << mpvmin << " - " << mpvmax << endl;

//Now we can go over the loop
TreeRec->GetEntry(0);
TreeMC->GetEntry(0);

cout << "Total Number of events / Mainscale Factor = " << TreeRec->GetEntries()/MainLoopScaleFactor << endl;

//Using i to loop over every event in the tree
//for(unsigned int i = 0; i < 100; i+=MainLoopScaleFactor){
//for(unsigned int i = 60500; i < 90600; i+=MainLoopScaleFactor){
for(unsigned int i = 0; i < TreeRec->GetEntries()/MainLoopScaleFactor; i++){ //This is not the "correct" way to do this, but it's probably fine. Should be skipping M each time, but that seems to be really slow!!
    TreeRec->GetEntry(i);
    TreeMC->GetEntry(i);

	//Cuts are implemented in this chunk:

	if( ((int)i % (int)ceil(TreeRec->GetEntries()/(MainLoopScaleFactor*10))) == 0){
	    cout << "Event number " << i << endl;
		cout << Event->GetActiveReconstruction() << endl;
	}

	//Lead with MC weighting, EACH EVENT has a weighting based on the ground flux of muons, angular and beta weighted!
        double RateScale = 0;
        if( MC_Weight ){
            //Theoretically save on compute time if put weighting stuff at last possible place.
            double AcceptanceScale;
            if (TreeSimulationParameter != nullptr){
                //acceptance scaling factor based on beta of the primary
                AcceptanceScale = MainLoopScaleFactor*StartingPlaneAcceptance/(BinWidthFactor*HPrimaryBeta->GetBinContent(HPrimaryBeta->FindBin(Event->GetPrimaryBetaGenerated())));

                if(HPrimaryBeta->GetBinContent(HPrimaryBeta->FindBin(Event->GetPrimaryBetaGenerated())) == 0) AcceptanceScale = 0;
                //Find the bin associated with the generated beta in the vector, see if it doesn't exist?
            }else AcceptanceScale = 1; //Set it = to 1 if there's a nullptr? That's a surprise...
            //cout << "AcceptanceScale = " << AcceptanceScale << endl;

       	    int AngularRegion = -1;
            for(unsigned int a = 0; a < CosZenithCut.size(); a++) if(Event->GetPrimaryMomentumDirectionGenerated().CosTheta() < CosZenithCut.at(a).first && Event->GetPrimaryMomentumDirectionGenerated().CosTheta() > CosZenithCut.at(a).second) AngularRegion = a;
            //This is just checking which "bin" the generated cos(theta) is in
            if(AngularRegion < 0) continue; //Don't bother if angular region wasn't found
            RateScale = FluxScaleFactor*AcceptanceScale*GMuonTotalFluxUnscaled.at(AngularRegion)->Eval(Event->GetPrimaryBetaGenerated());
            //cout << "Event is " << i << " Beta is " << Event->GetPrimaryBetaGenerated() << " Cos(theta) is " << Event->GetPrimaryMomentumDirectionGenerated().CosTheta() << " RateScale?? = " << RateScale << endl;
        }
        Total_Weighted = Total_Weighted + RateScale;


	CTrackRec* pt = Event->GetPrimaryTrack();
	uint pt_index = 0;
    for( ; pt_index < Event->GetNTracks(); pt_index++) if( Event->GetTrack(pt_index)->IsPrimary() ) break;

	//Note downwards beta enforced by Event->GetPrimaryBeta() (should be positive) multiplied by Event->GetPrimaryMomentumDirection()[2] (z trajectory of particle)
	//if(pt != nullptr && fabs(Event->GetPrimaryBeta()) >  betacut && fabs(Event->GetPrimaryBeta()) <  betahigh ){
	if(pt != nullptr && fabs(Event->GetPrimaryBetaGenerated()) >  betacut && fabs(Event->GetPrimaryBetaGenerated()) <  betahigh ){
	//if((pt != nullptr && pt->GetChi2()/pt->GetNdof()) < 3.2 && -fabs(Event->GetPrimaryMomentumDirectionGenerated().CosTheta()) > -coslow && -fabs(Event->GetPrimaryMomentumDirectionGenerated().CosTheta()) < -coshigh && Event->GetPrimaryBetaGenerated()*MCEvent->GetPrimaryMomentumDirection()[2] < 0 && fabs(Event->GetPrimaryBetaGenerated()) >  betacut && fabs(Event->GetPrimaryBetaGenerated()) <  betahigh ){
	//if(pt != nullptr && -fabs(Event->GetPrimaryMomentumDirection().CosTheta()) > -coslow && -fabs(Event->GetPrimaryMomentumDirection().CosTheta()) < -coshigh && Event->GetPrimaryBeta()*MCEvent->GetPrimaryMomentumDirection()[2] < 0 && fabs(Event->GetPrimaryBeta()) >  betacut && fabs(Event->GetPrimaryBeta()) <  betahigh ){
	//if(pt != nullptr && -fabs(Event->GetPrimaryMomentumDirectionGenerated().CosTheta()) > -coslow && -fabs(Event->GetPrimaryMomentumDirectionGenerated().CosTheta()) < -coshigh && Event->GetPrimaryBetaGenerated()*Event->GetPrimaryMomentumDirectionGenerated()[2] < 0 && fabs(Event->GetPrimaryBetaGenerated()) >  betacut && fabs(Event->GetPrimaryBetaGenerated()) <  betahigh ){

        //MC weighting over
        //cout << "Event is " << i << " Beta is " << Event->GetPrimaryBetaGenerated() << " Cos(theta) is " << Event->GetPrimaryMomentumDirectionGenerated().CosTheta() << " RateScale?? = " << RateScale << endl;
	    MCpassB++;
		MCpassB_Weighted = MCpassB_Weighted + RateScale;

		//-----------EVENT LEVEL CUT APPLIED
		//cout << endl << "Event " << i << endl;
		//cout << "MC info " << endl;
		//cout << "Number of tracks " << MCEvent->GetNTracks() << endl;

		//Might make sense to put this into KYtools as the MC Event framework
		int pid = MCEvent->GetTrack(0)->GetPdg(); //Identity of the primary track should be this
	    vector<int> MCSpec;
		vector<unsigned int> MCVolid;
		vector<double> MCEdep;
		bool MCST = 0; //MC Single track flag

		    for(uint t = 0; t < MCEvent->GetNTracks();t++){ //Each particle has its own track
                //cout << "Track is " << t << endl;
                for(uint isig=0; isig<MCEvent->GetTrack(t)->GetEnergyDeposition().size(); isig++){ //Iterate over the energy depositions of the track
                    unsigned int VolumeId  = MCEvent->GetTrack(t)->GetVolumeId(isig); //Check the volumeId of each of the hits in the iteration
                    int index = -1; //Note the index = -1 right now

                    for (int k = 0; k < MCVolid.size(); k++) { //Iterate over the VolumeId vector
                        if (MCVolid[k] == VolumeId) { //If volumeid is in the vector already, add the energy deposition to the vector counting edeps there
                            MCEdep[k] = MCEdep[k] + MCEvent->GetTrack(t)->GetEnergyDeposition(isig);
                            index = k; //Set index != -1 to skip over the addition of elements to the vector in the next part of the loop.
                        }
                    }

                    if(index == -1){ //If VolumeId is not in the vector, add it to the vector and add the energy deposition to the energy desposition vector
                        //cout << "VID not found, adding" << endl;
                        MCVolid.push_back(VolumeId);
                        MCSpec.push_back(MCEvent->GetTrack(t)->GetPdg()); //This is one sticky point, how to deal with two different particles having edeps in the same location, unsure right now.
                        //Actually maybe fine; if muon ejects an electron which deposits energy into the same strip, why bother separating those edeps when the instrument would see them as the same?
                        MCEdep.push_back(MCEvent->GetTrack(t)->GetEnergyDeposition(isig));
                    }

                } //End loop over the hits in an MC truth track
            } //End loop over tracks

            int TrueUMBflag = 0;
            int TrueCORflag = 0;
            int TrueCBEtopflag = 0;
            int TrueCBEbotflag = 0;
            int TrueCBEsideflag = 0;

            //cout << endl;

        /*
        for(int k = 0; k < MCVolid.size(); k++){
            if(MCEdep[k] > 0.4)cout << "MC hit: " << MCEdep[k] << " at " << MCVolid[k] << " by " << MCSpec[k] << endl;
        }*/

        for (int k = 0; k < MCVolid.size(); k++) { //Check to see if there's a significant hit from a primary particle (need to choose)
            if(MCEdep[k] > 0.4 && MCSpec[k] == pid){ MCST = 1;} //Checking to make sure all significant hits from the primary track
            if(MCEdep[k] > 0.4 && MCSpec[k] != pid){MCST = 0; break;} //cout << "NOT single track MC Event" << endl;
        }

        if(MCST){ //Checking for single tracks
            MCmust++; //Counter for MC true single tracks
            MCmust_Weighted = MCmust_Weighted + RateScale;
            //cout << "Event " << i << " MC Single track! " << endl;

            for (int k = 0; k < MCVolid.size(); k++) {
                if(MCEdep[k] > 0.4){
                    //if(volspec(MCVolid[k],0,3) == 111) cout << "MC hit: " << MCEdep[k] << " at " << MCVolid[k] << " by " << MCSpec[k] << endl;

                    if(volspec(MCVolid[k],0,3) == 100)TrueUMBflag++;
                    if(volspec(MCVolid[k],0,3) == 110){ TrueCBEtopflag++; /*cout << "VolumeId " << MCVolid[k] << " CBEtop hit! " << volspec(MCVolid[k],3,4) << endl; HCBEtop->Fill(volspec(MCVolid[k],3,4));*/ }
                    if(volspec(MCVolid[k],0,3) == 111){ TrueCBEbotflag++;  /*HCBEbot->Fill(volspec(MCVolid[k],3,4));*/ }
                    if(volspec(MCVolid[k],0,3) == 112 || volspec(MCVolid[k],0,3) == 113 || volspec(MCVolid[k],0,3) == 114 || volspec(MCVolid[k],0,3) == 115 || volspec(MCVolid[k],0,3) == 116)TrueCBEsideflag++;
                    if(volspec(MCVolid[k],0,3) == 102 || volspec(MCVolid[k],0,3) == 103 || volspec(MCVolid[k],0,3) == 104 || volspec(MCVolid[k],0,3) == 105 || volspec(MCVolid[k],0,3) == 106)TrueCORflag++;
                }
                //cout << "VID " << MCVolid[k] << " total edep " << MCEdep[k] << endl;
            } //This seems to be a fine way of determining total edeps in the instrument from MC

            //Flags checked!
                if(TrueCBEsideflag == 0 && TrueCORflag == 0 ){ //No side TOF hits
                    //cout << "No sides True Event: " << i << endl;
                    //cout << "Beta is " << Event->GetPrimaryBetaGenerated() << " Angle is " << Event->GetPrimaryMomentumDirectionGenerated().CosTheta() << " RateScale = " << RateScale << endl;
                    MCnosides++; //cout << "Event: " << i << " NO SIDES HIT!" << endl;
                    MCnosides_Weighted = MCnosides_Weighted + RateScale;
                    if( (TrueUMBflag>0) && (TrueCBEtopflag>0) && (TrueCBEbotflag>0) ) {
                        MCyesumb_ctop_cbot++;
                        MCyesumb_ctop_cbot_Weighted = MCyesumb_ctop_cbot_Weighted + RateScale;
                        //myfile << i << "\t" << "Truth: YesUMBCtopCBot" << endl;
                        //cout << "Event: " << i << " UMB TOP BOT YES!" << endl;
                    }
                    if( (TrueUMBflag>0) && (TrueCBEtopflag>0) && (TrueCBEbotflag==0) ){
                        MCyesumb_ctop_nocbot++;
                        MCyesumb_ctop_nocbot_Weighted = MCyesumb_ctop_nocbot_Weighted + RateScale;
                        //myfile << i << "\t" << "Truth: YesUMBCtop NoCBot" << endl;
                        //cout << "Event: " << i << " UMB TOP NO BOT!" << endl;
                    }
                    /*if( (TrueUMBflag>0) && (TrueCBEtopflag==0) && (TrueCBEbotflag>0) ) {
                        MCyesumb_cbot_nocbot++;
                        cout << "Event: " << i << " UMB BOT YES TOP NO!" << endl;
                    }*/
                }

        } // End MC single track trigger

        if(Event->GetNTracks() == 1){
            must++;
            must_Weighted = must_Weighted + RateScale;

            int UMBflag = 0;
            int CORflag = 0;
            int CBEtopflag = 0;
            int CBEbotflag = 0;
            int CBEsideflag = 0;

            //cout << endl;
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

                 if(CBEsideflag == 0 && CORflag == 0 ){ //No side TOF hits
                     nosides++; //cout << "Event: " << i << " NO SIDES HIT!" << endl;
                     nosides_Weighted = nosides_Weighted + RateScale;
                     if( (UMBflag>0) && (CBEtopflag>0) && (CBEbotflag>0) ) {
                         yesumb_ctop_cbot++;
                         yesumb_ctop_cbot_Weighted = yesumb_ctop_cbot_Weighted + RateScale;
                         //myfile << i << "\t" << "Rec: YesUMBCtopCBot" << endl;
                         //cout << "Event: " << i << " UMB TOP BOT YES!" << endl;
                     }
                     if( (UMBflag>0) && (CBEtopflag>0) && (CBEbotflag==0) ){
                         yesumb_ctop_nocbot++;
                         yesumb_ctop_nocbot_Weighted = yesumb_ctop_nocbot_Weighted + RateScale;
                         //myfile << i << "\t" << "Rec: YesUMBCtop NoCBot" << endl;
                         //cout << "Event: " << i << " UMB TOP NO BOT!" << endl;
                     }
                 } //No sides


			//-----------EVENT LEVEL CUTS END

            } //Closed bracket for reconstructed single track

		} //Closed bracket for event level cut

}  //Closed bracket for iteration through tree events, move on to the next event i

myfile << "Total number of events: " << (float)TreeRec->GetEntries() << endl;
myfile.open(out_path + "MCTruthCuts_V2.txt",std::ios::app);
myfile << "Total number of events: " << (float)TreeRec->GetEntries() << endl;

myfile << "MC Total Events in Beta Range " << MCpassB << endl;
myfile << fixed << setprecision(2) << 100*(float)MCpassB/(float)TreeRec->GetEntries() <<"%" <<  endl;
myfile << "MC Single Track Events " << MCmust << endl;
myfile << fixed << setprecision(2) << 100*(float)MCmust/(float)MCpassB << "%" << endl;
myfile << "MC Single Track Satisfied, No sides " << MCnosides << endl;
myfile << fixed << setprecision(2) << 100*(float)MCnosides/(float)MCmust << "%" << endl;
myfile << "MC Single Track Satisfied, No sides, UMB, CT, CB " << MCyesumb_ctop_cbot << endl;
myfile << fixed << setprecision(2) << 100*(float)MCyesumb_ctop_cbot/(float)MCnosides <<"%" <<  endl;
myfile << "MC Single Track Satisfied, No sides or CB, UMB, CT " << MCyesumb_ctop_nocbot << endl;
myfile << fixed << setprecision(2) << 100*(float)MCyesumb_ctop_nocbot/(float)MCnosides <<"%" <<  endl;

myfile << "MC Total Events in Beta Range " << MCpassB << endl;
myfile << fixed << setprecision(2) << 100*(float)MCpassB/(float)TreeRec->GetEntries() <<"%" <<  endl;
myfile << "Reco Single Track Events " << must << endl;
myfile << fixed << setprecision(2) << 100*(float)must/(float)MCpassB << "%" << endl;
myfile << "Reco Single Track Satisfied, No sides " << nosides << endl;
myfile << fixed << setprecision(2) << 100*(float)nosides/(float)must << "%" << endl;
myfile << "Reco Single Track Satisfied, No sides, UMB, CT, CB " << yesumb_ctop_cbot << endl;
myfile << fixed << setprecision(2) << 100*(float)yesumb_ctop_cbot/(float)nosides <<"%" <<  endl;
myfile << "Reco Single Track Satisfied, No sides or CB, UMB, CT " << yesumb_ctop_nocbot << endl;
myfile << fixed << setprecision(2) << 100*(float)yesumb_ctop_nocbot/(float)nosides <<"%" <<  endl;

myfile << "Total number of events Weighted: " << Total_Weighted << endl;
myfile << "MC Weighted Total Events in Beta Range " << MCpassB_Weighted << endl;
myfile << fixed << setprecision(2) << 100*(float)MCpassB_Weighted/(float)Total_Weighted << "%" << endl;
myfile << "MC Weighted Single Track Events " << MCmust_Weighted << endl;
myfile << fixed << setprecision(2) << 100*(float)MCmust_Weighted/(float)MCpassB_Weighted << "%" << endl;
myfile << "MC Weighted Single Track Satisfied, No sides " << MCnosides_Weighted << endl;
myfile << fixed << setprecision(2) << 100*(float)MCnosides_Weighted/(float)MCmust_Weighted << "%" << endl;
myfile << "MC Weighted Single Track Satisfied, No sides, UMB, CT, CB " << MCyesumb_ctop_cbot_Weighted << endl;
myfile << fixed << setprecision(2) << 100*(float)MCyesumb_ctop_cbot_Weighted/(float)MCnosides_Weighted <<"%" <<  endl;
myfile << "MC Weighted Single Track Satisfied, No sides or CB, UMB, CT " << MCyesumb_ctop_nocbot_Weighted << endl;
myfile << fixed << setprecision(2) << 100*(float)MCyesumb_ctop_nocbot_Weighted/(float)MCnosides_Weighted <<"%" <<  endl;

myfile << "MC Weighted Total Events in Beta Range " << MCpassB_Weighted << endl;
myfile << "Reco Weighted Single Track Events " << must_Weighted << endl;
myfile << fixed << setprecision(2) << 100*(float)must_Weighted/(float)MCpassB_Weighted << "%" << endl;
myfile << "Reco Weighted Single Track Satisfied, No sides " << nosides_Weighted << endl;
myfile << fixed << setprecision(2) << 100*(float)nosides_Weighted/(float)must_Weighted << "%" << endl;
myfile << "Reco Weighted Single Track Satisfied, No sides, UMB, CT, CB " << yesumb_ctop_cbot_Weighted << endl;
myfile << fixed << setprecision(2) << 100*(float)yesumb_ctop_cbot_Weighted/(float)nosides_Weighted <<"%" <<  endl;
myfile << "Reco Weighted Single Track Satisfied, No sides or CB, UMB, CT " << yesumb_ctop_nocbot_Weighted << endl;
myfile << fixed << setprecision(2) << 100*(float)yesumb_ctop_nocbot_Weighted/(float)nosides_Weighted <<"%" <<  endl;

myfile.close();

//histplot1f("c1", HCBEtop, "CBE top Paddle Occupancy Plot","Paddle Number","NEvents", out_path + "CBE1dOccu");
//histplot1f("c2", HCBEbot, "CBE bot Paddle Occupancy Plot","Paddle Number","NEvents", out_path + "CBE1dOccu");

myfile.close();

cout << endl << "I am done" << endl;

return 1;

}
