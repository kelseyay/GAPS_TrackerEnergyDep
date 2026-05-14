//How to use: ./MCCutEff -i /home/kelsey/simulations/simdat/mu/v.2.1.2/mu-_gaps_triggerlevel1_FTFP_BERT_1744342800_rec -o test -w 1
//Super not thrilled with how slow this is...
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
double TofCutLow = 0.4;
//double TrackerCut = 0.3; //Threshold for an energy deposition to be considered a hit

//Prepare textile for saving values
std::ofstream myfile;
myfile.open(out_path + "MCTruthCuts.txt");
myfile << TString::Format( "Filename : %s", reco_path.c_str() )  << endl;
myfile << TString::Format( "Beta High : %f", betahigh) << endl;
myfile << TString::Format( "Beta Low : %f", betacut) << endl << endl;
myfile.close();

const int NCuts = 7; //Start simple!
string MCTrue_cutnames[NCuts+1] = {
    " Total NEvents ",
    " Total Events in Generated Beta Range ",
    " Single Track Events in Beta Range ",
    " All above, No sides ",
    " No Sides, YES UMB ",
    " No Sides, YES UMB, YES CT ",
    " No Sides, YES UMB, CT, YES CB ",
    " No Sides, YES UMB, CT, NO CB ",
};

string Weighted_MCTrue_cutnames[NCuts+1];
string MCReco_cutnames[NCuts+1];
string Weighted_MCReco_cutnames[NCuts+1];

for(int k = 0; k < NCuts + 1; k++){
    MCReco_cutnames[k] = "MC Reco: " + MCTrue_cutnames[k];
    Weighted_MCTrue_cutnames[k] = "Weighted MC True" + MCTrue_cutnames[k];
    Weighted_MCReco_cutnames[k] = "Weighted MC Reco" + MCTrue_cutnames[k];
    MCTrue_cutnames[k] = "MC True: " + MCTrue_cutnames[k];
}

int MCTrue_cuts[NCuts+1] = {};
int MCReco_cuts[NCuts+1] = {};
float Weighted_MCTrue_cuts[NCuts+1] = {};
float Weighted_MCReco_cuts[NCuts+1] = {};

MCTrue_cuts[0] = TreeRec->GetEntries(); //Number of events (not a cut)
MCReco_cuts[0] = TreeRec->GetEntries();

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
for(unsigned int i = 0; i < TreeRec->GetEntries()/MainLoopScaleFactor; i++){ //This is not the "correct" way to do this, but it's probably fine. Should be skipping M each time, but that seems to be really slow!!
    TreeRec->GetEntry(i);
    TreeMC->GetEntry(i);

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
            AcceptanceScale = MainLoopScaleFactor*StartingPlaneAcceptance/(BinWidthFactor*HPrimaryBeta->GetBinContent(HPrimaryBeta->FindBin(MCEvent->GetPrimaryBeta())));

            if(HPrimaryBeta->GetBinContent(HPrimaryBeta->FindBin(MCEvent->GetPrimaryBeta())) == 0) AcceptanceScale = 0;
            //Find the bin associated with the generated beta in the vector, see if it doesn't exist?
         }else AcceptanceScale = 1; //Set it = to 1 if there's a nullptr? That's a surprise...
        //cout << "AcceptanceScale = " << AcceptanceScale << endl;

        int AngularRegion = -1;
        for(unsigned int a = 0; a < CosZenithCut.size(); a++) if(MCEvent->GetPrimaryMomentumDirection().CosTheta() < CosZenithCut.at(a).first && MCEvent->GetPrimaryMomentumDirection().CosTheta() > CosZenithCut.at(a).second) AngularRegion = a;
        //This is just checking which "bin" the generated cos(theta) is in
        if(AngularRegion < 0) continue; //Don't bother if angular region wasn't found
        RateScale = FluxScaleFactor*AcceptanceScale*GMuonTotalFluxUnscaled.at(AngularRegion)->Eval(MCEvent->GetPrimaryBeta());
        //cout << "Event is " << i << " Beta is " << MCEvent->GetPrimaryBeta() << " Cos(theta) is " << MCEvent->GetPrimaryMomentumDirection().CosTheta() << " RateScale?? = " << RateScale << endl;
    }
    //Total_Weighted = Total_Weighted + RateScale;
    //cout << "RateScale = " << RateScale << endl;
    Weighted_MCTrue_cuts[0] = Weighted_MCTrue_cuts[0] + RateScale;
    //cout << " Weighted_MCTrue_cuts[0] = " << Weighted_MCTrue_cuts[0] << endl;
    Weighted_MCReco_cuts[0] = Weighted_MCReco_cuts[0] + RateScale;


	if(fabs(MCEvent->GetPrimaryBeta()) >  betacut && fabs(MCEvent->GetPrimaryBeta()) <  betahigh ){
	    MCTrue_cuts[1]++;
		Weighted_MCTrue_cuts[1] = Weighted_MCTrue_cuts[1] + RateScale;
		MCReco_cuts[1]++; //Single Track
        Weighted_MCReco_cuts[1] = Weighted_MCReco_cuts[1] + RateScale;

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

        for (int k = 0; k < MCVolid.size(); k++) { //Check to see if there's a significant hit from a primary particle (need to choose)
            if(MCEdep[k] > 0.4 && MCSpec[k] == pid){ MCST = 1;} //Checking to make sure all significant hits from the primary track
            if(MCEdep[k] > 0.4 && MCSpec[k] != pid){MCST = 0; break;} //cout << "NOT single track MC Event" << endl;
        }
        if(!MCST){cout << "Event " << i << " MC Truth Non Single Track!" << endl;}

        if(MCST){ //Checking for single tracks
            MCTrue_cuts[2]++; //Single Tracks in beta range
            Weighted_MCTrue_cuts[2] = Weighted_MCTrue_cuts[2] + RateScale;
            cout << "Event " << i << " MC Truth Single track! " << endl;

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
                    MCTrue_cuts[3]++; //MC No Sides
                    Weighted_MCTrue_cuts[3] = Weighted_MCTrue_cuts[3] + RateScale;
                    if(TrueUMBflag > 0){
                        MCTrue_cuts[4]++; //MC No sides, YES UMB
                        Weighted_MCTrue_cuts[4] = Weighted_MCTrue_cuts[4] + RateScale;
                        if(TrueCBEtopflag > 0){
                            MCTrue_cuts[5]++; //No sides, YES UMB, YES CT
                            Weighted_MCTrue_cuts[5] = Weighted_MCTrue_cuts[5] + RateScale;
                            if(TrueCBEbotflag>0){
                                MCTrue_cuts[6]++; //No sides, YES UMB, YES CT, YES CB
                                Weighted_MCTrue_cuts[6] = Weighted_MCTrue_cuts[6] + RateScale;
                            }else{
                                MCTrue_cuts[7]++; //No sides, YES UMB, YES CT, NO CB
                                Weighted_MCTrue_cuts[7] = Weighted_MCTrue_cuts[7] + RateScale;
                            }
                        }
                    }
                }

        } // End MC single track trigger

       	CTrackRec* pt = Event->GetPrimaryTrack();
        uint pt_index = 0;
        for( ; pt_index < Event->GetNTracks(); pt_index++) if( Event->GetTrack(pt_index)->IsPrimary() ) break;

        if(pt != nullptr && Event->GetNTracks() != 1){cout << "Event " << i << " MC Reco Non Single Track!" << endl;}
        //Start cuts on reconstructed data.
        if(pt != nullptr && Event->GetNTracks() == 1){
            cout << "Event " << i << " MC Reco Single Track!" << endl;
            MCReco_cuts[2]++; //Single Track
            Weighted_MCReco_cuts[2] = Weighted_MCReco_cuts[2] + RateScale;

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
                 MCReco_cuts[3]++; //No sides
                 Weighted_MCReco_cuts[3] = Weighted_MCReco_cuts[3] + RateScale;
                 if(UMBflag > 0){
                     MCReco_cuts[4]++; //No sides, YES UMB
                     Weighted_MCReco_cuts[4] = Weighted_MCReco_cuts[4] + RateScale;
                     if(CBEtopflag > 0){
                         MCReco_cuts[5]++; //No sides YES UMB, YES CBEtop
                         Weighted_MCReco_cuts[5] = Weighted_MCReco_cuts[5] + RateScale;
                         if((CBEbotflag>0)){
                             MCReco_cuts[6]++; //No sides YES UMB, YES CBEtop
                             Weighted_MCReco_cuts[6] = Weighted_MCReco_cuts[6] + RateScale;
                         }else{
                             MCReco_cuts[7]++; //No sides YES UMB, YES CBEtop
                             Weighted_MCReco_cuts[7] = Weighted_MCReco_cuts[7] + RateScale;
                         }
                     }
                 }
             } //No sides

			//-----------EVENT LEVEL CUTS END

            } //Closed bracket for reconstructed single track

		} //Closed bracket for event level cut

}  //Closed bracket for iteration through tree events, move on to the next event i

//myfile << "Total number of events: " << (float)TreeRec->GetEntries() << endl;
myfile.open(out_path + "MCTruthCuts.txt",std::ios::app);
//myfile << "Total number of events: " << (float)TreeRec->GetEntries() << endl;

for(int k = 0; k < NCuts+1; k++){
    myfile << MCTrue_cutnames[k] << MCTrue_cuts[k] << endl;
    if(k > 0 && k!= NCuts) myfile << fixed << setprecision(2) << 100*(float)MCTrue_cuts[k]/(float)MCTrue_cuts[k-1] <<"%" <<  endl; //NEntries doesn't need to divide by anything
    if(k == NCuts) myfile << fixed << setprecision(2) << 100*(float)MCTrue_cuts[k]/(float)MCTrue_cuts[k-2] <<"%" <<  endl; //Final cut is a percentage of two above
}

myfile << endl;

for(int k = 0; k < NCuts+1; k++){
    myfile << Weighted_MCTrue_cutnames[k] << Weighted_MCTrue_cuts[k] << endl;
    if(k > 0 && k!= NCuts) myfile << fixed << setprecision(2) << 100*(float)Weighted_MCTrue_cuts[k]/(float)Weighted_MCTrue_cuts[k-1] <<"%" <<  endl; //NEntries doesn't need to divide by anything
    if(k == NCuts) myfile << fixed << setprecision(2) << 100*(float)Weighted_MCTrue_cuts[k]/(float)Weighted_MCTrue_cuts[k-2] <<"%" <<  endl; //Final cut is a percentage of two above
}

myfile << endl;

for(int k = 0; k < NCuts+1; k++){
    myfile << MCReco_cutnames[k] << MCReco_cuts[k] << endl;
    if(k > 0 && k!= NCuts) myfile << fixed << setprecision(2) << 100*(float)MCReco_cuts[k]/(float)MCReco_cuts[k-1] <<"%" <<  endl; //NEntries doesn't need to divide by anything
    if(k == NCuts) myfile << fixed << setprecision(2) << 100*(float)MCReco_cuts[k]/(float)MCReco_cuts[k-2] <<"%" <<  endl; //Final cut is a percentage of two above
}

myfile << endl;

for(int k = 0; k < NCuts+1; k++){
    myfile << Weighted_MCReco_cutnames[k] << Weighted_MCReco_cuts[k] << endl;
    if(k > 0 && k!= NCuts) myfile << fixed << setprecision(2) << 100*(float)Weighted_MCReco_cuts[k]/(float)Weighted_MCReco_cuts[k-1] <<"%" <<  endl; //NEntries doesn't need to divide by anything
    if(k == NCuts) myfile << fixed << setprecision(2) << 100*(float)Weighted_MCReco_cuts[k]/(float)Weighted_MCReco_cuts[k-2] <<"%" <<  endl; //Final cut is a percentage of two above
}

myfile << endl;


myfile.close();

//histplot1f("c1", HCBEtop, "CBE top Paddle Occupancy Plot","Paddle Number","NEvents", out_path + "CBE1dOccu");
//histplot1f("c2", HCBEbot, "CBE bot Paddle Occupancy Plot","Paddle Number","NEvents", out_path + "CBE1dOccu");

myfile.close();

cout << endl << "I am done" << endl;

return 1;

}
