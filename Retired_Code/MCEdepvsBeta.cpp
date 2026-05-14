//To use: ./MCEdepvB -i /home/kelsey/simulations/simdat/mu/v.3.0.0/triggerlevel1/mu-_gaps_triggerlevel1_FTFP_BERT_1757850432_rec -o test2 -l 0.6 -u 1 -b 20 -w 1

using namespace std;

#include "KYtools.C"

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
//To run, do ./OneHist -i /home/kelsey/simulations/simdat/mu/v.2.1.2/mu -o test

int main(int argc, char *argv[]){

GOptionParser* parser = GOptionParser::GetInstance();
parser->AddProgramDescription("MC truth energy depositions in beta bins");
parser->AddCommandLineOption<string>("in_path", "path to instrument data files", "./*", "i");
parser->AddCommandLineOption<string>("out_file", "name of output root file", "", "o");
parser->AddCommandLineOption<double>("beta_low", "low Beta Cut",0.2,"l");
parser->AddCommandLineOption<double>("beta_high", "upper Beta Cut",1,"u");
parser->AddCommandLineOption<int>("bbins", "number of beta bins",4,"b");
parser->AddCommandLineOption<int>("cbins", "number of cos(theta) bins",2,"t");
parser->AddCommandLineOption<bool>("MC_Weighting", "Monte Carlo Weighting on or off?",0,"w");
parser->ParseCommandLine(argc, argv);
parser->Parse();

string reco_path = parser->GetOption<string>("in_path");
string out_path = parser->GetOption<string>("out_file");
double betacut = parser->GetOption<double>("beta_low");
double betahigh = parser->GetOption<double>("beta_high");
int bbins = parser->GetOption<int>("bbins");
int cbins = parser->GetOption<int>("cbins");
bool MC_Weight = parser->GetOption<bool>("MC_Weighting");
cout << "MC_weight " << MC_Weight << endl;
if(betacut <= 0 || betacut >=1){ betacut = 0.2; cout << "Error with low beta choice. Setting Beta low to 0.2" << endl; }
cout << "beta cut = " << betacut << endl;

if(betahigh <= 0 || betahigh >=2 || betahigh < betacut){ betahigh = 1; cout << "Error with high/low beta choice! Setting Beta upper to 1" << endl; }
cout << "beta high = " << betahigh << endl;

double bwid = (betahigh-betacut)/bbins;

cout << reco_path << endl;
//sprintf(FilenameRoot,"%s/%s*.root",argv[1], argv[2]);
//cout << argv[1] << endl;

char FilenameRoot[400];
sprintf(FilenameRoot,"%s*.root",reco_path.c_str());
//sprintf(FilenameRoot,reco_path.c_str()); //210 simu data on my computer!
cout << FilenameRoot << endl;

if(out_path != "" && out_path[out_path.length()-1] != '/' ){ cout <<  "out path no slash!" << endl; out_path = out_path + '/'; }

//Prepare Rec Event
CEventRec* Event = new CEventRec(); //New reconstructed event
TChain * TreeRec = new TChain("TreeRec"); //New TreeRec Tchain object (this is new to me)
TreeRec->SetBranchAddress("Rec", &Event); //Set the branch address using Event (defined above)
TreeRec->Add(FilenameRoot);

//Prepare MC event
CEventMc* MCEvent = new CEventMc(); //New reconstructed event
TChain * TreeMC = new TChain("TreeMc"); //New TreeMC Tchain object (this is new to me)
TreeMC->SetBranchAddress("Mc", &MCEvent); //Set the branch address using Event (defined above)
TreeMC->Add(FilenameRoot);

TChain*TreeSimulationParameter = new TChain("SimulationParameterTree");
TreeSimulationParameter->Add(FilenameRoot);
GSimulationParameter * Parameter = new GSimulationParameter;
TreeSimulationParameter->SetBranchAddress("SimulationParameter", &Parameter);
TreeSimulationParameter->GetEntry(0);


int MainLoopScaleFactor = 1; //Set this number to scale the step size. Larger means runs faster and fewer events
double TrackerCut = 0.4; //Threshold for an energy deposition to be considered a hit
double TofCut = 0.1;

double xlow = 0.1; //Low range for histogram MeV
double xhigh = 10; //High range for histogram MeV

double coshigh = 1; //0.995; //0.92 //0.54 is the highest angle that can hit UMB, CBEtop, CBEbot
double coslow = 0; //0.62 //0.8
double cwid = (coshigh-coslow)/cbins;

double fitlow = 0.5;
double fithigh = 3.5;
const Int_t NBins = 50;

//Full tracker histogram range
double mpvmin = 0.66;
double mpvmax = 0.75;

//Prepare textile for saving values
std::ofstream myfile;
myfile.open(out_path + "Edep_vs_Beta.txt");
myfile << TString::Format( "Filename : %s", reco_path.c_str() )  << endl;
myfile << TString::Format( "Beta High : %f", betahigh) << endl;
myfile << TString::Format( "Beta Low : %f", betacut) << endl;
myfile << TString::Format( "Beta Bins : %d", bbins) << endl;
myfile << "Beta = [" ; //This bracket thing should be a function!! woo

TH1F * htkr[bbins];
TH1F * htof[bbins];
TH1F * htof_weight[bbins];
TH1F * htof_cos[cbins];
TH1F * htkr_cos[cbins];

//record_bins( myfile, "Beta", bbins, betacut, bwid);

for(int i = 0; i < bbins; i++){
    cout << "low " << betacut+bwid*i << endl;
    cout << "high "<< betacut+bwid*(i+1) << endl;
    myfile << betacut+bwid*(i+0.5);
    if(i != bbins - 1) myfile << ",";
    htkr[i] = new TH1F (("htkr"+to_string(i)).c_str(), ("Edep l Beta " + to_string(betacut+bwid*i) + " - " + to_string(betacut+bwid*(i+1)) ).c_str(), NBins, xlow,xhigh);
    htof[i] = new TH1F (("htof"+to_string(i)).c_str(), ("Edep l Beta " + to_string(betacut+bwid*i) + " - " + to_string(betacut+bwid*(i+1)) ).c_str(), NBins, xlow,xhigh);
}

myfile << "]" << endl;

myfile << "Cos(theta) = [" ;

for(int i = 0; i < cbins; i++){
    myfile << coslow+cwid*(i+0.5);
    if(i != cbins - 1) myfile << ",";
    htof_cos[i] = new TH1F (("htof_cos"+to_string(i)).c_str(), ("MIP Edep l Cos " + to_string(coslow+cwid*i) + " - " + to_string(coslow+cwid*(i+1)) ).c_str(), NBins, xlow,xhigh);
    htkr_cos[i] = new TH1F (("htkr_cos"+to_string(i)).c_str(), ("MIP Edep l Cos " + to_string(coslow+cwid*i) + " - " + to_string(coslow+cwid*(i+1)) ).c_str(), NBins, xlow,xhigh);
}

myfile << "]" << endl;

//TH2F(name, title, nbinsx, xlow, xup, nbinsy, ylow, yup)
auto h2dbetaTKR = new TH2F("h2dbetaTKR","Energy x Cos(theta) Distribution vs Beta" ,bbins,betacut,betahigh,NBins, xlow,10);
auto h2dbetaTOF = new TH2F("h2dbetaTOF","Energy x Cos(theta) Distribution vs Beta" ,bbins,betacut,betahigh,NBins, xlow,10);
auto h2dbetaTOF_Weight = new TH2F("h2dbetaTOF_Weight","Energy x Cos(theta) Distribution vs Beta" ,bbins,betacut,betahigh,NBins, xlow,10);
auto h2dbetaTKR_Weight = new TH2F("h2dbetaTOF_Weight","Energy x Cos(theta) Distribution vs Beta" ,bbins,betacut,betahigh,NBins, xlow,10);
auto h2dcosTOF = new TH2F("h2dcosTOF","Energy x Cos(theta) Distribution vs Cos(Theta)" ,cbins,coslow,coshigh,NBins, xlow,10);
auto h2dcosTKR = new TH2F("h2dcosTKR","Energy x Cos(theta) Distribution vs Cos(Theta)" ,cbins,coslow,coshigh,NBins, xlow,10);
auto h2dcosTOF_Weight = new TH2F("h2dcosTOF_Weight","Energy x Cos(theta) Distribution vs Cos(Theta)" ,cbins,coslow,coshigh,NBins, xlow,10);
auto h2dcosTKR_Weight = new TH2F("h2dcos_TKR_Weight","Energy x Cos(theta) Distribution vs Cos(Theta)" ,cbins,coslow,coshigh,NBins, xlow,10);

//Prepare cuts:
map<int, unsigned int> TofIndexVolumeIdMap;

//For Plotting purposes
ca::GPlottingTools Plotting;
char text[400]; //This variable is used later to name the plots

//How many entries:
cout << "Total Number of events / Mainscale Factor = " << TreeRec->GetEntries()/MainLoopScaleFactor << endl;
cout << "Beta high " << betahigh << endl;
cout << "Beta cut " << betacut << endl;
cout << "MPV range on 2D Full Tracker will be " << mpvmin << " - " << mpvmax << endl;

bool print = 0;

//Now we can go over the loop
TreeRec->GetEntry(0);
TreeMC->GetEntry(0);

cout << "Total Number of events / Mainscale Factor = " << TreeRec->GetEntries()/MainLoopScaleFactor << endl;

//MC Weighting things below:

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

//End MC Weighting things! //Stopping point here, have not yet implemented the weighting


//Using i to loop over every event in the tree
//for(unsigned int i = 80; i < 100; i+=MainLoopScaleFactor){
for(unsigned int i = 0; i < TreeRec->GetEntries(); i+=MainLoopScaleFactor){
    TreeRec->GetEntry(i);
    TreeMC->GetEntry(i);

    if( ((int)i % (int)ceil(TreeRec->GetEntries()/10)) == 0){
	    cout << "Event number " << i << endl;
	}

    //cout << "Event is " << i << endl;
    //cout << "Beta is " << MCEvent->GetPrimaryBeta() << " Cos Theta " << MCEvent->GetPrimaryMomentumDirection().CosTheta() << endl;


	CTrackRec* pt = Event->GetPrimaryTrack();
	uint pt_index = 0;
      	        for( ; pt_index < Event->GetNTracks(); pt_index++) if( Event->GetTrack(pt_index)->IsPrimary() ) break;
    //Still don't want to do anything with weirt nullptr tracks in the reconstruction

	if(pt != nullptr && -fabs(MCEvent->GetPrimaryMomentumDirection().CosTheta()) < -fabs(coslow) && -fabs(MCEvent->GetPrimaryMomentumDirection().CosTheta()) > -fabs(coshigh) && MCEvent->GetPrimaryBeta()*MCEvent->GetPrimaryMomentumDirection()[2] < 0 && fabs(MCEvent->GetPrimaryBeta()) >  betacut && fabs(MCEvent->GetPrimaryBeta()) <  betahigh ){
	    //-----------EVENT LEVEL CUTS START

		int pid = MCEvent->GetTrack(0)->GetPdg(); //Identity of the primary track should be this
	    vector<int> MCSpec;
		vector<unsigned int> MCVolid;
		vector<double> MCEdep;
		bool MCST = 0; //MC Single track flag
		double MaxNPEdep = 0; //Maximum energy deposition not from the primary species

        for(uint t = 0; t < MCEvent->GetNTracks();t++){ //Each particle has its own track
            if(print) cout << "Track is " << t << endl;
            for(uint isig=0; isig<MCEvent->GetTrack(t)->GetEnergyDeposition().size(); isig++){ //Iterate over the energy depositions of the track
                unsigned int VolumeId  = MCEvent->GetTrack(t)->GetVolumeId(isig); //Check the volumeId of each of the hits in the iteration
                int index = -1; //Note the index = -1 right now
                if(print)  cout << "Hit " << isig << " Edep " << MCEvent->GetTrack(t)->GetEnergyDeposition(isig) << " VolumeId " << VolumeId << " PID " << MCEvent->GetTrack(t)->GetPdg() << endl;

                for (int k = 0; k < MCVolid.size(); k++) { //Iterate over the VolumeId vector
                    if (MCVolid[k] == VolumeId) { //If volumeid is in the vector already, add the energy deposition to the vector counting edeps there
                        MCEdep[k] = MCEdep[k] + MCEvent->GetTrack(t)->GetEnergyDeposition(isig);
                        index = k; //Set index != -1 to skip over the addition of elements to the vector in the next part of the loop.
                    }
                }

                if(index == -1){ //If VolumeId is not in the vector, add it to the vector and add the energy deposition to the energy desposition vector
                    //if(print)  cout << "VID not found, adding" << endl;
                    MCVolid.push_back(VolumeId);
                    MCSpec.push_back(MCEvent->GetTrack(t)->GetPdg()); //This is one sticky point, how to deal with two different particles having edeps in the same location, unsure right now.
                    //Actually maybe fine; if muon ejects an electron which deposits energy into the same strip, why bother separating those edeps when the instrument would see them as the same?
                    //See Event 5 of mu-_gaps_triggerlevel1_FTFP_BERT_1757857268_rec as an example. Muon deposits energy and electron deposits energy in the same strip, should be fine.
                    MCEdep.push_back(MCEvent->GetTrack(t)->GetEnergyDeposition(isig));
                }

            } //End loop over the hits in an MC truth track
        } //End loop over tracks

        for (int k = 0; k < MCVolid.size(); k++) { //Check to see if there's a significant hit from a primary particle (need to choose)
            if(print) cout << "Hit in event is " << k << " Volid " << MCVolid[k] << " Edep " << MCEdep[k] << " Species " << MCSpec[k] << endl;
            if(MCEdep[k] > 0.4 && MCSpec[k] == pid){ MCST = 1;} //Checking to make sure all significant hits from the primary track
            if(MCEdep[k] > 0.4 && MCSpec[k] != pid){MCST = 0; break;} //cout << "NOT single track MC Event" << endl;
        }

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
            //cout << "RateScale?? = " << RateScale << endl;
        }
        //cout << "RateScale = " << RateScale << endl;

        if(MCST){
            int TrueUMBflag = 0;
            int TrueCORflag = 0;
            int TrueCBEtopflag = 0;
            int TrueCBEbotflag = 0;
            int TrueCBEsideflag = 0;
            double beta = MCEvent->GetPrimaryBeta();
            double costheta = fabs(MCEvent->GetPrimaryMomentumDirection().CosTheta());
            int bbin = (int)floor((beta-betacut)/bwid);
            int cbin = (int)floor((costheta-coslow)/cwid);

            //First iteration over events to check for TOF hits
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
            } //End for loop for iterating over the hit series for TOF cuts

            if(TrueUMBflag > 0 && TrueUMBflag < 3 && TrueCBEtopflag > 0 && TrueCBEtopflag < 3 ){ //Require one or two hits in the TOF section
                if(print) cout << "Event is " << i << " Single track and TOF cut passed! " << endl;

                for (int k = 0; k < MCVolid.size(); k++) { //Iterating over the hit series of a single track event (MC)
                    if(MCEdep[k] > 0.4){
                        unsigned int VolumeId = MCVolid[k];

                        if(GGeometryObject::IsTofVolume(VolumeId) && MCEdep[k] > TofCut){
                            if(volspec(VolumeId,2,1) == 0 || volspec(VolumeId,2,1) == 1){
                                htof[bbin]->Fill( MCEdep[k]*fabs(costheta) );
                                h2dbetaTOF->Fill( beta ,MCEdep[k]*fabs(costheta) );
                                if( MC_Weight ) h2dbetaTOF_Weight->Fill( beta ,MCEdep[k]*fabs(costheta),RateScale );
                                if(beta > 0.8 && beta < 1.0){
                                    htof_cos[cbin]->Fill( MCEdep[k]*fabs(costheta) );
                                    h2dcosTOF->Fill( costheta, MCEdep[k]*fabs(costheta) );
                                    if( MC_Weight ) h2dcosTOF_Weight->Fill( costheta, MCEdep[k]*fabs(costheta),RateScale );
                                }

                                if(print) cout << "Flat paddle hit!! " << VolumeId << endl;
                            }else{
                                htof[bbin]->Fill( MCEdep[k]*sqrt(1-pow(costheta,2)) );
                                h2dbetaTOF->Fill( beta ,MCEdep[k]*sqrt(1-pow(costheta,2)) );
                                if( MC_Weight ) h2dbetaTOF_Weight->Fill( beta ,MCEdep[k]*sqrt(1-pow(costheta,2)),RateScale );
                                if(beta > 0.8 && beta < 1.0){
                                    htof_cos[cbin]->Fill( MCEdep[k]*sqrt(1-pow(costheta,2)) );
                                    h2dcosTOF->Fill( costheta, MCEdep[k]*sqrt(1-pow(costheta,2)) );
                                    if( MC_Weight ) h2dcosTOF_Weight->Fill( costheta, MCEdep[k]*sqrt(1-pow(costheta,2)), RateScale );
                                }
                                if(print) cout << "Vertical paddle hit!! " << endl;
                            }
                        }

                        if(GGeometryObject::IsTrackerVolume(VolumeId) && MCEdep[k] > TrackerCut){
                            htkr[bbin]->Fill( MCEdep[k]*fabs(costheta) );
                            h2dbetaTKR->Fill( beta ,MCEdep[k]*fabs(costheta) );
                            if( MC_Weight ) h2dbetaTKR_Weight->Fill( beta, MCEdep[k]*fabs(costheta),RateScale );
                            if(beta > 0.8 && beta < 1.0){
                                htkr_cos[cbin]->Fill( MCEdep[k]*fabs(costheta) );
                                h2dcosTKR->Fill( costheta, MCEdep[k]*fabs(costheta) );
                                if( MC_Weight ) h2dcosTKR_Weight->Fill( costheta, MCEdep[k]*fabs(costheta),RateScale );
                            }
                        }

                    }//Closed bracket for 0.4 MeV requirement
                }//Closed bracket requirement for iteration over hit series

            } //TOF cuts histogram filling done


        } //Closed bracket for single track MC truth requirement

		//-----------EVENT LEVEL CUTS END
	} //Closed bracket for event level cut

}  //Closed bracket for iteration through tree events, move on to the next event i



//Histogram section
//--------------------------------------

//Scan the 2D histogram for its max values in each beta bin.

TGraph *pts = new TGraph(); //Points to plot along with the histogram.
TGraph *errs = new TGraph(); //Points to plot along with the histogram.
TGraph *pts_TOF = new TGraph(); //Points to plot along with the histogram.
TGraph *errs_TOF = new TGraph(); //Points to plot along with the histogram.
TGraph *cos_pts_TOF = new TGraph(); //Points to plot along with the histogram.
TGraph *cos_errs_TOF = new TGraph(); //Points to plot along with the histogram.
TGraph *cos_pts_TKR = new TGraph(); //Points to plot along with the histogram.
TGraph *cos_errs_TKR = new TGraph(); //Points to plot along with the histogram.


//Good
//cout << "GetNBinsX histogram " <<  h2dbetaTKR->GetNbinsX() << " GetNBinsY histogram " << h2dbetaTKR->GetNbinsY() << endl;

label_2Dhisto(h2dbetaTKR, pts, errs);
label_2Dhisto(h2dbetaTOF, pts_TOF, errs_TOF);
label_2Dhisto(h2dcosTOF, cos_pts_TOF, cos_errs_TOF);
label_2Dhisto(h2dcosTKR, cos_pts_TKR, cos_errs_TKR);

histplot2f_pts_errs("c0", h2dbetaTKR,pts,errs, "MC True: TKR Energy x Cos(theta) Distribution vs Beta", "Beta" , "Energy x Cos(theta)", "NEntries", out_path + "TKR2D_vs_Beta" );
histplot2f_pts_errs("c1", h2dbetaTKR_Weight,pts,errs,"MC True: TKR Energy x Cos(theta) Distribution vs Beta", "Beta" , "Energy x Cos(theta)", "NEntries", out_path + "TKR2D_vs_Beta_Weight" );
histplot2f_pts_errs("c2", h2dbetaTOF,pts_TOF,errs_TOF,"MC True: TOF Energy x Sin/Cos(theta) Distribution vs Beta", "Beta" , "Angle Corrected Energy", "NEntries", out_path + "TOF2D_vs_Beta" );
histplot2f_pts_errs("c3", h2dbetaTOF_Weight,pts_TOF,errs_TOF,"MC True: TOF Energy x Sin/Cos(theta) Distribution vs Beta", "Beta" , "Angle Corrected Energy", "NEntries", out_path + "TOF2D_vs_Beta_Weight" );
histplot2f_pts_errs("c4", h2dcosTOF,cos_pts_TOF,cos_errs_TOF,"MC True: MIP TOF Energy x Sin/Cos(theta) Distribution vs Cos(Theta)", "Cos(Theta)" , "Angle Corrected Energy", "NEntries", out_path + "TOF2D_vs_Cos" );
histplot2f_pts_errs("c5", h2dcosTKR,cos_pts_TKR,cos_errs_TKR,"MC True: MIP TKR Energy x Cos(theta) Distribution vs Cos(Theta)", "Cos(Theta)" , "Energy x Cos(Theta)", "NEntries", out_path + "TKR2D_vs_Cos" );
histplot2f_pts_errs("c6", h2dcosTOF_Weight,cos_pts_TOF,cos_errs_TOF,"MC True: MIP TOF Energy x Sin/Cos(theta) Distribution vs Cos(Theta)", "Cos(Theta)" , "Angle Corrected Energy", "NEntries", out_path + "TOF2D_vs_Cos_Weight" );
histplot2f_pts_errs("c7", h2dcosTKR_Weight,cos_pts_TKR,cos_errs_TKR,"MC True: MIP TKR Energy x Cos(theta) Distribution vs Cos(Theta)", "Cos(Theta)" , "Energy x Cos(Theta)", "NEntries", out_path + "TKR2D_vs_Cos_Weight" );


myfile << "EdepTKR = [" ;

for(int i = 0; i < bbins; i++){
    //cout << "Mean of tkr hist " << i << " = " << htkr[i]->GetMean() << endl;
    myfile << htkr[i]->GetMean();
    if(i != bbins - 1) myfile << ",";
    histplot1f(("ctkr" + to_string(i)).c_str(), htkr[i], ( "Edep Beta " + to_string(betacut+bwid*i) + " - " + to_string(betacut+bwid*(i+1)) ).c_str(),"Energy Deposition Angle Corrected (Sin or Cos)","NEvents", out_path + ("Fulltkr" + to_string(i)).c_str() );
}

myfile << "]" << endl;

myfile << "EdepTOF = [" ;

for(int i = 0; i < bbins; i++){
    //cout << "Mean of tkr hist " << i << " = " << htkr[i]->GetMean() << endl;
    myfile << htof[i]->GetMean();
    if(i != bbins - 1) myfile << ",";
    histplot1f(("ctof" + to_string(i)).c_str(), htof[i], ( "Edep Beta " + to_string(betacut+bwid*i) + " - " + to_string(betacut+bwid*(i+1)) ).c_str(),"Energy Deposition Angle Corrected (Sin or Cos)","NEvents", out_path + ("Fulltof" + to_string(i)).c_str() );
}
myfile << "]" << endl;

for(int i = 0; i < cbins; i++){
    histplot1f(("ctof_cos" + to_string(i)).c_str(), htof_cos[i], ( "MIP Edep Cos " + to_string(coslow+cwid*i) + " - " + to_string(coslow+cwid*(i+1)) ).c_str(),"Cos(theta) Energy Deposition Angle Corrected (Sin or Cos)","NEvents", out_path + ("Fulltof_Cos" + to_string(i)).c_str() );
    histplot1f(("ctkr_cos" + to_string(i)).c_str(), htkr_cos[i], ( "MIP Edep Cos " + to_string(coslow+cwid*i) + " - " + to_string(coslow+cwid*(i+1)) ).c_str(),"Cos(theta) Energy Deposition Angle Corrected (Sin or Cos)","NEvents", out_path + ("Fulltkr_Cos" + to_string(i)).c_str() );
}

//Histogram for NEntries at a strip level

//--------------------------------------

myfile.close();

cout << endl << "I am done" << endl;

return 1;

}
