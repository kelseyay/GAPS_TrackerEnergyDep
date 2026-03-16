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
myfile.open(out_path + "MCCharge.txt");
myfile << TString::Format( "Filename : %s", reco_path.c_str() )  << endl;
myfile << TString::Format( "Beta High : %f", betahigh) << endl;
myfile << TString::Format( "Beta Low : %f", betacut) << endl;
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
double fitlow = 0.5;
double fithigh = 3.5;
double acutlow = 1.5;
const Int_t NBins = 50;

//Full tracker histogram range
double mpvmin = 0.66;
double mpvmax = 0.75;

//From fitted peaks from flight! Check over!
double TofAngleCorrectedMip = 0.82;
double TrackerAngleCorrectedMip = 0.57;

int passctr = 0;
int strkctr = 0;

int alphactr = 0;
int pctr = 0;
int muctr = 0;
int yesid_actr = 0;
int yesid_pctr = 0;

float Ztof = 5.574;
float Atof = 10.3;
float Ltof = 0.635;
float rtof = 1.032;
float Gmctof = (1/1.1646);
float Gtofnew = (1/0.71);

float Ztkr = 14;
float Atkr = 28.3;
float Ltkr = 0.22;
float rtkr = 2.33;
float Gmctkr = (1/0.928);
float Gtkrnew = (1/0.83);
//Prepare cuts:
map<int, unsigned int> TofIndexVolumeIdMap;

//For Plotting purposes
ca::GPlottingTools Plotting;
char text[400]; //This variable is used later to name the plots

TH1D * HTruncatedMeanZevent = Plotting.DefineTH1D("HTruncatedMeanZevent",100, 0, 10, "sqrt(sqrt(truncated mean E))nergy deposition downgoing MIP [MIP]", "entries", 0.5, 1e4);
TH1D * HChargeMip = Plotting.DefineTH1D("HChargeMip",200, 0, 6, "particle charge for downgoing MIP", "entries", 0.5, 1e4);

//TH2D * HGenB_vs_GenZ = new TH2D("HGenB_vs_GenZ","Gen_Beta * Gen_Z vs Gen_Beta",50,betacut - 0.1, betahigh + 0.1, 50, 0.5*betacut -0.1 , 1.5*betahigh*2 + 0.1 );
TH2D * HRecB_vs_RecBTrunM = new TH2D("HRecB_vs_RecBTrunM","Rec_Beta * Tr_Mean vs Rec_Beta",50,betacut - 0.1, betahigh + 0.1, 50, 0.5*betacut -0.1 , 1.5*betahigh*2 + 0.1 );
TH2D * HGenB_vs_GenBTrunM = new TH2D("HGenB_vs_GenB2TrunM2","Gen_Beta * Tr_Mean vs Gen_Beta",50, betacut - 0.1, betahigh + 0.1, 50, 0.5*betacut -0.1 , 1.5*betahigh*2 + 0.1 ) ;
TH2D * HRecB_vs_GenZ = new TH2D("HRecB_vs_GenZ","Gen_Z vs Rec_Beta",50, betacut - 0.1, betahigh + 0.1,50,  0.5*betacut -0.1 , 1.5*betahigh*2 + 0.1);

TH2D * HRecB_vs_GenB= new TH2D("HRecB_vs_GenB","Rec_B vs Gen_Beta",50, betacut - 0.1, betahigh + 0.1,50, betacut - 0.1,betahigh + 0.1);
TH2D * HMisIDRecB_vs_GenB= new TH2D("HMisIDRecB_vs_GenB","Rec_B vs Gen_Beta",50, betacut - 0.1, betahigh + 0.1,50, betacut - 0.5,betahigh + 0.5);
TH2D * HTrunM_vs_GenB= new TH2D("HTrunM_vs_GenB","Tr_Mean vs Gen_Beta",50, betacut - 0.1, betahigh + 0.1,50,  0.5 , 3.5);


//TH1D * hedep;
//hedep = new TH1D ("h0", ("Edep l Beta " + to_string(betacut) + " - " + to_string(betahigh) ).c_str(), NBins, xlow,xhigh);


//How many entries:
cout << "Total Number of events / Mainscale Factor = " << TreeRec->GetEntries()/MainLoopScaleFactor << endl;
cout << "Beta high " << betahigh << endl;
cout << "Beta cut " << betacut << endl;
cout << "MPV range on 2D Full Tracker will be " << mpvmin << " - " << mpvmax << endl;

//Now we can go over the loop
TreeRec->GetEntry(0);
TreeMC->GetEntry(0);

cout << "Total Number of events / Mainscale Factor = " << TreeRec->GetEntries()/MainLoopScaleFactor << endl;

//Using i to loop over every event in the tree
//for(unsigned int i = 0; i < 10; i+=MainLoopScaleFactor){
//for(unsigned int i = 0; i < TreeRec->GetEntries(); i+=MainLoopScaleFactor){
for(unsigned int i = 0; i < TreeRec->GetEntries()/MainLoopScaleFactor; i++){ //This is not the "correct" way to do this, but it's probably fine. Should be skipping M each time, but that seems to be really slow!!
    TreeRec->GetEntry(i);
    TreeMC->GetEntry(i);

	//Cuts are implemented in this chunk:
	if(Event->GetNTracks() == 1){  //First select the single track event
	    strkctr++;
		bool Umbflag = 0;
		bool CBEtopflag = 0;
		bool CBEbotflag = 0;

		if( ((int)i % (int)ceil(TreeRec->GetEntries()/(MainLoopScaleFactor*10))) == 0){
		    cout << "Event number " << i << endl;
		}

		//cout << "Even is " << i << endl;

		CTrackRec* pt = Event->GetPrimaryTrack();
		uint pt_index = 0;
      	        for( ; pt_index < Event->GetNTracks(); pt_index++) if( Event->GetTrack(pt_index)->IsPrimary() ) break;

		//Note downwards beta enforced by Event->GetPrimaryBeta() (should be positive) multiplied by Event->GetPrimaryMomentumDirection()[2] (z trajectory of particle)
		if((pt != nullptr && pt->GetChi2()/pt->GetNdof()) < 3.2 && -fabs(MCEvent->GetPrimaryMomentumDirection().CosTheta()) > -coslow && -fabs(MCEvent->GetPrimaryMomentumDirection().CosTheta()) < -coshigh && MCEvent->GetPrimaryBeta()*MCEvent->GetPrimaryMomentumDirection()[2] < 0 && fabs(MCEvent->GetPrimaryBeta()) >  betacut && fabs(MCEvent->GetPrimaryBeta()) <  betahigh ){
		//if(pt != nullptr && -fabs(Event->GetPrimaryMomentumDirection().CosTheta()) > -coslow && -fabs(Event->GetPrimaryMomentumDirection().CosTheta()) < -coshigh && Event->GetPrimaryBeta()*MCEvent->GetPrimaryMomentumDirection()[2] < 0 && fabs(Event->GetPrimaryBeta()) >  betacut && fabs(Event->GetPrimaryBeta()) <  betahigh ){
		//if(pt != nullptr && -fabs(Event->GetPrimaryMomentumDirectionGenerated().CosTheta()) > -coslow && -fabs(Event->GetPrimaryMomentumDirectionGenerated().CosTheta()) < -coshigh && Event->GetPrimaryBetaGenerated()*Event->GetPrimaryMomentumDirectionGenerated()[2] < 0 && fabs(Event->GetPrimaryBetaGenerated()) >  betacut && fabs(Event->GetPrimaryBetaGenerated()) <  betahigh ){

		    //-----------EVENT LEVEL CUT APPLIED

			//First iteration over event for flags and Zevent vector filling
			vector<double> Zevent;
			for(uint isig=0; isig<Event->GetTrack(0)->GetEnergyDeposition().size(); isig++){
                unsigned int VolumeId  = Event->GetTrack(0)->GetVolumeId(isig); //Check the VolumeId of the event
                double Bgen = MCEvent->GetPrimaryBeta();
                double Brec = Event->GetPrimaryBeta();


                if(GGeometryObject::IsTofVolume(VolumeId) && Event->GetTrack(0)->GetEnergyDeposition(isig) > TofCutLow){
                    //A hit in the COR or CBE_sides needs to be multiplied by sin(theta) instead of cos(theta)
                    if(volspec(VolumeId,2,1) == 0 || volspec(VolumeId,2,1) == 1){
                        //cout << Event->GetTrack(0)->GetEnergyDeposition(isig)*fabs(Event->GetPrimaryMomentumDirection().CosTheta())/TofAngleCorrectedMip << endl;
                        Zevent.push_back( Gtofnew*sqrt( (Event->GetTrack(0)->GetEnergyDeposition(isig)*fabs(Event->GetPrimaryMomentumDirection().CosTheta()))*Gmctof*pow(Brec,2)*(Atof/Ztof)*(1/(rtof*Ltof*0.307)) / (log( 1.022*pow(Brec,2) / ((1 - pow(Brec,2))*0.000016*pow(Ztof,0.9) ) ) - pow(Brec,2)) ) );
                        //cout << "FLAT PADDLE HIT" << endl;
                    } else{
                        //cout << Event->GetTrack(0)->GetEnergyDeposition(isig)*sqrt(1-pow(Event->GetPrimaryMomentumDirection().CosTheta(),2))/TofAngleCorrectedMip << endl;
                        Zevent.push_back( Gtofnew*sqrt( (Event->GetTrack(0)->GetEnergyDeposition(isig)*fabs(sqrt(1-pow(Event->GetPrimaryMomentumDirection().CosTheta(),2)))) *Gmctof*pow(Brec,2)*(Atof/Ztof)*(1/(rtof*Ltof*0.307)) / (log( 1.022*pow(Brec,2) / ((1 - pow(Brec,2))*0.000016*pow(Ztof,0.9) ) ) - pow(Brec,2)) ));
                        //cout << "VERTICAL PADDLE HIT" << endl;
                    }
                }


                if(GGeometryObject::IsTrackerVolume(VolumeId) && Event->GetTrack(0)->GetEnergyDeposition(isig) > TrackerCut){
                    //cout << Event->GetTrack(0)->GetEnergyDeposition(isig)*fabs(Event->GetPrimaryMomentumDirection().CosTheta())/TrackerAngleCorrectedMip << endl;
                    Zevent.push_back(Gtkrnew* sqrt(  (Event->GetTrack(0)->GetEnergyDeposition(isig)*fabs(Event->GetPrimaryMomentumDirection().CosTheta()))*Gmctkr*pow(Brec,2)*(Atkr/Ztkr)*(1/(rtkr*Ltkr*0.307)) / (log( 1.022*pow(Brec,2) / ((1 - pow(Brec,2))*0.000016*pow(Ztkr,0.9) ) ) - pow(Brec,2)) ));
                    //cout << "TRACKER HIT" << endl;
                }

                if(volspec(VolumeId,0,3) == 100){ Umbflag = 1; } // cout << "UMB hit!" <<endl ;
                //if(volspec(VolumeId,0,2) == 10 && volspec(VolumeId,2,1) != 0 ){ } //cout << "COR hit! " << endl;
                if(volspec(VolumeId,0,3) == 110) { CBEtopflag = 1; }// cout << "CBE top hit!" << endl;
                if(volspec(VolumeId,0,3) == 111) { CBEbotflag = 1; }// cout << "CBE bot hit!" << endl;

			}

			//Do we want to only run this on certain tracks? Yeah probably. Can remove the TOF flags. Just run on whatever lol.
			if(/*Umbflag && CBEtopflag && CBEbotflag && */ 1){

				//Calculate the truncated mean:
				//Prepare calculation :o
				//Want to take lowest energy deposits to avoid Landau long tail!

				std::sort(Zevent.begin(), Zevent.begin()+Zevent.size());
				double TruncatedMeanZevent = 0;
				double CtrTruncatedMeanZevent = 0;

				if(Zevent.size() == 0){ TruncatedMeanZevent = 0;
				}else{

				    for(unsigned int isig = 0; isig < double(Zevent.size())/2; isig++){
						TruncatedMeanZevent += Zevent.at(isig);
						CtrTruncatedMeanZevent++;
					}

					if(CtrTruncatedMeanZevent == 0){ TruncatedMeanZevent = 0;
					}else{ TruncatedMeanZevent /= CtrTruncatedMeanZevent;}

				}

				if(Zevent.size() == 1) TruncatedMeanZevent = Zevent.at(0);

				HTruncatedMeanZevent->Fill(TruncatedMeanZevent);
				HChargeMip->Fill(sqrt(TruncatedMeanZevent));
				//cout << "sqrt(truncated mean E) = " << TruncatedMeanZevent << endl;
				HRecB_vs_RecBTrunM->Fill(Event->GetPrimaryBeta(), (Event->GetPrimaryBeta()) * sqrt(TruncatedMeanZevent) );
				HGenB_vs_GenBTrunM->Fill(MCEvent->GetPrimaryBeta(), (MCEvent->GetPrimaryBeta()) * sqrt(TruncatedMeanZevent) );
				HTrunM_vs_GenB->Fill(MCEvent->GetPrimaryBeta(),sqrt(TruncatedMeanZevent));
				//No need for another iteration over the events

				//cout << "Event number " << i << " passes the cuts!" << endl;
				//cout << "Particle species " << MCEvent->GetTrack(0)->GetPdg() << endl;
				HRecB_vs_GenB->Fill(MCEvent->GetPrimaryBeta(),Event->GetPrimaryBeta());
				if(MCEvent->GetPrimaryPdg() == 1000020040){
				    alphactr++;
					//HGenB_vs_GenZ->Fill( MCEvent->GetPrimaryBeta(), 2*(MCEvent->GetPrimaryBeta())  );
					HRecB_vs_GenZ->Fill( Event->GetPrimaryBeta(), 2*(Event->GetPrimaryBeta() ) );
					if(sqrt(TruncatedMeanZevent)*Event->GetPrimaryBeta() > acutlow){
					    yesid_actr++; }else{
	                    HMisIDRecB_vs_GenB->Fill(MCEvent->GetPrimaryBeta(),Event->GetPrimaryBeta());
						//cout << "Wrong! Calculated Alpha Charge " << sqrt(TruncatedMeanZevent)*Event->GetPrimaryBeta() << endl;
						//cout << "Rec_B = " << Event->GetPrimaryBeta() << endl;
						//cout << "Gen_B = " << MCEvent->GetPrimaryBeta() << endl;
					}
				}else{
				    //cout << "Non-alpha event?! " << MCEvent->GetTrack(0)->GetPdg() << endl;
					//cout << "Event # " << i << endl;
				}


				if(MCEvent->GetPrimaryPdg() == 2212){
				    pctr++;
					//HGenB_vs_GenZ->Fill( MCEvent->GetPrimaryBeta(), 1*(MCEvent->GetPrimaryBeta()) );
					if(sqrt(TruncatedMeanZevent)*Event->GetPrimaryBeta() > acutlow){
					}
					HRecB_vs_GenZ->Fill( Event->GetPrimaryBeta(), 1*(Event->GetPrimaryBeta()) );
					if(sqrt(TruncatedMeanZevent)*Event->GetPrimaryBeta() < acutlow){
					    yesid_pctr++; }else{
					    HMisIDRecB_vs_GenB->Fill( MCEvent->GetPrimaryBeta(), Event->GetPrimaryBeta() );
					    //cout << "Wrong! Calculated Proton Charge " << sqrt(TruncatedMeanZevent)*Event->GetPrimaryBeta() << endl;
						//cout << "Rec_B = " << Event->GetPrimaryBeta() << endl;
						//cout << "Gen_B = " << MCEvent->GetPrimaryBeta() << endl;
					}
				}else{
				    //cout << "Non-proton event?! " << MCEvent->GetTrack(0)->GetPdg() << endl;
					//cout << "Event # " << i << endl;
					//cout << "Non-proton event?! " << endl;
					//for(int k = 0; k < MCEvent->GetNTracks(); k++){
					//    cout << "Track " << k << " Particle " << MCEvent->GetTrack(k)->GetPdg() << endl;
					//}

				}

				if(MCEvent->GetTrack(0)->GetPdg() == 13){
				    muctr++;
					//HGenB_vs_GenZ->Fill( MCEvent->GetPrimaryBeta(), 1*(MCEvent->GetPrimaryBeta()) );
					HRecB_vs_GenZ->Fill( Event->GetPrimaryBeta(), 1*(Event->GetPrimaryBeta()) );
				}

				/*
				for(uint isig=0; isig<Event->GetTrack(0)->GetEnergyDeposition().size(); isig++){
					unsigned int VolumeId  = Event->GetTrack(0)->GetVolumeId(isig);
					if(GGeometryObject::IsTrackerVolume(VolumeId) && Event->GetTrack(0)->GetEnergyDeposition(isig) > TrackerCut){
							hedep->Fill(  (Event->GetTrack(0)->GetEnergyDeposition(isig)*fabs(Event->GetPrimaryMomentumDirection().CosTheta()))   );
					} //Closed bracket for Tracker volume and tracker cutoff
				} //Closed bracket for iteration over event with TOF cuts
				*/

			} //Closed bracket for if statement for cuts


			//-----------EVENT LEVEL CUTS END


		} //Closed bracket for event level cut

	} //Closed bracket for single track cut

}  //Closed bracket for iteration through tree events, move on to the next event i


myfile.open(out_path + "MCCharge.txt",std::ios::app);
myfile << "Total Events/Mainscale Factor " << TreeRec->GetEntries()/MainLoopScaleFactor << endl;
myfile << "Single Track Events " << strkctr << endl;
myfile << "Events that pass cuts " << passctr << endl;

myfile << "# Alphas pass cuts " << alphactr << endl;
myfile << "# Alphas correctly identified " << yesid_actr << endl;
myfile << "Proportion alphas correctly ID-ed " << setprecision(2) << (float)yesid_actr/(float)alphactr << endl;
myfile << "# protons pass cuts " << pctr << endl;
myfile << "# protons correctly identified " << yesid_pctr << endl;
myfile << "Proportion protons correctly ID-ed " << setprecision(2)  << (float)yesid_pctr/(float)pctr << endl;
myfile.close();


//Histogram section
//--------------------------------------

histplot1d("c1", HChargeMip, "MIP Charge for "+to_string(alphactr)+" alphas, "+to_string(pctr)+" protons, "+to_string(muctr)+" mu","Charge","NEvents", out_path + "Both");
//histplot2d("c2", HGenB_vs_GenZ, "Gen_Z versus Gen_B","Generated Beta","Generated Z","NEntries", out_path + "test2D");
histplot2d("c3", HRecB_vs_RecBTrunM, "Rec_B * sqrt(Tr_Mean) versus Rec_B","Reconstructed Beta","sqrt(truncated mean E) * Reconstructed B","NEntries", out_path + "BothRrec2D");
histplot2d("c4", HGenB_vs_GenBTrunM, "Gen_B * sqrt(Tr_Mean) versus Gen_B","Generated Beta","sqrt(truncated mean E) * Generated B","NEntries",out_path + "BothGenBTrunM2D");
//histplot2d("c5", HRecB_vs_GenZ, "Rec_B * Gen_Z versus Rec_B","Reconstructed Beta","Generated Z","NEntries", out_path + "testGenZ2D");
histplot2d("c6", HTrunM_vs_GenB, "sqrt(Tr_Mean) versus Gen_B","Generated Beta","sqrt(truncated mean E)","NEntries", out_path + "BothGenBTrunM");
histplot2d("c7",HRecB_vs_GenB,"Rec_B versus Gen_B","Generated Beta", "Reconstructed Beta","NEntries", out_path + "BothgenBRecB" );
//histplot2d("c7",HMisIDRecB_vs_GenB,"Rec_B versus Gen_B Mis-ID Particles","Generated Beta", "Reconstructed Beta","NEntries", out_path + "testMisIDgenBRecB" );


cout << endl << "I am done" << endl;

return 1;

}
