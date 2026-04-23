//To use: ./RecEdepvB -i /home/kelsey/simulations/simdat/mu/v.3.0.0/triggerlevel1/mu-_gaps_triggerlevel1_FTFP_BERT_1757850432_rec -o test2 -l 0.6 -u 1 -b 20

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
parser->AddCommandLineOption<int>("TRG", "Which trigger?",0,"r");
parser->AddCommandLineOption<bool>("MC_Weighting", "Monte Carlo Weighting on or off?",0,"w");
parser->ParseCommandLine(argc, argv);
parser->Parse();

string reco_path = parser->GetOption<string>("in_path");
string out_path = parser->GetOption<string>("out_file");
double betacut = parser->GetOption<double>("beta_low");
double betahigh = parser->GetOption<double>("beta_high");
int bbins = parser->GetOption<int>("bbins");
int cbins = parser->GetOption<int>("cbins");

int TRG = parser->GetOption<int>("TRG");
bool MC_Weight = parser->GetOption<bool>("MC_Weighting");

if(betacut <= 0 || betacut >=1){ betacut = 0.2; cout << "Error with low beta choice. Setting Beta low to 0.2" << endl; }
cout << "beta cut = " << betacut << endl;

if(betahigh <= 0 || betahigh >=2 || betahigh < betacut){ betahigh = 1; cout << "Error with high/low beta choice! Setting Beta upper to 1" << endl; }
cout << "beta high = " << betahigh << endl;

double bwid = (betahigh-betacut)/bbins;
cout << reco_path << endl;

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

int MainLoopScaleFactor = 1; //Set this number to scale the step size. Larger means runs faster and fewer events
double TrackerCut = 0.4; //Threshold for an energy deposition to be considered a hit
double TofCut = 0.1;

double xlow = 0.1; //Low range for histogram MeV
double xhigh = 2.5; //High range for histogram MeV

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
TH1F * htof_cos[cbins];
TH1F * htkr_cos[cbins];

//MC Weighting things below:

//This object is managing the plotting
ca::GPlottingTools Plotting;

double FluxScaleFactor = 71.1552/32.058750*0.25;

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

//End MC Weighting things! //Stopping point here, have not yet implemented the weighting


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
auto h2dcosTOF = new TH2F("h2dcosTOF","Energy x Cos(theta) Distribution vs Cos(Theta)" ,cbins,coslow,coshigh,NBins, xlow,10);
auto h2dcosTKR = new TH2F("h2dcosTKR","Energy x Cos(theta) Distribution vs Cos(Theta)" ,cbins,coslow,coshigh,NBins, xlow,10);

//Prepare cuts:
map<int, unsigned int> TofIndexVolumeIdMap;

char text[400]; //This variable is used later to name the plots

//How many entries:
cout << "Total Number of events / Mainscale Factor = " << TreeRec->GetEntries()/MainLoopScaleFactor << endl;
cout << "Beta high " << betahigh << endl;
cout << "Beta cut " << betacut << endl;
cout << "MPV range on 2D Full Tracker will be " << mpvmin << " - " << mpvmax << endl;

bool print = 0;

//Now we can go over the loop
TreeRec->GetEntry(0);

cout << "Total Number of events / Mainscale Factor = " << TreeRec->GetEntries()/MainLoopScaleFactor << endl;

//Using i to loop over every event in the tree
//for(unsigned int i = 80; i < 100; i+=MainLoopScaleFactor){
for(unsigned int i = 0; i < TreeRec->GetEntries(); i+=MainLoopScaleFactor){
    TreeRec->GetEntry(i);

    if( ((int)i % (int)ceil(TreeRec->GetEntries()/10)) == 0){
	    cout << "Event number " << i << endl;
	}

    //cout << "Event is " << i << endl;
    //cout << "Beta is " << MCEvent->GetPrimaryBeta() << " Cos Theta " << MCEvent->GetPrimaryMomentumDirection().CosTheta() << endl;


	CTrackRec* pt = Event->GetPrimaryTrack();
	uint pt_index = 0;
      	        for( ; pt_index < Event->GetNTracks(); pt_index++) if( Event->GetTrack(pt_index)->IsPrimary() ) break;
    //Still don't want to do anything with weirt nullptr tracks in the reconstruction

	if(pt != nullptr && ( (TRG == 0) || ((int)Event->GetTriggerSources().at(0) == TRG) ) && Event->GetPrimaryBeta()*Event->GetPrimaryMomentumDirection()[2] < 0 && -fabs(Event->GetPrimaryMomentumDirection().CosTheta()) < -fabs(coslow) && -fabs(Event->GetPrimaryMomentumDirection().CosTheta()) > -fabs(coshigh) && Event->GetPrimaryBeta()*Event->GetPrimaryMomentumDirection()[2] < 0 && fabs(Event->GetPrimaryBeta()) >  betacut && fabs(Event->GetPrimaryBeta()) <  betahigh ){

	//-----------EVENT LEVEL CUTS START

        //if(MCST == 1){cout << "Single Track Event!" << endl;}else{cout << "Not ST Event!" << endl;}

        if(Event->GetNTracks() == 1){  //First select the single track event
            int UMBflag = 0;
            int CORflag = 0;
            int CBEtopflag = 0;
            int CBEbotflag = 0;
            int CBEsideflag = 0;
            double beta = Event->GetPrimaryBeta();
            double costheta = fabs(Event->GetPrimaryMomentumDirection().CosTheta());
            int bbin = (int)floor((beta-betacut)/bwid);
            int cbin = (int)floor((costheta-coslow)/cwid);

            //First iteration over events to check for TOF hits
            for(uint isig=0; isig<Event->GetTrack(0)->GetEnergyDeposition().size(); isig++){
                unsigned int VolumeId = Event->GetTrack(0)->GetVolumeId(isig);
                if(volspec(VolumeId,0,3) == 100)UMBflag++;
                if(volspec(VolumeId,0,3) == 110){ CBEtopflag++; /*cout << "VolumeId " << VolumeId << " CBEtop hit! " << volspec(VolumeId,3,4) << endl; HCBEtop->Fill(volspec(VolumeId,3,4));*/ }
                if(volspec(VolumeId,0,3) == 111){ CBEbotflag++;  /*HCBEbot->Fill(volspec(VolumeId,3,4));*/ }
                if(volspec(VolumeId,0,3) == 112 || volspec(VolumeId,0,3) == 113 || volspec(VolumeId,0,3) == 114 || volspec(VolumeId,0,3) == 115 || volspec(VolumeId,0,3) == 116)CBEsideflag++;
                if(volspec(VolumeId,0,3) == 102 || volspec(VolumeId,0,3) == 103 || volspec(VolumeId,0,3) == 104 || volspec(VolumeId,0,3) == 105 || volspec(VolumeId,0,3) == 106)CORflag++;
            }

            if(UMBflag > 0 && UMBflag < 3 && CBEtopflag > 0 && CBEtopflag < 3 && (pt->GetChi2()/pt->GetNdof()) < 3.2 ){ //Require one or two hits in the TOF section
                if(print) cout << "Event is " << i << " Single track and TOF cut passed! " << endl;

                for(uint isig=0; isig<Event->GetTrack(0)->GetEnergyDeposition().size(); isig++){
                    unsigned int VolumeId = Event->GetTrack(0)->GetVolumeId(isig);

                    if(GGeometryObject::IsTofVolume(VolumeId) && Event->GetTrack(0)->GetEnergyDeposition(isig) > TofCut){
                        if(volspec(VolumeId,2,1) == 0 || volspec(VolumeId,2,1) == 1){
                            htof[bbin]->Fill(  Event->GetTrack(0)->GetEnergyDeposition(isig)*fabs(costheta) );
                            h2dbetaTOF->Fill( beta , Event->GetTrack(0)->GetEnergyDeposition(isig)*fabs(costheta) );
                            if(beta > 0.8 && beta < 1.0){
                                htof_cos[cbin]->Fill(  Event->GetTrack(0)->GetEnergyDeposition(isig)*fabs(costheta) );
                                h2dcosTOF->Fill( costheta,  Event->GetTrack(0)->GetEnergyDeposition(isig)*fabs(costheta) );
                            }

                            if(print) cout << "Flat paddle hit!! " << VolumeId << endl;
                        }else{
                            htof[bbin]->Fill( Event->GetTrack(0)->GetEnergyDeposition(isig)*sqrt(1-pow(costheta,2)) );
                            h2dbetaTOF->Fill( beta, Event->GetTrack(0)->GetEnergyDeposition(isig)*sqrt(1-pow(costheta,2)) );
                             if(beta > 0.8 && beta < 1.0){
                                htof_cos[cbin]->Fill( Event->GetTrack(0)->GetEnergyDeposition(isig)*sqrt(1-pow(costheta,2)) );
                                h2dcosTOF->Fill( costheta, Event->GetTrack(0)->GetEnergyDeposition(isig)*sqrt(1-pow(costheta,2)) );
                            }
                            if(print) cout << "Vertical paddle hit!! " << endl;
                        }
                    }

                    if(GGeometryObject::IsTrackerVolume(VolumeId) && Event->GetTrack(0)->GetEnergyDeposition(isig) > TrackerCut){
                        htkr[bbin]->Fill( Event->GetTrack(0)->GetEnergyDeposition(isig)*fabs(costheta) );
                        h2dbetaTKR->Fill( beta ,Event->GetTrack(0)->GetEnergyDeposition(isig)*fabs(costheta) );
                        if(beta > 0.8 && beta < 1.0){
                            htkr_cos[cbin]->Fill( Event->GetTrack(0)->GetEnergyDeposition(isig)*fabs(costheta) );
                            h2dcosTKR->Fill( costheta, Event->GetTrack(0)->GetEnergyDeposition(isig)*fabs(costheta) );
                        }
                    }

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

histplot2f_pts_errs("c1", h2dbetaTKR,pts,errs, "Rec: TKR Energy x Cos(theta) Distribution vs Beta", "Beta" , "Energy x Cos(theta)", "NEntries", out_path + "TKR2D_vs_Beta_Rec" );
histplot2f_pts_errs("c0", h2dbetaTOF,pts_TOF,errs_TOF,"Rec: TOF Energy x Sin/Cos(theta) Distribution vs Beta", "Beta" , "Angle Corrected Energy", "NEntries", out_path + "TOF2D_vs_Beta_Rec" );
histplot2f_pts_errs("c2", h2dcosTOF,cos_pts_TOF,cos_errs_TOF,"Rec: MIP TOF Energy x Sin/Cos(theta) Distribution vs Cos(Theta)", "Cos(Theta)" , "Angle Corrected Energy", "NEntries", out_path + "TOF2D_vs_Cos_Rec" );
histplot2f_pts_errs("c3", h2dcosTKR,cos_pts_TKR,cos_errs_TKR,"Rec: MIP TKR Energy x Sin/Cos(theta) Distribution vs Cos(Theta)", "Cos(Theta)" , "Energy x Cos(Theta)", "NEntries", out_path + "TKR2D_vs_Cos_Rec" );


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
