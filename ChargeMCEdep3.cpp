//FIX ME TO FLOOR THE ODD NHITS!!! I goofed ;__; Kind of a small error, but non negligible for sure

//How to use: ./DZedep3 -i /data1/nextcloud/cra_data/data/2025/production/v26.01/reconstructed/flight/251226/starlink251226_17 -s 0 -t 1 -f 1 -k 1 -r 2 -o test/
//Also use: ./DZedep3 -i /home/kelsey/simulations/simdat/flight/251226/26.01/starlink251226_1 -s 0 -t 1 -f 1 -k 1 -r 2 -o test/

//This NEW AND IMPROVED script will use the MPV method for calculating Z with the truncated NHits

using namespace std;

#include "KYtools.C"
#include <boost/math/special_functions/lambert_w.hpp>
//Headers included in KYtools are included here, so this is nice.

//FIXME: does this work on mac?
#include <sys/stat.h>

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
double TrackerCut = 0.4; //Threshold for an energy deposition to be considered a hit
int NHitsmin = 3;

double xlow = 0.1; //Low range for histogram MeV
double xhigh = 3; //High range for histogram MeV

double coshigh = 0.54; //0.995; //0.92 //0.54 is the highest angle that can hit UMB, CBEtop, CBEbot
double coslow = 1; //0.62 //0.8
double fitlow = 0.5;
double fithigh = 3.5;
double acutlow = 1.5;
const Int_t NBins = 50;

//Full tracker histogram range
double mpvmin = 0.66;
double mpvmax = 0.75;

bool TKR = parser->GetOption<bool>("TKR");
bool TF = parser->GetOption<bool>("TF");
int TRG = parser->GetOption<int>("TRG");

double tk = parser->GetOption<double>("tkr_factor");
double tf = parser->GetOption<double>("tof_factor");

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
float Gtofnew = tf; //1/0.9; //(1/0.71);
//cout << "sqrt(Gtofnew) = " << sqrt(tf) << endl;

float Ztkr = 14;
float Atkr = 28.3;
float Ltkr = 0.23; //0.249;
float rtkr = 2.33;
float Gtkrnew = tk; //1; // (1/0.83);
//cout << "sqrt(Gtkrnew) = " << sqrt(tk) << endl;

//Prepare cuts:
map<int, unsigned int> TofIndexVolumeIdMap;

//For Plotting purposes
ca::GPlottingTools Plotting;
char text[400]; //This variable is used later to name the plots

TH1D * HTruncatedMeanEnergyDepositionMip = Plotting.DefineTH1D("HTruncatedMeanEnergyDepositionMip",100, 0, 10, "sqrt(sqrt(truncated mean E))nergy deposition downgoing MIP [MIP]", "entries", 0.5, 1e4);
TH1D * HChargeMip = Plotting.DefineTH1D("HChargeMip",200, 0.1, 4, "particle charge for downgoing MIP", "entries", 0.5, 1e4);
TH1D * HTofMult = Plotting.DefineTH1D("HTofMult",200, 0.1, 2, "Multiplicative Factor TOF", "entries", 0.5, 1e4);
TH1D * HTkrMult = Plotting.DefineTH1D("HTkrMult",200, 0.1, 2, "Multiplicative Factor TKR", "entries", 0.5, 1e4);

//TH2D * HGenB_vs_GenZ = new TH2D("HGenB_vs_GenZ","Gen_Beta * Gen_Z vs Gen_Beta",50,betacut - 0.1, betahigh + 0.1, 50, 0.5*betacut -0.1 , 1.5*betahigh*2 + 0.1 );
TH2D * HRecB_vs_CalcZ = new TH2D("HRecB_vs_CalcZ","Rec_Beta * Tr_Mean vs Rec_Beta",50,betacut - 0.1, betahigh + 0.1, 50, 0.1 , 4);
//TH2D * HTrunM_vs_RecB= new TH2D("HTrunM_vs_RecB","Tr_Mean vs Rec_Beta",50, betacut - 0.1, betahigh + 0.1,50,  0.5 , 3.5);

//TH2D * HRecB_vs_Cal = new TH2D("HRecB_vs_RecBTrunM","Rec_Beta * Tr_Mean vs Rec_Beta",50,betacut - 0.1, betahigh + 0.1, 50, 0.5*betacut -0.1 , 1.5*betahigh*2 + 0.1 );

TH1D * hedep; //Full tracker energy deposition plot
hedep = new TH1D ("h0", ("Edep l Beta " + to_string(betacut) + " - " + to_string(betahigh) ).c_str(), NBins, xlow,xhigh);


//How many entries:
cout << "Total Number of events / Mainscale Factor = " << TreeRec->GetEntries()/MainLoopScaleFactor << endl;
cout << "Beta high " << betahigh << endl;
cout << "Beta cut " << betacut << endl;
cout << "MPV range on 2D Full Tracker will be " << mpvmin << " - " << mpvmax << endl;

//Now we can go over the loop
TreeRec->GetEntry(0);

cout << "Total Number of events / Mainscale Factor = " << TreeRec->GetEntries()/MainLoopScaleFactor << endl;

//Using i to loop over every event in the tree
//for(unsigned int i = 0; i < 1000; i+=MainLoopScaleFactor){
for(unsigned int i = 0; i < TreeRec->GetEntries(); i+=MainLoopScaleFactor){
//for(unsigned int i = 0; i < TreeRec->GetEntries()/MainLoopScaleFactor; i++){ //This is not the "correct" way to do this, but it's probably fine. Should be skipping M each time, but that seems to be really slow!!
    TreeRec->GetEntry(i);
    TreeMC->GetEntry(i);

    if( ((int)i % (int)ceil(TreeRec->GetEntries()/(MainLoopScaleFactor*10))) == 0){
		    cout << "Event number " << i << endl;
	}

	//Cuts are implemented in this chunk:
	if(Event->GetNTracks() == 1){  //First select the single track event
	    strkctr++;
		bool Umbflag = 0;
		bool CBEtopflag = 0;
		bool CBEbotflag = 0;

		CTrackRec* pt = Event->GetPrimaryTrack();
		uint pt_index = 0;
      	        for( ; pt_index < Event->GetNTracks(); pt_index++) if( Event->GetTrack(pt_index)->IsPrimary() ) break;


        //cout << "TRG = " << TRG << " and Event GTS = " << (int)Event->GetTriggerSources().at(0) << endl;

        //cout << ((int)Event->GetTriggerSources().at(0) == TRG) << endl;

		//Note downwards beta enforced by Event->GetPrimaryBeta() (should be positive) multiplied by Event->GetPrimaryMomentumDirection()[2] (z trajectory of particle)
		//if(pt != nullptr && (int)Event->GetTriggerSources().at(0) == TRG && (pt->GetChi2()/pt->GetNdof()) < 3.2 && -fabs(Event->GetPrimaryMomentumDirection().CosTheta()) > -coslow && -fabs(Event->GetPrimaryMomentumDirection().CosTheta()) < -coshigh && Event->GetPrimaryBeta()*Event->GetPrimaryMomentumDirection()[2] < 0 && fabs(Event->GetPrimaryBeta()) > betacut && fabs(Event->GetPrimaryBeta()) <  betahigh ){
		if( pt != nullptr && ( (TRG == 0) || ((int)Event->GetTriggerSources().at(0) == TRG) ) && pt->GetChi2()/pt->GetNdof() < 3.2 && -fabs(Event->GetPrimaryMomentumDirection().CosTheta()) > -coslow && -fabs(Event->GetPrimaryMomentumDirection().CosTheta()) < -coshigh && Event->GetPrimaryBeta()*Event->GetPrimaryMomentumDirection()[2] < 0 && fabs(MCEvent->GetPrimaryBeta()) > betacut && fabs(MCEvent->GetPrimaryBeta()) <  betahigh ){
		//if(pt != nullptr && -fabs(Event->GetPrimaryMomentumDirection().CosTheta()) > -coslow && -fabs(Event->GetPrimaryMomentumDirection().CosTheta()) < -coshigh && Event->GetPrimaryBeta()*MCEvent->GetPrimaryMomentumDirection()[2] < 0 && fabs(Event->GetPrimaryBeta()) >  betacut && fabs(Event->GetPrimaryBeta()) <  betahigh ){
		//if(pt != nullptr && -fabs(Event->GetPrimaryMomentumDirectionGenerated().CosTheta()) > -coslow && -fabs(Event->GetPrimaryMomentumDirectionGenerated().CosTheta()) < -coshigh && Event->GetPrimaryBetaGenerated()*Event->GetPrimaryMomentumDirectionGenerated()[2] < 0 && fabs(Event->GetPrimaryBetaGenerated()) >  betacut && fabs(Event->GetPrimaryBetaGenerated()) <  betahigh ){

			//-----------EVENT LEVEL CUT APPLIED

			//cout << "Event: " << i << endl;
			//First iteration over event for flags and EnergyDepositionMip vector filling
			vector<double> Zevent;
			double Bgen = MCEvent->GetPrimaryBeta();
			//double gamma = sqrt( 1/(1 - pow(Bgen,2)) );
			//cout << "Bgen " << Bgen << endl;

			for(uint isig=0; isig<Event->GetTrack(0)->GetEnergyDeposition().size(); isig++){
                unsigned int VolumeId  = Event->GetTrack(0)->GetVolumeId(isig); //Check the VolumeId of the event

                if(TF && GGeometryObject::IsTofVolume(VolumeId) && Event->GetTrack(0)->GetEnergyDeposition(isig) > TofCutLow){
                    //A hit in the COR or CBE_sides needs to be multiplied by sin(theta) instead of cos(theta)
                    double ion = 0.000016*pow(Ztof,0.9);
                    double eps_noz = 0.15355*(Ztof/Atof)*(rtof*Ltof)*(1/pow(Bgen,2));
                    double const_terms = 0.2 - pow(Bgen,2) + log( 1.022*pow(Bgen,2)/( ion*(1 - pow(Bgen,2))) )  ;

                    if(volspec(VolumeId,2,1) == 0 || volspec(VolumeId,2,1) == 1){
                        //This isn't actually the truncated mean. It's the energy deposition, the truncating happens in the Zevent phase
                        double Tr_Mean = Gtofnew*Event->GetTrack(0)->GetEnergyDeposition(isig)*fabs(Event->GetPrimaryMomentumDirection().CosTheta()) ;
                        Zevent.push_back( sqrt( (ion/eps_noz)*exp( boost::math::lambert_w0( (Tr_Mean/ion)*exp(const_terms) ) - const_terms  ) ) );

                        //cout << "z from MPV form Tr Mean = " << sqrt( (ion/eps_noz)*exp( boost::math::lambert_w0( (Tr_Mean/ion)*exp(const_terms) ) - const_terms  ) ) << endl;
                        //Zevent.push_back( sqrt((Gtofnew*Event->GetTrack(0)->GetEnergyDeposition(isig)*fabs(Event->GetPrimaryMomentumDirection().CosTheta()))*pow(Bgen,2)*(Atof/Ztof)*(1/(rtof*Ltof*0.307)) / (log( 1.022*pow(Bgen,2) / ((1 - pow(Bgen,2))*0.000016*pow(Ztof,0.9) ) ) - pow(Bgen,2)) ) );
                        //cout << "FLAT PADDLE HIT" << endl;
                    } else{
                        double Tr_Mean = Gtofnew*Event->GetTrack(0)->GetEnergyDeposition(isig)*fabs( sqrt(1-pow(Event->GetPrimaryMomentumDirection().CosTheta(),2)) ) ;
                        Zevent.push_back( sqrt( (ion/eps_noz)*exp( boost::math::lambert_w0( (Tr_Mean/ion)*exp(const_terms) ) - const_terms  ) ) );

                        //cout << Event->GetTrack(0)->GetEnergyDeposition(isig)*sqrt(1-pow(Event->GetPrimaryMomentumDirection().CosTheta(),2))/TofAngleCorrectedMip << endl;
                        //Zevent.push_back( sqrt( (Gtofnew*Event->GetTrack(0)->GetEnergyDeposition(isig)*fabs(sqrt(1-pow(Event->GetPrimaryMomentumDirection().CosTheta(),2)))) *pow(Bgen,2)*(Atof/Ztof)*(1/(rtof*Ltof*0.307)) / (log( 1.022*pow(Bgen,2) / ((1 - pow(Bgen,2))*0.000016*pow(Ztof,0.9) ) ) - pow(Bgen,2)) ));
                        //cout << "VERTICAL PADDLE HIT" << endl;
                    }
                }

                if(TKR && GGeometryObject::IsTrackerVolume(VolumeId) && Event->GetTrack(0)->GetEnergyDeposition(isig) > TrackerCut){
                    //A hit in the COR or CBE_sides needs to be multiplied by sin(theta) instead of cos(theta)
                    double ion = 0.000016*pow(Ztkr,0.9);
                    double eps_noz = 0.15355*(Ztkr/Atkr)*(rtkr*Ltkr)*(1/pow(Bgen,2));
                    double const_terms = 0.2 - pow(Bgen,2) + log( 1.022*pow(Bgen,2)/( ion*(1 - pow(Bgen,2))) )  ;
                    double Tr_Mean = Gtkrnew*Event->GetTrack(0)->GetEnergyDeposition(isig)*fabs(Event->GetPrimaryMomentumDirection().CosTheta()) ;
                    Zevent.push_back( sqrt( (ion/eps_noz)*exp( boost::math::lambert_w0( (Tr_Mean/ion)*exp(const_terms) ) - const_terms  ) ) );
                    //cout << "z from MPV form Tr Mean = " << sqrt( (ion/eps_noz)*exp( boost::math::lambert_w0( (Tr_Mean/ion)*exp(const_terms) ) - const_terms  ) ) << endl;

                    //cout << Event->GetTrack(0)->GetEnergyDeposition(isig)*fabs(Event->GetPrimaryMomentumDirection().CosTheta())/TrackerAngleCorrectedMip << endl;
                    //cout << "My factor TKR: " << sqrt( Gtkrnew*(Atkr/Ztkr)*(1/(rtkr*Ltkr*0.307)) / (log( 1.022*pow(Bgen,2) / ((1 - pow(Bgen,2))*0.000016*pow(Ztkr,0.9) ) ) - pow(Bgen,2)) ) << endl;
                    //cout << "TRACKER HIT" << endl;
                }

                if(volspec(VolumeId,0,3) == 100){ Umbflag = 1; } // cout << "UMB hit!" <<endl ;
                //if(volspec(VolumeId,0,2) == 10 && volspec(VolumeId,2,1) != 0 ){ } //cout << "COR hit! " << endl;
                if(volspec(VolumeId,0,3) == 110) { CBEtopflag = 1; }// cout << "CBE top hit!" << endl;
                if(volspec(VolumeId,0,3) == 111) { CBEbotflag = 1; }// cout << "CBE bot hit!" << endl;

			}

			//Do we want to only run this on certain tracks? Yeah probably. Can remove the TOF flags. Just run on whatever lol.
			if(Umbflag && CBEtopflag /*&& CBEbotflag && 1*/ ){
			    //cout << "Size: " << Zevent.size() << endl;

			    std::sort(Zevent.begin(), Zevent.begin()+Zevent.size());
				double TrZ = 0;
				double CtrTrZ = 0;

				if(Zevent.size() < NHitsmin){
				    TrZ = 0;
				}else{
			        for(unsigned int isig = 0; isig < double(Zevent.size())/2; isig++){
						TrZ += Zevent.at(isig);
						CtrTrZ++;
					}
					if(CtrTrZ == 0){ TrZ = 0;
					}else{ TrZ /= CtrTrZ;}
				}
				//cout << "Combined TrZ " << TrZ << endl;
				HChargeMip->Fill(TrZ);
				HRecB_vs_CalcZ->Fill(Event->GetPrimaryBeta(),TrZ);

				for(uint isig=0; isig<Event->GetTrack(0)->GetEnergyDeposition().size(); isig++){
				    unsigned int VolumeId  = Event->GetTrack(0)->GetVolumeId(isig);
					if(GGeometryObject::IsTrackerVolume(VolumeId) && Event->GetTrack(0)->GetEnergyDeposition(isig) > TrackerCut){
						hedep->Fill(  (Event->GetTrack(0)->GetEnergyDeposition(isig)*fabs(Event->GetPrimaryMomentumDirection().CosTheta()))   );
					} //Closed bracket for Tracker volume and tracker cutoff
				} //Closed bracket for iteration over event with TOF cuts

			    //cout << "Cuts passed " << endl;

			} //Closed bracket for if statement for TOF cuts

			//-----------EVENT LEVEL CUTS END


		} //Closed bracket for event level cut (beta, cos, pt exists)

	} //Closed bracket for single track cut

}  //Closed bracket for iteration through tree events, move on to the next event i


myfile.open(out_path + "MCCharge.txt",std::ios::app);
myfile << "Total Events/Mainscale Factor " << TreeRec->GetEntries()/MainLoopScaleFactor << endl;
myfile.close();


//Histogram section
//--------------------------------------

HRecB_vs_CalcZ->SaveAs( (out_path + "HRecB_vs_CalcZ.root").c_str() );
histplot1d("c1", HChargeMip, ("Charge Distribution Rec Beta " + to_string(betacut) + " - " + to_string(betahigh) ).c_str(),"Charge","NEvents", out_path + "Both");
histplot2d("c2", HRecB_vs_CalcZ, "Z_calc versus Gen_B","Generated Beta","Z_calc","NEntries", out_path + "BothBgen2D");

histplot1d("c5", hedep, "Tracker Energy Deposition for "+to_string(betacut)+" - "+to_string(betahigh),"Energy Deposit x Cos(theta)","NEvents", out_path + "Hedep");
cout << "Hedep Max Bin Center = " << hedep->GetBinCenter(hedep->GetMaximumBin()) << endl;


//histplot2d("c6", HTrunM_vs_RecB, "sqrt(Tr_Mean) versus Rec_B","Reconstructed Beta","sqrt(truncated mean E)","NEntries", out_path + "BothRecBTrunM");
//cout << "Max bin is " << HChargeMip->GetMaximumBin() << /*" factor should be " << 1/(HChargeMip->GetMaximumBin()) << */ endl;
//cout << "Max bin NEntries is " << HChargeMip->GetBinContent(HChargeMip->GetMaximumBin()) << endl;
cout << "TOF: " << TF << " TKR: " << TKR << endl;
cout << "Max bin Center is  " << HChargeMip->GetBinCenter(HChargeMip->GetMaximumBin()) << " factor should be " << pow(1/HChargeMip->GetBinCenter(HChargeMip->GetMaximumBin()),2) << endl;

cout << endl << "I am done" << endl;

return 1;

}
