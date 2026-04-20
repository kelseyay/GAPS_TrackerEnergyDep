//To run, do ./MCOneDet -i /home/kelsey/simulations/simdat/mu/v.3.0.0/triggerlevel1/mu-_gaps_triggerlevel1_FTFP_BERT_17
//HEY I'M MESSING AROUND HERE DON'T USE THIS CARELESSLY!
//I HAVE ADDED A FACTOR TO MC VERSION OF THIS CODE!

using namespace std;
#include "KYtools.C"

#include <TStyle.h>
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
#include "langaufun.C"

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


int main(int argc, char *argv[]){

GOptionParser* parser = GOptionParser::GetInstance();
parser->AddProgramDescription("Minimal Reproducable Example for Extracing Data from Reco Data");
parser->AddCommandLineOption<string>("in_path", "path to instrument data files", "./*", "i");
parser->AddCommandLineOption<string>("out_file", "name of output root file", "out.root", "o");
parser->ParseCommandLine(argc, argv);
parser->Parse();

string reco_path = parser->GetOption<string>("in_path");
string out_path = parser->GetOption<string>("out_file");
cout << reco_path << endl;
//sprintf(FilenameRoot,"%s/%s*.root",argv[1], argv[2]);
//cout << argv[1] << endl;

int MainLoopScaleFactor = 1;

char FilenameRoot[400];
sprintf(FilenameRoot,"%s*.root",reco_path.c_str());
//sprintf(FilenameRoot,reco_path.c_str()); //210 simu data on my computer!
cout << FilenameRoot << endl;

double TrackerCut = 0.3; //Threshold for an energy deposition to be considered a hit

double xlow = 0.3; //Low range for histogram MeV
double xhigh = 4; //High range for histogram MeV

double coshigh = 0.54; //0.995; //0.92 //0.54 is the highest angle that can hit UMB, CBEtop, CBEbot
double coslow = 1; //0.62 //0.8
double fitlow = 0.4;
double fithigh = 4;
const Int_t NBins = 50;
double betacut = 0.8; //Currently we're only doing a beta > 0 cutoff for real data. Beta > 0.8 recommended for sim

//Number of layers and strips
const int nstrips = 32;
const int nmods = 6;
const int nrows = 6;
const int nlayers = 7;
int dt[4] = {3,4,1,2};

int lyr[nlayers];
int rw[nrows];
int md[nmods];
int strps[nstrips];

for(int l = 0; l < nlayers; l++){lyr[l] = l;}
for(int r = 0; r < nrows; r++){rw[r] = r;}
for(int k = 0; k < nmods; k++){md[k] = k;}
for(int s = 0; s < nstrips; s++){strps[s] = s;}


//The Chosen One TM
const int chl = 0;
const int chr = 0;
const int chm = 0;
const int chd = 2;
const int chslow = 24;
const int chshi = 31;
int numParameter = 4;

//Just need one histogram and one fit for one strip
TH1F * hdet = new TH1F (TString::Format("h0_l%ir%im%id%i",chl,chr,chm,chd), ("Edep l" + to_string(chl) + "r" + to_string(chr) + "m" + to_string(chm) + "d" + to_string(chd)).c_str(), NBins, xlow,xhigh);
TF1 * gdet = new TF1("f_landau_gauss",langaufun,fitlow,fithigh,numParameter);
gdet->SetParameter(0, 0.1);
gdet->SetParameter(1, 0);
gdet->SetParameter(2, 1000);
gdet->SetParameter(3, 1);

//Prepare reconstronstructed event
CEventRec* Event = new CEventRec(); //New reconstructed event
TChain * TreeRec = new TChain("TreeRec"); //New TreeRec Tchain object (this is new to me)
TreeRec->SetBranchAddress("Rec", &Event); //Set the branch address using Event (defined above)
TreeRec->Add(FilenameRoot);

//Prepare textile for The Chosen One
std::ofstream chosenfile;
chosenfile.open("TheChosenDet.txt");
chosenfile << TString::Format( "The Chosen Det : l%ir%im%id%i ", chl,chr,chm,chd )  << endl;
chosenfile << TString::Format( "Filename : %s", FilenameRoot )  << endl;
chosenfile << TString::Format( "Beta Cut : %f", betacut) << endl;
chosenfile << TString::Format( "Total Entries : %i", static_cast<int>(TreeRec->GetEntries()/MainLoopScaleFactor)) << endl;
chosenfile << "Event \t LRMS \t EDep \t Cos(theta) \t Edep*Cos(theta)" << endl;
chosenfile.close();

//For mapping:
map<int, unsigned int> TofIndexVolumeIdMap;

//How many entries:
cout << "Total Number of events / Mainscale Factor = " << TreeRec->GetEntries()/MainLoopScaleFactor << endl;

//Using i to loop over every event in the tree
for(unsigned int i = 0; i < TreeRec->GetEntries(); i+=MainLoopScaleFactor){
	TreeRec->GetEntry(i);

	if( ((int)i % (int)ceil(TreeRec->GetEntries()/10)) == 0){
        cout << "Event number " << i << endl;
    }

	//Cuts are implemented in this chunk:
	if(Event->GetNTracks() == 1){  //First select the single track event
		bool Umbflag = 0;
		bool CBEtopflag = 0;
		bool CBEbotflag = 0;

		CTrackRec* pt = Event->GetPrimaryTrack();
		uint pt_index = 0;
   	        for( ; pt_index < Event->GetNTracks(); pt_index++) if( Event->GetTrack(pt_index)->IsPrimary() ) break;

		//Note downwards beta enforced by Event->GetPrimaryBeta() (should be positive) multiplied by Event->GetPrimaryMomentumDirection()[2] (z trajectory of particle)
		if(pt != nullptr && -fabs(Event->GetPrimaryMomentumDirectionGenerated().CosTheta()) > -coslow && -fabs(Event->GetPrimaryMomentumDirectionGenerated().CosTheta()) < -coshigh && Event->GetPrimaryBetaGenerated()*Event->GetPrimaryMomentumDirectionGenerated()[2] < 0 && fabs(Event->GetPrimaryBetaGenerated()) >  betacut){
			//cout << "Event is " << i << endl;

			//-----------EVENT LEVEL CUT APPLIED

			for(uint isig=0; isig<Event->GetTrack(0)->GetEnergyDeposition().size(); isig++){
                                unsigned int VolumeId  = Event->GetTrack(0)->GetVolumeId(isig); //Check the VolumeId of the event
                                if(volspec(VolumeId,0,3) == 100){ Umbflag = 1;} // cout << "UMB hit!" <<endl ;
                                if(volspec(VolumeId,0,3) == 110) {CBEtopflag = 1;}// cout << "CBE top hit!" << endl;
                                if(volspec(VolumeId,0,3) == 111) {CBEbotflag = 1;}// cout << "CBE bot hit!" << endl;
                        }

			if(Umbflag && CBEtopflag /*&& CBEbotflag*/ && (pt->GetChi2()/pt->GetNdof()) < 3.2 ){
				//cout << "Event number " << i << " passes the cuts!" << endl;
				for(uint isig=0; isig<Event->GetTrack(0)->GetEnergyDeposition().size(); isig++){
					unsigned int VolumeId  = Event->GetTrack(0)->GetVolumeId(isig);
					if(GGeometryObject::IsTrackerVolume(VolumeId) && Event->GetTrack(0)->GetEnergyDeposition(isig) > TrackerCut){

						int layer = GGeometryObject::GetTrackerLayer(VolumeId);

						int sdmod = GGeometryObject::GetLayerModule(VolumeId);
						int det = GGeometryObject::GetModuleDetector(VolumeId);
						int sdstrip = GGeometryObject::GetDetectorStrip(VolumeId);

						//cout << "SD: lmd = " << layer << sdmod << det << endl;

						int row = getrow(layer,sdmod);
						int mod = getmod(layer,sdmod);
						int strip = getch(layer, det, sdstrip);

						//cout << "Edep * Cos(theta) " << Event->GetTrack(0)->GetEnergyDeposition(isig)*fabs(Event->GetPrimaryMomentumDirection().CosTheta() << endl;
						//h[layer][row][mod][strip]->Fill(  (Event->GetTrack(0)->GetEnergyDeposition(isig)*fabs(Event->GetPrimaryMomentumDirection().CosTheta()))   );
						//hnentries->Fill(row*32+strip,layer*6+mod);

						if(layer == chl && row == chr && mod == chm && strip >= chslow && strip <= chshi){

							hdet->Fill((0.23/0.25)*Event->GetTrack(0)->GetEnergyDeposition(isig)*fabs(Event->GetPrimaryMomentumDirection().CosTheta()));

							chosenfile.open("TheChosenDet.txt",std::ios::app);
						chosenfile << TString::Format( "%i \t %i%i%i%i \t %f \t %f \t %f",i,layer,row,mod,strip,Event->GetTrack(0)->GetEnergyDeposition(isig),fabs(Event->GetPrimaryMomentumDirection().CosTheta()), (Event->GetTrack(0)->GetEnergyDeposition(isig)*fabs(Event->GetPrimaryMomentumDirection().CosTheta()))   ) << endl;
							chosenfile.close();
						}

					} //Closed bracket for Tracker volume and tracker cutoff

				} //Closed bracket for iteration over event with TOF cuts

			} //Closed bracket for if statement for cuts

			//-----------EVENT LEVEL CUTS END


		} //Closed bracket for event level cut

	} //Closed bracket for single track cut

} //Closed bracket for iteration through tree events, move on to the next event i



//Histogram section
//--------------------------------------


TCanvas * EdepCompare = new TCanvas("EdepCompare", "EdepCompare", 200, 10, 900, 900);
EdepCompare->SetLeftMargin(0.11);
EdepCompare->SetRightMargin(0.04);
EdepCompare->SetTopMargin(0.04);
TLegend* LegEdepCompare = new TLegend(0.5, 0.75, 0.95, 0.95);
LegEdepCompare->SetFillColor(0);

hdet->SetLineColor(1);
hdet->GetXaxis()->SetTitle("Energy Deposition of Hit (MeV)");
hdet->SaveAs("h300detprefit");
hdet->SaveAs("h300detprefit.root");

hdet->GetYaxis()->SetTitle("Number of Events");
hdet->Fit(gdet,"R");
gPad->SetGridx(1);
gPad->SetGridy(1);
gPad->SetLogy(1);
gStyle->SetTitleW(0.9);
gStyle->SetOptFit();
hdet->Draw();
hdet->SaveAs("h300det");
gdet->SaveAs("g300det");
hdet->SaveAs("h300det.root");
gdet->SaveAs("g300det.root");

string title = "l" + to_string(chl) + "r" + to_string(chr) + "m" + to_string(chm) + "d" + to_string(chd);
char name[400];
sprintf(name, "%s.root",title.c_str());
EdepCompare->SaveAs(name);
sprintf(name, "%s.png",title.c_str());
EdepCompare->SaveAs(name);
sprintf(name, "%s.pdf",title.c_str());
EdepCompare->SaveAs(name);



//--------------------------------------

cout << endl << "I am done" << endl;



return 1;

}
