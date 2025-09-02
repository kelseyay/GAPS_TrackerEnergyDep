using namespace std;

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
//using Crane::Calibration;
//To run, do ./EdepHist -i /home/kelsey/simulations/simdat/simnew/mu-_gaps_triggerlevel1_FTFP_BERT_1744342800_rec.root
//This function will take in a VolumeID number and then two other numbers that specify which part of the VolumeID you want to interrogate
//so ideally volspec(VolumeID,0,3) will give you the first three numbers of VolumeID which can tell you which CBE part it is.
//How to use, volspec(12345,1,4) outputs 234 as an integer for your comparing needs
int volspec(int volnum,int a, int b){
        stringstream ss;
        ss << volnum;
        return atoi(ss.str().substr(a, b).c_str());
}

//The two functions below aim to take the SimpleDet mod from 0-35 and translate that into a row and module
int getmod(int sdlayer, int sdmod){
        if(sdlayer % 2 == 0){
                return 5 - (sdmod % 6);
        }else{
                return 5 - (floor(sdmod/6));
        }
}

int getrow(int sdlayer, int sdmod){
        if(sdlayer % 2 == 0){
                return floor(sdmod/6);
        }else{
                return sdmod % 6;
        }
}

//This function aims to take layer, sddet, and sdstrip and translate that into a channel 0-31
int getch(int sdlayer, int sddet, int sdstrp){
	if(sdlayer % 2 == 0){
		if(sddet == 0){return sdstrp + 24;}
		if(sddet == 1){return 23-sdstrp;}
		if(sddet == 2){return sdstrp;}
		if(sddet == 3){return 15-sdstrp;}
		else{return -1;}
	}
	if(sdlayer % 2 != 0){
                if(sddet == 0){return 7 - sdstrp;}
                if(sddet == 1){return 31 - sdstrp;}
                if(sddet == 2){return sdstrp + 8;}
                if(sddet == 3){return sdstrp + 16;}
		else{return -1;}
	}
	else{return -1;}
}

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

char FilenameRoot[400];
sprintf(FilenameRoot,"%s*.root",reco_path.c_str());
//sprintf(FilenameRoot,reco_path.c_str()); //210 simu data on my computer!
cout << "File name is : " << FilenameRoot << endl;

CEventRec* Event = new CEventRec(); //New reconstructed event
TChain * TreeRec = new TChain("TreeRec"); //New TreeRec Tchain object (this is new to me)
TreeRec->SetBranchAddress("Rec", &Event); //Set the branch address using Event (defined above)
TreeRec->Add(FilenameRoot);

int MainLoopScaleFactor = 1; //Set this number to scale the step size. Larger means runs faster and fewer events
double TrackerCut = 0.3; //Threshold for an energy deposition to be considered a hit

double xlow = 0.3; //Low range for histogram MeV
double xhigh = 5; //High range for histogram MeV

double coshigh = 0.54; //0.995; //0.92 //0.54 is the highest angle that can hit UMB, CBEtop, CBEbot
double coslow = 1; //0.62 //0.8
double fitlow = 0.5;
double fithigh = 3.5;
const Int_t NBins = 50;
double betacut = 0.8; //Currently we're only doing a beta > 0 cutoff for real data. Beta > 0.8 recommended for sim

//Full tracker histogram range
double mpvmin = 0.66;
double mpvmax = 0.75;

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
for(int m = 0; m < nmods; m++){md[m] = m;}
for(int s = 0; s < nstrips; s++){strps[s] = s;}


//Need a histogram and fitting function for every strip
TH1F * h[nlayers][nrows][nmods][nstrips];

//auto hcol21 = new TH2F("hcol21","MPV Full Tracker",nrows*nstrips,0,nrows*nstrips,nlayers*nmods,0,nlayers*nmods);
auto hnentries = new TH2F("hnentries","Full Tracker Strip-Level NHits",nrows*nstrips,0,nrows*nstrips,nlayers*nmods,0,nlayers*nmods);

for(int l = 0;l<nlayers;l++){
        for(int r = 0;r<nrows;r++){
                for(int m = 0; m < nmods; m++){
                        for(int s = 0; s < nstrips; s++){
                                h[l][r][m][s] = new TH1F (TString::Format("h0_l%ir%im%is%i",l,r,m,s), ("Edep l" + to_string(lyr[l]) + "r" + to_string(rw[r]) + "m" + to_string(md[m]) + "s" + to_string(strps[s])).c_str(), NBins, xlow,xhigh);
                        }
                }
        }
}

//Prepare textile for saving values
std::ofstream myfile;
myfile.open("EdepList.txt");
myfile << TString::Format( "Filename : %s", reco_path.c_str() )  << endl;
myfile << TString::Format( "Beta Cut : %f", betacut) << endl;
myfile << TString::Format( "Total Entries : %i", static_cast<int>(TreeRec->GetEntries()/MainLoopScaleFactor)) << endl;
myfile << "Layer \t Row \t Mod \t Strip \t NEntries" << endl;
myfile.close();

//Prepare cuts:
map<int, unsigned int> TofIndexVolumeIdMap;

//For Plotting purposes
ca::GPlottingTools Plotting;
char text[400]; //This variable is used later to name the plots

//How many entries:
cout << "Total Number of events / Mainscale Factor = " << TreeRec->GetEntries()/MainLoopScaleFactor << endl;
cout << "Beta cut " << betacut << endl;
cout << "MPV range on 2D Full Tracker will be " << mpvmin << " - " << mpvmax << endl;

//Now we can go over the loop
TreeRec->GetEntry(0);

cout << "Total Number of events / Mainscale Factor = " << TreeRec->GetEntries()/MainLoopScaleFactor << endl;

//Using i to loop over every event in the tree
for(unsigned int i = 0; i < TreeRec->GetEntries(); i+=MainLoopScaleFactor){
    TreeRec->GetEntry(i);

	//Cuts are implemented in this chunk:
	if(Event->GetNTracks() == 1){  //First select the single track event
		bool Umbflag = 0;
		bool CBEtopflag = 0;
		bool CBEbotflag = 0;

		if( ((int)i % (int)ceil(TreeRec->GetEntries()/10)) == 0){
		    cout << "Event number " << i << endl;
		}

		CTrackRec* pt = Event->GetPrimaryTrack();
		uint pt_index = 0;
      	        for( ; pt_index < Event->GetNTracks(); pt_index++) if( Event->GetTrack(pt_index)->IsPrimary() ) break;

		//Note downwards beta enforced by Event->GetPrimaryBeta() (should be positive) multiplied by Event->GetPrimaryMomentumDirection()[2] (z trajectory of particle)
		if(pt != nullptr && -fabs(Event->GetPrimaryMomentumDirection().CosTheta()) > -1 && -fabs(Event->GetPrimaryMomentumDirection().CosTheta()) < 0 && Event->GetPrimaryBeta()*Event->GetPrimaryMomentumDirection()[2] < 0 && fabs(Event->GetPrimaryBeta()) >  0){

			if(/*Umbflag && CBEtopflag && CBEbotflag && (pt->GetChi2()/pt->GetNdof()) < 3.2 */ 1){
				//cout << "Event number " << i << " passes the cuts!" << endl;
				for(uint isig=0; isig<Event->GetTrack(0)->GetEnergyDeposition().size(); isig++){
					unsigned int VolumeId  = Event->GetTrack(0)->GetVolumeId(isig);
					if(GGeometryObject::IsTrackerVolume(VolumeId) && Event->GetTrack(0)->GetEnergyDeposition(isig) > TrackerCut){

						int layer = GGeometryObject::GetTrackerLayer(VolumeId);

						int sdmod = GGeometryObject::GetLayerModule(VolumeId);
						int det = GGeometryObject::GetModuleDetector(VolumeId);
						int sdstrip = GGeometryObject::GetDetectorStrip(VolumeId);

						int row = getrow(layer,sdmod);
						int mod = getmod(layer,sdmod);
						int strip = getch(layer, det, sdstrip);

						//cout << "lrms" << layer << row << mod << strip << endl;
						//cout << "Edep * Cos(theta) " << Event->GetTrack(0)->GetEnergyDeposition(isig)*fabs(Event->GetPrimaryMomentumDirection().CosTheta() << endl;
						if(layer < nlayers){ //This line prevents a segfault in the case of wanting to do fewer layers than the whole tracker
							h[layer][row][mod][strip]->Fill(  (Event->GetTrack(0)->GetEnergyDeposition(isig)*fabs(Event->GetPrimaryMomentumDirection().CosTheta()))   );
							hnentries->Fill(row*32+strip,layer*6+mod);
						}

					} //Closed bracket for Tracker volume and tracker cutoff

				} //Closed bracket for iteration over event with TOF cuts

			} //Closed bracket for if statement for cuts


			//-----------EVENT LEVEL CUTS END


		} //Closed bracket for event level cut

	} //Closed bracket for single track cut

}  //Closed bracket for iteration through tree events, move on to the next event i



//Histogram section
//--------------------------------------

myfile.open("EdepList.txt",std::ios::app);
TList *hlist = new TList();
TH1F *hnotes;
hnotes = new TH1F (TString::Format("Bcut%f",betacut),reco_path.c_str(),1,1,1);
hlist->Add(hnotes);

for(int l = 0; l<nlayers;l++){
        for(int r = 0; r<nrows;r++){
                for(int m = 0; m < nmods; m++){
			for(int s = 0; s < nstrips;s++){
				hlist->Add(h[l][r][m][s]);
				myfile << (TString::Format(    "%i \t %i \t %i \t %i \t %i \n",l,r,m,s, static_cast<int>(h[l][r][m][s]->GetEntries())   ));
			}
		}
        }
}

TFile *f = new TFile("histlist.root","RECREATE");
hlist->Write();
f->Close();


//Histogram for NEntries at a strip level
TCanvas * c2 = new TCanvas("c2", "c2", 200, 10, 900, 900);
c2->SetLeftMargin(0.1);
c2->SetRightMargin(0.16);
c2->SetTopMargin(0.1);
c2->SetBottomMargin(0.1);
hnentries->SetBit(TH1::kNoStats);
hnentries->GetXaxis()->SetTitle("row(0-6)*32 + strp (0-31)");
hnentries->GetYaxis()->SetTitle("layer(0-7)*6 + mod(0-6)");
hnentries->GetZaxis()->SetTitle("Number of Hits");
hnentries->Draw("COLZ");

//hnentries->SaveAs("hnentries.root");
//hnentries->SaveAs("hnentries");

string title = "HistFullTrackerNEntries";
char histname[400];
sprintf(histname, "%s.root",title.c_str());
c2->SaveAs(histname);
sprintf(histname, "%s.png",title.c_str());
c2->SaveAs(histname);

//--------------------------------------

cout << endl << "I am done" << endl;

return 1;

}
