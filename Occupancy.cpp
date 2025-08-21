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
//To run, do ./Occupancy -i /home/kelsey/simulations/simdat/simnew/mu-_gaps_triggerlevel1_FTFP_BERT_1744342800_rec.root
//This function will take in a VolumeID number and then two other numbers that specify which part of the VolumeID you want to interrogate
//so ideally volspec(VolumeID,0,3) will give you the first three numbers of VolumeID which can tell you which CBE part it is.
//How to use, volspec(12345,1,4) outputs 234 as an integer for your comparing needs
int volspec(int volnum,int a, int b){
        stringstream ss;
        ss << volnum;
        return atoi(ss.str().substr(a, b).c_str());
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

CEventRec* Event = new CEventRec(); //New reconstructed event
TChain * TreeRec = new TChain("TreeRec"); //New TreeRec Tchain object (this is new to me)
TreeRec->SetBranchAddress("Rec", &Event); //Set the branch address using Event (defined above)
TreeRec->Add(reco_path.c_str());

int MainLoopScaleFactor = 1; //Set this number to scale the step size. Larger means runs faster and fewer events
double TrackerCut = 0.4; //Threshold for an energy deposition to be considered a hit
double TofCutLow = 0;
double xlow = 0.3; //Hist range
double xhigh = 6; //Hist range
double coshigh = 0; //0.92; //0.995; //0.92 //0.54 is the highest angle that can hit UMB, CBEtop, CBEbot
double coslow = 1; //0.62 //0.8
double eventbetacut = 0.2; //Cut that is applied to all events
double betacut = 0.9; //Separation of quandrants in the beta plot
const Int_t NBins = 50;

//Prepare cuts:
map<int, unsigned int> TofIndexVolumeIdMap;

//For Plotting purposes
ca::GPlottingTools Plotting;
char text[400]; //This variable is used later to name the plots

//All of the plots are declared here
//("Title",Number of bins,xmin,xmax,"xlabel","ylabel",ymin,ymax)

//2D Histos
//Currently just two histograms, one for Umbrella, one for CBE
TH2D* HTofUmbOccu = Plotting.DefineTH2D("HTofUmbOccu", 30, -2000, 2000, 30, -2000, 2000, "rec. hit position x [mm]", "rec. hit position y [mm]", "events", 0.5, 1e3);
TH2D* HTofCBEtopOccu = Plotting.DefineTH2D("HTofCBEtopOccu", 30, -2000, 2000, 30, -2000, 2000, "rec. hit position x [mm]", "rec. hit position y [mm]", "events", 0.5, 1e3);

//Now we can go over the loop
TreeRec->GetEntry(0);

cout << "Total Number of events / Mainscale Factor = " << TreeRec->GetEntries()/MainLoopScaleFactor << endl;

//Using i to loop over every event in the tree
for(unsigned int i = 0; i < TreeRec->GetEntries(); i+=MainLoopScaleFactor){
        TreeRec->GetEntry(i);

        //Cuts are implemented in this chunk:
        if(Event->GetNTracks() == 1){  //First select the single track event

                //cout << "Single Track Event " << endl;
                bool Umbflag = 0;
                bool CBEtopflag = 0;
                bool CBEbotflag = 0;

                CTrackRec* pt = Event->GetPrimaryTrack();
                uint pt_index = 0;
                for( ; pt_index < Event->GetNTracks(); pt_index++) if( Event->GetTrack(pt_index)->IsPrimary() ) break;
        //Apply event-level cuts
                if(pt != nullptr && -fabs(Event->GetPrimaryMomentumDirection().CosTheta()) > -coslow && -fabs(Event->GetPrimaryMomentumDirection().CosTheta()) < -coshigh && Event->GetPrimaryBeta()*Event->GetPrimaryMomentumDirection()[2] < 0 && fabs(Event->GetPrimaryBeta()) > eventbetacut){

                        //cout << endl << "Event is " << i << endl;

                        //First loop over the event for your needed flags and variables
                        for(unsigned int isig = 0; isig < Event->GetTrack(0)->GetEnergyDeposition().size(); isig++){
                                unsigned int VolumeId  = Event->GetTrack(0)->GetVolumeId(isig); //Event->GetVolumeId().at(isig); //Check the VolumeId of the event
                                if(volspec(VolumeId,0,3) == 100){ Umbflag = 1;} // cout << "UMB hit!" <<endl ;
                                if(volspec(VolumeId,0,3) == 110) {CBEtopflag = 1;}// cout << "CBE top hit!" << endl;
                                if(volspec(VolumeId,0,3) == 111) {CBEbotflag = 1;}// cout << "CBE bot hit!" << endl;
                        }

                        //If the desired flags are checked, proceed to fill the relevant histograms.
                        if(Umbflag && CBEtopflag && CBEbotflag && (pt->GetChi2()/pt->GetNdof()) < 3.2){
                                //cout << "Event is " << i << endl;

                                //Loop over the events on track 0
                                for(unsigned int k = 0; k < Event->GetTrack(0)->GetEnergyDeposition().size(); k++){
                                unsigned int VolumeId = Event->GetTrack(0)->GetVolumeId(k); //For each hit, check the volume ID
                                        //Check if the volume is a TOF volume and make sure you pass the low energy criteria
                                        if(GGeometryObject::IsTofVolume(VolumeId) && Event->GetTrack(0)->GetEnergyDeposition(k) > TofCutLow){
                                                if(volspec(VolumeId,0,3) == 100){ //If hit is Umb, fill the Umb Occu plot
                                                        HTofUmbOccu->Fill(Event->GetTrack(0)->GetPosition(k).X()+Event->GetTrack(0)->GetPositionResidual(k).X(), Event->GetTrack(0)->GetPosition(k).Y()+Event->GetTrack(0)->GetPositionResidual(k).Y());
                                                }
                                                if(volspec(VolumeId,0,3) == 110){ //If hit is CBEtop, fill the CBEtop Occu plot
                                                        HTofCBEtopOccu->Fill(Event->GetTrack(0)->GetPosition(k).X()+Event->GetTrack(0)->GetPositionResidual(k).X(), Event->GetTrack(0)->GetPosition(k).Y()+Event->GetTrack(0)->GetPositionResidual(k).Y());
                                                }

                                        }
                                }
                        }

                }
        }
}

TCanvas * c1 = new TCanvas("c1", "c1", 200, 10, 900, 900);
c1->SetLeftMargin(0.15);
c1->SetRightMargin(0.17);
c1->SetTopMargin(0.11);
c1->SetBottomMargin(0.11);
HTofUmbOccu->SetTitle("Umbrella Occupancy Plot");
HTofUmbOccu->GetXaxis()->SetTitle("X Location Umb Hit");
HTofUmbOccu->GetYaxis()->SetTitle("Y Location Umb Hit");
HTofUmbOccu->GetYaxis()->SetTitleOffset(2);
HTofUmbOccu->GetZaxis()->SetTitle("Number of Entries");

HTofUmbOccu->Draw("COLZ");
//gStyle->SetTitleAlign(33);
//gStyle->SetTitleY(.99);
//gStyle->SetTitleX(.80);
gPad->SetLogz();

char histname[400];
string TofUmbTitle = "TofUmbOccu";
sprintf(histname, "%s.root",TofUmbTitle.c_str());
c1->SaveAs(histname);
sprintf(histname, "%s.png",TofUmbTitle.c_str());
c1->SaveAs(histname);
sprintf(histname, "%s.pdf",TofUmbTitle.c_str());
c1->SaveAs(histname);

TCanvas * c2 = new TCanvas("c2", "c2", 200, 10, 900, 900);
c2->SetLeftMargin(0.15);
c2->SetRightMargin(0.17);
c2->SetTopMargin(0.11);
c2->SetBottomMargin(0.11);
HTofCBEtopOccu->SetTitle("CBE Top Occupancy Plot");
HTofCBEtopOccu->GetXaxis()->SetTitle("X Location CBE Top Hit");
HTofCBEtopOccu->GetYaxis()->SetTitle("Y Location CBE Top Hit");
HTofCBEtopOccu->GetYaxis()->SetTitleOffset(2);
HTofCBEtopOccu->GetZaxis()->SetTitle("Number of Entries");
HTofCBEtopOccu->Draw("COLZ");
gPad->SetLogz();

string TofCBEtopTitle = "TofCBEtopOccu";
sprintf(histname, "%s.root",TofCBEtopTitle.c_str());
c2->SaveAs(histname);
sprintf(histname, "%s.png",TofCBEtopTitle.c_str());
c2->SaveAs(histname);
sprintf(histname, "%s.pdf",TofCBEtopTitle.c_str());
c2->SaveAs(histname);



cout << endl << "I am done" << endl;
return 1;

}
