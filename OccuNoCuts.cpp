//To run, do ./OccuNoCuts -i /home/kelsey/simulations/simdat/ground/251204/ethernet251204_0 -o test

using namespace std;

#include "KYtools.C"

using namespace Crane::Analysis;
namespace ca = Crane::Analysis;
namespace cl = Crane::Common;
//using Crane::Calibration;

int main(int argc, char *argv[]){

GOptionParser* parser = GOptionParser::GetInstance();
parser->AddProgramDescription("Occupancy plot with no cuts");
parser->AddCommandLineOption<string>("in_path", "path to instrument data files", "./*", "i");
parser->AddCommandLineOption<string>("out_file", "name of output root file", "", "o");
parser->AddCommandLineOption<int>("zsc", "Scale factor z axis", 10, "s");
parser->AddCommandLineOption<int>("MainloopScale", "Main loop scale factor",1,"m");
parser->ParseCommandLine(argc, argv);
parser->Parse();

string reco_path = parser->GetOption<string>("in_path");
string out_path = parser->GetOption<string>("out_file");
int MainLoopScaleFactor = parser->GetOption<int>("MainloopScale");
int zsc = parser->GetOption<int>("zsc");

cout << "out path " << out_path << endl;
if(out_path != "" && out_path[out_path.length()-1] != '/' ){ cout <<  "out path no slash!" << endl; out_path = out_path + '/'; }

cout << reco_path << endl;
char FilenameRoot[400];
sprintf(FilenameRoot,"%s*.root",reco_path.c_str());
cout << FilenameRoot << endl;

CEventRec* Event = new CEventRec(); //New reconstructed event
TChain * TreeRec = new TChain("TreeRec"); //New TreeRec Tchain object (this is new to me)
TreeRec->SetBranchAddress("Rec", &Event); //Set the branch address using Event (defined above)
TreeRec->Add(FilenameRoot);

double TrackerCut = 0.3; //Threshold for an energy deposition to be considered a hit
double TofCutLow = 0.1;

//Prepare cuts:
map<int, unsigned int> TofIndexVolumeIdMap;

//For Plotting purposes
ca::GPlottingTools Plotting;
char text[400]; //This variable is used later to name the plots

//Tracker dimensions
int nlayers = 7;
int nrows = 6;
int nmods = 6;
int nstrips = 32;

//All of the plots are declared here
//("Title",Number of bins,xmin,xmax,"xlabel","ylabel",ymin,ymax)

//2D Histos
//Currently just two histograms, one for Umbrella, one for CBE
TH2D* HTofUMBOccu = Plotting.DefineTH2D("HTofUMBOccu", 168, -2000, 2000, 168, -2000, 2000, "rec. hit position x [mm]", "rec. hit position y [mm]", "events", 2, TreeRec->GetEntries()/(MainLoopScaleFactor*zsc));
TH2D* HTofCBEtopOccu = Plotting.DefineTH2D("HTofCBEtopOccu", 25, -937.5, 937.5, 25, -937.5, 937.5, "rec. hit position x [mm]", "rec. hit position y [mm]", "events", 10, TreeRec->GetEntries()/(MainLoopScaleFactor*zsc));
TH2D* HTofCBEbotOccu = Plotting.DefineTH2D("HTofCBEbotOccu", 25, -937.5, 937.5, 25, -937.5, 937.5, "rec. hit position x [mm]", "rec. hit position y [mm]", "events", 10, TreeRec->GetEntries()/(MainLoopScaleFactor*zsc));

//auto hcol21 = new TH2F("hcol21","MPV Full Tracker",nrows*nstrips,0,nrows*nstrips,nlayers*nmods,0,nlayers*nmods);
auto hnentries = new TH2F("hnentries","Full Tracker Strip-Level NHits",nrows*nstrips,0,nrows*nstrips,nlayers*nmods,0,nlayers*nmods);

//Now we can go over the loop
TreeRec->GetEntry(0);

cout << "Total Number of events / Mainscale Factor = " << TreeRec->GetEntries()/MainLoopScaleFactor << endl;

//Using i to loop over every event in the tree
//for(unsigned int i = 0; i < 100; i+=MainLoopScaleFactor){
for(unsigned int i = 0; i < TreeRec->GetEntries(); i+=MainLoopScaleFactor){
        TreeRec->GetEntry(i);

        if( ((int)i % (int)ceil(TreeRec->GetEntries()/10) == 0) ){
		    cout << "Event number " << i << endl;
		}

        //cout << "Event is " << i << endl;
        for(uint isig=0; isig<Event->GetVolumeId().size(); isig++){
            int VolumeId = Event->GetVolumeId().at(isig);
            //cout << "Volume Id " << VolumeId << endl;

            if(GGeometryObject::IsTofVolume(VolumeId) && Event->GetHitSeries().at(isig).GetTotalEnergyDeposition() > TofCutLow){
                if(volspec(VolumeId,0,3) == 100){ //If hit is Umb, fill the Umb Occu plot
                    HTofUMBOccu->Fill(Event->GetHitSeries().at(isig).GetPosition().X(), Event->GetHitSeries().at(isig).GetPosition().Y());
                    //HTofUmbOccu->Fill(Event->GetTrack(0)->GetPosition(k).X()+Event->GetTrack(0)->GetPositionResidual(k).X(), Event->GetTrack(0)->GetPosition(k).Y()+Event->GetTrack(0)->GetPositionResidual(k).Y());
                }
                if(volspec(VolumeId,0,3) == 110){ //If hit is CBEtop, fill the CBEtop Occu plot
                    HTofCBEtopOccu->Fill(Event->GetHitSeries().at(isig).GetPosition().X(), Event->GetHitSeries().at(isig).GetPosition().Y());
                    //HTofCBEtopOccu->Fill(Event->GetTrack(0)->GetPosition(k).X()+Event->GetTrack(0)->GetPositionResidual(k).X(), Event->GetTrack(0)->GetPosition(k).Y()+Event->GetTrack(0)->GetPositionResidual(k).Y());
                }
                if(volspec(VolumeId,0,3) == 111){ //If the hit is CBEbot, fill the Cbebot Occu plot
                    HTofCBEbotOccu->Fill(Event->GetHitSeries().at(isig).GetPosition().X(), Event->GetHitSeries().at(isig).GetPosition().Y());
                }
            }

            if(GGeometryObject::IsTrackerVolume(VolumeId) && Event->GetHitSeries().at(isig).GetTotalEnergyDeposition() > TrackerCut){

				int layer = GGeometryObject::GetTrackerLayer(VolumeId);
				int sdmod = GGeometryObject::GetLayerModule(VolumeId);
				int det = GGeometryObject::GetModuleDetector(VolumeId);
				int sdstrip = GGeometryObject::GetDetectorStrip(VolumeId);

				int row = getrow(layer,sdmod);
				int mod = getmod(layer,sdmod);
				int strip = getch(layer, det, sdstrip);

				if(layer < nlayers){ //This line prevents a segfault in the case of wanting to do fewer layers than the whole tracker
					hnentries->Fill(row*32+strip,layer*6+mod);
				}

			} //Closed bracket for Tracker volume and tracker cutoff



            //cout << "Hit " << isig << " VolumrId " << Event->GetVolumeId().at(isig) << " Edep " << Event->GetHitSeries().at(isig).GetTotalEnergyDeposition() << endl;
        }

        //No cuts at all, just loop over the hit series.

}

//void histplot2d(string ctitle, TH2D* h1, string title, string xtitle, string ytitle, string ztitle, string savename){


//Save all histograms so they can be reopened.
HTofUMBOccu->SaveAs( (out_path + "HTofUMBOccu.root").c_str() );
HTofCBEtopOccu->SaveAs( (out_path + "HTofCBEtopOccu.root").c_str() );
HTofCBEbotOccu->SaveAs( (out_path + "HTofCBEbotOccu.root").c_str() );
hnentries->SaveAs( (out_path + "NhitsTracker.root").c_str() );

histplot2d("c1",HTofUMBOccu,"Umbrella Occupancy Plot","X Location Umb Hit","Y Location Umb Hit","NEntries",out_path+"TofUmbOccu");
histplot2d("c2",HTofCBEtopOccu,"CBE Top Occupancy Plot","X Location Umb Hit","Y Location Umb Hit","NEntries",out_path+"TofCBEtopOccu");
histplot2d("c3",HTofCBEbotOccu,"CBE Bot Occupancy Plot","X Location Umb Hit","Y Location Umb Hit","NEntries",out_path+"TofCBEbotOccu");
histplot2f("c4",hnentries,"Full Tracker Strip Level NHits","row(0-5)*32 + strip(0-31)","layer(0-5)*6 + mod(0-5)","NEntries",out_path+"TrackerEntries");

cout << endl << "I am done" << endl;
return 1;

}
