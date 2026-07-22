//To run, do /home/kelsey/GAPS_TrackerEnergyDep/build/HoursOccu -i /home/kelsey/simulations/simdat/flight/251221/FPSI/starlink251221_0

//After this project exploded into what it is now, I really should have used an array of histograms, maybe one for each section, like:
//  HCOR[4] rip but it's fine.

using namespace std;

#include "KYtools.C"

using namespace Crane::Analysis;
namespace ca = Crane::Analysis;
namespace cl = Crane::Common;
//using Crane::Calibration;

void format_hist2d(TH2D* h1, string title, string xtitle, string ytitle, string ztitle, int zsc){
    h1->GetZaxis()->SetRangeUser(2,h1->GetEntries()/zsc);
    h1->SetTitle(title.c_str());
    h1->SetBit(TH1::kNoStats);
    h1->GetXaxis()->SetTitle(xtitle.c_str());
    h1->GetYaxis()->SetTitle(ytitle.c_str());
    h1->GetZaxis()->SetTitle(ztitle.c_str());
}

void format_histplot2f(TH2F* h1, string title, string xtitle, string ytitle, string ztitle){
    h1->SetTitle(title.c_str());
    h1->SetBit(TH1::kNoStats);
    h1->GetXaxis()->SetTitle(xtitle.c_str());
    h1->GetYaxis()->SetTitle(ytitle.c_str());
    h1->GetZaxis()->SetTitle(ztitle.c_str());

}

void format_pad(TCanvas * c1, int i){
    c1->GetPad(i)->cd();
    c1->GetPad(i)->SetLeftMargin(0.14);
    c1->GetPad(i)->SetRightMargin(0.16);
    c1->GetPad(i)->SetTopMargin(0.1);
    c1->GetPad(i)->SetBottomMargin(0.1);

    gPad->SetLogz();
}

int main(int argc, char *argv[]){

GOptionParser* parser = GOptionParser::GetInstance();
parser->AddProgramDescription("Occupancy plot with no cuts");
parser->AddCommandLineOption<string>("in_path", "path to instrument data files", "./*", "i");
parser->AddCommandLineOption<string>("out_file", "name of output root file", "", "o");
parser->AddCommandLineOption<string>("name_pfx", "additional title to files", "", "n");
parser->AddCommandLineOption<int>("zsc", "Scale factor z axis", 10, "s");
parser->AddCommandLineOption<int>("MainloopScale", "Main loop scale factor",1,"m");
parser->ParseCommandLine(argc, argv);
parser->Parse();

string reco_path = parser->GetOption<string>("in_path");
string out_path = parser->GetOption<string>("out_file");
string pfx = parser->GetOption<string>("name_pfx");
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

int Npaddles = 12;
int Hlen = Npaddles*3;

//2D Histos
//Currently just two histograms, one for Umbrella, one for CBE
TH2D* HTofUMBOccu = Plotting.DefineTH2D("HTofUMBOccu", 168, -2000, 2000, 168, -2000, 2000, "rec. hit position x [mm]", "rec. hit position y [mm]", "events", 2, 10000);
TH2D* HTofCBEtopOccu = Plotting.DefineTH2D("HTofCBEtopOccu", 25, -937.5, 937.5, 25, -937.5, 937.5, "rec. hit position x [mm]", "rec. hit position y [mm]", "events", 10, 10000);
TH2D* HTofCBEbotOccu = Plotting.DefineTH2D("HTofCBEbotOccu", 25, -937.5, 937.5, 25, -937.5, 937.5, "rec. hit position x [mm]", "rec. hit position y [mm]", "events", 10, 10000);
TH2D* HTofCBEtopOccu_test = Plotting.DefineTH2D("HTofCBEtopOccu_test", Hlen, 0, Hlen,  Hlen, 0, Hlen, "CBE Top Panel X Span", "CBE Top Panel Y Span", "events", 10, 10000);

TH2D* HTofCOR_XOccu = Plotting.DefineTH2D("HTofCOR_XOccu", 50, -1200, 1200, 49, -175, 1540, "rec. hit position y [mm]", "rec. hit position z [mm]", "events", 10, 10000);
TH2D* HTofCOR_min_XOccu = Plotting.DefineTH2D("HTofCOR_min_XOccu", 50, -1200, 1200, 49, -175, 1540, "rec. hit position y [mm]", "rec. hit position z [mm]", "events", 10, 10000);
TH2D* HTofCOR_YOccu = Plotting.DefineTH2D("HTofCOR_YOccu", 50, -1200, 1200, 49, -175, 1540, "rec. hit position x [mm]", "rec. hit position z [mm]", "events", 10, 10000);
TH2D* HTofCOR_min_YOccu = Plotting.DefineTH2D("HTofCOR_min_YOccu", 50, -1200, 1200, 49, -175, 1540, "rec. hit position x [mm]", "rec. hit position z [mm]", "events", 10, 10000);

TH2D* HTofCBE_XOccu = Plotting.DefineTH2D("HTofCBE_XOccu", 20, -800, 800, 36, -120, 1300, "rec. hit position y [mm]", "rec. hit position z [mm]", "events", 10, 10000);
TH2D* HTofCBE_min_XOccu = Plotting.DefineTH2D("HTofCBE_min_XOccu", 20, -800, 800, 36, -120, 1300, "rec. hit position y [mm]", "rec. hit position z [mm]", "events", 10, 10000);
TH2D* HTofCBE_YOccu = Plotting.DefineTH2D("HTofCBE_YOccu", 20, -800, 800, 36, -120, 1300,  "rec. hit position x [mm]", "rec. hit position z [mm]", "events", 10, 10000);
TH2D* HTofCBE_min_YOccu = Plotting.DefineTH2D("HTofCBE_min_YOccu", 20, -800, 800, 36, -120, 1300,  "rec. hit position x [mm]", "rec. hit position z [mm]", "events", 10, 10000);

//auto hcol21 = new TH2F("hcol21","MPV Full Tracker",nrows*nstrips,0,nrows*nstrips,nlayers*nmods,0,nlayers*nmods);
auto hnentries = new TH2F("hnentries","Full Tracker Strip-Level NHits",nrows*nstrips,0,nrows*nstrips,nlayers*nmods,0,nlayers*nmods);

//Prepare textile for saving values
std::ofstream myfile;
myfile.open(out_path + "OccuNoCuts.txt");
myfile << TString::Format( "Filename : %s", reco_path.c_str() )  << endl;
myfile << "NEvents = " << TreeRec->GetEntries()/MainLoopScaleFactor << endl;
myfile.close();

//Now we can go over the loop
TreeRec->GetEntry(0);
TTimeStamp T(Event->GetEventTime());
cout << "First event time " << Event->GetEventTime() << endl;
cout << "First event time formatted " << T.AsString("s") << endl;
string time_previous = T.AsString("s");

TreeRec->GetEntry(TreeRec->GetEntries()-1);
TTimeStamp T2(Event->GetEventTime());
//T = Event->GetEventTime();
cout << "Last event time " << Event->GetEventTime() << endl;
cout << "Last event time formatted " << T2.AsString("s") << endl;

int nhours = floor( (T2-T) / 3600);
int hour_flag = 0; //When this flag gets tripped, hour_flag++, wait for the next hour before saving all of the variables.

cout << "Number of hours " << nhours << endl;
cout << "Total Number of events / Mainscale Factor = " << TreeRec->GetEntries()/MainLoopScaleFactor << endl;

TCanvas * c_flat[nhours+1];
TCanvas * c_sides[nhours+1];

TString oo = out_path+pfx+"_HG_Occu_Hour_by_Hour.pdf";

//Using i to loop over every event in the tree
//for(unsigned int i = 0; i < 100; i+=MainLoopScaleFactor){
for(unsigned int i = 0; i < TreeRec->GetEntries(); i+=MainLoopScaleFactor){
        TreeRec->GetEntry(i);
        T2 = Event->GetEventTime();
        int Hour = floor( (T2-T) / 3600);

        if( ((int)i % (int)ceil(TreeRec->GetEntries()/10) == 0) ){
		    cout << "Event number " << i << endl;
		}

        if(Hour != hour_flag){
            cout << "hour_flag " << hour_flag << endl;
            TreeRec->GetEntry(i-10);
            TTimeStamp T_hourend(Event->GetEventTime());
            string timestamp_end = T_hourend.AsString("s");
            TreeRec->GetEntry(i);
            string time_thisstep = T2.AsString("s");
            string time = time_previous + " to " + timestamp_end.substr(timestamp_end.length() - 8);
            time_previous = time_thisstep;
            cout << time << endl;

            format_histplot2f(hnentries,"Tracker " + time,"row(0-5)*32 + strip(0-31)","layer(0-5)*6 + mod(0-5)","NEntries");
            format_hist2d(HTofUMBOccu,"UMB " + time,"X Location [cm]","Y Location [cm]","NEntries",1000);
            format_hist2d(HTofCBEtopOccu,"CBEtop " + time,"X Location [cm]","Y Location [cm]","NEntries",100);
            format_hist2d(HTofCBEbotOccu,"CBEbot " + time,"X Location [cm]","Y Location [cm]","NEntries",100);
            format_hist2d(HTofCOR_XOccu,"COR +X " + time,"Y Location [cm]","Z Location [cm]","NEntries",100);
            format_hist2d(HTofCOR_min_XOccu,"COR -X " + time,"Y Location [cm]","Z Location [cm]","NEntries",100);
            format_hist2d(HTofCOR_YOccu,"COR +Y " + time,"X Location [cm]","Z Location [cm]","NEntries",100);
            format_hist2d(HTofCOR_min_YOccu,"COR -Y " + time,"X Location [cm]","Z Location [cm]","NEntries",100);

            format_hist2d(HTofCBE_XOccu,"CBE +X " + time,"Y Location [cm]","Z Location [cm]","NEntries",100);
            format_hist2d(HTofCBE_min_XOccu,"CBE -X " + time,"Y Location [cm]","Z Location [cm]","NEntries",100);
            format_hist2d(HTofCBE_YOccu,"CBE +Y " + time,"X Location [cm]","Z Location [cm]","NEntries",100);
            format_hist2d(HTofCBE_min_YOccu,"CBE -Y " + time,"X Location [cm]","Z Location [cm]","NEntries",100);

            //format_hist2f(hnentries,"Full Tracker Strip Level NHits","row(0-5)*32 + strip(0-31)","layer(0-5)*6 + mod(0-5)","NEntries");

            //new TCanvas(Form("c%d", i), Form("Canvas %d", i), 200 + i*50, 200 + i*50, 600, 400);
            c_flat[hour_flag] = new TCanvas(Form("c_flat%d", hour_flag),Form("UMB CBEtop CBEbot TKR%d", hour_flag),600,600); //Four flat structures on one canvas
            c_sides[hour_flag] = new TCanvas(Form("c_sides%d", hour_flag),Form("COR CBE sides%d", hour_flag),1200,600); //Side TOF panels take up 8 (4x2) pieces of one canvas


            c_flat[hour_flag]->cd();
            c_flat[hour_flag]->Divide(2,2);

            format_pad(c_flat[hour_flag],1);
            HTofUMBOccu->Draw("COLZ");
            format_pad(c_flat[hour_flag],2);
            hnentries->Draw("COLZ");
            format_pad(c_flat[hour_flag],3);
            HTofCBEtopOccu->Draw("COLZ");
            format_pad(c_flat[hour_flag],4);
            HTofCBEbotOccu->Draw("COLZ");

            c_sides[hour_flag]->cd();
            c_sides[hour_flag]->Divide(4,2);

            format_pad(c_sides[hour_flag],1);
            HTofCOR_XOccu->Draw("COLZ");
            format_pad(c_sides[hour_flag],2);
            HTofCOR_min_XOccu->Draw("COLZ");

            format_pad(c_sides[hour_flag],3);
            HTofCOR_YOccu->Draw("COLZ");
            format_pad(c_sides[hour_flag],4);
            HTofCOR_min_YOccu->Draw("COLZ");
            format_pad(c_sides[hour_flag],5);
            HTofCBE_XOccu->Draw("COLZ");
            format_pad(c_sides[hour_flag],6);
            HTofCBE_min_XOccu->Draw("COLZ");
            format_pad(c_sides[hour_flag],7);
            HTofCBE_YOccu->Draw("COLZ");
            format_pad(c_sides[hour_flag],8);
            HTofCBE_min_YOccu->Draw("COLZ");

            if(hour_flag == 0){c_flat[hour_flag]->Print(Form("%s[",oo.Data())); c_flat[hour_flag]->Print(Form("%s",oo.Data())); c_sides[hour_flag]->Print(Form("%s",oo.Data()));
            }else{c_flat[hour_flag]->Print(Form("%s",oo.Data()));c_sides[hour_flag]->Print(Form("%s",oo.Data()));}

            HTofUMBOccu->Reset("ICESM");
            hnentries->Reset("ICESM");
            HTofCBEtopOccu->Reset("ICESM");
            HTofCBEbotOccu->Reset("ICESM");
            HTofCOR_XOccu->Reset("ICESM");
            HTofCOR_min_XOccu->Reset("ICESM");
            HTofCOR_YOccu->Reset("ICESM");
            HTofCOR_min_YOccu->Reset("ICESM");
            HTofCBE_XOccu->Reset("ICESM");
            HTofCBE_min_XOccu->Reset("ICESM");
            HTofCBE_YOccu->Reset("ICESM");
            HTofCBE_min_YOccu->Reset("ICESM");

            hour_flag++;

        }

        //cout << "Event is " << i << endl;

       	for(unsigned int k = 0; k < Event->GetTriggerVolumeId().size(); k++){
            unsigned int VolumeId = Event->GetTriggerVolumeId().at(k);
            if(volspec(VolumeId,0,3) == 110){ //If hit is CBEtop, fill the CBEtop Occu plot
            for(int k = 0; k < Hlen; k++){
                HTofCBEtopOccu_test->Fill(Hlen - (3*volspec(VolumeId,5,2)+2), k, 1/(float)Hlen);
            }
            }
        }

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

                //CBE filling:
                if(volspec(VolumeId,0,3) == 114){ //If the hit is COR +Y, fill the COR +Y Occu plot
                    HTofCBE_YOccu->Fill(Event->GetHitSeries().at(isig).GetPosition().X(), Event->GetHitSeries().at(isig).GetPosition().Z());
                    //cout << "Event is " << i << " COR +Y Hit! " << endl;
                }

                if(volspec(VolumeId,0,3) == 115){ //If the hit is COR -X, fill the COR -X Occu plot
                    HTofCBE_min_YOccu->Fill(Event->GetHitSeries().at(isig).GetPosition().X(), Event->GetHitSeries().at(isig).GetPosition().Z());
                    //cout << "Event is " << i << " COR +Y Hit! " << endl;
                }

                if(volspec(VolumeId,0,3) == 112){ //If the hit is COR +X, fill the COR +X Occu plot
                    HTofCBE_XOccu->Fill(Event->GetHitSeries().at(isig).GetPosition().Y(), Event->GetHitSeries().at(isig).GetPosition().Z());
                    //cout << "Event is " << i << " COR +X Hit! " << endl;
                }
                if(volspec(VolumeId,0,4) == 1160 || volspec(VolumeId,0,4) == 1161){ //If the hit is COR +X, fill the COR +X Occu plot
                    HTofCBE_XOccu->Fill(Event->GetHitSeries().at(isig).GetPosition().Y(), Event->GetHitSeries().at(isig).GetPosition().Z());
                    //cout << "Event is " << i << " COR +X 3PP Hit! " << endl;
                }
                if(volspec(VolumeId,0,3) == 113){ //If the hit is COR -X, fill the COR -X Occu plot
                    HTofCBE_min_XOccu->Fill(Event->GetHitSeries().at(isig).GetPosition().Y(), Event->GetHitSeries().at(isig).GetPosition().Z());
                    //cout << "Event is " << i << " COR +X Hit! " << endl;
                }
                if(volspec(VolumeId,0,4) == 1162 || volspec(VolumeId,0,4) == 1163){ //If the hit is COR +X, fill the COR +X Occu plot
                    HTofCBE_min_XOccu->Fill(Event->GetHitSeries().at(isig).GetPosition().Y(), Event->GetHitSeries().at(isig).GetPosition().Z());
                    //cout << "Event is " << i << " COR +X 3PP Hit! " << endl;
                }


                //COR filling:
                if(volspec(VolumeId,0,3) == 104){ //If the hit is COR +Y, fill the COR +Y Occu plot
                    HTofCOR_YOccu->Fill(Event->GetHitSeries().at(isig).GetPosition().X(), Event->GetHitSeries().at(isig).GetPosition().Z());
                    //cout << "Event is " << i << " COR +Y Hit! " << endl;
                }

                if(volspec(VolumeId,0,3) == 105){ //If the hit is COR -X, fill the COR -X Occu plot
                    HTofCOR_min_YOccu->Fill(Event->GetHitSeries().at(isig).GetPosition().X(), Event->GetHitSeries().at(isig).GetPosition().Z());
                    //cout << "Event is " << i << " COR +Y Hit! " << endl;
                }

                if(volspec(VolumeId,0,3) == 102){ //If the hit is COR +X, fill the COR +X Occu plot
                    HTofCOR_XOccu->Fill(Event->GetHitSeries().at(isig).GetPosition().Y(), Event->GetHitSeries().at(isig).GetPosition().Z());
                    //cout << "Event is " << i << " COR +X Hit! " << endl;
                }
                if(volspec(VolumeId,0,4) == 1060 || volspec(VolumeId,0,4) == 1061){ //If the hit is COR +X, fill the COR +X Occu plot
                    HTofCOR_XOccu->Fill(Event->GetHitSeries().at(isig).GetPosition().Y(), Event->GetHitSeries().at(isig).GetPosition().Z());
                    //cout << "Event is " << i << " COR +X 3PP Hit! " << endl;
                }

                if(volspec(VolumeId,0,3) == 103){ //If the hit is COR -X, fill the COR -X Occu plot
                    HTofCOR_min_XOccu->Fill(Event->GetHitSeries().at(isig).GetPosition().Y(), Event->GetHitSeries().at(isig).GetPosition().Z());
                    //cout << "Event is " << i << " COR +X Hit! " << endl;
                }
                if(volspec(VolumeId,0,4) == 1062 || volspec(VolumeId,0,4) == 1063){ //If the hit is COR +X, fill the COR +X Occu plot
                    HTofCOR_min_XOccu->Fill(Event->GetHitSeries().at(isig).GetPosition().Y(), Event->GetHitSeries().at(isig).GetPosition().Z());
                    //cout << "Event is " << i << " COR +X 3PP Hit! " << endl;
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



TreeRec->GetEntry(TreeRec->GetEntries()-1);
T2 = Event->GetEventTime();

string time_thisstep = T2.AsString("s");
string time = time_previous + " to " + time_thisstep.substr(time_thisstep.length() - 8);

cout << time << endl;
format_histplot2f(hnentries,"Tracker " + time,"row(0-5)*32 + strip(0-31)","layer(0-5)*6 + mod(0-5)","NEntries");
format_hist2d(HTofUMBOccu,"UMB " + time,"X Location [cm]","Y Location [cm]","NEntries",1000);
format_hist2d(HTofCBEtopOccu,"CBEtop " + time,"X Location [cm]","Y Location [cm]","NEntries",100);
format_hist2d(HTofCBEbotOccu,"CBEbot " + time,"X Location [cm]","Y Location [cm]","NEntries",100);
format_hist2d(HTofCOR_XOccu,"COR +X " + time,"Y Location [cm]","Z Location [cm]","NEntries",100);
format_hist2d(HTofCOR_min_XOccu,"COR -X " + time,"Y Location [cm]","Z Location [cm]","NEntries",100);
format_hist2d(HTofCOR_YOccu,"COR +Y " + time,"X Location [cm]","Z Location [cm]","NEntries",100);
format_hist2d(HTofCOR_min_YOccu,"COR -Y " + time,"X Location [cm]","Z Location [cm]","NEntries",100);

format_hist2d(HTofCBE_XOccu,"CBE +X " + time,"Y Location [cm]","Z Location [cm]","NEntries",100);
format_hist2d(HTofCBE_min_XOccu,"CBE -X " + time,"Y Location [cm]","Z Location [cm]","NEntries",100);
format_hist2d(HTofCBE_YOccu,"CBE +Y " + time,"X Location [cm]","Z Location [cm]","NEntries",100);
format_hist2d(HTofCBE_min_YOccu,"CBE -Y " + time,"X Location [cm]","Z Location [cm]","NEntries",100);
//format_hist2f(hnentries,"Full Tracker Strip Level NHits","row(0-5)*32 + strip(0-31)","layer(0-5)*6 + mod(0-5)","NEntries");

//new TCanvas(Form("c%d", i), Form("Canvas %d", i), 200 + i*50, 200 + i*50, 600, 400);
c_flat[nhours] = new TCanvas(Form("c_flat%d", nhours),Form("UMB CBEtop CBEbot TKR%d", nhours),600,600); //Four flat structures on one canvas
c_sides[nhours] = new TCanvas(Form("c_sides%d", nhours),Form("COR CBE sides%d", nhours),1200,600); //Side TOF panels take up 8 (4x2) pieces of one canvas


c_flat[nhours]->cd();
c_flat[nhours]->Divide(2,2);

format_pad(c_flat[nhours],1);
HTofUMBOccu->Draw("COLZ");
format_pad(c_flat[nhours],2);
hnentries->Draw("COLZ");
format_pad(c_flat[nhours],3);
HTofCBEtopOccu->Draw("COLZ");
format_pad(c_flat[nhours],4);
HTofCBEbotOccu->Draw("COLZ");

c_sides[nhours]->cd();
c_sides[nhours]->Divide(4,2);

format_pad(c_sides[nhours],1);
HTofCOR_XOccu->Draw("COLZ");
format_pad(c_sides[nhours],2);
HTofCOR_min_XOccu->Draw("COLZ");

format_pad(c_sides[nhours],3);
HTofCOR_YOccu->Draw("COLZ");
format_pad(c_sides[nhours],4);
HTofCOR_min_YOccu->Draw("COLZ");
format_pad(c_sides[nhours],5);
HTofCBE_XOccu->Draw("COLZ");
format_pad(c_sides[nhours],6);
HTofCBE_min_XOccu->Draw("COLZ");
format_pad(c_sides[nhours],7);
HTofCBE_YOccu->Draw("COLZ");
format_pad(c_sides[nhours],8);
HTofCBE_min_YOccu->Draw("COLZ");

//Last histogram in the pdf
c_flat[nhours]->Print(Form("%s",oo.Data())); c_sides[nhours]->Print(Form("%s",oo.Data())); c_flat[nhours]->Print(Form("%s]",oo.Data()));

myfile.close();

cout << endl << "I am done" << endl;
return 1;

}
