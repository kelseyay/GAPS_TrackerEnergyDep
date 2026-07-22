using namespace std;

#include "KYtools.C"

#include "GDataEvent.hh"
#include "GGeometryTools.hh"
#include "GDataPoint.hh"
#include "GDataTrack.hh"
#include "GDataVertex.hh"

using namespace Crane::Analysis;
namespace ca = Crane::Analysis;
namespace cl = Crane::Common;
//using Crane::Calibration;

bool is_selected(const CEventRec* event){
    return true;
}

int main(int argc, char *argv[]){

GOptionParser* parser = GOptionParser::GetInstance();
parser->AddProgramDescription("Minimal Reproducable Example for Extracing Data from Reco Data");
parser->AddCommandLineOption<string>("in_path", "path to instrument data files", "./*", "i");
parser->AddCommandLineOption<bool>("save", "save the special events to a root file?",0,"s");
parser->AddCommandLineOption<string>("out_file", "name of output path", "", "o");
parser->AddCommandLineOption<string>("sts_root_name", "name of output root file", "sts_root_test.root", "n");
parser->ParseCommandLine(argc, argv);
parser->Parse();

bool SAVE = parser->GetOption<bool>("save");

string out_path = parser->GetOption<string>("out_file");
cout << "out path: " << out_path << endl;
if(out_path != "" && out_path[out_path.length()-1] != '/' ){ cout <<  "out path no slash! Adding! " << endl; out_path = out_path + '/'; }

string reco_path = parser->GetOption<string>("in_path");

string sts_file_name = parser->GetOption<string>("sts_root_name");
if(sts_file_name.compare(sts_file_name.length()-5,sts_file_name.length(),".root") != 0){ cout << "NO .root at the end of the root name! Adding!" << endl; sts_file_name = sts_file_name + ".root"; }

cout << reco_path << endl;
if(reco_path.compare(reco_path.length()-5,reco_path.length(),".root") == 0){ cout << ".root at the end of the reco path! Deleting!" << endl; reco_path = reco_path.substr(0,reco_path.length()-5); }

char FilenameRoot[400];
sprintf(FilenameRoot,"%s*.root",reco_path.c_str());
cout << FilenameRoot << endl;

TChain* events = new TChain("TreeRec");
CEventRec* reco_event = new CEventRec;
events->SetBranchAddress("Rec", &reco_event);
events->Add(FilenameRoot);

TTree* selected_events = new TTree("TreeRec","TreeRec");
selected_events->Branch("Rec",&reco_event);

for(uint i=0; i<10; i++){
//for(uint i=0; i<events->GetEntries(); i++){
    events->GetEntry(i);
    cout << "Event is " << i << endl;
    selected_events->Fill();
    /*
    if(is_selected(reco_event)){
        selected_events->Fill();
    }*/
}

TFile good_evts_file((out_path+ sts_file_name).c_str(), "recreate");
events->GetEntry(0);
cout << "source root file " << events->GetCurrentFile()->GetName() << endl;
//cout << "source root file Attempt Get GGeometry " << events->GetCurrentFile()->Get("GGeometry") << endl;
//cout << "source root file Attempt Get GOptions " << events->GetCurrentFile()->Get("GOptions") << endl;
//cout << "source root file Attempt Get nonsense " << events->GetCurrentFile()->Get("sdfaskjdf") << endl;

TObject* geo_tree = events->GetCurrentFile()->Get("GGeometry");
geo_tree->Write("GGeometry"); //This does work!
//TObject* goptions_tree = events->GetCurrentFile()->Get("GOptions"); //Well this doesn't work!!!
//goptions_tree->Write("GOptions"); //Why doesn't this work :') 2D viewer and everything else still works so maybe can move on.

selected_events->Write();
good_evts_file.Close();

cout << "I am done" << endl;

}
