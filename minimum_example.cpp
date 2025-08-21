#include "GOptionParser.hh"
#include "CEventRec.hh"
#include "GGeometry.hh"
#include "TChain.h"
#include "TH1.h"
#include "TFile.h"

using std::string;

int main(int argc, char *argv[]){
    GOptionParser* parser = GOptionParser::GetInstance();
    parser->AddProgramDescription("Minimal Reproducable Example for Extracing Data from Reco Data");
    parser->AddCommandLineOption<string>("in_path", "path to instrument data files", "./*", "i");
    parser->AddCommandLineOption<string>("out_file", "name of output root file", "out.root", "o");
    parser->ParseCommandLine(argc, argv);
    parser->Parse();

    string reco_path = parser->GetOption<string>("in_path");
    string out_path = parser->GetOption<string>("out_file");

    TChain* Reco_Events = new TChain("TreeRec");
    CEventRec* Reco_Event = new CEventRec;
    Reco_Events->SetBranchAddress("Rec", &Reco_Event);

    Reco_Events->Add(reco_path.c_str());

    TH1D* hits_by_layer = new TH1D("layers", "Layer Numbers", 10, -0.5, 10.5);

    int vol_id;

    for(uint i=0; i<Reco_Events->GetEntries(); i++){
        Reco_Events->GetEntry(i);
        for(GRecoHit hit:Reco_Event->GetHitSeries()){
            vol_id = hit.GetVolumeId();
            if(GGeometryObject::IsTrackerVolume(vol_id)){
                hits_by_layer->Fill(GGeometryObject::GetTrackerLayer(vol_id));
            }
        }
    }

    TFile out_file(out_path.c_str(), "recreate");
    out_file.cd();
    hits_by_layer->Write();
    out_file.Close();

    return 1;
}
