#include "MCiniWriter.hh"

MCINIWriter::MCINIWriter() {
    InitTree();
}

MCINIWriter::~MCINIWriter() noexcept {
    std::string fileName;
    if ( runData.fileName.empty()) {fileName = "GRATE_"+runData.SysA+runData.SysB+"_"+std::to_string(runData.KinEnPerNucl/GeV)+"_GeV_"+std::to_string(runData.iterations)+"_events";}
    else {fileName = runData.fileName;}
    std::string  fileType = "root";
    std::string fileFullName = fileName+"_mcini_"+"."+fileType;
    std::unique_ptr<TFile> file( TFile::Open(fileFullName.c_str(), "RECREATE", fileName.c_str(), 9) );

    file->WriteObject(tMCini.get(), "events");
    file->WriteObject(GetURun(tMCini.get(), runData).get(), "run");
    G4cout << "\n----> Data were written into the file in MCini format " << fileFullName<< G4endl;
}

void MCINIWriter::operator()(AAMCCEvent *ev, AAMCCrun *run, aamcc::NucleonVector *nucleons) {
    if(run != nullptr) runData = (*run);
    if(ev != nullptr && nucleons != nullptr){ FillTree(ev);}
}
