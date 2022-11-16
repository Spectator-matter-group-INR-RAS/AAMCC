#include "AAMCCWriter.hh"

AAMCCWriter::AAMCCWriter() {
    InitTrees(); InitRunTree();
}

AAMCCWriter::~AAMCCWriter() {
    FillRunTree(runData);
    std::string fileName;
    if ( runData.fileName.empty()) {fileName = "GRATE_"+runData.SysA+runData.SysB+"_"+std::to_string(runData.KinEnPerNucl/GeV)+"_GeV_"+std::to_string(runData.iterations)+"_events";}
    else {fileName = runData.fileName;}
    std::string  fileType = "root";
    std::string fileFullName = fileName+"."+fileType;
    std::unique_ptr<TFile> file( TFile::Open(fileFullName.c_str(), "RECREATE", fileName.c_str(), 9) );

    file->WriteObject(tGlauber.get(), "Glauber");
    file->WriteObject(tFermiMom.get(), "FermiMomentum");
    file->WriteObject(tClusters.get(), "MST-Clusters");
    file->WriteObject(tRun.get(), "Conditions");
    G4cout << "\n----> Data were written into the file " << runData.fileName+".root" << G4endl;
}

void AAMCCWriter::operator()(AAMCCEvent *ev, AAMCCrun *run, aamcc::NucleonVector *nucleons) {
    if(run != nullptr) runData = (*run);
    if(callflag){InitHisto(runData); callflag = false;}
    if(ev != nullptr && nucleons != nullptr){FillTrees((*ev)); FillHisto((*ev), (*nucleons));}
}
