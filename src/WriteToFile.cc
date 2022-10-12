#include "WriteToFile.hh"

std::once_flag flag;

TFile* rFile;
TTree* tGlauber;
TTree* tRun;
TTree* tClusters;
TTree* tFermiMom;
TH1D*  histo[20];
TH2D*  histo2[10];

int calls = 0;

AAMCCEvent event;

void InitFile(AAMCCrun run){
    std::string fileName;
    if ( run.fileName.empty()) {fileName = "GRATE_"+run.SysA+run.SysB+"_"+std::to_string(run.KinEnPerNucl/GeV)+"_GeV_"+std::to_string(run.iterations)+"_events";}
    else {fileName = run.fileName;}
    std::string  fileType = "root";
    std::string fileFullName = fileName+"."+fileType;
    Int_t compressionFactor = 9;
    rFile = new TFile(fileFullName.c_str(), "RECREATE", fileName.c_str(),compressionFactor);
}

void InitTrees(AAMCCrun run){
    tGlauber = new TTree("Glauber","Events from glauber modeling");
    tRun = new TTree("Conditions","preconditions for modeling");
    tClusters = new TTree("MST-Clusters","TTree to store clusters");
    tFermiMom = new TTree("FermiMomentum", "Fermi momentum");

    tGlauber->Branch("id", &event.id, "id/i");
    tGlauber->Branch("A_on_A", "std::vector" ,&event.MassOnSideA);
    tGlauber->Branch("A_on_B", "std::vector" ,&event.MassOnSideB);
    tGlauber->Branch("Z_on_A", "std::vector" ,&event.ChargeOnSideA);
    tGlauber->Branch("Z_on_B", "std::vector" ,&event.ChargeOnSideB);
    tGlauber->Branch("Nhard", &event.Nhard, "Nhard/I");
    tGlauber->Branch("Ncoll", &event.Ncoll, "Ncoll/I");
    tGlauber->Branch("Ncollpp", &event.Ncollpp, "Ncollpp/I");
    tGlauber->Branch("Ncollpn", &event.Ncollpn, "Ncollpn/I");
    tGlauber->Branch("Ncollnn", &event.Ncollnn, "Ncollnn/I");
    tGlauber->Branch("Npart", &event.Npart, "Npart/I");
    tGlauber->Branch("NpartA", &event.NpartA, "NpartA/I");
    tGlauber->Branch("NpartB", &event.NpartB, "NpartB/I");
    if(true /*writePseudorapidity*/){
        tGlauber->Branch("pseudorapidity_on_A", "std::vector", &event.pseudorapidity_A);
        tGlauber->Branch("pseudorapidity_on_B", "std::vector", &event.pseudorapidity_B);
    }
    if(true/*WriteMomentum*/){
        tGlauber->Branch("pX_on_A", "std::vector" ,&event.pXonSideA,128000,1);
        tGlauber->Branch("pY_on_A", "std::vector" ,&event.pYonSideA,128000,1);
        tGlauber->Branch("pZ_on_A", "std::vector" ,&event.pZonSideA,128000,1);
        tGlauber->Branch("pX_on_B", "std::vector" ,&event.pXonSideB,128000,1);
        tGlauber->Branch("pY_on_B", "std::vector" ,&event.pYonSideB,128000,1);
        tGlauber->Branch("pZ_on_B", "std::vector" ,&event.pZonSideB,128000,1);
    }
    tGlauber->Branch("impact_parameter", &event.b, "impact_parameter/f");
    tGlauber->Branch("PhiRotA", &event.PhiRotA, "PhiRotA/f");
    tGlauber->Branch("ThetaRotA", &event.ThetaRotA, "ThetaRotA/f");
    tGlauber->Branch("PhiRotB", &event.PhiRotB, "PhiRotB/f");
    tGlauber->Branch("ThetaRotB", &event.ThetaRotB, "ThetaRotB/f");
    tGlauber->Branch("Ecc", &event.Ecc, "Ecc[10]/f");
    tGlauber->Branch("id", &event.id, "id/i");
    tGlauber->Branch("A_on_A", "std::vector" ,&event.MassOnSideA);
    tGlauber->Branch("A_on_B", "std::vector" ,&event.MassOnSideB);
    tGlauber->Branch("Z_on_A", "std::vector" ,&event.ChargeOnSideA);
    tGlauber->Branch("Z_on_B", "std::vector" ,&event.ChargeOnSideB);
    tGlauber->Branch("Nhard", &event.Nhard, "Nhard/I");
    tGlauber->Branch("Ncoll", &event.Ncoll, "Ncoll/I");
    tGlauber->Branch("Ncollpp", &event.Ncollpp, "Ncollpp/I");
    tGlauber->Branch("Ncollpn", &event.Ncollpn, "Ncollpn/I");
    tGlauber->Branch("Ncollnn", &event.Ncollnn, "Ncollnn/I");
    tGlauber->Branch("Npart", &event.Npart, "Npart/I");
    tGlauber->Branch("NpartA", &event.NpartA, "NpartA/I");
    tGlauber->Branch("NpartB", &event.NpartB, "NpartB/I");

    tClusters->Branch("Aa_cl", "std::vector" ,&event.A_cl);
    tClusters->Branch("Za_cl", "std::vector" ,&event.Z_cl);
    tClusters->Branch("d", &event.d_MstA ,"d/d");
    tClusters->Branch("Clust_num_a", &event.ClustNumA ,"Clust_num/I");
    tClusters->Branch("Ab_cl", "std::vector" ,&event.Ab_cl);
    tClusters->Branch("Zb_cl", "std::vector" ,&event.Zb_cl);
    tClusters->Branch("d_b", &event.d_MstB ,"d/d");
    tClusters->Branch("Clust_num_b", &event.ClustNumB ,"Clust_num_b/I");

    tClusters->Branch("Aa_cl", "std::vector" ,&event.A_cl);
    tClusters->Branch("Za_cl", "std::vector" ,&event.Z_cl);
    tClusters->Branch("d", &event.d_MstA ,"d/d");
    tClusters->Branch("Clust_num_a", &event.ClustNumA ,"Clust_num/I");
    tClusters->Branch("Ab_cl", "std::vector" ,&event.Ab_cl);
    tClusters->Branch("Zb_cl", "std::vector" ,&event.Zb_cl);
    tClusters->Branch("d_b", &event.d_MstB ,"d/d");
    tClusters->Branch("Clust_num_b", &event.ClustNumB ,"Clust_num_b/I");


    tGlauber->Branch("impact_parameter", &event.b, "impact_parameter/f");

    tGlauber->Branch("PhiRotA", &event.PhiRotA, "PhiRotA/f");
    tGlauber->Branch("ThetaRotA", &event.ThetaRotA, "ThetaRotA/f");
    tGlauber->Branch("PhiRotB", &event.PhiRotB, "PhiRotB/f");
    tGlauber->Branch("ThetaRotB", &event.ThetaRotB, "ThetaRotB/f");
    tGlauber->Branch("Ecc", &event.Ecc, "Ecc[10]/f");

    tGlauber->Branch("Ex_En_per_nucleon", &event.ExEnA, "Ex_En_per_nucleon/f");

    tFermiMom->Branch("Fermi_momentum_x_side_A", &event.FermiMomA_x, "Fermi_momentumA_x/d");
    tFermiMom->Branch("Fermi_momentum_y_side_A", &event.FermiMomA_y, "Fermi_momentumA_y/d");
    tFermiMom->Branch("Fermi_momentum_z_side_A", &event.FermiMomA_z, "Fermi_momentumA_y/d");
    tFermiMom->Branch("Fermi_momentum_x_side_B", &event.FermiMomB_x, "Fermi_momentumB_x/d");
    tFermiMom->Branch("Fermi_momentum_y_side_B", &event.FermiMomB_y, "Fermi_momentumB_y/d");
    tFermiMom->Branch("Fermi_momentum_z_side_B", &event.FermiMomB_z, "Fermi_momentumB_y/d");

}

void InitRunTree(AAMCCrun run){
    //needs a realisation
}

void InitHisto(AAMCCrun run){
    histo[0] =  new TH1D("Charge distruibution for side B"," ;Z;entries",run.ZinitB+1,-0.5, run.ZinitB+0.5);
    histo[1] =  new TH1D("M distr"," ;M;entries",100, -0.5, 100+0.5);
    //histo[2] =  new TH1D("pz for neutrons, side A",";pz;",1e+3,InCond->GetPzA()/MeV/InCond->GetSourceA() - 20e+3, InCond->GetPzA()/MeV/InCond->GetSourceA() + 20e+3);
    //histo[3] =  new TH1D("pz for protons, side A"," ;pz;",1e+3, InCond->GetPzA()/MeV/InCond->GetSourceA() - 20e+3, InCond->GetPzA()/MeV/InCond->GetSourceA() + 20e+3);
    //histo[4] =  new TH1D("pz (2 < Z < 20), side A"," ;pz;",1e+3, InCond->GetPzA()/MeV/20 - 500e+3, InCond->GetPzA()/MeV/20 + 800e+3);
    //histo[5] =  new TH1D("pz (Z > 20), side A"," ;pz;",1e+3, InCond->GetPzA()/MeV/2 - 1.5e+6, InCond->GetPzA()/MeV/2 + 2e+6);
    histo[6] =  new TH1D("Charge distruibution"," ;Z;entries",run.ZinitA+1,-0.5, run.ZinitA+0.5);
    histo[7] =  new TH1D("Mass distribution", " ;A,entries",run.AinitA,0.5, run.AinitA+0.5);
    double CritDist = 5;
    histo2[1] = new TH2D("Ex En distribution, side A"," ;E* / A;A_{pf}/A",300, 0, 15, run.AinitA+1, 0, 1);
    histo2[2] = new TH2D("Mass and Charge distribution"," ;Z;A",run.ZinitA+1, -0.5, run.ZinitA+0.5, run.AinitA+1, -0.5, run.AinitA+0.5);
    histo2[3] = new TH2D("px vs py for neutrons", ";px;py", 100,-200,200,100,-200,200);
    histo2[4] = new TH2D("px vs py for protons", ";px;py", 100,-200,200,100,-200,200);
    histo2[5] = new TH2D("px vs py for IMF", ";px;py", 100,-200,200,100,-200,200);
    histo2[6] = new TH2D("px vs py for heavy fragments", ";px;py", 100,-200,200,100,-200,200);
    histo2[7] = new TH2D("Ex En VS d, side A"," ;E* / A;d",300, 0, 15, 280, 0, CritDist + 0.5);
    histo2[8] = new TH2D("fermi px vs fermi py for deltaA = 1", ";px;py", 50,-200,200,100,-200,200);
    int Rmax = 20;
    int Num_ent = 1200;
    histo[8] = new TH1D ("Neutron distribution A", ";R;entries", Num_ent, 0, Rmax);
    histo[9] = new TH1D ("Proton distribution A", ";R;entries", Num_ent, 0, Rmax);
    histo[10] = new TH1D("Neutron distribution B", ";R;entries", Num_ent, 0, Rmax);
    histo[11] = new TH1D("Proton distribution B", ";R;entries", Num_ent, 0, Rmax);
}

void Initialize(AAMCCrun run){
    InitFile(run);
    InitTrees(run);
    InitHisto(run);
    InitRunTree(run);
};

void FillTrees(AAMCCEvent ev){
    event = ev;
    tGlauber->Fill();
    tClusters->Fill();
    tFermiMom->Fill();
    //tRun->Fill(); //should be activated after reslisation of InitRunTree();
}

void FillHisto(AAMCCEvent ev, NucleonVector nucleons){
    //needs a realisation
}

void WriteToFile(AAMCCEvent* ev, AAMCCrun run, NucleonVector* nucleons){
    std::call_once(flag,Initialize, run);
    FillTrees((*ev));
    FillHisto((*ev), (*nucleons));
    calls++;
    if(calls == run.iterations) {rFile->Write();
    G4cout << "\n----> Data were written into the file " << run.fileName+".root" << G4endl;}
};