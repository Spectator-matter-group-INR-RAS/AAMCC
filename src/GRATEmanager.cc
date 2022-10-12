#include "GRATEmanager.hh"

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h" 
#include "TMath.h"


GRATEmanager::GRATEmanager()
  : sourceZ(0), sourceA(0), KinEn(-1.), SqrtSnn(-1.), XsectNN(-1.), lowLimitExEn( 0.), upperLimitExEn( 100.), binsExEn(1), eventsPerBin(1), StatisticsLabel(-1), iterations(1), wM(0), wP(0), NucleusInputLabel(0), IsCollider(0), upperLimitB(-1), CritDist(2.7), angle(0), DeExModel("")
{  
  std::cout << "######### Abrasion-Ablation model using Glauber Monte Carlo and Geant4" <<std::endl;
  while(!NucleusInputLabel){
    std::cout << "Please enter colliding nucleus name (side A). U, U2, Pb, Pbrw, Pbpn, Pbpnrw, Au, Aurw, Au2, Au2rw, Xe, Ag, Br, Cu, Ca2 (with SRC), Ar, Al, O, O2 (with SRC), Oho (for HO param.), C, He4, He3, H3, d, p is available : ";
    std::cin >> SysA;
    NucleusInputLabel = InCond->SetSysA(SysA);
    sourceA = InCond->GetSourceA();
    sourceZ = InCond->GetSourceZ();
  }


  NucleusInputLabel = 0;

  while(!NucleusInputLabel){
    std::cout << "Please enter colliding nucleus name (side B). U, U2, Pb, Pbrw, Pbpn, Pbpnrw, Au, Aurw, Au2, Au2rw, Xe, Ag, Br, Cu, Ca2 (with SRC), Ar, Al, O, O2 (with SRC), Oho (for HO param.), C, He4, He3, H3, d, p is available : ";
    std::cin >> SysB;
    NucleusInputLabel = InCond->SetSysB(SysB);
    sourceAb = InCond->GetSourceAb();
    sourceZb = InCond->GetSourceZb();
  }

  std::cout<<"Input lower limit for impact parameter in fm (MB if negative) : ";
  std::cin >> lowLimitB;

  if(lowLimitB > -0.000001){
    std::cout<<"Input upper limit for impact parameter in fm (MB if negative) : ";
    std::cin >> upperLimitB;
  }
  else if(0){
    lowLimitB = 0;
    upperLimitB = 20;
  }

  std::cout<<"Do you want to calculate collisions for collider or for fixed target geometry (1 for collider, 0 for fixed target) : ";
  std::cin >> IsCollider;

  while ( KinEn<280.0*MeV/GeV && SqrtSnn <0.) {
    if (!IsCollider) {
      std::cout << "Please enter kinetic enegy of projectile nucleus (per nucleon in GeV) : ";
      std::cin >> KinEn;
      if(KinEn<280.0*MeV/GeV){
        std::cout << "AAMCC works at kinetic energies above 280A MeV, please input higher energy!" << std::endl;
      }
    }
    else {
      std::cout << "Please enter s^1/2 of colliding nuclei (per nucleon in GeV) : ";
      std::cin >> SqrtSnn;
      G4double NuclMass = nucleonAverMass/CLHEP::GeV;
      G4double CollKinCheck = (SqrtSnn/2.0 - NuclMass);
      G4double KinEnAtFixTargetCheck = (2.0*(CollKinCheck + NuclMass)*(CollKinCheck + NuclMass)/(NuclMass*NuclMass) - 1.0)*NuclMass;
      if(KinEnAtFixTargetCheck < 280.0*MeV/GeV){
        std::cout << "AAMCC works at kinetic energies above 280A MeV, please input higher energy!" << std::endl;
        SqrtSnn = -1.0;
      }
    }
  }

  InCond->SetCollider(IsCollider);
  if(IsCollider){
    InCond->SetKinematics(SqrtSnn);
    SqrtSnn = InCond->GetSqrtSnn();
    KinEn = InCond->GetKinEnergy();
  }
  else{
    InCond->SetKinematics(KinEn);
    SqrtSnn = InCond->GetSqrtSnn();
    KinEn = InCond->GetKinEnergy();
  }

  sourceA=InCond->GetSourceA();
  sourceAb=InCond->GetSourceAb();
  lowLimitExEnB = lowLimitExEn;
  upperLimitExEnB = upperLimitExEn;
  lowLimitExEn *= sourceA;
  lowLimitExEn *= MeV;
  upperLimitExEn *=sourceA;
  upperLimitExEn *=MeV;
  lowLimitExEnB *= sourceAb;
  lowLimitExEnB *= MeV;
  upperLimitExEnB *=sourceAb;
  upperLimitExEnB *=MeV;

  while ( (StatisticsLabel<0) || (StatisticsLabel>7) || (upperLimitExEn<lowLimitExEn) ) {
    std::cout << "Please choose the level density function to be used: 1 - Ericson, 2 - Gaimard-Schmidt, 3 - ALADIN parametrization, 4 - Hybrid of 1 and 3, 7 - Fit of Hybrid : ";
    std::cin >> StatisticsLabel;
  }

  while(CritDist < 0) {
    std::cout<<"Please enter critical distance (in fm) : ";
    std::cin >> CritDist;
  }

  while(DeExModel.empty()){
      std::cout<<"Choose a model for fragment deexcitation. G4, ABLAXX, AAMCC or MIX (random mix of G4, ABLAXX and AAMCC) options are available: ";
      std::cin>>DeExModel;
  }

  std::cout<<"Write coordinates of nucleons in the text file or not (one event)? (1 - yes, 0 - no): ";
  std::cin >> InFileOrNot;


  if(!InFileOrNot){
    std::cout<<"Write momentum of each fragment?  (1 - yes, 0 - no) ";
    std::cin>>wM;

    std::cout<<"Write pseudorapidity of each fragments?  (1 - yes, 0 - no) ";
    std::cin>>wP;
  }


  while ( ((iterations<2) || (iterations>10000000)) && !InFileOrNot ) {
    std::cout<<"Please enter number of events to be generated: ";
    std::cin >> iterations;
  }

  std::cout << "Please enter the file name to write histograms (.root will be supplied): ";
  std::cin >> fileName; 

  runData.ZinitA = sourceZ; runData.ZinitB = sourceZb; runData.AinitA = sourceA; runData.AinitB = sourceAb; runData.SysA = SysA; runData.SysB = SysB; runData.fileName = fileName;
  runData.KinEnPerNucl = KinEn; runData.isCollider = IsCollider; runData.iterations = iterations;

  //for (G4int j=0; j<20; j++) histo[j] = 0;
  //for (G4int l=0; l<10; l++) histo2[l] = 0;
}


GRATEmanager::~GRATEmanager()
{
}

void GRATEmanager::BookHisto()
{
  // Open a file to keep histograms inside it
  if ( fileName.empty()) fileName = "GRATE_"+SysA+SysB+"_"+std::to_string(KinEn/GeV)+"_GeV_"+std::to_string(iterations)+"_events";
  fileType = "root";
  fileFullName = fileName+"."+fileType;
  compressionFactor = 9;
  fFile = new TFile(fileFullName, "RECREATE", fileName, compressionFactor);
  //Creating a Trees
  Glauber = new TTree("Glauber","Events from glauber modeling");
  modelingCo = new TTree("Conditions","preconditions for modeling");
  Clusters = new TTree("MST-Clusters","TTree to store clusters");
  FermiMom = new TTree("FermiMomentum", "Fermi momentum");
  // Book all histograms there ...
  histo[0] =  new TH1D("Charge distruibution for side B"," ;Z;entries",sourceZb+1,-0.5, sourceZb+0.5); 

  histo[1] =  new TH1D("M distr"," ;M;entries",100, -0.5, 100+0.5);

  histo[2] =  new TH1D("pz for neutrons, side A",";pz;",1e+3,InCond->GetPzA()/MeV/InCond->GetSourceA() - 20e+3, InCond->GetPzA()/MeV/InCond->GetSourceA() + 20e+3);
  histo[3] =  new TH1D("pz for protons, side A"," ;pz;",1e+3, InCond->GetPzA()/MeV/InCond->GetSourceA() - 20e+3, InCond->GetPzA()/MeV/InCond->GetSourceA() + 20e+3);
  histo[4] =  new TH1D("pz (2 < Z < 20), side A"," ;pz;",1e+3, InCond->GetPzA()/MeV/20 - 500e+3, InCond->GetPzA()/MeV/20 + 800e+3);
  histo[5] =  new TH1D("pz (Z > 20), side A"," ;pz;",1e+3, InCond->GetPzA()/MeV/2 - 1.5e+6, InCond->GetPzA()/MeV/2 + 2e+6);


  histo[6] =  new TH1D("Charge distruibution"," ;Z;entries",sourceZ+1,-0.5, sourceZ+0.5);

  histo[7] =  new TH1D("Mass distribution", " ;A,entries",sourceA,0.5, sourceA+0.5);

  histo2[1] = new TH2D("Ex En distribution, side A"," ;E*/A;A_{pf}/A",300, 0, 15, sourceA+1, 0, 1);

  histo2[2] = new TH2D("Mass and Charge distribution"," ;Z;A",sourceZ+1, -0.5, sourceZ+0.5, sourceA+1, -0.5, sourceA+0.5);

  histo2[3] = new TH2D("px vs py for neutrons", ";px;py", 100,-200,200,100,-200,200);
  histo2[4] = new TH2D("px vs py for protons", ";px;py", 100,-200,200,100,-200,200);
  histo2[5] = new TH2D("px vs py for IMF", ";px;py", 100,-200,200,100,-200,200);
  histo2[6] = new TH2D("px vs py for heavy fragments", ";px;py", 100,-200,200,100,-200,200);
  histo2[7] = new TH2D("Ex En VS d, side A"," ;E*/A;d",300, 0, 15, 280, 0, CritDist + 0.5);
  histo2[8] = new TH2D("fermi px vs fermi py for deltaA = 1", ";px;py", 50,-200,200,100,-200,200);

  G4int Rmax = 20;
  G4int Num_ent = 1200;
  histo[8] = new TH1D ("Neutron distribution A", ";R;entries", Num_ent, 0, Rmax);
  histo[9] = new TH1D ("Proton distribution A", ";R;entries", Num_ent, 0, Rmax);
  histo[10] = new TH1D("Neutron distribution B", ";R;entries", Num_ent, 0, Rmax);
  histo[11] = new TH1D("Proton distribution B", ";R;entries", Num_ent, 0, Rmax);

  this->InitTree();

  G4cout << "Histograms will be written to " << fileFullName << G4endl;
}


void GRATEmanager::CalcXsectNN()
{
  G4double shadowing = 41.5/70; //according to Eskola K.J. et al. PHYSICAL REVIEW LETTERS 125, 212301 (2020)
  G4double KinEnAtFixTarget = 0;
  if(IsCollider){KinEnAtFixTarget = (2.0*(KinEn + nucleonAverMass*G4double(sourceA))*(KinEn + nucleonAverMass*G4double(sourceA))/(nucleonAverMass*G4double(sourceA)*nucleonAverMass*G4double(sourceA)) - 1.0)*nucleonAverMass*G4double(sourceA);}
  else{KinEnAtFixTarget = KinEn;}

  if(KinEnAtFixTarget/G4double(sourceA) < 425*GeV){
    G4double Tkin[2] = {0};
    G4double xsect[2] = {-1};
    std::string filepath(__FILE__);
    std::string filename(basename(__FILE__));
    filepath.erase(filepath.length() - filename.length(), filename.length());
    filepath += "bystricky.dat";
    XsectFile.open(filepath.c_str());
    while (Tkin[0]*GeV < KinEnAtFixTarget/G4double(sourceA)) {
      Tkin[0] = Tkin[1];	
      xsect[0] = xsect[1];
      XsectFile >> Tkin[1] >> xsect[1];
      if (!XsectFile.good()) break;
    }
    G4double a = (xsect[1]-xsect[0])/(Tkin[1]-Tkin[0]);
    G4double b = xsect[1] - a*Tkin[1];
    XsectNN = a*KinEnAtFixTarget/(G4double(sourceA)*GeV)+b;
  }
  else{
    G4double S = SqrtSnn*SqrtSnn;
    XsectNN = 25.0+0.146*pow(log(S/(GeV*GeV)),2);
  }
  //XsectNN *= shadowing;
}

void GRATEmanager::CalcNucleonDensity(TObjArray* nucleons_pre, G4double b)
{
  G4float X_pre, Y_pre, Z_pre, R_pre;
  for (G4int iter = 0; iter < nucleons_pre->GetEntries(); iter++) {
    TGlauNucleon* nucleon_pre = (TGlauNucleon*)(nucleons_pre->At(iter));
    if (nucleon_pre->IsInNucleusA()){X_pre = nucleon_pre->GetX() + b/2;}
    if (nucleon_pre->IsInNucleusB()){X_pre = nucleon_pre->GetX() - b/2;}
    Y_pre = nucleon_pre->GetY();
    Z_pre = nucleon_pre->GetZ();
    R_pre = std::sqrt(pow(X_pre,2) + pow(Y_pre,2) + pow(Z_pre,2));
    if (nucleon_pre->IsProton() && nucleon_pre->IsInNucleusA()) { (*this).GetHisto(9)->Fill(R_pre); }
    if (nucleon_pre->IsNeutron() && nucleon_pre->IsInNucleusA()) { (*this).GetHisto(8)->Fill(R_pre); }
    if (nucleon_pre->IsProton() && nucleon_pre->IsInNucleusB()) { (*this).GetHisto(11)->Fill(R_pre); }
    if (nucleon_pre->IsNeutron() && nucleon_pre->IsInNucleusB()) { (*this).GetHisto(10)->Fill(R_pre); }
  }
}

void GRATEmanager::CleanHisto()
{
  fFile->Write();
  G4cout << "\n----> Histograms were written into the file " << fileFullName << G4endl;
  delete fFile;
}

void GRATEmanager::FillConditionsTree(G4double Xsect)
{
  G4double XsectTot = 0;
  G4double KineticEnergy = 0;
  G4double pzA = 0;
  G4double pzB = 0;
  G4double Mass_on_A = 0;
  G4double Mass_on_B = 0;
  G4double Charge_on_A = 0;
  G4double Charge_on_B = 0;

  modelingCo->Branch("Xsect_total", &XsectTot,"Xsect_total/d");
  modelingCo->Branch("Kinetic_energy_per_nucleon_of_projectile_in_GeV", &KineticEnergy,"Kinetic_energy_of_per_nucleon_projectile_in_GeV/d");
  modelingCo->Branch("SqrtS_nn_in_GeV", &SqrtSnn,"SqrtS_nn_in_GeV/d");
  modelingCo->Branch("pZ_in_MeV_on_A", &pzA,"pZ_in_MeV_on_A/d");
  modelingCo->Branch("pZ_in_MeV_on_B", &pzB,"pZ_in_MeV_on_B/d");
  modelingCo->Branch("Mass_on_A", &Mass_on_A,"Mass_on_A/d");
  modelingCo->Branch("Mass_on_B", &Mass_on_B,"Mass_on_B/d");
  modelingCo->Branch("Charge_on_A", &Charge_on_A,"Charge_on_A/d");
  modelingCo->Branch("Charge_on_B", &Charge_on_B,"Charge_on_B/d");

  XsectTot = Xsect;
  SqrtSnn = InCond->GetSqrtSnn()/GeV;
  KineticEnergy = InCond->GetKinEnergy()/(sourceA*GeV);
  pzA = InCond->GetPzA()/MeV;
  pzB = InCond->GetPzB()/MeV;
  Mass_on_A = InCond->GetSourceA();
  Mass_on_B = InCond->GetSourceAb();
  Charge_on_A = InCond->GetSourceZ();
  Charge_on_B = InCond->GetSourceZb();

  modelingCo->Fill();
}

void GRATEmanager::WriteNucleonsCoordinatesInFile(GMSTClusterVector clusters_to_excit_A, GMSTClusterVector clusters_to_excit_B, G4double b)
{
  ofstream Event("Event.txt");

  Event<<"Clusters on side A: \n";

  for(G4int i = 0; i < clusters_to_excit_A.size(); ++i) {
    Event<<"Nucleons from "<<i+1<<"-th cluster\n";
    for(G4int iCoord = 0; iCoord<(clusters_to_excit_A[i]).GetCoordinates().size(); iCoord++){
      Event<<"X = "<<(clusters_to_excit_A[i]).GetCoordinates().at(iCoord).X()<<" Y = "<<(clusters_to_excit_A[i]).GetCoordinates().at(iCoord).Y()<<" Z = "<<(clusters_to_excit_A[i]).GetCoordinates().at(iCoord).Z()<<"\n";
    }
  }

  Event<<"\n";  Event<<"Clusters on side B: \n";

  for(G4int i = 0; i < clusters_to_excit_B.size(); ++i) {
    Event<<"Nucleons from "<<i+1<<"-th cluster\n";
    for(G4int iCoord = 0; iCoord<(clusters_to_excit_B[i]).GetCoordinates().size(); iCoord++){
      Event<<"X = "<<(clusters_to_excit_B[i]).GetCoordinates().at(iCoord).X()<<" Y = "<<(clusters_to_excit_B[i]).GetCoordinates().at(iCoord).Y()<<" Z = "<<(clusters_to_excit_B[i]).GetCoordinates().at(iCoord).Z()<<"\n";
    }
  }
  Event<<"\nb = "<<b<<" fm \n";
}

void GRATEmanager::FillEventTree(AAMCCEvent *ev_in) {
    event = (*ev_in);
    this->GetTree()->Fill();
    this->GetTreeMST()->Fill();
    this->GetTreeFermiMom()->Fill();
}

void GRATEmanager::InitTree() {
    this->GetTree()->Branch("id", &event.id, "id/i");
    this->GetTree()->Branch("A_on_A", "std::vector" ,&event.MassOnSideA);
    this->GetTree()->Branch("A_on_B", "std::vector" ,&event.MassOnSideB);
    this->GetTree()->Branch("Z_on_A", "std::vector" ,&event.ChargeOnSideA);
    this->GetTree()->Branch("Z_on_B", "std::vector" ,&event.ChargeOnSideB);
    this->GetTree()->Branch("Nhard", &event.Nhard, "Nhard/I");
    this->GetTree()->Branch("Ncoll", &event.Ncoll, "Ncoll/I");
    this->GetTree()->Branch("Ncollpp", &event.Ncollpp, "Ncollpp/I");
    this->GetTree()->Branch("Ncollpn", &event.Ncollpn, "Ncollpn/I");
    this->GetTree()->Branch("Ncollnn", &event.Ncollnn, "Ncollnn/I");
    this->GetTree()->Branch("Npart", &event.Npart, "Npart/I");
    this->GetTree()->Branch("NpartA", &event.NpartA, "NpartA/I");
    this->GetTree()->Branch("NpartB", &event.NpartB, "NpartB/I");

    this->GetTreeMST()->Branch("Aa_cl", "std::vector" ,&event.A_cl);
    this->GetTreeMST()->Branch("Za_cl", "std::vector" ,&event.Z_cl);
    this->GetTreeMST()->Branch("d", &event.d_MstA ,"d/d");
    this->GetTreeMST()->Branch("Clust_num_a", &event.ClustNumA ,"Clust_num/I");
    this->GetTreeMST()->Branch("Ab_cl", "std::vector" ,&event.Ab_cl);
    this->GetTreeMST()->Branch("Zb_cl", "std::vector" ,&event.Zb_cl);
    this->GetTreeMST()->Branch("d_b", &event.d_MstB ,"d/d");
    this->GetTreeMST()->Branch("Clust_num_b", &event.ClustNumB ,"Clust_num_b/I");

    if(this->WritePseudorapidity()){
        this->GetTree()->Branch("pseudorapidity_on_A", "std::vector", &event.pseudorapidity_A);
        this->GetTree()->Branch("pseudorapidity_on_B", "std::vector", &event.pseudorapidity_B);
    }
    if(this->WriteMomentum()){
        this->GetTree()->Branch("pX_on_A", "std::vector" ,&event.pXonSideA,128000,1);
        this->GetTree()->Branch("pY_on_A", "std::vector" ,&event.pYonSideA,128000,1);
        this->GetTree()->Branch("pZ_on_A", "std::vector" ,&event.pZonSideA,128000,1);
        this->GetTree()->Branch("pX_on_B", "std::vector" ,&event.pXonSideB,128000,1);
        this->GetTree()->Branch("pY_on_B", "std::vector" ,&event.pYonSideB,128000,1);
        this->GetTree()->Branch("pZ_on_B", "std::vector" ,&event.pZonSideB,128000,1);
    }

    this->GetTree()->Branch("impact_parameter", &event.b, "impact_parameter/f");

    this->GetTree()->Branch("PhiRotA", &event.PhiRotA, "PhiRotA/f");
    this->GetTree()->Branch("ThetaRotA", &event.ThetaRotA, "ThetaRotA/f");
    this->GetTree()->Branch("PhiRotB", &event.PhiRotB, "PhiRotB/f");
    this->GetTree()->Branch("ThetaRotB", &event.ThetaRotB, "ThetaRotB/f");
    this->GetTree()->Branch("Ecc", &event.Ecc, "Ecc[10]/f");
    this->GetTree()->Branch("id", &event.id, "id/i");
    this->GetTree()->Branch("A_on_A", "std::vector" ,&event.MassOnSideA);
    this->GetTree()->Branch("A_on_B", "std::vector" ,&event.MassOnSideB);
    this->GetTree()->Branch("Z_on_A", "std::vector" ,&event.ChargeOnSideA);
    this->GetTree()->Branch("Z_on_B", "std::vector" ,&event.ChargeOnSideB);
    this->GetTree()->Branch("Nhard", &event.Nhard, "Nhard/I");
    this->GetTree()->Branch("Ncoll", &event.Ncoll, "Ncoll/I");
    this->GetTree()->Branch("Ncollpp", &event.Ncollpp, "Ncollpp/I");
    this->GetTree()->Branch("Ncollpn", &event.Ncollpn, "Ncollpn/I");
    this->GetTree()->Branch("Ncollnn", &event.Ncollnn, "Ncollnn/I");
    this->GetTree()->Branch("Npart", &event.Npart, "Npart/I");
    this->GetTree()->Branch("NpartA", &event.NpartA, "NpartA/I");
    this->GetTree()->Branch("NpartB", &event.NpartB, "NpartB/I");

    this->GetTreeMST()->Branch("Aa_cl", "std::vector" ,&event.A_cl);
    this->GetTreeMST()->Branch("Za_cl", "std::vector" ,&event.Z_cl);
    this->GetTreeMST()->Branch("d", &event.d_MstA ,"d/d");
    this->GetTreeMST()->Branch("Clust_num_a", &event.ClustNumA ,"Clust_num/I");
    this->GetTreeMST()->Branch("Ab_cl", "std::vector" ,&event.Ab_cl);
    this->GetTreeMST()->Branch("Zb_cl", "std::vector" ,&event.Zb_cl);
    this->GetTreeMST()->Branch("d_b", &event.d_MstB ,"d/d");
    this->GetTreeMST()->Branch("Clust_num_b", &event.ClustNumB ,"Clust_num_b/I");

    if(this->WritePseudorapidity()){
        this->GetTree()->Branch("pseudorapidity_on_A", "std::vector", &event.pseudorapidity_A);
        this->GetTree()->Branch("pseudorapidity_on_B", "std::vector", &event.pseudorapidity_B);
    }
    if(this->WriteMomentum()){
        this->GetTree()->Branch("pX_on_A", "std::vector" ,&event.pXonSideA,128000,1);
        this->GetTree()->Branch("pY_on_A", "std::vector" ,&event.pYonSideA,128000,1);
        this->GetTree()->Branch("pZ_on_A", "std::vector" ,&event.pZonSideA,128000,1);
        this->GetTree()->Branch("pX_on_B", "std::vector" ,&event.pXonSideB,128000,1);
        this->GetTree()->Branch("pY_on_B", "std::vector" ,&event.pYonSideB,128000,1);
        this->GetTree()->Branch("pZ_on_B", "std::vector" ,&event.pZonSideB,128000,1);
    }

    this->GetTree()->Branch("impact_parameter", &event.b, "impact_parameter/f");

    this->GetTree()->Branch("PhiRotA", &event.PhiRotA, "PhiRotA/f");
    this->GetTree()->Branch("ThetaRotA", &event.ThetaRotA, "ThetaRotA/f");
    this->GetTree()->Branch("PhiRotB", &event.PhiRotB, "PhiRotB/f");
    this->GetTree()->Branch("ThetaRotB", &event.ThetaRotB, "ThetaRotB/f");
    this->GetTree()->Branch("Ecc", &event.Ecc, "Ecc[10]/f");

    this->GetTree()->Branch("Ex_En_per_nucleon", &event.ExEnA, "Ex_En_per_nucleon/f");

    this->GetTreeFermiMom()->Branch("Fermi_momentum_x_side_A", &event.FermiMomA_x, "Fermi_momentumA_x/d");
    this->GetTreeFermiMom()->Branch("Fermi_momentum_y_side_A", &event.FermiMomA_y, "Fermi_momentumA_y/d");
    this->GetTreeFermiMom()->Branch("Fermi_momentum_z_side_A", &event.FermiMomA_z, "Fermi_momentumA_y/d");
    this->GetTreeFermiMom()->Branch("Fermi_momentum_x_side_B", &event.FermiMomB_x, "Fermi_momentumB_x/d");
    this->GetTreeFermiMom()->Branch("Fermi_momentum_y_side_B", &event.FermiMomB_y, "Fermi_momentumB_y/d");
    this->GetTreeFermiMom()->Branch("Fermi_momentum_z_side_B", &event.FermiMomB_z, "Fermi_momentumB_y/d");

    this->GetTree()->Branch("Ex_En_per_nucleon", &event.ExEnA, "Ex_En_per_nucleon/f");

    this->GetTreeFermiMom()->Branch("Fermi_momentum_x_side_A", &event.FermiMomA_x, "Fermi_momentumA_x/d");
    this->GetTreeFermiMom()->Branch("Fermi_momentum_y_side_A", &event.FermiMomA_y, "Fermi_momentumA_y/d");
    this->GetTreeFermiMom()->Branch("Fermi_momentum_z_side_A", &event.FermiMomA_z, "Fermi_momentumA_y/d");
    this->GetTreeFermiMom()->Branch("Fermi_momentum_x_side_B", &event.FermiMomB_x, "Fermi_momentumB_x/d");
    this->GetTreeFermiMom()->Branch("Fermi_momentum_y_side_B", &event.FermiMomB_y, "Fermi_momentumB_y/d");
    this->GetTreeFermiMom()->Branch("Fermi_momentum_z_side_B", &event.FermiMomB_z, "Fermi_momentumB_y/d");
}

void GRATEmanager::FillHisto(AAMCCEvent ev) {
}

void GRATEmanager::ToFile(AAMCCEvent* event, NucleonVector* nucleons, void (*toFile)(AAMCCEvent*, AAMCCrun, NucleonVector*)) {
 toFile(event, runData, nucleons);
}
