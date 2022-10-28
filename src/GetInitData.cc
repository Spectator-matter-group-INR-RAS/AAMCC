#include "GetInitData.hh"


AAMCCrun getTheRunData(TString inputFileName="test.root")
{
	AAMCCrun runData;
	TFile *fIn = new TFile(inputFileName, "read");
	URun *run = (URun *)fIn->Get("run");

	// TTree *tree = (TTree *)fIn->Get("events");
	// EventInitialState *iniState{nullptr};
	// tree->SetBranchAddress("iniState", &iniState);

	runData.AinitA = run->GetAProj();
	std::cout << runData.AinitA << endl;
	runData.ZinitA = run->GetZProj();
	std::cout << runData.ZinitA << endl;
	runData.AinitB = run->GetATarg();
	runData.ZinitB = run->GetZTarg();

	runData.pzB = run->GetPProj();
	runData.pzA = run->GetPTarg();

	runData.SqrtSnn = run->GetNNSqrtS();
	std::cout << runData.SqrtSnn << endl;
	runData.XsectTot = run->GetSigma();

	runData.iterations = run->GetNEvents();
	std::cout << runData.iterations << endl;
	

	runData.isCollider = true;

	runData.isQMD = true;

	return runData;
}