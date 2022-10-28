#include "GetInitData.hh"
// R__LOAD_LIBRARY(libMcIniData.so)

AAMCCrun getTheRunData(TString inputFileName="test.root")
{
	using namespace mcini;
	AAMCCrun runData;
	TFile *fIn = new TFile(inputFileName, "read");
	URun *run = (URun *)fIn->Get("run");

	// TTree *tree = (TTree *)fIn->Get("events");
	// EventInitialState *iniState{nullptr};
	// tree->SetBranchAddress("iniState", &iniState);

	runData.AinitA = run->GetAProj();
	std::cout << runData.AinitA;
	runData.ZinitA = run->GetZProj();
	runData.AinitB = run->GetATarg();
	runData.ZinitB = run->GetZTarg();

	runData.pzB = run->GetPProj();
	runData.pzA = run->GetPTarg();

	runData.SqrtSnn = run->GetNNSqrtS();
	runData.XsectTot = run->GetSigma();

	runData.iterations = run->GetNEvents();

	runData.isCollider = true;

	runData.isQMD = true;

	return runData;
}