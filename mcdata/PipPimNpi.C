#define PipPimNpi_cxx
#include "PipPimNpi.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TF1.h"
#include "TGraphErrors.h"

void PipPimNpi::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L PipPimNpi.C
//      Root > PipPimNpi t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

	TCanvas *c1 = new TCanvas();

   Long64_t nentries = fChain->GetEntriesFast();
/*
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
   }
*/
	
	TFile *f = new TFile("presolution.root","recreate");

	const int Npart = 40;
	double pcut[Npart+1];
	for (int i=0; i<Npart+1; i++){
	  pcut[i] = 0.0 + (2.0-0.0)/Npart*i;
	}
	TH1D *hpsigma = new TH1D("hpsigma","hpsigma",200,-0.05,0.05);
	gStyle->SetOptFit(1111);
	gStyle->SetFitFormat("5.6g");
	int Nevt;
	int np=0;
	double pmean[Npart];
	double pmeer[Npart];
	double psigma[Npart];
	double psiger[Npart];

	char cuts[1000];
	for (int i=0; i<Npart; i++){
		sprintf(cuts,"mc_pipp>%f & mc_pipp<%f",pcut[i],pcut[i+1]);
		hpsigma->Reset();
   		Nevt = fChain->Draw("(pipp-mc_pipp)/mc_pipp>>hpsigma",cuts);
		if (Nevt < 100 ) continue;
		hpsigma->Fit("gaus");
		sprintf(cuts,"p_in_cuts_%02d",i);
		hpsigma->SetTitle(cuts);
		hpsigma->Write();

		pmean[np]  = (pcut[i+1]+pcut[i])/2.;
		pmeer[np]  = (pcut[i+1]-pcut[i])/2.;
		psigma[np] = hpsigma->GetFunction("gaus")->GetParameter(2);
		psiger[np] = hpsigma->GetFunction("gaus")->GetParError(2);
		np ++;
	} 

	TGraphErrors *graph = new TGraphErrors(np,pmean,psigma,pmeer,psiger);
	graph->SetMarkerStyle(5);
	graph->SetMarkerColor(2);
	graph->SetName("resolution");

	graph->Write();
	f->Close();

}
