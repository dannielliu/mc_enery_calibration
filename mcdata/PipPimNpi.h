//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Mar 23 21:33:40 2015 by ROOT version 5.34/24
// from TTree Jpsi3pi/ks N-Tuple example
// found on file: Jpsi_3pi.root
//////////////////////////////////////////////////////////

#ifndef PipPimNpi_h
#define PipPimNpi_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class PipPimNpi {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Double_t        rec_truth_mkkpi0;
   Int_t           run;
   Int_t           rec;
   Int_t           evttag;
   Int_t           indexmc;
   Int_t           pdgid[100];   //[indexmc]
   Int_t           motheridx[100];   //[indexmc]
   Int_t           Nchrg;
   Int_t           Nneu;
   Int_t           NCharge;
   Double_t        mdc_rvxy0;
   Double_t        mdc_rvz0;
   Double_t        mdc_rvthe0;
   Double_t        mdc_rvphi0;
   Double_t        dthe;
   Double_t        dphi;
   Double_t        dang;
   Double_t        eraw;
   Double_t        mc_pippx;
   Double_t        mc_pippy;
   Double_t        mc_pippz;
   Double_t        mc_pipp;
   Double_t        mc_pipe;
   Double_t        mc_pimpx;
   Double_t        mc_pimpy;
   Double_t        mc_pimpz;
   Double_t        mc_pimp;
   Double_t        mc_pime;
   Double_t        pippx;
   Double_t        pippy;
   Double_t        pippz;
   Double_t        pipp;
   Double_t        pipe;
   Double_t        pimpx;
   Double_t        pimpy;
   Double_t        pimpz;
   Double_t        pimp;
   Double_t        pime;
   Double_t        mpippim;

   // List of branches
   TBranch        *b_rec_truth_mkkpi0;   //!
   TBranch        *b_run;   //!
   TBranch        *b_rec;   //!
   TBranch        *b_evttag;   //!
   TBranch        *b_indexmc;   //!
   TBranch        *b_pdgid;   //!
   TBranch        *b_motheridx;   //!
   TBranch        *b_Nchrg;   //!
   TBranch        *b_Nneu;   //!
   TBranch        *b_NCharge;   //!
   TBranch        *b_mdc_rvxy0;   //!
   TBranch        *b_mdc_rvz0;   //!
   TBranch        *b_mdc_rvthe0;   //!
   TBranch        *b_mdc_rvphi0;   //!
   TBranch        *b_dthe;   //!
   TBranch        *b_dphi;   //!
   TBranch        *b_dang;   //!
   TBranch        *b_eraw;   //!
   TBranch        *b_mc_pippx;   //!
   TBranch        *b_mc_pippy;   //!
   TBranch        *b_mc_pippz;   //!
   TBranch        *b_mc_pipp;   //!
   TBranch        *b_mc_pipe;   //!
   TBranch        *b_mc_pimpx;   //!
   TBranch        *b_mc_pimpy;   //!
   TBranch        *b_mc_pimpz;   //!
   TBranch        *b_mc_pimp;   //!
   TBranch        *b_mc_pime;   //!
   TBranch        *b_pippx;   //!
   TBranch        *b_pippy;   //!
   TBranch        *b_pippz;   //!
   TBranch        *b_pipp;   //!
   TBranch        *b_pipe;   //!
   TBranch        *b_pimpx;   //!
   TBranch        *b_pimpy;   //!
   TBranch        *b_pimpz;   //!
   TBranch        *b_pimp;   //!
   TBranch        *b_pime;   //!
   TBranch        *b_mpippim;   //!

   PipPimNpi(TTree *tree=0);
   virtual ~PipPimNpi();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef PipPimNpi_cxx
PipPimNpi::PipPimNpi(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Jpsi_3pi.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("Jpsi_3pi.root");
      }
      f->GetObject("Jpsi3pi",tree);

   }
   Init(tree);
}

PipPimNpi::~PipPimNpi()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t PipPimNpi::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t PipPimNpi::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void PipPimNpi::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("rec_truth_mkkpi0", &rec_truth_mkkpi0, &b_rec_truth_mkkpi0);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("rec", &rec, &b_rec);
   fChain->SetBranchAddress("evttag", &evttag, &b_evttag);
   fChain->SetBranchAddress("indexmc", &indexmc, &b_indexmc);
   fChain->SetBranchAddress("pdgid", pdgid, &b_pdgid);
   fChain->SetBranchAddress("motheridx", motheridx, &b_motheridx);
   fChain->SetBranchAddress("Nchrg", &Nchrg, &b_Nchrg);
   fChain->SetBranchAddress("Nneu", &Nneu, &b_Nneu);
   fChain->SetBranchAddress("NCharge", &NCharge, &b_NCharge);
   fChain->SetBranchAddress("mdc_rvxy0", &mdc_rvxy0, &b_mdc_rvxy0);
   fChain->SetBranchAddress("mdc_rvz0", &mdc_rvz0, &b_mdc_rvz0);
   fChain->SetBranchAddress("mdc_rvthe0", &mdc_rvthe0, &b_mdc_rvthe0);
   fChain->SetBranchAddress("mdc_rvphi0", &mdc_rvphi0, &b_mdc_rvphi0);
   fChain->SetBranchAddress("dthe", &dthe, &b_dthe);
   fChain->SetBranchAddress("dphi", &dphi, &b_dphi);
   fChain->SetBranchAddress("dang", &dang, &b_dang);
   fChain->SetBranchAddress("eraw", &eraw, &b_eraw);
   fChain->SetBranchAddress("mc_pippx", &mc_pippx, &b_mc_pippx);
   fChain->SetBranchAddress("mc_pippy", &mc_pippy, &b_mc_pippy);
   fChain->SetBranchAddress("mc_pippz", &mc_pippz, &b_mc_pippz);
   fChain->SetBranchAddress("mc_pipp", &mc_pipp, &b_mc_pipp);
   fChain->SetBranchAddress("mc_pipe", &mc_pipe, &b_mc_pipe);
   fChain->SetBranchAddress("mc_pimpx", &mc_pimpx, &b_mc_pimpx);
   fChain->SetBranchAddress("mc_pimpy", &mc_pimpy, &b_mc_pimpy);
   fChain->SetBranchAddress("mc_pimpz", &mc_pimpz, &b_mc_pimpz);
   fChain->SetBranchAddress("mc_pimp", &mc_pimp, &b_mc_pimp);
   fChain->SetBranchAddress("mc_pime", &mc_pime, &b_mc_pime);
   fChain->SetBranchAddress("pippx", &pippx, &b_pippx);
   fChain->SetBranchAddress("pippy", &pippy, &b_pippy);
   fChain->SetBranchAddress("pippz", &pippz, &b_pippz);
   fChain->SetBranchAddress("pipp", &pipp, &b_pipp);
   fChain->SetBranchAddress("pipe", &pipe, &b_pipe);
   fChain->SetBranchAddress("pimpx", &pimpx, &b_pimpx);
   fChain->SetBranchAddress("pimpy", &pimpy, &b_pimpy);
   fChain->SetBranchAddress("pimpz", &pimpz, &b_pimpz);
   fChain->SetBranchAddress("pimp", &pimp, &b_pimp);
   fChain->SetBranchAddress("pime", &pime, &b_pime);
   fChain->SetBranchAddress("mpippim", &mpippim, &b_mpippim);
   Notify();
}

Bool_t PipPimNpi::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void PipPimNpi::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t PipPimNpi::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef PipPimNpi_cxx
