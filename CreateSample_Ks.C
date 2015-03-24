
#ifndef Ks0Alg_h
#define Ks0Alg_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

//struct Event;

class Ks0Alg {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Double_t        rec_truth_Mks;
   Double_t        rec_truth_Pks;
   Double_t        rec_truth_Eks;
   Int_t           indexmc;
   Int_t           pdgid[1];   //[indexmc]
   Int_t           motheridx[1];   //[indexmc]
   Int_t           npip;
   Int_t           npim;
   Int_t           nKsCan;
   Double_t        m_chisq0[1];   //[nKsCan]
   Double_t        m_chisqvtxnd[1];   //[nKsCan]
   Double_t        m_decayL[1];   //[nKsCan]
   Double_t        m_decayLerr[1];   //[nKsCan]
   Double_t        m_ctau[1];   //[nKsCan]
   Int_t           Ks_id[1];   //[nKsCan]
   Double_t        Mpippim[1];   //[nKsCan]
   Double_t        Ppippim[1];   //[nKsCan]
   Double_t        Epippim[1];   //[nKsCan]
   Double_t        Thepippim[1];   //[nKsCan]
   Double_t        Phipippim[1];   //[nKsCan]
   Double_t        Ppip[1];   //[nKsCan]
   Double_t        Ppim[1];   //[nKsCan]
   Double_t        Thepip[1];   //[nKsCan]
   Double_t        Thepim[1];   //[nKsCan]
   Double_t        Phipip[1];   //[nKsCan]
   Double_t        Phipim[1];   //[nKsCan]
   Double_t        pippx[1];   //[npip]
   Double_t        pippy[1];   //[npip]
   Double_t        pippz[1];   //[npip]
   Double_t        pipe[1];   //[npip]
   Double_t        pimpx[1];   //[npim]
   Double_t        pimpy[1];   //[npim]
   Double_t        pimpz[1];   //[npim]
   Double_t        pime[1];   //[npim]

   // List of branches
   TBranch        *b_rec_truth_Mks;   //!
   TBranch        *b_rec_truth_Pks;   //!
   TBranch        *b_rec_truth_Eks;   //!
   TBranch        *b_indexmc;   //!
   TBranch        *b_pdgid;   //!
   TBranch        *b_motheridx;   //!
   TBranch        *b_npip;   //!
   TBranch        *b_npim;   //!
   TBranch        *b_nKsCan;   //!
   TBranch        *b_m_chisq0;   //!
   TBranch        *b_m_chisqvtxnd;   //!
   TBranch        *b_m_decayL;   //!
   TBranch        *b_m_decayLerr;   //!
   TBranch        *b_m_ctau;   //!
   TBranch        *b_Ks_id;   //!
   TBranch        *b_Mpippim;   //!
   TBranch        *b_Ppippim;   //!
   TBranch        *b_Epippim;   //!
   TBranch        *b_Thepippim;   //!
   TBranch        *b_Phipippim;   //!
   TBranch        *b_Ppip;   //!
   TBranch        *b_Ppim;   //!
   TBranch        *b_Thepip;   //!
   TBranch        *b_Thepim;   //!
   TBranch        *b_Phipip;   //!
   TBranch        *b_Phipim;   //!
   TBranch        *b_pippx;   //!
   TBranch        *b_pippy;   //!
   TBranch        *b_pippz;   //!
   TBranch        *b_pipe;   //!
   TBranch        *b_pimpx;   //!
   TBranch        *b_pimpy;   //!
   TBranch        *b_pimpz;   //!
   TBranch        *b_pime;   //!

   Ks0Alg(TTree *tree=0);
   virtual ~Ks0Alg();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   
};

#endif





#include "TRandom.h"
#include <iostream>
#include <fstream>
#include "TMath.h"
#include <string>
#include <map>
class Ks0Alg;

int main(int argc,char** argv)
{
	TFile *f = new TFile("Ks_toymc.root","recreate");
	TTree *tree = new TTree("Ks_info","Ks_info");
	Ks0Alg ks;
	ks.Init(tree);


  ofstream sample;
  double args[10];
  sample.open("sample.txt");
  double cmspsig,cmspsig1,cmspsig2,cmspbck,cmspbck1,cmspbck2,msig,mbck;
  double labpsig,labpsig1,labpsig2,labpbck,labpbck1,labpbck2;
  double cmstheta1,cmstheta2,cmsphi1,cmsphi2;
  double labtheta1,labtheta2,labphi1,labphi2;
  double betac,gammac;
  double cmsbeta;
  double pi=TMath::Pi();
  double twopi=2*pi;
  //double mpeak=1.019455;
  //double mpeak=3.0969;
  //double sigma=0.0025;
  //double mD0 = 1.86484;
  //double mk=0.493677;
  double mKs = 0.497614;
  double mpi = 0.13957;
  double &sigNo=args[1]=100000;
  double &backNo=args[2]=0;
  double &mpeak=args[3]=mKs;
  double &sigma=args[4]=0.005;
  double &cmsp=args[5]=1.1;
  double &runmode=args[6]=0;
  double &factor=args[7]=1;

  // deal with input arguments
  std::map<std::string,int> ArgMap;
  ArgMap.insert(std::make_pair("-signal",1));
  ArgMap.insert(std::make_pair("-s",1));
  ArgMap.insert(std::make_pair("-background",2));
  ArgMap.insert(std::make_pair("-b",2));
  ArgMap.insert(std::make_pair("-peak",3));
  ArgMap.insert(std::make_pair("-m",3));
  ArgMap.insert(std::make_pair("-sigma",4));
  ArgMap.insert(std::make_pair("-w",4));
  ArgMap.insert(std::make_pair("-cmsp",5));
  ArgMap.insert(std::make_pair("-pc",5));
  ArgMap.insert(std::make_pair("-runmode",6));
  ArgMap.insert(std::make_pair("-k",7));
  ArgMap.insert(std::make_pair("-f",7));
  ArgMap.insert(std::make_pair("-factor",7));
  // the following map just for showing information
  std::map<std::string,int> InfoMap;
  InfoMap["SignalNumber    "] = 1;
  InfoMap["BackgroundNumber"] = 2;
  InfoMap["InvariantMass   "] = 3;
  InfoMap["Sigma           "] = 4;
  InfoMap["CMSMomentum     "] = 5;
  InfoMap["RunMode         "] = 6;
  InfoMap["Factor          "] = 7;
  //InfoMap["SignalNumber"]=1;

  std::cout<<"argc is "<<argc<<"\n";
  for (int i=1;i<argc;i++){
    std::cout<<"arg "<<i<<" is "<<argv[i]<<"\n";
    std::string tmpstr=argv[i];
    if(tmpstr.find("-")==0 && argc>i){
      std::cout<<"setting arg: "<<&tmpstr[1]<<"\n";
      args[ArgMap[tmpstr]] = atof(argv[i+1]);
      i++;
    }
  }
  for (std::map<std::string,int>::iterator iter=InfoMap.begin();iter!=InfoMap.end();iter++){
    std::cout<<iter->first<<" : "<<args[iter->second]<<"\n";
  } 

  //double mk=0.000511;
  double c = 299792458;
  gRandom->SetSeed(time(0));
  //signal part,
  if (runmode>0) { 
    for(int i=0;i <sigNo;i++){
      msig=gRandom->Gaus(mpeak,sigma)/factor;
      sample<<"\n"<<msig;
    }
  }
  else {
    //std::cout<<sigNo<<std::endl;
    for (int i=0;i<sigNo;i++){
      double alpha, beta;
	  double theta;
      //labpsig1 = gRandom->Uniform(0.2, 2.0);
	  //labpsig2 = gRandom->Uniform(0.2, 2.0);
      labpsig1 = gRandom->Gaus(1.1, 0.20);
	  labpsig2 = gRandom->Gaus(1.1, 0.20);
	  double costheta = gRandom->Uniform(-1,1);
      theta = acos(costheta);
      //double cosalpha = gRandom->Uniform(0,1);
	  //alpha = acos(cosalpha);
	  // ~~~~~~~~~
	  // calculate angle: alpha and beta

	//double cosAB = (2*sqrt(labpsig1*labpsig1+mpi*mpi)*sqrt(labpsig2*labpsig2+mpi*mpi)+2*mpi*mpi-mpeak*mpeak)/(2*labpsig1*labpsig2);
	//if ( fabs(cosAB)>1 ){
	//  i--;
	//  continue;
	//}
	//double consA = labpsig2*labpsig2*(1-cosAB*cosAB)/(labpsig1*labpsig1+labpsig2*labpsig2+2*labpsig1*labpsig2*cosAB);
	//alpha = asin(sqrt(consA));
	//beta = asin(labpsig1/labpsig2*sqrt(consA));
	  
	  double totE,E1,E2;
	  E1 = sqrt(labpsig1*labpsig1+mpi*mpi);
	  E2 = sqrt(labpsig2*labpsig2+mpi*mpi);
	  totE = E1+E2;
	  double totPsq;
	  totPsq = totE*totE - mpeak*mpeak;
	  double cosal = (totPsq + labpsig1*labpsig1 - labpsig2*labpsig2)/(2*sqrt(totPsq)*labpsig1);
      if (fabs(cosal)<1) alpha = acos(cosal);
	  else {
	    i--;
	    continue;
	  }
	  beta = asin(labpsig1*sin(alpha)/labpsig2);

	  //beta = acos(cosAB) - alpha;
	  labtheta1 = theta + alpha;
	  labtheta2 = theta - beta;
	  labphi1 = gRandom->Uniform(0,pi);
	  labphi2 = labphi1;
	  if (theta + alpha > pi) labphi1 = labphi1+pi;
	  if (theta < beta ) labphi2 = labphi2 + pi;
	  //labphi2 = pi + labphi1;
	  //if (labphi2 > twopi) labphi2 -= twopi;

	  // smear momentum in lab system
	  labpsig1 = gRandom->Gaus(labpsig1, 0.01*labpsig1);
	  labpsig2 = gRandom->Gaus(labpsig2, 0.01*labpsig2);
	  //labpsig1 = gRandom->Uniform(labpsig1-0.01*labpsig1, labpsig1+0.01*labpsig1);
	  //labpsig2 = gRandom->Uniform(labpsig2-0.01*labpsig2, labpsig2+0.01*labpsig2);


      double p1[3],p2[3];
      p1[0] = labpsig1*sin(labtheta1)*cos(labphi1)/factor;
      p1[1] = labpsig1*sin(labtheta1)*sin(labphi1)/factor;
      p1[2] = labpsig1*cos(labtheta1)/factor;
      p2[0] = labpsig2*sin(labtheta2)*cos(labphi2)/factor;
      p2[1] = labpsig2*sin(labtheta2)*sin(labphi2)/factor;
      p2[2] = labpsig2*cos(labtheta2)/factor;
    //p1[0] = labpsig1*sin(labtheta1)*cos(0)/factor;
    //p1[1] = labpsig1*sin(labtheta1)*sin(0)/factor;
    //p1[2] = labpsig1*cos(labtheta1)/factor;
    //p2[0] = labpsig2*sin(labtheta2)*cos(pi)/factor;
    //p2[1] = labpsig2*sin(labtheta2)*sin(pi)/factor;
    //p2[2] = labpsig2*cos(labtheta2)/factor;
      //psig1 = gRandom->Gaus(1.,0.5);
      //psig2 = TMath::Sqrt(mpeak*mpeak - msig*msig)-psig1;
      sample<<"\n"<<p1[0]<<"\t"<<p1[1]<<"\t"<<p1[2]<<"\t";
      sample<<"\t"<<p2[0]<<"\t"<<p2[1]<<"\t"<<p2[2]<<"\t";
 
	  ks.pippx[0] = p1[0];
	  ks.pippy[0] = p1[1];
	  ks.pippz[0] = p1[2];
	  ks.pimpx[0] = p2[0];
	  ks.pimpy[0] = p2[1];
	  ks.pimpz[0] = p2[2];
	  ks.Ks_id[0] = 1;
	  ks.npip     = 1;
	  ks.npim     = 1;
	  ks.nKsCan   = 1;
	  ks.Mpippim[0]=mpeak;
	  ks.m_decayL[0]=3.0;
	  ks.m_decayLerr[0] = 1.0;

      //ks.pim4 = msig;
	  tree->Fill();
    }
	tree->Write();
    //background part,
    sample<<"\n";
    
    for(int i=0;i<backNo;i++){
      //initial
      msig    = gRandom->Uniform(mpeak-10*sigma,mpeak+30*sigma);
      //cmspsig = gRandom->Exp(0.5);
      cmspsig = gRandom->Gaus(0,3);
      betac   = cmspsig/TMath::Sqrt(msig*msig+cmspsig*cmspsig);
      gammac  = TMath::Sqrt(msig*msig+cmspsig*cmspsig)/msig;
      //split
      //cms
      cmspsig1 = TMath::Sqrt(msig*msig/4.0 - mpi*mpi);
      cmspsig2 = cmspsig1;
      cmstheta1 = gRandom->Uniform(0,pi);
      cmsphi1   = gRandom->Uniform(0,twopi);
      cmstheta2 = pi - cmstheta1;
      cmsphi2   = cmsphi1 + pi;
      // change to lab
      cmsbeta   = cmspsig1 / (msig/2.);
      if(cos(cmstheta1)+betac/cmsbeta==0){
        labtheta1 = pi/2;
      }
      else
        labtheta1 = atan(sin(cmstheta1)/(gammac*(cos(cmstheta1)+betac/cmsbeta)));
      if(labtheta1<0) labtheta1 += pi;
      labpsig1 = cmspsig1*sin(cmstheta1)/sin(labtheta1);
      labphi1  = cmsphi1;
      if(cos(cmstheta2)+betac/cmsbeta==0){
        labtheta2 = pi/2;
      }else
        labtheta2 = atan(sin(cmstheta2)/(gammac*(cos(cmstheta2)+betac/cmsbeta)));
      if(labtheta2<0) labtheta2 += pi;
      labpsig2 = cmspsig1*sin(cmstheta2)/sin(labtheta2);
      labphi2  = cmsphi2;
      double p1[3],p2[3];
      p1[0] = labpsig1*sin(labtheta1)*cos(labphi1)/factor;
      p1[1] = labpsig1*sin(labtheta1)*sin(labphi1)/factor;
      p1[2] = labpsig1*cos(labtheta1)/factor;
      p2[0] = labpsig2*sin(labtheta2)*cos(labphi2)/factor;
      p2[1] = labpsig2*sin(labtheta2)*sin(labphi2)/factor;
      p2[2] = labpsig2*cos(labtheta2)/factor;
      //psig1 = gRandom->Gaus(1.,0.5);
      //psig2 = TMath::Sqrt(mpeak*mpeak - msig*msig)-psig1;
      sample<<"\n"<<p1[0]<<"\t"<<p1[1]<<"\t"<<p1[2]<<"\t";
      sample<<"\t"<<p2[0]<<"\t"<<p2[1]<<"\t"<<p2[2]<<"\t";
    //  mbck=gRandom->Uniform(mpeak-30*sigma,mpeak-30*sigma);
    //  pbck1 = gRandom->Uniform(0,2.0);
    //  pbck2 = TMath::Sqrt(mbck);
    //  sample<<"\n"<<pbck;
    }
  } 
  sample.close();
  return 0;
}


//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Feb  6 09:56:41 2015 by ROOT version 5.24/00b
// from TTree Ks_info/track N-Tuple example
// found on file: Ksto2pi_583_140226_0035862.root
//////////////////////////////////////////////////////////
#define Ks0Alg_cxx
//#include "Ks0Alg.h"
#include <TH2.h>

void Ks0Alg::Loop(){}

#ifdef Ks0Alg_cxx
Ks0Alg::Ks0Alg(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   
   Init(tree);
}

Ks0Alg::~Ks0Alg()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Ks0Alg::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Ks0Alg::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void Ks0Alg::Init(TTree *tree)
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

   fChain->Branch("rec_truth_Mks"	, &rec_truth_Mks, "rec_truth_Mks/D"	);
   fChain->Branch("rec_truth_Pks"	, &rec_truth_Pks, "rec_truth_Pks/D"	);
   fChain->Branch("rec_truth_Eks"	, &rec_truth_Eks, "rec_truth_Eks/D"	);
   fChain->Branch("indexmc"			, &indexmc		, "indexmc/I"		);
   fChain->Branch("pdgid"			, pdgid			, "pdgid[1]/I"		);
   fChain->Branch("motheridx"		, motheridx		, "motheridx[1]/I"	);
   fChain->Branch("npip"			, &npip			, "npip/I"			);
   fChain->Branch("npim"			, &npim			, "npim/I"			);
   fChain->Branch("nKsCan"			, &nKsCan		, "nKsCan/I"		);
   fChain->Branch("m_chisq0"		, m_chisq0		, "m_chisq0[1]/D"		);
   fChain->Branch("m_chisqvtxnd"	, m_chisqvtxnd	, "m_chisqvtxnd[1]/D"	);
   fChain->Branch("m_decayL"		, m_decayL		, "m_decayL[1]/D"		);
   fChain->Branch("m_decayLerr"		, m_decayLerr	, "m_decayLerr[1]/D"	);
   fChain->Branch("m_ctau"			, m_ctau		, "m_ctau[1]/D"		);
   fChain->Branch("Ks_id"			, Ks_id			, "Ks_id[1]/I"		);
   fChain->Branch("Mpippim"			, Mpippim		, "Mpippim[1]/D"	);
   fChain->Branch("Ppippim"			, Ppippim		, "Ppippim[1]/D"	);
   fChain->Branch("Epippim"			, Epippim		, "Epippim[1]/D"	);
   fChain->Branch("Thepippim"		, Thepippim		, "Thepippim[1]/D"	);
   fChain->Branch("Phipippim"		, Phipippim		, "Phipippim[1]/D"	);
   fChain->Branch("Ppip"			, Ppip			, "Ppip[1]/D"		);
   fChain->Branch("Ppim"			, Ppim			, "Ppim[1]/D"		);
   fChain->Branch("Thepip"			, Thepip		, "Thepip[1]/D"		);
   fChain->Branch("Thepim"			, Thepim		, "Thepim[1]/D"		);
   fChain->Branch("Phipip"			, Phipip		, "Phipip[1]/D"		);
   fChain->Branch("Phipim"			, Phipim		, "Phipim[1]/D"		);
   fChain->Branch("pippx"			, pippx			, "pippx[1]/D"		);
   fChain->Branch("pippy"			, pippy			, "pippy[1]/D"		);
   fChain->Branch("pippz"			, pippz			, "pippz[1]/D"		);
   fChain->Branch("pipe"			, pipe			, "pipe[1]/D"		);
   fChain->Branch("pimpx"			, pimpx			, "pimpx[1]/D"		);
   fChain->Branch("pimpy"			, pimpy			, "pimpy[1]/D"		);
   fChain->Branch("pimpz"			, pimpz			, "pimpz[1]/D"		);
   fChain->Branch("pime"			, pime			, "pime[1]/D"		);
   Notify();
}

Bool_t Ks0Alg::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Ks0Alg::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Ks0Alg::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Ks0Alg_cxx


