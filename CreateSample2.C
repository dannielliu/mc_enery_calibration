
#ifndef gepep_kk_h
#define gepep_kk_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class gepep_kk {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           run;
   Int_t           rec;
   Int_t           evttag;
   Int_t           indexmc;
   Int_t           pdgid[100];   //[indexmc]
   Int_t           motheridx[100];   //[indexmc]
   Int_t           ngch;
   Int_t           ncharg;
   Int_t           nneu;
   Double_t        kappx;
   Double_t        kappy;
   Double_t        kappz;
   Double_t        kape;
   Double_t        kampx;
   Double_t        kampy;
   Double_t        kampz;
   Double_t        kame;
   Double_t        mphi;
   Double_t        kkm4;
   
   
   gepep_kk(TTree *tree=0);
   virtual ~gepep_kk();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   //virtual void     Loop();
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

int main(int argc,char** argv)
{
	TFile *f = new TFile("kk_toymc.root","recreate");
	TTree *tree = new TTree("gepep_kk","gepep_kk");
	gepep_kk kk;
	kk.Init(tree);


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
  double mD0 = 1.86484;
  double &sigNo=args[1]=100000;
  double &backNo=args[2]=0;
  double &mpeak=args[3]=mD0;
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

  double mk=0.493677;
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
      labpsig1 = gRandom->Gaus(1.5, 0.15);
	  labpsig2 = gRandom->Gaus(1.5, 0.15);
	  double costheta = gRandom->Uniform(-1,1);
      theta = acos(costheta);
      //double cosalpha = gRandom->Uniform(0,1);
	  //alpha = acos(cosalpha);
	  // ~~~~~~~~~
	  // calculate angle: alpha and beta

	//double cosAB = (2*sqrt(labpsig1*labpsig1+mk*mk)*sqrt(labpsig2*labpsig2+mk*mk)+2*mk*mk-mpeak*mpeak)/(2*labpsig1*labpsig2);
	//if ( fabs(cosAB)>1 ){
	//  i--;
	//  continue;
	//}
	//double consA = labpsig2*labpsig2*(1-cosAB*cosAB)/(labpsig1*labpsig1+labpsig2*labpsig2+2*labpsig1*labpsig2*cosAB);
	//alpha = asin(sqrt(consA));
	//beta = asin(labpsig1/labpsig2*sqrt(consA));
	  
	  double totE,E1,E2;
	  E1 = sqrt(labpsig1*labpsig1+mk*mk);
	  E2 = sqrt(labpsig2*labpsig2+mk*mk);
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
	  if (theta + alpha > pi) labphi1 += pi;
	  if (theta < beta) labphi2 += pi;

	  // smear momentum in lab system
	  labpsig1 = gRandom->Gaus(labpsig1, 0.01*labpsig1);
	  labpsig2 = gRandom->Gaus(labpsig2, 0.01*labpsig2);
	  //labpsig1 = gRandom->Uniform(labpsig1-0.01*labpsig1, labpsig1+0.01*labpsig1);
	  //labpsig2 = gRandom->Uniform(labpsig2-0.01*labpsig2, labpsig2+0.01*labpsig2);


      double p1[3],p2[3];
    //p1[0] = labpsig1*sin(labtheta1)*cos(labphi1)/factor;
    //p1[1] = labpsig1*sin(labtheta1)*sin(labphi1)/factor;
    //p1[2] = labpsig1*cos(labtheta1)/factor;
    //p2[0] = labpsig2*sin(labtheta2)*cos(labphi2)/factor;
    //p2[1] = labpsig2*sin(labtheta2)*sin(labphi2)/factor;
    //p2[2] = labpsig2*cos(labtheta2)/factor;
      p1[0] = labpsig1*sin(labtheta1)*cos(0)/factor;
      p1[1] = labpsig1*sin(labtheta1)*sin(0)/factor;
      p1[2] = labpsig1*cos(labtheta1)/factor;
      p2[0] = labpsig2*sin(labtheta2)*cos(0)/factor;
      p2[1] = labpsig2*sin(labtheta2)*sin(0)/factor;
      p2[2] = labpsig2*cos(labtheta2)/factor;
      //psig1 = gRandom->Gaus(1.,0.5);
      //psig2 = TMath::Sqrt(mpeak*mpeak - msig*msig)-psig1;
      sample<<"\n"<<p1[0]<<"\t"<<p1[1]<<"\t"<<p1[2]<<"\t";
      sample<<"\t"<<p2[0]<<"\t"<<p2[1]<<"\t"<<p2[2]<<"\t";
 
	  kk.kappx = p1[0];
	  kk.kappy = p1[1];
	  kk.kappz = p1[2];
	  kk.kampx = p2[0];
	  kk.kampy = p2[1];
	  kk.kampz = p2[2];
      kk.kkm4 = msig;
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
      cmspsig1 = TMath::Sqrt(msig*msig/4.0 - mk*mk);
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
      if (labtheta2<0) labtheta2 += pi;
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



#define gepep_kk_cxx
//#include <TH2.h>
//#include <TH3.h>
//#include <TStyle.h>
//#include <TCanvas.h>
//#include "function.h"
//#include "TF1.h"
//#include "TPaveText.h"
//#include "TGraphErrors.h"
//#include <fstream>
//#include "RooFit.h"
//#include "RooFitResult.h"
//#include "RooRealVar.h"
//#include "RooGaussian.h"
//#include "RooChebychev.h"
//#include "RooCBShape.h"
//#include "RooExponential.h"
//#include "RooPolynomial.h"
//#include "RooBreitWigner.h"
////#include "RooLandau.h"
//#include "RooDataHist.h"
//#include "RooDataSet.h"
//#include "RooAddPdf.h"
//#include "RooArgList.h"
//#include "RooPlot.h"
//#include <iostream>
//extern std::string outputdir;
//using RooFit::Title;
//using RooFit::Components;
//using RooFit::LineStyle;
//using RooFit::LineColor;
//using RooFit::Range;
//using namespace RooFit;
/*
void gepep_kk::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L gepep_kk.C
//      Root > gepep_kk t
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

   Long64_t nentries = fChain->GetEntriesFast();

   std::cout<<"Toral entry is "<<nentries<<std::endl;
   int nBins=200;
   double factorstart=0.99;
   double mk=0.493677;
   // D0 -> K K
   double philow=1.82;
   double phiup=1.90;
   double peakvalue=1.86484;// mphi
   // phi -> K K
   //double philow=0.995;
   //double phiup=1.045;
   //double peakvalue=1.019455;// mphi
   int pointNo=10;
   double factors[pointNo];
   double factorserr[pointNo];
   double deltapeaks[pointNo];
   double deltapeakserr[pointNo];
   double factor=factorstart;
   double factorstep=(1.-factor)*2/pointNo;
   
   // try to use roofit
   RooRealVar x("x","energy",peakvalue,philow,phiup,"GeV");
   //RooRealVar x("x","energy",peakvalue,1.015,1.025,"GeV");
   RooRealVar mean("mean","mean of gaussian",peakvalue,philow,phiup);
   RooRealVar mean2("mean2","mean of gaussian2",peakvalue,philow,phiup);
   //RooRealVar sigma("sigma","width of gaussian",0.0023,0.0015,0.0040);//phi version
   //RooRealVar sigma2("sigma2","width of gaussian",0.005,0.003,0.008);
   RooRealVar sigma("sigma","width of gaussian",0.007,0.003,0.0075);//D0 version
   RooRealVar sigma2("sigma2","width of gaussian",0.02,0.008,0.05);
   RooGaussian gaus("gaus","gauss(x,m,s)",x,mean,sigma);
   RooGaussian gaus2("gaus2","gauss(x,m,s)",x,mean,sigma2);
   
   RooRealVar co1("co1","coefficient #1",0,-0.5,0.5);
   RooRealVar co2("co2","coefficient #2",0,-0.01,0.5);
   RooChebychev bkg("bkg","background",x,RooArgList(co1,co2));
   
   RooRealVar alpha("alpha","#alpha",1.16,-5,5);
   RooRealVar nnn("nnn","n",50,1,200);
   RooCBShape cbshape("cbshape","crystal ball",x,mean,sigma,alpha,nnn);

   RooBreitWigner brewig("brewig","brewig",x,mean,sigma);
   //RooLandau landau("landau","landau",x,mean,sigma);

   RooRealVar tao("tao","#tao",-200.,-1000.,-100.);
   RooExponential expo("expo","expo",x,tao);

   RooRealVar a0("a0","coefficient #0",100,100,100000);
   RooRealVar a1("a1","coefficient #1",-50,-100000,100000);
   RooPolynomial ground("ground","ground",x,RooArgList(a0,a1));
   
   RooRealVar signal("signal"," ",1200,10,1000000);//event number
   RooRealVar signal2("signal2"," ",1200,0,1000000);//event number
   RooRealVar background("background"," ",200,0,1000000);
   RooRealVar background2("background2"," ",200,0,1000000);
   
   RooPlot *xframe;
   //RooDataHist *data_k;
   RooDataSet *dataset;
   RooAddPdf *sum;
   
   //RooDataSet *dataset = new RooDataSet("dataset","data",dataraw,x);
   RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR); // set out put message level of roofit

   int Npar;
   char fname[1000];
   ofstream ofpar;
   sprintf(fname,"%s/pars.txt",outputdir.c_str());
   ofpar.open(fname,std::ios::app);
   ofstream purepar;
   sprintf(fname,"%s/parspure.txt",outputdir.c_str());
   purepar.open(fname,std::ios::app);
   sprintf(fname,"%s/plot_kk.root",outputdir.c_str());
   TFile *f = new TFile(fname,"RECREATE");
   TTree *dataraw=new TTree("dataraw","dataraw");
   double mass;
   dataraw->Branch("x",&mass,"x/D");
   TTree *vars = new TTree("vars","vars");
   double phi1,phi2;
   double costheta1,costheta2;
   double p1,p2;
   vars->Branch("phi1",&phi1,"phi1/D");
   vars->Branch("phi2",&phi2,"phi2/D");
   vars->Branch("costheta1",&costheta1,"costheta1/D");
   vars->Branch("costheta2",&costheta2,"costheta2/D");
   vars->Branch("p1",&p1,"p1/D");
   vars->Branch("p2",&p2,"p2/D");
   vars->Branch("mass",&mass,"mass/D");


   TF1 *facfit = new TF1("facfit",line2,0.9,1.1,2);
   TH1D *h1    = new TH1D("h1","2 kaon invariant mass",nBins,philow,phiup);
   TCanvas *c1 = new TCanvas("","",800,600);

   const int Npart=20;
   double m0=peakvalue;
   double sigma_m=0.0024;//0.0024 for phi,
   double width = 10.*sigma_m;
   double mparticle=0.493677;
   std::vector<int> partmap;
   std::vector<std::pair<int,double> > facmap;

   int realsize=0;
   double partid[Npart];
   double parter[Npart];
   double corfac[Npart];
   double corerr[Npart];
   
   double pcut[Npart+1];
   double facv[Npart];
   double facev[Npart];


   //pcut[0]=0;
   //pcut[1] = 2.0;
   pcut[0] =0.0 ;    facv[0] =1.0;       facev[0] =1.0;  
   pcut[1] =0.10;    facv[1] =1.0     ;  facev[1] =1.0        ;  
   pcut[2] =0.20;    facv[2] =1.0     ;  facev[2] =1.0        ;
   pcut[3] =0.30;    facv[3] =1.0     ;  facev[3] =1.0        ;
   pcut[4] =0.40;    facv[4] =1.0     ;  facev[4] =1.0        ;
   pcut[5] =0.50;    facv[5] =1.00338 ;  facev[5] =1.0        ;
   pcut[6] =0.60;    facv[6] =1.00484 ;  facev[6] =1.0        ;
   pcut[7] =0.70;    facv[7] =1.0024  ;  facev[7] =1.0        ;
   pcut[8] =0.80;    facv[8] =0.999958;  facev[8] =1.0        ;
   pcut[9] =0.90;    facv[9] =0.997556;  facev[9] =1.0        ;
   pcut[10]=1.00;    facv[10]=0.997925;  facev[10]=1.0        ;
   pcut[11]=1.10;    facv[11]=0.999218;  facev[11]=1.0        ;
   pcut[12]=1.20;    facv[12]=1.00027 ;  facev[12]=1.0        ;
   pcut[13]=1.30;    facv[13]=1.0     ;  facev[13]=1.0        ;
   pcut[14]=1.40;    facv[14]=1.0     ;  facev[14]=1.0        ;
   pcut[15]=1.50;    facv[15]=1.0     ;  facev[15]=1.0        ;
   pcut[16]=1.60;    facv[16]=1.0     ;  facev[16]=1.0        ;
   pcut[17]=1.70;    facv[17]=1.0     ;  facev[17]=1.0        ;
   pcut[18]=1.80;    facv[18]=1.0;       facev[18]=1.0;  
   pcut[19]=1.90;    facv[19]=1.0;       facev[19]=1.0;  
   pcut[20]=2.00;           

   char name[100];
   TH1D *hmD0[Npart];
   for (int partj=0;partj<Npart;partj++){
     sprintf(name,"mass_part%0d",partj);
     hmD0[partj] = new TH1D(name,name,100,philow,phiup);
   }
 
   // loop data
   Event evt(mk);
   std::vector<Event> evts[Npart];
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
     
	    evt.SetVal(kappx,kappy,kappz,kampx,kampy,kampz);
        mass = evt.InvMass();
        p1 = evt.GetP1();
        p2 = evt.GetP2();
        costheta1 = evt.GetCostheta1();
		costheta2 = evt.GetCostheta2();
        phi1 = evt.GetPhi1();
        phi2 = evt.GetPhi2();
        if (mass>philow && mass<phiup){
          vars->Fill();
        }
        if (costheta1 > 0.8 && costheta2 <-0.8) continue;
		if (p1<0.5 || p1>1.3) continue;

        //if ( partj>=Npart || partj<0 ) continue;
        for (int partj=0;partj<Npart;partj++){
          if (p2<pcut[partj] || p2>pcut[partj+1]) continue;
          //if (p2<pcut[partj] || p2>pcut[partj+1]) continue;
          if (mass>philow-0.02 && mass<phiup+0.02){
            hmD0[partj]->Fill(mass);
			evts[partj].push_back(evt);
	      }
          break;
		}
   
   }
   std::cout<<"fffffffffffffffff"<<std::endl;
   vars->Write();
   for (int partj=0;partj<Npart;partj++){
     hmD0[partj]->Write();
     if (hmD0[partj]->GetEntries() > 100
       &&hmD0[partj]->GetMaximumBin()> 30 //check the the peak position,
       &&hmD0[partj]->GetMaximumBin()< 70 
       ){
       partmap.push_back(partj);
     }
   }

   // ~~~~~~~~ draw end

   for (int loopj=0;loopj<partmap.size();loopj++){
     int partj=partmap.at(loopj);

   // for saving the fit result
   double factori = 1.0;
   int fittimes = 0;
   factor = factorstart;
   for (int i=0;i<pointNo;i++){
      //xframe->Clear();
      xframe = x.frame(Title("fit kaon"));

      //h1->Reset();
      dataraw->Reset();
      std::cout<<"factor is "<<factor<<std::endl;
      for (Long64_t jentry=0; jentry<evts[partj].size();jentry++) {
		  p1 = evts[partj].at(jentry).GetP1();
		  p2 = evts[partj].at(jentry).GetP2();
          for (int i=0;i<Npart;i++){
            if (p1>=pcut[i]&&p1<pcut[i+1]){
              factori = facv[i];
              break;
            }
          }
          mass = evts[partj].at(jentry).InvMass(factori,factor);
          if (mass>philow && mass<phiup){
            dataraw->Fill();
          }
         // if (Cut(ientry) < 0) continue;
      }
      //dataraw->Write();

      dataset = new RooDataSet("dataset","data",dataraw,x);
      char tmpchr[100];
      sprintf(tmpchr,"data_k_%02d",fittimes);
      //data_k = new RooDataHist(tmpchr,"data_k",x,h1);
      sum = new RooAddPdf("sum","sum",RooArgList(gaus,gaus2,bkg),RooArgList(signal,signal2,background));
      Npar=8;
      //sum = new RooAddPdf("sum","sum",RooArgList(gaus,ground),RooArgList(signal,background2));
      mean.setVal(peakvalue+0.06*(factor-1.0));
      sum->fitTo(*dataset,Range(philow,phiup));
      dataset->plotOn(xframe);
      sum->plotOn(xframe,Components(gaus),LineStyle(2),LineColor(2));
      sum->plotOn(xframe,Components(gaus2),LineStyle(2),LineColor(2));
      sum->plotOn(xframe,Components(ground),LineStyle(3),LineColor(3));
      sum->plotOn(xframe);
      xframe->Draw();
      TPaveText *pt = new TPaveText(0.12,0.50,0.5,0.90,"BRNDC");
      pt->SetBorderSize(0);
      pt->SetFillStyle(4000);
      pt->SetTextAlign(12);
      pt->SetTextFont(42);
      pt->SetTextSize(0.035);
      sprintf(tmpchr,"#mu_{1} = %1.6f #pm %1.6f",mean.getVal(),mean.getError());
      pt->AddText(tmpchr);
      sprintf(tmpchr,"#sigma_{1} = %1.6f #pm %1.6f",sigma.getVal(),sigma.getError());
      pt->AddText(tmpchr);
      sprintf(tmpchr,"#sigma_{2} = %1.6f #pm %1.6f",sigma2.getVal(),sigma2.getError());
      pt->AddText(tmpchr);
      sprintf(tmpchr,"signal1 = %.2f #pm %.2f",signal.getVal(),signal.getError());
      pt->AddText(tmpchr);
      sprintf(tmpchr,"signal2 = %.2f #pm %.2f",signal2.getVal(),signal2.getError());
      pt->AddText(tmpchr);
      sprintf(tmpchr,"backNo = %.2f #pm %.2f",background.getVal(),background.getError());
      pt->AddText(tmpchr);
      sprintf(tmpchr,"#chi^{2}/(100-%d) = %5.6f",Npar,xframe->chiSquare(Npar));
      pt->AddText(tmpchr);
      pt->Draw();
      sprintf(name,"part%d_fitFor_%dth_time",partj,fittimes);
      c1->SetName(name);
      c1->Write();
      //delete data_k;
      delete dataset;
      delete xframe;
      delete sum;

      //sprintf(tmpchr,"data_k_%d.eps",fittimes);
      //h1->Draw();
      //c1->Print(tmpchr);
      
      // save pars
      factors[i]=factor;
      factorserr[i]=0;
      deltapeaks[i] = mean.getVal() - peakvalue;
      deltapeakserr[i] = mean.getError();

      fittimes++;
      factor += factorstep;
   }
   
   TGraphErrors *graph1 = new TGraphErrors(pointNo,factors,deltapeaks,factorserr,deltapeakserr);
   graph1->SetTitle("delta peak");
   graph1->SetMarkerStyle(5);
   graph1->Draw("AP");
   gStyle->SetOptFit(1111);
   facfit->SetParameters(1,0.3);
   facfit->SetParNames("factor","slope");
   graph1->Fit(facfit,"","",factors[0],factors[pointNo-1]);
   //factor1=facfit->GetParameter(0);
   //factor1err=facfit->GetParError(0);
   factor = facfit->GetParameter(0);
   sprintf(name,"factors_kk_part%d",partj);
   graph1->SetName(name);
   graph1->Write();

   // draw the best fitting
      xframe = x.frame(Title("fit kaon"));
      dataraw->Reset();
      std::cout<<"factor is "<<factor<<std::endl;
      for (Long64_t jentry=0; jentry<evts[partj].size();jentry++) {
		 p1 = evts[partj].at(jentry).GetP1();
		 p2 = evts[partj].at(jentry).GetP2();
         for (int i=0;i<Npart;i++){
            if (p1>=pcut[i]&&p1<pcut[i+1]){
              factori = facv[i];
              break;
            }
         }
         mass = evts[partj].at(jentry).InvMass(factori,factor);
         if (mass>philow && mass<phiup){
           dataraw->Fill();
         }
      }
      //dataraw->Write();

      char tmpchr[100];
      sprintf(tmpchr,"data_k");
      dataset = new RooDataSet("dataset","data",dataraw,x);
      sum = new RooAddPdf("sum","sum",RooArgList(gaus,gaus2,bkg),RooArgList(signal,signal2,background));
      Npar=8;
      //sum = new RooAddPdf("sum","sum",RooArgList(gaus,ground),RooArgList(signal,background2));
      mean.setVal(peakvalue+0.06*(factor-1.0));
      sum->fitTo(*dataset,Range(philow,phiup));
      dataset->plotOn(xframe);
      sum->plotOn(xframe,Components(gaus),LineStyle(2),LineColor(2));
      sum->plotOn(xframe,Components(gaus2),LineStyle(2),LineColor(2));
      //sum->plotOn(xframe,Components(expo),LineStyle(3),LineColor(3));
      sum->plotOn(xframe,Components(ground),LineStyle(3),LineColor(3));
      sum->plotOn(xframe);
      xframe->Draw();
      TPaveText *pt = new TPaveText(0.12,0.50,0.5,0.90,"BRNDC");
      pt->SetBorderSize(0);
      pt->SetFillStyle(4000);
      pt->SetTextAlign(12);
      pt->SetTextFont(42);
      pt->SetTextSize(0.035);
      sprintf(tmpchr,"#mu_{1} = %1.6f #pm %1.6f",mean.getVal(),mean.getError());
      pt->AddText(tmpchr);
      sprintf(tmpchr,"#sigma_{1} = %1.6f #pm %1.6f",sigma.getVal(),sigma.getError());
      pt->AddText(tmpchr);
      sprintf(tmpchr,"#sigma_{2} = %1.6f #pm %1.6f",sigma2.getVal(),sigma2.getError());
      pt->AddText(tmpchr);
      sprintf(tmpchr,"signal1 = %.2f #pm %.2f",signal.getVal(),signal.getError());
      pt->AddText(tmpchr);
      sprintf(tmpchr,"signal2 = %.2f #pm %.2f",signal2.getVal(),signal2.getError());
      pt->AddText(tmpchr);
      sprintf(tmpchr,"backNo = %.2f #pm %.2f",background.getVal(),background.getError());
      pt->AddText(tmpchr);
      sprintf(tmpchr,"#chi^{2}/(100-%d) = %5.6f",Npar,xframe->chiSquare(Npar));
      pt->AddText(tmpchr);
      double factor4err=TMath::Sqrt(TMath::Power(mean.getError()/facfit->GetParameter(1),2) + TMath::Power(facfit->GetParError(0),2));
      sprintf(name,"factor = %.6f #pm %.6f",factor,factor4err);
      pt->AddText(name);

      pt->Draw();
      sprintf(name,"part%d_final_fit",partj);
      c1->SetName(name);
      c1->Write();      
      
      partid[loopj] = pcut[partj]+(pcut[partj+1]-pcut[partj])/2;
      parter[loopj] = 0;
      corfac[loopj] = factor;
      corerr[loopj] = factor4err;
   

	  //delete data_k;
      delete dataset;
      delete xframe;
      delete sum;
   }

  realsize = partmap.size();
   
  for (int i=0;i<realsize;i++){
    ofpar<<"p="<<partid[i]<<"\tfactor: "<<corfac[i]<<"\t +/- \t"<< corerr[i]<<std::endl;
    purepar<<partid[i]<<"\t"<<corfac[i]<<"\t"<< corerr[i]<<std::endl;
  }
   
   
   f->Close();
   return;
}
*/
#ifdef gepep_kk_cxx
gepep_kk::gepep_kk(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
 //if (tree == 0) {
 //   TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("data/RValue_kk_3850.root");
 //   if (!f || !f->IsOpen()) {
 //      f = new TFile("data/RValue_kk_3850.root");
 //   }
 //   f->GetObject("gepep_kk",tree);

 //}
   Init(tree);
}

gepep_kk::~gepep_kk()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t gepep_kk::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t gepep_kk::LoadTree(Long64_t entry)
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

void gepep_kk::Init(TTree *tree)
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

   fChain->Branch("run", &run,"run/I");
   fChain->Branch("rec", &rec,"rec/I");
   fChain->Branch("evttag", &evttag,   "evttag/I");
   fChain->Branch("indexmc", &indexmc, "indexmc/I");
   fChain->Branch("pdgid"   , pdgid,   "pdgid[100]/I");
   fChain->Branch("motheridx", motheridx, "motheridx[100]/I");
   fChain->Branch("ngch", &ngch, "ngch/I");
   fChain->Branch("ncharg", &ncharg, "ncharg/I");
   fChain->Branch("nneu", &nneu, "nneu/I");
   fChain->Branch("kappx", &kappx, "kappx/D");
   fChain->Branch("kappy", &kappy, "kappy/D");
   fChain->Branch("kappz", &kappz, "kappz/D");
   fChain->Branch("kape" , &kape,  "kape/D" );
   fChain->Branch("kampx", &kampx, "kampx/D");
   fChain->Branch("kampy", &kampy, "kampy/D");
   fChain->Branch("kampz", &kampz, "kampz/D");
   fChain->Branch("kame" , &kame , "kame/D" );
   fChain->Branch("mphi" , &mphi , "mphi/D" );
   fChain->Branch("kkm4" , &kkm4 , "kkm4/D" );
   Notify();
}

Bool_t gepep_kk::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void gepep_kk::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t gepep_kk::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef gepep_kk_cxx
