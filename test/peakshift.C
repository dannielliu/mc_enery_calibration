
namespace PeakEstimate{
//private:
  // par[0]: p distribution slope
  // par[1]: gaus mean;
  // par[2]: gaus sigma
  // par[3]: dm(p) E1
  // E1 = par[3];
  // E2 = par[4];
  // p1 = par[5];
  // p2 = par[6];
  // costheta = par[7];
  //double pdfpar[3];
  //double pspar[2];
  //double dmpar[5];

//public:
  //void SetPdfPar(double a, double b, double c)
  //{}
  //void SetPsPar(double mean, double sigma)
  //{
  //  pspar[0] = mean;
	//pspar[1] = sigma;
  //}
  double mresonace = 1.86484;
  double mfinal    = 0.493677;
  double ppdf(double *x, double *par)
  {
    //return par[0]*x[0];
	double res = TMath::Gaus(x[0],par[0],par[1],true);
    //std::cout<<"in ppdf, p1 is "<<x[0]<<", mean is "<<par[0]<<", sigma is "<< par[1] <<", value is "<<res<< std::endl;
	return res;
  }
  double psmear(double *x, double *par)
  { 
    double sigma = 0.01*x[0];
    double res = TMath::Gaus(par[0],x[0],sigma,true);
    //std::cout<<"in psmear, p1 "<<par[0]<<", mean is "<<x[0]<<", sigma is "<< sigma <<", value is "<<res << std::endl;
    return res;
  }
  double deltam(double *x, double *par)
  {
    double m,E1,E2,p1,p2,costheta;
    //double m=1.86484;
	//double mk=0.493677;
    m = mresonace;
	double mk = mfinal;
	p1 = par[0];
    p2 = par[1];
	double p1r = x[0]; // p1 real position, smear to par[0]
	E1 = sqrt(mk*mk + p1r*p1r);
	E2 = sqrt(mk*mk + p2*p2);
    costheta = ((E1+E2)*(E1+E2)-m*m-(p1r*p1r + p2*p2))/(2*p1r*p2);
    //std::cout<<"E1: "<<E1<<"\tE2: "<<E2 <<"\tp1: "<<p1<<"\tp2: "<<p2<<"\tcostheta: "<<costheta<< std::endl;
	if (fabs(costheta)>1) return -1;

    double delta;
    double dp1 = p1 - x[0];
    delta = 1/m*((E1+E2)*p1*dp1/E1-(p1*dp1+p2*dp1*costheta));
    //std::cout<<"in delta m, p1 is "<<p1<<", dp1 is "<<dp1<<", delta m is "<< delta << std::endl;

	return delta;
  }
  double multifgdelta(double *x, double *par)
  {
    //std::cout<<"in multifgdelta"<<std::endl;
	//std::cout<<deltam(x,&par[4])<<std::endl;
    if (deltam(x,&par[4]) <-0.9) return -1;
    return ppdf(x,&par[0])*psmear(x,&par[2])*deltam(x,&par[4]);
  }
  double multifg(double *x, double *par)
  {
    return ppdf(x,&par[0])*psmear(x,&par[2]);
  }
  double dmatp(double *x, double *par)
  {
    //std::cout<<"in dmatp"<<std::endl;
    double dm;

	double pars[6];
	pars[0] = par[0];//p shape
	pars[1] = par[1];
	pars[2] = x[0];// set smear function
	//pars[3] = 0.01*x[0];
	pars[4] = x[0];// set dm
	pars[5] = par[5];

    double p = x[0];
    double sigmap = par[3];
    TF1 fgd("fgd",PeakEstimate::multifgdelta,0,5,6);
    TF1 fg("fg",PeakEstimate::multifg,0,5,4);
    fgd.SetParameters(pars);
	//fgd.SetParameter(2,x[0]);
    //fgd.SetParameter(4,x[0]);//set p1 as a parameter for delta m
    
	fg.SetParameters(pars);
	//std::cout<<multifgdelta(x,par)<<std::endl;
    
	if ( fgd(x[0])<-0.9 ) return -1;
	if ( fgd(x[0]-5*sigmap)<-0.9 ) return -1;
	if ( fgd(x[0]+5*sigmap)<-0.9 ) return -1;

	//int np = 1000;
	//double *xx = new double[np];
	//double *w = new double[np];
	//fgd.CalcGaussLegendreSamplingPoints(np,xx,w,1e-15);
	//double integral1 = fgd.IntegralFast(np,xx,w,p-5*sigmap,p+5*sigmap);
    double integral1 = fgd.Integral(p-5*sigmap,p+5*sigmap);
	//std::cout<<"fgd integral is "<< integral1<< std::endl;
	//fg.CalcGaussLegendreSamplingPoints(np,xx,w,1e-15);
	//double integral2 = fg.IntegralFast(np,xx,w,p-5*sigmap,p+5*sigmap);
	double integral2 = fg.Integral(p-5*sigmap,p+5*sigmap);
	//std::cout<<"fg  integral is "<< integral2<< std::endl;
	//dm = fgd.Integral(p-5*sigmap,p+5*sigmap)/fg.Integral(p-5*sigmap,p+5*sigmap);
    //delete []xx;
	//delete []w;
	return integral1/integral2;
  }
  double multifdm(double *x, double *par)
  {
    //std::cout<<"in multifdm"<<std::endl;
    return ppdf(x,par)*dmatp(x,&par[0]);
  }

  double peakshift(double plow=0, double pup=0)
  {
    // par[0]: p dis mean
    // par[1]: p dis sigma
    // par[2]: gaus mean;
    // par[3]: gaus sigma
    // par[4]: dm(p) E1
    // p1 = par[5];
    // p2 = par[6];
    double par[6] ;
    double mk = 0.493677;
    double &pdismean = par[0] = 1.5;
    double &pdissigma= par[1] = 0.15;
    double &mean  = par[2] = 1.0;
    double &sigma = par[3] = 0.01*mean;
    double &p1    = par[4] = 1.0;
    double &p2    = par[5] = 1.0;
    //double &e1    = par[3] = sqrt(p1*p1+mk*mk);
    //double &e2    = par[4] = sqrt(p2*p2+mk*mk);
 
    TF1 fdm("fdm", PeakEstimate::multifdm,0,2,6);
    TF1 pdis("pdis",PeakEstimate::ppdf,0,2,2);
    //TF1 dm("dm",PeakEstimate::dmatp,0,2,6);
    int Npart = 20;
    for (int i=0;i<Npart;i++){
      p1 = 0.5+2./Npart*i;
      plow = p1-0.05;
      pup = p1+0.05;
      mean = p1;
      sigma = 0.01*mean;
      fdm.SetParameters(par);
      pdis.SetParameters(par);
      //dm.SetParameters(par);
      //std::cout<<"dm at "<<p1<<"GeV is "<<std::endl; 
      //std::cout<<dm(p1)<<std::endl;
        std::cout<<"p1 is "<<p1<<" "<<plow<<" "<<pup<<std::endl;
        std::cout<<"peak shift calculation ..."<<std::endl;
 
      //int np = 1000;
      //double *xx = new double[np];
      //double *w  = new double[np];
      //double ps = fdm.IntegralFast(np,xx,w,plow,pup)/pdis.IntegralFast(np,xx,w,plow,pup);
        double ps = fdm.Integral(plow,pup)/pdis.Integral(plow,pup);
        std::cout<<"peak shift is "<< ps <<std::endl;
        std::cout<<"peak       is "<< ps + 1.86484 <<"\n"<<std::endl;
    }
    TCanvas *c1 = new TCanvas;
    pdis.Draw();
    c1->Print("pdis.pdf");
    return 0;
  }
};
