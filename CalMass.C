void CalMass(double dp1 = 0.0, double dp2 = 0.0)
{
  double p1 = 1.0;
  double p2 = 1.0;
  double p1sq = p1*p1;
  double p2sq = p2*p2;
  double E1 = sqrt(p1sq+TMath::Power(0.493677,2));
  double E2 = sqrt(p2sq+TMath::Power(0.493677,2));
  double m = sqrt((E1+E2)*(E1+E2) - (p1sq+p2sq+2*p1*p2*(-0.25)));
  std::cout<<"m is "<<m <<std::endl;
  dp1 = p1*dp1;
  dp2 = p2*dp2;
  double dm = ((E1+E2)*(p1*dp1/E1+p2*dp2/E2)-(p1*dp1+p2*dp2+p1*dp2*(-0.25)+p2*dp1*(-0.25)))/m;
  std::cout<<"delta m is "<<dm <<std::endl;
}
/*
double CalInvMass(double m1, double p1,double m2, double p2, double costheta, double f1=1.0, double f2=1.0)
{
  double p1,p2,px,py,pz;
  double e1,e2;
  double minv;
  double f1,f2;
  
  p1=TMath::Sqrt(px1*px1+py1*py1+pz1*pz1);
  p2=TMath::Sqrt(px2*px2+py2*py2+pz2*pz2);
  if (n==0) {f1=1; f2=1;}
  else if(n==1) {f1=x[0]; f2=x[0];}
  else if(n==-2) {f1=x[0]; f2=x[1];}
  else {
    int tmpindex;
    tmpindex=(int)((p1-par[0])/(par[1]-par[0])*n);
    if (tmpindex<0) tmpindex =0;
    if (tmpindex>n-1) tmpindex=n-1;
    f1 = x[tmpindex];
    tmpindex=(int)((p2-par[0])/(par[1]-par[0])*n);
    if (tmpindex<0) tmpindex =0;
    if (tmpindex>n-1) tmpindex=n-1;
    f2 = x[tmpindex];
  }
  e1=TMath::Sqrt(m1*m1+f1*f1*p1*p1);
  e2=TMath::Sqrt(m2*m2+f2*f2*p2*p2);
  px=f1*px1+f2*px2;
  py=f1*py1+f2*py2;
  pz=f1*pz1+f2*pz2;
  minv = TMath::Sqrt((e1+e2)*(e1+e2)-(px*px+py*py+pz*pz));
  
  return minv;
}
*/

