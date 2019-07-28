#include "cuadratura.hh"
#include "dif_finitas.hh"

double asdf(double a,double b){
  return sqrt(1-a*a-b*b);
}
double x(double a){
  return sqrt(1-a*a);
}
double xx(double a){return -sqrt(1-a*a);}

int main(int argc,char* argv[]){
  double e=exp(1);
  double a=gauss_leg2(asdf,-1.,1.,xx,x,20);
  cout<<2*a<<endl;

  return 0;

}
