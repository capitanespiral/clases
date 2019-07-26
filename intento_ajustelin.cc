#include "ajustelineal.h"
#include "random.h"
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

double f(double a){
  return 0.3*a*a+7;
}
double g(double a){
  return a*a;
}

double t(double a){
  return a;
}

double r(double a){return 1;}

int main(){
  vector<double> v={1,4,4,7,6,2,4,0,2,9,7,1,8,1};
  matriz<double> asd(v,2);
  cout<<asd<<endl;
  /* matriz<double> m=dat_func(sin,0,15,300);
  env_archivo("original.dat",m);
  m=noise(m,0.5,1);
  env_archivo("noise.dat",m);
  matriz<double> n=rec_archivo("noise.dat",2,2.);*/
  ajus_lin poto(asd);
  vector<double (*)(double)> w={t,g};
  asd=poto.min_cuad(w);
  cout<<asd<<endl;
  /*n=evalua_coef(w,m,0,15,300);
    env_archivo("ajuste.dat",n);*/
  
  return 0;
}
