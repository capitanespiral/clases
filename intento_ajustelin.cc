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
  matriz<double> m=dat_func(sin,0,15,300);
  env_archivo("original.dat",m);
  m=noise(m,0.5,1);
  env_archivo("noise.dat",m);
  matriz<double> n=rec_archivo("noise.dat",2,2.);
  ajus_lin poto(n);
  vector<double (*)(double)> w={sin,cos};
  m=poto.min_cuad(w);
  cout<<m<<endl;
  n=evalua_coef(w,m,0,15,300);
  env_archivo("ajuste.dat",n);
  
  return 0;
}
