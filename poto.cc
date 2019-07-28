#include "ceros.hh"
#include "dif_finitas.hh"

//funcion x+1,y-x,z-2x
matriz<double> f(const matriz<double> &m){
  matriz<double> res(2,1);
  res(0,0)=0.5*m(0,0)*m(0,0)*m(0,0)+m(1,0)*m(1,0)*m(1,0);res(1,0)=3*m(0,0)+5*m(1,0);
  return res;
}

//primera fila
double f11(const matriz<double> &m){
  return 1;
}

double f12(const matriz<double> &m){
  return 0;
}

double f21(const matriz<double> &m){
  return -1;
}

double f31(const matriz<double> &m){
  return -2;
}


int main(int argc,char* argv[]){
  vector<double (*)(const matriz<double> &)> v={f11,f12,f12,f21,f11,f12,f31,f12,f11};
  vector<double> w={1,1};
  matriz<double> h(w,1);
  matriz<double (*)(const matriz<double> &)> m(3,v);
  //matriz<double> n=newton_n(f,m,h);
  // cout<<n<<endl;
  matriz<double> n=newton_n(f,h);
  cout<<n<<endl;
  return 0;

}
