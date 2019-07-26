#include "cuadratura.h"

double asd(double a){return 3*a*a+2;}

int main(int argc,char* argv[]){
  double a=montecarlo3d(-1.,1.,atoi(argv[1]));
  cout<<a<<endl;
  return 0;

}
