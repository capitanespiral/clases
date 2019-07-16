#include "random.h"
#include <fstream>
#include <cmath>

double f(double x){
	return sin(x)+1;
}

int main(int argc,char* argv[]){
	double min=0;
	double max=9;
	rdom andale;
	ofstream archivo("orale.dat");
	matriz<double> res=distr(f,0,2,-10,10,2000);
	for(int i=0;i<res.fila();++i){
		archivo<<res(i,0)<<' '<<f(res(i,0))<<endl;
}
return 0;
}
