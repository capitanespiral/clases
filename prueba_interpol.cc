#include "interpol.h"
#include <iostream>
#include <fstream>
using namespace std;
 
int main(int argc, char* argv[]){
	vector<double> v;
	double x;
	ifstream archivo(argv[1]);
	while(archivo>>x) v.push_back(x);
	matriz<double> m(v,2);
	cout<<m<<endl;
	i_pol orale(m);
	matriz<poli<double>> p=orale.spline_p();
	cout<<p<<endl;
	matriz<double> potito=puntos(p,orale,10000);
	file("splinep",potito);
	return 0;
}
