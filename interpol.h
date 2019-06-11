#ifndef INTERPOL_H
#define INTERPOL_H

#include "matriz.h"
#include "poli.h"
#include <fstream>
#include <string>

using namespace std;

//Clase de interpoladores, lo no explicado se explica en el cc.
class i_pol{
 private:
  matriz<double> datos;//matriz de los datos
  int cant;//cantidad de pares ordenados
 public:
  i_pol(const matriz<double> &);
  poli<double> lagrange();
  matriz<poli<double>> spline();
  matriz<poli<double>> spline_p();
  friend matriz<double> puntos(const matriz<poli<double>> &,i_pol,int);
  friend matriz<double> puntos(poli<double>,i_pol,int);
};

void file(string,const matriz<double> &);


#endif
