#ifndef AJUSTLIN_H
#define AJUSTLIN_H

#include "matriz.h"

using namespace std;

class ajus_lin{
 private:
  matriz<double> datos;//matriz de los datos
  int cant;//cantidad de pares ordenados
 public:
  ajus_lin(const matriz<double> &);
  matriz<double> min_cuad(const vector<double (*)(double)> &);//entrega coeficientes
};

//entrega pts evaluados en la base, vector de funciones, matriz de coefs, de donde hasta a donde y cuantos puntos.
matriz<double> evalua_coef(const vector<double (*)(double)> &,const matriz<double> &,double,double,int);


#endif
