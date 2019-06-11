#ifndef ODEINT_H
#define ODEINT_H

#include "matriz.h"
#include <cmath>
#include <fstream>
#include <vector>
#include <string>
using namespace std;

//A considerar, vel_velver y pefrl solo funcionan con aceleraciones NO dependientes de la velocidad. ¿Los otros? NI IDEA.
//Corrido de pts en euler-cromer?
//Definir version para solo una edo? por ahora es para edos acopladas
//rk4 adaptativo
//darle besitos a la polilla (all kind of besitos) 1313


//Clase de integradores de ecuaciones diferenciales
//Lo no explicado se explica en el .cc
class odeint{
 private:
  double h,t;//paso (al menos inicial en caso adaptativo) y tiempo inicial
  matriz<double> x,v,q;//Posicion inicial (funcion), velocidad inicial (su primera derivada), masa o carga (o lo que pueda afectar al potencial)
 public:
  odeint(double,double,const matriz<double> &,const matriz<double> &,const matriz<double> &);
  odeint(double,double,const matriz<double> &,const matriz<double> &); 
  matriz<matriz<double>> euler_cromer(matriz<double> (*g)(double,const matriz<double> &,const matriz<double> &,const matriz<double> &),double=0,char='a');
  matriz<matriz<double>> rk4(matriz <double> (*g)(double,const matriz<double> &,const matriz<double> &,const matriz<double> &),double=0,char='a');
  //rk4 para una derivada, unica diferencia es que la funcion recibe una matriz menos.
  matriz<matriz<double>> rk4(matriz <double> (*g)(double,const matriz<double> &,const matriz<double> &),double=0,char='a');
  matriz<matriz<double>> rk4_adap(matriz<double> (*g)(double,const matriz<double> &,const matriz<double> &,const matriz<double> &),double=0,char='a');
  matriz<matriz<double>> vel_velver(matriz <double> (*g)(double,const matriz<double> &,const matriz<double> &,const matriz<double> &),double=0,char='a');
  matriz<matriz<double>> pefrl(matriz<double> (*g)(double,const matriz<double> &,const matriz<double> &,const matriz<double> &),double=0,char='a');
};

int cuentapasos(double,double,double);

matriz<double> desarma(const matriz<matriz<double>> &,int);

matriz<double> desarma1(const matriz<matriz<double>> &,int);

void file(string,const matriz<matriz<double>> &,int,char='a');

#endif
