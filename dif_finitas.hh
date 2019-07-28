#ifndef DIF_FINITAS_H
#define DIF_FINITAS_H

#include "matriz.hh"

using namespace std;

//Forward diferencia finita R a R
template <class T>
T diff_forw(T (*f)(T),T a,double h=1e-8){
  return ((*f)(a+h)-(*f)(a))/h;
}

//Backward diff R a R
template <class T>
T diff_back(T (*f)(T),T a,double h=1e-8){
  return ((*f)(a)-(*f)(a-h))/h;
}

//Central diff R a R
template<class T>
T diff_cen(T (*f)(T),T a,double h=1e-8){
  return ((*f)(a+h)-(*f)(a-h))/(2*h);
}

//Central diff general (derivada parcial), se espera el a evaluar como columna. Funcion de n a m, valor donde evaluo, respecto a que coordenada, y que coordenada de la funcion.
template<class T>
T diff_cen(matriz<T> (*f)(const matriz<T> &),const matriz<T> &a,int var,int varf,double h=1e-8){
  matriz<T> temp(a);temp(var,0)=temp(var,0)+h;
  T acum=(*f)(temp)(varf,0);
  temp(var,0)=temp(var,0)-2*h;
  acum=acum-(*f)(temp)(varf,0);
  return acum/(2*h);
}

//Jacobiana de una funcion, asumo vectores columna
template<class T>
matriz<T> jacob(matriz<T> (*f)(const matriz<T> &),const matriz<T> &a,double h=1e-8){
  matriz<T> temp=(*f)(a);
  matriz<T> res(temp.fila(),a.fila());
  for(int i=0;i<res.fila();++i){
    for(int j=0;j<res.colu();++j){
      res(i,j)=diff_cen(f,a,j,i,h);
    }
  }
  return res;
}


#endif 
