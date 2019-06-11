#ifndef DUALES_H
#define DUALES_H

#include <iostream>
#include <cmath>

using namespace std;

//Clase de duales con templates
template <class W>
class dual{
 private:
  double r,e;//Miembros privados, parte real y parte dual
 public:
  dual():r(0),e(1){} //constructor por defecto
  dual(W x,W y=1):r(x),e(y){} //constructor que recibe la parte real y la dual (de no recibir dual la hace 1) 	
  dual(const dual &d){//constructor de copia
    r=d.r;
    e=d.e;
  }
  dual & operator=(const dual& d){//asignador de igual (deep copy)
    if(this!=&d){
      r=d.r;
      e=d.e;
    }
    return *this;
  }
  W operator[](int i) const{//Devolver parte real o dual para ver
    if(i==0) return r;
    else return e;
  } 
  W & operator[](int i){//Lo mismo para sobreescribir
    if(i==0) return r;
    else return e;
  }
};

////////////////ARITMETICA///////////////////////7

//SUMA
template <class W>
dual<W> operator+(const dual<W> &x,const dual<W> &y){
  dual<W> d(x[0]+y[0],x[1]+y[1]);
  return d;
}

//RESTA
template <class W>
dual<W> operator-(const dual<W> &x,const dual<W> &y){
  dual<W> d(x[0]-y[0],x[1]-y[1]);
  return d;
}

//MULTIPLICACION
template <class W>
dual<W> operator*(const dual<W> &x,const dual<W> &y){
  dual<W> d(x[0]*y[0],x[1]*y[0]+x[0]*y[1]);
  return d;
}

//DIVISION
template <class W>
dual<W> operator/(const dual<W> &x,const dual<W> &y){
  dual<W> d(x[0]/y[0],(x[1]*y[0]-x[0]*y[1])/(y[0]*y[0]));
  return d;
}

/////////////////////////FUNCIONES MAS COMPLEJAS//////////////////////

//SENO
template <class W>
dual<W> sin(const dual<W> &d){
  dual<W> f(sin(d[0]),cos(d[0])*d[1]);
  return f;
}

//SENO HIPERBOLICO
template <class W>
dual<W> sinh(const dual<W> &d){
  dual<W> f(sinh(d[0]),cosh(d[0])*d[1]);
  return f;
}

//COSENO
template <class W>
dual<W> cos(const dual<W> &d){
  dual<W> f(cos(d[0]),-sin(d[0])*d[1]);
  return f;
}

//COSENO HIPERBOLICO
template <class W>
dual<W> cosh(const dual<W> &d){
  dual<W> f(cosh(d[0]),sinh(d[0])*d[1]);
  return f;
}

//LOGARITMO
template <class W>
dual<W> log(const dual<W> &d){
  dual<W> f(log(d[0]),d[1]/d[0]);
  return f;
}

//EXPONENCIAL
template <class W>
dual<W> exp(const dual<W> &d){
  dual<W> f(exp(d[0]),exp(d[0])*d[1]);
  return f;
}

/////////////////////////VARIOS////////////////////////

//IMPRESION
template <class W>
ostream & operator<<(ostream &os,const dual<W> &d)
{
  os<<"{"<<d[0]<<" , "<<d[1]<<"}";
  return os;
}

#endif
