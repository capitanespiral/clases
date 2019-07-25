#ifndef COMPLEJOS_
#define COMPLEJOS_

#include <iostream>
#include <cmath>
using namespace std;

template <class T>
class complejo{
 private:
  T r;
  T i;
 public:
  complejo():r(0),i(0){}
  complejo(T,T);
  complejo(T);
  complejo(T,T,char);
  complejo(const complejo &);
  complejo& operator=(const complejo &);
  T operator[](int) const;
  T& operator[](int);
  complejo & operator++();
  complejo & operator++(int);
  complejo operator-() const &;
  complejo conj() const &;
  T mod() const &;
  complejo inv() const &;
  T arg() const &;
  operator double() const {return (*this).mod();}
};

template <class T>
complejo<T>::complejo(T a,T b){
  r=a;
  i=b;
}

template <class T>
complejo<T>::complejo(T a){
  r=a;
  i=0;
}

//constructor recibiendo r y arg (a y b)
template <class T>
complejo<T>::complejo(T a,T b, char){
  r=a*cos(b);
  i=a*sin(b);
  if(abs(r)<1e-15) r=0;
  if(abs(i)<1e-15) i=0;
}


template <class T>
complejo<T>::complejo(const complejo<T> &c){
  r=c.r;
  i=c.i;
}

template <class T>
complejo<T>& complejo<T>::operator=(const complejo<T> &c){
  if(this!=&c){
    r=c.r;
    i=c.i;
  }
  return (*this);
}
  
template <class T> 
T complejo<T>::operator[](int h) const{
  if(h==0) return r;
  else return i;
}

template <class T>
T& complejo<T>::operator[](int h){
  if(h==0) return r;
  else return i;
}

template <class T>
complejo<T> & complejo<T>::operator++(){
  r=r+1;
  return *this;
}

template <class T>
complejo<T> & complejo<T>::operator++(int){
  i=i+1;
  return *this;
}

//argumento (angulo) del complejo
template <class T>
T complejo<T>::arg() const &{
  if(r>0){
    if(i>=0) return atan(i/r);
    else return atan(i/r)+2*M_PI;
  }
  else if(r<0) return atan(i/r)+M_PI;
  else {
    if(i>0) return M_PI/2;
    else if(i<0) return 3*M_PI/2;
    else {cout<<"Mal definido el angulo para el complejo nulo"<<endl;return 0;}
  }
}

////aritmetica basica////
//suma
template <class T>
complejo<T> operator+(const complejo<T> &z,const complejo<T> &w){
  complejo<T> res(z[0]+w[0],z[1]+w[1]);
  return res;
}

//poner el menos
template <class T>
complejo<T> complejo<T>::operator-() const &{
  complejo<T> res(-1*((*this)[0]),-1*((*this)[1]));
  return res;
}

//restar
template <class T>
complejo<T> operator-(const complejo<T> &z,const complejo<T> &w){
  return z+(-w);
}

//conjugado
template <class T>
complejo<T> complejo<T>::conj() const &{
  complejo<T> res((*this)[0],-1*((*this)[1]));
  return res;
}

//Multiplicacion
template <class T>
complejo<T> operator*(const complejo<T> &z,const complejo<T> &w){
  complejo<T> res(z[0]*w[0]-z[1]*w[1],z[0]*w[1]+z[1]*w[0]);
  return res;
}

//Modulo
template <class T>
T complejo<T>::mod() const &{
  return sqrt((*this)[0]*(*this)[0]+(*this)[1]*(*this)[1]);
}

//Multiplicacion por escalar
template <class T,class U>
complejo<T> operator*(U a,const complejo<T> &z){
  complejo<T> res(T(a)*z[0],T(a)*z[1]);
  return res;
}

template <class T,class U>
complejo<T> operator*(const complejo<T> &z, U a){
  complejo<T> res(T(a)*z[0],T(a)*z[1]);
  return res;
}

//Division por escalar
template<class T,class U>
complejo<T> operator/(const complejo<T> &z,U a){
  return z*(1/T(a));
}

//Inverso de un complejo
template <class T>
complejo<T> complejo<T>::inv() const &{
  complejo<T> res=((*this).conj())/((*this).mod()*(*this).mod());
  return res;
}

//division de complejos
template <class T>
complejo<T> operator/(const complejo<T> &z,const complejo<T> &w){
  return z*(w.inv());
}

//Raiz de complejo, recibe complejo, indice de la raiz y cual entregar
template <class T>
complejo<T> raiz(const complejo<T> &c,int a,int b){
  if (b<=a){
    complejo<T> res(pow(c.mod(),1/double(a)),(c.arg()+(b-1)*2*M_PI)/a,'w');
    return res;
  }
  else{
    cout<<"No tiene sentido raiz pedida"<<endl;
    complejo<T> res;
    return res;
  }
}

//impresion
template <class T>
ostream & operator<<(ostream &os, const complejo<T> &z){
  os<<'('<<z[0]<<" , "<<z[1]<<')';
  return os;
}


//////operadores de comparacion///////
template <class T>
bool operator<(const complejo<T> &z,const complejo<T> &w){
  return z.mod()<w.mod();
}

template <class T>
bool operator>(const complejo<T> &z,const complejo<T> &w){
  return z.mod()>w.mod();
}

template <class T>
bool operator==(const complejo<T> &z,const complejo<T> &w){
  return z.mod()==w.mod();
}

template <class T>
bool operator!=(const complejo<T> &z,const complejo<T> &w){
  return z.mod()!=w.mod();
}
#endif
