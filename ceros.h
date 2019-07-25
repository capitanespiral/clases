#ifndef CEROS_H
#define CEROS_H

#include "poli.h"
#include "matriz.h"
#include "complejos.h"
#include <cmath>

typedef unsigned int uint;
//Funciones para sacar ceros de funciones

//metodo de newton, recibe funcion, derivada, semilla y error.
template <class T>
T newton(T (*f)(T),T (*g)(T),T seed=T(1),double error=1e-10){
  T err=T(1);
  T xn=seed;
  T xn1;
  int i=0;
  while(err>error && err!=0){
    xn1=xn-(*f)(xn)/(*g)(xn);
    err=abs(xn-xn1);
    xn=xn1;
    ++i;
    if(i>100){cout<<"No convergió"<<endl;break;}
  }
  return xn;
}

//newton con polinomios (es necesario dar una semilla con parte imaginaria si quiero raiz compleja)
template <class T,class U>
T newton(const poli<U> &p,T seed=T(1),double error=1e-10){
  T err(1);
  T xn=seed;
  T xn1;
  poli<U> q=p.deriv();
  int i=0;
  while(err>error && err!=0){
    xn1=xn-(eval(p,xn)/eval(q,xn));
    err=abs(xn-xn1);
    xn=xn1;
    ++i;
    if(i>1000) {cout<<"No convergió"<<endl;break;}
  }
  return xn;
}

//biseccion, recibo funcion, intervalo [a,b] y error
template <class T>
T biseccion(T (*f)(T),T a=T(-10),T b=T(10),double error=1e-10){
  if((*f)(a)*(*f)(b)<0){
    T err(1);
    T medio;
    while(err!=0 && abs(err)>error){
      medio=(a+b)/2.;
      err=(*f)(medio);
      if(err*(*f)(a)<0) b=medio;
      else a=medio;
    }
    return medio;
  } 
  else{
    cout<<"Mal definición de intervalo"<<endl;
    return 0;
  }
}

//biseccion para polinomios
template <class T,class U>
T biseccion(const poli<U> &p,T a=T(-10),T b=T(10),double error=1e-10){
  if(eval(p,a)*eval(p,b)<0){
    T err(1);
    T medio;
    while(err!=0 && abs(err)>error){
      medio=(a+b)/2;
      err=p(medio);
      if(err*eval(p,a)<0) b=medio;
      else a=medio;
    }
    return medio;
  } 
  else{
    cout<<"Mal definición de intervalo"<<endl;
    return 0;
  }
}

//durand (solo funcionara si T es complejo)
template <class T,class U>
matriz<T> durand(const poli<U> &p,matriz<T> &m,double error=1e-10){
  complejo<double> uno(1,0);
  if(m.fila()==1 && m.colu()==1){
    matriz<complejo<double>> start(p.grad(),1,uno);
    for(int i=0;i<p.grad();++i) start(i,0)=raiz(uno,p.grad(),i+1);
    m=start;
  }
  complejo<double> unoo(uno);
  poli<double> q=p/p[p.grad()];
  double err;
  bool flag=true;
  int k=0;
  while(true){
    for(int i=0;i<m.fila();++i){
      for(int j=0;j<m.fila();++j){
	if(i!=j) uno=uno*(m(i,0)-m(j,0));
      }
      m(i,0)=m(i,0)-(eval(q,m(i,0))/uno);
      uno=unoo;
    }
    for(int i=0;i<m.fila();++i){
      err=eval(q,m(i,0)).mod();
      if(err>error){flag=false;break;}
    }
    if(flag) return m;
    flag=true;
    ++k;
    if(k>1000){cout<<"No convergió"<<endl;break;}
  }
  return m;
}

//Decidir intervalo (específico para legendre), fijo a y devuelvo un b tentativo (que cumpla la diferencia de signos). Trozo son los trozos en que corto en el intervalo entregado (inicialmente). 
template <class T,class U>
T int_legendre(const poli<U> &p,T a, T b,int trozo = 3){
  //Si sirve, de inmediato lo devuelvo
  if(eval(p,a)*eval(p,b)<0) return b;
  else{
    double temp;
    bool flag=true;
    int j=0;
    //Voy cortando el intervalo y evaluando en cada corte.
    while(flag && j<100){
      for(int i=trozo;i>=1;--i){
	temp=(i*(b-a)/double(trozo)+a);
	if(eval(p,a)*eval(p,temp)<0) {flag=!flag;break;}
      }
      //Si ninguno sirvio realizo cortes mas pequeños
      ++j;trozo+=j;
    }
    //Medida de seguridad si nunca funciono
    if(j>=100) cout<<"Cuidado, intervalo talvez no sirva"<<endl;
    return temp;
  }
}

//Encontrar todos las raices de un polinomio modo newton y modo biseccion, usare a de semilla si es metodo de newton. Biseccion no está bien hecho para complejos.
template <class T,class U>
matriz<T> raices_pol(const poli<U> &p,char m,T a=T(-1),T b=T(1),double error=1e-10){
  T llenar = T(1);
  matriz<T> res(p.grad(),1,llenar);
  if(m =='n'){
    poli<T> aux(1,llenar);
    poli<T> temp=cambia_tipo(p,llenar);
    for(int i=0;i<p.grad();++i){
      llenar=newton(temp,a,error);
      res(i,0)=llenar;
      aux[0]=-1*llenar;
      temp=temp/aux;
    }
    return res;
  }
  else if(m == 'b'){
    poli<T> temp(p);
    poli<T> aux(1,llenar);
    for(int i=0;i<p.grad();++i){
      if(temp.grad()%2!=0){
	a=-1;b=1;
      }
      else{
	a=-1;b=int_legendre(temp,a,double(1));
      }
      llenar=biseccion(temp,a,b,error);
      res(i,0)=llenar;
      aux[0]=-1*llenar;
      temp=temp/aux;
    }
    return res;
  }
  else{cout<<"Opcion no existente"<<endl;return res;}
}

#endif

