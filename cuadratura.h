#ifndef CUADRATURA_H
#define CUADRATURA_H

#include "poli.h"
#include "ceros.h"
#include "matriz.h"
#include "random.h"
using namespace std;

//integracion del trapezoide, recibe la funcion, de donde hasta donde integrar (a,b) y en cuanto subdividir el intervalo
template <class T>
T trapezoide(T (*f)(T),T a,T b,int inicial=100000,int paso=10000,double error=1e-6){
  T res=T(0);
  T res1=T(1);
  int N=inicial;
  while(abs(res-res1)>error && N<1000000){
    res1=res;
    res=0;
    for(int i=0;i<N;++i) res=res+(*f)(a+i*(b-a)/T(N));
    res=((b-a)/T(N))*((*f)(a)*0.5+(*f)(b)*0.5+res);
    N+=paso;
  }
  return res;
}

//Integracion gauss-legendre, recibe lo mismo pero el N es ahora que puntos tomar. Tambien valor auxiliar por el valor extra "S" que tiene la funcion pedida
template <class T>
T gauss_leg(T (*f)(T,T),T a, T b,int N,T aux){
  poli<T> p;
  p=p.legendre(N);
  matriz<T> raices=raices_pol(p,'n',T(-1));
  for(int i=0;i<raices.fila();++i)
    if(abs(raices(i,0))<1e-15) raices(i,0)=0;
  p=p.deriv();
  matriz<T> pesos(N,1);
  for(int i=0;i<pesos.fila();++i)
    pesos(i,0)=2/((1-raices(i,0)*raices(i,0))*p(raices(i,0))*p(raices(i,0)));
  T res=0;
  for(int i=0;i<N;++i)
    res=res+pesos(i,0)*(*f)((b-a)*0.5*raices(i,0)+(a+b)*0.5,aux);
  return (b-a)*0.5*res;
}

//Metodo de romberg, recibe funcion, valores donde integrar, pasos maximos y error esperado
template <class T>
T romberg(T (*f)(T),T a,T b,int max=20,double error = 1e-8){
  T acum=0;
  T hm=1;
  bool flag=true;
  int I;
  matriz<T> m(max,max);
  T hn=b-a;
  m(0,0)=0.5*hn*((*f)(a)+(*f)(b));
  for(int i=1;i<max;++i){
    hn=hn/2.;
    for(int k=1;k<=pow(2,i-1);++k) acum=acum+(*f)(a+(2*k-1)*hn);
    m(i,0)=0.5*m(i-1,0)+hn*acum;
    for(int j=1;j<=i;++j){
      hm=hm*4;
      m(i,j)=m(i,j-1)+(m(i,j-1)-m(i-1,j-1))/(hm-1);
    }
    acum=0;hm=1;
    if(abs(m(i,i)-m(i,i-1))<error) {flag=!flag;I=i;break;}
  }
  if(flag) {cout<<"No se alcanzó el error deseado"<<endl;return m(0,0);}
  return m(I,I);
}

//Integracion montecarlo 2-D, densidad cuantos puntos en una unidad
//Talvez encontrar max o min con difs finitas
//Pa ajustarlo a cualquier problema, cambiar el penultimo if y la funcion que recibe
template<class T>
T montecarlo(T (*f)(T),T a,T b,int puntos,double densidad=500){
  matriz<T> m((b-a)*densidad,1);
  rdom random;
  bool neg=false;
  for(int i=0;i<m.fila();++i){
    m(i,0)=(*f)(((b-a)*i/densidad)+a);
    if(m(i,0)<0) neg=true;
  }
  if(neg) cout<<"Cuidado, valores menores que cero"<<endl;
  T maxx=max(m);
  T x,y;int cont=0;
  for(int i=0;i<puntos;++i){
    x=random.drand(b,a);y=random.drand(maxx+1,0);
    if(y<=(*f)(x)) cont+=1;
    if(i%1000==0) random.upd();
  }
  return cont*(b-a)*(maxx+1)/puntos;
}

//Le dare 2 puntos, defino un cubo
template<class T>
T montecarlo3d(T min,T max,int puntos){
  rdom random;
  T x,y,z;
  int cont=0;
  for(int i=0;i<puntos;++i){
    x=random.drand(max,min);y=random.drand(max,min);z=random.drand(max,min);
    //Condiciones pa contar o no
    if(x*x+y*y+z*z<=1) cont+=1;
    if(i%1000==0) random.upd();
  }
  return cont*(max-min)*(max-min)*(max-min)/puntos;
}

#endif
