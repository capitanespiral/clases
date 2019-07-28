#ifndef CUADRATURA_H
#define CUADRATURA_H

#include "poli.hh"
#include "ceros.hh"
#include "matriz.hh"
#include "random.hh"
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

//Gauss-legendre para dos integrales, con todos los limites definidos (integrando respecto a "y" primero, osea segundo parametro el c d y primer parametro a b.)
template <class T>
T gauss_leg2(T (*f)(T,T),T a,T b,T c,T d,int N){
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
  for(int i=0;i<N;++i){
    for(int j=0;j<N;++j){
      res=res+pesos(i,0)*pesos(j,0)*(*f)((b-a)*0.5*raices(i,0)+(a+b)*0.5,(d-c)*0.5*raices(j,0)+(c+d)*0.5);
    }
  }
  return 0.25*(b-a)*(d-c)*res;
}


//Doble integral cuando los limites de la interior son funciones. Por defecto, se integra y primero (o sea, el segundo valor de f), de manera que las funciones c y d estan en funcion de x.
template <class T>
T gauss_leg2(T (*f)(T,T),T a,T b,T (*c)(T),T (*d)(T),int N){
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
  for(int i=0;i<N;++i){
    res=res+pesos(i,0)*gauss_aux(f,c,d,pesos,raices,0.5*((b-a)*raices(i,0)+a+b));
  }
  return 0.5*(b-a)*res;
}

//Funcion auxiliar gauss_2 general
template <class T>
T gauss_aux(T (*f)(T,T),T (*c)(T),T (*d)(T),const matriz<T> &pesos,const matriz<T> &raices,T x){
  T acum=T(0);
  for(int i=0;i<pesos.fila();++i){
    acum=acum+pesos(i,0)*f(x,0.5*((d(x)-c(x))*raices(i,0)+c(x)+d(x)));
  }
  return 0.5*(d(x)-c(x))*acum;
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
  rdom random;
  T max=T(0);T cualquiera;
  bool neg=false;
  for(int i=0;i<(b-a)*densidad;++i){
    cualquiera=(*f)(((b-a)*i/densidad)+a);
    if(cualquiera>max) max=cualquiera;
    if(cualquiera<0) neg=true;
  }
  if(neg) cout<<"Cuidado, valores menores que cero"<<endl;
  T x,y;int cont=0;
  for(int i=0;i<puntos;++i){
    x=random.drand(b,a);y=random.drand(max+1,0);
    if(y<=(*f)(x)) cont+=1;
    if(i%1000==0) random.upd();
  }
  return cont*(b-a)*(max+1)/puntos;
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
