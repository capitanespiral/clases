#include "racionales.hh"
#include <cmath>
using namespace std;

///////////CONSTRUCTORES/////////////

//Constructor por defecto (0/1)
racional::racional():num(0),den(1){}

//Constructor que recibe el numerador y el denominador
racional::racional(long a,long b){
  //Primero saco el maximo comun divisor, para que siempre se usen fracciones irreducibles.
  long c=mcd(a,b);
  a/=c;b/=c;
  num=a;den=b;
}

//genero el racional mas cercano al double (funciona bien hasta dos decimales)(alguna forma de regular en error en función de la cantidad de decimales entregados, (como determinar esta?)
racional::racional(double a){
  long pot=1;
  long b=long(a);
  while((abs(b-a))!=0){
    a*=10;
    pot*=10;
    b=long(a);   
  }
  long c=mcd(long(a),pot);
  a/=c;pot/=c;
  num=a;den=pot;
}

//constructor de copia
racional::racional(const racional &a){
  num=a.num;
  den=a.den;
}

/////////////////////METODOS////////////////////////
//Asignacion con igual (deep copy)
racional & racional::operator=(const racional &a){
  if(this!=&a){
    num=a.num;
    den=a.den;
  }
  return *this;
}

//Para ver el numerador o el denominador
long racional::operator[](int a) const{
  if(a==0) return num;
  else return den;
}

//Para cambiarlo
long & racional::operator[](int a){
  if(a==0) return num;
  else return den;
}

//Entrego su forma decimal
double racional::dec(){
  return double(num)/den;
}

//Calculo su potencia (racional elevado a)
racional racional::pow(int a){
  if(a==0) return racional(1,1);
  else{
    long n=1;
    long d=1;
    for(int i=0;i<abs(a);++i){
      n*=num;
      d*=den;
    }
    if(a>0) return racional(n,d);
    else return racional(d,n);
  }
}

//////////////////////ARITMETICA//////////////////////

//suma
racional operator+(const racional &a,const racional &b){
  racional r(a[0]*b[1]+b[0]*a[1],a[1]*b[1]);
  return r;
}

//multiplicacion
racional operator*(const racional &a,const racional &b){
  racional r(a[0]*b[0],a[1]*b[1]);
  return r;
}

//resta
racional operator-(const racional &a,const racional &b){
  racional r(a[0]*b[1]-b[0]*a[1],a[1]*b[1]);
  return r;
}

//division
racional operator/(const racional &a,const racional &b){
  racional r(a[0]*b[1],a[1]*b[0]);
  return r;
}

///////////////////FUNCIONES VARIAS//////////////////////

//Impresion, de ser negativo siempre el menos en el numerador
ostream & operator<<(ostream &os,const racional &a){
  if(a[0]*a[1]<0)
    os<<-abs(a[0])<<'/'<<abs(a[1]);
  else
    os<<abs(a[0])<<'/'<<abs(a[1]);
  return os;
}

//Maximo comun divisor.
long mcd(long a, long b)
{
  if(b==0)
    return a;
  else
    return mcd(b,a%b);
}
