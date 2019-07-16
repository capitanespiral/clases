#ifndef POLI_
#define POLI_

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>

using namespace std;
typedef unsigned int uint;

//Clase de polinomios con templates
template <class T>
class poli{
 private://miembros privados como el grado y un vector de los coeficientes
  uint grado;
  vector<T> coef;//guardados de menor a mayor grado
 public://estos se explicaran en su definicion
  poli();
  poli(uint, vector<T>);
  poli(uint,T=0,T=1);
  poli(const poli&);
  poli& operator=(const poli&);
  T operator[](int) const;
  T& operator[](int);
  uint grad() const {return grado;} //devuelvo grado
  poli operator+(const poli&);
  poli operator-();
  poli operator-(poli);
  void b_grad();
  T operator()(T);
  poli deriv();
  T deriv(T);
  poli integ(T=0);
  T integ(T,T);
  poli legendre(uint);
};

/////////////////////////CONSTRUCTORES////////////////////////////

//constructor por defecto, polinomio nulo.
template <class T>
poli<T>::poli(){
  grado=0;
  vector<T> v={0};
  coef=v;
}

//Constructor que recibe grado y vector con coeficientes. Se esperan los elementos de mayor a menor grado, pero se almacenan de menor a mayor grado. NO CONSIDERA QUE ENTREGUE GRADO MENOR A CANTIDAD DE ELEMENTOS DEL VECTOR!!!!!
template <class T>
poli<T>::poli(uint i, vector<T> v){//grado de mayor a menor en el vector entregado
  grado=i;
  //Si se entrega vector nulo, guardo un cero.
  if(v.size()==0)
    v.push_back(0);
  //Si entregue menos elementos que lo que el grado exige, guardo los ceros necesarios.
  if(i+1>v.size()){
    for(uint j=0;j<=i-v.size();++j)
      coef.push_back(0);
  }
  //Guardo el resto de los elementos del vector en coef.
  for(uint j=1;j<=v.size();++j)
    coef.push_back(v[v.size()-j]);
  //Verifico que el grado mayor no tenga coeficiente cero
  if(coef.back()==0 && coef.size()>1){
    //Lo itero hasta encontrar el coeficiente de mayor grado no nulo
    for(uint j=0;j<coef.size();++j){
      if(coef.back()==0 && coef.size()>1){
	int s=coef.size()-1;
	coef.erase(coef.begin()+s);
	grado-=1;
      }
      else
	break;
    }
  }
}

//Constructor que recibe el grado, y permite rellenar todos los coeficientes con un valor 'b' o solo el de mayor grado y todos los demas cero (primera opcion por default).
template <class T>
poli<T>::poli(uint a,T b,T c){
  if(c==1){
    grado=a;
    vector<T> v(a+1,b);
    coef=v;
  }
  else{
    grado=a;
    vector<T> v(a+1,0);
    v.push_back(b);
    coef=v;
  }
}


//Constructor de copia
template <class T>
poli<T>::poli(const poli<T> &p){
  grado=p.grado;
  coef=p.coef;
}
 
//////////////METODOS/FUNCIONES INTERNAS/////////////////

//Operador de asignacion (deep copy)
template <class T>
poli<T>& poli<T>:: operator=(const poli<T> &p){
  if(this!=&p){
    grado=p.grado;
    coef=p.coef;
  }
  return *this;
}

//Para ver los coeficientes, de menor a mayor.   
template <class T>
T poli<T>:: operator[](int i) const {
  return coef[i];
}

//Para cambiarlos de ser necesario.
template <class T>
T& poli<T>:: operator[](int i){
	return coef[i];
}

//Suma de polinomios, tres casos, si tienen misma dimension, si uno es mayor o si el otro lo es.
template <class T>
poli<T> poli<T>::operator+(const poli<T> &p){
  if(this->grado==p.grado){
    vector<T> v;
    for(uint i=0;i<=p.grado;++i)
      v.push_back((*this)[i]+p[i]);
    reverse(v.begin(),v.end());
    poli<T> h(v.size()-1,v);
    return h;
  }
  else if(this->grado>p.grado){
    vector<T> v;
    for(uint i=0;i<=p.grado;++i)
      v.push_back((*this)[i]+p[i]);
    for(uint i=p.grado+1;i<=this->grado;++i)
      v.push_back((*this)[i]);
    reverse(v.begin(),v.end());
    poli<T> h(v.size()-1,v);
    return h;
  }
  else{
    vector<T> v;
    for(uint i=0;i<=this->grado;++i)
      v.push_back((*this)[i]+p[i]);
    for(uint i=this->grado+1;i<=p.grado;++i)
      v.push_back(p[i]);
    reverse(v.begin(),v.end()); 
    poli<T> h(v.size()-1,v);
    return h;
  }
}

//Cambiar el signo del polinomio
template <class T>
poli<T> poli<T>::operator-(){
  poli<T> h(this->grado);
  for(uint i=0;i<=this->grado;++i)
    h[i]=-(*this)[i];
  return h;
}

//Resta entre polinomios (por eso se sobrecargo el operador anterior)
template <class T>
poli<T> poli<T>::operator-(poli<T> p){
  return (*this)+(-p);
}

//Funcion que baja un grado del polinomio (eliminando ese coeficiente tambien). Se usa para facilitar la division de polinomios. 
template <class T>
void poli<T>::b_grad(){
  coef.pop_back();
  grado-=1;
}

//Evaluacion del polinomio en algun valor 'a'
template <class T>
T poli<T>::operator()(T a){
  T acum=0;
  for(uint i=0;i<=this->grado;++i)
    acum+=(*this)[i]*pow(a,i);
  return acum;
}

//Derivada, entrega nuevo polinomio
template <class T>
poli<T> poli<T>::deriv(){
  if(grado==0){
    vector<T> v={0};
    poli<T> p(0,v);
    return p;
  }
  else{
    vector<T> v;
    for(uint i=1;i<=grado;++i){
      v.push_back(coef[i]*i);
    }
    reverse(v.begin(),v.end());
    poli<T> p(grado-1,v);
    return p;
  }
}

//Derivada evaluada en algun valor 'a'
template <class T>
T poli<T>::deriv(T a){
  poli<T> p=this->deriv();
  return p(a);
}

//Integral, devuelve funcion y permite escoger constante de integracion (0 por defecto).
template <class T>
poli<T> poli<T>::integ(T a){
  vector<T> v;
  for(uint i=0;i<=grado+1;++i){
    if(i==0)
      v.push_back(a);
    else
      v.push_back(coef[i-1]/(i)); 
  }
  reverse(v.begin(),v.end());
  poli<T> p(grado+1,v);
  return p;
}

//Integral evaluada entre el punto a y b.
template <class T>
T poli<T>::integ(T a,T b){
  poli<T> p=this->integ();
  return p(b)-p(a);
}

//////////////////////////////////ARITMETICA NO INTERNA////////////////////////////////////////////////

//Multiplicacion entre polinomios
template <class T>
poli<T> operator*(const poli<T> &p, const poli<T> &q){
  vector<T> v(p.grad()+q.grad()+1,0);
  for(uint i=0;i<=p.grad();++i){
    for(uint j=0;j<=q.grad();++j){
      v[i+j]+=p[i]*q[j];
    }
  }
  reverse(v.begin(),v.end());
  poli<T> h(p.grad()+q.grad(),v);
  return h;
}

//Multiplicacion entre polinomio y escalar
template <class T,class U>
poli<T> operator*(const poli<T> &p,U a){
  vector<T> v(p.grad()+1,0);
  for(uint i=0;i<=p.grad();++i)
    v[i]+=T(a)*p[i];
  reverse(v.begin(),v.end());
  poli<T> h(p.grad(),v);
  return h;
}

//Lo mismo para conmutatividad
template <class T,class U>
  poli<T> operator*(U a,const poli<T> &p){
  vector<T> v(p.grad()+1,0);
  for(uint i=0;i<=p.grad();++i)
    v[i]+=T(a)*p[i];
  reverse(v.begin(),v.end());
  poli<T> h(p.grad(),v);
  return h;
}

//Division entre polinomios, p/q, no se usa directamente, se accede a traves de los operadores / y %. Por defecto entrega cuociente.
template <class T>
poli<T> div(const poli<T> &p,const poli<T>&q, char c='q'){
  vector<T> cuo;
  if(q.grad()>p.grad()){
    cout<<"El divisor no puede tener mayor grado que el dividendo"<<endl;
    poli<T> p;
    return p;
  }
  poli<T> aux1(p);poli<T> aux2(q);
  for(uint i=0;i<=p.grad()-q.grad();++i){
    if(aux1.grad()+i!=p.grad()){
      cuo.push_back(0);continue;
    }
    T a = aux1[aux1.grad()]/q[q.grad()];
    cuo.push_back(a);
    poli<T> bla(aux1.grad()-aux2.grad(),a,0);
    aux1.b_grad();aux2.b_grad();
    aux1=aux1-(aux2*bla);
    aux2=q;
  }
  if(c=='q'){
    poli<T> Cuo(cuo.size()-1,cuo);
    return Cuo;
  }
  else return aux1;
}

//Entrega cuociente de p/q
template <class T>
poli<T> operator/(const poli<T> &p,const poli<T> &q){
  return div(p,q);
}

//Entrega resto de la division entre p y q
template <class T>
poli<T> operator%(const poli<T> &p,const poli<T> &q){
  return div(p,q,'r');
}

//////////////////////////////FUNCIONES VARIAS////////////////////////

//Impresion de polinomios, se cuida no imprimir ceros (a menos que sea el unico elemento del polinomio), no imprimir + antes de un menos e imprimirlos de mayor a menos.
template <class T>
ostream& operator<<(ostream &os,const poli<T> &p){
  for(int i=p.grad();i>=0;--i){
    if(i==p.grad() && i!=0)
      os<<p[i]<<"x^"<<i;
    else if(p[i]<0){
      if(i!=0)
	os<<p[i]<<"x^"<<i;
      else
	os<<p[i];
    }
    else if(p[i]>0){
      if(i!=0)
	os<<'+'<<p[i]<<"x^"<<i;
      else
	os<<'+'<<p[i];	
    }
    else if(p[i]==0 && p.grad()==0)
      os<<0;
  }
  return os;
}

template <class T>
poli<T> poli<T>::legendre(uint n){
  poli<T> temp(2,1);temp[1]=0;temp[0]*=-1;
  poli<T> p(0,1);
  for(uint i=1;i<=n;++i) p=p*temp;
  for(uint i=1;i<=n;++i) p=(1/(2*double(i)))*p.deriv();
  return p;
}
#endif
