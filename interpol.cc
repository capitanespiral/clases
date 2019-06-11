#include "interpol.h"

using namespace std;

//Constructor de i_pol, recibe la matriz con los datos, se asume primera fila como x, segunda como y, y cantidad de filas como cantidad de datos
i_pol::i_pol(const matriz<double> &m):datos(m),cant(m.fila()){}

//Interpolacion de lagrange, retorna el polinomio de la interpolacion
poli<double> i_pol::lagrange(){
  poli<double> res,p,temp;//polinomios para el calculo
  poli<double> lineal(1,1);//polinomio lineal
  poli<double> uno(0,1);//polinomio cualquiera (parte como uno)
  //primer bucle para cada uno de los polinomios (uno por dato)
  for(int i=0;i<cant;++i){
    temp=uno;
    //segundo bucle que construye cada uno de estos
    for(int j=0;j<cant;++j){
      if(i==j) continue;
      p=lineal;
      p[0]=-1*datos(j,0);p=p*(1.0/(datos(i,0)-datos(j,0)));
      temp=temp*p;
    }
    //Lo voy acumulando en res, el polinomio final.
    res=res+temp*datos(i,1);
  }
  return res;
}

//Intepolacion con splines periodicos, se entrega una matriz con una columna y n filas, conteniendo en cada una uno de los polinomios del spline.
//Por ahora funciona con 4 o mas puntos.Siempre de menor a mayor x. Solo splines naturales.
//n+1 pts -> n polinomios de splines -> n+1 M's -> matriz -m0 y -mn queda de n-1xn-1.
matriz<poli<double>> i_pol::spline_p(){
  //coef guarda la matriz tridiagonal, coef1 los "resultados" al otro lado de la ecuacion matricial.Los h y los k, variables que aparecen en ambas.
  matriz<double> coef(cant-2,cant-2);
  vector<double> coef1;double hi,hi1,ki,ki1,h1,k1,hn,kn;
  //Guardo la primera y la ultima fila que son casos singulares
  //Primera fila
  hi=datos(1,0)-datos(0,0);hi1=datos(2,0)-datos(1,0);
  ki=datos(1,1)-datos(0,1);ki1=datos(2,1)-datos(1,1);
  hn=datos(cant-1,0)-datos(cant-2,0);
  kn=datos(cant-1,1)-datos(cant-2,1);
  coef(0,0)=((hi+hi1)/3)-hi*hi/(12*(hi+hn));
  coef(0,1)=hi1/6;coef(0,cant-3)=-1*hi*hn/(12*(hi+hn));
  coef1.push_back(ki1/hi1-ki/hi+kn*hi/(2*hn*(hi+hn))-ki/(2*(hi+hn)));
 
  //Ciclo que guarda los casicoeficientes de al medio
  for(int i=1;i<cant-3;++i){
    hi=datos(i+1,0)-datos(i,0);hi1=datos(i+2,0)-datos(i+1,0);
    ki=datos(i+1,1)-datos(i,1);ki1=datos(i+2,1)-datos(i+1,1);
    coef(i,i-1)=hi/(hi+hi1);coef(i,i)=2;coef(i,i+1)=hi1/(hi+hi1);
    coef1.push_back(6*(ki1/hi1-ki/hi)/(hi+hi1));
  }
   //ultima fila
  hi=datos(cant-2,0)-datos(cant-3,0);hi1=datos(cant-1,0)-datos(cant-2,0);
  ki=datos(cant-2,1)-datos(cant-3,1);ki1=datos(cant-1,1)-datos(cant-2,1);
  h1=datos(1,0)-datos(0,0);
  k1=datos(1,1)-datos(0,1);
  coef(cant-3,0)=-1*h1*hi1/(12*(h1+hi1));
  coef(cant-3,cant-4)=hi/6;
  coef(cant-3,cant-3)=(hi+hi1)/3-hi1*hi1/(12*(h1+hi1));
  coef1.push_back(ki1/hi1-ki/hi+ki1/(2*(h1+hi1))-hi1*k1/(2*h1*(h1+hi1)));
  
  //Ahora sacamos las M (seg derivada evaluada), falta M0 y Mn, asi que los "cateo"
  matriz<double> m;
  coef=ANW(coef,coef1);
  m(0,0)=(3/(hi1+h1))*(k1/h1-ki1/hi1-h1*coef(0,0)/6-hi1*coef(0,cant-3)/6);
  coef=cat(m,coef);coef=cat(coef,m);
  //Armamos polinomios y los guardamos en res.
  poli<double> p,temp,temp1;matriz<poli<double>> res(cant-1,1,p);
  poli<double> lineal(1,1);
  for(int i=1;i<cant;++i){
    temp=-1*lineal;temp1=lineal;
    hi=datos(i,0)-datos(i-1,0);
    ki=(datos(i-1,1)-(coef(0,i-1)*hi*hi)/6)/hi;
    ki1=(datos(i,1)-(coef(0,i)*hi*hi)/6)/hi;
    temp[0]=datos(i,0);temp1[0]=-1*datos(i-1,0);
    p=(coef(0,i-1)/(6*hi))*temp*temp*temp+(coef(0,i)/(6*hi))*temp1*temp1*temp1+ki*temp+ki1*temp1;
    res(i-1,0)=p;
  }
  return res;
}


//splines naturales
matriz<poli<double>> i_pol::spline(){
  //coef guarda la matriz tridiagonal, coef1 los "resultados" al otro lado de la ecuacion matricial.Los h y los k, variables que aparecen en ambas.
  matriz<double> coef(cant-2,cant-2);
  vector<double> coef1;double hi,hi1,ki,ki1;
  //Guardo la primera y la ultima fila que son casos singulares
  //Primera fila
  hi=datos(1,0)-datos(0,0);hi1=datos(2,0)-datos(1,0);
  ki=datos(1,1)-datos(0,1);ki1=datos(2,1)-datos(1,1);
  coef(0,0)=2;coef(0,1)=hi1/(hi1+hi);coef1.push_back(6*(ki1/hi1-ki/hi)/(hi+hi1));
  //Ciclo que guarda los casicoeficientes de al medio
  for(int i=1;i<cant-3;++i){
    hi=datos(i+1,0)-datos(i,0);hi1=datos(i+2,0)-datos(i+1,0);
    ki=datos(i+1,1)-datos(i,1);ki1=datos(i+2,1)-datos(i+1,1);
    coef(i,i-1)=hi/(hi+hi1);coef(i,i)=2;coef(i,i+1)=hi1/(hi+hi1);
    coef1.push_back(6*(ki1/hi1-ki/hi)/(hi+hi1));
  }
  //ultima fila
  hi=datos(cant-2,0)-datos(cant-3,0);hi1=datos(cant-1,0)-datos(cant-2,0);
  ki=datos(cant-2,1)-datos(cant-3,1);ki1=datos(cant-1,1)-datos(cant-2,1);
  coef(cant-3,cant-4)=hi/(hi+hi1);coef(cant-3,cant-3)=2;
  coef1.push_back(6*(ki1/hi1-ki/hi)/(hi+hi1));

  //Ahora sacamos las M (seg derivada evaluada), falta M0 y Mn que son cero (para caso natural), asi que los "cateo"
  matriz<double> m;
  coef=thomas(coef,coef1);
  coef=cat(m,coef);coef=cat(coef,m);
  //Armamos polinomios y los guardamos en res.
  poli<double> p,temp,temp1;matriz<poli<double>> res(cant-1,1,p);
  poli<double> lineal(1,1);
  for(int i=1;i<cant;++i){
    temp=-1*lineal;temp1=lineal;
    hi=datos(i,0)-datos(i-1,0);
    ki=(datos(i-1,1)-(coef(0,i-1)*hi*hi)/6)/hi;
    ki1=(datos(i,1)-(coef(0,i)*hi*hi)/6)/hi;
    temp[0]=datos(i,0);temp1[0]=-1*datos(i-1,0);
    p=(coef(0,i-1)/(6*hi))*temp*temp*temp+(coef(0,i)/(6*hi))*temp1*temp1*temp1+ki*temp+ki1*temp1;
    res(i-1,0)=p;
  }
  return res;
}

//Funcion que crea matriz con datos evaluados segun interpolacion de splines.
//Se asume primer y ultimo punto los correspondientes a i_pol, a como totalidad de puntos (no pueden ser menores a la cantidad de puntos iniciales claramente.
matriz<double> puntos(const matriz<poli<double>> &m,i_pol n,int a){
  matriz<double> res(a,2);double x,y;int j=0;
  for(int i=0;i<a;++i){
    x=(n.datos(0,0)*((a-1-double(i))/(a-1))+n.datos(n.cant-1,0)*(double(i)/(a-1)));//particion de a pts tomando en cuenta el primero y el ultimo.
    //Aqui cambio el polinomio con que se evalua si es que cambie de intervalo
    if(x>n.datos(j+1,0)) ++j;
    y=m(j,0)(x);
    res(i,0)=x;res(i,1)=y;
  }
  return res;
}

//Funcion que hace lo mismo que la anterior pero para el caso de lagrange
matriz<double> puntos(poli<double> p, i_pol n,int a){
  matriz<double> res(a,2);double x,y;
  for(int i=0;i<a;++i){
    x=(n.datos(0,0)*((a-1-double(i))/(a-1))+n.datos(n.cant-1,0)*(double(i)/(a-1)));
    y=p(x);
    res(i,0)=x;res(i,1)=y;
  }
  return res;
}

//Funcion que manda la info a un archivo. Recibe un string para usarlo de nombre con extension .dat
void file(string s,const matriz<double> &m){
  ofstream archivo(s+".dat");
  for(int i=0;i<m.fila();++i){
    archivo<<m(i,0)<<' '<<m(i,1)<<endl;
  }
  archivo.close();
}
