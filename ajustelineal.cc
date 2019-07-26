#include "ajustelineal.h"

ajus_lin::ajus_lin(const matriz<double> &m){
  datos=m;
  cant=m.fila();
}

matriz<double> ajus_lin::min_cuad(const vector<double (*)(double)> &v){
  matriz<double> temp(cant,v.size());
  for(int i=0;i<temp.fila();++i){
    for(int j=0;j<temp.colu();++j){
      temp(i,j)=v[j](datos(i,0));
    }
  }
  matriz<double> temptrans=temp.tras();
  temp=temptrans*temp;
  temptrans=temptrans*datos.colu(1);
  temp=sist_ec_lin(temp,temptrans);
  return temp;
}

matriz<double> evalua_coef(const vector<double (*)(double)> &v,const matriz<double> &m,double start,double fin,int puntos){
  matriz<double> res(puntos,2);double acum=0;double step=(fin-start)/puntos;
  for(int i=0;i<res.fila();++i){
    res(i,0)=start;
    for(int j=0;j<v.size();++j){
      acum=acum+m(j,0)*v[j](res(i,0));
    }
    res(i,1)=acum;acum=0;start+=step;
  }
  return res;
}
