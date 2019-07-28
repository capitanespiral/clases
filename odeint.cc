#include "odeint.hh"

//Defino uint (mas comodo)
typedef unsigned int uint;

//Ctes de PEFRL
const double XI=0.1786178958448091;
const double LAMBDA=-0.2123418310626054;
const double JI=-0.06626458266981849;

//Constructor de odeint, recibe tiempo inicial, paso, posiciones iniciales, velocidades iniciales, y variable "comodin" inicial (como masas o cargas).
//Las matrices se esperan como -> cada fila un cuerpo, cada columna nueva coordenada.
odeint::odeint(double a,double b,const matriz<double> &xx,const matriz<double> &vv,const matriz<double> &qq){
  t=a;h=b;x=xx,v=vv;q=qq;
}

//Para una derivada
odeint::odeint(double a,double b,const matriz<double> &xx,const matriz<double> &qq){
  t=a;h=b;x=xx;q=qq;
}

//Funcion que cuenta pasos totales desde ti hasta tf avanzando h en cada paso.
int cuentapasos(double ti,double tf,double h){
  int res=1;
  while(ti<=tf){
    ti+=h;
    res+=1;
  }
  return res;
}

//Metodo de euler-cromer. Recibe funcion (aceleracion - segunda derivada) tiempo final y un indicador para guardar todos los datos o solo el ultimo. Esta misma recibe tiempo, posiciones, velocidades y la variable "comodin"
matriz<matriz<double>> odeint::euler_cromer(matriz<double> (*g)(double,const matriz<double> &,const matriz<double> &,const matriz<double> &),double tf, char c){
  //Si no se define tf, se corre un paso.
  if(tf==0) tf=t;
  //Creo variables importantes para el algoritmo
  matriz<double> m,x_c,v_c;
  matriz<double> a=(*g)(t,x,v,q);
  int colum=1;
  double t_c;
  vector<double> t_i;
  //Si solo quiero un resultado para un punto o si quiero todos hasta cierto punto
  if(c!='u')
    //Cantidad de columnas para la matriz de rptas
    colum=cuentapasos(t,tf,h);
  //Creo la matriz de matrices que almacenara las rptas
  matriz<matriz<double>> res(4,colum,m);//Filas t,x,v,a. Columnas cada paso (x,v,a son de todos los cuerpos con todas las coordenadas requeridas, por eso matriz de matrices)
  if(c!='u'){//Si quiero todos los resultados guardo el inicial primero
    t_i={t};matriz<double> T(1,t_i);
    res(0,0)=T;res(1,0)=x;res(2,0)=v;res(3,0)=a;
  }
  //Comienza bucle de ODEINT
  //Reciclo colum, y copio los valores iniciales
  t_c=t;x_c=x;v_c=v;
  colum=1;
  //mientras el tiempo es menor que el punto pedido
  while(t_c<=tf){
    //Velocidad
    v_c=v_c+h*a;
    //Posicion
    x_c=x_c+h*v_c;
    //Avanzo un paso y recalculo la aceleracion
    t_c+=h;a=(*g)(t_c,x_c,v_c,q);
    //Si quiero todos los pasos guardo este
    if(c!='u'){
      t_i={t_c};matriz<double> T(1,t_i);
      res(0,colum)=T;res(1,colum)=x_c;res(2,colum)=v_c;res(3,colum)=a;
      colum+=1;
    }
    //De querer solo un punto guardo ese (el final)
    if(c=='u' && t_c>tf){
      t_i={t_c};matriz<double> T(1,t_i);
      res(0,0)=T;res(1,0)=x_c;res(2,0)=v_c;res(3,0)=a;
    } 
  }
  //Devuelvo resultado
  return res;
}

//Metodo rk4-Recibe lo mismo que euler-cromer, y casi todo es igual, mencionare lo distinto.
matriz<matriz<double>> odeint::rk4(matriz<double> (*g)(double,const matriz<double> &,const matriz<double> &,const matriz<double> &),double tf, char c){
  if(tf==0) tf=t;
  vector<double> t_i;
  matriz<double> k1,k_1,k2,k_2,k3,k_3,k4,k_4,x_c,v_c,m;//ki sobre la primera derivada y k_i sobre la segunda
  int colum=1;
  double t_c;matriz<double> a=(*g)(t,x,v,q);
  if(c!='u')
    colum=cuentapasos(t,tf,h);
  matriz<matriz<double>> res(4,colum,m);
  if(c!='u'){
    t_i={t};matriz<double> T(1,t_i);
    res(0,0)=T;res(1,0)=x;res(2,0)=v;res(3,0)=a;
  }
  t_c=t;x_c=x;v_c=v;
  colum=1;
  while(t_c<=tf){
    //Primero calculamos los ki y los k_i
    k1=v_c;k_1=a;
    k2=v_c+0.5*h*k_1;
    k_2=(*g)(t_c+0.5*h,x_c+0.5*h*k1,v_c+0.5*h*k_1,q);
    k3=v_c+0.5*h*k_2;
    k_3=(*g)(t_c+0.5*h,x_c+0.5*h*k2,v_c+0.5*h*k_2,q);
    k4=v_c+h*k_3;k_4=(*g)(t_c+h,x_c+h*k3,v_c+h*k_3,q);
    //Calculo las nuevas posiciones y velocidades
    x_c=x_c+(h/6.0)*(k1+2*k2+2*k3+k4);
    v_c=v_c+(h/6.0)*(k_1+2*k_2+2*k_3+k_4);
    //Avanzo un paso y calculo la aceleracion
    t_c+=h;a=(*g)(t_c,x_c,v_c,q);
    if(c!='u'){
      t_i={t_c};matriz<double> T(1,t_i);
      res(0,colum)=T;res(1,colum)=x_c;res(2,colum)=v_c;res(3,colum)=a;
      colum+=1;
    }
    if(c=='u' && t_c>tf){
      t_i={t_c};matriz<double> T(1,t_i);
      res(0,0)=T;res(1,0)=x_c;res(2,0)=v_c;res(3,0)=a;
    } 
  }
  return res;
}

//rk4 para una derivada
matriz<matriz<double>> odeint::rk4(matriz <double> (*g)(double,const matriz<double> &,const matriz<double> &),double tf,char c){
  if(tf==0) tf=t;
  vector<double> t_i;
  matriz<double> k1,k2,k3,k4,x_c,m;//ki sobre la primera derivada y k_i sobre la segunda
  int colum=1;
  double t_c;matriz<double> V=(*g)(t,x,q);
  if(c!='u')
    colum=cuentapasos(t,tf,h);
  matriz<matriz<double>> res(3,colum,m);
  if(c!='u'){
    t_i={t};matriz<double> T(1,t_i);
    res(0,0)=T;res(1,0)=x;res(2,0)=V;
  }
  t_c=t;x_c=x;
  colum=1;
  while(t_c<=tf){
    //Primero calculamos los ki y los k_i
    k1=V;k2=(*g)(t+0.5*h,x_c+0.5*h*k1,q);
    k3=(*g)(t+0.5*h,x_c+0.5*h*k2,q);k4=(*g)(t+h,x_c+h*k3,q);
    //Calculo las nuevas posiciones
    x_c=x_c+(h/6.0)*(k1+2*k2+2*k3+k4);
    //Avanzo un paso y calculo la velocidad
    t_c+=h;V=(*g)(t_c,x_c,q);
    if(c!='u'){
      t_i={t_c};matriz<double> T(1,t_i);
      res(0,colum)=T;res(1,colum)=x_c;res(2,colum)=V;
      colum+=1;
    }
    if(c=='u' && t_c>tf){
      t_i={t_c};matriz<double> T(1,t_i);
      res(0,0)=T;res(1,0)=x_c;res(2,0)=V;
    } 
  }
  return res;
}

//Metodo velocity-velver (tambien recibe lo mismo), todo excepto las tres primeras lineas del while es identico a euler-cromer.
matriz<matriz<double>> odeint::vel_velver(matriz <double> (*g)(double,const matriz<double> &,const matriz<double> &,const matriz<double> &),double tf,char c){
  if(tf==0) tf=t;
  vector<double> t_i;
  int colum=1;
  double t_c;matriz<double> x_c,v_c,x_t,m;
  matriz<double> a=(*g)(t,x,v,q);
  if(c!='u')
    colum=cuentapasos(t,tf,h);
  matriz<matriz<double>> res(4,colum,m);
  if(c!='u'){
    t_i={t};matriz<double> T(1,t_i);
    res(0,0)=T;res(1,0)=x;res(2,0)=v;res(3,0)=a;
  }
  t_c=t;x_c=x;v_c=v;
  colum=1; 
  while(t_c<=tf){
    x_t=x_c+h*v_c+0.5*h*h*a;
    v_c=v_c+0.5*h*((*g)(t_c,x_c,v_c,q)+(*g)(t_c+h,x_t,v_c,q));
    x_c=x_t;t_c+=h;a=(*g)(t_c,x_c,v_c,q);
    if(c!='u'){
      t_i={t_c};matriz<double> T(1,t_i);
      res(0,colum)=T;res(1,colum)=x_c;res(2,colum)=v_c;res(3,colum)=a;
      colum+=1;
    }
    if(c=='u' && t_c>tf){
      t_i={t_c};matriz<double> T(1,t_i);
      res(0,0)=T;res(1,0)=x_c;res(2,0)=v_c;res(3,0)=a;
    } 
  }
  return res;
}

//Metodo pefrl. Tambien recibe lo mismo y es muy parecido, solo mencionare las diferencias
matriz<matriz<double>> odeint::pefrl(matriz<double> (*g)(double,const matriz<double> &,const matriz<double> &,const matriz<double> &),double tf,char c){
  if(tf==0) tf=t;
  vector<double> t_i;
  int colum=1;
  double t_c;matriz<double> x_c,v_c,m;
  matriz<double> a=(*g)(t,x,v,q);
  if(c!='u')
    colum=cuentapasos(t,tf,h);
  matriz<matriz<double>> res(4,colum,m);
  if(c!='u'){
    t_i={t};matriz<double> T(1,t_i);
    res(0,0)=T;res(1,0)=x;res(2,0)=v;res(3,0)=a;
  }
  t_c=t;x_c=x;v_c=v;
  colum=1; 
  while(t_c<=tf){
    //El metodo utiliza las tres ctes definidas a principio del archivo
    x_c=x_c+::XI*h*v_c;
    v_c=v_c+(1-2*::LAMBDA)*(h/2.0)*(*g)(t_c,x_c,v_c,q);
    x_c=x_c+::JI*h*v_c;
    v_c=v_c+::LAMBDA*h*(*g)(t_c,x_c,v_c,q);
    x_c=x_c+(1-2*(::JI+::XI))*h*v_c;
    v_c=v_c+::LAMBDA*h*(*g)(t_c,x_c,v_c,q);
    x_c=x_c+::JI*h*v_c;
    v_c=v_c+(1-2*::LAMBDA)*(h/2.0)*(*g)(t_c,x_c,v_c,q);
    x_c=x_c+::XI*h*v_c;
    t_c+=h;a=(*g)(t_c,x_c,v_c,q);
    if(c!='u'){
      t_i={t_c};matriz<double> T(1,t_i);
      res(0,colum)=T;res(1,colum)=x_c;res(2,colum)=v_c;res(3,colum)=a;
      colum+=1;
    }
    if(c=='u' && t_c>tf){
      t_i={t_c};matriz<double> T(1,t_i);
      res(0,0)=T;res(1,0)=x_c;res(2,0)=v_c;res(3,0)=a;
    } 
  }
  return res;
}

matriz<matriz<double>> odeint::rk4_adap(matriz<double> (*g)(double,const matriz<double> &,const matriz<double> &,const matriz<double> &),double tf,char c){
  if(tf==0) tf=t;
  //err_u de usuario.
  double paso=h;double err_u=1.0e-3;double s1=0.9;double s2=2;
  vector<double> t_i;
  matriz<double> k1,k_1,k2,k_2,k3,k_3,k4,k_4,x_c1,v_c1,x_c2,v_c2,x_aux,v_aux,m;//ki sobre la primera derivada y k_i sobre la segunda
  //err_l error local, paso nuevo
  double t_c,t_aux,err_l,nue_paso;matriz<double> a=(*g)(t,x,v,q);
  matriz<matriz<double>> res(4,1,m);matriz<matriz<double>> temp=res;
  if(c!='u'){
    t_i={t};matriz<double> T(1,t_i);
    res(0,0)=T;res(1,0)=x;res(2,0)=v;res(3,0)=a;
  }
  t_c=t;t_aux=t;x_c1=x;v_c1=v;x_c2=x;v_c2=v;
  while(t_c<=tf){
    //Primero calculamos los ki y los k_i PASO GRANDE
    x_aux=x_c1;v_aux=v_c1;
    k1=v_c1;k_1=a;
    k2=v_c1+0.5*paso*k_1;
    k_2=(*g)(t_c+0.5*paso,x_c1+0.5*paso*k1,v_c1+0.5*paso*k_1,q);
    k3=v_c1+0.5*paso*k_2;
    k_3=(*g)(t_c+0.5*paso,x_c1+0.5*paso*k2,v_c1+0.5*paso*k_2,q);
    k4=v_c1+paso*k_3;k_4=(*g)(t_c+paso,x_c1+paso*k3,v_c1+paso*k_3,q);
    //Calculo las nuevas posiciones y velocidades
    x_c1=x_c1+(paso/6.0)*(k1+2*k2+2*k3+k4);
    v_c1=v_c1+(paso/6.0)*(k_1+2*k_2+2*k_3+k_4);
    //Comienza PASO CHICO
    for(int i=0;i<2;++i){
      k1=v_c2;k_1=(*g)(t_aux,x_c2,v_c2,q);
      k2=v_c2+0.5*0.5*paso*k_1;
      k_2=(*g)(t_aux+0.5*0.5*paso,x_c2+0.5*0.5*paso*k1,v_c2+0.5*0.5*paso*k_1,q);
      k3=v_c2+0.5*0.5*paso*k_2;
      k_3=(*g)(t_aux+0.5*0.5*paso,x_c2+0.5*0.5*paso*k2,v_c2+0.5*0.5*paso*k_2,q);
      k4=v_c2+0.5*paso*k_3;k_4=(*g)(t_aux+0.5*paso,x_c2+0.5*paso*k3,v_c2+0.5*paso*k_3,q);
      x_c2=x_c2+(0.5*paso/6)*(k1+2*k2+2*k3+k4);
      v_c2=v_c2+(0.5*paso/6)*(k_1+2*k_2+2*k_3+k_4);
      t_aux+=0.5*paso;
    }
    t_aux=t_c;t_c+=paso;
    //Analizo el error
    err_l=max(mod(x_c1-x_c2));
    //Si es mayor calculo un nuevo paso y repito el ciclo para el mismo momento.
    nue_paso=paso*pow(err_u/err_l,0.2);
    if(nue_paso*s1>s2*paso) nue_paso=s2*paso;
    else if(nue_paso*s1<paso/s2) nue_paso=paso/s2;
    else nue_paso=s1*nue_paso;
    paso=nue_paso;
    if(err_l>err_u){
      x_c1=x_aux;x_c2=x_aux;v_c1=v_aux;v_c2=v_aux;t_c=t_aux;
      continue;
    }
    a=(*g)(t_c,x_c1,v_c1,q);
    if(c!='u'){
      t_i={t_c};matriz<double> T(1,t_i);
      temp(0,0)=T;temp(1,0)=x_c1;temp(2,0)=v_c1;temp(3,0)=a;
      res=cat(res,temp);
    }
    if(c=='u' && t_c>tf){
      t_i={t_c};matriz<double> T(1,t_i);
      res(0,0)=T;res(1,0)=x_c1;res(2,0)=v_c1;res(3,0)=a;
    } 
  }
  return res;
}


 

//Funcion que "desarma" las variables de la matriz de matrices resultante en solo una matriz (esta actua sobre las variables de una dimension y con solo una, por defecto, siempre el tiempo al menos)
//Recibe la matriz de matrices resultante y la fila que uno desea.
matriz<double> desarma(const matriz<matriz<double>> &m,int a){
  matriz<double> n(m(a,0).fila(),m.colu());
  for(int j=0;j<m(a,0).fila();++j){
    for(int i=0;i<m.colu();++i){
      n(j,i)=(m(a,i))(j,0);
    }
  }
  return n;
}

//Hace lo mismo pero de forma general, cuando cada elemento de la matriz de la matrices puede ser una matriz de dimensiones arbitrarias. Recibe lo mismo que la anterior.
matriz<double> desarma1(const matriz<matriz<double>> &m,int a){
  matriz<double> n(m(a,0).fila(),m.colu()*m(a,0).colu());
  for(int i=0;i<m(a,0).fila();++i){
    for(int j=0;j<m.colu();++j){
      for(int k=0;k<m(a,0).colu();++k){
	n(i,j*m(a,0).colu()+k)=(m(a,j))(i,k);
      }
    }
  }
  return n;
}

//Funcion que crea archivos para cada cuerpo involucrado, con la primera columna el tiempo y la segunda la fila que pidamos de la matriz de matrices. Recibe s (nombre del archivo sin extension, .dat por defecto) la matriz de matrices y la fila que deseemos.
void file(string s,const matriz<matriz<double>> &m,int a,char b){
  //variables importantes para el algoritmo
  string part;
  matriz<double> n,t;
  //almaceno el tiempo en "t"
  if(b=='t')
    t=desarma(m,0);
  //Si tiene mas de una dimension uso desarma1 y "marco la marca"
  if(m(a,0).colu()>1){
    n=desarma1(m,a);
  }
  //De no ser el caso, uso desarma normal.
  else
    n=desarma(m,a);
  //Primer ciclo que crea un archivo por cuerpo interesado.
  for(int k=0;k<m(a,0).fila();++k){
    part=to_string(k+1);
    ofstream f(s+'_'+part+".dat");
    //Segundo ciclo que guarda el tiempo
    for(int i=0;i<m.colu();++i){
      if(b=='t')
	f<<t(0,i)<<' ';
      //Si la marca fue marcada, guardo cada coordenada separada por un espacio,y luego el fin de linea.
      for(int j=0;j<m(a,0).colu();++j)
	f<<n(k,i*m(a,0).colu()+j)<<' ';
      f<<endl;
    }
    //Cierro el archivo, y se repite el ciclo de haber mas de un cuerpo.
    f.close();
  }
}
