#ifndef RACIONALES_
#define RACIONALES_
#include <iostream>

using namespace std;

//Clase de racionales, lo no explicado se explica en el cc.
class racional{
 private:
  long num,den; //Partes privadas como numerador y denominador
 public:
  racional(); 
  racional(long,long); 
  racional(double); 
  racional(const racional&); 
  racional & operator=(const racional&); 
  long operator[](int) const; 
  long & operator[](int); 
  double dec(); 
  racional pow(int); 
};

long mcd(long,long);
racional operator+(const racional&,const racional&);
racional operator-(const racional&,const racional&);
racional operator*(const racional&,const racional&);
racional operator/(const racional&,const racional&);
ostream & operator<<(ostream&,const racional &);

#endif
