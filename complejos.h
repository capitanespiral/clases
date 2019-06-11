#ifndef COMPLEJOS_
#define COMPLEJOS_

#include <iostream>
using namespace std;

template <class T>

class comple{
 private:
  T r;
  T i;
 public:
  comple():r(0),i(0){}
  comple(T a,T b){
    r=a;
    i=b;
  }
  ~comple(){}
  comple(const comple &c){
    r=c.r;
    i=c.i;
  }
  comple& operator=(const comple &c){
    if(this!=&c){
      r=c.r;
      i=c.i;
    }
  }
  T operator[](int h) const{
    if(h==0) return r;
    else return i;
  }
  T& operator[](int h){
    if(h==0) return r;
    else return i;
    }
  comple & operator++(){
    r=r+1;
    return *this;
  }
  comple & operator++(int){
    i=i+1;
    return *this;
  }
};
#endif
