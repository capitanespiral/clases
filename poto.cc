#include "matriz.h"

int main(){
  vector<double> v={1,7,1,3,8,1,2,5,4};
  vector<double> w={1,2,3,7};
  matriz<double> m(3,v);
  matriz<double> h(w,1);
  cout<<m<<endl;
  cout<<m.triang(true)<<endl;
cout<<m.det(true)<<endl;
  return 0;

}
