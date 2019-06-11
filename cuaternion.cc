#include "cuaternion.h"

cuaternion::cuaternion():r(0),i(0),j(0),k(0)
{}

cuaternion::cuaternion(double x, double ii, double jj, double kk):r(x),i(ii),j(jj),k(kk)
{}

cuaternion::~cuaternion(){}

cuaternion::cuaternion(const cuaternion &q):r(q.r),i(q.i),j(q.j),k(q.k)
{
}

cuaternion &cuaternion::operator=(const cuaternion &d){
	if(this!=&d){
		r=d.r;
		i=d.i;
		j=d.j;
		k=d.k;	
	}
	return *this;	
}

cuaternion cuaternion::inv(){
	cuaternion d=this->conj();
	for(int m=0;m<4;++m){
		d[m]/=this->mod()*(this->mod());
	}
	return d;
}

cuaternion cuaternion::conj(){
	return cuaternion(this->r,-this->i,-this->j,-this->k);
}

double cuaternion::mod(){
	return sqrt(this->r*(this->r)+this->i*(this->i)+this->j*(this->j)+this->k*(this->k));
}

cuaternion operator+(const cuaternion&q1, const cuaternion&q2){
	return cuaternion(q1[0]+q2[0],q1[1]+q2[1],q1[2]+q2[2],q1[3]+q2[3]);
}

cuaternion operator-(const cuaternion&q1, const cuaternion&q2){
	return cuaternion(q1[0]-q2[0],q1[1]-q2[1],q1[2]-q2[2],q1[3]-q2[3]);
}

cuaternion operator*(const cuaternion&q1, const cuaternion&q2){
	return cuaternion(q1[0]*q2[0]-q1[1]*q2[1]-q1[2]*q2[2]-q1[3]*q2[3],q1[0]*q2[1]+q1[1]*q2[0]+q1[2]*q2[3]-q1[3]*q2[2],q1[0]*q2[2]-q1[1]*q2[3]+q1[2]*q2[0]+q1[3]*q2[1],q1[0]*q2[3]+q1[1]*q2[2]-q1[2]*q2[1]+q1[3]*q2[0]);
}

cuaternion operator/(const cuaternion&q1, const cuaternion&q2){
	cuaternion h=q2;
	cuaternion d=q1*(h.inv());
	return d;
}

double cuaternion::operator[](int m) const 
{
	if(m==0)
		return r;
	else if(m==1)
		return i;
	else if(m==2)
		return j;
	else if(m==3)
		return k;
}

double &cuaternion::operator[](int m) {
	if(m==0)
		return r;
	else if(m==1)
		return i;
	else if(m==2)
		return j;
	else if(m==3)
		return k;
}