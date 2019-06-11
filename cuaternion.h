#ifndef CUATERNION_H
#define CUATERNION_H

#include <iostream>
#include <cmath>

using namespace std;

class cuaternion{
	protected: 
		double r,i,j,k;
	
	public:
		cuaternion();
		cuaternion(double,double=0,double=0,double=0);
		~cuaternion();
		cuaternion(const cuaternion&);//shallow copy, pasar por valor
		cuaternion & operator=(const cuaternion&);//deep copy
		double & operator[](int);
		double operator[](int) const;
		double mod();
		cuaternion conj();
		cuaternion inv();
		friend ostream & operator<<(ostream &os, const cuaternion &q)  //<< insercion, >> extraccion
		{
			//os<<"[["<<q.xR<<","<<qouble x, double i, double j, double k):xR(x),xI(i),xJ(j),xK(k).xI<<","<<q.xJ<<","<<q.xK<<"]]";
			return q.print(os);		
		}
		virtual ostream & print(ostream &os) const
			{
				os<<"[["<<this->r<<","<<this->i<<","<<this->j<<","<<this->k<<"]]";
				return os;
			}
};

cuaternion operator*(const cuaternion&,const cuaternion&);
cuaternion operator+(const cuaternion&,const cuaternion&);
cuaternion operator-(const cuaternion&,const cuaternion&);
cuaternion operator/(const cuaternion&,const cuaternion&);
#endif
