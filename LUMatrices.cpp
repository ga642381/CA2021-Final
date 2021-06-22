#include "classes.h"

LowerMatrix::LowerMatrix(int n) {
    this->n=n;
    this->dim=n*n;
    this->diagonal.assign(dim,1);
    this->tridiagonal.assign(dim-1,0);
    this->identity.assign(dim-n,0);
}

LowerMatrix::~LowerMatrix() {
    vector<double>().swap(diagonal);
    vector<double>().swap(tridiagonal);
    vector<double>().swap(identity);
}

int LowerMatrix::Size() {
    return dim;
}

void LowerMatrix::Set(int i,int j,double value) {
    if(i==j) {
        diagonal[i]=value;
    } else if(j==(i-1)) {
        tridiagonal[j]=value;
    } else if(j==(i-n)) {
        identity[j]=value;
    }
}

double LowerMatrix::Get(int i,int j) {
    if(i==j) {
        return diagonal[i];
    } else if(j==(i-1)) {
        return tridiagonal[j];
    } else if(j==(i-n)) {
        return identity[j];
    } else {
        return 0.0;
    }
}

std::vector<double> LowerMatrix::operator*(const std::vector<double>& z) {
    dim=z.size();
    std::vector<double> tmp(z);
    for(int i=0;i<dim;i++) {
        if(i>=n) {
            tmp[i]-=tmp[i-n]*Get(i,i-n);
        }
        if(i%n!=0) {
            tmp[i]-=tmp[i-1]*Get(i,i-1);
        }
        tmp[i]/=Get(i,i);
    }
    return tmp;
}

UpperMatrix::UpperMatrix(int n) : LowerMatrix(n) {
   this->n=n;
   this->dim=n*n;
   this->diagonal.assign(dim,0);
   this->tridiagonal.assign(dim-1,0);
   this->identity.assign(dim-n,0);
}

UpperMatrix::~UpperMatrix() {
   vector<double>().swap(diagonal);
   vector<double>().swap(tridiagonal);
   vector<double>().swap(identity);
}

double UpperMatrix::Get(int i,int j) {
   return LowerMatrix::Get(j,i);
}

void UpperMatrix::Set(int i,int j,double value) {
   return LowerMatrix::Set(j,i,value);
}

std::vector<double> UpperMatrix::operator*(const std::vector<double>& z) {
    std::vector<double> tmp(z);
    for(int i=dim-1;i>=0;i--) {
        tmp[i]-=tmp[i]*Get(i,i);
        if(i<=(dim-n)) {
            tmp[i]-=tmp[i+n]*Get(i,i+n);
        }
        if(i%n!=0) {
            tmp[i]-=tmp[i+1]*Get(i,i+1);
        }
        tmp[i]/=Get(i,i);
    }
    return tmp;
}