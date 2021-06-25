#include "classes.h"

PoissonMatrix::PoissonMatrix(int n) {
   this->n=n;
   this->dim=n*n;
}


PoissonMatrix::~PoissonMatrix() {
    vector<vector<int> >().swap(HashMatrix);
}

int PoissonMatrix::Size() {
   return dim;
}

vector<double> PoissonMatrix::operator*(const vector<double>& x) {
    int dim=x.size(),n=sqrt(dim);
    vector<double> tmp(dim,0);
    for(int i=0;i<dim;i++) {
        tmp[i]+=x[i]*4.0*pow(n+1,2);
        if(i<(dim-n)) {
            tmp[i]+=x[i+n]*-1.0*pow(n+1,2);
            tmp[i+n]+=x[i]*-1.0*pow(n+1,2);
        }
        if(i%n!=0) {
            tmp[i]+=x[i-1]*-1.0*pow(n+1,2);
            tmp[i-1]+=x[i]*-1.0*pow(n+1,2);
        }
    }
    return tmp;
}
