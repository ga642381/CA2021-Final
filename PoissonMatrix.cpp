#include "classes.h"

PoissonMatrix::PoissonMatrix(int n) {
   this->n=n;
   this->dim=n*n;
   this->diagonal=4.0*pow((n+1),2);
   this->tridiagonal=-1.0*pow((n+1),2);
   this->identity=-1.0*pow((n+1),2);
}

void PoissonMatrix::InitHashMatrix() {
    vector<int> tmp(5,-1);
    for(int i=0;i<this->dim;i++) {
        if(i<this->n) {
            if(i==0) {
                tmp[0]=i;
                tmp[1]=i+1;
                tmp[2]=i+this->n;
                tmp[3]=-1;
                tmp[4]=-1;
            } else if(i%this->n!=this->n-1 && i!=0) {
                tmp[0]=i-1;
                tmp[1]=i;
                tmp[2]=i+1;
                tmp[3]=i+this->n;
                tmp[4]=-1;
            } else {
                tmp[0]=i-1;
                tmp[1]=i;
                tmp[2]=i+this->n;
                tmp[3]=-1;
                tmp[4]=-1;
            }
        } else if(i>=this->n && i<this->dim-this->n) {
            if(i%this->n==0) {
                tmp[0]=i-this->n;
                tmp[1]=i;
                tmp[2]=i+1;
                tmp[3]=i+this->n;
                tmp[4]=-1;
            } else if(i%this->n==this->n-1) {
                tmp[0]=i-this->n;
                tmp[1]=i-1;
                tmp[2]=i;
                tmp[3]=i+this->n;
                tmp[4]=-1;
            } else {
                tmp[0]=i-this->n;
                tmp[1]=i-1;
                tmp[2]=i;
                tmp[3]=i+1;
                tmp[4]=i+this->n;
            }
        } else if(i>=this->dim-this->n) {
            if(i==this->dim-1) {
                tmp[0]=i-this->n;
                tmp[1]=i-1;
                tmp[2]=i;
                tmp[3]=-1;
                tmp[4]=-1;
            } else if(i%this->n!=0 && i!=this->dim-1) {
                tmp[0]=i-this->n;
                tmp[1]=i-1;
                tmp[2]=i;
                tmp[3]=i+1;
                tmp[4]=-1;
            } else {
                tmp[0]=i-this->n;
                tmp[1]=i;
                tmp[2]=i+1;
                tmp[3]=-1;
                tmp[4]=-1;
            }
        }
        HashMatrix.push_back(tmp);
    }
}

PoissonMatrix::~PoissonMatrix() {
    vector<vector<int> >().swap(HashMatrix);
}

int PoissonMatrix::Size() {
   return dim;
}

double PoissonMatrix::Get(int i,int j) {
   if(i==j) {
      return diagonal;
   } else if((j==(i+1) && i%n!=(n-1)) || (j==(i-1) && i%n!=0)) {
      return tridiagonal;
   } else if((j==(i+n)) || (j==(i-n))) {
      return identity;
   } else {
      return 0.0;
   }
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