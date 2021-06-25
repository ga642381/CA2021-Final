#include "classes.h"

Algorithms::Algorithms(int n)
{
    this->n = n;
    this->dim = n * n;
    this->h = 1.0 / (double)(n + 1);
    this->Vcounter = 0;
}

Algorithms::~Algorithms() {}

int Algorithms::MultiGridMethod(Matrix& A,vector<double> &x, const vector<double> &b, const vector<double> &solved, string method)
{
    // MultiGrid Method
    // To be implemented
    // input should be Matrix
    int steps = 0, dim = x.size();
    vector<double> r(dim);
    r = x - solved;
    int numberOfGrids = 2, VW = 0;
    if (method == "Two-Grid")
    {
        numberOfGrids = 1;
        VW = 0;
    }
    
    if (method == "V-Cycle")
    {
        numberOfGrids = 2; 
        VW = 0;
    }

    
    if (method == "W-Cycle")
    {
        numberOfGrids = 2;
        VW = 1;
    }
    
    double TOL = pow(10, -3) * (r | r);
    while (TOL <= (r | r)) 
    {
        x = Cycle(A, x, b, numberOfGrids, VW, solved);
        r = x - solved;
        steps++;
    }
    return steps;
}

vector<double> Algorithms::Cycle(Matrix& A, vector<double>& x, const vector<double>& b, int lambda, int theta, const vector<double>& solved) {
    int dim = x.size(), n = sqrt(dim), N2h = (n + 1) / 2 - 1, dim2h = pow(N2h, 2);
    vector<double> r(dim, 0), E(dim, 0), sol(dim, 0), r2h(dim2h, 0), E2h(dim2h, 0), sol2h(dim2h, 0);
    if (this->Vcounter == lambda) {
        //JacobiMethod(A, x, b, solved);
        //SORMethod(x, b, solved);
        JacobiRelaxation(A, x, b, 200);
        return x;
    }
    else {
        this->Vcounter++;
        JacobiRelaxation(A, x, b, 3);
        r = b - A * x;
        sol = solved;
        Restriction(r, r2h, n);
        Restriction(sol, sol2h, n);
        E2h = Cycle(A, E2h, r2h, lambda, theta, sol2h);
        if (theta == 1) E2h = Cycle(A, E2h, r2h, lambda, theta, sol2h);
        Interpolation(E2h, E, n);
        x += E;
        JacobiRelaxation(A, x, b, 3);
        this->Vcounter--;
        return x;
    }
}

void Algorithms::JacobiMethod(Matrix& A,vector<double>& x,const vector<double>& b,const vector<double>& solved) {
    vector<double> tmp;
    double sum;
    int n=sqrt(x.size()),dim=n*n,steps=0;
    double h=1.0/(double)(n+1);
    vector<double> r(x.size());
    r=x-solved;
    double TOL=(r|r)*pow(10,-3);
    while(TOL<=(r|r)) {
        cout << (r|r)/TOL << endl;
        tmp=x;
        for(int i=0;i<n*n;i++) {
            sum=0.0;
            if(i>=n && i<dim-n) sum+=tmp[i-n]+tmp[i+n];
            if(i<n) sum+=tmp[i+n];
            if(i>=dim-n) sum+=tmp[i-n];
            if(i%n!=0) sum+=tmp[i-1];
            if(i%n!=n-1) sum+=tmp[i+1];
            x[i]=1.0/(4.0*pow(1.0/h,2))*(b[i]+pow(1.0/h,2)*sum);
            r[i]=x[i]-solved[i];
        }
        steps++;
    }
    //return steps;
}

void Algorithms::JacobiRelaxation(Matrix& A, vector<double>& x, const vector<double>& b, int steps) {   //a,x,b,3
    vector<double> tmp;
    tmp = x;
    double sum, omega = 4.0 / 5.0;
    int n = sqrt(x.size()), dim = n * n;
    for (int k = 0;k < steps;k++) {
        for (int i = 0;i < n * n;i++) {
            sum = 0.0;
            if (i >= n && i < dim - n) sum += tmp[i - n] + tmp[i + n];
            if (i < n) sum += tmp[i + n];
            if (i >= dim - n) sum += tmp[i - n];
            if (i % n != 0) sum += tmp[i - 1];
            if (i % n != n - 1) sum += tmp[i + 1];
            x[i] = omega * 1.0 / (4.0 * pow(1.0 / h, 2)) * (b[i] + pow(1.0 / h, 2) * sum) + tmp[i] * (1 - omega);
            tmp[i] = x[i];
        }
    }
}

void Algorithms::Restriction(const vector<double>& r,vector<double>& r2h,int n) {
    for(int i=1,l=0,k=0;i<=n;i++) {
        for(int j=1;j<=n;j++,k++) {
            if(i%2==0 && j%2==0) {
                r2h[l]=1.0/16.0*(4.0*r[k]+2.0*(r[k-1]+r[k+1]+r[k-n]+r[k+n])+r[k+n-1]+r[k+n+1]+r[k-n-1]+r[k-n+1]);
                l++;
            }
        }
    }
}


void Algorithms::Interpolation(const vector<double>& E2h, vector<double>& E, int n) {
    for (int i = 1, k = 0, l = 0;i <= n;i++) {
        for (int j = 1;j <= n;j++, k++) {
            if (i % 2 == 0 && j % 2 == 0) {
                E[k] = E2h[l];
                l++;
            }
        }
    }
    for (int i = 1, k = 0;i <= n;i++) {
        for (int j = 1;j <= n;j++, k++) {
            if (i % 2 == 0 && j % 2 != 0) {
                if (j != 1 && j != n) E[k] = 1.0 / 2.0 * (E[k - 1] + E[k + 1]);
                if (j == 1) E[k] = 1.0 / 1.0 * (E[k + 1]);
                if (j == n) E[k] = 1.0 / 1.0 * (E[k - 1]);
            }
            if (i % 2 != 0 && j % 2 == 0) {
                if (i != 1 && i != n) E[k] = 1.0 / 2.0 * (E[k - n] + E[k + n]);
                if (i == 1) E[k] = 1.0 / 1.0 * (E[k + n]);
                if (i == n) E[k] = 1.0 / 1.0 * (E[k - n]);
            }
            if (i % 2 != 0 && j % 2 != 0) {
                if (i != 1 && j != 1 && i != n && j != n) E[k] = 1.0 / 4.0 * (E[k + n - 1] + E[k + n + 1] + E[k - n + 1] + E[k - n - 1]);
                if (i == 1 && j == 1) E[k] = 1.0 / 1.0 * (E[k + n + 1]);
                if (i == n && j == n) E[k] = 1.0 / 1.0 * (E[k - n - 1]);
                if (i == 1 && j == n) E[k] = 1.0 / 1.0 * (E[k + n - 1]);
                if (i == n && j == 1) E[k] = 1.0 / 1.0 * (E[k - n + 1]);
                if (i == 1 && j != 1 && j != n) E[k] = 1.0 / 2.0 * (E[k + n + 1] + E[k + n - 1]);
                if (i == n && j != 1 && j != n) E[k] = 1.0 / 2.0 * (E[k - n - 1] + E[k - n + 1]);
                if (j == 1 && i != 1 && i != n) E[k] = 1.0 / 2.0 * (E[k + n + 1] + E[k - n + 1]);
                if (j == n && i != 1 && i != n) E[k] = 1.0 / 2.0 * (E[k + n - 1] + E[k - n - 1]);
            }

        }
    }
}

int Algorithms::SORMethod(vector<double> &x, const vector<double> &b, const vector<double> &solved)
{

    double sum;
    double Pi = 3.141592654;
    double omega = 2 / (1 + sqrt(1 - pow(cos(Pi * h), 2)));

    int n = sqrt(x.size());
    int steps = 0;
    double h = 1.0 / (double)(n + 1);

    // residual
    vector<double> r(x.size());
    r = x - solved;
    double TOL = (r | r) * pow(10, -3);

    // update loop
    while (TOL <= (r | r))
    {
        // iterate each cell
        for (int i = 0; i < n * n; i++)
        {
            sum = 0.0;
            if (i >= n && i < dim - n)
                sum += x[i - n] + x[i + n];
            if (i < n)
                sum += x[i + n];
            if (i >= dim - n)
                sum += x[i - n];
            if (i % n != 0)
                sum += x[i - 1];
            if (i % n != n - 1)
                sum += x[i + 1];
            x[i] = omega * 1.0 / (4.0 * pow(1.0 / h, 2)) * (b[i] + pow(1.0 / h, 2) * sum) + x[i] * (1 - omega);
            r[i] = x[i] - solved[i];
        }
        steps++;
    }
    return steps;
}



