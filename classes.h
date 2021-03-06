#ifndef CLASSES_H
#define CLASSES_H

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
// #include <cstdio>
#include <stdio.h>
#include <vector>
#include <assert.h>
#include <time.h>
#include <string>
#include <math.h>
#include <omp.h>

using namespace std;

vector<double> operator-(const vector<double> &, const vector<double> &);
void operator+=(vector<double> &, const vector<double> &);
vector<double> operator*(double, vector<double>);
double operator|(const vector<double> &, const vector<double> &);
double operator*(const vector<double> &, const vector<double> &);
void print_vector(const vector<double> &);

class Matrix
{
public:
    vector<vector<int>> HashMatrix;
    virtual int Size() = 0;
    virtual vector<double> operator*(const vector<double> &) = 0;
};

class WriteableMatrix : public Matrix
{
public:
    virtual void Set(int, int, double) = 0;
};

class PoissonMatrix : public Matrix
{
private:
    int dim;
    int n;

public:
    PoissonMatrix(int);
    ~PoissonMatrix();
    int Size();
    vector<double> operator*(const vector<double> &);
};

/*================================Matrix end here===================================================================*/
class Vectors
{
public:
    virtual double Get(int) = 0;
    virtual int Size() = 0;
    void WriteToFile();
    double f(double, double, int);
    double g(double, double, int);
};

class Boundary : public Vectors
{
private:
    int dim;
    int n;
    double h;
    int k;

public:
    Boundary(int, int);
    ~Boundary();
    vector<double> b;
    vector<double> solved;
    double Get(int);
    int Size();
};

class Startvector : public Vectors
{
private:
    int dim;
    int n;
    double value;

public:
    Startvector(int, double);
    ~Startvector();
    vector<double> x;
    double Get(int);
    int Size();
};

class Algorithms
{
private:
    int n;
    int num_threads;
    int dim;
    double h;
    int Vcounter;

public:
    int update_step;
    vector<double> error_array;
    Algorithms(int, int);
    ~Algorithms();
    int SORMethod(vector<double> &, const vector<double> &, const vector<double> &, int fixed_step);
    int MultiGridMethod(Matrix &, vector<double> &, const vector<double> &, const vector<double> &, string, int fixed_step);
    void JacobiMethod(Matrix &, vector<double> &, const vector<double> &, const vector<double> &);
    vector<double> Cycle(Matrix &, vector<double> &, const vector<double> &, int, int);
    void JacobiRelaxation(Matrix &, vector<double> &, const vector<double> &, int);
    void SORRelaxation(Matrix &, vector<double> &, const vector<double> &, int);
    void Prolongation(const vector<double> &, vector<double> &, int);
    void Restriction(const vector<double> &, vector<double> &, int);
};

#endif
