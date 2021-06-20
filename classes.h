#ifndef CLASSES_H
#define CLASSES_H

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <cstdio>
#include <vector>
#include <assert.h>
#include <time.h>
#include <string>

using namespace std;

vector<double> operator-(const vector<double> &, const vector<double> &);
void operator+=(vector<double> &, const vector<double> &);
vector<double> operator*(double, vector<double>);
double operator|(const vector<double> &, const vector<double> &);
double operator*(const vector<double> &, const vector<double> &);

/*TODO Matrix WriteableMatrix PoissonMatrix UpperMatrix LowerMatrix*/

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
    int dim;
    double h;
    int Vcounter;

public:
    Algorithms(int);
    ~Algorithms();
    int SORMethod(vector<double> &, const vector<double> &, const vector<double> &);
    int MultiGridMethod(vector<double> &, const vector<double> &, const vector<double> &, string);
};

#endif