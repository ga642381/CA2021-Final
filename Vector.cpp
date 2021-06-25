#include "classes.h"

double Vectors::f(double x, double y, int k)
{
	
    double val;
    if (k == 1)
        val = -4.0;
    if (k == 2)
        val = 0.0;
    if (k == 3)
        val = 2*(6*pow(x,2)-6*x+1)*pow(y,2)*pow(y-1,2) + 2*(6*pow(y,2)-6*y+1)*pow(x,2)*pow(x-1,2);
    if (k == 4)
        val = 2*x*(x-1) + 2*y*(y-1);
    return val;
}

double Vectors::g(double x, double y, int k)
{
    double val;
    if (k == 1)
        val = pow(x, 2) + pow(y, 2);
    if (k == 2)
        val = 1.0;
    if (k == 3)
        val = pow(x,2)*pow(x-1,2)*pow(y,2)*pow(y-1,2);
    if (k == 4)
        val = x*(x-1)*y*(y-1);
    return val;
}

void Vectors::WriteToFile()
{
    FILE *file;
    file = fopen("./plot/dat/plot.dat", "w");
    if (file == NULL)
        printf("<td colspan=2>ERROR: Could not open file!</td>");
    else
    {
        for (int i = 0, k = 0; i <= sqrt(Size()); i++)
        {
            for (int j = 0; j <= sqrt(Size()); j++)
            {
                if (i == 0)
                {
                    continue;
                }
                else if (i != 0)
                {
                    if (j == 0)
                    {
                        continue;
                    }
                    else if (j != 0)
                    {
                        fprintf(file, "%f %f %f\n", (double)j / (double)(sqrt(Size())), (double)i / (double)(sqrt(Size())), Get(k));
                        k++;
                    }
                }
            }
            fprintf(file, "\n");
        }
    }
    fclose(file);
}

Boundary::Boundary(int n, int k)
{
    this->n = n;
    this->dim = n * n;
    this->h = 1.0 / (double)(n + 1);
    this->k = k;
    this->b.resize(dim);
    this->solved.resize(dim);

    for (int i = 1, k = 0; i <= n; i++)
    {
        for (int j = 1; j <= n; j++, k++)
        {
            solved[k] = g(j * h, i * h, this->k);
            b[k] = f(i * h, j * h, this->k);
            if (i == 1)
                b[k] += pow(1.0 / h, 2) * g(j * h, 0, this->k);
            if (i == n)
                b[k] += pow(1.0 / h, 2) * g(j * h, 1, this->k);
            if (j == 1)
                b[k] += pow(1.0 / h, 2) * g(0, i * h, this->k);
            if (j == n)
                b[k] += pow(1.0 / h, 2) * g(1, i * h, this->k);
        }
    }
}

Boundary::~Boundary()
{
    vector<double>().swap(b);
}

double Boundary::Get(int i)
{
    return this->b[i];
}

int Boundary::Size()
{
    return b.size();
}

Startvector::Startvector(int n, double value)
{
    this->n = n;
    this->dim = n * n;
    this->value = value;
    this->x.assign(dim, value);
}

Startvector::~Startvector()
{
    vector<double>().swap(x);
}

double Startvector::Get(int i)
{
    return this->x[i];
}

int Startvector::Size()
{
    return x.size();
}

/* === Vector Operators ===*/
vector<double> operator-(const vector<double> &lhs, const vector<double> &rhs)
{
    vector<double> tmp(lhs);
    for (int i = 0; i < (int)lhs.size(); i++)
    {
        tmp[i] = lhs[i] - rhs[i];
    }
    return tmp;
}

void operator+=(vector<double> &lhs, const vector<double> &rhs)
{
    for (int i = 0; i < (int)rhs.size(); i++)
    {
        lhs[i] += rhs[i];
    }
}

vector<double> operator*(double x, vector<double> rhs)
{
    vector<double> tmp(rhs);
    for (int i = 0; i < (int)rhs.size(); i++)
    {
        tmp[i] = x * rhs[i];
    }
    return tmp;
}

double operator|(const std::vector<double> &x, const std::vector<double> &y)
{
    double norm = 0.0;
    for (int i = 0; i < (int)x.size(); i++)
    {
        norm += x[i] * y[i];
    }
    return sqrt(norm);
}

double operator*(const std::vector<double> &x, const std::vector<double> &y)
{
    double ip = 0.0;
    for (int i = 0; i < (int)x.size(); i++)
        ip += x[i] * y[i];
    return ip;
}
