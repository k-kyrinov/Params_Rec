#ifndef MATRIX_H
#define MATRIX_H
#include <cmath>
#include <fstream>
#include <iostream>

class Matrix{
private:
    double **A;
    int m;
    int n;
    void Create(){
    A = new double*[m];
    for (int i = 0; i < m; ++i)
    A[i] = new double[n];
    zeros();
    }
    void zeros(){
        for (int i = 0; i < m; i++){
            for (int j = 0; j < n; j++){
                A[i][j] = 0;
            }
        }
    }
public:
    Matrix() :  A(nullptr), m(0), n(0)
    {}
    Matrix(int m1, int n1)
    {
    m = m1;
    n = n1;
    Create();
    }
    Matrix(const Matrix& _M)    //конструктор копирования
    {
    m = _M.m;
    n = _M.n;
    Create();
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            A[i][j] = _M.A[i][j];
    }

    double **Get_arr(){
        return A;
    }
    void fill_arr();
    void show_arr();
    double sum();
    double max();
    void add_arr(Matrix a, Matrix b);
    void clear();
    int Get_m();
    int Get_n();
    double Get_el(int i, int j);
    void Set_el(int i, int j, double el);
    void Del();
    double shpur();
    double Dispersion(Matrix a);
    double Mean(Matrix a);
    void merge(Matrix a, Matrix b);
    void unit();
    void diag(Matrix a);
    void Add_to_el(int i, int j, double el);
    void Mult_on_el(int i, int j, double w);
    void reverse();
    void abs();
    void pow(int degree);
    void LU(char x);
    double det();
    Matrix difference(Matrix a, Matrix b);
    Matrix mult_arr(Matrix a, Matrix b);

    Matrix operator=(const Matrix& _M)
    {
        Del();
        m = _M.m;
        n = _M.n;
        Create();
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                A[i][j] = _M.A[i][j];
        return *this;                   //скрытый указатель (содержит адресс object)
        }

    Matrix T(const Matrix& _M){
        if(_M.m == _M.n){
            double temp{};
            for (int i = 0; i < m; ++i){
                for (int j = 0; j < n; ++j){
                    if(i!=j && i < j){
                        temp = _M.A[i][j];
                        A[i][j] = _M.A[j][i];
                        A[j][i] = temp;
                    }
                }
            }
        }
        else{
            Del();
            m = _M.n;
            n = _M.m;
            Create();
            for(int i = 0; i < _M.m; ++i){
                for(int j = 0; j < _M.n; ++j){
                    A[j][i] = _M.A[i][j];
                }
            }
        }
        return *this;
    }

//    ~Matrix()
//    {
//      if (n > 0)
//      {
//        for (int i = 0; i < m; i++)
//          delete[] A[i];
//      }
//      if (m > 0)delete[] A;
//      }
};

#endif // MATRIX_H
