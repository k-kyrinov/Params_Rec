#include "matrix.h"

void Matrix::fill_arr()
{
    for (int i = 0; i < m; ++i){
        for (int j = 0; j < n; ++j){
        //    std::cin >> A[i][j];
    }
  }
}

void Matrix::show_arr()
{
    for (int i = 0; i < m; ++i){
        if(i!=0)std::cout << '\n';
        for (int j = 0; j < n; ++j){
            std::cout << A[i][j] << " ";
        }
    }
}

void Matrix::add_arr(Matrix a, Matrix b)
{
    for (int i = 0; i < m; i++){
        for (int j = 0; j < n; j++){
            A[i][j] = a.A[i][j] + b.A[i][j];
        }
    }
}

void Matrix::abs(){
    for(int i = 0; i < m; ++i)
        for(int j = 0; j < n; ++j)
            A[i][j] = std::abs(A[i][j]);
}

void Matrix::pow(int degree){
    double delta = 1e-5;
    for(int i = 0; i < m; ++i){
        for(int j = 0; j < n; ++j){
            if(degree < 0 and A[i][j] < delta) A[i][j] = std::pow(delta, degree);
            else A[i][j] = std::pow(A[i][j],degree);
        }
    }
}

Matrix Matrix::mult_arr(Matrix a, Matrix b)
{
    if(a.n == b.m){
        Del();
        m = a.m;
        n = b.n;
        Create();
        for (int i = 0; i < a.m; i++){
            for (int j = 0; j < b.n; j++){
                //A[i][j] = 0;
                for (int k = 0; k < a.n; k++)
                    A[i][j] += a.A[i][k] * b.A[k][j];
            }
        }
    }
    else return a;
}

Matrix Matrix::difference(Matrix a, Matrix b)
{
    if(a.m == b.m and a.n == b.n){
        Del();
        m = a.m;
        n = a.n;
        Create();
        for(int i = 0; i < m; ++i){
            for(int j = 0; j < n; ++j){
                A[i][j] = a.A[i][j] - b.A[i][j];
            }
        }
    }
    else return a;
}

int Matrix::Get_m()
{
    return m;
};

double Matrix::sum()
{
    double sum{};
    for(int i = 0; i < m; ++i){
        for(int j = 0; j < n; ++j)
            sum+=A[i][j];
    }
    return sum;
}

double Matrix::max(){
    double max{};
    for(int i = 0; i < m; ++i)
        for(int j = 0; j < n; ++j)
            if(A[i][j]>max) max = A[i][j];
    return max;
}

int Matrix::Get_n()
{
    return n;
};

void Matrix::Del()
{
   for (int i = 0; i < m; i++)
   delete[] A[i];
   delete[] A;
}

double Matrix::shpur()
{
    double sh{};
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < n; ++j)
            if(i==j)sh+=A[i][j];
    return sh;
}

double Matrix::Dispersion(Matrix a){              // дисперсия
    double Summ{},Summ2{}, size{};
        for (int i = 0; i < m; ++i){
            if(a.A[i][0]>0) size+=1;
            Summ += a.A[i][0];
            Summ2 += a.A[i][0]*a.A[i][0];
        }
        return (1.0/(size-1))*Summ2 - (1.0/(size*(size-1)))*Summ*Summ;
}

double Matrix::Mean(Matrix a){
    double Summ{};
    for (int i = 0;i < m; i++) Summ+=a.A[i][0];
    return Summ/m;
}

void Matrix::merge(Matrix a, Matrix b)
{
    Del();
    n = a.n + b.n;
    m = a.m;
    Create();
    for (int i = 0; i < m; i++){
        for (int j = 0; j < n; j++){
        if(j<a.n) A[i][j] = a.A[i][j];
        else A[i][j] = b.A[i][j-a.n];
        }
    }
}

double Matrix::Get_el(int i, int j)
{
    if(i <= m-1 and j <= n-1) return A[i][j];
}

void Matrix::Set_el(int i, int j, double el)
{
    if(i <= m-1 and j <= n-1) A[i][j] = el;
}

void Matrix::unit()
{
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            if(i==j) A[i][j] = 1;
            else A[i][j] = 0;
        }
    }
}

void Matrix::diag(Matrix a)  // create diagonal matrix from array
{
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            if(i==j and a.m > a.n) {
                A[i][j] = a.A[i][0];
            }
            else if (i==j and a.m < a.n){
                A[i][j] = a.A[0][i];
            }
            else A[i][j] = 0;
        }
    }
}

void Matrix::Add_to_el(int i, int j, double el)
{
    A[i][j]+=el;
}

void Matrix::Mult_on_el(int i, int j, double w)
{
/*    for(int i = 0; i < m; ++i){
        for(int j = 0; j < n; ++j){
            if(A[i][j]!=0) A[i][j]*=w;
        }
    }*/
    A[i][j]*=w;
}

void Matrix::reverse()
{
    double **I;
    I = new double*[m];
    for (int i = 0; i < m; ++i)
    I[i] = new double[m];
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            if(i==j) I[i][j] = 1;
            else I[i][j] = 0;
        }
    }

    double tem{};
    tem = A[0][0];
    for (int kk = 0; kk < m; kk++) {
        A[0][kk]/=tem;
        I[0][kk]/=tem;
    }

    int k{};
    double temp{};

    // Прямой ход

    for (int s = 1; s < m; s++) {
        if(A[s][s]!=0){
        for (int i = s; i < m; i++) {
            ++k;
            temp = A[i][s-1];
            for (int j = 0; j < m; j++) {
                A[i][j] = A[i][j] - A[i-k][j]*temp;
                I[i][j] = I[i][j] - I[i-k][j]*temp;
            }
        }
        if(A[s][s]!=0){
            tem = A[s][s];
            for (int kk = 0; kk < m; kk++) {
                A[s][kk]/=tem;
                I[s][kk]/=tem;
            }
        }
        k=0;
        }
    }

    // Обратный ход

    for (int s = m-2; s >= 0 ; s--) {
        for (int i = s; i >= 0; i--) {
            ++k;
            temp = A[i][s+1];
            for (int j = m-1; j >= 0; j--) {
                A[i][j] = A[i][j] - A[i+k][j]*temp;
                I[i][j] = I[i][j] - I[i+k][j]*temp;
            }
        }
        k=0;
      }

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            A[i][j] = I[i][j];

        }
    }
}

void Matrix::LU(char x)
{
    double **L, **U;
    L = new double*[n];
    U = new double*[n];
    for (int i = 0; i < n; i++){
    L[i] = new double[n];
    U[i] = new double[n];
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            U[i][j] = 0;
            L[i][j] = 0;
        }
    }
    for (int j = 0; j < n; j++) {
        U[0][j] = A[0][j];
        if(U[0][0]!=0) L[j][0] = A[j][0]/U[0][0];
        //else qDebug() << "Division on zero" << '\n'; //std::cout << "Division on zero" << std::endl;
    }
    double pr_summ1{}, pr_summ2{};
    int i{};
    for (i = 1; i < n; i++) {
        for (int j = i; j < n; j++) {
        for (int k = 0; k <= i-1; k++) {
            pr_summ1+=L[i][k]*U[k][j];
            pr_summ2+=L[j][k]*U[k][i];
            }
        U[i][j] = A[i][j] - pr_summ1;
        if(U[i][i]!=0)L[j][i] = (A[j][i] - pr_summ2)/U[i][i];
        pr_summ1 = pr_summ2 = 0;
        //std::cout << U[i][j] << "   " << L[i][j] << std::endl;
        //std::cout << pr_summ1 << "   " << pr_summ2 << std::endl;
        }
    }
    if(x=='U')
    {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A[i][j] = U[i][j];
        }
    }
    }
    if(x=='L')
    {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A[i][j] = L[i][j];
        }
    }
    }
}

double Matrix::det()
{
    double determ = 1;
    this->LU('U');
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if(i==j)determ*=A[i][j];
        }
    }
    return determ;
}

void Matrix::clear()
{
    for(int i = 0; i < m; ++i)
        for(int j = 0; j < n; ++j)
            A[i][j] = 0;
}
