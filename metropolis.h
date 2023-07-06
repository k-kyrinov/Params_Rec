#ifndef METROPOLIS_H
#define METROPOLIS_H
#include <random>
#include <vector>
#include <fstream>
#include <algorithm>
#include <iterator>

class Metropolis        // Hyperparameters for finding maximum  - Number of iterations (N) and step value (eps)
{
public:
    Metropolis(std::vector<double> &x1, std::vector<double> &y1, std::vector<int> &xt1, std::vector<double> &yt, std::string &s1) : x(x1), y(y1), t(xt1), y_t(yt)
    {
        std::ifstream conf(s1.c_str()); //config_rect config
        double x_d{}, y_d{};
        while(conf >> x_d && conf >> y_d){
            x_det.push_back(x_d);
            y_det.push_back(y_d);
        }
        if(x_det.size()==y_det.size()) Num_det = x_det.size();
    }
    ~Metropolis()
    {
    x.clear();
    y.clear();
    t.clear();
    y_t.clear();
    }

    std::vector<double> &x;
    std::vector<double> &y;     // TRUE values
    std::vector<int> &t;
    std::vector<double> &y_t;
    void start_init();
    void find_min(int N);
    void find_grad(int N);
    double num_grad(int &i);
    void find_minf(int N);
    void arr_dir();
    double calc_del_r();
    double calc_psi();
    double nkg(int &i);
    double nkg2(int &i);
    double Cors_nkg(int &i);
    double Lagutin(int &i);
    double Lagutin_mod(int &i);
    double Kascade(int &i);
    double Greisen(int &i);
    double Uchakin(int &i);
    // ------------ Tibet attitude
    double nkg_t(int &i);
    double Tien_Shan(int &i);
    double Argo(int &i);
    double LHAASO(int &i);

    double ne_calc();
    double L();
    double L2();
    void print_params();
    std::vector<double> par;
    std::vector<double> angl;

private:
    double S = 0.36;
    int threshold = 7;                      // treshold of trigged detectors
    double pi = 3.14159265358979323846;
    double c = 0.2998; // m/nsec
    std::random_device rd;
    int Num_det{}; /*= 16;*/
    std::vector<double> x_det;
    std::vector<double> par_prom;
    std::vector<double> y_det;
    std::vector<double> eps {5, 5, 0.1, 1000}; // eps declaration {x, y, s, Ne}; 50k
    std::vector<double> eps_g {1, 1, 0.1, 5000}; // eps declaration {x, y, s, Ne}; 50k
    //std::string s; // filename
};

#endif // METROPOLIS_H
