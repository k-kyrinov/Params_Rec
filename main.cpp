#include "functions.h"
#include "metropolis.h"
#include <random>
#include <iostream>
#include <future>
#include <thread>
#include <ctime>
#include <mutex>
#include <sstream>

const int Num_det = 16;
const int number_of_threads = 1;

std::vector <double> del_R;
std::vector <double> del_psi;
std::vector<std::thread> thread_pool;

std::mutex mtx;
std::ofstream out("C:\\Qt\\progects\\MC_EAS\\pred.txt");
std::ofstream out1("C:\\Qt\\progects\\MC_EAS\\mod.txt");
std::string s;

void loop3(int j, std::vector<double> &xx, std::vector<double> &yy, std::vector<int> &xx_t ,std::vector<double> &yy_t)
{
    while (j < xx.size()/Num_det){
    std::vector<double> x1(&xx[Num_det*j], &xx[Num_det*(j+1)]), y1(&yy[2*j], &yy[2*(j+1)]), yt1(&yy_t[2*j], &yy_t[2*(j+1)]);
    std::vector<int> xt1(&xx_t[Num_det*j], &xx_t[Num_det*(j+1)]);
    Metropolis M(x1, y1, xt1 ,yt1, s);
    M.start_init();
    M.arr_dir();
    M.find_min(10000); //800
    //M.find_grad(1000);
//    M.find_minf(10000);
//    if(!(M.par[0] > 9.5 or M.par[0] < -9.5 or M.par[1] > 9.5 or M.par[1] < -9.5) and M.par[2] > 0.2 and M.par[3] > 1 and M.par[2] < 2.0  /*and M.angl[0] < 40 and (M.angl[0]!=-1 && M.angl[1]!=-1 )*/){      // отбор внутри границ установки
    del_R.push_back(M.calc_del_r());
    if(M.calc_psi()!=-1) del_psi.push_back(M.calc_psi());
    mtx.lock();
    out << M.par[0] << ','<< M.par[1] << ','<< M.par[2] << ',' << M.par[3] << ',' << M.ne_calc() << ',' << M.calc_psi() << ',' << M.calc_del_r() << '\n';
    out1 << M.y[0] << ',' << M.y[1] << ',' << M.angl[0] << ',' << M.angl[1] << ',' << M.y_t[0] << ',' << M.y_t[1] << ',' << M.calc_psi() << ',' << M.calc_del_r()  << '\n';
    //std::cout << M.par[0] << ' ' << M.par[1] << ' ' << M.y[0] << ' ' << M.y[1] << ' ' <<  M.calc_del_r() << '\n';
    mtx.unlock();
//    }
    std::cout << j << ' ' << M.calc_del_r() << ' ' << /*M.angl[0] << ' ' << M.y_t[0] << ' ' <<*/ M.calc_psi()  << '\n';
//    std::cout << "0.72 quantile: " << quantile(del_R, 0.72) << '\n';
    j+=number_of_threads;
    }
}

int main(int argc, char **argv)
{
//    if(argc != 3){
//    std::cout << "Please choose config and data file" << '\n';
//    return 0;
//    }
    std::fstream f1("C:\\Qt\\progects\\MC_EAS\\X1_new.txt");    //"C:\\Qt\\progects\\MC_EAS\\input_data.txt"
    // std::fstream f1("C:\\Qt\\progects\\MC_EAS\\input_data.txt");    //"C:\\Qt\\progects\\MC_EAS\\input_data.txt"
    s = "C:\\Qt\\progects\\MC_EAS\\conf_16_INR.txt"; //std::string(argv[2]);
    //s = "C:\\Qt\\progects\\MC_EAS\\conf_16_HZS.txt"; //std::string(argv[2]);

    std::vector<double> xx, yy, yy_t;
    std::vector<int> xx_t;
    std::string line;
    std::vector<std::string> parsedRow2;
    int count{};
    while(std::getline(f1,line))
    {
        std::stringstream lineStream(line);
        std::string cell;
        std::vector<std::string> parsedRow;
        while(std::getline(lineStream,cell, ',') ) parsedRow.push_back(cell);
        for(int i = 0; i < 2*Num_det + 4; ++i)
        parsedRow2.push_back(parsedRow[i]);
        ++count;
    }


    for(int i = 0; i < parsedRow2.size(); ++i){
        if(i%(2*Num_det + 4) ==0){
            for(int k = 0; k < 2*Num_det + 4; ++k){
                if(k < Num_det) xx.push_back(stod(parsedRow2[k+i]));
                else if(k < 2*Num_det and k >= Num_det)  xx_t.push_back(stod(parsedRow2[k+i]));
                else if (k >= 2*Num_det and k < 2*Num_det+2 ) yy.push_back(stod(parsedRow2[k+i]));
                else yy_t.push_back(stod(parsedRow2[k+i]));
            }
        }
    }

    for(int i = 0; i < number_of_threads; ++i) thread_pool.emplace_back(loop3, i, std::ref(xx), std::ref(yy), std::ref(xx_t), std::ref(yy_t));
    for(int i = 0; i < number_of_threads; ++i){
        std::thread& entry = thread_pool[i];
        entry.join();
    }

    std::cout << "Number of events: " << count << '\n';
    std::cout << "Events in border: " << del_R.size() << ' ' << static_cast<double>(del_R.size())/count*100 << '\n'; //count
    std::cout << "Distance error: " << '\n';
    std::cout << "Mean "  << MX(del_R) << ' ' << " Mean Square Deviation " << std::sqrt(DX(del_R)) << '\n';
    std::cout << "0.72 quantile: " << quantile(del_R, 0.72) << '\n';

    print(del_R, "C:\\Qt\\progects\\MC_EAS\\del_R.txt");

    std::cout << "Events direc rec: " << del_psi.size() << ' ' << static_cast<double>(del_psi.size())/count*100 << '\n';
    std::cout << "psi error: " << '\n';
    std::cout << "Mean "  << MX(del_psi) << ' ' << " Mean Square Deviation " << std::sqrt(DX(del_psi)) << '\n';
    std::cout << "0.72 quantile: " << quantile(del_psi, 0.72) << '\n';
    //print(del_psi, "C:\\Qt\\progects\\MC_EAS\\psi.txt");

    std::vector<double> den;
    den = fun_dence(del_R);
    print(den, "C:\\Qt\\progects\\MC_EAS\\dence.txt");
    return 0;
}
