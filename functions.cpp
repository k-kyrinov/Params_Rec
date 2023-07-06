#include"functions.h"

double MX(std::vector<double> &f){                       // Expectation value
    double Summ{};
    for (int i = 0;i < f.size(); ++i) Summ+=f[i];
    return Summ/f.size();
}

double DX(std::vector<double> &f){               // Dispersion
    double Summ{},Summ2{};
    for (int i = 0;i < f.size(); ++i){
        Summ += f[i];
        Summ2 += f[i]*f[i];
    }
    return (1.0/(f.size()-1))*Summ2 - (1.0/(f.size()*(f.size()-1)))*Summ*Summ;
}

double quantile(std::vector<double> &f, double alpha){
    std::sort(f.begin(), f.end());
    const std::size_t pos = alpha * std::distance(f.begin(), f.end());
    return f[pos];
}

std::vector<double> fun_dence(std::vector<double> &f){
    std::vector<double> f_res(50);
    for(int i = 0; i < 50; ++i){
        for(int j = 0; j < f.size(); ++j)
            if(f[j]<i+1 and f[j]>i) f_res[i]+=1;
    }
    for(int i = 0; i < f_res.size(); ++i)
        //f_res[i]/=std::accumulate(f_res.begin(), f_res.end(), 0);
        f_res[i]/=f.size();
    return f_res;
}

void print(std::vector<double> &f, std::string p){
    //std::ofstream out("C:\\Qt\\progects\\markov\\out.txt");
    std::ofstream file_del(p.c_str());
    for(int i = 0; i < f.size(); ++i) file_del << f[i] << '\n';
}
