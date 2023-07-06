#ifndef FUNCTIONS_H
#define FUNCTIONS_H
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>
#include <fstream>

double MX(std::vector<double> &f);   

double DX(std::vector<double> &f);     

std::vector<double> fun_dence(std::vector<double> &f);

double quantile(std::vector<double> &f, double alpha);

void print(std::vector<double> &f, std::string p);

#endif // FUNCTIONS_H
