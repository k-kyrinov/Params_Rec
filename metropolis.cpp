#include "metropolis.h"
#include "functions.h"
#include "matrix.h"
#include <iostream>
#include <mutex>

double x_a[51] = {1.97376115143293e-001,2.39192422628323e-001,2.98565164093573e-001,3.58271550162608e-001,4.29917884243487e-001,
                  5.15891890127776e-001,6.19058783209113e-001,7.42856758173555e-001,8.91411572102221e-001,1.06967404164358e+000,
                  1.28358503655918e+000,1.54027346830543e+000,1.84829387192396e+000,2.21791149902107e+000,2.66144442299609e+000,
                  3.19367405769951e+000,3.83233777068348e+000,4.59872000813587e+000,5.51836163164124e+000,6.62191110650248e+000,
                  7.94614590877754e+000,9.53519819098440e+000,1.14420255536611e+001,1.37301759384948e+001,1.64759054607864e+001,
                  1.97707197612597e+001,2.37244235716564e+001,2.84687801255600e+001,3.41619023699168e+001,4.09935223211034e+001,
                  4.91913141748989e+001,5.90284818976854e+001,7.08328641669701e+001,8.49978600973166e+001,1.01995539868228e+002,
                  1.22392377185738e+002,1.46868127886072e+002,1.76238483839769e+002,2.11482257132292e+002,2.53774000475579e+002,
                  3.04523150975701e+002,3.66021651167416e+002,4.39217832914120e+002,5.27051621494212e+002,6.32450212407458e+002,
                  7.58926175087446e+002,9.10694514656579e+002,1.09281314342600e+003,1.31135144356828e+003,1.57359253857203e+003,
                  1.88827601448448e+003};
double y_a[51] = {2.82407407407407e-001,2.74074074074074e-001,2.63477366255144e-001,2.54732510288066e-001,2.44444444444444e-001,
                  2.34362139917696e-001,2.23971193415638e-001,2.11316872427984e-001,1.99279835390947e-001,1.86934156378601e-001,
                  1.72222222222222e-001,1.58539094650206e-001,1.44444444444444e-001,1.29115226337449e-001,1.14814814814815e-001,
                  1.00925925925926e-001,8.61111111111111e-002,7.31481481481482e-002,6.12139917695473e-002,4.97942386831276e-002,
                  4.03292181069959e-002,3.24074074074074e-002,2.61316872427984e-002,2.20164609053498e-002,2.04732510288066e-002,
                  2.08847736625514e-002,2.36625514403292e-002,2.91152263374486e-002,3.73456790123457e-002,4.71193415637860e-002,
                  5.92592592592593e-002,7.40740740740741e-002,8.95061728395062e-002,1.06378600823045e-001,1.26131687242798e-001,
                  1.44958847736626e-001,1.64403292181070e-001,1.85699588477366e-001,2.05144032921811e-001,2.24485596707819e-001,
                  2.44753086419753e-001,2.62551440329218e-001,2.79526748971193e-001,2.96913580246914e-001,3.11213991769547e-001,
                  3.25205761316872e-001,3.39711934156379e-001,3.51851851851852e-001,3.63683127572016e-001,3.75205761316872e-001,
                  3.85288065843621e-001};

double Metropolis::nkg(int &i){
    double rm = 80.0;
    double delt = 0.1;
    double r = std::sqrt(std::pow(par[0]-x_det[i], 2) + std::pow(par[1]-y_det[i], 2)) /*+ 0.1*/;
    //if(r < 0.1) r = 0.1;
    double alpha_r{};
    for(int i = 0; i < 51-1; ++i){
        if(r>=x_a[i] && r<x_a[i+1]){
        alpha_r = (y_a[i]-y_a[i+1])/(x_a[i]-x_a[i+1])*r + y_a[i]-x_a[i]*(y_a[i]-y_a[i+1])/(x_a[i]-x_a[i+1]);
       }
    }
    return  (1/rm/rm)*std::tgamma(4.5-par[2])/(2*pi*tgamma(par[2])*tgamma(4.5-2*par[2]))                   // age with corresctions
            *(std::pow((r/rm),(par[2]+alpha_r-2.)))*std::pow((1+r/rm),(par[2]+alpha_r-4.5));

//    return (1/rm/rm)*std::tgamma(4.5-par[2])/(2*pi*tgamma(par[2])*tgamma(4.5-2*par[2]))                   // standard NKG
//           *(std::pow((r/rm),(par[2]-2.)))*std::pow((1+r/rm),(par[2]-4.5));
//    return (1/rm/rm)*(0.366)*par[2]*par[2]*std::pow((2.07-par[2]),1.25)                                     // Modified nkg
//           *(std::pow((r/rm),(par[2]-2.)))*std::pow((1+r/rm),(par[2]-4.5));

}

double Metropolis::nkg2(int &i){
    double rm = 80.0;
    double delt = 0.1;
    double r = std::sqrt(std::pow(par[0]-x_det[i], 2) + std::pow(par[1]-y_det[i], 2)) /*+ 0.1*/;
    if(r < 0.1) r = 0.1;
    double alpha_r{};
    for(int i = 0; i < 51-1; ++i){
        if(r>=x_a[i] && r<x_a[i+1]){
        alpha_r = (y_a[i]-y_a[i+1])/(x_a[i]-x_a[i+1])*r + y_a[i]-x_a[i]*(y_a[i]-y_a[i+1])/(x_a[i]-x_a[i+1]);
       }
    }
    return  (1/rm/rm)*std::tgamma(4.5-par[2])/(2*pi*tgamma(par[2])*tgamma(4.5-par[2]))                   // age with corresctions
            *(std::pow((r/rm),(par[2]+alpha_r-2.)))*std::pow((1+r/rm),(par[2]+alpha_r-4.5));

//    return (1/rm/rm)*std::tgamma(4.5-par[2])/(2*pi*tgamma(par[2])*tgamma(4.5-par[2]))                   // standard NKG
//           *(std::pow((r/rm),(par[2]-2.)))*std::pow((1+r/rm),(par[2]-4.5));
//    return (1/rm/rm)*(0.366)*par[2]*par[2]*std::pow((2.07-par[2]),1.25)*std::tgamma(4.5-par[2])/(2*pi*tgamma(par[2])*tgamma(4.5-2*par[2])) // Modified
//           *(std::pow(((r+delt)/rm),(par[2]-2.)))*std::pow((1+(r+delt)/rm),(par[2]-4.5));

}

double Metropolis::Cors_nkg(int &i){
    double rm = 80.0;
    double delt = 0.1;
    double r = std::sqrt(std::pow(par[0]-x_det[i], 2) + std::pow(par[1]-y_det[i], 2));
    double sm = 0.78-0.21*par[2];

    //if(r < 0.1) r = 0.1;
    return (1/rm/rm/sm/sm)*std::tgamma(4.5-par[2])/(2*pi*tgamma(par[2])*tgamma(4.5-2*par[2]))           // Lagutin
           *(std::pow(((r+delt)/rm/sm),(par[2]-2.)))*std::pow((1+(r+delt)/rm/sm),(par[2]-4.5));
//    return (1/rm/rm/sm/sm)*std::tgamma(4.5-par[2])/(2*pi*tgamma(par[2])*tgamma(4.5-2*par[2]))
//           *(std::pow(((r+delt)/rm/sm/sm),(par[2]-2.)))*std::pow((1+(r+delt)/rm/sm/sm),(par[2]-4.5));
}

double Metropolis::Lagutin(int &i){
    double R_ms = 153.;
    double r = std::sqrt(std::pow(par[0]-x_det[i], 2) + std::pow(par[1]-y_det[i], 2));
    if(r < 0.1) r = 0.1;
    return (0.28/R_ms/R_ms)*std::pow((r/R_ms),(-1.2))*std::pow((1+r/R_ms),(-3.33))*std::pow((1.+std::pow((1./10./R_ms),2)),(-0.6));
}

double Metropolis::Lagutin_mod(int &i){
    int rm = 80.;
    double m = (0.78-0.21*par[2]);
    double r = std::sqrt(std::pow(par[0]-x_det[i], 2) + std::pow(par[1]-y_det[i], 2));
    r/=m;
    //if(r < 0.1) r = 0.1;
    return std::pow(m,-2)*(1/rm/rm)*std::tgamma(4.5-par[2])/(2*pi*std::tgamma(par[2])*std::tgamma(4.5-par[2]))*
            std::pow((r/rm),(par[2]-2))*std::pow((1+r/rm),(par[2]-4.5));
}

double Metropolis::Uchakin(int &i){
    int rm = 80.;
    double a = (4/par[2])*std::exp(0.915*(par[2]-1));
    double b = 0.15*1/(1+par[2]);
    double r = std::sqrt(std::pow(par[0]-x_det[i], 2) + std::pow(par[1]-y_det[i], 2));
    double C1 = std::pow(a, par[2]/b)/(2*pi*(std::tgamma(par[2]/b)+(4*std::tgamma(par[2]+1)/b)/(par[2]*std::pow(a,1/b))));
    return C1*std::pow(r,(par[2]-2)/rm)*(1+4/par[2]*r/rm)*std::exp(-a*std::pow(r,b)/rm);
}

//def Kascade(Ne, r):
//    rm = 80
//    return (1/rm/rm)*gamma(3.6-s)/(2*np.pi*gamma(s-1.5+2)*gamma(1.5+3.6-2*s-2))*(r/rm)**(s-2)*(1+r/rm)**(s-4.5)

double Metropolis::Kascade(int &i){
    double rm = 80.0;
    double delt = 0.1;
    double r = std::sqrt(std::pow(par[0]-x_det[i], 2) + std::pow(par[1]-y_det[i], 2));
    //if(r < 0.1) r = 0.1;
//    return (1/rm/rm)*std::tgamma(4.5-par[2])/(2*pi*tgamma(par[2])*tgamma(4.5-par[2]))
//           *(std::pow((r/rm),(par[2]-2.)))*std::pow((1+r/rm),(par[2]-4.5));
    return (1/rm/rm)*std::tgamma(3.6-par[2])/(2*pi*tgamma(par[2]+0.5)*tgamma(3.1-2*par[2]))
           *(std::pow((r/rm),(par[2]-2.)))*std::pow((1+r/rm),(par[2]-4.5));
}

double Metropolis::Greisen(int &i){
    double rm = 80;
    double r = std::sqrt(std::pow(par[0]-x_det[i], 2) + std::pow(par[1]-y_det[i], 2));
    if(r < 0.1) r = 0.1;
    return 0.4/rm/rm*std::pow((rm/r),(0.75))*std::pow((rm/(r+rm)),(3.25))*(1+(r/11.4/rm));
}

// Tibet attitude
double Metropolis::Tien_Shan(int &i){
    double rm = 136.;
    double r = std::sqrt(std::pow(par[0]-x_det[i], 2) + std::pow(par[1]-y_det[i], 2));
    return (1/rm/rm)*0.366*par[2]*par[2]*std::pow((2.07-par[2]),1.25)*std::pow((r/rm),par[2]-2)*std::pow((1+r/rm),par[2]-4.5);
}

double Metropolis::nkg_t(int &i){
    double rm = 136.;
    double delt = 0.1;
    double r = std::sqrt(std::pow(par[0]-x_det[i], 2) + std::pow(par[1]-y_det[i], 2)) /*+ 0.1*/;
    return (1/rm/rm)*std::tgamma(4.5-par[2])/(2*pi*tgamma(par[2])*tgamma(4.5-2*par[2]))                   // standard NKG
           *(std::pow((r/rm),(par[2]-2.)))*std::pow((1+r/rm),(par[2]-4.5));
}

double Metropolis::LHAASO(int &i){
    double rm = 136.;
    double r = std::sqrt(std::pow(par[0]-x_det[i], 2) + std::pow(par[1]-y_det[i], 2)) /*+ 0.1*/;
    return (1/rm/rm)*std::tgamma(4.5-par[2])/(2*pi*tgamma(par[2]-0.5)*tgamma(5-2*par[2]))
           *(std::pow((r/rm),(par[2]-2.5)))*std::pow((1+r/rm),(par[2]-4.5));
}


// ------------------------------

double Metropolis::L(){
    double F{}, N_i{};
    for(int i = 0; i < Num_det; ++i){
        //double ang;
        //if(std::cos(angl[0]*pi/180) < 40 and std::cos(angl[0]*pi/180)>=0) ang = std::cos(angl[0]*pi/180) ;
        //else  ang = 0;
        //if(std::sqrt(std::pow(par[0]-x_det[i], 2) + std::pow(par[1]-y_det[i], 2))>=5){
        N_i = par[3]*nkg(i)*S;
        //N_i = par[3]*nkg_t(i)*S;
        //N_i = par[3]*LHAASO(i)*S;
        //N_i = par[3]*nkg(i)*S*ang;
        //N_i = par[3]*Cors_nkg(i)*S;
        //if(x[i]< 15) F += x[i]/**ang*/* std::log(N_i) - N_i; //25
        //F += std::log(1/std::sqrt(2*pi*N_i)) - std::pow((x[i]/**ang*/-N_i)/(2*std::sqrt(N_i)),2);
        F += std::log(1/std::sqrt(2*pi*1)) - std::pow((x[i]/**ang*/-N_i)/(2*std::sqrt(1)),2);
        //else F += std::log(1/std::sqrt(2*pi*N_i)) - std::pow((x[i]/**ang*/-N_i)/(2*std::sqrt(N_i)),2);
        }
    //}
    return F;
}

double Metropolis::L2(){
    double F{}, N_i{};
    for(int i = 0; i < Num_det; ++i){
        //double ang;
        //if(std::cos(angl[0]*pi/180) < 40 and std::cos(angl[0]*pi/180)>=0) ang = std::cos(angl[0]*pi/180) ;
        //else  ang = 0;
        //N_i = par[3]*nkg2(i)*S;
        N_i = par[3]*nkg(i)*S;
        //N_i = par[3]*Cors_nkg(i)*S;
        /*if(x[i]< 15)*/ F += x[i]/**ang*/* std::log(N_i) - N_i; //25
        //else F += std::log(1/std::sqrt(2*pi*N_i)) - std::pow((x[i]/**ang*/-N_i)/(2*std::sqrt(N_i)),2);
    }
    return F;
}

void Metropolis::start_init(){      // Some zero initialization
        std::mt19937 mt(rd());
        double Sum_x{}, Sum_y{}, Q_all{}, Ne_0{};
        for(int i = 0; i < Num_det; ++i){
            Sum_x += x[i]*x_det[i];
            Sum_y += x[i]*y_det[i];
            Q_all += x[i];              // x - n (in MIP)
        }
        par.push_back(Sum_x/Q_all);
        par.push_back(Sum_y/Q_all);
        double s = 1.0;
        par.push_back(s);
        Sum_x = 0;
        for(int i = 0; i < Num_det; ++i)
            Sum_x+=S*nkg(i);
            //Sum_x+=S*LHAASO(i);
            //Sum_x+=S*Cors_nkg(i);
        Ne_0 = Q_all/Sum_x;
        par.push_back(Ne_0);
}

double Metropolis::ne_calc(){
    double Sum_x{},Q_all{};
    for(int i = 0; i < Num_det; ++i) Q_all += x[i];
    for(int i = 0; i < Num_det; ++i)
        Sum_x+=S*nkg(i);
        //Sum_x+=S*LHAASO(i);
        //Sum_x+=S*Cors_nkg(i);
    return Q_all/Sum_x;
}

double Metropolis::calc_del_r(){
    return std::sqrt(std::pow(par[0] - y[0],2)+std::pow(par[1] - y[1],2));
}

std::ofstream iter_("C:\\Qt\\progects\\MC_EAS\\iter.txt");
std::mutex g_lock;
void Metropolis::find_min(int N){
    int j{};
    int flag = 0;
    int flag2{};
    std::mt19937 mt(rd());
    int num_par{};
    std::uniform_real_distribution<double> distribution{0, 1};
    std::uniform_int_distribution<> it(0, 3);
    double gamma{};
    double step{};
    par_prom = par;
    double loss{}, loss2{};
    while (j < N){
        loss = L();
        num_par = it(mt);
        gamma = distribution(mt);
        step = (-1+2*gamma)*eps[num_par];
        par[num_par] +=step;
        loss2 = L();
        if (loss2 <= loss) {
            par[num_par] -=step;
        }
        if(std::abs(loss2-loss) < 1e-9 and flag2!=0) break;
        if(/*std::abs(loss2-loss) < 1e-5 and*/ j > N/2 and flag!=1){
            for(int i = 0; i <= 3; ++i) eps[i]/=2;//eps[i]*=0.1;
            flag = 1;
        }
        if(/*std::abs(loss2-loss) < 1e-7 and */ j > 2*N/3 and  flag2!=1 and flag!=0){
            for(int i = 0; i < 2; ++i) eps[i]/=2;
            flag2 = 1;
        }
        j+=1;
    }
}

void Metropolis::find_minf(int N){
    int j{};

    int flag = 0;
    int flag2{};
    std::mt19937 mt(rd());
    int num_par{};
    std::uniform_real_distribution<double> distribution{0, 1};
    std::uniform_int_distribution<> it(0, 2);
    if(eps[0] < 1.){
        for(int i = 0; i < 2; ++i) eps[i]*=10;
    }
    double gamma{};
    double step{};
    par_prom = par;
    double loss{}, loss2{};
    while (j < N){
        loss = L2();
        num_par = it(mt);
        //if(num_par==2) num_par=it(mt);

        gamma = distribution(mt);
        step = (-1+2*gamma)*eps[num_par];
        par[num_par] +=step;
        loss2 = L2();
        if (loss2 <= loss) par[num_par] -=step;
        if(std::abs(loss2-loss) < 1e-5 and flag!=1){
            for(int i = 0; i < 2; ++i) eps[i]*=0.1;
            flag = 1;
        }
        //if(j%(N/10)==0) std::cout << j << ' ' << par[0] << ' ' << par[1] << ' ' << par[2] << ' ' << loss << '\n';
        //if(loss <= 1e-3) break;
        j+=1;
    }
}

std::ofstream f1("C:\\Qt\\progects\\markov\\tot.txt");
void Metropolis::arr_dir(){
    int count{}, u{}, min{};
    double tetha{}, phi{};
    double alpha{}, betta{}, C{};
    int *tim = &t[0];
    //double *ni = &x[0];
    double *ni = new double[Num_det];
    //std::vector<double> ni;

    std::copy(std::begin(x), std::end(x), ni);
    Matrix M(3,3), Omega(3,1), X(3,1);
    Matrix x(Num_det, 3), x_T(Num_det, 3), y(Num_det, 1), W(Num_det, Num_det), z(1,1), _B(3,1), B(3,1);
    double tot{};
    int k{};
    double max_ni{};

    // Normalization
    for (int i = 0; i < Num_det; ++i) {
        if(tim[i]!=-1) {
            min = tim[i];
            break;
        }
    }
    for (int i = 0; i < Num_det; ++i) {
        if(tim[i] < min && tim[i]!=-1) min = tim[i];
    }
    for (int i = 0; i < Num_det; ++i) {
        if(tim[i]!=-1) tim[i]-=min;
    }

    // delete outliers
    for (int j = 0; j < Num_det; ++j ) {
        if(tim[j] > 1000) tim[j] = -1;
        if(tim[j]>=0) count+=1;
    }

    Matrix yy(16,1), xx(16,1);
    if(count > threshold){
        for(int j = 0; j < Num_det; ++j){
            if(tim[j]!=-1){
//                ++u;
                  x.Add_to_el(j, 0, x_det[j]);
                  x.Add_to_el(j, 1, y_det[j]);
                  x.Add_to_el(j, 2, 1);
                  y.Add_to_el(j, 0, c*tim[j]);

                  xx.Add_to_el(j,0,x_det[j]);
                  yy.Add_to_el(j,0,c*tim[j]);

/*                M.Add_to_el(0,0,x_det[j]*x_det[j]);
                M.Add_to_el(0,1,x_det[j]*y_det[j]);
                M.Add_to_el(0,2,x_det[j]);
                M.Add_to_el(1,1,y_det[j]*y_det[j]);
                M.Add_to_el(1,2,y_det[j]);
                M.Add_to_el(1,0,x_det[j]*y_det[j]);
                M.Add_to_el(2,0,x_det[j]);
                M.Add_to_el(2,1,y_det[j]);
                M.Set_el(2,2,u);
                Omega.Add_to_el(0,0,c*tim[j]*x_det[j]);
                Omega.Add_to_el(1,0,c*tim[j]*y_det[j]);
                Omega.Add_to_el(2,0,c*tim[j]);*/
            }
        }
//        std::cout << "X" << '\n';
//        x.show_arr();
//        std::cout << "X" << '\n';
//        y.show_arr();
//        std::cout << "y" << '\n';

        x_T = x;
        x_T.T(x);
        W.unit();               // ones

//        for(int i = 0; i < Num_det; ++i){  //f(ri)
//            if(min_r >= 0.1) ri[i]/=min_r;
//            else ri[i]/=0.1;
//            W.Mult_on_el(i,i,ri[i]);
//        }


//        for(int i = 0; i < Num_det; ++i){         //f(ni)
//            if(ni[i] < 1.) ni[i] = 0;
//            else if (ni[i] > 1. and ni[i] <=5) ni[i] = (ni[i]-0.5)/4.5;
//            else ni[i] = 1;
//            W.Mult_on_el(i,i,ni[i]);
//        }
        M.mult_arr(x_T, W);
        M.mult_arr(M,x);
        M.reverse();
        Omega.mult_arr(x_T, W);
        Omega.mult_arr(Omega, y);
        X.mult_arr(M, Omega);   // tetha

        _B = X;
        _B.mult_arr(x, X);
        z.difference(y, _B);
        //tot = z.sum();
        //double sigma2 = tot*tot/(count - 3);

       // z.abs();
       // tot = z.max(); //z.sum();

        //f1 << tot << '\n';
        //f1 << _B.Get_m() << ' ' << _B.Get_n() << ' ' << y.Get_m() << ' ' << y.Get_n() << '\n';
        double sigma2{};
        for(int i = 0; i < 16; ++i){
            //f1 << _B.Get_el(i,0) << ' ' << y.Get_el(i,0) << ' ' << z.Get_el(i,0) << ' ' << _B.Get_el(i,0) - y.Get_el(i,0) <<'\n';
            sigma2 += std::pow(_B.Get_el(i,0) - y.Get_el(i,0),2);
            //f1 << z.Get_el(i,0) << '\n';
        }
        sigma2/=(count-3);

        // ------IRLS ----------

//        while(k<15){
//            _B = X;
//            B = X;
//            _B.mult_arr(x, X);
//            z.difference(y, _B);
//            z.abs();
//            z.pow(-1);
//            //z.show_arr();
//            W.diag(z);
//            // W.show_arr();
//            M.mult_arr(x_T, W);
//            M.mult_arr(M,x);
//            M.reverse();
//            Omega.mult_arr(x_T, W);
//            Omega.mult_arr(Omega, y);
//            X.mult_arr(M, Omega);
//            z.difference(X, B);
//            z.abs();
//            tot = z.sum();
//          //  std::cout << tot << '\n';
//            ++k;
//        }

        alpha = X.Get_el(0,0);
        betta = X.Get_el(1,0);
        C = sqrt(1-alpha*alpha-betta*betta);


//        if(count >= 3){
//        tot = std::sqrt((1/(count - 2))*(yy.Dispersion(yy)/xx.Dispersion(xx)-alpha*alpha));
//        f1 << alpha << ' ' << xx.Dispersion(xx) << '\n';
//        }
        xx.clear();
        yy.clear();

        tetha = acos(C)*180.0/M_PI;
        phi = atan2(betta, alpha)*180.0/M_PI;
        if(phi < 0) phi+=360.0;

        if(tetha>7){
            angl.insert(angl.end(), {tetha, phi}); //out << tetha << "," << phi << "," << count << std::endl;
//            z.difference(y, _B);
//            double ch = z.sum();
//            double zn{};
//            for(int i = 0; i < 16; ++i){
//                if(y.Get_el(i,0)!=-1){
//                    zn+=y.Get_el(i,0)-y.Mean(yy);
//                }
//            }
//            double R2 = 1.0 - ch*ch/zn/zn;
//            f1 << R2 << ' ' << ch << ' ' << zn << '\n';
            f1 << alpha << ' ' << betta << ' ' << C << ' ' << sigma2*M.Get_el(0,0) << ' ' <<  sigma2*M.Get_el(1,1) << ' ' << sigma2*M.Get_el(2,2) << ' ' << sigma2*M.Get_el(2,2)/C*100  /*<< ' ' << sigma2*/ << '\n';
        }
        else {
        angl.insert(angl.end(), {-1, -1}); //out << tetha << "," << -1 << "," << count << std::endl;
        f1 << -1 << '\n';
        }
        //std::cout << tetha << "," << phi << "," << count  << std::endl;
        count = k = 0;

        x.clear();
        X.clear();
        _B.clear();
        z.clear();
        x_T.clear();
        W.clear();
        y.clear();
        M.clear();
        Omega.clear();

        x.Del();
        X.Del();
        _B.Del();
        z.Del();
        x_T.Del();
        W.Del();
        y.Del();
        M.Del();
        Omega.Del();
    }
    else{
        angl.insert(angl.end(), {-1, -1}); //out << -1 << "," << -1 << "," << count << std::endl;
        count = 0;
        f1 << -1 << '\n';
    }
}

void Metropolis::print_params(){
    if(par.size()==4)  std::cout << "EAS Parameters: " << par[0] << ' ' << par[1] << ' '<< par[2] << ' '<< par[3] << '\n';
    if(angl.size()==2) std::cout << "Arrival Direction: " << angl[0] << ' ' << angl[1] << ' ' << "Real: " << y_t[0] << ' ' << y_t[1] << '\n';
}

double Metropolis::calc_psi(){
    if(angl[0]!=-1 and angl[1]!=-1){
    double t_v = angl[0]*pi/180., p_v = angl[1]*pi/180., t = y_t[0]*pi/180., p = y_t[1]*pi/180.;
    return std::acos(std::sin(t_v)*std::cos(p_v)*std::sin(t)*std::cos(p)+
    std::sin(t_v)*std::sin(p_v)*std::sin(t)*std::sin(p)+
    std::cos(t_v)*std::cos(t))*180./pi;
    }
    else return -1;
}
