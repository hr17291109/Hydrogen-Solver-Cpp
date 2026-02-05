#include <cstdlib>
#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <fstream>
#include <numeric>

using namespace std;

const double me = 0.51099895069;
const double mp = 938.27208943;
const double alpha = 7.2973525643e-3;
const double a0 = 0.529177210544*pow(10,5);
const double hbarc = 197.3263721;
const double c = 299792458;
double mu = me*mp/(me+mp);

double k2(double r, double l, double E);
void chi(vector<double> &r_list, vector<double> &chi_list, int i, double l, double E, double h);

double k2(double r, double l, double E){
    double k;
    double ra = r*a0;
    if (r == 0.0){
        k = 0.0;
    } else {
        k = l*(l+1)/(ra*ra) - 2*mu*alpha/(hbarc*ra) - 2*mu*E/(hbarc*hbarc);
    }
    
    return k;
}

void chi(vector<double> &r_list, vector<double> &chi_list, int i, double l, double E, double h){
    double c1;
    double h1 = h*a0;
    c1 = ((2.0+(5.0/6.0)*k2(r_list[i-1], l, E)*h1*h1)*chi_list[i-1] - (1.0 - (1.0/12.0)*k2(r_list[i-2], l, E)*h1*h1)*chi_list[i-2]) / (1.0 - (1.0/12.0)*k2(r_list[i], l, E)*h1*h1);
    //c1 = ((2.0+(5.0/6.0)*k2(r_list[i-1], l, E))*chi_list[i-1] - (1.0 - (1.0/12.0)*k2(r_list[i-2], l, E))*chi_list[i-2]) / (1.0 - (1.0/12.0)*k2(r_list[i], l, E));
    chi_list[i] = c1;
    //cout << "k: " << k2(r_list[i], l, E) << endl;
}


int main(){
    double E = -13.5959e-6;
    // double E = -3.4e-6;
    double h = 0.01;
    int l = 0;
    double r;
    double r_int = 0.0;
    double r_max = 8;
    //double r_max2 = 20;
    double chi_int = 0.0;
    double dxdr_int = 0.0;
    double dxdr_tmp = 0.0;

    int N;
    N = (r_max-r_int)/h;
    
    vector<double> r_list(N);
    vector<double> chi_list(N);
    vector<double> dxdr_list(N);
    vector<double> E_list(100000);
    vector<double> chi_last_list(100000);
    vector<double> r_best(N);
    vector<double> chi_best(N);

    // cout << "k: " << k2(r_list[0], l, E) << endl;
    // cout << "k: " << k2(r_list[1], l, E) << endl;

    int ii = 0;
    double Ebest = 0;
    for (double E=-13.6; E < -13.5; E+=0.00001) {
        r_list[0] = r_int;
        r_list[1] = r_int + h;
        chi_list[0] = chi_int;
        chi_list[1] = pow(r_list[1]*a0, l+1);
        dxdr_list[0] = dxdr_int;
        dxdr_list[1] = dxdr_list[0] + k2(r_list[0], l, E)*chi_list[0]*h;
        r = r_list[1];
        double Etmp =  E * 1e-6;
        E_list[ii] = Etmp;
        for (int i = 2; i < N; i++) {
            r = r + h;
            r_list[i] = r;
            chi(r_list, chi_list, i, l, Etmp, h);
            dxdr_list[i] = dxdr_list[i-1] + k2(r_list[i-1], l, E)*chi_list[i-1]*h;
        }
        chi_last_list[ii] = 1 / (chi_list[N-1]*chi_list[N-1]);
        //cout << chi_list[N-1] << endl;
        if ((ii == 0) || (1 / (chi_list[N-1]*chi_list[N-1]) > 1 / (chi_best[N-1]*chi_best[N-1]))){
            r_best = r_list;
            chi_best = chi_list;
            Ebest = E;
            //cout << Ebest << endl;
        }
        ii++;
    }
    cout << Ebest << endl;
    cout << mu << endl;
    cout << -alpha*alpha*mu*0.5 << endl;

    //vector<double> a = {1,2,3};
    //double ipcheck = std::inner_product(a.begin(), a.end(), a.begin(), 0.0);
    //cout << ipcheck << endl;

    double ip_chi1 = std::inner_product(chi_best.begin(), chi_best.end(), chi_best.begin(), 0.0) * h;
    cout << ip_chi1 << endl;
    
    double norm_check = 0;

    // cout << chi_best[1] << endl;
    // cout << chi_best[2] << endl;
    // cout << chi_best[N-1] << endl;

    for (int i = 0; i < N; i++){
        chi_best[i] /= sqrt(ip_chi1);
        norm_check += chi_best[i]*chi_best[i]*h;
    }

    cout << norm_check << endl;

    vector<double> r_list2(N);
    vector<double> chi_list2(N);
    vector<double> dxdr_list2(N);
    E = Ebest;

    r_int = r_max;
    r_list2[0] = r_int;
    r_list2[1] = r_int - h;
    chi_list2[0] = chi_int;
    chi_list2[1] = exp(-sqrt(-2*mu*E/(hbarc*hbarc))*r_list2[1]*a0);
    dxdr_list2[0] = dxdr_int;
    dxdr_list2[1] = dxdr_list2[0] + k2(r_list2[0], l, E)*chi_list2[0]*h;

    r = r_list2[1];

    for (int i = 2; i < N; i++) {
        r = r - h;
        r_list2[i] = r;
        chi(r_list2, chi_list2, i, l, E, h);
        dxdr_list2[i] = dxdr_list2[i-1] + k2(r_list2[i-1], l, E)*chi_list2[i-1]*h;
    }

    double ip_chi2 = std::inner_product(chi_list2.begin(), chi_list2.end(), chi_list2.begin(), 0.0) * h;
    cout << ip_chi2 << endl;

    for (int i = 0; i < N; i++){
        chi_list2[i] /= sqrt(ip_chi2);
    }

    for (int i = 0; i < N; i++){
        double c = dxdr_list[i]/chi_list[i] - dxdr_list[i-1]/chi_list[i-1];
        if(c < 10e-7){
            //cout << i << endl;
        }
    }

    //cout << chi_list[2] << endl;

    std::ofstream ofile("Numerov001_l0.dat");
    std::ofstream ofile2("Numerov001_l0_2.dat");
    std::ofstream ofile1_5("Numerov001_l0_1.dat");
    std::ofstream ofile3("Esearch_l0.dat");
    ofile << "r" << " " << "chi" << endl;
    ofile2 << "r" << " " << "chi" << endl;
    ofile2 << "E" << " " << "1/chi2" << endl;
    ofile << std::setprecision(15);
    ofile1_5 << std::setprecision(15);
    ofile2 << std::setprecision(15);
    ofile3 << std::setprecision(15);
    for (int i = 0; i < N; i++){
        ofile << r_best[i] << " " << chi_best[i]*chi_best[i] << " " << dxdr_list[i] << endl;
        ofile2 << r_list2[i] << " " << chi_list2[i] << " " << dxdr_list[i] << endl;
        if (i < N/2) {
            ofile1_5 << r_best[i] << " " << chi_best[i]*chi_best[i] << " " << dxdr_list[i] << endl;
        } else {
            ofile1_5 << r_best[i] << " " << chi_best[i]*chi_best[i] << " " << dxdr_list[i] << endl;                    
        }

    }
    for (int i = 0; i < 100000; i++){
        ofile3 << E_list[i] << " " << chi_last_list[i] << endl;
    }
    ofile.close();
    ofile2.close();
    ofile3.close();

    return 0;

}
