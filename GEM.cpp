#include <cstdlib>
#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <fstream>
#include <numeric>
#include <Eigen/Dense>

using namespace std;

const double me = 0.51099895069;
const double mp = 938.27208943;
const double alpha = 7.2973525643e-3;
const double a0 = 0.529177210544*pow(10,5);
const double hbarc = 197.3263721;
const double c = 299792458;
double mu = me*mp/(me+mp);
const double PI   = 3.14159265358979323846;
const double Z = 1.0; 

double Knm(double rn, double rm) {
    return 3*sqrt(PI)*pow(rn*rm, 3) / (4*mu*hbarc*pow(rn*rn+rm*rm, 2.5));
}

double Vnm(double rn, double rm) {
    return -0.5*Z*alpha*rn*rn*rm*rm / ((rn*rn+rm*rm)*hbarc*hbarc); 
}

double Nnm(double rn, double rm) {
    return 0.25*sqrt(PI)*pow(rn*rm, 3) / (pow(rn*rn+rm*rm, 1.5)*hbarc*hbarc*hbarc);
}

int main() {
    double E;
    double r_int = a0*0.01;
    double rn;
    double rm;
    int N = 30;
    double a = 1.4;
    double best_E;
    double max_inv_det = -1.0;
    double inv_det = -1.0;

    vector<double> rn_list(N);
    vector<double> rm_list(N);
    vector<double> E_list(N*20000);
    vector<double> det_list(N*20000);

    Eigen::MatrixXd Mat_A(N, N);
    Eigen::MatrixXd Mat_K(N, N);
    Eigen::MatrixXd Mat_V(N, N);
    Eigen::MatrixXd Mat_N(N, N);
    Eigen::MatrixXd N_inv(N, N);
    
    double ii = 0;
    for (double E=-13.7; E < -13.5; E+=0.0001){
        double Etmp =  E * 1e-6;
        for (int n = 0; n < N; n++) {
            for (int m=0; m < N; m++){
                rn=pow(a,n)*r_int;
                //cout << "rn: " << rn << endl;
                rm=pow(a,m) * r_int;
                Mat_K(n, m) = Knm(rn, rm);
                Mat_V(n, m) = Vnm(rn, rm);
                Mat_N(n, m) = Nnm(rn, rm);
                Mat_A(n, m) = Knm(rn, rm) + Vnm(rn, rm) - Nnm(rn, rm)*Etmp;
            }
        }
        E_list[ii]=Etmp;
        double det = Mat_A.determinant();
        det_list[ii] = 1 / abs(det);

        inv_det = 1 / abs(det);
        

        if (inv_det > max_inv_det) {
            max_inv_det = inv_det;
            best_E = Etmp;
        }
        // cout << "det: " << det << endl;
        ii++;
    }

    // cout << "Mat_K: " << Mat_K << endl;
    // cout << "Mat_V: " << Mat_V << endl;
    // cout << "Mat_N: " << Mat_N << endl;

    cout << "best_E: " << best_E << endl;
    cout << "E_exact: " << -mu*alpha*alpha/2 << endl;

    std::ofstream ofile("GEM_energy.dat");
    ofile << "E" << " " << "1/det" << endl;
    ofile << std::setprecision(15);
    for (int i = 0; i < N*20; i++) {
        ofile << E_list[i] << " " << det_list[i] << endl;
    }

    double E_approx_MeV = -13.6057 * 1e-6;
    //best_E = -13.59872 * 1e-6;

    N_inv = Mat_N.inverse();

    Eigen::MatrixXd Mat_H_eff = N_inv * (Mat_K + Mat_V);
    Eigen::MatrixXd Mat_H = Mat_K + Mat_V;

    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(N, N);
    Eigen::MatrixXd Mat_A_shift = Mat_H_eff - best_E * I;

    Eigen::MatrixXd A_inv = Mat_A_shift.inverse();

    Eigen::VectorXd c_vec = Eigen::VectorXd::Random(N);
    Eigen::VectorXd c_vec_prev = c_vec;

    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(Mat_H, Mat_N);
    Eigen::VectorXd c_vec2 = eigensolver.eigenvectors().col(0).real();
    //Eigen::VectorXd c_vec2 = eigensolver.eigenvectors();

    int max_iter = 100;
    for(int iter=0; iter<max_iter; ++iter) {
        c_vec = A_inv * c_vec;

        double norm = c_vec.transpose() * Mat_N * c_vec;
        c_vec = c_vec / sqrt(abs(norm));

        double diff = (c_vec - c_vec_prev).norm();
        if(diff < 1e-10) {
            break;
        }
        c_vec_prev = c_vec;
    }

    for(int iter=0; iter<max_iter; ++iter) {
        c_vec2 = A_inv * c_vec2;

        double norm = c_vec2.transpose() * Mat_N * c_vec2;
        c_vec2 = c_vec2 / sqrt(abs(norm));
        cout << "norm: " << norm << endl;

        double diff = (c_vec2 - c_vec_prev).norm();
        if(diff < 1e-10) {
            break;
        }
        c_vec_prev = c_vec2;
    }

    cout << "c: " << c_vec(0) << endl;

    std::ofstream wf_file("wavefunction.dat");
    if (!wf_file) {
        cerr << "Error: Cannot open output file." << endl;
        return 1;
    }
    
    wf_file << "# r [fm]  Psi(r)  Exact(1s)" << endl;

    double r_max_plot = 8.0; 
    double dr = r_max_plot / 100.0;

    for (double r = 0.0; r <= r_max_plot; r += dr) {
        double psi_gem = 0.0;
        double psi_gem2 = 0.0;
        double ra = r*a0;
        for (int n = 0; n < N; n++) {
            
            double rn = pow(a, n) * r_int;
            psi_gem += c_vec(n) * exp( - (ra*ra) / (rn*rn) );
            psi_gem2 += c_vec2(n) * exp( - (ra*ra) / (rn*rn) );
        }

        double P_gem = psi_gem * psi_gem * ra * ra / pow(hbarc, 3);
        double P_gem2 = psi_gem2 * psi_gem2 * ra * ra / pow(hbarc, 3);

        //double psi_exact_shape = exp(-r);
        //double P_exact = ra*ra * exp(-ra) / a0;
        double P_exact = 4.0 * r * r * exp(-2.0 * r);
        //cout << "psi_gem: " << psi_gem << endl;
        wf_file << r << " " << P_gem*a0 << " " << P_gem2*a0 << " " << P_exact << endl;
    }

    wf_file.close();

}