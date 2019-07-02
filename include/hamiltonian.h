#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include <vector>
#include <algorithm>
#include <stdexcept>
#include <cmath>
#include <armadillo>

# define GAP 0.000001
# define Pi 3.141592653589793
# define E  2.718281828459
# define STH 0.6931471805599453
# define GR 0.618033989

typedef double (*hmpointer) (arma::uword, arma::uword, arma::uword, std::vector<double>);// vec(0) hopping, vec(1) characteristic amplitude of potential
typedef double (*lenpointer) (std::vector<double> hm_param); // function return localization length

double avco(arma::vec & state); // calculate the position expectation value given wave function
std::vector<int> get_no_pos(std::vector<double> hm_param);
std::vector<int> get_qp_pos(std::vector<double> hm_param);
std::vector<int> get_tridiag_pos(std::vector<double> hm_param);

typedef std::vector<int> (*ranpospointer) (std::vector<double>);

//double random_h (arma::uword i, arma::uword j, arma::uword size, std::vector<double> param);
double qp_h(arma::uword i, arma::uword j, arma::uword size, std::vector<double> param);
double two_band_h(arma::uword i, arma::uword j, arma::uword size, std::vector<double> param);
double tbath_h(arma::uword i, arma::uword j, arma::uword size, std::vector<double> param);
double tridiag_h(arma::uword i, arma::uword j, arma::uword size, std::vector<double> param);

double qp_len(std::vector<double> hm_param);// suitable for case with mixed random and qp potentials
double two_band_len(std::vector<double> hm_param);
double tridiag_len(std::vector<double> hm_param);

class Hamiltonian
{
public:
    std::vector<double> lambda;
    std::vector<double> position;
    int n;
    double length;
    Hamiltonian(int size);
    void state_init(hmpointer hm, std::vector<double> hm_param, lenpointer lenpt, ranpospointer random_pos);
};

#endif //HAMILTONIAN_H
