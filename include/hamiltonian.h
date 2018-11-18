#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include <vector>
#include <algorithm>
#include <armadillo>
#include <cmath>

# define GAP 0.000001
# define Pi 3.141592653589793
# define E  2.718281828459
# define STH 0.6931471805599453
# define GR 0.618033989

typedef double (*hmpointer) (arma::uword, arma::uword, arma::uword, std::vector<double>);// vec(0) hopping, vec(1) characteristic amplitude of potential
typedef double (*lenpointer) (std::vector<double> hm_param); // function return localization length

double avco(arma::vec & state); // calculate the position expectation value given wave function
std::vector<int> get_no_pos();
std::vector<int> get_qp_pos();

double random_h (arma::uword i, arma::uword j, arma::uword size, std::vector<double> param);
double linear_h (arma::uword i, arma::uword j, arma::uword size, std::vector<double> param);
double qp_h(arma::uword i, arma::uword j, arma::uword size, std::vector<double> param);

double random_len(std::vector<double> hm_param);
double qp_len(std::vector<double> hm_param);// suitable for case with mixed random and qp potentials

class Hamiltonian
{
public:
    std::vector<double> lambda;
    std::vector<double> position;
    int n;
    double length;
    Hamiltonian(int size);
    void state_init(hmpointer hm, std::vector<double> hm_param, lenpointer lenpt, std::vector<int> random_pos = get_no_pos());
};

#endif //HAMILTONIAN_H
