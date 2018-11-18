#ifndef RSRG_H
#define RSRG_H

#include <vector>
#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric> // for accumulate function compatibility
#include "fio.h"
#include "hamiltonian.h"
#include "set.h"


bool const LOG = false; // the switch for detailed debug logs



bool comparelts(const std::vector<int>& p1,const  std::vector<int>& p2);

double jack_knife(std::vector<double> & data, double tmean);

class Model // the model and data for RG
{
    public:

        int n;
        int ntotal;
        std::vector<int> cluster;
        std::vector<double> lambda;
        std::vector< std::vector<double> > deltae;
        std::vector< std::vector<double> > gamma;
        std::vector< std::vector<int> > sites;

        Model(int size);

        void rginit(Hamiltonian & h, double V);
        void rgstep();
        void fixedpoint();
        double entangle(int cut=-1);
        double maxlength(int benchmark=-1);
        double averagelength(int benchmark=-1);


};


typedef double (Model::*measurement) (int);


std::vector<measurement> get_measures();
std::vector<int> get_param_inpath();

std::vector<double> generate_measures(hmpointer hm, std::vector<double> hm_param, lenpointer lenpt,
                         int size, double V, int repeat, std::vector< measurement > measures,
                         std::vector<int> random_pos, bool dist, std::string, std::vector<int> param_path);
/*the main function of the package*/


#endif //RSRG_H
