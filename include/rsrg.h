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


std::vector<double> generate_measures(hmpointer hm, std::vector<double> hm_param, lenpointer lenpt,
                         int size, double V, int repeat, std::vector< measurement > measures,
                         ranpospointer random_pos, bool dist, std::string pathheader, std::vector<int> param_path);

/* the main function of the package
hm: the hamiltonian function pointer for single-particle hamiltonian matrix elements
hm_param: vector of Hamiltonian paramters
lenpt: function pointer for localization length function
size: system size
V: typical interaction value
repeat: numbers of configurations
measures: vector of function pointer for different observables
random_pos: vector of 0 and 1, when 1 occurs, it means in that position the paramter of Hamiltonian is random within window given by hm_param
dist: boolean for whether output detail observables for each configuration
pathheader: the path header of dist output file
param_path: vector indicating which position of hm_param should occur in the dist file name (indicate by 1 otherwise 0)
*/


#endif //RSRG_H
