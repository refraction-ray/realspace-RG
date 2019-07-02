#ifndef FIO_H
#define FIO_H

#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <vector>
#include <unistd.h>
#include "hamiltonian.h"

#define INPUT 4096

std::vector<int> get_param_inpath();
std::vector<int> get_param_inpath_pr();
std::vector<int> get_param_inpath_tb();


std::vector<int> readin (const char*  filename, std::vector<double> & data);
void arg_parser(char* inputpath, char* outputpath, bool& dist, char* distoutputpath,
        hmpointer& hmf, lenpointer& lenf, ranpospointer& random_pos, std::vector<int>& param_path,
        int argc,char** argv );

template <class T>
void print_1d(std::vector<T> &aa)
{
    int bb=(int)(aa.size());
    for(int i=0;i<bb;i++)
    {
        std::cout<<aa[i]<<" ";
    }
    std::cout<<"\n";
};

template <class T>
void print_2d(std::vector< std::vector<T> > &aa)
{
    int bb=(int)(aa.size());
    for(int i=0;i<bb;i++)
    {
        print_1d(aa[i]);
    }
    std::cout<<"\n";
};

template <class T>
void fprint_1d(std::vector<T> &aa, std::ofstream & fout)
{

    int bb = (int)(aa.size());
    for (int i=0;i<bb;i++)
    {
        fout<<aa[i]<<" ";
    }
    fout<<"\n";

}

template <class T>
void fprint_2d(std::vector< std::vector<T> > &aa, char* path)
{
    std::ofstream fout;
    fout.open(path);
    int bb = (int)(aa.size());
    for (int i=0;i<bb;i++)
    {
        fprint_1d(aa[i], fout);
    }
    fout.close();
}
#endif //FIO_H
