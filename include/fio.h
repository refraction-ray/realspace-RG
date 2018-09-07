#ifndef FIO_H
#define FIO_H

#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <vector>
#include <unistd.h>
#include "hamiltonian.h"

#define INPUT 1024

std::vector<int> readin (const char*  filename, std::vector<double> & data);
void arg_parser(char* inputpath, char* outputpath, lenpointer lenf, int argc,char** argv );

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