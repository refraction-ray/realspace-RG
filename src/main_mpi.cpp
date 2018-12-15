#include <mpi.h>
#include "rsrg.h"

int main(int argc, char** argv)
{
    using namespace std;

    char inputpath[INPUT], outputpath[INPUT], distoutputpath[INPUT];
    lenpointer lenf;
    bool dist = false;
    hmpointer hmf;
    std::vector<int> random_pos;
    arg_parser(inputpath,outputpath,dist,distoutputpath,hmf,lenf,random_pos,argc,argv);

    string dist_header(distoutputpath);
    vector<double> data(INPUT);
    vector<int> signal(2);

    signal=readin(inputpath, data);
    if(signal[0]==0)
    {
        return 1;
    }


    int rank;
    int world;
    vector< vector<double> > res;
    int si = signal[0]/signal[1]+2*get_measures().size();



    MPI_Status state;
    MPI_Init(NULL, NULL);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world);

    res.resize(world-1, vector<double>(si) );

    if(rank==0)
    {
            for(int i=1;i<world;i++)
            {
                for(int j=0;j<si;j++)
                {
                    double num;
                    MPI_Recv(&num, 1, MPI_DOUBLE, i, j, MPI_COMM_WORLD, &state);
                    res[i-1][j] = num;
                }


            }

            fprint_2d(res,outputpath);
        }

    else
    {
        if(signal[1]<=rank-1)
        {
            cout<<"no input data for "<<rank<<"\n";
            for(int i=0;i<si;i++)
            {
                double num=0;
                MPI_Send(&num,1,MPI_DOUBLE,0,i,MPI_COMM_WORLD);
            }

        }
        else
        {
            int r = rank-1;
            srand((unsigned)time(0)+r);
            vector<double> a(signal[0]/signal[1]-3);

            for(int j=0;j<signal[0]/signal[1]-3;j++)
            {
                a[j] = data[r*signal[0]/signal[1]+j];
            }
            vector<double> restemp;

            stringstream ss;
            ss<<"r"<<rank<<"_";
            dist_header = dist_header+ss.str();
            restemp = generate_measures(hmf, a, lenf, data[(r+1)*signal[0]/signal[1]-3], data[(r+1)*signal[0]/signal[1]-2],
                                  data[(r+1)*signal[0]/signal[1]-1], get_measures(), random_pos,
                                  dist, dist_header, get_param_inpath());

            a.push_back(data[(r+1)*signal[0]/signal[1]-3]);
            a.push_back(data[(r+1)*signal[0]/signal[1]-2]);
            a.push_back(data[(r+1)*signal[0]/signal[1]-1]);
            a.insert(a.end(),restemp.begin(),restemp.end());

            print_1d(a);

            for (int i=0;i<si;i++)
            {
                double num = a[i];
                MPI_Send(&num, 1, MPI_DOUBLE, 0, i, MPI_COMM_WORLD);
            }


        }
    }

    MPI_Finalize();
    return 0;
}
