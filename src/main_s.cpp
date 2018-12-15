#include "rsrg.h"



int main(int argc, char** argv)
{
    using namespace std;
    srand((unsigned)time(0));
    char inputpath[INPUT], outputpath[INPUT], distoutputpath[INPUT];
    bool dist = false;
    lenpointer lenf;
    hmpointer hmf;
    std::vector<int> random_pos;

    arg_parser(inputpath,outputpath,dist,distoutputpath,hmf,lenf,random_pos,argc,argv);

    string dist_header(distoutputpath);
    vector<double> data(INPUT);
    vector<int> signal(2);

    signal = readin(inputpath, data);
    if(signal[0]==0)
    {
        return 1;
    }
    vector< vector<double> > res;



    for(int i=0;i<signal[1];i++)
    {
        vector<double> a(signal[0]/signal[1]-3);
        for(int j=0;j<signal[0]/signal[1]-3;j++)
        {
            a[j] = data[i*signal[0]/signal[1]+j];

        }

        vector<double> restemp;

        restemp = generate_measures(hmf, a, lenf, data[(i+1)*signal[0]/signal[1]-3], data[(i+1)*signal[0]/signal[1]-2],
                data[(i+1)*signal[0]/signal[1]-1], get_measures(), random_pos, dist, dist_header, get_param_inpath());

        a.push_back(data[(i+1)*signal[0]/signal[1]-3]);
        a.push_back(data[(i+1)*signal[0]/signal[1]-2]);
        a.push_back(data[(i+1)*signal[0]/signal[1]-1]);
        a.insert(a.end(),restemp.begin(),restemp.end());
        res.push_back(a);
        print_1d(res[i]);

    }

    fprint_2d(res, outputpath);

    return 0;
}