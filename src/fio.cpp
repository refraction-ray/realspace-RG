#include "fio.h"


std::vector<int> readin (const char*  filename, std::vector<double>& data)
{
    using namespace std;
    ifstream fin;
    string line;
    double item;
    data.resize(INPUT);
    int i=0; int l=0;
    vector<int> result(2);
    fin.open(filename);
    if(fin.fail())
    {
        fin.close();

        return result;
    }
    while(getline(fin,line))
    {
        istringstream items(line);
        l++;
        while(items>>item)
        {
            data[i]=item;
            i++;
        }
    }
    fin.close();

    result[0] = i;// number of items
    result[1] = l;// number of lines
    return result;
}


std::vector<int> get_param_inpath()
{
    std::vector<int> result(6,0);
    result[1]=1;
    result[2]=1;
    result[5]=1;
    return result;

}

std::vector<int> get_param_inpath_pr() // used for two band model with all parameters in file path, only worked with -d
{
    std::vector<int> result(4,1);
    return result;

}

std::vector<int> get_param_inpath_tb()
{
    const int arr[] = {0,1,1,0,1,1,1};
    std::vector<int> result (arr, arr + sizeof(arr) / sizeof(arr[0]) );
    return result;
}


void arg_parser(char* inputpath, char* outputpath, bool& dist, char* distoutputpath, hmpointer& hmf,
        lenpointer& lenf, std::vector<int>& random_pos, std::vector<int>& param_path, int argc,char** argv )
{
    const char *optString = "i:o:d:rqpt";
    int opt=0;
    opt = getopt( argc, argv, optString );
    strcpy(inputpath, "input.txt");
    strcpy(outputpath,"output.txt");
    strcpy(distoutputpath,"");
    lenf = qp_len;
    hmf = qp_h;
    random_pos = get_qp_pos();
    param_path = get_param_inpath();

    while( opt != -1 ) {
        switch( opt ) {
            case 'i':
                strcpy(inputpath, optarg);
                break;

            case 'o':
                strcpy(outputpath, optarg);
                break;

            case 'd':
                strcpy(distoutputpath, optarg);
                dist = true;
                break;

            case 'r':
                lenf = random_len;
                break;

            case 'q':
                lenf = qp_len;
                break;

            case 'p':
                random_pos = get_no_pos();
                hmf = two_band_h;
                lenf = two_band_len;
                param_path = get_param_inpath_pr();
                break;

            case 't':
                hmf = tbath_h;
                lenf = qp_len;
                param_path = get_param_inpath_tb();

            default:
                /* You won't actually get here. */
                break;
        }

        opt = getopt( argc, argv, optString );
    }
}