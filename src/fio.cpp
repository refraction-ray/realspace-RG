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


void arg_parser(char* inputpath, char* outputpath, bool& dist, char* distoutputpath, hmpointer& hmf,
        lenpointer& lenf, std::vector<int>& random_pos, int argc,char** argv )
{
    const char *optString = "i:o:d:rqp";
    int opt=0;
    opt = getopt( argc, argv, optString );
    strcpy(inputpath, "input.txt");
    strcpy(outputpath,"output.txt");
    strcpy(distoutputpath,"");
    lenf = qp_len;
    hmf = qp_h;
    random_pos = get_qp_pos();

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
                break;

            default:
                /* You won't actually get here. */
                break;
        }

        opt = getopt( argc, argv, optString );
    }
}