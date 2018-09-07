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


void arg_parser(char* inputpath, char* outputpath, lenpointer lenf, int argc,char** argv )
{
    const char *optString = "i:o:rq";
    int opt=0;
    opt = getopt( argc, argv, optString );
    strcpy(inputpath, "input.txt");
    strcpy(outputpath,"output.txt");
    lenf = qp_len;

    while( opt != -1 ) {
        switch( opt ) {
            case 'i':
                strcpy(inputpath, optarg); /* true */
                break;

            case 'o':
                strcpy(outputpath, optarg);
                break;

            case 'r':
                lenf = random_len;
                break;

            case 'q':
                lenf = qp_len;
                break;

            default:
                /* You won't actually get here. */
                break;
        }

        opt = getopt( argc, argv, optString );
    }
}