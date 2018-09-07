#include "rsrg.h"


using namespace std;


Model::Model(int size)
{
    n=size;
    ntotal = size;
    cluster.resize(n);
    lambda.resize(n);
    deltae.resize(n,vector<double>(n));
    gamma.resize(n,vector<double>(n));
    sites.resize(n,vector<int>(n,-1));
};


void Model::rginit(Hamiltonian &hamiltonian, double V)
{
    lambda = hamiltonian.lambda;

    if(LOG == true)
    {
        cout<<"initial lambda status\n";
        print_1d(lambda);
    }

    for(int i=0;i<n;i++)
    {
        cluster[i]=1;
    };

    for (int j = 0; j <n; j++)
    {
        sites[j][0]=j;
        for (int k = 0; k < n; k++)
        {
            deltae[j][k] = fabs(lambda[j]-lambda[k]);
            gamma[j][k] = V*(pow(E,-fabs(hamiltonian.position[j]-hamiltonian.position[k])/hamiltonian.length));
                   // +V*(pow(E,-(n-fabs(hamiltonian.position[j]-hamiltonian.position[k]))/hamiltonian.length));

        }

    };

};


void Model::rgstep()
{

    Graph g(n);
    for(int i=0;i<n-1;i++)
    {
        for(int j=i+1;j<n;j++)
        {
            if (gamma[i][j]-deltae[i][j]>0)
            {
                g.addEdge(i, j);
            }
        };
    };
    vector<int>* gp=g.connectedComponents(g.V);
    vector<vector<int> > np;
    for (int i = 0; i < n; i++)
    {
        if( gp[i].size()==0 )
        {
            break;
        }
        else
        {
            vector<int> tmp(gp[i].size()+1);
            tmp[0]=(gp[i].size());
            for (int j=0;j<gp[i].size();j++)
            {
                tmp[j+1]=(gp[i][j]);
            }
            np.push_back(tmp);
        }
    }
    delete [] gp;



    int new_n = np.size();
    sort(np.begin(),np.end(),comparelts);
    vector< vector<int> > new_sites;
    new_sites.resize(new_n, vector<int>(ntotal,-1));

    if(LOG==true)
    {
        cout<<"np status\n";
        print_2d(np);
    }

    for(int j=0;j<new_n;j++)
    {
        int start=0;
        for(int k=0;k<np[j][0];k++)
        {
            for(int l=0;l < sites[np[j][k+1]].size();l++)
            {
                if (sites[np[j][k+1]][l]>=0)
                {

                    new_sites[j][start] = sites[np[j][k+1]][l];
                    start++;
                }
                else
                {
                    break;
                }


            }

        }
    };

    if(LOG==true)
    {
        cout<<"the cluster status\n";
        print_2d(new_sites);
    }

    vector<double> lam_lampart(new_n);
    for(int j=0;j<new_n;j++)
    {
        vector<double> lamsquare;

        for(int kk=1;kk<=np[j][0];kk++)
        {
            lamsquare.push_back(pow(lambda.at(np[j][kk]),2));
        };

        lam_lampart[j]=(accumulate(lamsquare.begin(),lamsquare.end(), 0.0));
    };

    vector<double> lam_gampart(new_n);
    for(int j=0;j<new_n;j++)
    {
        vector <double> gammatemp;
        if(np[j][0]==1)
        {
            gammatemp.push_back(0);
        }
        else
        {
            for(int kk=0;kk<np[j][0]-1;kk++)
            {
                for(int mm=kk+1;mm<np[j][0];mm++)
                {
                    gammatemp.push_back(pow(gamma[np[j][kk+1]][np[j][mm+1]],2));
                }
            };
        }
        lam_gampart[j]=(accumulate(gammatemp.begin(),gammatemp.end(), 0.0));
    };

    vector<double> new_lambda(new_n);

    for(int j=0;j<new_n;j++)
    {
        new_lambda[j] = (sqrt(      lam_lampart[j]+lam_gampart[j]      ));
    };

    vector<int> new_cluster(new_n);
    for(int j=0;j<new_n;j++)
    {
        vector<int> clutemp;
        for(int mm=0;mm<np[j][0];mm++)
        {
            clutemp.push_back(cluster.at(np[j][mm+1]));
        }
        new_cluster[j] = (accumulate(clutemp.begin(), clutemp.end(), 0));
    };

    vector<double> del(new_n);
    for(int mm=0;mm<new_n;mm++)
    {

        del[mm] = new_lambda[mm]/(pow(2,new_cluster[mm])-1);

    };

    double ran,tan;
    vector< vector<double> > new_deltae;
    new_deltae.resize(new_n,vector<double>(new_n));

    for(int mm=0;mm<new_n;mm++)
    {
        for(int nn=0;nn<new_n;nn++)
        {
            tan=max(del[mm]-new_lambda[nn],del[nn]-new_lambda[mm]);
            if(tan<0)
            {
                ran = del[nn]*del[mm]/(min(new_lambda[mm],new_lambda[nn]));
            }
            else
            {
                ran = max(tan,min(del[mm],del[nn]));
            }

            new_deltae[mm][nn] = ran;
        };
    };

    vector< vector<double> > new_gamma;
    new_gamma.resize(new_n,vector<double>(new_n));

    for(int mm=0;mm<new_n;mm++)
    {
        for(int nn=0;nn<new_n;nn++)
        {
            if(np[mm][0]==1&&np[nn][0]==1)
            {

                new_gamma[mm][nn]=0;

            }
            else
            {

                double gamax=0; int locx=0,locy=0;
                for(int rr=0;rr<np[mm][0];rr++)
                {
                    for(int qq=0;qq<np[nn][0];qq++)
                    {
                        double comp = gamma[np[mm][rr+1]][np[nn][qq+1]];
                        if(comp>gamax)
                        {
                            gamax=comp;
                            locx=rr+1;
                            locy=qq+1;
                        }

                    }
                }

                new_gamma[mm][nn] =
                        gamax*pow(E,-(new_cluster[mm]+new_cluster[nn]-cluster[np[mm][locx]]-cluster[np[nn][locy]])*STH/2);


            }
        }
    }

    n = new_n;
    cluster = new_cluster;
    lambda = new_lambda;
    deltae = new_deltae;
    gamma = new_gamma;
    sites = new_sites;

    if(LOG == true)
    {
        cout<<"new DeltaE matrix\n";
        print_2d(deltae);
        cout<<"new Gamma matrix\n";
        print_2d(gamma);
    }

};


void Model::fixedpoint ()
{
    int rem;
    while(1)
    {
        rem = n;
        this->rgstep();
//        cout<<n<<" ";
        if(n==rem)
        {
            break;
        }
    }
};

double Model::maxlength (int benchmark)
{
    double maxlen;
    maxlen=*max_element(cluster.begin(), cluster.end());
    return (maxlen+(double)benchmark)/(double)ntotal;
}

double Model::averagelength (int benchmark)
{
    double avlen;
    avlen=(double)ntotal/(double)n;
    return (avlen+(double)benchmark)/(double)ntotal;
}

double Model::entangle (int cut)
{
    vector<double> piece;

    if(cut == -1)
    {
        cut=ntotal/2;
    }

    double p=0,q=0,totals=0;
    for(int i=0;i<n;i++)
    {
        p=0;q=0;
        for(int j=0;j<sites[i].size();j++)
        {
            if (sites[i][j]<0)
            {
                break;
            }
            if(sites[i][j]<cut)
            {
                p++;
            }
            else
            {
                q++;
            }
        }
        piece.push_back(min(p,q)*STH);

    }
    totals=accumulate(piece.begin(),piece.end(),0.0);

    return totals/(double)cut/STH;
}

vector<measurement> get_measures()
{
    std::vector<measurement> default_measures(3);
    default_measures[0] = &Model::averagelength;
    default_measures[1] = &Model::maxlength;
    default_measures[2] = &Model::entangle;
    return default_measures;
}

vector<double> generate_measures(hmpointer hm, vector<double> hm_param, lenpointer lenpt,
        int size, double V, int repeat, vector< measurement > measures, vector<int> random_pos)
{
    Hamiltonian hamiltonian(size);
    vector< vector<double> > results;
    int n_measures = measures.size();
    results.resize(n_measures, vector<double>(repeat));
    vector<double> average_results(2*n_measures);
    for(int i=0;i<repeat;i++)
    {
        if (LOG == true)
        {
            cout<<"\n"<<"--------------\n"<<i<<" round start\n";
        }

        hamiltonian.state_init(hm,hm_param,lenpt,random_pos);
        Model model(size);
        model.rginit(hamiltonian, V);
        model.fixedpoint();

        for(int j=0;j<n_measures;j++)
        {
            results[j][i] = (model.*(measures[j]))(-1);
        }


    }

    for(int i =0;i<n_measures;i++)
    {
        average_results[2*i] = accumulate(results[i].begin(),results[i].end(),0.00)/(double)repeat;
        average_results[2*i+1] = jack_knife(results[i], average_results[2*i]);
    }

    return average_results;
}

bool comparelts (const vector<int>& p1,const  vector<int>& p2)
{
    int cut=min(p1.size(),p2.size());
    for(int i=0;i<cut;i++)
    {
        if(p1[i]!=p2[i])
        {
            return( p1[i]>p2[i]);
        }
    }
    return ( p1[0]>p2[0] );
}

double jack_knife(vector<double> & data, double tmean)
{

    int si = data.size();
    double sum = tmean*si;
    double var =0;
    double mean;
    for (int i=0;i<si;i++)
    {
        mean = (sum-data[i])/(double)(si-1);
        var = var+pow(mean-tmean,2.0);
    }
    var = sqrt(var*(si-1)/si);
    return var;
}