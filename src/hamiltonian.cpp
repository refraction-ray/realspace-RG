#include "hamiltonian.h"

using namespace arma;


double avco (vec & state)
{
    double si=state.n_elem;
    vec lp=square(abs(state));
    rowvec rp=linspace(0,si-1,si).t();
    return as_scalar(rp*lp);
}

/*

//simple random hamiltonian and localization length is deprecated, please use the qp_h for mixed potential instead

double random_h (uword i, uword j, uword size, std::vector<double> param)
{
    double hop = param[0];
    double amp = param[1];
    double hop2;
    if (param.size()>2)
    {
         hop2 = param[2];
    }
    else
    {
        hop2 = 0;
    }
    double ele = 0;
    switch (std::abs((int)(i-j)))
    {
        case 0:
        {
            ele = amp*rand()/(double)(RAND_MAX);
            break;
        }

        case 1:
        {
          ele=hop;
          break;
        }

        case 2:
        {
            ele=hop2;
            break;
        }
    }
    return ele;
}

double random_len(std::vector<double> hm_param)
{
    return 2/log(1+hm_param[5]*hm_param[5]+GAP);
}

*/

double qp_h(uword i, uword j, uword size, std::vector<double> param) //hamiltonian with both random and quasi potential
{
    double hop = param[0];
    double amp = param[1];
    double wavevector = param[2];
    double phase = param[3];
    double hop2;
    if (param.size()>4)
    {
         hop2 = param[4];
    }
    else
    {
         hop2 = 0;
    }
    double ramp;
    if (param.size()>5)
    {
        ramp = param[5];
    }
    else
    {
        ramp=0;
    }

    double ele=0;
    switch (std::abs((int)(i-j)))
    {
        case 0:
        {
            ele=amp*cos(wavevector*(double)i+phase);
            if (ramp >= 0)
            {
                ele = ele + ramp*rand()/(double)(RAND_MAX);
            }
            else
            {
                ele = ele + ramp*(rand()/(double)(RAND_MAX)-0.5);
            }

            break;
        }

        case 1:
        {
            ele=hop;
            break;
        }

        case 2:
        {
            ele=hop2;
            break;
        }
    }
    return ele;

}

double qp_len(std::vector<double> hm_param)
{
    double ramp = 0;

    if (hm_param.size()>5)
    {
        ramp = hm_param[5];
    }
    else
    {
        ramp = 0;
    }

    double r_ll = 2/log(1+ramp*ramp+GAP);
    if (hm_param[1] <= 2+GAP)
    {
        return r_ll;
    }

    double qp_ll = 1/log(1+hm_param[1]/2+GAP);

    if (ramp == 0)
    {
        return qp_ll;
    }

    return qp_ll*r_ll/(qp_ll+r_ll);
}


double two_band_h (uword i, uword j, uword size, std::vector<double> param)
{
    if ( i != j)
    {
        return 0;
    }
    double band1 = param[0];
    double  band2 = param[1];
    double gap1 = param[2];
    double criteria = param[3];//eg. the probability to stay on the same band, 0.005
    static int lasttime = 0;
    double p_transfer = rand()/(double)(RAND_MAX);
    double p_pos = rand()/(double)(RAND_MAX);
    if (p_transfer>criteria)
    {
        if (lasttime == 0)
        {
            lasttime = 1;
        }
        else
        {
            lasttime = 0;
        }
    }

    double ele;

    if (lasttime == 1)
    {
        ele = band1+gap1+band2*p_pos;
    }
    else
    {
        ele = band1*p_pos;
    }

    return ele;
}

double two_band_len(std::vector<double> hm_param)
{
    double weff = hm_param[0]+hm_param[1]+hm_param[2];
    return 2/log(1+weff*weff+GAP);
}

double tbath_h(uword i, uword j, uword size, std::vector<double> param)
{
    // 0 NN-hopping 1 qp-potential-amplitude 2 qp-potential-wavevector 3 qp-potential-sample-window
    // 4 NNN-hopping 5 random-onsite-potential 6 bath-size
   int bathsize = (int) param[6];
   if ( i<bathsize || j<bathsize )
   {
        if (i==j){ return GAP*rand()/(double)(RAND_MAX);}
        else {return 0;}
   }
   else
   {
        return qp_h(i, j, size, param);
   }

}



std::vector<int> get_no_pos() // the model has no random sample parameter like phi in qp
{
    std::vector<int> pos(1);
    return pos;
}

std::vector<int> get_qp_pos()
{
    std::vector<int> pos(4,0);
    pos[3]=1;
    return pos;

}

void Hamiltonian::state_init(hmpointer hm, std::vector<double> hm_param, lenpointer lenpt,
        std::vector<int> random_pos )
{
    mat ha(n,n);

    for (int i=0;i<random_pos.size();i++)
    {
        if( random_pos[i] == 1)
        {
            hm_param[i] = hm_param[i]*rand()/(double)(RAND_MAX);
        }
        else if( random_pos[i] == -1)
        {
            hm_param[i] = hm_param[i]*(2*rand()/(double)(RAND_MAX)-1);
        }
    }

    for(uword i=0;i<n;i++)
    {
        for(uword j=0;j<n;j++)
        {
            ha(i,j)=hm(i,j,n,hm_param);
        }
    }

    vec spec,mi;
    mat wavf;
    eig_sym(spec, wavf, ha);

    rowvec ss(n);
    for(uword i=0;i<n;i++)
    {
        mi = wavf.col(i);
        ss(i) = avco(mi);
    }
    urowvec dd = stable_sort_index(ss).t();
    double lower_bound = min(spec);
    for(uword i=0;i<n;i++)
    {
        lambda[i] = ( spec(dd(i))-lower_bound+GAP*hm_param[1] );
        position[i] = ( ss(dd(i)) );
    }

    length = lenpt(hm_param);
}



Hamiltonian::Hamiltonian(int size)
{
    this->n=size;
    lambda.resize(n);
    position.resize(n);
};

