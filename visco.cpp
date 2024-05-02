#include <limits>
#include <map>
#include <list>
#include <vector>
#include <random>
#include <cmath>
#include <iostream>
#include <fstream>
#include <set>
#include <sstream>
using namespace std;

inline int positive_modulo(int i, int n) {
    return (i % n + n) % n;
}



class viscoelastic
{
    public:
        int Lx;
        int Ly;
        int N;

        double* f_th;
        double* f; //viscoelastic part
        double* g;
        default_random_engine rnd_eng;
        weibull_distribution<double> thresholds_distr;
        double k0,k1,k2;

        int* activeSites;
        int* activeSitesNew;
        int* avalancheProfile;
        bool* isActiveNew;

        int get_nn(int site, int nn_idx)
        {
            int r = site / Ly;
            int c = site % Ly;
            int nn_site = -1;
            if(nn_idx==0)
            {
                nn_site = r*Ly + positive_modulo(c-1, Ly);
            }
            else if(nn_idx==1)
            {
                nn_site = r*Ly + positive_modulo(c+1, Ly);
            }
            else if(nn_idx==2)
            {
                nn_site = positive_modulo(r-1, Ly)*Ly + c;
            }
            else
            {
                nn_site = positive_modulo(r+1, Ly)*Ly + c;
            }
            return nn_site;
        }

        bool consistency_check()
        {
            for(int i = 0; i < N; i++)
            {
                double delta = f_th[i] - f[i] - g[i];
                if(delta < 0)
                {
                    cout << "incoherence (instabilities): " << i << " " << delta << endl;
                    return false;
                }
            } 
            return true;
        }


        viscoelastic(int n_rows, int n_cols, double param_k0, double param_k1, double param_k2, int seed=0)
        {
            k0 = param_k0;
            k1 = param_k1;
            k2 = param_k2;
            rnd_eng = default_random_engine();
            rnd_eng.seed(seed);
            thresholds_distr = weibull_distribution<double>(2.0, 1.0);
            Lx = n_rows;
            Ly = n_cols;
            N = Lx*Ly;
            f_th = new double[N];
            f = new double[N];
            g = new double[N];
            activeSites = new int[N];
            activeSitesNew = new int[N];
            avalancheProfile = new int[N];
            isActiveNew = new bool[N];
            for(int i = 0; i < N; i++)
            {
                activeSites[i] = -1;
                activeSitesNew[i] = -1;
                avalancheProfile[i] = 0;
                isActiveNew[i] = false;
                f_th[i] = thresholds_distr(rnd_eng);
                f[i]  = 0.0;
                g[i] = 0.0;   
            }
        }

        // system size trigger
        
        
        void avalanche(int origin, int& S_tot)
        {
            int nActiveSites = 1;
            activeSites[0] = origin;
            S_tot = 0;
            int nSitesTouched = 0;
            while(nActiveSites > 0)
            {   
                S_tot  = S_tot + nActiveSites;
                for(int k = 0; k < nActiveSites; k++)
                {
                    int site = activeSites[k];
                    f[site] -= 4*k2;
                    g[site] -= (4*k1+k0);
                    f_th[site] = thresholds_distr(rnd_eng);
                    avalancheProfile[site] += 1;
                    for(int nn_idx = 0; nn_idx < 4; nn_idx++)
                    {
                        int nn_site = get_nn(site, nn_idx);
                        f[nn_site] += k2;
                        g[nn_site] += k1;
                    }
                }
                
                int nActiveSitesNew = 0;
                for(int k = 0; k < nActiveSites; k++)
                {
                    int site = activeSites[k];
                    double delta = f_th[site] - f[site] - g[site];
                    if(delta < 0)
                    {
                        if(isActiveNew[site] == false)
                        {
                            isActiveNew[site] = true;
                            activeSitesNew[nActiveSitesNew] = site;
                            nActiveSitesNew++;
                        }
                    }                    
                    for(int nn_idx = 0; nn_idx < 4; nn_idx++)
                    {
                        int nn_site = get_nn(site, nn_idx);
                        double nn_delta = f_th[nn_site] - f[nn_site] - g[nn_site];
                        if(nn_delta < 0)
                        {
                            if(isActiveNew[nn_site] == false)
                            {
                                isActiveNew[nn_site] = true;
                                activeSitesNew[nActiveSitesNew] = nn_site;
                                nActiveSitesNew++;
                            }
                        }
                    }
                
                    activeSites[k] = -1; //we passed this site, reset
                }
                //cleanup
                for(int k = 0; k < nActiveSitesNew; k++)
                {
                    int site_to_reset = activeSitesNew[k];
                    isActiveNew[site_to_reset] = false;
                }
                swap(activeSites, activeSitesNew);
                nActiveSites = nActiveSitesNew;      
            }
            //consistency check
            //for(int i = 0; i < N; i++)
            //{
            //    if(f_th[i] < f[i] + g[i])
            //    {
            //        cout << "site still unstable" << i << " " << (f_th[i] - f[i] - g[i]) << endl;
            //    }
            //} 
        }
        
};

int main()
{
    int L = 512;
    double k0 = 0.01;
    double k1 = 0.4;
    double k2 = 0.3;
    int seed = 1;
    viscoelastic model = viscoelastic(L, L, k0, k1, k2,seed);
    
    stringstream filename;
    filename << "visco/recording_L=" << L << "_";
    filename << "k0=" << setprecision(5) << k0 << "_";
    filename << "k1=" << setprecision(5) << k1 << "_";
    filename << "k2=" << setprecision(5) << k2;
    ofstream snapshots;
    snapshots.open(filename.str() + ".dat", ios::out | ios::binary);
    ofstream recording;
    recording.open(filename.str() + ".txt");
    int sequence_counter = -1;
    int event_counter = 0;
    int min_num_events = 4000000;
    int end_of_avalanche_flag = -1;
    while(true)
    {
        int origin = -1;
        double trigger_amount = 0;
        bool is_relaxing;
        model.trigger(origin, is_relaxing, trigger_amount);
        //cout << "Origin " << origin << endl;
        int S_tot = 0;

       

        if(is_relaxing == false and event_counter >= min_num_events) //stop the simualtion when we reached a given amount of events AND the last sequence ended
        {
            break;
        }

        recording << (origin / L) << " " << (origin % L) << " ";
        recording << is_relaxing << " ";

        if(is_relaxing)
        {
            recording << std::setprecision(7) << (-log(trigger_amount)) << " ";
        }
        else
        {
            //sequence just ended, starting a new one
            sequence_counter += 1;
            recording << trigger_amount << " ";
        }

        



        //map<int, int> avalanche_history;
        //model.propagate_avalanche(origin, S_tot, avalanche_history);
        model.avalanche(origin, S_tot);


        for(int i =0; i < model.N; i++)
        {
            if(model.avalancheProfile[i]>0)
            {
                snapshots.write((char*)&i, sizeof(int));
                snapshots.write((char*)&model.avalancheProfile[i], sizeof(int));
            }
        }
        snapshots.write((char*)&end_of_avalanche_flag, sizeof(int));
        //snapshots.write((char*)(-1), sizeof(int));

        recording << S_tot << endl;
        //cout << "S: " << S_tot << endl;

        //if(model.consistency_check() == false)
        //{
        //    cout << trigger_amount << endl;
        //    cout << "Inconsistent: " << t << endl;
        //}
        //cout << "..." << endl;

        //avalanche_history.clear();
        //for(auto amount : avalanche_history)
        //{
        //    cout << (amount.first / model.Ly) << " " << (amount.first % model.Ly) << " "<< amount.second << endl;
        //}

        event_counter += 1;
    }
    recording.close();
    snapshots.close();
}