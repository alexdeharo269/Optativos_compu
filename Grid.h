
#include <vector>
#include<iostream>
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <chrono>
#include <random>
using namespace std; 

class Grid{
    private:
        using index = vector<int>::size_type; // Define index as the size type of vector<int>
        vector<vector<int>> grid;
        index size;
        mt19937_64 generator;
        uniform_real_distribution<double> r_distribution;

        public:
            // fonstructor
            Grid(index n, unsigned seed) : size(n), generator(seed), r_distribution(0.0, 1.0)
        {
            grid.resize(n, vector<int>(n));
        }

        
        double magnetization(int domain,float init_mg) const
        // Magnetización del grid
        // Domains: completa -->3, abajo-->1, arriba-->2
        // Es const al final para poder acceder a ella dentro de otras funciones.
        {
            index plus =static_cast<index>( floor((init_mg)*static_cast<float>(size)));
            index l_lim; 
            index u_lim; 
            if (domain==0){l_lim=0;u_lim=plus;}else if(domain==1){l_lim=plus;u_lim=size;}else if(domain==3){l_lim=0;u_lim=size;}
            double sum = 0;
            for (index i = 0 ; i <size ; i++)
            {
                for (index j =l_lim; j <  u_lim  ; j++)
                {
                    sum += grid[static_cast<index>(i)][static_cast<index>(j)];
                }
            }
            return sum / static_cast<double>( size*(u_lim-l_lim));
        }
        
        struct domains
        {
            float ur;
            float ul;
            float dr;
            float dl;
        };
        
        struct region{
            double up;
            double down;
        };

        double initialize(const string &init_options,float init_mg)
        // Queremos iniciarla con magnetización inicial nula. Varios enfoques:
        // uqr, opposite sign up, pongo el signo contrario en el j+1 cuando inicialize un j;
        // ubd, until bound, magnetizo hasta 0 con un cierto error.
        {
            if (init_options == "ubd")
            {
                for (index i = 0; i < size; i++)
                {
                    grid[i][size - 1] = -1;
                    grid[i][0] = 1;
                    for (index j = 1; j < size - 1; j++)
                    {
                        double s;
                        s = r_distribution(generator);
                        if (s < 0.5)
                        {
                            grid[i][j] = 1;
                        }
                        else
                        {
                            grid[i][j] = -1;
                        }
                    }
                }
                return magnetization(3,init_mg);
                
            }
            else if (init_options == "uqr")
            {
                for (index i = 0; i < size; i++)
                {
                    grid[i][size - 1] = -1;
                    grid[i][0] = 1;
                    if(i%2==0)
                    {
                        for (index j = 1; j < size - 1; j++)
                        {
                            if (j % 2 == 0)
                            {
                                grid[i][j] = -1;
                            }
                            else
                            {
                                grid[i][j] = 1;
                            }
                        }
                    }
                    else{
                        for (index j = 1; j < size - 1; j++)
                        {
                            if (j % 2 == 0)
                            {
                                grid[i][j] = 1;
                            }
                            else
                            {
                                grid[i][j] = -1;
                            }
                        }
                    }
                }
                return 0.0f;
            }
            else if(init_options=="mgt"){
                //Vamos a poner dominios uniformes
                index plus_bound=static_cast<index>( floor(init_mg*static_cast<double> (size)));
                for (index i=0; i<size;i++){
                    for (index j=0; j<plus_bound; j++){
                        grid[i][j]=1;
                    }
                    for (index j = plus_bound; j < size; j++){
                        grid[i][j]=-1;
                    }
                }
            }
            return 0;
        }

        double energy(){ //Energy of the configuration
            double e=0;
            index i,j;
            int left, right;  //Para poder calcular los vecinos izquierdo y derecho teniendo en cuenta CC continuas en x


            j=0;
            e+=grid[0][1]+1;  
            for(i=1;i<size-1;i++){ 
                e+=grid[i][1]+2;
            }
            e+=grid[size-1][1]+1;

            for(j=1;j<size-1;j++){i=0;
                left=grid[size-1][j];
                e+=grid[i][j]*(grid[i+1][j]+left+grid[i][j-1]+grid[i][j+1]);    //Ahoramismo se esta saliendo para j
                for(i=1; i<size-1;i++){
                    left=grid[i-1][j]; right=grid[i+1][j];
                    e+=grid[i][j]*(left+right+grid[i][j-1]+grid[i][j+1]);
                }i=size-1;                                                          //Se podría comprobar que aquí i llega hasta size-1.
                right=grid[0][j];
                e+=grid[i][j]*(grid[i][j-1]+grid[i][j+1]+grid[i-1][j]+right);               
            }
            j=size-1;
            
            e-=grid[0][size-2]-1;
            for(i=1;i<size-1;i++){
                e-=grid[i][size-2]-2;
            }
            e-=grid[size-1][size-2]-1;

            e=-0.5*e;

            return e;
        }

        double deltaE(int x1, int y1, int x2, int y2,int vecino){
            double deltaE;
            index i1 = static_cast<index>(x1);
            index i2 = static_cast<index>(x2);
            index j1 = static_cast<index>(y1);
            index j2 = static_cast<index>(y2);

            //Condiciones de contorno continuas
            int l1, l2, r1, r2; //left1, left 2...
            if(i1==0){l1=grid[size-1][j1];r1=grid[1][j1];}
            else if(i1==size-1){r1=grid[0][j1];l1=grid[size-2][j1];}
            else{l1=grid[i1-1][j1];r1=grid[i1+1][j1];}

            if(i2==0){l2=grid[size-1][j2];r2=grid[1][j2];}
            else if(i2==size-1){r2=grid[0][j2];l2=grid[size-2][j2];}
            else{l2=grid[i2-1][j2];r2=grid[i2+1][j2];}

            if(vecino==1){
                deltaE=2*(grid[i1][j1]*(r1+grid[i1][j1+1]+grid[i1][j1-1])+grid[i2][j2]*(l2+grid[i2][j2+1]+grid[i2][j2-1]));
            }

            if(vecino==2){
                deltaE=2*(grid[i1][j1]*(l1+r1+grid[i1][j1-1])+grid[i2][j2]*(l2+r2+grid[i2][j2+1]));
            }
            if(vecino==3){
                deltaE=2*(grid[i1][j1]*(l1+grid[i1][j1+1]+grid[i1][j1-1])+grid[i2][j2]*(r2+grid[i2][j2+1]+grid[i2][j2-1]));
            }
            if(vecino==4){
                deltaE=2*(grid[i1][j1]*(l1+r1+grid[i1][j1+1])+grid[i2][j2]*(l2+r2+grid[i2][j2-1]));
            }
            deltaE=0.5*deltaE;


            return deltaE;
        }

        
        int get(int i, int j){
            return grid[static_cast<index>(i)][static_cast<index>(j)];
        }
        
        void set(int i, int j, int val){
            grid[static_cast<index>(i)][static_cast<index>(j)]=val;
        }

        region magnet_regions(float init_mgt){
            region hemisferio;
            hemisferio.up = magnetization(0, init_mgt);
            hemisferio.down = magnetization(1, init_mgt);
            return hemisferio;
        }
};
