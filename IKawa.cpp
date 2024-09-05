
// g++ programa.cpp -o Pprograma nohup ./Program
/*DUDAS:

*/
#include <vector>
#include<iostream>
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <chrono>
#include <random>
#include "Grid.h"
using namespace std;
using index = vector<int>::size_type;

const int n = 16;            
const int t_max = 500;             
const int t_size=50;
float max_asymetry=0.01f;
float init_mgt=0.5f;

double p_val(int i1, int i2, int j1, int j2, float temperature,int vecino);
void print(FILE *out);
int i_val(int vecino,int i1);
int j_val(int vecino,int j1);

unsigned seed1 = 1946402705; // semilla generador de números aleatorios
mt19937_64 generator(seed1); // generador  mt19937_64

//Tenemos varias distribuciones de numeros aleatorios
uniform_int_distribution<int> i_distribution(0, n - 1);     //Distribución random para los i
uniform_int_distribution<int> j_distribution(1, n - 2);     //Distribución random para los j (respetar extremos en y fijos)
uniform_int_distribution<int> vecino_distribution(1,4);
uniform_real_distribution<double> r_distribution(0., 1.);   // initialize the distribution r_distribution




//Constructor
Grid grid(n, seed1);



int main()
{
    float temperature_vals[t_size];
    float temperature;
    for(int i=0;i<t_size;i++){
        temperature_vals[i]=0.01f+static_cast<float>(i)*0.1f;
    }
    char filename[31]; // Restringido a 31 caracteres, entonces solo caben en el nombre temperaturas menores de 10 y con 1 cifra decimal
    

    /*Apartado 2
    Vamos a calcular la magnetización por dominios para cada temperatura:
    El archivo de salida contiene en cada columna el promedio sobre los pasos montecarlo de la magnetización para cada dominio.
    Distitas filas corresponden a distintas temperaturas. Realmente es calcular sobre las mitades superior e inferior.
    */
    FILE *magnetization_data;
    magnetization_data = fopen("./apartados/magnet_domains.dat", "w");

    /*Apartado 4
    Calculamos la energía por partícula en función de la temperatura. Esto es equivalente a calcular grid.energy()/n^2 al final de cada 
    paso Montecarlo (aka una medida) y promediarlo sobre t_max.
    */
    FILE *epp;
    epp=fopen("./apartados/energy_per_particle.dat","w");

    /*Apartado 5
    Calculamos la densidad de partículas positivas en el grid, en la dirección y
    */
   FILE *rho;
   rho=fopen("./apartados/rho.dat","w");

    /*Apartado 6:
    Calor específico a partir de las fluctuaciones de energía. Teniendo en cuenta que <E> es la suma sobre todas las medidas 
    del hamiltoniano y divididas por el nº de medidas (t_max), se ve que <E^2> es la suma de las energías al cuadrado.
    Estos apartados, al ser simulaciones Monte Carlo, creo que en principio no pasa nada si en vez de empeze en tiempo=0
    empiezo en tiempo=100, las diferencias deberían ser mayores. Luego ya es calcular la temperatura critica
    */

    FILE *CV_dat;
    CV_dat=fopen("./apartados/specific_heat.dat","w");
    
    /*Apartado 7:
    Calculo de la susceptibilidad a partir de las fluctuaciones de magnetización en cada dominio. ¿No falta volver a 
    dividir por t_max en la formula? Mismo planteamiento que en el apartado 7
    */

    FILE *Chi_dat;
    Chi_dat = fopen("./apartados/chi.dat", "w");


    /*Apartado 8: 
    Quiero partir de una magnetización no nula
    */

    for (int t_iter = 0; t_iter < t_size; t_iter++)
    {
        // Inicializo el grid
        
        for (int tries=0; tries<10;tries ++){
            double asymetry=grid.initialize("mgt",init_mgt);
            if (abs(asymetry)<max_asymetry){
                break;                         
            }
            else{fprintf(stdout,"\n Intento %i, asimetría: %f",tries+1,asymetry);}
        }


        //fprintf(stdout, "\n \n \n Magnet 0: %f", grid.magnetization(3,init_mgt));

        temperature = temperature_vals[t_iter];
        sprintf(filename, "./ising_data/ising_dataT%.1f.dat", temperature); // si no funciona %i probar %d
        //fprintf(stdout, "\n Archivo: %s",filename);

        FILE *flips_data = fopen(filename, "w");

        if (flips_data == NULL)
        {
            perror("Error opening file: "); // Esto se asegura de que si hay algun error en la temperatura pase a la siguiente.
            t_iter++;
        }
        
        // Inicializo los resultados de la magnetización por dominios
        Grid::region result;
        result.up=0;result.down=0;
        int vecino;
        int i1,i2,j1,j2;
        double p;
        double ji;
        double energy; double m_up; double m_down;
        double sum_m_up=0; double sum_m_down=0; 
        double sum_epp=0;
        double sum_eppsrd=0;
        double sum_msqrd_up=0; double sum_msqrd_down=0;
        vector<double>density(n,0.0);
        print(flips_data);
        for (int t = 0; t < t_max; t++) // va en pasos Monte Carlo de n^2 intentos de spin flip
        {
            for (int try_ = 0; try_ < pow(n, 2); try_++)
            {
                //Dos parejas
                i1 = i_distribution(generator);
                j1=j_distribution(generator);
                vecino=vecino_distribution(generator);

                i2=i_val(vecino,i1);
                j2=j_val(vecino,j1);
                 //No vamos a calcular si el espín es el mismo.
                p = p_val(i1, i2, j1, j2, temperature, vecino);
                // Cambio de sitio (o no)
                ji = r_distribution(generator);
                if (ji < p) // No está mal, hay que tener en cuenta que al ser grid global realmente ya he cambiado los siios al calucular Eprime
                {
                    int aux;
                    aux = grid.get(i1, j1);
                    grid.set(i1, j1, grid.get(i2, j2));
                    grid.set(i2, j2, aux);
                }
                
            }
            if((t_iter+1)%10==0){
                fprintf(flips_data, "\n");
                print(flips_data);
            }
            
            //Aqui no lo estoy haciendo con el valor absoluto, aunque luego lo pruebo.
            energy=grid.energy();
            sum_epp+=energy;
            sum_eppsrd+=energy*energy;

            result = grid.magnet_regions(init_mgt);
            m_up = result.up;
            m_down = result.down;
            sum_m_up += m_up;
            sum_m_down += m_down;


            sum_msqrd_up+=pow(m_up,2); sum_msqrd_down+=pow(m_down,2);
            //fprintf(stdout, "eep: %lf \n", sum_epp);

            for (int j = 0; j < n; j++)
            {
                for (index i = 0; i < static_cast<index>(n); i++)
                {
                    density[i] += (grid.get(j,static_cast<int>(i)) + 1) / 2;
                }
            }
        }
        //Densidad definidad como (s+1)/2=0 o 1.
        fprintf(rho, "%f ", temperature);
        for(index j=0; j<static_cast<index>(n);j++){
            fprintf(rho, "%f ", density[j]/(n*t_max));
        }
        fprintf(rho,"\n");
        //fprintf(magnetization_data, "%lf", grid.energy());
        fprintf(magnetization_data,"%.1f, %f, %f \n",
        temperature, sum_m_up/t_max, sum_m_down/t_max);
        
        fprintf(epp,"%.1f, %lf, %lf \n",temperature,sum_epp/(n*n*t_max ),sqrt(sum_eppsrd-sum_epp*sum_epp/t_max)/sum_epp);
        
        fprintf(CV_dat,"%.1f, %lf \n",temperature,(sum_eppsrd/t_max-sum_epp*sum_epp/(t_max*t_max))/(n*n*temperature*temperature));
        
        fprintf(Chi_dat, "%lf, %lf, %lf  \n", temperature,
                n*n/2*(sum_msqrd_up/t_max - sum_m_up*sum_m_up/(t_max*t_max))/(temperature), // Aquí me la estoy jugando poniendo el abs a la magnetización, formula del apartado 7
                n*n/2*(sum_msqrd_down/t_max - sum_m_down*sum_m_down/(t_max*t_max)) /(temperature));
        
        fclose(flips_data);
    }
    fclose(magnetization_data);
    fclose(epp);
    fclose(CV_dat);
    fclose(Chi_dat);

    return 0;
}

int i_val(int vecino,int i1){    //tanto en esta función como en la otra estoy cogiendo numeros del 1 al 4 como los posibles 
//vecinos de un cubo 3x3 al que le he quitado el centro.
    if(i1==0){
        if(vecino==1){
            return n-1;
        }
        else if(vecino==3){
            return 1;
        }
        else{return i1;}
    }
    else if(i1==n-1){
        if(vecino==1){
            return n-2;
        }
        if(vecino==3){
            return 0;
        }
        else{return i1;}
    }
    else{
        if(vecino==1){
        return i1-1;
    }
    else if (vecino==3)
    {
        return i1+1;
    }
    else {return i1;}}
}
int j_val(int vecino,int j1){
    if(j1==1){
        if(vecino!=2){return j1;}
        //de esta forma i2,j2 es igual a i1, j1 y no hay intercambio de espin, pero la probabilidad de elección se mantiene
        //en el caso de que vecino==4
        else{return 2;}
    }
    else if (j1==n-2){
        if(vecino!=4){return j1;}//igual, si el vecino es 2, se mantiene igual 
        else{return n-3;}
    }
    else{
        if(vecino==2){
        return j1+1;
    }
    else if (vecino==4)
    {
        return j1-1;
    }
    else {return j1;}
    }
}

double p_val(int i1, int i2, int j1, int j2, float temperature, int vecino)
{
    double p;
    double deltaE;
    //Energría de la configuración
    //E=grid.energy();
    //Cambio por su vecino aleatorio y vuelvo a calcular la energía
    //aux = grid.get(i1, j1);
    //grid.set(i1, j1, grid.get(i2, j2));
    //grid.set(i2, j2, aux);
    //Eprime=grid.energy();
    if(grid.get(i1,j1)==grid.get(i2,j2)){
        p=0;
    }
    else{
        deltaE = grid.deltaE(i1, j1, i2, j2,vecino);
        // Calculo la probabiliad de cambio de posición como
        p = exp(-(deltaE) / temperature);
    }
    
    if (p >1){p=1;};
    return p;
}


void print(FILE *out)   //Solo para imprimir el grid
{
    for (int i = 0; i < n; i++)
    {
        fprintf(out, "%i", grid.get(i,0));
        for (int j = 1; j < n; j++)
        {
            fprintf(out, ",%i  ", grid.get(i,j));
        }
        fprintf(out, "\n");
    }
    fflush(out);
};


