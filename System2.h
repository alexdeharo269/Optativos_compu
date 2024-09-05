#include<stdio.h>
#include<math.h>
#include <iostream>
#include <fstream>
#include <chrono>
#include<random>
#include<cmath>
#include<vector>
#include<algorithm>

using namespace std;
using index_t = vector<int>::size_type; // Define index_t as the size type of vector<int>

class Planet {
    public:
        double dtS;
        double phi;
        double vr;
        double vphi;
        double mass;
        double radius;
        double U=0; //energia interna.
        int rock=0; //
        bool real=true;
        float years=0; //nº de vueltas al sol.
        int interacciones=0;
        double periodo=0;
        double exc;
        double avgdtS;

        // Constructor
        Planet()=default;

        //Paso a cartesianas 
        double x(){
            double x;
            x=dtS*cos(phi);
            return x;
        }
        double y()
        {
            double y;
            y = dtS * sin(phi);
            return y;
        }
        
};

class System {
    private:
        int numPlanetas;
        float D_max,D_min, rock_density;
        double vrmod,radmult;
        double asteroidmult;
        mt19937_64 generator;
        uniform_real_distribution<double> real_distribution;
        uniform_real_distribution<double> rad_distribution;
        uniform_real_distribution<double> vr_distribution;
        const double pi = M_PI;

    public:
        vector<Planet> plts;
        int interacciones=0;
        int absorbidos=0;
        int supervivientes;
        int n_old=0;

        // Constructor
        System(int nplt, unsigned seed, float Dmax,float Dmin, float rock_den, double vr_mod, double rad_mult, float asteroid_mult)
            : numPlanetas(nplt), D_max(Dmax),D_min(Dmin), rock_density(rock_den), vrmod(vr_mod), radmult(rad_mult),
              asteroidmult(asteroid_mult), generator(seed),
              real_distribution(0.0, 1.0), rad_distribution(D_min / D_max, 1.0), // pongo otra distribución para el radio ya que no pueden tener radio 0.
              // el menor radio será la órbita de Mercurio.
              vr_distribution(-1.0, 1.0),
              plts(static_cast<index_t>(numPlanetas))
        {
        }

        // Método para generar posiciones aleatorias para los planetas
        void Parametrosiniciales()
        {
            double p=0;
            // Generar posiciones aleatorias dentro del rango especificado
            for (index_t i = 0; i < static_cast<index_t>(numPlanetas); ++i)
            {
                // Multiplicamos por una función de densidad que da cuenta de la probabiliad de que
                // en el disco de formación del SS los planetas estén a cierta distancia.
                p = rad_distribution(generator);
                plts[i].dtS=pow(D_min/D_max,1-p);
                plts[i].dtS *=D_max;
                //plts[i].dtS = rad_distribution(generator) * D_max;

                plts[i].phi = real_distribution(generator) * 2 * pi;

                plts[i].vr = vrmod*vr_distribution(generator); // PASO 2
                plts[i].radius=0.1;
                //Para la velocidad angular vamos a calcular una velocidad angular (radianes/59.1días) en función de la distancia
                //Usamos la velocidad de una trayectoria circular:
                
                plts[i].mass=0.0014/numPlanetas;
                plts[i].vphi = sqrt(1/ plts[i].dtS);
                plts[i].radius=0.85*4.264e-5*1000/numPlanetas*radmult;
                if(real_distribution(generator)<=rock_density){plts[i].rock=1;}

                
            }
        }

        void Mover(double h, FILE* data,int iter_max){
            
            
            //Algoritmo de Verlet
            index_t nplt = static_cast<index_t>(numPlanetas);
            vector<double> ax(nplt,0.0), ay(nplt,0.0),wx(nplt,0.0),wy(nplt,0.0);

            vector<double> energy(nplt, 0.0);
            vector<double> y_tminus1(nplt);   // Aqui guardo la y del instante previo para calcular el periodo.
            vector<double> aux(nplt, 0.0);    // Variable auxiliar para el tiempo (ultima vez que pasó el planeta por y=0).
            vector<double> vueltas_sin_iter(nplt,0.0);   //Variable auxiliar para el tiempo que lleva un planeta sin interaccionar.
            vector<int> tiempo_sin_iter(nplt,0);
            vector<int> aux_iter(nplt,0);
            vector<double> dist(nplt,0.0);
            vector<double> exc(nplt,0.0);

            // Initialize positions and velocities
            vector<vector<double>> pic(nplt, vector<double>(4, 0.0));
            
            //Paso a cartesianas
            for (index_t i = 0; i < static_cast<index_t>(nplt); ++i){
                pic[i][0] = plts[i].x();
                pic[i][1] = plts[i].y();
                pic[i][2] = -(plts[i].vphi * plts[i].y() + plts[i].vr * plts[i].x()) / plts[i].dtS; // vx = -vphi * y
                pic[i][3] = (plts[i].vphi * plts[i].x() + plts[i].vr * plts[i].y()) / plts[i].dtS;  // vy = vphi * x

                ax[i] = - plts[i].x() / (plts[i].dtS * plts[i].dtS * plts[i].dtS);
                ay[i] = - plts[i].y() / (plts[i].dtS * plts[i].dtS * plts[i].dtS);
            }
            //Algoritmo de Verlet
            for(double t=0;t<iter_max;t+=h){
                for (index_t i = 0; i < static_cast<index_t>(nplt); ++i)
                {
                    if (plts[i].real)
                    {
                        y_tminus1[i] = pic[i][1]; // Guardo la y anterior para calcular el periodo

                        pic[i][0] = pic[i][0] + h * pic[i][2] + pow(h, 2) / 2 * (ax[i]); // r(t+h)=r(t)+hv(t)+h^2/2a(t)
                        pic[i][1] = pic[i][1] + h * pic[i][3] + pow(h, 2) / 2 * (ay[i]); //

                        wx[i] = pic[i][2] + h / 2 * ax[i]; // variables auxialiares para poder guardar ax[t]
                        wy[i] = pic[i][3] + h / 2 * ay[i];

                        ax[i] = -pic[i][0] / sqrt(pow(pow(pic[i][0], 2) + pow(pic[i][1], 2), 3));
                        ay[i] = -pic[i][1] / sqrt(pow(pow(pic[i][0], 2) + pow(pic[i][1], 2), 3));

                        pic[i][2] = wx[i] + h / 2 * ax[i]; // v(t+h)=w(t)+h/2a(t+h)
                        pic[i][3] = wy[i] + h / 2 * ay[i];
                    }
                    else
                    {
                        pic[i][0] = 0.0;
                        pic[i][1] = 0.0;
                    }
                }
                
                for (index_t i = 0; i < static_cast<index_t>(nplt); ++i)
                {
                    plts[i].dtS = sqrt(pic[i][0] * pic[i][0] + pic[i][1] * pic[i][1]);
                    plts[i].phi=atan2(pic[i][1],pic[i][0]);

                    plts[i].vr = (pic[i][0] * pic[i][2]+pic[i][1]*pic[i][3]) / sqrt(pic[i][0] * pic[i][0] + pic[i][1] * pic[i][1]);
                    plts[i].vphi = (pic[i][0] * pic[i][3] - pic[i][1] * pic[i][2]) / sqrt(pic[i][0] * pic[i][0] + pic[i][1] * pic[i][1]);
                    //Printear en planets_data
                    //Imprimir cada 10*h
                    if(static_cast<int>(round(t))%5==0){
                        fprintf(data, "%i,   %lf, %lf, %lf, %i\n", plts[i].real, pic[i][0], pic[i][1], plts[i].radius, plts[i].rock);
                    }
                    
                }
                fprintf(data, "\n");
                

                    Interacciones();

                /*Calculo el periódo en años terrestres de los planetas.
                Consideramos que podemos hablar de periodo cuando el planeta no ha interaccionado en la última órbita
                y esta ya está formada. Con este bucle calculo tan sólo el último periodo instantáneo, por lo que estoy
                interesado en la media de los periódos desde que el planeta no interacciona.
                */
                
                for (index_t i = 0; i < nplt; i++)
                {
                    if((plts[i].real)&(plts[i].interacciones==aux_iter[i])){
                        if (signbit(y_tminus1[i]) != signbit(pic[i][1])) // signbit(): bit del signo de un float (de cmath)
                        {
                            vueltas_sin_iter[i] += (t - aux[i]) * 59.1 / 365.25;
                            aux[i] = t;
                            plts[i].years+=0.5f;
                            
                        }
                        tiempo_sin_iter[i]++;
                        dist[i]+=plts[i].dtS;

                        exc[i] += sqrt(1 + 2*plts[i].dtS *plts[i].dtS* plts[i].vphi *plts[i].vphi* 
                        (0.5*(plts[i].vr * plts[i].vr + plts[i].vphi * plts[i].vphi)-1/plts[i].dtS));
                    }
                    else{
                        aux_iter[i]=plts[i].interacciones;
                        dist[i]=0; exc[i]=0;
                        vueltas_sin_iter[i]=0;
                        tiempo_sin_iter[i]=0;
                        plts[i].years=0;
                    }
                }
                
                
            }

            
            for(index_t i=0; i<nplt; i++){
                if(plts[i].years>=0.99){     //Para no tener en cuenta los planetas que no llegan a dar una vuelta completa
                    plts[i].periodo = vueltas_sin_iter[i] / plts[i].years;
                    plts[i].avgdtS=dist[i]/tiempo_sin_iter[i];
                    plts[i].exc=exc[i]/tiempo_sin_iter[i];
                    n_old++;
                    }
            }
            supervivientes=numPlanetas-interacciones-absorbidos;
        }

        void Sizes(FILE *energies, FILE *sizes, FILE *multiverse, double h){
            vector<vector<double>> plt(static_cast<index_t>(supervivientes),vector<double>(7,0.0));
            int smalls=0; int mediums=0; int larges=0; int very_larges=0;
            int old_small=0; int old_mediumns=0;int old_larges=0;int old_very_larges=0; //para tener en cuenta los que tienen más de un año
            index_t i,j;

            
            index_t real_i=0;
            for (i = 0; i < static_cast<index_t>(numPlanetas); ++i){
                if(plts[i].real){
                    plt[real_i][1]=plts[i].radius;
                    plt[real_i][2]=plts[i].periodo;
                    plt[real_i][3]=plts[i].exc;
                    plt[real_i][4]=plts[i].avgdtS;
                    plt[real_i][5]=plts[i].U;
                    plt[real_i][6]=static_cast<double>(plts[i].rock);
                    real_i++;
                }
            }
             // Ordenar según la segunda columna, el radio.
             sort(plt.begin(), plt.end(), [](const vector<double> &a, const vector<double> &b)
                  { return a[1] < b[1]; });

             index_t q1_index = static_cast<index_t>(supervivientes / 4);
             index_t q2_index = static_cast<index_t>(supervivientes / 2);
             index_t q3_index = static_cast<index_t>(3 * supervivientes / 4);

             double q1 = plt[q1_index][1];
             double q2 = plt[q2_index][1];
             double q3 = plt[q3_index][1];

             for (i = 0; i < static_cast<index_t>(supervivientes); ++i) {
             if (plt[i][1] <= q1) {plt[i][0] = 1.0;smalls++;}
             else if (plt[i][1] <= q2) {plt[i][0] = 2.0;mediums++;}
             else if (plt[i][1] <= q3) {plt[i][0] = 3.0;larges++;}
             else {plt[i][0] = 4.0;very_larges++;}
             }

             // Sumas de las columnas de plt para calcular los valores medios de cada propiedad para cada tipo de planeta
             vector<double> sum_small(7, 0.0), sum_medium(7, 0.0), sum_large(7, 0.0), sum_very_large(7, 0.0);

            for (i = 0; i < static_cast<index_t>(supervivientes); ++i)
            {
                switch (static_cast<int>(plt[i][0]))
                {
                case 1: // small
                    for (j = 1; j < 7; ++j){
                        if((j==2||j==3||j==4)&&(plt[i][2]!=0.0)){
                            sum_small[j] += plt[i][j];
                            if(j==2){old_small++;}
                        }else{sum_small[j] += plt[i][j];}
                        }
                    break;
                case 2: // medium
                    for (j = 1; j < 7; ++j){
                        if((j==2||j==3||j==4)&&(plt[i][2]!=0.0)){
                            sum_medium[j] += plt[i][j];
                            if(j==2){old_mediumns++;}
                        }
                        else{sum_medium[j] += plt[i][j];}}
                    break;
                case 3: // large
                    for (j = 1; j < 7; ++j){
                        if((j==2||j==3||j==4)&&(plt[i][2]!=0.0)){
                            sum_large[j] += plt[i][j];
                            if(j==2){old_larges++;}
                        }
                        else{sum_large[j] += plt[i][j];}}
                    break;
                case 4: // very large
                    for (j = 1; j < 7; ++j){
                        if((j==2||j==3||j==4)&&(plt[i][2]!=0.0)){
                            sum_very_large[j] += plt[i][j];
                            if(j==2){old_very_larges++;}
                        }
                        else{sum_very_large[j] += plt[i][j];}}
                     break;
                }
            }

            // Calculate averages
            auto calc_avg = [](const vector<double> &sum, int count,int old_count)
            {
                vector<double> avg(7, 0.0);
                if (count > 0)
                {
                    for (index_t k = 1; k < 7; ++k)
                    {
                        if (k == 2 || k == 3 || k == 4) {
                           avg[k] = sum[k] / old_count;
                        } else {
                            avg[k] = sum[k] / count;
                        }
                    }
                }
                avg[6] = (sum[6] / count) * 100; // Rock composition as a percentage
                return avg;
            };

            vector<double> avg_small = calc_avg(sum_small, smalls,old_small);
            vector<double> avg_medium = calc_avg(sum_medium, mediums,old_mediumns);
            vector<double> avg_large = calc_avg(sum_large, larges,old_larges);
            vector<double> avg_very_large = calc_avg(sum_very_large, very_larges,old_very_larges);

            // Print to file
            //fprintf(sizes, "%lf, %lf  ", vrmod, radmult, )
            fprintf(sizes, "Categoría & Totales(Estables) & Radio &    Periodo &   Exc &    Órbita(UA) &    U &   Rocoso(\\%)\\\\\n");

            fprintf(sizes, "Pequeños&      %i(%i)&  %.2lf& %.1lf& %.2lf& %.2lf& %e& %.1lf\\\\\n", smalls, old_small, avg_small[1], avg_small[2], avg_small[3], avg_small[4], avg_small[5], avg_small[6]);
            fprintf(sizes, "Medianos&     %i(%i)&  %.2lf& %.1lf& %.2lf& %.2lf& %e& %.1lf\\\\\n", mediums, old_mediumns, avg_medium[1], avg_medium[2], avg_medium[3], avg_medium[4], avg_medium[5], avg_medium[6]);
            fprintf(sizes, "Grandes&     %i(%i)&  %.2lf& %.1lf& %.2lf& %.2lf& %e& %.1lf\\\\\n", larges, old_larges, avg_large[1], avg_large[2], avg_large[3], avg_large[4], avg_large[5], avg_large[6]);
            fprintf(sizes, "Muy grandes& %i(%i)&  %.2lf& %.1lf& %.2lf& %.2lf& %e& %.1lf\\\\\n", very_larges, old_very_larges, avg_very_large[1], avg_very_large[2], avg_very_large[3], avg_very_large[4], avg_very_large[5], avg_very_large[6]);
            fprintf(sizes, "Iniciales & Interacciones & Absorciones & Supervivientes(Estables) & Volumen \\\\\n %i &     %i &     %i   &    %i(%i)   &   %.2lf",
                    numPlanetas, interacciones, absorbidos, supervivientes, n_old, static_cast<double>(numPlanetas) / supervivientes);


            fprintf(energies, "i,  Size,    Rad        T       Exc       avgdts      U      Rock\n");

            for (i = 0; i < static_cast<index_t>(supervivientes); ++i){
               fprintf(energies, "%i, %lf, %lf, %lf, %lf, %lf, %lf, %lf\n", static_cast<int>(i), plt[i][0], plt[i][1], plt[i][2],
                       plt[i][3], plt[i][4],
                       plt[i][5], plt[i][6] );
            }
            
            fprintf(multiverse, "\n%.1lf %.4lf %.1lf", h, vrmod, radmult);
            fprintf(multiverse, " %i %i %i %i %i %.3lf", numPlanetas, interacciones, absorbidos, supervivientes, n_old, static_cast<double>(numPlanetas) / supervivientes);
            fprintf(multiverse, " %i %i %e %e %e %e %e %.1lf", smalls,old_small, avg_small[1], avg_small[2], avg_small[3], avg_small[4], avg_small[5], avg_small[6]);
            fprintf(multiverse, " %i %i %e %e %e %e %e %.1lf", mediums,old_mediumns, avg_medium[1], avg_medium[2], avg_medium[3], avg_medium[4], avg_medium[5], avg_medium[6]);
            fprintf(multiverse, " %i %i %e %e %e %e %e %.1lf", larges, old_larges, avg_large[1], avg_large[2], avg_large[3], avg_large[4], avg_large[5], avg_large[6]);
            fprintf(multiverse, " %i %i %e %e %e %e %e %.1lf", very_larges, old_very_larges, avg_very_large[1], avg_very_large[2], avg_very_large[3], avg_very_large[4], avg_very_large[5], avg_very_large[6]);
            
        }
        void Colision(index_t i, index_t j)
        {

                 plts[i].radius = plts[i].radius * pow((plts[i].mass + plts[j].mass) / plts[i].mass, 1.0 / 3.0);
                 plts[i].vphi = (plts[i].vphi * plts[i].mass + plts[j].vphi * plts[j].mass) / (plts[j].mass + plts[i].mass);
                 plts[i].vr = (plts[i].mass * plts[i].vr + plts[j].mass * plts[j].vr) / (plts[j].mass + plts[i].mass);
                 plts[i].mass += plts[j].mass;
                 // Energía interna como valor absoluto de la diferencia de enegía cinética mas la energía interna del otro planeta.
                 plts[i].U += 0.5 * abs(plts[i].mass * (plts[i].vr * plts[i].vr + plts[i].vphi * plts[i].vphi) -
                                        plts[j].mass * (plts[j].vr * plts[j].vr + plts[j].vphi * plts[j].vphi)) +
                              plts[j].U;
                 plts[j].U = 0;
                 plts[j].mass = 0;
                 plts[j].real = false;
                 plts[j].radius = 0.0;

                 //

        }

        void Interacciones(){
            index_t i,j;
            double d;
            for(i=0; i<static_cast<index_t>(numPlanetas);i++){
                if(plts[i].real){
                    if (plts[i].rock == 1)
                    {
                        for (j = 0; j < i; j++)
                        {
                            if((plts[j].real)&(plts[j].rock==1)){
                                d = sqrt(plts[i].dtS * plts[i].dtS + plts[j].dtS * plts[j].dtS - 2 * plts[i].dtS * plts[j].dtS * cos(plts[j].phi - plts[i].phi));

                                if (d < (plts[i].radius + plts[j].radius))
                                {
                                    interacciones++;
                                    plts[i].interacciones++;
                                    Colision(i, j);
                                }
                            }
                        }
                    }
                    if ((plts[i].rock == 0) && (plts[i].dtS >= 2.7 * asteroidmult))
                    {
                        for (j = 0; j < i; j++)
                        {
                            if((plts[j].real)&(plts[j].rock==0)){
                                d = sqrt(plts[i].dtS * plts[i].dtS + plts[j].dtS * plts[j].dtS - 2 * plts[i].dtS * plts[j].dtS * cos(plts[j].phi - plts[i].phi));

                                if (d < (plts[i].radius + plts[j].radius))
                                {
                                    interacciones++;
                                    plts[i].interacciones++;
                                    Colision(i, j);
                                }
                            }
                        }
                    }
                    //Consideramos aquí interacciones con el Sol, planetas absorbidos por éste. 
                    //El radio solar en UA es aproximadamente 0.005
                    if(plts[i].dtS<=0.005){
                        plts[i].real=false;
                        absorbidos++;
                    }
                }
            }

        }

};