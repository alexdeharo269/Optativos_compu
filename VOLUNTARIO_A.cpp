#include<stdio.h>
#include<math.h>
#include <iostream>
#include <fstream>
#include<random>
#include<cmath>
#include<vector>
#include "System2.h"
#include<filesystem>

//Radio solar (absorbidos) =0.005
//Los planetas no pueden aparecer antes que mercurio
//Imprimo las cosas al fichero cada 10*h

const int nplt = 1000;      // esto va a ser el numero inicial de planetesimales.
const int iter_max = 10*round(165*365/59.1); // nº máximo de iteraciones, lo pongo para que un planeta tan lejano como neptuno de al menos 10 vueltas

using namespace std;
string formatDouble(double value, int precision);
//Para ahorrar computación lo hago con una seed y con un h y luego lo cambio individualmente para ver las diferencias con un mismo valor.



int main(){
    FILE *multiverse = fopen("resumen.dat", "w");

    vector<unsigned int> seed_values = {11546458}; //Valores seed
    vector<double> h_values = {0.5};                           //Posibles valores h
    vector<double> vr_values = {0.001, 0.001, 0.01, 0.1 ,1.0,10.0};   //Posibles valores vr aleatoria
    vector<double> radii_values = {10, 100, 1000, 10000,100000};                 //Posibles valores radio

    for (unsigned seed : seed_values)
    {
        for (double h : h_values)
        {
            for (double vr : vr_values)
            {
                for (double radii : radii_values)
                {
                    // Construct unique filenames based on parameter combinations
                    string file_suffix = "_seed_" + to_string(seed) + "_h_" + formatDouble(h,1) +
                                         "_vr_" + formatDouble(vr,3) + "_radii_" + formatDouble(radii,1) + ".dat";

                    fprintf(stdout, "%s\n", file_suffix.c_str());
                    
                    // Construct file paths
                    string data_filename = "./data/planets_data" + file_suffix;
                    string energies_filename = "./energies/planets_energies" + file_suffix;
                    string sizes_filename = "./sizes/planets_size" + file_suffix;

                    // Open files
                    FILE *data = fopen(data_filename.c_str(), "w");
                    FILE *energies = fopen(energies_filename.c_str(), "w");
                    FILE *sizes = fopen(sizes_filename.c_str(), "w");
                    // Initialize and run system
                    System system(nplt, seed, 30.0f, 0.39f, 0.1f, vr,radii , 1.0f);

                    system.Parametrosiniciales();
                    system.Mover(h, data, iter_max);
                    system.Sizes(energies, sizes,multiverse, h);


                    // Close files
                    fclose(data);
                    fclose(energies);
                    fclose(sizes);
                }
            }
        }
    }
    fclose(multiverse);
    return 0;
}

string formatDouble(double value, int precision)
{
    stringstream stream;
    stream << fixed << setprecision(precision) << value;
    return stream.str();
}