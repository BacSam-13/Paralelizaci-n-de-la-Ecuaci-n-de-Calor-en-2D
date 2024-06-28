//Baruc Samuel Cabrera García

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
//#include <mpi.h>
#include <chrono>
#include <fstream>
#include <iostream>
#include <memory>
#include <mutex>
#include <thread>

using namespace std;
//using namespace std::chrono;

//Compilar y Ejecutar
/*
g++ Calor_2D_sec.cpp -o Calor_2D_sec
sbatch slurm_script_sec

*/

//Nota personal: Que el .log este en blanco, no significa que el codigo fallo. Se tarda en imprimir la respuesta

int main(){
    double start = clock();

    //Inicialización de valores
    int kx = 1.0; //Termino difusivo en x       ! Termino difusivo en x
    int ky = 1.0; //Termino difusivo en y        ! Termino difusivo en y

    int xI = 0.0; //Inicio dominio de x
    int xF = 1.0; //Fin dominio de x
    int yI = 0.0; //Inicio dominio de y
    int yF = 1.0; //Fin dominio de y

    int tI = 0.0; //Tiempo inicial
    int tF = 1.0; //Tiempo final
    int Nt  = 200000; //Numero de tiempos a considerar
    int NxG = 200; //Numero de puntos en x a considerar
    int NyG = 200; //Numero de puntos en y a considerar

    //Calculo de constantes
    double Dx = (xF-xI)/(NxG-1); //Distancia entre puntos en x
    double Dy = (yF-yI)/(NyG-1); //Distancia entre puntos en y
    double Dt = (tF-tI)/(Nt-1); //Disntancia entre tiempos
    double rx = kx*(Dt/(Dx*Dx)); //Constante necesaria
    double ry = ky*(Dt/(Dy*Dy)); //Constante necesaria

    //Condicion de error
    if(rx + ry < 1/2){
        cout<<"ERROR"<<endl;
        return 0;
    }

    /*
    Sea A una matriz de dimension mxn, entonces se tiene que
    A[i][j] = B[i*n + j]
    Donde B es un vector de tamaño m*n

    En este caso, hacemos la convección de que la casilla i,j de la matriz malla (sentido cartesiano) 
    esta dada por [i*Nyg + j]
    i recorrera el eje x, y j el eje y
    */

    //Declaramos las "matrices" para las iteraciones
    // De forma iteraitva, mu_old sera n, y mu_new sera n+1
    // Usaremos puntoeros inteligentes para facilitar ciertos procesos (Eliminacion automatica)
    //double *mu_old = (double *)malloc((NxG*NyG) * sizeof(double));
    //double *mu_new = (double *)malloc((NxG*NyG) * sizeof(double));
    shared_ptr<double> MU_OLD = shared_ptr<double>(new double[NxG*NyG], default_delete<double[]>());
    shared_ptr<double> MU_NEW = shared_ptr<double>(new double[NxG*NyG], default_delete<double[]>());

    double *mu_old = MU_OLD.get();
    double *mu_new = MU_NEW.get();

    //Inicializacion de la malla vieja
    /*
    Notemos que este tipo de recorrido lo que hace es recorrer 
    columnas de izquierda a derecha, de abajo hacia arriba
    i-th columna, j-th fila
    */
    for(int i = 0; i < NxG; i++){
        for(int j = 0; j < NyG; j++){
            mu_old[i*NyG + j] = sin(( xI + (i)*Dx) + (yI + (j)*Dy));
            mu_old[i*NyG + j] *= mu_old[i*NyG + j];
        }
    }

    //Condiciones de frontera en la malla nueva
    for(int i = 0; i < NxG; i++){
        mu_new[i*NyG + 0] = 1;//Frontera sur
        mu_new[i*NyG + (NyG-1)] = 1;//Frontera norte
    }

    for(int j = 0; j < NyG; j++){
        mu_new[0*NyG + j] = 1;//Frontera oeste
        mu_new[(NxG-1)*NyG + j] = 1;//Frontera este
    }


    //Iteraciones
    for(int t = 0; t < Nt; t++){
        for(int i = 1; i < NxG-1; i++){
            for(int j = 1; j < NyG-1; j++){
                double a1 = mu_old[i*NyG + j];
                double a2 = rx*(mu_old[(i+1)*NyG + j] + 2*mu_old[i*NyG + j] + mu_old[(i-1)*NyG + j]);
                double a3 = ry*(mu_old[i*NyG + j+1] + 2*mu_old[i*NyG + j] + mu_old[i*NyG + j-1]);
                mu_new[i*NyG + j] = a1 + a2 + a3;
            }
        }

        //Debido a que calculamos la nueva matriz, es momento de actualizar la vieja
        for(int i = 0; i < NxG; i++){
            for(int j = 0; j < NyG; j++){
                mu_old[i*NyG + j] = mu_new[i*NyG + j];
            }
        }
    }

    //Calculo de Tiempo Final
    double end = clock();
    double total = (end - start)/CLOCKS_PER_SEC;
    cout<<"Tiempo total -> "<< total<<endl;

    //free(mu_new);
    //free(mu_old);

    return 0;
}