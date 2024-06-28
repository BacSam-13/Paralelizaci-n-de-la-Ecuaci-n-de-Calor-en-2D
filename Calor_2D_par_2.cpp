//Baruc Samuel Cabrera García

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>
#include <chrono>
#include <fstream>
#include <memory>

using namespace std;
using namespace std::chrono;

//Compilar y Ejecutar
/*
mpic++ Calor_2D_par_2.cpp -o Calor_2D_par_2
sbatch slurm_script_par_2

*/


int main(){
    double start = clock();

    int numtasks, taskid;

    MPI_Init(NULL, NULL); //Inicializacion
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks); //Numero de procesos
    MPI_Comm_rank(MPI_COMM_WORLD, &taskid); // Identificador de cada proceso

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

    //Variables para SEND_RECV
    MPI_Status statut;
    MPI_Comm comm;
    int etiqueta = 13;

    //Condicion de error
    if(rx + ry < 1/2){
        cout<<"ERROR"<<endl;
        return 0;
    }


    //Calculo y Declaracion de variables para procesos

    //Bloques vecinos
    int vec[2];
    vec[0] = MPI_PROC_NULL;
    vec[1] = MPI_PROC_NULL;

    //Variables para la separación porbloques
    int Ncol = NxG/numtasks;//Numero de columnas para este proceso
    int i_I = taskid*Ncol;//A partir de que columna empezaremos
    int i_F;//Columna donde terminaremos
    /*
    Recordemos que en el caso en que el numero de columnas no sea divisible entre el numero de procesos.
    Se tiene que considerar el ultimo proceso con un tamaño distinto
    */
    if(taskid == numtasks - 1){
        i_F = NxG;
    }else{
        i_F = i_I + Ncol;//Columna en la que terminaremos
    }

    Ncol = i_F - i_I;//Debido al if anterior, debemos volver a actualizar


    //Divisiones de los dominios
    int Nx, Ny;
    if (numtasks == 1){
        Nx = NxG;
    }else if (taskid == 0){
        Nx = Ncol;
    }else if(taskid == numtasks-1){
        Nx = NxG - Ncol*taskid + 2;
    }else{
        Nx = Ncol + 2;
    }
    Ny = NyG;
    
    //Indice global para la inicializacion de la malla
    int *index_global = (int *)malloc((Nx) * sizeof(int));
    if (numtasks == 1){
        for(int i = 0; i < Nx; i++){
            index_global[i] = i;
        }
    }else{
        if (taskid == 0){
            for(int i = 0; i < Nx; i++){
                index_global[i] = i;
            }
        }else{
            for(int i = 0; i < Nx; i++){
                index_global[i] = taskid*(Ncol)-1 + (i-1);
            }
        }
    }
    double *xG = (double *)malloc((NxG) * sizeof(double));
    double *x = (double *)malloc((Nx) * sizeof(double));
    double *yG = (double *)malloc((NyG) * sizeof(double));
    double *y = (double *)malloc((Ny) * sizeof(double));

    for(int i = 0; i < NxG; i++){
        xG[i] = xI + (i-1)*Dx;
    }
    for(int j = 0; j < NyG; j++){
        yG[j] = yI + (j-1)*Dy;
    }

    for(int i = 0; i < Nx; i++){
        x[i] = xG[index_global[i]];
    }

    for(int j = 0; j < Ny; j++){
        y[j] = yG[j];
    }



    //Declaramos las "matrices" para las iteraciones
    // De forma iteraitva, mu_old sera n, y mu_new sera n+1.
    double *mu_old = (double *)malloc((Nx*Ny) * sizeof(double));
    double *mu_new = (double *)malloc((Nx*Ny) * sizeof(double));

    //Inicializacion de la malla vieja
    /*
    Notemos que este tipo de recorrido lo que hace es recorrer 
    columnas de izquierda a derecha, de abajo hacia arriba
    i-th columna, j-th fila
    */
    for(int i = 0; i < Nx; i++){
        for(int j = 0; j < Ny; j++){
            mu_old[i*Ny + j] = sin(x[i] + y[j]);
            mu_old[i*Ny + j] *= mu_old[i*Ny + j];
            mu_new[i*Ny + j] = mu_old[i*Ny + j];
        }
    }

    //Condiciones de frontera en la malla nueva
    for(int j = 0; j < Ny; j++){
        mu_new[0*Ny + j] = 1; //Oeste
        mu_new[(Nx-1)*Ny + j] = 1; //Este
    }
    for(int i = 0; i < Nx; i++){
        mu_new[i*Ny + 0] = 1; //Sur
        mu_new[i*Ny + (Ny-1)] = 1; //Norte
    }


    //Iteraciones
    for(int t = 0; t < Nt; t++){

        //Calculo de la nueva iteracion, considerando los bloques de este proceso
        for(int i = 1; i < Nx-1; i++){
            for(int j = 1; j < Ny-1; j++){
                double a1 = mu_old[i*NyG + j];
                double a2 = rx*(mu_old[(i+1)*Ny + j] + 2*mu_old[i*Ny + j] + mu_old[(i-1)*Ny + j]);
                double a3 = ry*(mu_old[i*Ny + j+1] + 2*mu_old[i*Ny + j] + mu_old[i*Ny + j-1]);
                mu_new[i*Ny + j] = a1 + a2 + a3;
            }
        }

        //Comunicacion
        if(numtasks > 1){
            //Enviamos el borde izquierdo, y recibimos un borde izquierda de la derecha.
			MPI_Sendrecv(&mu_new[1*Ny + 0], Ny, MPI_DOUBLE, vec[0], etiqueta, &mu_new[(Nx-1)*Ny + 0], Ny, MPI_DOUBLE, vec[1], etiqueta, comm, &statut);
			
			//Enviamos el borde derecho, y recibimos un borde derecho de la izquierda.
			MPI_Sendrecv(&mu_new[((Nx-1)-1)*Ny + 0], Ny, MPI_DOUBLE, vec[1], etiqueta, &mu_new[0*Ny + 0], Ny, MPI_DOUBLE, vec[0], etiqueta, comm, &statut);
		}


        //Debido a que calculamos la nueva matriz, es momento de actualizar la vieja
        for(int i = 0; i < Nx; i++){
            for(int j = 0; j < Ny; j++){
                mu_old[i*NyG + j] = mu_new[i*NyG + j];
            }
        }        


    }

    //Calculo de tiempo final
    double end = clock();
    
    //Liberamos memoria usada
    free(x);
    free(y);
    free(xG);
    free(yG);
    free(mu_old);
    free(mu_new);
    free(index_global);

    //Finalizamos MPI
    MPI_Finalize();

    //Imprimimos tiempos finales
    double total = (end - start)/CLOCKS_PER_SEC;
    cout<<"Tiempo total -> "<<total<<endl;
}