#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#define n 1
#define D 1.0
#define tEND 10
#define xEND M_PI
#define lam .4
#define dx (M_PI/10)
#define dt (lam/D*dx*dx)

#pragma warning(disable:4996)
int xSize = xEND/dx +1;
int tSize = (double)tEND/dt +1;
char fileName[100] = "Explicit1Dheat.csv";
char fileNameTimeSepCSV[100];
char fileNameBlockCSV[100];
char fileNameBlockActualCSV[100];
char errorFile[101];

double F( double t ,double x){ return 0; }
void initializeU( double u[xSize][tSize] );
void boundaryConditions( double u[tSize][xSize],int);
void fillRow( double[tSize][xSize], int );
void printAll( double u[tSize][xSize] );
double exact( double, double);

int main()
{
    double U[tSize][xSize];

    while (true){

        strcpy(fileNameBlockCSV, fileName);
        strcpy(errorFile, fileName);
        strcpy(fileNameBlockActualCSV, fileName);
        strcat(fileNameBlockCSV, "BLOCK.csv");
        strcat(errorFile, "ERROR.csv");
        strcat(fileNameBlockActualCSV, "BLOCKACTUAL.csv");

        initializeU( U );

        int i;
        for( i=1; i<tSize; i++ ){
            boundaryConditions( U, i );

            fillRow( U, i );
        }


        printAll( U );
        if ( lam>=0.5 ){
            printf("\n\nThe algorithm is unstable. Lambda = %1.5lf", lam);
        }
        printf("\n\nLambda: %lf\n", lam);

        printf("\nWould you like to continue again? (y/n) ");
        char continu;
        scanf(" %c", &continu);
        if (continu == 'n'){
            return 0;
        }
    }

    return 0;
}

void printAll( double u[tSize][xSize] ){
    int i, j;

    FILE *ef = fopen(errorFile, "w");
    FILE *bf = fopen(fileNameBlockCSV, "w");
    FILE *abf = fopen(fileNameBlockActualCSV, "w");

    for( i=0; i<tSize; i++ ){

        for( j=0; j<xSize; j++ ){
            if( j==0 )printf("ESTIMATE:\t");
            printf("\t%lf", u[i][j] );
            fprintf(bf,"%lf,", u[i][j]);
        }
        fprintf(bf,"\n");

        for( j=0; j<xSize; j++ ){
            if( j==0 )printf("\nACTUAL:\t\t");
            printf("\t%lf", exact(i*dt, j*dx) );
            fprintf(abf, "%lf,", exact(i*dt, j*dx) );
        }
        fprintf(abf,"\n");

        for( j=0; j<xSize; j++ ){
            if( j==0 )printf("\nPercent ERROR:\t");
            if( j!=0 && j!=xSize-1)
                fprintf(ef, "%lf,", (exact(i*dt,1*dx)-u[i][1])/exact(i*dt,1*dx)*100 );
        }
        fprintf(ef, "\n");

        printf("\n\n\n\n");
    }

    fclose(ef);
    fclose(abf);
    fclose(bf);

}
void initializeU( double u[tSize][xSize] ){
    int i;
    for( i=0; i<xSize; i++ ){
        u[0][i] = sin(i*dx*n);
    }
}
void boundaryConditions( double u[tSize][xSize], int tRow ){
    u[tRow][0] = 0;
    u[tRow][xSize-1] = 0;
}
void fillRow( double u[tSize][xSize], int tRow){
    int i;
    for( i=1; i<xSize-1; i++ ){
        u[tRow][i] = u[tRow-1][i] + lam*( u[tRow-1][i-1] - 2*u[tRow-1][i] + u[tRow-1][i+1] ) + dt*F(i*dt, i*dx);
    }
}
double exact( double t, double x ){ return exp(-n*n*t)*sin(n*x); }
