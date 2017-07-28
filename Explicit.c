
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#define string.h
#define n 1
#define D 1.0
#define tMAX 10
#define xMAX M_PI // does not work in VS 2017
#define lam .4 //unstable after 0.5
#define dx (M_PI/10)
#define dt (lam/D*dx*dx)

int xstep = xMAX / dx + 1;
int tstep = (double)tMAX / dt + 1;
char fileName[100] = "Exp1D.csv";
char fileNameBlockCSV[100];
char fileNameBlockActualCSV[100];
char errorFile[101];

double F(double, double);
void U0(double u[xstep][tstep]);
void Boundary(double u[tstep][xstep], int tRow);
void fill(double[tstep][xstep], int);
void print(double u[tstep][xstep]);
double exact(double, double);

int main(void)
{
	double U[tstep][xstep];

	U0(U);
	int i;
	for (i = 1; i<tstep; i++)
	{
		Boundary(U, i);
		fill(U, i);
	}

	print(U);
	if (lam >= 0.5)
	{
		printf("\n\ algorithm is unstable because Lambda is = %1.5lf", lam);
	}

	return 0;
}
void print(double u[tstep][xstep])
{
	int i, j;


	FILE *error = fopen("error.csv", "w");
	FILE *estimate = fopen(fileName, "w");
	FILE *actual = fopen("actual.csv", "w");

	for (i = 0; i<tstep; i++)
	{

		for (j = 0; j<xstep; j++)
		{
			printf("\t%lf", u[i][j]);
			fprintf(estimate, "%lf,", u[i][j]);
		}
		fprintf(estimate, "\n");

		for (j = 0; j<xstep; j++)
		{
			printf("\t%lf", exact(i*dt, j*dx));
			fprintf(actual, "%lf,", exact(i*dt, j*dx));
		}
		fprintf(actual, "\n");

		for (j = 0; j<xstep; j++)
		{
			printf("\t%lf", (exact(i*dt, j*dx) - u[i][j]) / exact(i*dt, j*dx) * 100);
			if (j != 0 && j != xstep - 1)
				fprintf(error, "%lf,", (exact(i*dt, 1 * dx) - u[i][1]) / exact(i*dt, 1 * dx) * 100);
		}
		fprintf(error, "\n");

		printf("\n\n\n\n");
	}

	fclose(error);
	fclose(actual);
	fclose(estimate);

}


void U0(double u[tstep][xstep])
{
	int i;
	for (i = 0; i<xstep; i++)
	{
		u[0][i] = sin(i*dx*n);
	}
}
void Boundary(double u[tstep][xstep], int tRow)
{
	u[tRow][0] = 0;
	u[tRow][xstep - 1] = 0;
}
void fill(double u[tstep][xstep], int tRow)
{
	int i;
	for (i = 1; i<xstep - 1; i++)
	{
		u[tRow][i] = u[tRow - 1][i] + lam*(u[tRow - 1][i - 1] - 2 * u[tRow - 1][i] + u[tRow - 1][i + 1]) + dt*F(i*dt, i*dx);
	}
}
double exact(double t, double x)
{
	return exp(-n*n*t)*sin(n*x);
}
double F(double t, double x)
{
	return 0;
}
