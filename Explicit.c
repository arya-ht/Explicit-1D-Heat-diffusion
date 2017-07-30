#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#define n 1
#define D 1.0
#define tMAX 10
#define xMAX M_PI // does not work in VS 2017
#define lam .4 //unstable after 0.5
#define dx (M_PI/10)
#define dt (lam/D*dx*dx)

int xstep = xMAX / dx + 1;
int tstep = (double)tMAX / dt + 1;

double F(double, double);
void U0(double u[tstep][xstep]);
void Boundary(double u[tstep][xstep], int tRow);
void fill(double u[tstep][xstep], int);
void print(double u[tstep][xstep]);
double exact(double, double);
void Explicit();

int main(void) {
	Explicit();
	return 0;
}
void print(double u[tstep][xstep]) {
	int i, j;
	FILE *error = fopen("error.csv", "w"), *estimate = fopen("Exp1D.csv", "w"),
		*actual = fopen("actual.csv", "w"), *t = fopen("t.csv", "w"), *x = fopen("x.csv", "w");
	double error_n = 0.0;
	for (i = 0; i<tstep; i++) {
		printf("Estimate: \n");
		for (j = 0; j<xstep; j++) {
			fprintf(x, "%d", j);      ///x axis
			printf("%lf, ", u[i][j]);
			fprintf(estimate, "%lf ,", u[i][j]);
		}
		fprintf(estimate, "\n");
		printf("\nExact: \n");
		for (j = 0; j<xstep; j++) {
			printf("%lf, ", exact(i*dt, j*dx));
			fprintf(actual, "%lf,", exact(i*dt, j*dx));
		}
		fprintf(actual, "\n");
		printf("\nError: \n");

		for (j = 0; j<xstep; j++) {

			error_n = fabs((exact(i*dt, j*dx) - u[i][j]) / exact(i*dt, j*dx)) * 100;

			if (j == xstep - 1 || j == 0) { // because error becomes infinity (0/0).
				error_n = 0.0;
			}
			printf("%lf%% ", error_n);
			if (j != 0 && j != xstep - 1)
				fprintf(error, "%lf,", fabs((exact(i*dt, 1 * dx) - u[i][1]) / exact(i*dt, 1 * dx)) * 100);
		}
		fprintf(error, "\n");
		fprintf(t, "%d", i); ///y axis
		printf("\n===================================================================================================================\n");
	}

	fclose(error);
	fclose(actual);
	fclose(estimate);
	fclose(t);
	fclose(x);

}
void U0(double u[tstep][xstep]) {
	int i;
	for (i = 0; i<xstep; i++) {
		u[0][i] = sin(i*dx*n);
	}
}
void Boundary(double u[tstep][xstep], int tRow) {
	u[tRow][0] = 0;
	u[tRow][xstep - 1] = 0;
}
void fill(double u[tstep][xstep], int tRow) {
	int i;
	for (i = 1; i<xstep - 1; i++) {
		u[tRow][i] = u[tRow - 1][i] + lam*(u[tRow - 1][i - 1] - 2 * u[tRow - 1][i] + u[tRow - 1][i + 1]) + dt*F(i*dt, i*dx);
	}
}
double exact(double t, double x) {
	return exp(-n*n*t)*sin(n*x);
}
double F(double t, double x) {
	return 0;
}
void Explicit() {
	double U[tstep][xstep];

	U0(U);
	int i;
	for (i = 1; i<tstep; i++) {
		Boundary(U, i);
		fill(U, i);
	}

	print(U);
	if (lam >= 0.5) {
		printf("\n\a Explicit is unstable because Lambda is = %1.5lf", lam); //alert with a beep!
	}
}
