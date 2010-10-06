//  runHydro.cpp
//  Created by Joshua Vredevoogd on 8/30/10.

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <string>
#define PRECISION "%0.6g"

int main(int argc, char * const argv[]){
/*
	double val = 0.07;
	FILE *f = fopen("pow007.dat","w");
	for (int i=0;i<30;i++)
		fprintf(f,"%d %0.9g\n",i,pow(val,i));
*/

	int power = 30;
	FILE *f = fopen("pow30.dat","w");
	for (double val=0.01;val<=0.5;val+=0.01)
		fprintf(f,"%0.3g %0.9g\n",val,pow(val,power));
}