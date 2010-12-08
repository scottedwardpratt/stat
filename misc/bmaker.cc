#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <complex>

const double PI=3.14159265358979323844;
const double HBARC=197.3269602;

using namespace std;

int main(){
	int NB=10,ib;
	double bb2,oldbb2,b;
	//double bb2max=283.296; // = sigma_tot/pi, sigmatot=8.9 barns
	double bb2max=216.45; // = sigma_tot/pi, sigmatot=6.8 barns
	oldbb2=0.0;
	for(ib=0;ib<NB;ib++){
		bb2=(double(ib+1)/double(NB))*bb2max;
		b=sqrt(0.5*(bb2+oldbb2));
		oldbb2=bb2;
		printf("%3.0f-%3.0f%%: %g\n",100.0*ib/double(NB),100.0*(ib+1)/double(NB),b);
	}

  return 0;
}


