#include <fstream>
#include <vector>

using namespace std;

int main(int argc, char* argv[]){
	srand((unsigned)time(0));
	double MinVals[8] = {0.6, 0.6, 2.4, 0.5, 0.4, 0.03, 0.2, 0.0};
	double MaxVals[8] = {1.1, 1.0, 3.3, 2.0, 1.2, 0.25, 0.8, 100.0};
	
	ofstream outputfile;
	
	outputfile.open("params.dat" , fstream::out);
	
	if(outputfile){
		
	}
	for(int i = 0; i < 5; i++){
		for(int index = 0; index< 7; index++){
			//From stack overflow
			double randfrac = (double)rand()/(double)RAND_MAX;
			double temp = (randfrac * (MaxVals[index]-MinVals[index])) + MinVals[index];
			outputfile << temp << " ";
		}
		outputfile << 1.0 << endl;
	}
}