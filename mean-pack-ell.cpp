#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <cstring>
using namespace std;
#include "random.h"
#include "ellipsoid.h"

int main (int argc, char *argv[]) {
	int No_parts, No_species, n = 0;
	double Din, Pactual, pack = 0;
	vector Box;
	cout.setf(ios::fixed, ios::floatfield);
	while (--argc > 1) {
		ifstream pResltx (*++argv);
		if (!pResltx) {
			cout << "Cannot open file " << *argv << endl;
			exit(1);
		}
		pResltx.read((char *)&No_parts, sizeof(int));
		pResltx.read((char *)&Box, sizeof(vector));
		pResltx.read((char *)&Pactual, sizeof(double));
		pResltx.read((char *)&Din, sizeof(double));
		pResltx.read((char *)&No_species, sizeof(int));
		n++;
		pack += Pactual;
		cout << *argv << " - " << Pactual << endl;
		pResltx.close();
	}
	cout << "--------------------------------------" << endl;
	cout << pack / n << endl;
	ofstream myfile;
  	myfile.open (*++argv);
  	myfile << pack / n;
  	myfile.close();
	return 0;
}


