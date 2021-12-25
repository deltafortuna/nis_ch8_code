#include <algorithm>
#include <string>
#include <iterator>
#include <bitset>
#include <random>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <sstream>
#include <regex>

using namespace std;

#include "matrix.h"
#include "params.h"
#include "params.cc" // inclusion causes parameter values to be read
#include "allele.h"
#include "individual.h"
#include "population.h"
#include "metapopulation.h"

int main(int argc, char *argv[]) {
	mt19937 engine(time(0));  //initialize the random engine
	Population::e  = engine;
	mt19937 engine2(time(0) * 2);
	Metapopulation::f = engine2;

	cout << popsize << " is popsize" << endl;
	cout << "MIGRATION matrix:" << endl;
	mig.print_matrix();

	bool simulate = true;

	while(simulate) {
		Metapopulation meta;
		for (int i =0; i < runlength; ++i) {
			if (i % 1 == 0) { cout << "gen " << i << endl;}
			bool test = meta.reproduce_and_migrate(i);
			if (! test) break;
			if (i == runlength-1) simulate=false;
		}
		meta.close_output_files();
	}
	return 0;
}

// static variables
mt19937 Population::e;
mt19937 Metapopulation::f;
