#include <random>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <iostream>

using namespace std;

int main(int argc, char *argv[]) {

	// command line arguments
	double N = atof(argv[1]);
	string suffix = argv[2];
	double s = atof(argv[3]);
	double h = atof(argv[4]);

	map<int, int> pop;
	vector<double> fitness;
	fitness.push_back(1.-s);
	fitness.push_back(1.-(h*s));
	fitness.push_back(1.);

	default_random_engine engine(time(0));  //initialize the random engine
	map<int, vector<double>> data;
	uniform_int_distribution<int> randind(0,N-1);
	uniform_real_distribution<double> randd(0.,1.);

	bool trigger = true;
	int gen;
// int numtries = 0;    

	while(trigger) {
// numtries++;
		// populate the population with N-1 homozygotes for ancestral allele
		// ... and 1 heterozygote for the derived, beneficial allele
		for (int i=0; i<N-1; i++)
			pop[i] = 0;
		pop[N-1] = 1;
		gen = 0;
		double p = 1/(2*N);

		while (p != 0. && p != 1.) {
			gen++;
			vector<int> survivor_indices;
			map<int, int> nextpop;
			for (int i=0; i<N; i++) {
				if (randd(engine) <= fitness[pop[i]])
					survivor_indices.push_back(i);
			}
			uniform_int_distribution<int> randind(0,survivor_indices.size());

			for (int i=0; i<N; i++) { // constant population size
				double offspring = 0.;
				for (int j=0; j<2; j++) {
					switch (pop[survivor_indices[ randind(engine) ] ] ) {
						case 0: break;
						case 1: if (randd(engine) <= 0.5) { offspring+=1.; }
								break;
						case 2: offspring+=1.;
					}
				}
				nextpop[i] = offspring;
			}
			pop = nextpop;
			p = 0.;
			for (int i = 0; i<N; i++) p += pop[i];
			p /= (2*N);
			data[gen].push_back(p);
			double q = 1. - p;
			double popfit = p*(p*fitness[2]+q*fitness[1]) + q*(p*fitness[1]+q*fitness[0]); // mean population fitness
			data[gen].push_back(popfit);
// cout << p << endl;
		}
		if (p == 1.) trigger = false;
		else data.clear();
// cout << "****************************" << endl;
	}

// cout << "Number of tries: " << numtries << endl;

	// output data held in the map, data
	string fname = "stochastic_selection_resultsN" + suffix;
	ofstream output;
	output.open(fname.c_str());
	output << "gen\tfrequency\tpopfit" << endl;
	for (int i=1; i<data.size(); i++)
		output << i << "\t" << data[i][0] << "\t" << data[i][1] << endl;
	output.close();

	return 0;
}
