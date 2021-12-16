#include <random>
#include <fstream>
#include <map>
#include <vector>
#include <algorithm>
#include "matrix.h"
using namespace std;

vector<int> genos_from_haps(const vector<int> &haps);
vector<int> reproduce(const vector<vector<int> > &parents, const double &r, const vector<double> &prerand);

int main(int argc, char *argv[]) {

	// command line arguments
	double r = atof(argv[1]); // recombination rate r on the range (0.,0.5)
	double s_A = atof(argv[2]); // selection coefficients
	double s_B = atof(argv[3]);
	double t_A = atof(argv[4]);
	double t_B = atof(argv[5]);
	double h0 = atof(argv[6]); // initial haplotype COUNTS of ... A-B
	double h1 = atof(argv[7]); // ... A-b
	double h2 = atof(argv[8]); // ... a-B
	double h3 = atof(argv[9]); // ... a-b
	double eAB = atof(argv[10]); // epistatic parameter for AB/AB
	double eAb = atof(argv[11]); // ... Ab/Ab
	double eaB = atof(argv[12]); // ... aB/aB
	double eab = atof(argv[13]); // ... ab/ab
	int gens = atoi(argv[14]);

	uniform_real_distribution<double> randomnum(0.,1.);
	mt19937 engine(time(0));

	double p_A, p_B, w0, w1, w2, w3, D; // allele frequencies, marginal fitnesses and D

	// calculate popuation size and fitness matrix
	int numhaplotypes = h0+h1+h2+h3;
	int N = numhaplotypes/2;
	vector<double> hapcounts = {h0, h1, h2, h3};
	vector<double> fitness;
	fitness.push_back(1-s_A -s_B-eAB);
	fitness.push_back(1-s_B);
	fitness.push_back(1-t_A-s_B-eAb);
	fitness.push_back(1-s_A);
	fitness.push_back(1);
	fitness.push_back(1-t_A);
	fitness.push_back(1-s_A-t_B-eaB);
	fitness.push_back(1-t_B);
	fitness.push_back(1-t_A-t_B-eab);
	double* f = &fitness[0];
	Matrix<double> fit(3, 3, f); //fitness matrix

	//print fitness matrix
	for (int i=0; i<3; ++i) {
		for (int j=0; j<3; ++j)
			cout << fit[i][j] <<"\t" ;
		cout << endl;
	}

	// create vector of avaialble haplotypes
	vector<int> haplotypes;
	for (int i=0; i<4; ++i)
		for (int j=0; j<hapcounts[i]; ++j)
			haplotypes.push_back(i);

	// randomly pair haplotypes to individuals
	map<int, vector<int> > individuals;
	random_shuffle(haplotypes.begin(), haplotypes.end());
	for (int i=0; i<numhaplotypes; i+=2) {
		individuals[i/2].push_back(haplotypes[i]);
		individuals[i/2].push_back(haplotypes[i+1]);
		individuals[i/2].push_back(0);
		individuals[i/2].push_back(0);
	}

	// add A and B locus genotypes of each individual as [2] and [3] entries
	for (int i=0; i<N; ++i){
		vector<int> haps = {individuals[i][0], individuals[i][1]};
		vector<int> genos = genos_from_haps(haps);
		individuals[i][2] += genos[0];
		individuals[i][3] += genos[1];
	}

	ofstream datafile("twolocus_selection_data");
	datafile << "gen\th1\th2\th3\th4\tp1\tp2\tD"<< endl;

	for (int gen=0; gen<gens; ++gen) {
		// calc and store current haplotype frequencies
		vector<double> haplotypeFreqs;
		for (int i=0; i<4; ++i)
		 	haplotypeFreqs.push_back(hapcounts[i]/numhaplotypes);

		// calc and store current allele frequencies
		p_A = haplotypeFreqs[0] + haplotypeFreqs[1]; // A-B and A-b
		p_B = haplotypeFreqs[0] + haplotypeFreqs[2]; // A-B and a-B

		// calculate D using the A-B haplotype and print
		D = haplotypeFreqs[0] - p_A * p_B;
		datafile << gen << "\t" << haplotypeFreqs[0] << "\t" << haplotypeFreqs[1] << "\t" << haplotypeFreqs[2] << "\t" << haplotypeFreqs[3] << "\t" << p_A << "\t" << p_B << "\t" << D << endl;

		// cull individuals that don't make fitness test
		vector<int> remainder; // stores index of surviving individuals
		for (int i=0; i<N; ++i) {
			if (randomnum(engine) <= fit[ individuals[i][3] ][ individuals[i][2] ] ) //
				remainder.push_back(i);
		}

		uniform_int_distribution<int> randomint(0, remainder.size()-1);
		map<int, vector<int> > nextgen;

		for (int i=0; i<N; ++i) {
			vector<vector<int> > parents;
			parents.push_back( individuals[remainder[randomint(engine)]] );
			parents.push_back( individuals[remainder[randomint(engine)]] );
			vector<double> pr = {randomnum(engine), randomnum(engine), randomnum(engine)};
			nextgen[i] = reproduce(parents, r, pr);
		}

		individuals = nextgen;
		hapcounts = {0,0,0,0};
		for (int i=0; i<N; ++i) {
			hapcounts[individuals[i][0]]++;
			hapcounts[individuals[i][1]]++;
		}
	}
	datafile.close();
	return(0);
}

vector<int> genos_from_haps(const vector<int> &haps) {
	vector<int> vecky = {0,0};
	if (haps[0] == 2 || haps[0] == 3)
		vecky[0]++;
	if (haps[1] == 2 || haps[1] == 3)
		vecky[0]++;
	if (haps[0] == 1 || haps[0] == 3)
		vecky[1]++;
	if (haps[1] == 1 || haps[1] == 3)
		vecky[1]++;
	return(vecky);
}

vector<int> reproduce(const vector<vector<int> > &parents, const double &r, const vector<double> &prerand) {
	vector<int> v = {-1,-1,0,0};
	for (int i=0; i<2; ++i) {
		if (parents[i][0] + parents[i][1] == 3) { // then double heterozygote (AB/ab or Ab/aB)
										    // need to check for recombination ...
			if (prerand[2] <= r) {
				if (parents[i][0] == 0 || parents[i][1] == 0) {
					if (prerand[i] <= 0.5)
						v[i]=1;
					else
						v[i]=2;
				} else {
					if (prerand[i] <= 0.5)
						v[i]=0;
					else
						v[i]=3;
				}
			}
		}
		if(v[i]<0) {
			if (prerand[i] <= 0.5)
				v[i]=parents[i][0];
			else
				v[i]=parents[i][1];
		}
	}

	vector<int> haps = {v[0], v[1]};
	vector<int> genos = genos_from_haps(haps);
	v[2] = genos[0];
	v[3] = genos[1];
	return(v);
}
