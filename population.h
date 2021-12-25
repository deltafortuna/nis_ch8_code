#ifndef POPULATION_H
#define POPULATION_H

#include "params.h"  // provides access to global parameter values (extern)
#include "summarystats.h" // provides access to  summary statistic calculations/functions

class Population {

private:
double mu_sequence;
vector<Individual*> individuals;
map<int, Allele*> alleles; /// int key is position of the allele
uniform_int_distribution<int> randompos;
uniform_int_distribution<int> randomind;
uniform_real_distribution<double> randomnum;
poisson_distribution<int> randommut;
ofstream allele_file;
ofstream sumstat_file; // all sumstats printed here
double r_sequence;
poisson_distribution<int> randomrec;
int popn;
bool extant;
ofstream abf;
ofstream sinfo;
bool activeselection{}, standingvar{}, newvariant{}, nfdselection{}, negselection{};
int keypos=-999;
vector<double> fitness;
double nfd_s, nfd_h;
vector<int> purifying_sites;

void update_selected_freqAndFit(const int &gen) {
	double count = pop_schedule[popn][gen]*2; // start count at ALL alleles derived
	vector<double> genotype_freqs = {0.,0.,0.}; // holds frequency of A/A, A/a, a/a, where A is derived allele
	for (auto iter = individuals.begin(); iter != individuals.end(); ++iter) {
		int num_ancestral_alleles = (**iter).get_genotype(keypos);
		count -= num_ancestral_alleles;
		genotype_freqs[num_ancestral_alleles]++;
	}
	double p = count/(pop_schedule[popn][gen]*2); // frequency of derived allele
	double q = 1-p;
	for (int i=0; i<3; ++i)
		genotype_freqs[i] /= pop_schedule[popn][gen];
	sinfo << gen << "\t" << p << "\t" << genotype_freqs[0] << "\t" << genotype_freqs[1] << "\t" << genotype_freqs[2] << "\t";

	if (nfdselection) {
		/* fitness[0] = 1 - (nfd_s*p*p);
		fitness[1] = 1 - (2*nfd_s*nfd_h*p*q);
		fitness[2] = 1 - (nfd_s*q*q);*/

		double max_fit = -999;
		fitness[0] = 1 - nfd_h*genotype_freqs[1] + nfd_h*genotype_freqs[2];
		if (fitness[0] > max_fit) max_fit = fitness[0];
		fitness[1] = 1 - nfd_s*genotype_freqs[1];
		if (fitness[1] > max_fit) max_fit = fitness[1];
		fitness[2] = 1 - nfd_h*genotype_freqs[1] + nfd_h*genotype_freqs[0];
		if (fitness[2] > max_fit) max_fit = fitness[2];
		for (int i=0; i<3; ++i)
			fitness[i] /= max_fit;
	}
	sinfo << fitness[0] << "\t" << fitness[1] << "\t" << fitness[2] << endl;
}

bool update_alleles(const int &gen) {
	bool fixtest = true;
		// reset all counts to zero
	for (auto iter = alleles.begin(); iter != alleles.end(); ++iter)
		(*(iter->second)).set_count(0);
	map<int, int> new_allele_counts;
	for (auto iter = individuals.begin(); iter != individuals.end(); ++iter) {
		for (int i=0; i<2; ++i) {
			const vector<int> &a  = (**iter).get_seq(i);
			for (auto iter2 = a.begin(); iter2 != a.end(); ++iter2) {
				++new_allele_counts[*iter2];
			}
		}
	}

	for (auto iter3 = new_allele_counts.begin(); iter3 != new_allele_counts.end(); ++iter3)
		(*alleles[iter3->first]).set_count(iter3->second);

	// identify lost and fixed alleles and print to allele history file
	vector<int> to_remove;
	for (auto iter = alleles.begin(); iter != alleles.end(); ++iter) {
		int current_count = (*(iter->second)).get_count();
		if (current_count == 0) { // allele LOST from population
			to_remove.push_back(iter->first);  // first is position
			if (iter->first == keypos) fixtest = false;
			/*int birthgen = (*(iter->second)).get_birthgen();
			allele_file << iter->first << " " << birthgen << " " << gen - birthgen << " 0" << endl;*/
		}
		if (current_count == pop_schedule[popn][gen]*2) {  // then, derived allele FIXED
			if ( (activeselection && iter->first != keypos) || !activeselection) {
				to_remove.push_back(iter->first);
			//int birthgen = (*(iter->second)).get_birthgen();  // uncomment if using trackAlleleBirths
				for (auto iter2 = individuals.begin(); iter2 != individuals.end(); ++iter2)
					(**iter2).remove_fixed_allele(iter->first); // removed fixed allele's position from all individuals' sequences (currently stored in nextgen)
		//		allele_file << iter->first << " " << birthgen << " " << gen - birthgen << " 1" << endl;
			}
		}
	}

	// free memory associated w/ lost/fixed alleles and remove entry from alleles container
	for (auto iter = to_remove.begin(); iter != to_remove.end(); ++iter) {
		delete alleles[*iter];  // free memory from Allele object itself
		alleles.erase(*iter);  // erase alleles map entry corresponding to the deleted allele
	}
	return(fixtest);
}

void get_sample(int gen) {
	ofstream sequencefile;  //only used if gen % printhapfreq == 0
	string ofname = "deme" + to_string(popn) + "_" + to_string(gen);
	bool printhap = false;
	if (gen % printhapfreq == 0) printhap = true;
	if (printhap) sequencefile.open(ofname.c_str());
	vector<bitset<bitlength>> sample;
	map<int, int> allele_counts;
	int count = 0;
	int additional = sampsize;
	if (diploid_sample)
		additional /= 2;

	for (auto iter = individuals.begin(); iter != individuals.begin()+sampsize; ++iter) { // determines which alleles are present in sample
		vector<int> haplotype = (**iter).get_sequence(0);
		for (auto iter2 = haplotype.begin(); iter2 != haplotype.end(); ++iter2)
			++allele_counts[*iter2];
		if (diploid_sample) {
			vector<int> haplotype = (**iter).get_sequence(1);
			for (auto iter2 = haplotype.begin(); iter2 != haplotype.end(); ++iter2)
				++allele_counts[*iter2];
		}
	}

	vector<int> positions;
	for (auto iter = allele_counts.begin(); iter != allele_counts.end(); ++iter)
		positions.push_back(iter->first);

	if (printhap) { // print column headers
		for (auto iter = positions.begin(); iter != positions.end(); ++iter)
			sequencefile << "nt" << to_string(*iter) << " ";
		sequencefile << endl;
		for (auto iter = positions.begin(); iter != positions.end(); ++iter)
			sequencefile << (alleles[*iter]) -> get_originating_population() << " ";
		sequencefile << endl;
		for (auto iter = positions.begin(); iter != positions.end(); ++iter)
			sequencefile << (alleles[*iter]) -> get_birthgen() << " ";
		sequencefile << endl;
	}

	for (auto iter = individuals.begin(); iter != individuals.begin()+sampsize; ++iter) {  // creates haplotypes for each sequence in the sample and populates bitset
		for (int h=0; h<diploid_sample+1; ++h) {
			vector<int> haplotype = (**iter).get_sequence(h);
			sort(haplotype.begin(), haplotype.end());
			string hap;
			for (auto iter = allele_counts.begin(); iter != allele_counts.end(); ++iter)
				if ( binary_search (haplotype.begin(), haplotype.end(), iter->first)) {
					hap += "1";
					if (printhap) sequencefile << "1 ";
				} else {
					hap += "0";
					if (printhap) sequencefile << "0 ";
				}
			sample.push_back(bitset<bitlength> (hap));
			if (printhap) sequencefile << endl;
		}
	}

	int S = allele_counts.size();

	if (getWindowStats) {
		map<string, vector<double> > stats = get_windowStats(positions, sample, S);
		for (auto iter = stats.begin(); iter != stats.end(); ++iter) {
			sumstat_file << gen << " ";
			sumstat_file << iter->first << " ";
			for (auto iter2 = (iter->second).begin(); iter2 != (iter->second).end(); ++iter2)
				sumstat_file << *iter2 << " ";
			sumstat_file << endl;
		}
	} else {
		double pi = get_pi(sample);
		double watterson = get_watterson(S);
		double tajimasd = get_tajimas_d(pi, watterson, S);
		sumstat_file << gen << " " << pi << " " << watterson << " " << tajimasd << endl;
	}

	if (printhap) sequencefile.close();
}

vector<vector<int> > mutate(const vector<int> &parents, const int &gen) {
	vector<vector<int> > mutation_results;

	// determine which, if any, positions are mutated
	vector<int> mutnum{randommut(e)};
	mutnum.push_back(randommut(e));

	mutation_results.push_back({mutnum[0]});
	mutation_results.push_back({mutnum[1]});

	// determine which of the two homologs is transmitted by each parent
	mutation_results.push_back({(randomnum(e)<0.5) ? 0 : 1});
	mutation_results.push_back({(randomnum(e)<0.5) ? 0 : 1});

	// resolve any mutation(s) that did occur
	for (int i=0; i<2; ++i)  {
		for (int j = 0; j < mutnum[i]; ++j)  { // loop not entered if no mutation (i.e., mutnum[i] == 0)
			int position = randompos(e);
			if (alleles.find(position) == alleles.end()) { // new mutation to a derived allele in the population
				alleles.insert({position, new Allele(position, gen, popn)});
				if (trackAlleleBirths)  // new statement
					abf << "nt" << position << "\t" << gen << "\t" << popn << endl;
				mutation_results[i].push_back(position);
			} else { // mutation present in POPULATION; determine if derived allele found in the considered sequence
				vector<int> seq = (*(individuals[parents[i]])).get_sequence(mutation_results[i+2][0]);
				vector<int>::iterator p = find(seq.begin(), seq.end(), position);
				if (p != seq.end())   // back mutation
					mutation_results[i].push_back(position * -1);	 // negative position signals removal of allele by back mutation
				else
					mutation_results[i].push_back(position);
			}
 		}
	}
	return mutation_results;
}

vector<int> recombine() {
 vector<int> breakpoints;
 if (useRec) {  /// if false, empty vector passed to Individual constructor
	 int chiasmata = randomrec(e);
	 for (int i=0; i<chiasmata; ++i) {
		 if (useHotRec) {
			 double c = randomnum(e);
			 if (c < hotrecStart * recrate / r_sequence )
				 breakpoints.push_back( c * r_sequence / recrate );
			 else if (c < (hotrecStop*hotrecrate - hotrecStart*(hotrecrate-recrate)) / r_sequence)
				 breakpoints.push_back(  (c*r_sequence + hotrecStart*(hotrecrate - recrate))    /     hotrecrate);
			 else
				 breakpoints.push_back( (c*r_sequence - (hotrecStop-hotrecStart)*(hotrecrate-recrate)   )    / recrate );
		 } else {
			 breakpoints.push_back(randompos(e));
		 }
	 }
	 sort(breakpoints.begin(), breakpoints.end());
 }
 return breakpoints;
}



public:
void reproduce(int gen) {
	randomind.param(uniform_int_distribution<int>::param_type(0,pop_schedule[popn][gen]-1));
	int N = pop_schedule[popn][gen]; // N is number of individuals to produce
	vector<int> recode_parents;

	if (modelMigration) {
		int n_imm = 0;
		int n_emi = 0;
		for (int i=0; i<pop_num; ++i) {
			n_emi += mig[popn][i] * pop_schedule[i][gen];
			n_imm += mig[i][popn] * N;
		}
		N -= n_emi;
		N += n_imm;
	}

	if (activeselection)
		for (int i=0; i<individuals.size(); ++i)
			if(randomnum(e) <= fitness[(*individuals[i]).get_genotype(keypos)])
				recode_parents.push_back(i);
	if (negselection)
		for (int i=0; i<individuals.size(); ++i)
			if (randomnum(e) <= fitness[(*individuals[i]).get_negsel_genotype(purifying_sites)])
				recode_parents.push_back(i);
	if (activeselection || negselection)
		randomind.param(uniform_int_distribution<int>::param_type(0, recode_parents.size() - 1));

	#pragma omp parallel for
	for (int i=0; i< N; ++i) {
		vector<int> parents;
		if (activeselection || negselection) {
			parents.push_back(recode_parents[randomind(e)]);
			parents.push_back(recode_parents[randomind(e)]);
		} else {
			parents.push_back(randomind(e));
			parents.push_back(randomind(e));
		}
		individuals.push_back(new Individual(individuals[parents[0]], individuals[parents[1]], mutate(parents, gen), recombine()) );
	}

	// delete dynamically allocated individuasl of the last generation
	for (auto iter = individuals.begin(); iter != individuals.end() - N; ++iter)

		delete *iter;
	// remove orphaned pointers from individuals
	individuals.erase(individuals.begin(), individuals.end() - N);


	// update allele counts on sample generations
	if (gen % sampfreq == 0 && (!activeselection && !negselection))
		update_alleles(gen);

	if ( gen != 0 && gen % sampfreq == 0 ) {
		random_shuffle(individuals.begin(), individuals.end() ) ;
		get_sample(gen);
	}
	if (activeselection) update_selected_freqAndFit(gen);
}

void remove_emigrants(int Nm) {
	for (auto iter = individuals.begin(); iter != individuals.begin() + Nm; ++iter)
		delete *iter;
	individuals.erase(individuals.begin(), individuals.begin()+Nm);
}

vector<int> get_allele_positions() {
	vector<int> v;
	for (auto iter=alleles.begin(); iter!=alleles.end(); ++iter)
		v.push_back(iter->first);
	return v;
}

vector<int> get_allele_info(int s) {
	vector<int> v = {s};
	v.push_back( (*alleles[s]).get_birthgen() );
	v.push_back( (*alleles[s]).get_originating_population() );
	return v;
}

void insert_new_allele(vector<int> v) {
	alleles.insert( { v[0]  , new Allele( v[0], v[1], v[2] ) } );
}

bool sample(int gen) {
	bool fixtest = update_alleles(gen);
	if (!fixtest && sellocus[popn][4]==1) {return fixtest;}
	random_shuffle(individuals.begin(), individuals.end() ) ;
	get_sample(gen);
	return(fixtest);
}

void close_output_files () {
	allele_file.close();
	sumstat_file.close();
	if (trackAlleleBirths) abf.close();
	if (activeselection) sinfo.close();
}

vector<int> set_extant() {
	extant = 1;
	vector<int> i;
	if (splitgenesis[popn][0] > 0) {
		i.push_back(1);
		i.push_back(splitgenesis[popn][1]); // source population
		i.push_back(splitgenesis[popn][2]); // percent
	} else if (mergegenesis[popn][0] > 0 ) {
		i.push_back(2);
		i.push_back(mergegenesis[popn][1]); // first source population
		i.push_back(mergegenesis[popn][2]); // second source population
	} else
		i.push_back(0);
	return(i);
}

inline int get_popnum() { return popn;}
inline bool get_extant() {return extant;}
inline void set_extinct() {extant = 0;}
inline vector<vector<int>> get_sequences(int indnum) { return (*individuals[indnum]).get_sequences();}
inline void add_immigrant(vector<vector<int>> ses) {individuals.push_back( new Individual(ses) ); }
inline int get_current_popsize(int gen) {return pop_schedule[popn][gen];}

Population (int popnum, int eextant):popn(popnum), extant(eextant) { // added parameters to constructor
	// initialize random number distributions
	mu_sequence = seqlength * mutrate;
	randompos.param(uniform_int_distribution<int>::param_type(1,seqlength));
	randomind.param(uniform_int_distribution<int>::param_type(0,pop_schedule[popn][0] - 1)); // replaces ch3 line
	randomnum.param(uniform_real_distribution<double>::param_type(0.,1.));
 	randommut.param(poisson_distribution<int>::param_type(mu_sequence));

	individuals.reserve(popsize*10);

	if (possel[popn][0]!=0 || nfdsel[popn][0]!=0 || negsel[popn][0]!=0) {
			if (nfdsel[popn][0] != 0)
				nfdselection = true;
        if (negsel[popn][0] != 0) {
      		negselection = true;
          fitness = {1, 1-(negsel[popn][0]*negsel[popn][1]), 1-negsel[popn][0]};
          vector<int> possible_sites(negsel[popn][3]-negsel[popn][2]+1);
					iota(possible_sites.begin(), possible_sites.end(), negsel[popn][2]);
					random_shuffle(possible_sites.begin(), possible_sites.end());
					auto iter = possible_sites.begin();
					purifying_sites.assign(iter, iter+negsel[popn][4]);
					sort (purifying_sites.begin(), purifying_sites.end());
        } else
        	activeselection = true;

        if (activeselection) {
					keypos = sellocus[popn][0];
					if (sellocus[popn][1] > 1)
	        			standingvar=true;
					else
						newvariant=true;
					double s = possel[popn][0];
					double t = possel[popn][1];
					double h = possel[popn][2];
					nfd_s = nfdsel[popn][0];
					nfd_h = nfdsel[popn][1];
					if (t == 0) fitness = {1, 1-(h*s), 1-s}; // dominant or additive
					else fitness = {1-s, 1, 1-t}; // overdominant
        }
	}

	if (useRec) {
		if (useHotRec) {
			int hotspot_length = hotrecStop - hotrecStart + 1;
			r_sequence = ( hotspot_length * hotrecrate )  +  ((seqlength - hotspot_length) * recrate);
		} else
			r_sequence = seqlength * recrate;
		randomrec.param(poisson_distribution<int>::param_type(r_sequence));
	}

	vector<int> allele_positions;
	if (activeselection && !standingvar)
		alleles.insert( { keypos, new Allele(keypos, -1, popn) } );

	if (useMS[popn]) { // start population with MS generated variation
		cout << "using MS to initialize population " << popn << endl;
		system(mscommand[popn].c_str());
		ifstream ms_output("ms_output");
		string ms_line;
		regex query("positions");
		bool trigger = false;

		while(getline(ms_output, ms_line)) {
			if (regex_search(ms_line, query)) {
				trigger = true;
				istringstream iss(ms_line);
				string s;
				iss >> s; //skip the first subpart, which is "positions:"
				while (iss >> s) { // read decimal positions,
	                                   // convert to base pair position,
									// and create new allele at that position
					int position = seqlength * atof(s.c_str());
					allele_positions.push_back(position);
					if (standingvar || position != keypos)
						alleles.insert( { position  , new Allele(position,-1,popn) } );
				}
				continue;
			}

			if (trigger) { // allele positions determined;
				vector<int> s1, s2;
				for (int i=0; i < ms_line.length(); ++i)
					if (ms_line[i] == '1') {
						if (allele_positions[i] != keypos)
							s1.push_back(allele_positions[i]);
						if (standingvar)
							(*(alleles[allele_positions[i]])).increment_count();
					}
				getline(ms_output, ms_line);
				for (int i=0; i < ms_line.length(); ++i)
					if (ms_line[i] == '1') {
						if (allele_positions[i] != keypos)
							s2.push_back(allele_positions[i]);
						if (standingvar)
							(*(alleles[allele_positions[i]])).increment_count();
					}
				if (activeselection && !standingvar && newvariant) {
					s1.push_back(keypos);
					newvariant = false;
					sort(s1.begin(), s1.end());
				}
				vector<vector<int>> ses{s1,s2};
				individuals.push_back( new Individual(ses) );
			}
		}
	} else {
		for (int i=0; i<pop_schedule[popn][0]; ++i) { // replaces ch3 line
			vector<int> s1; vector<int> s2;
			vector<vector<int>> ses{s1,s2};
			individuals.push_back(  new Individual(ses)  );
		}
	}

	string ofname = "deme" + to_string(popn) + "_allele_births";
	if (trackAlleleBirths)
		abf.open(ofname.c_str());
	if (activeselection) {
		ofname = "deme" + to_string(popn) + "_selectioninfo";
		sinfo.open(ofname.c_str());
	}

	string fname = "allele_info";
	allele_file.open(fname.c_str());
	allele_file << "position birthgen lifespan extinct.fixed" << endl;

	fname = "sumstats" + to_string(popn);
	sumstat_file.open(fname.c_str());
	if (getWindowStats) {
		sumstat_file << "gen stat ";
		for (int w=0; w + windowSize <= seqlength; w += windowStep)
			sumstat_file << "w" << w+windowStep << " ";
		sumstat_file << endl;
	}	else
		sumstat_file << "gen pi watterson tajimasd" << endl;

	if (activeselection && standingvar) {   // determine keypos
		int low=keypos-floor(sellocus[popn][2]*seqlength);
		if (low < 0 ) low = 0;
		int high=keypos+floor(sellocus[popn][2]*seqlength);
		if (high >= seqlength) high = seqlength-1;
		vector<int> current_best{-999, 999999,-999};
		int acceptablediff = floor(pop_schedule[popn][0]*2*sellocus[popn][3]);
		int target_low = sellocus[popn][1] - acceptablediff;
		int target_high = sellocus[popn][1] + acceptablediff;
		for (auto iter=alleles.begin(); iter != alleles.end(); ++iter) {
			if (iter->first >= low && iter->first <= high) { // in range
				int ct = (*(iter->second)).get_count();
				if (ct >= target_low && ct <=target_high) {
					int diff = abs(sellocus[popn][1]-ct);
					if (diff < current_best[1] ) {
						current_best[0] = iter->first;
						current_best[1] = diff;
						current_best[2] = ct;
					}
				}
			} else {
				if (iter->first > high) break;
				else continue;
			}
		}
		if (current_best[0] != -999) {
			keypos = current_best[0];
			sinfo << "Standing variant identified at base pair " << keypos << " with a starting allele count at time 0 of " << current_best[2] << endl;
		} else {
			cout << "Did not find a suitable standing variant." << endl;
			throw("suitable standing variant not identified");
		}

		if (nfdselection)
 			update_selected_freqAndFit(0);
	}
	sinfo << "gen\tp\tf_AA\tf_Aa\tf_aa\tw_AA\tw_Aa\tw_aa" << endl;
}

static mt19937 e;

};

#endif
