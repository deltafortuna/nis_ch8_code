#ifndef METAPOPULATION_H
#define METAPOPULATION_H

class Metapopulation {

private:
vector<Population*> populations;
uniform_real_distribution<double> random01;

public:
	int reproduce_and_migrate(int gen) {
		// reproduction within all extant demes
        // AND check for birth/extinction of population
		bool fixtest = true;
		for (int i=0; i<pop_num; ++i) {
			if ((*populations[i]).get_extant()) {
				if (extinctgen[i] == gen)
					(*populations[i]).set_extinct();
				else
					(*populations[i]).reproduce(gen);
			}	else {
				if ( birthgen[i] == gen ) {
					vector<int> change = (*populations[i]).set_extant();
					map<int, int> alleles_to_add;
					if (change[0] == 1) { //split
						int N = (double) (change[2])  / 100 * (*populations[ change[1] ]).get_current_popsize(gen);
						for (int k = 0; k<N; ++k) {
							vector<vector<int>> v1 = (*populations[ change[1] ]).get_sequences(k);
							(*populations[i]).add_immigrant(  v1  );
							for (int m=0; m<2; ++m)
								for (auto iter=v1[m].begin(); iter!=v1[m].end(); ++iter)
									++alleles_to_add[*iter];
						}
						(*populations[ change[1] ]).remove_emigrants(N); // so that the other deme populated by the split doesn't receive same individuals
					}	else if (change[0] == 2) { //merger
						vector<int> N;
						N.push_back( (*populations[ change[1] ]).get_current_popsize(gen) );
						N.push_back( (*populations[ change[2] ]).get_current_popsize(gen) );
						for (auto q:N) {
							for (int k = 0; k<q; ++k) {
								vector<vector<int> > v1 = (*populations[ change[1] ]).get_sequences(k);
								(*populations[i]).add_immigrant( v1 );
								for (int m=0; m<2; ++m)
									for (auto iter=v1[m].begin(); iter!=v1[m].end(); ++iter)
										++alleles_to_add[*iter];
							}
						}
					}	 // else built from MS
					if (change[0] > 0) // deme not built from MS
							for (auto iter=alleles_to_add.begin(); iter != alleles_to_add.end(); ++iter) // populate alleles in new deme
								(*populations[i]).insert_new_allele( (*populations[ change[1] ]).get_allele_info(iter->first) );
					(*populations[i]).reproduce(gen);
				}
			}
		}

		// migration among all demes
		for (int i=0; i<pop_num; ++i) {
			if ( (*populations[i]).get_extant() ) {
				for (int j=0; j < pop_num; ++j) {
					if ( (*populations[j]).get_extant() ) {
						double Nm = mig[i][j] * pop_schedule[j][gen];
						if (Nm < 1) {
							if(random01(f) < Nm)
								Nm = 1;
							else
								Nm = 0;
						}
						else
							Nm = floor(Nm);

						for (int k=0; k<Nm; ++k) {
							vector<vector<int>> v1 = (*populations[i]).get_sequences(k);
							(*populations[j]).add_immigrant(  v1  );

							// check for new alleles introduced to population j
							for (int m = 0; m < 2; ++m) {
								vector<int> v2 = (*populations[j]).get_allele_positions();
								vector<int> diff;
								set_difference(v1[m].begin(), v1[m].end(), v2.begin(), v2.end(), inserter(diff, diff.begin()));
								for (auto q:diff)
										(*populations[j]).insert_new_allele( (*populations[i]).get_allele_info(q) ) ;
							}
						}
						(*populations[i]).remove_emigrants(Nm);
					} else continue;
				}
			} else continue;
		}
		if ( gen==0 || (gen+1) % sampfreq == 0) {
			for (int i=0; i<pop_num; ++i) {
				if ((*populations[i]).get_extant()) {
					if (possel[(*populations[i]).get_popnum()][0] != 0 && sellocus[(*populations[i]).get_popnum()][4] == 1)
						fixtest = (*populations[i]).sample(gen);
					else
						int dummy = (*populations[i]).sample(gen); //still update alleles
					if (! fixtest )
							break;
				}
			}
		}
		return(fixtest);
	}

	void close_output_files() {
		for (auto iter = populations.begin(); iter != populations.end(); ++iter)
			(*iter)->close_output_files();
	}

	Metapopulation()  {
		for (int i=0; i < pop_num; ++i) {
			if (birthgen[i] != 0 )
				populations.push_back( new Population(i, 0) );
			else
				populations.push_back( new Population(i, 1) );
		}
		random01.param(uniform_real_distribution<double>::param_type(0.,1.));
	}

	~Metapopulation() {}

	static mt19937 f;
};

#endif
