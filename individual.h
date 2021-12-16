#ifndef INDIVIDUAL_H
#define INDIVIDUAL_H

#include "allele.h"

class Individual {

private:
	vector<vector<int>> sequences;

	void remove_allele_by_position (int seqnum, int position) {
		auto pos = find(sequences[seqnum].begin(), sequences[seqnum].end(), position);
		if (pos != sequences[seqnum].end())  // ensures position in vector
			sequences[seqnum].erase(pos);
	}

	void resolve_crossover (vector<int>& breaks) {
		int numbreaks = breaks.size();
		vector<int>::iterator lower, upper;
		map<int, vector<int> > segments;
		vector<int> newvec;

		for (int seq = 0; seq < 2; ++seq) {
			lower = sequences[seq].begin();
			for (int i=0; i<numbreaks; ++i) {
				upper  = upper_bound(sequences[seq].begin(), sequences[seq].end(), breaks[i]);
				newvec.assign(lower, upper);
				lower = upper;
				segments[i + seq*(numbreaks+1)] = newvec;
			}
			newvec.assign(lower, sequences[seq].end());
			segments[numbreaks + seq*(numbreaks+1)] = newvec;
		}

		for(int i=0; i<= numbreaks; ++i) {
			if (i%2 == 0) {
				if (i == 0) {sequences[0] = segments[0]; sequences[1] = segments[numbreaks+1];}
				else {
					sequences[0].insert(sequences[0].end(), segments[i].begin(), segments[i].end());
					sequences[1].insert(sequences[1].end(), segments[i+numbreaks+1].begin(), segments[i+numbreaks+1].end());
				}
			} else{
				sequences[0].insert(sequences[0].end(), segments[i+numbreaks+1].begin(), segments[i+numbreaks+1].end());
				sequences[1].insert(sequences[1].end(), segments[i].begin(), segments[i].end());
			}
		}
	}

public:
	inline vector<int> get_sequence(int whichseq) { return sequences[whichseq]; }
	inline vector<vector<int> > get_sequences() { return sequences; }

	void remove_fixed_allele(int to_remove) {
		for (int i = 0; i<2; ++i) {
			vector<int>::iterator p = find(sequences[i].begin(), sequences[i].end(), to_remove);
			if (p != sequences[i].end()) // i.e., element not found
		    		sequences[i].erase(p);
		 }
	}

/*	vector<int> get_alleles() {
		vector<int> a = sequences[0];
		a.insert(a.end(), sequences[1].begin(), sequences[1].end());
		return(a);
	}
*/

	// generation 0 constructor
	Individual (vector<vector<int>> seqs): sequences(seqs) {
			;
	}

	// intra-simulation constructor
	Individual (Individual *p1, Individual *p2, vector<vector<int> > mutation_results, vector<int> breakpoints) {
		sequences.push_back((*p1).get_sequence(mutation_results[2][0]));
		sequences.push_back((*p2).get_sequence(mutation_results[3][0]));
		for (int i=0; i<2; ++i) {
			for (int j=1; j<mutation_results[i].size(); ++j)  {
				if (mutation_results[i][j] > 0) {
					sequences[i].push_back(mutation_results[i][j]);
					sort(sequences[i].begin(), sequences[i].end());
				} else
					remove_allele_by_position(i, -1 * mutation_results[i][j]);
			}
		}
		if (breakpoints.size() > 0)
			resolve_crossover(breakpoints);
	}

	~Individual() {} // destructor
};

#endif
