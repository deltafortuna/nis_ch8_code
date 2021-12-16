#ifndef SUMMARYSTATS_H
#define SUMMARYSTATS_H

#include "params.h"

// declarations
double get_pi (vector<bitset<bitlength>> &sample);
double get_watterson (int S);
double get_tajimas_d(double pi, double watterson, int S);
map<string, vector<double> > get_windowStats( vector<int> &positions, vector<bitset<bitlength>> &sample, int S );

// definitions
double get_pi (vector<bitset<bitlength>> &sample) {
	double sumdiffs = 0.;
	double numcomp = 0.;

	#pragma omp parallel for collapse(2)
	for (int i  = 0; i< sample.size() - 1; ++i) {
		for (int j = i+1; j<sample.size(); ++j) {
			sumdiffs += (sample[i] ^ sample[j]).count();
			numcomp+=1.;
		}
	}
	return (sumdiffs/numcomp);
}

double get_watterson (int S) {
	double denominator = 0.;
	for (double i=1.; i<sampsize; ++i)
		denominator += 1./i;
	return (S/denominator);
}

double get_tajimas_d (double pi, double watterson, int S) {
	double d = pi - watterson; // numerator
	double a1 = 0.;
	double a2 = 0.;
	for (double i=1.; i < sampsize; ++i) {
		a1 += 1./i;
		a2 += 1./(i*i);
	}
	double n = sampsize; // for easier expression
	double var = watterson * ( (n+1) / (3*(n-1)) - 1/a1 ) +
	    ( S * (S-1) * ( 1 / (a1*a1+a2)) *
		( (2*(n*n + n+3))/(9*n*(n-1)) - (n+2)/(a1*n) + a2/(a1*a1) ) );
	return(d / sqrt(var));
}

map<string, vector<double> > get_windowStats( vector<int> &positions, vector<bitset<bitlength>> &sample, int S ) {
	map<string, vector<double> >stats;
	map<int, vector<string> > haplotypes;
	map<int, vector<bitset<bitlength> > > windowBitSample;
	vector<double> windowS;
	for (int i = 0; i < sample.size(); ++i) {
		string a = sample[i].to_string();
		string b = a.substr(bitlength-S, S);
		vector<string> haps;
		int count = 0;
		for (int wl = 0; wl+windowSize<= seqlength; wl += windowStep) {
			int wu = wl + windowSize;
			vector<int>::iterator lowindex = upper_bound(positions.begin(), positions.end(), wl);
			vector<int>::iterator upindex = upper_bound(positions.begin(), positions.end(), wu);
			int l = lowindex - positions.begin();
			int u = upindex - positions.begin();
			haps.push_back(b.substr(l, u-l));
			bitset<bitlength> btemp (b.substr(l,u-l));
			windowBitSample[count].push_back(btemp);   // by window #
			count++;
			if (i == 0)
				windowS.push_back(u-l);
		}
		haplotypes[i] = haps;
	}
	stats["S"] = windowS;
	int windowCount = haplotypes[0].size();
	map<int, map<string, int> > haplotype_counts; // key is window number, internal map constists of string and its counts for this window
	for (int i = 0; i<windowCount; ++i) {
		map<string, int> counts_this_window;
		for (int j=0; j<sample.size(); ++j)
			counts_this_window[haplotypes[j][i]]++;
		stats["K"].push_back(counts_this_window.size());
		double pi, wat;
		pi = get_pi(windowBitSample[i]);
		wat = get_watterson(windowS[i]);
		stats["pi"].push_back(pi);
		stats["wat"].push_back(wat);
		stats["tajD"].push_back(get_tajimas_d(pi, wat, windowS[i]));
		haplotype_counts[i] = counts_this_window;    // use for rosenberg Set Statistic
	}
	return (stats);
}

#endif
