#ifndef ALLELE_H
#define ALLELE_H

class Allele {

private:
	int position;
	int birthgen;
	int count;
	int originating_population;

public:
	inline int get_count() { return count; }
	inline void set_count(int ccount) { count = ccount; }
	inline int get_position() { return position; }
	inline void set_position(int pposition) { position = pposition; cout << pposition << ": " << position << endl; }
	inline void increment_count() { ++count; }
	inline int get_birthgen() { return birthgen; }
	inline int get_originating_population() { return originating_population; }

	// constructor
	Allele (int pos, int gen, int op): position(pos), birthgen(gen), originating_population(op) {
		count = 0;
	}



};

#endif
