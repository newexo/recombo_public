#pragma once

#include "mmchain.h"
#include <algorithm>

class recomboFromFile{
private:
	int min_arc, max_arc, n_components;
	char read_mode;
	clkConformationAsList initialComp0, initialComp1;
	clkConformationBfacf3* knot;
	pseudorandom siteSelector;
	ifstream* in;
	ofstream* out;
	void do_recombo_knots();
	void do_recombo_links();
	bool read_comp_knots(ifstream* in);
	bool read_comp_links(ifstream* in);
public:
	recomboFromFile(int Min_arc, int Max_arc, char* Infile, char* Outfile, int n_components, char read_mode);
	void do_recombo();
	~recomboFromFile();
};

