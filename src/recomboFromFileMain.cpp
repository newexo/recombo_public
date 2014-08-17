#include "recomboFromFile.h"
#include <string.h>

void print_usage(){
	cout << "Usage:" << endl;
	cout << "recomboFromFile input_file output_file [operators/options]" << endl;
	cout << "recomboFromFile takes an input file of cubic lattice conformations " << 
		"written in Rob Schaerin's CUBE format and outputs a file of conformations " << 
		"after performing recombination with the user specified criteria." << endl << endl;
	cout << "recomboFromFile operators:" << endl;
	cout << "-min\tminimum arclength" << endl;
	cout << "-max\tmaximum arclength" << endl;
	cout << "-ncomp\tnumber of components" << endl;
	cout << "--m\tinput file format. use --m b for binary (default). use --m t for plain text." << endl;
}

int main(int argc, char* argv[]){
	int min_arc = 0, max_arc = 0, ncomp = 0; 
	char* infile = NULL, *outfile = NULL;
	char read_mode = 0;
	
	if (argc < 6){
		print_usage();
		return 0;
	}
	
	infile = argv[1];
	outfile = argv[2];

	for (int i=3; i < argc; i++){
		if (!strcmp(argv[i], "-min")){
			min_arc = atoi(argv[i+1]);
			i++;
		}
		else if (!strcmp(argv[i], "-max")){
			max_arc = atoi(argv[i+1]);
			i++;
		}
		else if (!strcmp(argv[i], "-ncomp")){
			ncomp = atoi(argv[i+1]);
			i++;
		}
		else if (!strcmp(argv[i], "--m")){
			read_mode = *(argv[i+1]);
			i++;
		}
		else{
			cout << "unrecognized operator/option. Terminating program...";
			return 0;
		}
	}

	//sanity check
	if (min_arc > max_arc){
		cout << "mininum arclength must be strictly <= maximum arclength. Terminating program...";
		return 0;
	}
	//set defauly read mode if none specified
	if (read_mode == 0){
		read_mode = 'b';
	}
	recomboFromFile recombo(min_arc, max_arc, infile, outfile, ncomp, read_mode);
	recombo.do_recombo();
	return 0;
}