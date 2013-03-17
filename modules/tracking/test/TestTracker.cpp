#include "interface.h"
#include "Input.h"
#include "Goldtest.h"
#include <string>
#include <stdio.h>
#include <sys/types.h>
#include <dirent.h>

#if(EVALRESULT)
#include "AnalyzeInput.h"
#endif

// usage is ./testtracker <inputdatafile> <outputdir> <output file basename> <number of files>

int main(int argc, char** argv) {
	
	char input_file[100];
	char output_path[100];
	char basename[100];
	int numfiles;

	if (argc > 1) {
		sprintf(input_file, "../../data/%s", argv[1]);
		printf("\nusing input file: %s \n", input_file);
	} 
	else{
		cerr << "No input file provided in call\nProgram exiting TestTracker.cpp)" << endl;
    		exit(1);
	}

	if (argc > 2) {
		sprintf(output_path, "../%s/", argv[2]);
		printf("using output path: %s\n\n", output_path);
	}
	else{
		cerr << "No output folder provided in call\nProgram exiting TestTracker.cpp)" << endl;
    		exit(1);
	}
	if (argc > 3) {
		sprintf(basename, "%s", argv[3]);
		printf("using output file basename: %s\n", basename);
	}
	else{
		cerr << "No base name provided for output files\nProgram exiting TestTracker.cpp)" << endl;
    		exit(1);
	}
	if (argc > 4) {
		numfiles = atoi(argv[4]);
		printf("number of files = %d\n",numfiles);
	}
	else{
		cerr << "No value provided for number of files\nProgram exiting TestTracker.cpp)" << endl;
    		exit(1);
	}
	
	//----Run Final Test----
	Input in;
	vector<Trajectory*> TotalTrajs;  
	
	for (unsigned int i = 0; i < numfiles; i++){
		char iFile[100];
		sprintf(iFile, "%s%s%d", output_path,basename,i);
		printf("%s\n",iFile);
		
		vector<Trajectory*> t = in.trajinput(iFile);

		for(vector<Trajectory*>::iterator it = t.begin(); it != t.end(); ++it)
			TotalTrajs.push_back(*it);
	}

	Goldtest final;
	vector<Trajectory*> GoldTrajs = in.trajinput(input_file);
	final.test(TotalTrajs, GoldTrajs, in.maxframes);	
	
	return 0;
}
