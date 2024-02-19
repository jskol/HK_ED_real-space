#include <unistd.h>
#include <sstream>

void read_cmd_line(int argc, char *argv[],
		double& U, // Interaction strenght
		double& t, // hopping
		double& tp, // tp- t^prime --> the "other" hopping in unit cell
		int& N, // N-number of sites
		bool& Hubbard, // Hubbard model - switch
		bool& calc_1p_spect, // calculate spectrum -switch 
		bool& calc_spin_spect, // calculate spin-spin excitations -switch
		bool& pbc, // use periodic boundary conditions -switch
		bool& two_p, // calculate 2-point correlator -switch
		bool& ss_corr // calculate spin-spin correlator -switch
		bool& SSH // calculate spin-spin correlator -switch
		)
{
	char c;
        while ((c = getopt(argc, argv, "t:p:U:N:HASPE")) != -1) {
		switch (c){
			case 'U':
				U=atof(optarg);
				break;
			case 't':
				t=atof(optarg);
				break;
			case 'p':
				tp=atof(optarg);
				SSH=true;
				break;
			case 'N':
				N=atoi(optarg);
				break;
			case 'H':
				Hubbard=true;
				break;
			case 'P':
				pbc=true;
				break;
			case 'E':
				two_p=true;
				ss_corr=true;
				break;
			case 'A':
				calc_1p_spect=true;
				break;
			case 'S':
				calc_spin_spect=true;
				break;
		}
	}
	/* Safety step to have a simple chian in case tp is not specified */
	if(!SSH){tp=t;}
}
