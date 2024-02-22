#include <unistd.h>
#include <sstream>

void read_cmd_line(int argc, char *argv[],
		struct Hamiltonian_params& params,
		struct measurments& flags
		/*
		
		double& U, // Interaction strenght
		double& t, // hopping
		double& tp, // tp- t^prime --> the "other" hopping in unit cell
		int& N, // N-number of sites
		bool& Hubbard, // Hubbard model - switch
		bool& calc_1p_spect, // calculate spectrum -switch 
		bool& calc_spin_spect, // calculate spin-spin excitations -switch
		bool& pbc, // use periodic boundary conditions -switch
		bool& two_p, // calculate 2-point correlator -switch
		bool& ss_corr, // calculate spin-spin correlator -switch
		bool& SSH // calculate spin-spin correlator -switch
		*/
		)
{
	char c;
        while ((c = getopt(argc, argv, "t:p:U:N:HASPEV:M:")) != -1) {
		switch (c){
			case 'U':
				params.interaction_U=atof(optarg);
				break;
			case 'V':
				params.el_field=atof(optarg);
				break;
			case 'M':
				params.mag_field=atof(optarg);
				break;
			case 't':
				params.hopping[0]=atof(optarg);
				break;
			case 'p':
				params.hopping.push_back(atof(optarg));
				break;
			case 'N':
				params.num_of_sites=atoi(optarg);
				break;
			case 'H':
				params.Hubbard=true;
				break;
			case 'P':
				params.pbc=true;
				break;
			case 'E':
				flags.two_p=true;
				flags.spin_spin_corr=true;
				break;
			case 'A':
				flags.single_p=true;
				break;
			case 'S':
				flags.spin_spect=true;
				break;
		}
	}
	std::cout<< "Calculating:\n";
	if(flags.single_p){ std::cout << "-> 1p spectra\n";}
	if(flags.spin_spect){std::cout << "-> spin correlation function\n";}	
	if(params.Hubbard){ std::cout <<" For Hubbard model ";}
	else{std::cout << " For Hatsugai-Kohmoto ";}
	if(params.pbc){std::cout<< " Using Periodic Boundary Condition (PBC) ";}
	if(flags.two_p){std::cout << " Calculating 2-point correlation function ";}
	if(flags.spin_spin_corr){std::cout << " Calculating spin-spin correlation function ";}
	std::cout<< " with:\n";
	std::cout<< " N= " << params.num_of_sites;
	std::cout<< " U= " << params.interaction_U;
	int num_of_t{0};
	for(const auto t: params.hopping){
		std::cout<< " t("<< num_of_t<< ")= " << t;
		num_of_t++;
	}
	std::cout << "Magnetic field M=" << params.mag_field;
	std::cout << "Voltage V=" << params.el_field;
	std::cout << std::endl;
}
