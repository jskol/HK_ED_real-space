#include <unistd.h>
#include <sstream>

void read_cmd_line(int argc, char *argv[],
		struct Hamiltonian_params& params,
		struct measurments& flags
		)
{
	char c;
        while ((c = getopt(argc, argv, "t:p:U:N:HASPEV:M:k:m:DR:")) != -1) {
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
			case 'k':
				params.k=atof(optarg);
				params.k_dep=true;
				break;

			case 'm':
				params.model=optarg;
				if(
					params.model == "Kane-Mele" ||
					params.model== "Haldane" ||
					params.model== "Graphene"					
					){
						params.cmplx=true;
					}
				break;

			case 'D':
				flags.electron_density=true;
				break;	
			case 'R':
				flags.retain_states=atoi(optarg);
				break;
		}
	}
	std::cout<< "Calculating "<< params.model <<" model :\n";
	if(flags.single_p){ std::cout << "-> 1p spectra (" << flags.retain_states << " max states kept)\n";}
	if(flags.spin_spect){std::cout << "-> spin correlation function\n";}	
	if(params.Hubbard){ std::cout <<" For Hubbard interaction ";}
	else{std::cout << " For Hatsugai-Kohmoto interaction";}
	if(params.cmplx) {std::cout << " Considering a complex Hamiltonian";}
	if(params.pbc){std::cout<< " Using Periodic Boundary Condition (PBC) ";}
	if(flags.two_p){std::cout << " Calculating 2-point correlation function ";}
	if(flags.spin_spin_corr){std::cout << " Calculating spin-spin correlation function ";}
	if(flags.electron_density){std::cout << " Calculating electron density ";}
	std::cout<< " with:\n";
	std::cout<< " N= " << params.num_of_sites << std::endl;
	std::cout<< " U= " << params.interaction_U << std::endl;
	int num_of_t{0};
	for(const auto t: params.hopping){
		std::cout<< " t("<< num_of_t<< ")= " << t;
		num_of_t++;
	}
	std::cout << std::endl;
	if(params.k_dep){
		std::cout << " k-resolved calculation k=" << params.k <<std::endl;
	}
	
	std::cout << " Magnetic field M=" << params.mag_field<< std::endl;
	std::cout << " Voltage V=" << params.el_field<< std::endl;
	std::cout << std::endl;
}
