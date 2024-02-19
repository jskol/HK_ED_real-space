#include <cstdlib>
#include <iostream>
#include <complex>
#include <vector>
#include <array>

#include <cassert>
#include <iostream>
#include <fstream>

#include <chrono>


typedef std::complex<double> cx_double;
typedef std::vector<cx_double> cx_vector;
std::array<std::string,2> spins_set{"up","do"};

#ifndef EIGEN_STUFF
#define EIGEN_STUFF
#include <Eigen/Dense>
#endif

// Globally define fermionic operators
#ifndef LIBCOMMUTE_STUFF
#define LIBCOMMUTE_STUFF
#include <libcommute/loperator/state_vector_eigen3.hpp>
#include <libcommute/libcommute.hpp>
using libcommute::static_indices::c_dag; // Create an electron
using libcommute::static_indices::c;     // Destroy an electron
using libcommute::static_indices::n;     // number operator
#endif

#include "read_params.hpp"
#include "gen_SSH_H.hpp"
#include "Gen_GF.hpp"
#include "spin_correlation_function.hpp"
#include "expectation_val.hpp"
#include "save_to_a_file.hpp"

int main(int argc, char* argv[]){
	//assuming real hoppings within the cluster
    libcommute::expression<double, int , std::string> H;
	/* Hamiltonina parameters*/
	int num_of_sites{0};
	double U{0};
	double t1{-1.},t2{-1.};	
	bool Hubbard{false};
	/* Flags */
	bool single_p{false};
	bool spin_spect{false};
	bool pbc{false};
	bool two_p{false};
	bool spin_spin_corr{false};
	
	read_cmd_line(argc,argv,U,t1,t2,num_of_sites,Hubbard,single_p,spin_spect,pbc,two_p, spin_spin_corr);	
	std::array<double,num_of_sublattices> t{t1,t2};
	std::cout<< "Calculating:\n";
	if(single_p){ std::cout << "-> 1p spectra\n";}
	if(spin_spect){std::cout << "-> spin correlation function\n";}	
	if(Hubbard){ std::cout <<" For Hubbard model ";}
	else{std::cout << " For Hatsugai-Kohmoto ";}
	if(pbc){std::cout<< " Using Periodic Boundary Condition (PBC) ";}
	if(two_p){std::cout << " Calculating 2-point correlation function ";}
	if(spin_spin_corr){std::cout << " Calculating spin-spin correlation function ";}
	std::cout<< " with:\n";
	std::cout<< " N= " << num_of_sites;
	std::cout<< " U= " << U;
	std::cout<< " t1= " << t1;
	std::cout<< " t2= " << t2;
	std::cout << std::endl;
	H=gen_Ham(t, num_of_sites,pbc);
	add_interaction(H,U,num_of_sites,Hubbard,pbc);
	std::cout << "Adding interaction U=" << U << std::endl;
	std::cout << H << std::endl;

	if(single_p){	
		auto tic = std::chrono::high_resolution_clock::now();
		std::vector<std::vector<double>> GF_init{GF_peaks<double>(H,spins_set[0],num_of_sites)};
		auto toc = std::chrono::high_resolution_clock::now();
		auto time=std::chrono::duration_cast<std::chrono::milliseconds>(toc-tic);
		std::cout << "Calculating GF took " << time.count() << " ms\n";

		std::cout << "Sorting GF-peaks: "; 
		sort_GF_peaks(GF_init);
		std::cout<< "DONE\n";

		std::cout << "Binnning peaks:";
		std::vector<std::vector<double>> GF=colllect(GF_init);
		std::cout<< "DONE\n";
		save_to_file(GF,"",U,t,num_of_sites,Hubbard,pbc);
		//save_to_file(GF_init,"_unbinned_",U,t,num_of_sites,Hubbard,pbc);
	}
	if(spin_spect){
		auto tic = std::chrono::high_resolution_clock::now();
		std::vector<std::vector<double>> spin_spect_init{spin_peaks<double>(H,spins_set[0],num_of_sites)};
		auto toc = std::chrono::high_resolution_clock::now();
		auto time=std::chrono::duration_cast<std::chrono::milliseconds>(toc-tic);
		std::cout << "Calculating of spin-spin took " << time.count() << " ms\n";

		std::cout << "Sorting spin-peaks: "; 
		sort_GF_peaks(spin_spect_init);
		std::cout<< "DONE\n";

		std::cout << "Binnning peaks:";
		std::vector<std::vector<double>> spin_spect=colllect(spin_spect_init);
		std::cout<< "DONE\n";
		save_to_file(spin_spect,"spin-spin",U,t,num_of_sites,Hubbard,pbc);
	}
	if(two_p){
		auto tic = std::chrono::high_resolution_clock::now();
		std::vector<std::vector<double>> two_p_corr{two_p_Correlator(H,spins_set[0],num_of_sites)};
		auto toc = std::chrono::high_resolution_clock::now();
		auto time=std::chrono::duration_cast<std::chrono::milliseconds>(toc-tic);
		std::cout << "Calculating of spin-spin took " << time.count() << " ms\n";
		save_to_file(two_p_corr,"_2_p_correlator_",U,t,num_of_sites,Hubbard,pbc);
	}
	if(spin_spin_corr){
		auto tic = std::chrono::high_resolution_clock::now();
		std::vector<std::vector<double>> SS_corr{spin_spin_Correlator(H,num_of_sites)};
		auto toc = std::chrono::high_resolution_clock::now();
		auto time=std::chrono::duration_cast<std::chrono::milliseconds>(toc-tic);
		std::cout << "Calculating of spin-spin took " << time.count() << " ms\n";
		save_to_file(SS_corr,"_s-s_correlator_",U,t,num_of_sites,Hubbard,pbc);
	}		
	return 0;
}
