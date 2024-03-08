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



#include "Hamiltonian_params.hpp"
#include "read_params.hpp"
#include "gen_Hamiltonian.hpp"
#include "Gen_GF.hpp"
#include "spin_correlation_function.hpp"
#include "expectation_val.hpp"
#include "save_to_a_file.hpp"

int main(int argc, char* argv[]){
	//assuming real hoppings within the cluster
    libcommute::expression<double, int , std::string> H;
	/* Hamiltonina parameters*/
	Hamiltonian_params H_params;
	/* Flags */
	measurments flags;
	
	/* Read parameters from the command line*/
	read_cmd_line(argc,argv,H_params,flags);	
	/* Initiante single particle Hamiltonian*/
	H=gen_Ham<double>(H_params);
	/* Add interactions */
	add_interaction(H,H_params);
	
	if(flags.single_p){	
		
		for(const auto& spin :spins_set){
			if(abs(H_params.mag_field) == 0. && spin == spins_set[1]){continue;} // Without magenetic field do only one spin
			std::cout << "Calculating GF for spin "<< spin << ": ";
			auto tic = std::chrono::high_resolution_clock::now();
			std::vector<std::vector<double>> GF_init{GF_peaks<double>(H,spin,H_params.num_of_sites)};
			auto toc = std::chrono::high_resolution_clock::now();
			auto time=std::chrono::duration_cast<std::chrono::milliseconds>(toc-tic);
			std::cout<<" DONE! Took " << time.count() << " ms\n";
			
			std::cout << "Sorting GF-peaks: "; 
			sort_GF_peaks(GF_init);
			std::cout<< "DONE\n";

			std::cout << "Binnning peaks:";
			std::vector<std::vector<double>> GF=colllect(GF_init);
			std::cout<< "DONE\n";
			std::string extension{"Aw_spin_"+spin};
			save_to_file(GF,extension,H_params);
		}

	}
	if(flags.spin_spect){
		for(const auto& spin :spins_set){
			if(abs(H_params.mag_field) == 0. && spin == spins_set[1]){continue;} // Without magenetic field do only one spin
			std::cout << "Calculating of spin-spin for spin " << spin << " : ";
			auto tic = std::chrono::high_resolution_clock::now();
			std::vector<std::vector<double>> spin_spect_init{spin_peaks<double>(H,spin,H_params.num_of_sites)};
			auto toc = std::chrono::high_resolution_clock::now();
			auto time=std::chrono::duration_cast<std::chrono::milliseconds>(toc-tic);
			std::cout<< "DONE!, took " << time.count() << " ms\n";

			std::cout << "Sorting spin-peaks: "; 
			sort_GF_peaks(spin_spect_init);
			std::cout<< "DONE\n";

			std::cout << "Binnning peaks:";
			std::vector<std::vector<double>> spin_spect=colllect(spin_spect_init);
			std::cout<< "DONE\n";
			std::string extension{ "spin-spin_spin_"+spin};
			save_to_file(spin_spect,extension,H_params);
		}
	}
	if(flags.two_p){
		for(const auto& spin :spins_set){
			if(abs(H_params.mag_field) == 0. && spin == spins_set[1]){continue;} // Without magenetic field do only one spin
			std::cout << "Calculating of 2-p correlator for spin " << spin << " : ";
			auto tic = std::chrono::high_resolution_clock::now();
			std::vector<std::vector<double>> two_p_corr{two_p_Correlator(H,spin,H_params.num_of_sites)};
			auto toc = std::chrono::high_resolution_clock::now();
			auto time=std::chrono::duration_cast<std::chrono::milliseconds>(toc-tic);
			std::cout << "DONE! took " << time.count() << " ms\n";
			std::string extension{"2_p_correlator_spin_"+spin};
			save_to_file(two_p_corr,extension,H_params);
		}
	}
	if(flags.spin_spin_corr){
		std::cout << "Calculating of spin-spin: ";
		auto tic = std::chrono::high_resolution_clock::now();
		std::vector<std::vector<double>> SS_corr{spin_spin_Correlator(H,H_params.num_of_sites)};
		auto toc = std::chrono::high_resolution_clock::now();
		auto time=std::chrono::duration_cast<std::chrono::milliseconds>(toc-tic);
		std::cout<< "DONE! took " << time.count() << " ms\n";
		save_to_file(SS_corr,"s-s_correlator_",H_params);
	}		
	if(flags.electron_density){
		auto tic = std::chrono::high_resolution_clock::now();
		std::vector<std::vector<double>> El_dens{Electron_Density(H,H_params.num_of_sites)};
		auto toc = std::chrono::high_resolution_clock::now();
		auto time=std::chrono::duration_cast<std::chrono::milliseconds>(toc-tic);
		std::cout << "Calculating of electron density took " << time.count() << " ms\n";
		save_to_file(El_dens,"_electron_density_",H_params);
	}
	return 0;
}
