#include<vector>


struct Hamiltonian_params{
    int num_of_sites {0};
    double interaction_U{0.}; // Interaction strenght 
    double mag_field{0.};
    std::vector<double> hopping{0.};
    bool k_dep{false}; // Do we have k-dependence ? 
    double k{0.};
    double el_field{0.};
    bool Hubbard{false}; //Do we want Hubbard interactions instead of H-K?
    bool pbc{false}; // What are the boundary conditions
    bool cmplx{false}; // flag for additional phase factors due to gauge transformation
    std::string model{"Chain"}; //Defult model is the "chain"
};

// Holds infromation on what do we want to calulate 
struct measurments{
    bool single_p{false};
	bool spin_spect{false};
	bool two_p{false};
	bool spin_spin_corr{false};
    bool electron_density{false};
    int retain_states{50};
};


