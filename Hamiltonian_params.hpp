#include<vector>

struct Hamiltonian_params{
    int num_of_sites {0};
    double interaction_U{0.};
    double mag_field{0.};
    std::vector<double> hopping{0.};
    bool k_dep{false};
    double k{0.};
    double el_field{0.};
    bool Hubbard{false};
    bool pbc{false};
    bool KM{false};
    std::string model{"Chain"};
};

struct measurments{
    bool single_p{false};
	bool spin_spect{false};
	bool two_p{false};
	bool spin_spin_corr{false};
    bool electron_density{false};
    int retain_states{50};
};