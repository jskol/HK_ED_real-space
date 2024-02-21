#include<vector>

struct Hamiltonian_params{
    int N{0};
    double U{0.};
    double mag_field{0.};
    std::vector<double> hopping{0.};
    double k{0.};
    double el_field{0.};
    bool Hubbard{false};
    bool pbc{false};
    bool k_dep{false};
};

struct measurments{
    bool single_p{false};
	bool spin_spect{false};
	bool two_p{false};
	bool spin_spin_corr{false};
};