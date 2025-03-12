class Hamiltonian{
    public:
        libcommute::expression<double,int,std::string> real;
        libcommute::expression<cx_double,int,std::string> cmplx;
        /* Contructor with interaction term*/
        Hamiltonian(Hamiltonian_params params){
            if(params.cmplx){
                std::cout << "Creating complex Hamiltonian\n";
                cmplx=gen_Ham<cx_double>(params); //initiate Hamiltonian with complex coeff's
                libcommute::static_indices::expr_complex<int,std::string> H_int{H_interaction_cmplx(params)};
                cmplx=H_int+cmplx;
            }
            else{
                std::cout << "Creating real Hamiltonian\n";
                real=gen_Ham<double>(params); //initiate Hamiltonian with only real coeff
                libcommute::static_indices::expr_real<int,std::string> H_int{H_interaction_real(params)};
                real = real +H_int;
            }
        }
};  