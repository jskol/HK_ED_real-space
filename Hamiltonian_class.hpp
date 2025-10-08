template<typename T>
class Hamiltonian{
    public:
        libcommute::expression<T,int,std::string> H;

        /* Contructor with interaction term*/
        Hamiltonian(Hamiltonian_params params){
                H=gen_Ham<T>(params); //initiate Hamiltonian
                libcommute::expression<T,int,std::string> H_int{H_interaction<T>(params)};
                H=H+H_int;
        }

};  