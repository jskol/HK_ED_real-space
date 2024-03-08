template<typename C>
libcommute::expression<C, int, std::string> gen_Graphene_Hamiltonian( 
    const struct Hamiltonian_params& params
)
{
    libcommute::expression<C, int , std::string> H;
    H.clear();
    /* Sanity check */
    assert(params.k_dep);
    assert(params.hopping.size()==2);
    assert(params.num_of_sites%2 ==0 ); /* Has to have even number of sites since its a two-sublattice model*/
    
    int num_of_sublattices{int(params.hopping.size())};
    for(auto spin :spins_set){
        for(int sub=0; sub< params.num_of_sites; sub += num_of_sublattices) {
            H += 2*params.hopping[0]*cos(params.k)*c_dag(sub,spin)*c(sub+1,spin);
            H += 2*params.hopping[0]*cos(params.k)*c_dag(sub+1,spin)*c(sub,spin);
            if(sub < params.num_of_sites-num_of_sublattices){
                H +=   params.hopping[0]*(c_dag(sub+1,spin)*c(sub+num_of_sublattices,spin) +c_dag(sub+num_of_sublattices,spin)*c(sub+1,spin));
            }
        }
        if(params.pbc){/*TODO*/;}
    }

    return H;
}