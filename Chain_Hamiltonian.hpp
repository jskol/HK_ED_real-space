template<typename C>
libcommute::expression<C, int, std::string> gen_Chain_Hamiltonian( 
    const struct Hamiltonian_params& params
)
{
    libcommute::expression<C, int , std::string> H;
    H.clear();

    int num_of_sublattices{int(params.hopping.size())};
    for(auto spin :spins_set){
      	for(int site=0; site< params.num_of_sites; site++){
        	if(site> 0){
                	H +=params.hopping[(site-1)%num_of_sublattices]*c_dag(site-1,spin)*c(site,spin);
            }
           	if(site < params.num_of_sites-1){
                	H += params.hopping[site%num_of_sublattices]*c_dag(site+1,spin)*c(site,spin);
            }
			if(params.k_dep){/*Adding k-dependence of a square lattice */
				H += 2.*params.hopping[0]*cos(params.k)*c_dag(site,spin)*c(site,spin);
			}
		}
		if(params.pbc){
				H += params.hopping[(params.num_of_sites-1)%num_of_sublattices]*( c_dag(0,spin)*c(params.num_of_sites-1,spin)+ c_dag(params.num_of_sites-1,spin)*c(0,spin) );
		} 
    }
    return H;
}