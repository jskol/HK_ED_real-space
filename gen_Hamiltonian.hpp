#include "Kane-Mele_Hamiltonian.hpp"
#include "Chain_Hamiltonian.hpp"

template<typename C>
libcommute::expression<C, int, std::string> gen_Ham( 
    const struct Hamiltonian_params& params
)
{
    
    libcommute::expression<C, int , std::string> H,H_empty;
    H_empty.clear(); //reference empty Hamiltonian for hermicity checks
    if(params.KM){
		H=gen_Kane_Mele_Hamiltonian<double>(params);
	}
	else{
		H=gen_Chain_Hamiltonian<double>(params);
	}

    for(auto spin :spins_set){
      	for(int site=0; site< params.num_of_sites; site++){
 
			if(abs(params.el_field )>0.){
				H += (params.el_field/(params.num_of_sites-1)*site- 0.5*params.el_field)*c_dag(site,spin)*c(site,spin);
			}
			if(abs(params.mag_field)>0.){
				double spin_sign=(spin == spins_set[0] ? 1. :-1.);
				H += 0.5*spin_sign*params.mag_field*c_dag(site,spin)*c(site,spin);
			}

		}
		 
    }

    bool hermicity{ ( (conj(H) - H) == H_empty )? true :false};
    assert(hermicity);
    assert((H != H_empty)); // Make sure H is non-empty
    std::cout << H << std::endl;
    return H;
}

template<typename T>
void add_interaction(
    libcommute::expression<T,int,std::string>& H, 
	const struct Hamiltonian_params& params
    )
{
    if(params.interaction_U==0){
        std::cout << "No interactions\n";
    }
    else{
        std::cout << "Adding terms created by interactions U=" << params.interaction_U << std::endl;
        libcommute::expression<double,int,std::string> H_int;
        if(params.Hubbard){
                std::cout<< "Doing Hubbard model\n";
                for(auto spin_ind=0; spin_ind< (int)spins_set.size(); spin_ind++){
                        for(auto site_up=0; site_up < params.num_of_sites; site_up++){
                                	H_int += 0.5*params.interaction_U
                                        *c_dag(site_up,spins_set[spin_ind])
                                        *c(site_up,spins_set[spin_ind])
                                        *c_dag(site_up,spins_set[(spin_ind+1)%(int)spins_set.size()])
                                        *c(site_up,spins_set[(spin_ind+1)%(int)spins_set.size()]);
                                
                                H_int += -0.5*params.interaction_U*c_dag(site_up,spins_set[spin_ind])*c(site_up,spins_set[spin_ind]);
                        }
                }
        }
        else{
                std::cout<< "Doing real-space HK model\n";
		if(!params.pbc){
			for(auto spin_ind=0; spin_ind< (int)spins_set.size(); spin_ind++){
				for(auto site_up=0; site_up < params.num_of_sites; site_up++){
					for(auto Delta=-site_up; Delta < params.num_of_sites-site_up; Delta++){
						for(auto site_do=0; site_do < params.num_of_sites; site_do++){
							if(site_do-Delta >= 0 && site_do-Delta < params.num_of_sites){
								H_int += (0.5*params.interaction_U/params.num_of_sites)
								    *c_dag(site_up,spins_set[spin_ind])
								    *c((site_up+Delta+params.num_of_sites)%params.num_of_sites,spins_set[spin_ind])
								    *c_dag((site_do+params.num_of_sites)%params.num_of_sites,spins_set[(spin_ind+1)%(int)spins_set.size()])
								    *c((site_do-Delta+params.num_of_sites)%params.num_of_sites,spins_set[(spin_ind+1)%(int)spins_set.size()]);
							}
						}
					}
					H_int += -0.5*params.interaction_U*c_dag(site_up,spins_set[spin_ind])*c(site_up,spins_set[spin_ind]);
				}
		    
			}
		}
		else{
                	std::cout<< "Assuming pbc\n";

			for(auto spin_ind=0; spin_ind< (int)spins_set.size(); spin_ind++){
			    for(auto site_up=0; site_up < params.num_of_sites; site_up++){
				for(auto Delta=-site_up; Delta < params.num_of_sites-site_up; Delta++){
				    for(auto site_do=0; site_do < params.num_of_sites; site_do++){
					H_int += (0.5*params.interaction_U/params.num_of_sites)
					    *c_dag(site_up,spins_set[spin_ind])
					    *c((site_up+Delta+params.num_of_sites)%params.num_of_sites,spins_set[spin_ind])
					    *c_dag((site_do+params.num_of_sites)%params.num_of_sites,spins_set[(spin_ind+1)%(int)spins_set.size()])
					    *c((site_do-Delta+params.num_of_sites)%params.num_of_sites,spins_set[(spin_ind+1)%(int)spins_set.size()]);
				    }
				}
				H_int += -0.5*params.interaction_U*c_dag(site_up,spins_set[spin_ind])*c(site_up,spins_set[spin_ind]);
			    }
			}
		}		
	}
        H += H_int;
    }
}
