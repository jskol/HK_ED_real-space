/*Model specific*/
const int num_of_sublattices{2};
/**/

template<typename C>
libcommute::expression<C, int, std::string> gen_Ham( 
    const std::array<C,num_of_sublattices>& hoppings, 
    const int num_of_sites,
    bool pbc
)
{
    
    libcommute::expression<C, int , std::string> H,H_empty;
    H_empty.clear(); //reference empty Hamiltonian for hermicity checks
    H.clear();

    for(auto spin :spins_set){
      	for(int site=0; site< num_of_sites; site++){
        	if(site> 0){
                	H += hoppings[(site-1)%num_of_sublattices]*c_dag(site-1,spin)*c(site,spin);
            	}
           	if(site < num_of_sites-1){
                	H += hoppings[site%num_of_sublattices]*c_dag(site+1,spin)*c(site,spin);
            	}
        }
	if(pbc){H +=hoppings[(num_of_sites-1)%num_of_sublattices]*( c_dag(0,spin)*c(num_of_sites-1,spin)+ c_dag(num_of_sites-1,spin)*c(0,spin) ); } 
	
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
    const double U,
    const int number_of_sites,
    bool Hubbard,
    bool pbc
    )
{
    if(U==0){
        std::cout << "No interactions\n";
    }
    else{
        std::cout << "Adding terms created by interactions U=" << U << std::endl;
        libcommute::expression<double,int,std::string> H_int;
        if(Hubbard){
                std::cout<< "Doing Hubbard model\n";
                for(auto spin_ind=0; spin_ind< (int)spins_set.size(); spin_ind++){
                        for(auto site_up=0; site_up < number_of_sites; site_up++){
                                	H_int += 0.5*U
                                        *c_dag(site_up,spins_set[spin_ind])
                                        *c(site_up,spins_set[spin_ind])
                                        *c_dag(site_up,spins_set[(spin_ind+1)%(int)spins_set.size()])
                                        *c(site_up,spins_set[(spin_ind+1)%(int)spins_set.size()]);
                                
                                H_int += -0.5*U*c_dag(site_up,spins_set[spin_ind])*c(site_up,spins_set[spin_ind]);
                        }
                }
        }
        else{
                std::cout<< "Doing real-space HK model\n";
		if(!pbc){
			for(auto spin_ind=0; spin_ind< (int)spins_set.size(); spin_ind++){
				for(auto site_up=0; site_up < number_of_sites; site_up++){
					for(auto Delta=-site_up; Delta < number_of_sites-site_up; Delta++){
						for(auto site_do=0; site_do < number_of_sites; site_do++){
							if(site_do-Delta >= 0 && site_do-Delta < number_of_sites){
								H_int += (0.5*U/number_of_sites)
								    *c_dag(site_up,spins_set[spin_ind])
								    *c((site_up+Delta+number_of_sites)%number_of_sites,spins_set[spin_ind])
								    *c_dag((site_do+number_of_sites)%number_of_sites,spins_set[(spin_ind+1)%(int)spins_set.size()])
								    *c((site_do-Delta+number_of_sites)%number_of_sites,spins_set[(spin_ind+1)%(int)spins_set.size()]);
							}
						}
					}
					H_int += -0.5*U*c_dag(site_up,spins_set[spin_ind])*c(site_up,spins_set[spin_ind]);
				}
		    
			}
		}
		else{
                	std::cout<< "Assuming pbc\n";

			for(auto spin_ind=0; spin_ind< (int)spins_set.size(); spin_ind++){
			    for(auto site_up=0; site_up < number_of_sites; site_up++){
				for(auto Delta=-site_up; Delta < number_of_sites-site_up; Delta++){
				    for(auto site_do=0; site_do < number_of_sites; site_do++){
					H_int += (0.5*U/number_of_sites)
					    *c_dag(site_up,spins_set[spin_ind])
					    *c((site_up+Delta+number_of_sites)%number_of_sites,spins_set[spin_ind])
					    *c_dag((site_do+number_of_sites)%number_of_sites,spins_set[(spin_ind+1)%(int)spins_set.size()])
					    *c((site_do-Delta+number_of_sites)%number_of_sites,spins_set[(spin_ind+1)%(int)spins_set.size()]);
				    }
				}
				H_int += -0.5*U*c_dag(site_up,spins_set[spin_ind])*c(site_up,spins_set[spin_ind]);
			    }
			}
		}		
	}
        H += H_int;
    }
}
