#include "Kane-Mele_Hamiltonian.hpp"
#include "Chain_Hamiltonian.hpp"
#include "Haldane_Hamiltonian.hpp"
#include "Graphene_Hamiltonian.hpp"

template<typename C>
libcommute::expression<C, int, std::string> gen_Ham( 
    const struct Hamiltonian_params& params
)
{
    
    libcommute::expression<C, int , std::string> H,H_empty;
    H_empty.clear(); //reference empty Hamiltonian for hermicity checks
    // For graphene-like lattices there is additional
	// Phase added on chosen sites to make the kinetic part
	// of the Hamiltonian real-> this changes the H-K which accumulates
	// phase in some terms, for that reason the Hamiltonian has to be
	// made complex to account for possible phase factors
	// after interactions are included 
	if(params.model=="Chain"){H=gen_Chain_Hamiltonian<C>(params);}
	else if (params.model=="Kane-Mele"){H=gen_Kane_Mele_Hamiltonian<C>(params);}
	else if (params.model=="Haldane"){ H=gen_Haldane_Hamiltonian<C>(params);}
	else if (params.model=="Graphene"){H=gen_Graphene_Hamiltonian<C>(params);}
	else{
		std::cout << "Model "<< params.model << " is not supported!\n";
		std::exit(EXIT_FAILURE);
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


libcommute::static_indices::expr_real<int,std::string> H_interaction_real(
	const struct Hamiltonian_params& params
    )
{
	cx_double one(1.,0.);
	libcommute::static_indices::expr_real<int,std::string> H_int;
    if(params.interaction_U==0){
        std::cout << "No interactions\n";
    }
    else{
        std::cout << "Adding terms created by interactions U=" << params.interaction_U << std::endl;
        if(params.Hubbard){
			std::cout<< "Doing Hubbard model\n";
			for(int spin_ind=0; spin_ind< (int)spins_set.size(); spin_ind++){
					for(int site_up=0; site_up < params.num_of_sites; site_up++){
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
			/* This part needs additional phasing 
				if the unit cell index (ind/2) is even
				the additional phase appearse at the 2nd site
				if the unit cell index is odd 
				the additional phase appearse at the 1st site 
			*/
			if(!params.pbc){
				for(auto spin_ind=0; spin_ind< (int)spins_set.size(); spin_ind++){
					for(int site_up=0; site_up < params.num_of_sites; site_up++){
						for(int Delta=-site_up; Delta < params.num_of_sites-site_up; Delta++){
							for(int site_do=0; site_do < params.num_of_sites; site_do++){
								if(site_do-Delta >= 0 && site_do-Delta < params.num_of_sites){

									int ind2{site_up+Delta};
									int ind4{site_do-Delta};
									H_int += (0.5*params.interaction_U/params.num_of_sites)
									*c_dag(site_up,spins_set[spin_ind])
									*c(ind2,spins_set[spin_ind])
									*c_dag(site_do,spins_set[(spin_ind+1)%(int)spins_set.size()])
									*c(ind4,spins_set[(spin_ind+1)%(int)spins_set.size()]);						
									
									/*	OLD IMPLEMENTATION !!!
										H_int += (0.5*params.interaction_U/params.num_of_sites)
										*c_dag(site_up,spins_set[spin_ind])
										*c((site_up+Delta+params.num_of_sites)%params.num_of_sites,spins_set[spin_ind])
										*c_dag((site_do+params.num_of_sites)%params.num_of_sites,spins_set[(spin_ind+1)%(int)spins_set.size()])
										*c((site_do-Delta+params.num_of_sites)%params.num_of_sites,spins_set[(spin_ind+1)%(int)spins_set.size()]);
									*/
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
						/* Add chemical potential */
						H_int += -0.5*params.interaction_U*c_dag(site_up,spins_set[spin_ind])*c(site_up,spins_set[spin_ind]);
					}
				}
			}		
		}
	}
	return H_int;
}


libcommute::static_indices::expr_complex<int,std::string> H_interaction_cmplx(
	const struct Hamiltonian_params& params
    )
{
	cx_double one(1.,0.);
	libcommute::static_indices::expr_complex<int,std::string> H_int;
    if(params.interaction_U==0){
        std::cout << "No interactions\n";
    }
    else{
        std::cout << "Adding terms created by interactions U=" << params.interaction_U << std::endl;
        if(params.Hubbard){
			std::cout<< "Doing Hubbard model\n";
			for(int spin_ind=0; spin_ind< (int)spins_set.size(); spin_ind++){
					for(int site_up=0; site_up < params.num_of_sites; site_up++){
						H_int = H_int + (0.5*one*params.interaction_U)
							*c_dag(site_up,spins_set[spin_ind])
							*c(site_up,spins_set[spin_ind])
							*c_dag(site_up,spins_set[(spin_ind+1)%(int)spins_set.size()])
							*c(site_up,spins_set[(spin_ind+1)%(int)spins_set.size()]);
						H_int = H_int -0.5*one*params.interaction_U*c_dag(site_up,spins_set[spin_ind])*c(site_up,spins_set[spin_ind]);
					}
			}
        }
        else{
            std::cout<< "Doing real-space HK model\n";
			/* This part needs additional phasing 
				if the unit cell index (ind/2) is even
				the additional phase appearse at the 2nd site
				if the unit cell index is odd 
				the additional phase appearse at the 1st site 
			*/
			if(!params.pbc){
				for(auto spin_ind=0; spin_ind< (int)spins_set.size(); spin_ind++){
					for(int site_up=0; site_up < params.num_of_sites; site_up++){
						for(int Delta=-site_up; Delta < params.num_of_sites-site_up; Delta++){
							for(int site_do=0; site_do < params.num_of_sites; site_do++){
								if(site_do-Delta >= 0 && site_do-Delta < params.num_of_sites){
									/* Calculating the overall phase */								
									cx_double phase(1.,0.);
									if( (((int)site_up / 2)%2 ==0 && (int)site_up%2 == 1) || 
									(((int)site_up / 2)%2 ==1 && (int)site_up%2 == 0) 
									){ phase = phase*exp(-I*0.5*params.k);}

									int ind2{site_up+Delta};
									if( (((int)ind2/2)%2 ==0 && (int)ind2%2 == 1) || 
									(((int)ind2/2)%2 ==1 && (int)ind2%2 == 0) 
									){ phase = phase*exp(I*0.5*params.k);}

									if( (((int)site_do / 2)%2 ==0 && (int)site_do%2 == 1) || 
									(((int)site_do / 2)%2 ==1 && (int)site_do%2 == 0) 
									){ phase = phase*exp(-I*0.5*params.k);}

									int ind4{site_do-Delta};
									if( (((int)ind4/2)%2 ==0 && (int)ind4%2 == 1) || 
									(((int)ind4/2)%2 ==1 && (int)ind4%2 == 0) 
									){ phase = phase*exp(I*0.5*params.k);}

									std::cout << "H_int for : (" << site_up;
									std::cout <<  " , " << ind2;
									std::cout <<  " , " << site_do;
									std::cout <<  " , " << ind4 ;
									std::cout <<  ") the phase factor is " << phase << "\n"; 
									H_int = H_int + (0.5*params.interaction_U/params.num_of_sites)*phase
									    *c_dag(site_up,spins_set[spin_ind])
									    *c((site_up+Delta+params.num_of_sites)%params.num_of_sites,spins_set[spin_ind])
									    *c_dag((site_do+params.num_of_sites)%params.num_of_sites,spins_set[(spin_ind+1)%(int)spins_set.size()])
									    *c((site_do-Delta+params.num_of_sites)%params.num_of_sites,spins_set[(spin_ind+1)%(int)spins_set.size()]);
								}
							}
						}
						H_int = H_int -0.5*one*params.interaction_U*c_dag(site_up,spins_set[spin_ind])*c(site_up,spins_set[spin_ind]);
					}
				}
			}
			else{
				std::cout<< "Assuming pbc\n";
				for(auto spin_ind=0; spin_ind< (int)spins_set.size(); spin_ind++){
					for(auto site_up=0; site_up < params.num_of_sites; site_up++){
						for(auto Delta=-site_up; Delta < params.num_of_sites-site_up; Delta++){
							for(auto site_do=0; site_do < params.num_of_sites; site_do++){
								cx_double phase(1.,0.);
								if( (((int)site_up / 2)%2 ==0 && (int)site_up%2 == 1) || 
								(((int)site_up / 2)%2 ==1 && (int)site_up%2 == 0) 
								){ phase = phase*exp(-I*0.5*params.k);}

								int ind2{site_up+Delta};
								if( (((int)ind2/2)%2 ==0 && (int)ind2%2 == 1) || 
								(((int)ind2/2)%2 ==1 && (int)ind2%2 == 0) 
								){ phase = phase*exp(I*0.5*params.k);}

								if( (((int)site_do / 2)%2 ==0 && (int)site_do%2 == 1) || 
								(((int)site_do / 2)%2 ==1 && (int)site_do%2 == 0) 
								){ phase = phase*exp(-I*0.5*params.k);}

								int ind4{site_do-Delta};
								if( (((int)ind4/2)%2 ==0 && (int)ind4%2 == 1) || 
								(((int)ind4/2)%2 ==1 && (int)ind4%2 == 0) 
								){ phase = phase*exp(I*0.5*params.k);}
								H_int =H_int +(0.5*params.interaction_U/params.num_of_sites)*one
									*c_dag(site_up,spins_set[spin_ind])
									*c(ind2,spins_set[spin_ind])
									*c_dag(site_do,spins_set[(spin_ind+1)%(int)spins_set.size()])
									*c(ind4,spins_set[(spin_ind+1)%(int)spins_set.size()]);

							}
						}
						/* Add chemical potential */
						H_int = H_int-0.5*one*params.interaction_U*c_dag(site_up,spins_set[spin_ind])*c(site_up,spins_set[spin_ind]);
					}
				}
			}		
		}
	}
	return H_int;
}
