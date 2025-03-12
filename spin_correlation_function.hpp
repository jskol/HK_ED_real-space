#include<vector>
#include <Eigen/Dense>
#include <libcommute/loperator/state_vector_eigen3.hpp>
#include <libcommute/libcommute.hpp>
using libcommute::static_indices::c_dag; // Create an electron
using libcommute::static_indices::c;     // Destroy an electron
using libcommute::static_indices::n;  

//#include "eigen_sys_lanczos.hpp"
//#include "gen_Hmat_in_subspace.hpp"

#include<chrono>
#include<omp.h>
#include<algorithm>



// reatin states defined in Gen_Gf
//int retain_states{40};

template<typename T>
std::vector<std::pair<double,std::vector<double>>> spin_peaks(
    const libcommute::expression<T,int, std::string>& H,
    const Hamiltonian_params& params,
    const std::string spin,
    const int retain_states
    )
{
    /* Get the location of the Ground state(s)*/  
    std::cout<< "\nLooking for the groundstate:";
    auto tic = std::chrono::high_resolution_clock::now();
	std::vector<std::pair<libcommute::sv_index_type,int>> gs_sub{get_GS_subspace(H)};
    auto toc=std::chrono::high_resolution_clock::now();
    auto time= std::chrono::duration_cast<std::chrono::milliseconds>(toc-tic);
    std::cout << "Finding the groundstate took " << time.count() << " ms\n";
    std::string opp_spin=(spin==spins_set[0] ? spins_set[1] : spins_set[0]);
    std::vector< std::pair< double,std::vector<double> > > res; // vector storing pairs of En-GS excitation energy and its probability
    std::vector<double> dens_vec(params.num_of_sites); // Vector storing the local electron density

    auto hs = libcommute::make_hilbert_space(H);
    auto Hop = libcommute::make_loperator(H, hs);
    
    auto sp =libcommute::space_partition(Hop, hs);
    std::vector<std::pair<libcommute::matrix_elements_map<double>,libcommute::matrix_elements_map<double>>> transitions;
    for(int site=0; site<params.num_of_sites;site++){
         auto c_dag_lop=libcommute::make_loperator(1.*c_dag(site,opp_spin)*c(site,spin),hs);
         auto c_lop=libcommute::make_loperator(1.*c_dag(site,spin)*c(site,opp_spin),hs);
         auto temp_map=sp.merge_subspaces(c_dag_lop,c_lop,hs,true);
         transitions.push_back(temp_map);
    }

    
    /*Calculate Z factor - T=0 limit assumed*/
    double Z{0};
    for(const auto& gs_sub_ind: gs_sub){ 
        Z += gs_sub_ind.second; //Add all the ground-state degeneracies
    }
    
    /* Get each Ground state*/
    std::cout << "Running over the GS's\n";
    std::pair<std::vector<double> ,std::vector< std::vector<T>>> GS_eigen_sys; //Will hold Lanczos results 
    for(const auto& gs_sub_ind: gs_sub){
    
        std::vector<libcommute::sv_index_type> basis_states_in_GS_subspace = sp.subspace_basis(gs_sub_ind.first);
        GS_eigen_sys.first.clear();
        GS_eigen_sys.second.clear();
        Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> GS_Hmat= gen_Hmat_in_subspace(H,sp,gs_sub_ind.first);
        std::pair<std::vector<double> ,std::vector< std::vector<T>>> GS_eigen_sys{eigen_sys_lanczos(GS_Hmat,gs_sub_ind.second)}; // Eigen energies and generic (type T) eigenvectors
        GS_Hmat.resize(0,0);
        double GS_en{GS_eigen_sys.first[0]};
    
        // Loop over the degenerate GS's within a certain subspace
        std::cout << "Looping over the lowest energy states within a subspace\n";
        for(int mult=0; mult< gs_sub_ind.second; mult++){
            //for(int mult=0; mult<1; mult++){
            std::vector<T> GS_in_H(hs.dim());
            // Create a GS vector
	        for(auto i=0; i< (int)basis_states_in_GS_subspace.size(); i++)
            {
                GS_in_H[basis_states_in_GS_subspace[i]] = GS_eigen_sys.second[mult][i];
	        }

            std::pair<std::vector<double> ,std::vector< std::vector<T>>> diag_n; /* Here will be the eigensystem in n-th subspace kept*/
    
            /* Calculate  |<n | c^+ |GS>|^2 only the electron-like excitations(E>E_GS) */

            /* Find subspaces conncted to the certain subspace with GS*/
            std::vector<libcommute::sv_index_type> connections;
            for (auto site=0; site<params.num_of_sites; site++){
                auto c_dag_lop=libcommute::make_loperator(1.*c_dag(site,opp_spin)*c(site,spin),hs);
                auto conn_c_dag=sp.find_connections(c_dag_lop,hs);
                /* Discard the elements not connecting to the groundstate */
                for(const auto& conn: conn_c_dag){
                    if(conn.first == gs_sub_ind.first){
                        connections.push_back(conn.second);
                    }
                }
            }
            /* Remove duplicates from the list */
            std::sort(connections.begin(),connections.end());
            connections.erase(std::unique(connections.begin(),connections.end()),connections.end());

            std::cout << "Groudstate subspace #"<<gs_sub_ind.first  << " is connected to subspaces:";
            for(const auto& conn : connections){std::cout << " #" << conn;}
            std::cout<< " through operator c_dag"<< std::endl;
            
            /* Start Calculations of the excitations */
            std::cout<< "Calculating weights for c_dag transitions from GS:";
            for(const auto& sp_it: connections){ /* Loop over target subspace of the excitation  C_dag*/
                
                /* Find eigen-states in a given subspace */
                std::vector<libcommute::sv_index_type> basis_states_in_subspace_n = sp.subspace_basis(sp_it);
                auto sp_dim = (int)basis_states_in_subspace_n.size();
                
                
                Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> Hmat=gen_Hmat_in_subspace(H,sp,sp_it);
                int n_states=(sp_dim < retain_states ? sp_dim : retain_states); /* Save only $retain_states lowest energy states */
                
                /* clear the container of Lanczos results*/
                diag_n.first.clear();
                diag_n.second.clear();
                
                /* Get eigenvectors of the excited subspace*/
                diag_n=eigen_sys_lanczos(Hmat,n_states);

                /*Empty the matrix Hmat*/
                Hmat.resize(0,0);

                /* Vector storing pairs of exciation energy and its probablity from GS to each excited state */
                std::vector<std::pair<double, std::vector<double>>> res_temp;
                res_temp.clear();
               
                /* Loop over all eigenstates of the excited state subspace (N+1) */
                for(auto eigen_states=0; eigen_states < (int)diag_n.first.size(); eigen_states++){
                
                    /*represent this eigenstate in the full H-space*/
                    /*Here the other eigenstate vector will be kept*/
                    std::vector<T> n_vec(hs.dim());
                    std::vector<T> res_vec(hs.dim());
                    for(int el=0; el< sp_dim; el++){
                        n_vec[basis_states_in_subspace_n[el]]=diag_n.second[eigen_states][el];
                    }

                    /*Clear the containers, duplicate res_arr for parallelization */
                    std::pair<double,std::vector<double>> res_arr; // stores energy with respect to GS and exitation element
                    res_arr.first=diag_n.first[eigen_states]-GS_en;
                    res_arr.second.clear();
                    for (int site=0; site<params.num_of_sites;site++){
                        /*Calculate the weigth of c_dag*/
                        auto c_dag_lop=libcommute::make_loperator(1.*c_dag(site,opp_spin)*c(site,spin),hs);
                        cx_double weight(0.,0.);
                        cx_double one(1.,0.);
                        c_dag_lop(GS_in_H, res_vec); /* c^+ acting in the full H-space*/
                        for(const auto& non_zero_el: basis_states_in_subspace_n){
                            weight = weight + std::conj(one*n_vec[non_zero_el])*res_vec[non_zero_el];
                        }
                        
                        res_arr.second.push_back(std::norm(weight)/Z);
                    }
                    /* Testing */
                    //std::cout << "\nweights=" << res_arr.first  << " ";
                    //for( const auto& els: res_arr.second){std::cout << els << " ";}
                    /* */
                    res_temp.push_back(res_arr);
                }

                for(const auto& flush: res_temp){res.push_back(flush);}
            }  
        
            std::cout << "\nDONE\n";

            for(auto site=0; site< params.num_of_sites; site++){
                auto n_op=libcommute::make_loperator(c_dag(site,spins_set[0])*c(site, spins_set[0])+c_dag(site,spins_set[1])*c(site, spins_set[1]),hs );
                std::vector<T> res_vec(hs.dim());
                n_op(GS_in_H,res_vec);
                double dens{0.};
                cx_double one(1.,0.);
                for(const auto& non_zero_el: basis_states_in_GS_subspace){
                    dens += std::real(std::conj(one*GS_in_H[non_zero_el])*res_vec[non_zero_el]);
                }
                dens_vec[site] += dens/Z;
            }
        }    

    /* Onto the next subspace with the groundstate*/
    }


    std::cout<< "Density distribution:\n";
    for(auto site=0; site< params.num_of_sites; site++){
        std::cout << "n_" << site << " = " << dens_vec[site] << "\n";
    }
    return res;
}
