#include <vector>

#include <Eigen/Dense>
#include <libcommute/loperator/state_vector_eigen3.hpp>
#include <libcommute/libcommute.hpp>
using libcommute::static_indices::c_dag; // Create an electron
using libcommute::static_indices::c;     // Destroy an electron
using libcommute::static_indices::n;  

#include "eigen_sys_lanczos.hpp"
#include "gen_Hmat_in_subspace.hpp"

#include<chrono>
#include<omp.h>
#include<algorithm>

int retain_states{50};

/*This function return vector of pairs (index of the invariant subspace, ground state degeneracy) */

template<typename T>
std::vector<std::pair<libcommute::sv_index_type,int>> get_GS_subspace(const libcommute::expression<T, int , std::string>& H)
{
    std::vector<std::pair<libcommute::sv_index_type,int>> res; /* Here will the final outcom go*/
    
    /* Set-up hilbert space to get non-zero elements of representation of H */
    /*                      as a linear operator                            */
    auto hs = libcommute::make_hilbert_space(H);
    auto Hop = libcommute::make_loperator(H, hs);
    auto sp = libcommute::space_partition(Hop, hs);
    double GS_en{100.};
    double deg_crit{1e-4};
    std::vector<std::tuple<libcommute::sv_index_type,double,int>> lowest_en(sp.n_subspaces()); /*Temporary storage of pairs serial number of a subspace and its lowest energy and the degeneracy,  */
    std::cout << "The calcuatlion will run on "  << omp_get_max_threads() << " threads\n";
    /* Start looping over the invariant subspaces*/
    #pragma omp parallel shared(lowest_en)
    {
        #pragma omp for
        for(libcommute::sv_index_type subspace=0; subspace< sp.n_subspaces(); subspace++){       
            
            Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> Hmat=gen_Hmat_in_subspace(H,sp,subspace);
            /*Get the lowest energy state and its energy --> Using Lanczos*/
            int GS_deg{1};
            double E_spacing{0};
            while(E_spacing < deg_crit){
                std::pair<std::vector<double> ,std::vector< std::vector<double>>> eigen_sys{eigen_sys_lanczos(Hmat,GS_deg+1)};
                if(eigen_sys.first[0]>GS_en+deg_crit){ break;} // skiping further  analysis of this subspace, its irrelevant
                
		        E_spacing=abs(eigen_sys.first[GS_deg]- eigen_sys.first[0]);
                if(E_spacing >deg_crit){
                    if(eigen_sys.first[0]-GS_en < deg_crit ){
                        GS_en=eigen_sys.first[0];
                    }
                    std::cout << "There is " << GS_deg << " states with the same energy\n";
                    std::tuple<libcommute::sv_index_type, double,int> res{subspace,eigen_sys.first[0],GS_deg};
                    lowest_en[subspace]=res;
                }
                else{
                    GS_deg++;
                }
            }
        }
    }
    for(const auto& els: lowest_en){
        if( abs(std::get<1>(els)-GS_en) <deg_crit){
            std::pair<libcommute::sv_index_type,int> temp{std::get<0>(els),std::get<2>(els)};
            res.push_back(temp);
        }
    }
    std::cout<< "Groundstate with energy E=" << GS_en << "is in subspace:\n"; 
    for(const auto& pos: res){ 
        std::cout << "#" << pos.first << " and it's " << pos.second-1 << " times degenerate,  ";
    }
    std::cout << "\n";
    return res;
}


/* Find degeneracy of the ground state !!*/


template<typename T>
std::vector<std::vector<double>> GF_peaks(
    const libcommute::expression<T,int, std::string>& H,
    const std::string spin,
    //const int site,
    const int sys_size)
{
    /* Get the location of the Ground state(s)*/  
    std::cout<< "Looking for the groundstate:";
    auto tic = std::chrono::high_resolution_clock::now();
	std::vector<std::pair<libcommute::sv_index_type,int>> gs_sub{get_GS_subspace(H)};
    auto toc=std::chrono::high_resolution_clock::now();
    auto time= std::chrono::duration_cast<std::chrono::milliseconds>(toc-tic);
    std::cout << "Finding the groundstate took " << time.count() << " ms\n";

    std::vector<std::vector<double>> res;
    std::vector<double> dens_vec(sys_size);

    auto hs = libcommute::make_hilbert_space(H);
    auto Hop = libcommute::make_loperator(H, hs);
    
    auto sp =libcommute::space_partition(Hop, hs);
    
    /*Calculate Z factor*/
    double Z{0};
    for(const auto& gs_sub_ind: gs_sub){
	Z += gs_sub_ind.second; //Add all the ground-state degeneracies
    }
    /* Get each Ground state*/
    std::cout << "Running over the GS's\n";
    for(const auto& gs_sub_ind: gs_sub){
        std::vector<libcommute::sv_index_type> basis_states_in_GS_subspace = sp.subspace_basis(gs_sub_ind.first);
        Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> GS_Hmat= gen_Hmat_in_subspace(H,sp,gs_sub_ind.first);
        std::pair<std::vector<double> ,std::vector< std::vector<double>>> GS_eigen_sys{eigen_sys_lanczos(GS_Hmat,gs_sub_ind.second)};
        double GS_en{GS_eigen_sys.first[0]};
    
    
        /*This will store temporary resulting vectors*/
        //std::vector<double> res_vec(hs.dim());
        // Loop over degenerate GS's within a certain subspace
        std::cout << "Looping over lowest energy states within a subspace\n";
        for(int mult=0; mult< gs_sub_ind.second; mult++){
        //for(int mult=0; mult<1; mult++){
            std::vector<double> GS_in_H(hs.dim(),0.),res_vec_target(hs.dim(),0.);
            // Create a GS vector
	    for(auto i=0; i< (int)basis_states_in_GS_subspace.size(); i++)
            {
                GS_in_H[basis_states_in_GS_subspace[i]] = GS_eigen_sys.second[mult][i];
	    }
            /* Loop over invariant subspaces and get elements of the eigen-system*/
            std::pair<std::vector<double> ,std::vector< std::vector<double>>> diag_n; /* Here will be the eigensystem in n-th subspace kept*/
    
            /* Calculate <GS| c |n><n | c^+ |GS> */
            std::vector<libcommute::sv_index_type> connections;
            for (auto site=0; site<sys_size;site++){
                auto c_dag_lop=libcommute::make_loperator(1.*c_dag(site,spin),hs);
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

            std::cout<< "Calculating weights for c_dag transitions from GS:";

        

            for(const auto& sp_it: connections){
                
                /* Find eigen-states in a given subspace */
                std::vector<libcommute::sv_index_type> basis_states_in_subspace_n = sp.subspace_basis(sp_it);
                auto sp_dim = (int)basis_states_in_subspace_n.size();
                Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> Hmat=gen_Hmat_in_subspace(H,sp,sp_it);
                int n_states=(sp_dim < retain_states ? sp_dim : retain_states); /* Save only $retain_states lowest energy states */
                
                /* clear the container */
                diag_n.first.clear();
                diag_n.second.clear();
                
                /* Get n_states eigenvectors */
                diag_n=eigen_sys_lanczos(Hmat,n_states);
                /*Store the results for parallelization*/
            
                
                std::vector<std::vector<double>> res_temp(diag_n.first.size());
                #pragma omp parallel shared(diag_n,res_temp)
                {
                    #pragma omp for
                    for(auto eigen_states=0; eigen_states < (int)diag_n.first.size(); eigen_states++){
                    
                        /*represent this eigenstate in the full H-space*/
                        /*Here the other eigen state vector will be kept*/
                        std::vector<T> n_vec(hs.dim(),0.);
                        std::vector<T> res_vec(hs.dim(),0.);
                        for(int el=0; el< sp_dim; el++){
                            n_vec[basis_states_in_subspace_n[el]]=diag_n.second[eigen_states][el];
                        }

                        /*Clear the containers, duplicate res_arr for parallelization */
                        std::vector<double> res_arr(sys_size+1);
                        res_arr[0]=diag_n.first[eigen_states]-GS_en;
                        for (auto site=0; site<sys_size;site++){
                            /*Calculate the weigth of c_dag*/
                            auto c_dag_lop=libcommute::make_loperator(1.*c_dag(site,spin),hs);
                            double weight{0.};
                            c_dag_lop(GS_in_H, res_vec); /* c^+ acting in the full H-space*/
                            for(const auto& non_zero_el: basis_states_in_subspace_n){
                                weight += n_vec[non_zero_el]*res_vec[non_zero_el];
                            }
                            res_arr[site+1] =(weight*weight)/Z;
                        }
                        res_temp[eigen_states]=res_arr;
                    }
                }

                for(const auto& flush: res_temp){res.push_back(flush);}
            }  
        
            std::cout << "DONE\n";

            /* Repeat for the swapped operator order --> second term in the commutator */
            /* <GS| c^+ |n><n| c |GS> */        
            /* Discard the elements not connecting to the groundstate */
            connections.clear();
            for (auto site=0; site<sys_size;site++){
                auto c_lop=libcommute::make_loperator(1.*c(site,spin),hs);
                auto conn_c=sp.find_connections(c_lop,hs);
                for(const auto& conn: conn_c){
                    if(conn.first == gs_sub_ind.first){
                        //std::cout <<conn.first << " " << conn.second << "\n";
                        connections.push_back(conn.second);
                    }
                }
            }

            /* Remove duplicates from the list */
            std::sort(connections.begin(),connections.end());
            connections.erase(std::unique(connections.begin(),connections.end()),connections.end());

            std::cout << "Groudstate subspace #"<<gs_sub_ind.first  << " is connected to subspaces:";
            for(const auto& conn : connections){std::cout << " #" << conn;}
            std::cout<< " through operator c"<< std::endl;

            std::cout<< "Calculating weights for c transitions from GS:";
            for(const auto& sp_it: connections){ 
                /* Find eigen-states in a given subspace */
                std::vector<libcommute::sv_index_type> basis_states_in_subspace_n = sp.subspace_basis(sp_it);
                auto sp_dim = (int)basis_states_in_subspace_n.size();
                Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> Hmat=gen_Hmat_in_subspace(H,sp,sp_it);
                int n_states=(sp_dim < retain_states ? sp_dim : retain_states); /* Save only $retain_states lowest energy states */

                /* clear the container */
                diag_n.first.clear();
                diag_n.second.clear();

                /*Get the n_states eigenvectors */
                diag_n=eigen_sys_lanczos(Hmat,n_states);
                /*Store the results for parallelization*/
                std::vector<std::vector<double>> res_temp(diag_n.first.size());
                #pragma omp parallel shared(diag_n,res_temp)
                {
                    #pragma omp for
                    for(auto eigen_states=0; eigen_states < (int)diag_n.first.size(); eigen_states++){
                        /*represent this eigenstate in the full H-space*/
                        /*Here the other eigen state vector will be kept*/
                        std::vector<double> res_vec(hs.dim());
                        std::vector<T> n_vec(hs.dim(),0.);
                        for(int el=0; el< sp_dim; el++){
                            n_vec[basis_states_in_subspace_n[el]]=diag_n.second[eigen_states][el];
                        }
                        /* Clear the containter*/
                        std::vector<double> res_arr(sys_size+1);
                        res_arr[0]=(GS_en-diag_n.first[eigen_states]);
                        for (auto site=0; site<sys_size;site++){
                            /* Calculate weights */
                            auto c_lop=libcommute::make_loperator(1.*c(site,spin),hs);
                            double weight{0.};
                            c_lop(GS_in_H,res_vec);
                            for(const auto& non_zero_el:basis_states_in_subspace_n){
                                weight += n_vec[non_zero_el]*res_vec[non_zero_el];
                            }
                            res_arr[site+1] =(weight*weight)/Z;
                        }
                        res_temp[eigen_states]=res_arr;
                    }
                }
                for(const auto& flush: res_temp){res.push_back(flush);}
            }
                    
            std::cout << "DONE\n";

            for(auto site=0; site< sys_size; site++){
                auto n_op=libcommute::make_loperator(c_dag(site,spins_set[0])*c(site, spins_set[0])+c_dag(site,spins_set[1])*c(site, spins_set[1]),hs );
                std::vector<double> res_vec(hs.dim());
                n_op(GS_in_H,res_vec);
                double dens{0.};
                for(const auto& non_zero_el: basis_states_in_GS_subspace){
                    dens += GS_in_H[non_zero_el]*res_vec[non_zero_el];
                }
                dens_vec[site] += dens/Z;
            }
        }    
    }
    std::cout<< "Density distribution:\n";
    for(auto site=0; site< sys_size; site++){
        std::cout << "n_" << site << " = " << dens_vec[site] << "\n";
    }
    return res;
}


void sort_GF_peaks(std::vector<std::vector<double>>& temp){
    std::sort(temp.begin(),temp.end(),
        [](const std::vector<double>& a, const std::vector<double>&b){return a[0]< b[0];} // Lambda to sort by the first el
        );
}
