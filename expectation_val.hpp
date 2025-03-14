/*
Routines for calculating the expectation values
(static quantities)
*/

#include<vector>

#include <libcommute/loperator/state_vector_eigen3.hpp>
#include <libcommute/libcommute.hpp>
using libcommute::static_indices::c_dag; // Create an electron
using libcommute::static_indices::c;     // Destroy an electron
using libcommute::static_indices::n;  

#include<chrono>
#include<omp.h>
#include<algorithm>

template<typename T>
std::vector<std::pair<double, std::vector<double>>> two_p_Correlator(
    const libcommute::expression<T,int, std::string>& H,
    std::string spin,
    const int sys_size
){
    std::vector<std::vector<double>> result(sys_size ,std::vector<double>(sys_size)); //inititate with a fixed size
    std::cout<< "Looking for the groundstate:";
    auto tic = std::chrono::high_resolution_clock::now(); //Timing
	std::vector<std::pair<libcommute::sv_index_type,int>> gs_sub{get_GS_subspace(H)};
    auto toc=std::chrono::high_resolution_clock::now();
    auto time= std::chrono::duration_cast<std::chrono::milliseconds>(toc-tic);
    std::cout << "Finding the groundstate took " << time.count() << " ms\n";

    auto hs = libcommute::make_hilbert_space(H);
    auto Hop = libcommute::make_loperator(H, hs);
    auto sp =libcommute::space_partition(Hop, hs);
    
    //Calculate Z factor-> Ground state degeneracy if T is (almost) zero
    double Z{0};
    for(const auto& gs_sub_ind: gs_sub){
	    Z += gs_sub_ind.second; //Add all the ground-state degeneracies
    }
     	   
    std::cout << "Running over the GS's\n";
    for(const auto& gs_sub_ind: gs_sub){
        std::vector<libcommute::sv_index_type> basis_states_in_GS_subspace = sp.subspace_basis(gs_sub_ind.first);
        Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> GS_Hmat= gen_Hmat_in_subspace(H,sp,gs_sub_ind.first);
        std::pair<std::vector<double> ,std::vector< std::vector<T>>> GS_eigen_sys{eigen_sys_lanczos(GS_Hmat,gs_sub_ind.second)};
        std::cout << "Looping over lowest energy states within a subspace\n";
        for(int mult=0; mult< gs_sub_ind.second; mult++){
            std::vector<T> GS_in_H(hs.dim()),res_vec_target(hs.dim());
            for(auto i=0; i< (int)basis_states_in_GS_subspace.size(); i++)
            {
                GS_in_H[basis_states_in_GS_subspace[i]] = GS_eigen_sys.second[mult][i];
            }
            for(int site_L=0; site_L< sys_size; site_L++){
                for(int site_R=0; site_R< sys_size; site_R++){
                    std::vector<T> res_vec(hs.dim());
                    auto lop=libcommute::make_loperator(c(site_L,spin)*c_dag(site_R,spin),hs);
                    lop(GS_in_H, res_vec);
                    double res_temp=0;
                    cx_double one(1.,0.);
                    for(const auto& non_zero_el: basis_states_in_GS_subspace){
                        res_temp += std::real(std::conj(one*GS_in_H[non_zero_el])*res_vec[non_zero_el]);
                    }
                    result[site_L][site_R] += res_temp/Z;
                }
            }

        }
    }
    /* Recast into format preffered by print_to_file*/
    std::vector< std::pair< double,std::vector<double> > > res_print;
    std::pair< double,std::vector<double> > temp;
    for( int row=0; row< sys_size; row++){
        temp.first=Z;
        temp.second.clear();
        for(int col=0; col < sys_size; col++){
            temp.second.push_back(result[row][col]);
        }
        res_print.push_back(temp);    
    }
    return res_print;
}


template<typename T>
std::vector< std::pair<double, std::vector<double> >> spin_spin_Correlator(
    const libcommute::expression<T,int, std::string>& H,
    const int sys_size
){
    std::vector<std::vector<double>> result(sys_size ,std::vector<double>(sys_size)); //inititate with a fixed size
    std::cout<< "Looking for the groundstate:";
    auto tic = std::chrono::high_resolution_clock::now();
	std::vector<std::pair<libcommute::sv_index_type,int>> gs_sub{get_GS_subspace(H)};
    auto toc=std::chrono::high_resolution_clock::now();
    auto time= std::chrono::duration_cast<std::chrono::milliseconds>(toc-tic);
    std::cout << "Finding the groundstate took " << time.count() << " ms\n";

    auto hs = libcommute::make_hilbert_space(H);
    auto Hop = libcommute::make_loperator(H, hs);
    auto sp =libcommute::space_partition(Hop, hs);
    
   /*Calculate Z factor*/
    double Z{0};
    for(const auto& gs_sub_ind: gs_sub){
	    Z += gs_sub_ind.second; //Add all the ground-state degeneracies
    }
    std::cout << "Running over the GS's\n";
    for(const auto& gs_sub_ind: gs_sub){
        std::vector<libcommute::sv_index_type> basis_states_in_GS_subspace = sp.subspace_basis(gs_sub_ind.first);
        Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> GS_Hmat= gen_Hmat_in_subspace(H,sp,gs_sub_ind.first);
        std::pair<std::vector<double> ,std::vector< std::vector<T>>> GS_eigen_sys{eigen_sys_lanczos(GS_Hmat,gs_sub_ind.second)};

        std::cout << "Looping over lowest energy states within a subspace\n";
        for(int mult=0; mult< gs_sub_ind.second; mult++){
            std::vector<T> GS_in_H(hs.dim()),res_vec_target(hs.dim());
            for(auto i=0; i< (int)basis_states_in_GS_subspace.size(); i++)
            {
                GS_in_H[basis_states_in_GS_subspace[i]] = GS_eigen_sys.second[mult][i];
            }
            for(int site_L=0; site_L< sys_size; site_L++){
                for(int site_R=0; site_R< sys_size; site_R++){
                    std::vector<T> res_vec(hs.dim(),0.);
                    libcommute::expression<double,int,std::string> spin_op;
                    spin_op.clear();
                    spin_op += 0.25*(c_dag(site_L,spins_set[0])*c(site_L,spins_set[0])-c_dag(site_L,spins_set[1])*c(site_L,spins_set[1]))*(c_dag(site_R,spins_set[0])*c(site_R,spins_set[0])-c_dag(site_R,spins_set[1])*c(site_R,spins_set[1])); //S_z(j)S_z(i)
                    spin_op += 0.5*(c_dag(site_L,spins_set[0])*c(site_L,spins_set[1])*c_dag(site_R,spins_set[1])*c(site_R,spins_set[0])); //S_+(j)S_-(i)
                    spin_op += 0.5*(c_dag(site_L,spins_set[1])*c(site_L,spins_set[0])*c_dag(site_R,spins_set[0])*c(site_R,spins_set[1])); // S_-(j)S_+(i)
                    auto lop=libcommute::make_loperator(spin_op,hs);
                    lop(GS_in_H, res_vec);
                    double res_temp=0;
                    cx_double one(1.,0.);
                    for(const auto& non_zero_el: basis_states_in_GS_subspace){
                        res_temp += std::real(std::conj(one*GS_in_H[non_zero_el])*res_vec[non_zero_el] );
                    }
                    result[site_L][site_R] += res_temp/Z;
                }
            }

        }
    }
    /* Recast into format preffered by print_to_file*/
    std::vector< std::pair< double,std::vector<double> > > res_print;
    std::pair< double,std::vector<double> > temp;
    for( int row=0; row< sys_size; row++){
        temp.first=Z;
        temp.second.clear();
        for(int col=0; col < sys_size; col++){
            temp.second.push_back(result[row][col]);
        }
        res_print.push_back(temp);    
    }
    return res_print;
}


template<typename T>
std::pair<double, std::vector<double>> Electron_Density(
    const libcommute::expression<T,int, std::string>& H,
    const int sys_size
){

    std::vector<std::vector<double>> result(sys_size ,std::vector<double>(sys_size)); //inititate with a fixed size
    std::cout<< "Looking for the groundstate:";
    auto tic = std::chrono::high_resolution_clock::now();
	std::vector<std::pair<libcommute::sv_index_type,int>> gs_sub{get_GS_subspace(H)};
    auto toc=std::chrono::high_resolution_clock::now();
    auto time= std::chrono::duration_cast<std::chrono::milliseconds>(toc-tic);
    std::cout << "Finding the groundstate took " << time.count() << " ms\n";

    auto hs = libcommute::make_hilbert_space(H);
    auto Hop = libcommute::make_loperator(H, hs);
    auto sp =libcommute::space_partition(Hop, hs);

   /*Calculate Z factor*/
    double Z{0};
    for(const auto& gs_sub_ind: gs_sub){
	     Z+= gs_sub_ind.second; //Add all the ground-state degeneracies
    }

    std::vector<double> dens_vec(sys_size);
    std::cout << "Running over the GS's\n";
    for(const auto& gs_sub_ind: gs_sub){
        std::vector<libcommute::sv_index_type> basis_states_in_GS_subspace = sp.subspace_basis(gs_sub_ind.first);
        Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> GS_Hmat= gen_Hmat_in_subspace(H,sp,gs_sub_ind.first);
        std::pair<std::vector<double> ,std::vector< std::vector<T>>> GS_eigen_sys{eigen_sys_lanczos(GS_Hmat,gs_sub_ind.second)};

        std::cout << "Looping over lowest energy states within a subspace\n";
        for(int mult=0; mult< gs_sub_ind.second; mult++){
            std::vector<T> GS_in_H(hs.dim()),res_vec_target(hs.dim());
            for(auto i=0; i< (int)basis_states_in_GS_subspace.size(); i++)
            {
                GS_in_H[basis_states_in_GS_subspace[i]] = GS_eigen_sys.second[mult][i];
            }
            for(auto site=0; site< sys_size; site++){
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
    }
    return std::make_pair(Z, dens_vec);
}