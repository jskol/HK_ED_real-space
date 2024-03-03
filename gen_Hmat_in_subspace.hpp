#include <Eigen/Dense>
#include <libcommute/loperator/state_vector_eigen3.hpp>
#include <libcommute/libcommute.hpp>
using libcommute::static_indices::c_dag; // Create an electron
using libcommute::static_indices::c;     // Destroy an electron
using libcommute::static_indices::n;  


template<typename T>
Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> gen_Hmat_in_subspace(
    const libcommute::expression<T, int , std::string>& H,
    const libcommute::space_partition& sp,
    const int subspace_num
){
    /* Set-up hilbert space to get non-zero elements of representation of H */
    /*                      as a linear operator                            */
    auto hs = libcommute::make_hilbert_space(H);
    auto Hop = libcommute::make_loperator(H, hs);
    
    std::vector<libcommute::sv_index_type> basis_states_in_subspace = sp.subspace_basis(subspace_num);
    libcommute::basis_mapper mapper(basis_states_in_subspace);
    auto sp_dim = (int)basis_states_in_subspace.size();
    std::vector<T> psi(sp_dim), phi(sp_dim);
    Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> Hmat(sp_dim,sp_dim);
    for(auto sp_el=0; sp_el < sp_dim; sp_el++){
        std::vector<double> psi(sp_dim), phi(sp_dim);
        psi[sp_el]=1;
        auto psi_view = mapper.make_const_view(psi);
        auto phi_view = mapper.make_view(phi);       
        Hop(psi_view,phi_view);           
        for(auto it=0; it< (int)phi.size(); it++){
            Hmat.col(sp_el)[it]=phi[it];
        }
    }

    return Hmat;

}
