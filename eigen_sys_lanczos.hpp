#include <Eigen/Dense>
#include <lambda_lanczos/lambda_lanczos.hpp>
#include <typeinfo>

/*
inline std::pair<std::vector<double>,std::vector<std::vector<double>>> eigen_sys_lanczos(
    const Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>& Hmat,
    const int num_of_eigstate
    )
{
    // Define matrix multiplication //
    auto mv_mul = [&](const std::vector<double>& in, std::vector<double>& out) {
        auto eigen_in = Eigen::Map<const Eigen::VectorXd>(&in[0], in.size());
        auto eigen_out = Eigen::Map<Eigen::VectorXd>(&out[0], out.size());

        eigen_out = Hmat * eigen_in; // Easy version
        //eigen_out.noalias() += Hmat * eigen_in; // Efficient version
    };

    lambda_lanczos::LambdaLanczos<double> engine(mv_mul, Hmat.cols(), false, num_of_eigstate); // Find 1 minimal eigenvalue
    std::vector<double> res1;
    std::vector<std::vector<double>> res2;
    engine.run(res1, res2);
    std::pair< std::vector<double>, std::vector<std::vector<double>> > res(res1,res2);
    return res;

}
*/

template<typename T>
inline std::pair<std::vector<double>,std::vector<std::vector<T>>> eigen_sys_lanczos(
    const Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>& Hmat,
    const int num_of_eigstate
    )
{
    /* Define matrix multiplication */
    auto mv_mul = [&](const std::vector<T>& in, std::vector<T>& out) {
        auto eigen_in = Eigen::Map<const Eigen::Matrix<T,Eigen::Dynamic,1>>(&in[0], in.size());
        auto eigen_out = Eigen::Map<Eigen::Matrix<T,Eigen::Dynamic,1>>(&out[0], out.size());

        eigen_out = Hmat * eigen_in; // Easy version
        //eigen_out.noalias() += Hmat * eigen_in; // Efficient version
    };

    lambda_lanczos::LambdaLanczos<T> engine(mv_mul, Hmat.cols(), false, num_of_eigstate); // Find 1 minimal eigenvalue
    std::vector<double> res1;
    std::vector<std::vector<T>> res2;
    engine.run(res1, res2);
    std::pair< std::vector<double>, std::vector<std::vector<T>> > res(res1,res2);
    return res;

}