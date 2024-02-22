#include<fstream>
#include<vector>
#include<array>

std::vector<std::vector<double>> colllect(const std::vector<std::vector<double>>& source){
    std::vector<std::vector<double>> res;
    std::vector<double> res_temp(source[0].size(),0.);
    res_temp[0]=source[0][0];
    for(const auto& el: source){
        if(abs(el[0]-res_temp[0]) <1e-8){
            for(auto i=1; i<(int)el.size();i++){
                res_temp[i] += el[i];
            }
        }
        else{
            res.push_back(res_temp);
            res_temp=el;
        }
    }
    res.push_back(res_temp);
    return res;
}


void save_to_file(
    const std::vector<std::vector<double>>& GFs,
    const std::string& prefix,
    const struct Hamiltonian_params& params   
    )
{
    std::ofstream file;
    std::string name{"finite_system_"};
    name += prefix;
    if(params.Hubbard){name += "_Hubbard";}
    if(params.pbc){name += "_PBC_";}
    name += "_N_"+std::to_string(params.num_of_sites);
    int num_of_t{0};
	for(const auto t: params.hopping){
        name += "_t("+std::to_string(num_of_t)+")_"+std::to_string(t);
    }
    name += "_U_"+std::to_string(params.interaction_U);
    name += ".dat";
    file.open(name);
    for (const auto& res: GFs){
        for (const auto& els: res){
            file << els << " ";
        }
        file << "\n";
    }
    file.close();
}
