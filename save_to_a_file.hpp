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
    const double U,
    const std::array<double,2> t,
    const int sys_size,
    bool Hubbard,
    bool pbc
    )
{
    std::ofstream file;
    std::string name{"finite_system_"};
    name += prefix;
    if(Hubbard){name += "_Hubbard";}
    if(pbc){name += "_PBC_";}
    name += "_N_"+std::to_string(sys_size);
    name += "_t_"+std::to_string(t[0])+"_"+std::to_string(t[1]);
    name += "_U_"+std::to_string(U);
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
