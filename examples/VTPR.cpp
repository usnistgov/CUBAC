#include "UNIFAQLibrary.h"
#include "UNIFAQ.h"

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

const double R = 8.3144598;

class VTPR{
private:
    UNIFAQ::UNIFAQMixture unifaq;
    std::vector<double> Tc, pc, omega, m_ii;
public:
    VTPR(UNIFAQ::UNIFAQMixture &unifaq) : unifaq(unifaq) {set_crits();};
    void set_mole_fractions(const std::vector<double> &z){ unifaq.set_mole_fractions(z); }
    void set_temperature(const double T){ unifaq.set_temperature(T); }
    void set_crits(){
        const std::vector<UNIFAQLibrary::Component> &comps = unifaq.get_components();
        for (std::vector<UNIFAQLibrary::Component>::const_iterator it = comps.begin(); it != comps.end(); ++it){
            Tc.push_back(it->Tc);
            pc.push_back(it->pc);
            omega.push_back(it->acentric);
            double o = omega.back();
            m_ii.push_back(0.37464 + 1.54226*o - 0.26992*o*o);
        }
    };
    /// The attractive part in cubic EOS
    double a_alpha(double T, std::size_t i){
        return 1 + m_ii[i]*(1-sqrt(T/Tc[i]));
    }
    /// Calculate the non-dimensionalized gE/RT term
    double gE_R_RT(){
        const std::vector<double> &z = unifaq.get_mole_fractions();
        double summer = 0;
        for (std::size_t i = 0; i < z.size(); ++i) {
            summer += z[i]*unifaq.ln_gamma_R(i);
        }
        return summer;
    }
    /// The co-volume for the i-th pure component
    double b_ii(std::size_t i){
        return 0.0778*R*Tc[i]/pc[i];
    }
    /// The attractive parameter for the i-th pure component
    double a_ii(std::size_t i){
        return 0.45724*pow(R*Tc[i], 2)/pc[i]*a_alpha(unifaq.get_temperature(), i);
    }
    double bm(){
        double _am,_bm; am_bm(_am, _bm); return _bm;
    }
    double am(){
        double _am,_bm; am_bm(_am, _bm); return _am;
    }
    double cm(){
        return 0;
    }
    /// Calculate both am and bm because am and bm are dependent on each other
    void am_bm(double &am, double &bm){
        const std::vector<double> &z = unifaq.get_mole_fractions();
        double summeram = 0, summerbm = 0;
        for (std::size_t i = 0; i < z.size(); ++i){
            summeram += z[i]*a_ii(i)/b_ii(i);
            for (std::size_t j = 0; j < z.size(); ++j){
                summerbm += z[i]*z[j]*pow((pow(b_ii(i), 0.75) + pow(b_ii(j),0.75))/2.0, 4.0/3.0);
            }
        }
        bm = summerbm;
        am = bm*(summeram + gE_R_RT()/(-0.53087));
    };
    /// Calculate the non-dimensionalized residual Helmholtz energy
    double alphar(double tau, double delta){
        double _am,_bm; am_bm(_am, _bm);
        double _cm = cm();
        double rhor = 1, Tr = 1;
        return -log(1-delta*rhor*(_bm-_cm)) - sqrt(2)*_am*tau/(4*R*Tr*_bm)*log( (1+delta*rhor*(_bm*(1+sqrt(2)+_cm))) / (1+delta*rhor*(_bm*(1-sqrt(2)+_cm))) );
    };
};

int main()
{
    std::ifstream t1("../dev/Horstmann_group_data.json"), 
                  t2("../dev/Horstmann_interaction_parameters.json"), 
                  t3("../dev/group_decompositions.json");
    std::stringstream b1, b2, b3;
    b1 << t1.rdbuf(); b2 << t2.rdbuf(); b3 << t3.rdbuf();
    double rrrr = 0;
    {
        using namespace UNIFAQLibrary;

        UNIFAQParameterLibrary lib;
        lib.populate(b1.str(), b2.str(), b3.str());
        UNIFAQ::UNIFAQMixture mix(lib);
        std::vector<std::string> names(2, "Ethane"); names[1] = "Ethanol";
        mix.set_components("name", names);
        mix.set_interaction_parameters();
        std::vector<double> z(2);
        
        VTPR vtpr(mix);

        FILE* fp = fopen("methane_water.txt","w");
        fprintf(fp, "z_CH4 T(K) gE(J/mol) ln(AC)_1 ln(AC)_2\n");

        for (double T = 300; T < 1000; T += 100) {
            for (double z0 = 0.047; z0 < 1; z0 += 1) {
                z[0] = z0; z[1] = 1 - z[0];
                vtpr.set_mole_fractions(z);
                vtpr.set_temperature(T);
                double ar = vtpr.alphar(1.2, 1.1);
                int rr =0;
            }
        }
        fclose(fp);
    }
}