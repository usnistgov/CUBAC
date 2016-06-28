#include "UNIFAQLibrary.h"
#include "UNIFAQ.h"

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

const double R = 8.3144598;

int main()
{
    std::ifstream t1("../dev/Horstmann_group_data.json"), t2("../dev/Horstmann_interaction_parameters.json"), t3("../dev/group_decompositions.json");
    std::stringstream b1, b2, b3;
    b1 << t1.rdbuf(); b2 << t2.rdbuf(); b3 << t3.rdbuf();
    double rrrr = 0;
    {
        using namespace UNIFAQLibrary;

        UNIFAQParameterLibrary lib;
        lib.populate(b1.str(), b2.str(), b3.str());

        UNIFAQ::UNIFAQMixture mix(lib);
        std::vector<std::string> names(2, "Acetone"); names[1] = "n-Pentane";
        mix.set_components("name", names);
        mix.set_interaction_parameters();
        std::vector<double> z(2);

        FILE* fp = fopen("methane_water.txt","w");
        fprintf(fp, "z_CH4 T(K) gE(J/mol) ln(AC)_1 ln(AC)_2\n");

        for (double T = 300; T < 1000; T += 100) {
            for (double z0 = 0.047; z0 < 1; z0 += 1) {
                z[0] = z0; z[1] = 1 - z[0];
                mix.set_mole_fractions(z);
                mix.set_temperature(T);
                double gE_over_RT = z[0]*log(mix.activity_coefficient(0)) + z[1]*log(mix.activity_coefficient(1));
                double gE = gE_over_RT*R*T;
                printf("%g %g %g %g %g\n", z0, T, gE, mix.activity_coefficient(0), mix.activity_coefficient(1));
                fprintf(fp, "%g %g %g %g %g\n", z0, T, gE, mix.activity_coefficient(0), mix.activity_coefficient(1));
            }
        }
        fclose(fp);
    }
}