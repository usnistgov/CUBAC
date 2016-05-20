#include "UNIFAQLibrary.h"
#include "UNIFAQ.h"

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <time.h>

int main()
{
    std::ifstream t1("../dev/Horstmann_group_data.json"), t2("../dev/Horstmann_interaction_parameters.json"), t3("../dev/Kang_decomps.json");
    std::stringstream b1, b2, b3;
    b1 << t1.rdbuf(); b2 << t2.rdbuf(); b3 << t3.rdbuf();
    double rrrr = 0;
    {
        using namespace UNIFAQLibrary;

        UNIFAQParameterLibrary lib;
        lib.populate(b1.str(), b2.str(), b3.str());

        double t1 = clock();
        long N = 10000;
        UNIFAQ::UNIFAQMixture mix(lib);
        std::vector<std::string> names(2, "Acetone"); names[1] = "n-Pentane";
        mix.set_components("name", names);
        mix.set_interaction_parameters();

        for (int ii = 0; ii < N; ++ii) {
            std::vector<double> z(2, 0.047); z[1] = 1 - z[0];
            mix.set_mole_fractions(z);
            mix.set_temperature(307);
            rrrr += mix.activity_coefficient(0);
        }
        double t2 = clock();
        std::cout << rrrr/N << " " << (t2 - t1) / ((double)CLOCKS_PER_SEC) / ((double)N)*1e6 << std::endl;
    }
}