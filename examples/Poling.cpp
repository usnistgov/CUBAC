#include "UNIFAQLibrary.h"
#include "UNIFAQ.h"

#include <string>
#include <fstream>
#include <sstream>

int main()
{
    std::ifstream t1("../dev/Horstmann_group_data.json"), 
                  t2("../dev/Horstmann_interaction_parameters.json"), 
                  t3("../dev/group_decompositions.json");
    std::stringstream b1, b2, b3;
    b1 << t1.rdbuf(); b2 << t2.rdbuf(); b3 << t3.rdbuf();
    {
        using namespace UNIFAQLibrary;
        
        UNIFAQParameterLibrary lib;
        lib.populate(b1.str(), b2.str(), b3.str());
        
        UNIFAQ::UNIFAQMixture mix(lib);
        std::vector<std::string> names(2, "Acetone"); names[1] = "n-Pentane"; 
        mix.set_components("name", names);

        mix.set_interaction_parameters();
        std::vector<double> z(2,0.047); z[1] = 1-z[0];
        mix.set_mole_fractions(z);
        mix.set_temperature(307);
        for (std::size_t i = 0; i < 2; ++i){
            printf("gamma_{%d}: %g\n", i, mix.activity_coefficient(i));
        }
    }
}