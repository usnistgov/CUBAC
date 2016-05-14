#include "UNIFAQLibrary.h"
#include "UNIFAQ.h"

#include <string>
#include <fstream>
#include <sstream>

int main()
{
    std::ifstream t1("../dev/Horstmann_group_data.json"), t2("../dev/Horstmann_interaction_parameters.json"), t3("../dev/Kang_decomps.json");
    std::stringstream b1, b2, b3;
    b1 << t1.rdbuf(); b2 << t2.rdbuf(); b3 << t3.rdbuf();
    {
        using namespace UNIFAQLibrary;
        
        UNIFAQParameterLibrary lib;
        lib.populate(b1.str(), b2.str(), b3.str());
        
        UNIFAQMixture mix(lib);
        std::vector<std::string> names(2, "n-Pentane"); names[1] = "Acetone"; 
        mix.set_components("name", names);

        mix.set_interaction_parameters();
    }
}