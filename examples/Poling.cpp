
#include "UNIFAQLibrary.h"

#include <string>
#include <fstream>
#include <sstream>

int main()
{
    std::ifstream t1("../dev/Horstmann_group_data.json"), t2("../dev/Horstmann_interaction_parameters.json");
    std::stringstream b1, b2;
    b1 << t1.rdbuf();
    b2 << t2.rdbuf();

    {
        using namespace UNIFAQLibrary;
        
        UNIFAQParameterLibrary lib;

        lib.populate(b1.str(), b2.str());
        int sgi = 12;
        Group g = lib.get_group(sgi);
        InteractionParameters ip = lib.get_interaction_parameters(1, 2);
        int rr = 0;
    }
}