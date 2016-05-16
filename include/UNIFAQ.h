#ifndef UNIFAQ_H_
#define UNIFAQ_H_

#include <map>

#include "UNIFAQLibrary.h"
namespace UNIFAQ
{
    class UNIFAQMixture
    {
    private:
        /// A const reference to the library of group and interaction parameters
        const UNIFAQLibrary::UNIFAQParameterLibrary &library;

        /// A map from (i, j) indices for subgroup, subgroup indices to the interaction parameters for this pair
        std::map<std::pair<int, int>, UNIFAQLibrary::InteractionParameters> interaction;

        /// A vector of unique groups in this mixture
        std::vector<UNIFAQLibrary::Group> unique_groups;
    
        std::vector<double> mole_fractions;
        std::vector<UNIFAQLibrary::Component> components;
    
    public:
        UNIFAQMixture(const UNIFAQLibrary::UNIFAQParameterLibrary &library) : library(library) {};

        /** 
        * \brief Set all the interaction parameters between groups
        *
        * \param subgroups A vector of the set of the unique Group forming the mixture - these 
        * permutations represent the set of posisble binary interactions
        */
        void set_interaction_parameters();

        /// Set the mole fractions of the components in the mixtures (not the groups)
        void set_mole_fractions(const std::vector<double> &z);

        /// Add a component with the defined groups defined by (count, sgi) pairs
        void add_component(const UNIFAQLibrary::Component &comp);
    
        void set_components(const std::string &identifier_type, std::vector<std::string> identifiers);
    };

} /* namespace UNIFAQ */

#endif