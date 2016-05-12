#ifndef UNIFAQ_H_
#define UNIFAQ_H_

#include <map>

#include "UNIFAQLibrary.h"

class UNIFAQMixture
{
private:
    /// A const reference to the library of group and interaction parameters
    const UNIFAQLibrary::UNIFAQParameterLibrary &library;

    /// A map from (i, j) indices for subgroup, subgroup indices to the interaction parameters for this pair
    std::map<std::pair<int, int>, UNIFAQLibrary::InteractionParameters> interaction;
    std::vector<double> mole_fractions;
    
public:
    UNIFAQMixture(const UNIFAQLibrary::UNIFAQParameterLibrary &library) : library(library) {};

    void set_interaction_parameters() {
        // TODO: generalize this
        std::pair< std::pair<int, int>, UNIFAQLibrary::InteractionParameters> m_pair(std::pair<int, int>(1, 2), library.get_interaction_parameters(1, 2));
        interaction.insert(m_pair);
    }

    /// Set the mole fractions of the components in the mixtures (not the groups)
    void set_mole_fractions(const std::vector<double> &z) {
        this->mole_fractions = z;
    }

};

#endif