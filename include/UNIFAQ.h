#ifndef UNIFAQ_H_
#define UNIFAQ_H_

#include <map>

#include "UNIFAQLibrary.h"

/// Structure containing data for the pure fluid in the mixture
struct ComponentData {
    std::map<std::size_t, double> X, theta, lnGamma;
};

namespace UNIFAQ
{
    class UNIFAQMixture
    {
    private:
        double m_T; ///< The temperature in K

        /// A const reference to the library of group and interaction parameters
        const UNIFAQLibrary::UNIFAQParameterLibrary &library;

        /// A map from (i, j) indices for subgroup, subgroup indices to the interaction parameters for this pair
        std::map<std::pair<int, int>, UNIFAQLibrary::InteractionParameters> interaction;

        /// A map from SGI to MGI
        std::map<std::size_t, std::size_t> m_sgi_to_mgi;

        /// A vector of unique groups in this mixture
        std::vector<UNIFAQLibrary::Group> unique_groups;
    
        std::vector<double> mole_fractions;

        std::vector<UNIFAQLibrary::Component> components;

        std::vector<ComponentData> pure_data;
    
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

        /// Set the mole fractions of the components in the mixtures (not the groups)
        void set_temperature(const double T);

        double Psi(std::size_t sgi1, std::size_t sgi2);

        double theta_pure(std::size_t i, std::size_t sgi);

        /// Add a component with the defined groups defined by (count, sgi) pairs
        void add_component(const UNIFAQLibrary::Component &comp);
    
        void set_components(const std::string &identifier_type, std::vector<std::string> identifiers);
    };

} /* namespace UNIFAQ */

#endif