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
    void set_interaction_parameters() {
        for (int i = 0; i < unique_groups.size(); ++i){
            for (int j = i+1; j < unique_groups.size(); ++j){
                int mgi1 = unique_groups[i].mgi, mgi2 = unique_groups[j].mgi;
                std::pair< std::pair<int, int>, UNIFAQLibrary::InteractionParameters> m_pair(std::pair<int, int>(mgi1, mgi2), library.get_interaction_parameters(mgi1, mgi2));
                interaction.insert(m_pair);
            }
        }
    }

    /// Set the mole fractions of the components in the mixtures (not the groups)
    void set_mole_fractions(const std::vector<double> &z) {
        this->mole_fractions = z;
    }

    /// Add a component with the defined groups defined by (count, sgi) pairs
    void add_component(const UNIFAQLibrary::Component &comp){
        components.push_back(comp);
        // Check if you also need to add group into list of unique groups
        for(std::vector<UNIFAQLibrary::ComponentGroup>::const_iterator it = comp.groups.begin(); it != comp.groups.end(); ++it){
            bool insert_into_unique = true;
            // if already in unique_groups, don't save it, go to next one
            for (std::vector<UNIFAQLibrary::Group>::const_iterator it2 = unique_groups.cbegin(); it2 !=unique_groups.end(); ++it2 ){
                if (it2->sgi == it->group.sgi){ insert_into_unique = false; break;}
            }
            if (insert_into_unique){ unique_groups.push_back(it->group); }
        }
    }
    
    void set_components(const std::string &identifier_type, std::vector<std::string> identifiers) {
        if (identifier_type == "name") {
            // Iterate over the provided names
            for (std::vector<std::string>::const_iterator it = identifiers.cbegin(); it != identifiers.cend(); ++it) {
                // Get and add the component
                UNIFAQLibrary::Component c = library.get_component("name", *it);
                add_component(c);
            }
        }
        else {
            throw std::exception("Cannot understand identifier_type");
        }
    }

};

#endif