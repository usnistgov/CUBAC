#ifndef UNIFAQ_H_
#define UNIFAQ_H_

#include <map>

#include "UNIFAQLibrary.h"

/// A structure containing a group (its count, index, etc.) for a subgroup forming a part of a component
struct ComponentGroup{
    int count;
    UNIFAQLibrary::Group group;
    ComponentGroup(const int count, const UNIFAQLibrary::Group group) : count(count), group(group) {};
};

/// A structure containing the groups and additional information for a component
struct Component
{
    std::vector<ComponentGroup> groups;
};

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
    std::vector<Component> components;
    
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
    void add_component(const std::vector<std::pair<int, UNIFAQLibrary::Group> > &groups){
        Component comp;
        for(std::vector<std::pair<int, UNIFAQLibrary::Group> >::const_iterator it = groups.begin(); it != groups.end(); ++it){
            UNIFAQLibrary::Group group = library.get_group(it->second.sgi);
            ComponentGroup cg(it->first, it->second);
            comp.groups.push_back(cg);
            bool insert_into_unique = true;
            // if already in unique_groups, don't save it, got to next one
            for (std::vector<UNIFAQLibrary::Group>::const_iterator it2 = unique_groups.cbegin(); it2 !=unique_groups.end(); ++it2 ){
                if (it2->sgi == it->second.sgi){ insert_into_unique = false; break;}
            }
            if (insert_into_unique){ unique_groups.push_back(it->second); }
        }
        components.push_back(comp);
    }

    /// TODO: use group decomposition list to populate 
    void set_components(const std::string &identifier_type, std::vector<int> identifiers){
        if (identifier_type == "PubChem"){
        }
        else if (identifier_type == "CAS"){
        }
        else{
            throw std::exception("Cannot understand identifier_type");
        }
    }
    /// TODO: use group decomposition list to populate 
    void set_components(const std::string &identifier_type, std::vector<std::string> identifiers) {
        if (identifier_type == "name") {
            for (std::vector<std::string>::const_iterator it = identifiers.cbegin(); it != identifiers.cend(); ++it){
                if (*it == "pentane"){
                    std::vector<std::pair<int, UNIFAQLibrary::Group> > pentane;
                    pentane.push_back(std::pair<int, UNIFAQLibrary::Group>(2, library.get_group(1)));
                    pentane.push_back(std::pair<int, UNIFAQLibrary::Group>(3, library.get_group(2)));
                    add_component(pentane);
                }
                else if (*it == "acetone"){
                    std::vector<std::pair<int, UNIFAQLibrary::Group> > acetone;
                    acetone.push_back(std::pair<int, UNIFAQLibrary::Group>(1, library.get_group(1)));
                    acetone.push_back(std::pair<int, UNIFAQLibrary::Group>(1, library.get_group(18)));
                    add_component(acetone);
                }
                else{
                    throw -1;
                }
            }
        }
        else {
            throw std::exception("Cannot understand identifier_type");
        }
    }

};

#endif