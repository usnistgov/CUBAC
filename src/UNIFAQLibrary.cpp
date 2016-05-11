#include "UNIFAQLibrary.h"

namespace UNIFAQLibrary{

    rapidjson::Document UNIFAQParameterLibrary::jsonize(std::string &s)
    {
        rapidjson::Document d;
        d.Parse<0>(s.c_str());
        if (d.HasParseError()) {
            throw - 1;
        }
        else {
            return d;
        }
    }
    void UNIFAQParameterLibrary::populate(rapidjson::Value &group_data, rapidjson::Value &interaction_data)
    {
        // Schema should have been used to validate the data already, so by this point we are can safely consume the data without checking ...
        for (rapidjson::Value::ValueIterator itr = group_data.Begin(); itr != group_data.End(); ++itr)
        {
            Group g;
            g.sgi = (*itr)["sgi"].GetInt();
            g.mgi = (*itr)["mgi"].GetInt();
            g.R_k = (*itr)["R_k"].GetDouble();
            g.Q_k = (*itr)["Q_k"].GetDouble();
            groups.push_back(g);
        }
        for (rapidjson::Value::ValueIterator itr = interaction_data.Begin(); itr != interaction_data.End(); ++itr)
        {
            InteractionParameters ip;
            ip.mgi1 = (*itr)["mgi1"].GetInt();
            ip.mgi2 = (*itr)["mgi2"].GetInt();
            ip.a_ij = (*itr)["a_ij"].GetDouble();
            ip.a_ji = (*itr)["a_ji"].GetDouble();
            ip.b_ij = (*itr)["b_ij"].GetDouble();
            ip.b_ji = (*itr)["b_ji"].GetDouble();
            ip.c_ij = (*itr)["c_ij"].GetDouble();
            ip.c_ji = (*itr)["c_ji"].GetDouble();
            interaction_parameters.push_back(ip);
        }
    }
    void UNIFAQParameterLibrary::populate(std::string &group_data, std::string &interaction_data)
    {
        rapidjson::Document group_JSON = jsonize(group_data);
        rapidjson::Document interaction_JSON = jsonize(interaction_data);
        populate(group_JSON, interaction_JSON);
    }
    Group UNIFAQParameterLibrary::get_group(int sgi) {
        for (std::vector<Group>::const_iterator it = groups.cbegin(); it != groups.cend(); ++it) {
            if (it->sgi == sgi) { return *it; }
        }
        throw std::exception("Could not find group");
    }
    InteractionParameters UNIFAQParameterLibrary::get_interaction_parameters(int mgi1, int mgi2) {
        for (std::vector<InteractionParameters>::const_iterator it = interaction_parameters.cbegin(); it != interaction_parameters.cend(); ++it) {
            if (it->mgi1 == mgi1 && it->mgi2 == mgi2) {
                // Correct order, return it
                return *it;
            }
            if (it->mgi2 == mgi1 && it->mgi1 == mgi2) {
                // Backwards, swap the parameters
                InteractionParameters ip = *it;
                ip.swap();
                return ip;
            }
        }
        throw std::exception("Could not find interaction pair");
    }

}; /* namespace UNIFAQLibrary */