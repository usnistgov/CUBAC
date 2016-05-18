#include "UNIFAQ.h"

void UNIFAQ::UNIFAQMixture::set_interaction_parameters() {
    for (int i = 0; i < unique_groups.size(); ++i) {
        for (int j = i + 1; j < unique_groups.size(); ++j) {
            int mgi1 = unique_groups[i].mgi, mgi2 = unique_groups[j].mgi;
            std::pair< std::pair<int, int>, UNIFAQLibrary::InteractionParameters> m_pair(std::pair<int, int>(mgi1, mgi2), library.get_interaction_parameters(mgi1, mgi2));
            interaction.insert(m_pair);
        }
    }
}

/// Set the mole fractions of the components in the mixtures (not the groups)
void UNIFAQ::UNIFAQMixture::set_mole_fractions(const std::vector<double> &z) {
    this->mole_fractions = z;
    std::size_t N = z.size();
    std::vector<double> r(N), q(N), l(N), phi(N), theta(N), ln_gamma_C(N);
    double summerzr = 0, summerzq = 0, summerzl = 0;
    for (std::size_t i = 0; i < z.size(); ++i) {
        double summerr = 0, summerq = 0;
        const UNIFAQLibrary::Component &c = components[i];
        for (std::size_t j = 0; j < c.groups.size(); ++j) {
            const UNIFAQLibrary::ComponentGroup &cg = c.groups[j];
            summerr += cg.count*cg.group.R_k;
            summerq += cg.count*cg.group.Q_k;
        }
        r[i] = summerr;
        q[i] = summerq;
        summerzr += z[i] * r[i];
        summerzq += z[i] * q[i];
    }
    for (std::size_t i = 0; i < z.size(); ++i) {
        phi[i] = z[i] * r[i] / summerzr;
        theta[i] = z[i] * q[i] / summerzq;
        l[i] = 10.0 / 2.0*(r[i] - q[i]) - (r[i] - 1);
        summerzl += z[i] * l[i];
    }
    for (std::size_t i = 0; i < z.size(); ++i) {
        ln_gamma_C[i] = log(phi[i]/z[i]) + 10.0/2.0*q[i] * log(theta[i]/phi[i]) + l[i] - phi[i]/z[i]*summerzl;
    }
    
    /// Calculate the parameters X and theta for the pure components, which does not depend on temperature
    for (std::size_t i = 0; i < z.size(); ++i){
        int totalgroups = 0;
        const UNIFAQLibrary::Component &c = components[i];
        ComponentData cd;
        double summerxq = 0;
        for (std::size_t j = 0; j < c.groups.size(); ++j) {
            const UNIFAQLibrary::ComponentGroup &cg = c.groups[j];
            double x = static_cast<double>(cg.count);
            double theta = static_cast<double>(cg.count*cg.group.Q_k);
            cd.X.insert( std::pair<int,double>(cg.group.sgi, x) );
            cd.theta.insert(std::pair<int, double>(cg.group.sgi, theta));
            totalgroups += cg.count;
            summerxq += x*cg.group.Q_k;
        }
        /// Now come back through and divide by the total # groups for this fluid
        for (std::map<std::size_t, double>::iterator it = cd.X.begin(); it != cd.X.end(); ++it) {
            it->second /= totalgroups;
            printf("X^(%d)_{%d}: %g\n", static_cast<int>(i + 1), static_cast<int>(it->first), it->second);
        }
        /// Now come back through and divide by the sum(X*Q) for this fluid
        for (std::map<std::size_t,double>::iterator it = cd.theta.begin(); it != cd.theta.end(); ++it){
            it->second /= summerxq;
            printf("theta^(%d)_{%d}: %g\n", static_cast<int>(i+1), static_cast<int>(it->first), it->second);
        }
        pure_data.push_back(cd);
    }
    for (std::size_t i = 0; i < z.size(); ++i) {
        printf("%g %g %g %g %g %g\n", l[i], phi[i], q[i], r[i], theta[i], ln_gamma_C[i]);
    }
}

double UNIFAQ::UNIFAQMixture::Psi(std::size_t sgi1, std::size_t sgi2) {
    std::size_t mgi1 = m_sgi_to_mgi.find(sgi1)->second;
    std::size_t mgi2 = m_sgi_to_mgi.find(sgi2)->second;
    std::map<std::pair<int, int>, UNIFAQLibrary::InteractionParameters>::iterator it = this->interaction.find(std::pair<int,int>(mgi1,mgi2));
    if (it != this->interaction.end()){
        return exp(-(it->second.a_ij + it->second.b_ij*this->m_T + it->second.c_ij*this->m_T*this->m_T)/this->m_T);
    }
    else{
        return 1;
    }
}

double UNIFAQ::UNIFAQMixture::theta_pure(std::size_t i, std::size_t sgi) {
    return pure_data[i].theta.find(sgi)->second;
}

void UNIFAQ::UNIFAQMixture::set_temperature(const double T){
    this->m_T = T;
    for (std::size_t i = 0; i < this->mole_fractions.size(); ++i) {
        const UNIFAQLibrary::Component &c = components[i];
        for (std::size_t k = 0; k < c.groups.size(); ++k) {
            double Q = c.groups[k].group.Q_k;
            int sgik = c.groups[k].group.sgi;
            double sum1 = 0;
            for (std::size_t m = 0; m < c.groups.size(); ++m) {
                int sgim = c.groups[m].group.sgi;
                sum1 += theta_pure(i, sgim)*Psi(sgim, sgik);
            }
            double s = 1 - log(sum1);
            for (std::size_t m = 0; m < c.groups.size(); ++m) {
                int sgim = c.groups[m].group.sgi;
                double sum2 = 0;
                for (std::size_t n = 0; n < c.groups.size(); ++n) {
                    int sgin = c.groups[n].group.sgi;
                    sum2 += theta_pure(i, sgin)*Psi(sgin, sgim);
                }
                s -= theta_pure(i, sgim)*Psi(sgik, sgim)/sum2;
            }
            ComponentData &cd = pure_data[i];
            cd.lnGamma.insert(std::pair<int, double>(sgik, Q*s));
            printf("ln(Gamma)^(%d)_{%d}: %g\n", static_cast<int>(i + 1), sgik, Q*s);
        }
    }
}

/// Add a component with the defined groups defined by (count, sgi) pairs
void UNIFAQ::UNIFAQMixture::add_component(const UNIFAQLibrary::Component &comp) {
    components.push_back(comp);
    // Check if you also need to add group into list of unique groups
    for (std::vector<UNIFAQLibrary::ComponentGroup>::const_iterator it = comp.groups.begin(); it != comp.groups.end(); ++it) {
        bool insert_into_unique = true;
        // if already in unique_groups, don't save it, go to next one
        for (std::vector<UNIFAQLibrary::Group>::const_iterator it2 = unique_groups.cbegin(); it2 != unique_groups.end(); ++it2) {
            if (it2->sgi == it->group.sgi) { insert_into_unique = false; break; }
        }
        if (insert_into_unique) { 
            unique_groups.push_back(it->group); 
            m_sgi_to_mgi.insert(std::pair<std::size_t, std::size_t>(it->group.sgi, it->group.mgi));
        }
    }
}

void UNIFAQ::UNIFAQMixture::set_components(const std::string &identifier_type, std::vector<std::string> identifiers) {
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