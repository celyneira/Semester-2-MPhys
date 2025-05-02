#ifndef G_MEDIUM_SIC_H
#define G_MEDIUM_SIC_H

#include <array>
#include <mutex>
#include <string>
#include <vector>

#include "Medium.hh"

namespace Garfield{

class MediumSiC : public Medium {
public:
    //Constructor
    MediumSiC();

    //destructor
    virtual ~MediumSiC() {}

    bool IsSemiconductor() const override { return true; }

    /// Set doping concentration [cm-3] and type ('i', 'n', 'p').
    void SetDoping(const char type, const double c);
    /// Retrieve doping concentration. 
    void GetDoping(char& type, double& c) const;
    void GetComponent(const unsigned int i, std::string& label, 
                    double& f) override;
        // Electron transport parameters
    bool ElectronVelocity(const double ex, const double ey, const double ez,
                            const double bx, const double by, const double bz,
                            double& vx, double& vy, double& vz) override;
    bool ElectronTownsend(const double ex, const double ey, const double ez,
                            const double bx, const double by, const double bz,
                            double& alpha) override;
    double ElectronMobility() override { return m_eMu; }
    double ElectronMobilityAnisotropic() { return m_eMu/0.83; }
    // Hole transport parameters
    bool HoleVelocity(const double ex, const double ey, const double ez,
                        const double bx, const double by, const double bz,
                        double& vx, double& vy, double& vz) override;
    bool HoleTownsend(const double ex, const double ey, const double ez,
                        const double bx, const double by, const double bz,
                        double& alpha) override;

    double HoleMobility() override { return m_hMu; }
    //low field mobility
    void SetLowFieldMobility(const double mue, const double muh);

    void SetHighFieldMobility();

    void SetSaturationVelocity(const double vsate, const double vsath);


    bool Initialise();



private:
    double m_bandGap = 3.26;

    // Doping
    char m_dopingType = 'i';
    // Doping concentration
    double m_cDop = 0.;

    //Lattice mobility
    double m_eMuLat = 0.945e-6;
    double m_hMuLat = 0.125e-6;

    //Low field mobility
    double m_eMu = 0.945e-6;
    double m_hMu = 0.125e-6;

    //High field mobility model parameters
    double alpha_sat = 1.2;
    double delta_sat = -0.44;
    double beta_sat = 1;

    //Saturation velocity
    double m_eVs = 2.2e-2;
    double m_hVs = 2.2e-2;

    //Low field parameters
    double e_mu_max_300 = 0.95e-6;
    double h_mu_max_300 = 0.125e-6;

    double e_mu_min_300 = 0.04e-6;
    double h_mu_min_300 = 0.0159e-6;

    double e_n_ref = 1.94e17;
    double h_n_ref = 1.76e19;

    double e_low_alpha = 0.61;
    double h_low_alpha = 0.34;

    double e_low_beta = -0.5;
    double h_low_beta = -0.5;

    double e_low_gamma = -2.40;
    double h_low_gamma = -2.15;

    //impact ionisation parameters
    double e_alpha_impact = 1.4686e6;
    double h_alpha_impact = 5.5222e6;
    double e_beta_impact = 1.2075e+07;
    double h_beta_impact = 1.2724e+07;

    //anisotropy parameters
    double e_c_mobility_ratio = 0.83; //mupar=muperp/0.83
    double h_c_mobility_ratio = 1;

    double e_c_velocity_ratio = 1.16;
    
    void UpdateSaturationVelocity();

    double ElectronMobility(const double e) const;
    double ElectronMobilityAnisotropic(const double e) const;
    double ElectronAlpha(const double e) const;

    double HoleMobility(const double e) const;
    double HoleAlpha(const double e) const;
    
    void UpdateDopingMobility();
    void UpdateHighFieldMobility();
    bool Update();

    bool m_hasUserMobility = false;
    bool m_hasUserSaturationVelocity = false;

    std::mutex m_mutex;












};
}
#endif
