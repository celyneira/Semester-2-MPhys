#include <cmath>
#include <iostream>

#include "Garfield/MediumSiliconCarbide.hh"
#include "Garfield/GarfieldConstants.hh"

namespace Garfield{

MediumSiC::MediumSiC() : Medium(){
    m_className = "MediumSiC";
    m_name = "SiC";

    SetTemperature(300.);
    SetDielectricConstant(9.66);
    Medium::SetAtomicNumber(14+6);
    Medium::SetAtomicWeight(40.11);
    Medium::SetMassDensity(3.211);

    m_driftable = true;
    m_ionisable = true;
    m_microscopic = false;

    m_w = 3.285;
    m_fano = 0.1;

    }

void MediumSiC::GetComponent(const unsigned int i, std::string& label, double & f){
  if (i==0){
    label  = "Si";
    f = 0.5;
  } else if (i == 1){
    label = "C";
     f = 0.5;
  } else {
     std::cerr << m_className << "::GetComponent: Index out of range.\n";
  }
  }


void MediumSiC::SetDoping(const char type, const double c) {
  if (toupper(type) == 'N') {
    m_dopingType = 'n';
    if (c > Small) {
      m_cDop = c;
    } else {
      std::cerr << m_className << "::SetDoping:\n"
                << "    Doping concentration must be greater than zero.\n"
                << "    Using default value for n-type silicon "
                << "(10^12 cm-3) instead.\n";
      m_cDop = 1.e12;
    }
  } else if (toupper(type) == 'P') {
    m_dopingType = 'p';
    if (c > Small) {
      m_cDop = c;
    } else {
      std::cerr << m_className << "::SetDoping:\n"
                << "    Doping concentration must be greater than zero.\n"
                << "    Using default value for p-type silicon "
                << "(10^18 cm-3) instead.\n";
      m_cDop = 1.e18;
    }
  } else if (toupper(type) == 'I') {
    m_dopingType = 'i';
    m_cDop = 0.;
  } else {
    std::cerr << m_className << "::SetDoping:\n"
              << "    Unknown dopant type (" << type << ").\n"
              << "    Available types are n, p and i (intrinsic).\n";
    return;
  }

  m_isChanged = true;
}
void MediumSiC::GetDoping(char& type, double& c) const {
  type = m_dopingType;
  c = m_cDop;
}


bool MediumSiC::ElectronVelocity(const double ex, const double ey,
                                const double ez, const double bx,
                                const double by, const double bz, double& vx,
                                double& vy, double& vz){
    
    vx = vy = vz = 0.;
     if (!Update()) return false;

       if (!m_eVelE.empty()) {
    // Interpolation in user table.
    return Medium::ElectronVelocity(ex, ey, ez, bx, by, bz, vx, vy, vz);
     }

     //Calculate the mobility

     const double e =  sqrt(ex * ex + ey * ey + ez * ez);

     const double mu = -ElectronMobility(e);
     const double mu_c_axis = -ElectronMobilityAnisotropic(e);

     vx = mu * ex;
     vy = mu* ey;
     vz = mu_c_axis * ez;
     
     return true;

                                  }

bool MediumSiC::ElectronTownsend(const double ex, const double ey,
                                     const double ez, const double bx,
                                     const double by, const double bz,
                                        double& alpha) {
    alpha = 0.;
    if (!Update()) return false;

    if (!m_eAlp.empty()) {
        // Interpolation in user table.
        return Medium::ElectronTownsend(ex, ey, ez, bx, by, bz, alpha);
    }

    const double e = sqrt(ex * ex + ey * ey + ez * ez);
    alpha = ElectronAlpha(e);
    return true;
    }

bool MediumSiC::HoleVelocity(const double ex, const double ey,
                                const double ez, const double bx,
                                const double by, const double bz, double& vx,
                                double& vy, double& vz){
    
    vx = vy = vz = 0.;
     if (!Update()) return false;

       if (!m_eVelE.empty()) {
    // Interpolation in user table.
    return Medium::HoleVelocity(ex, ey, ez, bx, by, bz, vx, vy, vz);
     }

     //Calculate the mobility

     const double e =  sqrt(ex * ex + ey * ey + ez * ez);

     const double mu = HoleMobility(e);


     vx = mu * ex;
     vy = mu* ey;
     vz = mu* ez;
     
     return true;

                                  }
bool MediumSiC::HoleTownsend(const double ex, const double ey,
                                 const double ez, const double bx,
                                 const double by, const double bz,
                                 double& alpha) {
  alpha = 0.;
  if (!Update()) return false;

  if (!m_hAlp.empty()) {
    // Interpolation in user table.
    return Medium::HoleTownsend(ex, ey, ez, bx, by, bz, alpha);
  }

  const double e = sqrt(ex * ex + ey * ey + ez * ez);
  alpha = HoleAlpha(e);
  return true;
}

void MediumSiC::SetLowFieldMobility(const double mue, const double muh){
      if (mue <= 0. || muh <= 0.) {
    std::cerr << m_className << "::SetLowFieldMobility:\n"
              << "    Mobility must be greater than zero.\n";
    return;

  }

  m_eMu = mue;
  m_hMu = muh;
  m_hasUserMobility = true;
  m_isChanged = true;

}

void MediumSiC::SetSaturationVelocity(const double vsate,
                                          const double vsath) {
  if (vsate <= 0. || vsath <= 0.) {
    std::cout << m_className << "::SetSaturationVelocity:\n"
              << "    Restoring default values.\n";
    m_hasUserSaturationVelocity = false;
  } else {
    m_eVs = vsate;
    m_hVs = vsath;
    m_hasUserSaturationVelocity = true;
  }

  m_isChanged = true;
}


bool MediumSiC::Initialise() {
  if (!m_isChanged) {
    if (m_debug) {
      std::cerr << m_className << "::Initialise: Nothing changed.\n";
    }
    return true;
  }
  return Update();
}

bool MediumSiC::Update() {
    if (!m_isChanged) return true;
    std::lock_guard<std::mutex> guard(m_mutex);


    if (!m_hasUserMobility){
        UpdateDopingMobility();
    }
    if (!m_hasUserSaturationVelocity) {
    UpdateSaturationVelocity();
  }


  m_isChanged = false;

  return true;
}

void MediumSiC::UpdateDopingMobility() {
      if (m_cDop < 1.e13) {
    m_eMu = m_eMuLat;
    m_hMu = m_hMuLat;
    return;
  }

    double t = m_temperature/300;
    
    double e_mu_max;
    double h_mu_max;

    double e_mu_min;
    double h_mu_min;


    const double e_coefficient_max = e_mu_max_300 * pow(t,e_low_gamma);
    const double h_coefficient_max = h_mu_max_300 * pow(t,h_low_gamma);

    const double e_coefficient_min = e_mu_min_300 * pow(t,e_low_beta);
    const double h_coefficient_min = h_mu_min_300 * pow(t,h_low_beta);
    
    e_mu_max = e_mu_max_300*e_coefficient_max;
    h_mu_max = h_mu_max_300*h_coefficient_max;

    e_mu_min = e_mu_min_300*e_coefficient_min;
    h_mu_min = h_mu_min_300*h_coefficient_min;

    m_eMu = e_mu_min + (e_mu_max - e_mu_min)/(1+pow((m_cDop / e_n_ref),e_low_alpha));
    m_hMu = h_mu_min + (h_mu_max - h_mu_min)/(1+pow((m_cDop / h_n_ref),h_low_alpha));
    
}

void MediumSiC::UpdateSaturationVelocity() {
    m_eVs = 2.2e-2 * pow(m_temperature/300,-0.44);
    m_hVs =  2.2e-2 * pow(m_temperature/300,-0.44);
}

double MediumSiC::ElectronMobility(const double e) const {
    if (e < Small) return 0;

    const double alpha = 1.2;
    
    double alpha_sat = alpha*m_temperature/300;
    double pow_1 = pow(((m_eMuLat*e)/m_eVs),alpha_sat);
    return m_eMu/pow((1+pow_1),(1/alpha_sat));

}

double MediumSiC::ElectronMobilityAnisotropic(const double e) const {
    if (e < Small) return 0;

    const double alpha = 1.2;
    
    double alpha_sat = alpha*m_temperature/300;
    double aniso_mu;
    aniso_mu = m_eMu/e_c_mobility_ratio;
    double pow_1 = pow(((aniso_mu*e)/(m_eVs/e_c_velocity_ratio)),alpha_sat);
    return aniso_mu/pow((1+pow_1),(1/alpha_sat));

}

double MediumSiC::HoleMobility(const double e) const {
    if (e < Small) return 0;

    const double alpha = 1.2;
    
    double alpha_sat = alpha*m_temperature/300;
    double pow_1 = pow(((m_hMuLat*e)/m_hVs),alpha_sat);
    return m_hMu/pow((1+pow_1),(1/alpha_sat));
    
}


double MediumSiC::ElectronAlpha(const double e) const{
    constexpr double hbaromega = 0.19;
    constexpr double T0 = 300;
    constexpr double kb=1.380649e-23;
    const double kb_eV = kb/1.602e-19;


    double e_a = 1.4686e+06;
    double e_b = 1.2075e+07;

    double gamma = tanh(hbaromega / (2 * kb_eV * T0)) / tanh(hbaromega / (2 * kb_eV * T0));
    
    double alpha_e = gamma * e_a * exp(-e_b * gamma / e);
  
    return alpha_e;
}

double MediumSiC::HoleAlpha(const double e) const {
    constexpr double hbaromega = 0.19;
    constexpr double T0 = 300;
    constexpr double kb=1.380649e-23;
    const double kb_eV = kb/1.602e-19;

    double h_a = 5.5222e6;
    double h_b = 1.2724e7;

    double gamma = tanh(hbaromega / (2 * kb_eV * T0)) / tanh(hbaromega / (2 * kb_eV * T0));
   
    double alpha_h = gamma * h_a * exp(-h_b * gamma / e);
    return alpha_h;
}



}