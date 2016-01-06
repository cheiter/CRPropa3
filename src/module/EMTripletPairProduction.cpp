#include "crpropa/module/EMTripletPairProduction.h"
#include "crpropa/Units.h"
#include "crpropa/ParticleID.h"
#include "crpropa/ParticleMass.h"
#include "crpropa/Random.h"

#include <fstream>
#include <limits>
#include <stdexcept>

namespace crpropa {

EMTripletPairProduction::EMTripletPairProduction(PhotonField photonField,
		bool haveElectrons, double limit) {
	setPhotonField(photonField);
	this->haveElectrons = haveElectrons;
	this->limit = limit;
}

void EMTripletPairProduction::setPhotonField(PhotonField photonField) {
	this->photonField = photonField;
	switch (photonField) {
	case CMB:
		setDescription("EMTripletPairProduction: CMB");
		initRate(getDataPath("EMTripletPairProduction_CMB.txt"));
		initCumulativeRate(getDataPath("EMTripletPairProduction_CDF_CMB.txt"));
		break;
	case IRB:  // default: Kneiske '04 IRB model
	case IRB_Kneiske04:
		setDescription("EMTripletPairProduction: IRB (Kneiske 2004)");
		initRate(getDataPath("EMTripletPairProduction_IRB_Kneiske04.txt"));
		initCumulativeRate(getDataPath("EMTripletPairProduction_CDF_IRB_Kneiske04.txt"));
		break;
	case IRB_Stecker05:
		setDescription("EMTripletPairProduction: IRB (Stecker 2005)");
		initRate(getDataPath("EMTripletPairProduction_IRB_Stecker05.txt"));
		initCumulativeRate(getDataPath("EMTripletPairProduction_CDF_IRB_Stecker05.txt"));
		break;
	case IRB_Franceschini08:
		setDescription("EMTripletPairProduction: IRB (Franceschini 2008)");
		initRate(getDataPath("EMTripletPairProduction_IRB_Franceschini08.txt"));
		initCumulativeRate(getDataPath("EMTripletPairProduction_CDF_IRB_Franceschini08.txt"));
		break;
	case IRB_Finke10:
		setDescription("EMTripletPairProduction: IRB (Finke 2010)");
		initRate(getDataPath("EMTripletPairProduction_IRB_Finke10.txt"));
		initCumulativeRate(getDataPath("EMTripletPairProduction_CDF_IRB_Finke10.txt"));
		break;
	case IRB_Dominguez11:
		setDescription("EMTripletPairProduction: IRB (Dominguez 2011)");
		initRate(getDataPath("EMTripletPairProduction_IRB_Dominguez11.txt"));
		initCumulativeRate(getDataPath("EMTripletPairProduction_CDF_IRB_Dominguez11.txt"));
		break;
	case IRB_Gilmore12:
		setDescription("EMTripletPairProduction: IRB (Gilmore 2012)");
		initRate(getDataPath("EMTripletPairProduction_IRB_Gilmore12.txt"));
		initCumulativeRate(getDataPath("EMTripletPairProduction_CDF_IRB_Gilmore12.txt"));
		break;
	case URB_Protheroe96:
		setDescription("EMTripletPairProduction: URB (Protheroe 1996)");
		initRate(getDataPath("EMTripletPairProduction_URB_Protheroe96.txt"));
		initCumulativeRate(getDataPath("EMTripletPairProduction_CDF_URB_Protheroe96.txt"));
		break;
	default:
		throw std::runtime_error(
				"EMTripletPairProduction: unknown photon background");
	}
}

void EMTripletPairProduction::setHaveElectrons(bool haveElectrons) {
	this->haveElectrons = haveElectrons;
}

void EMTripletPairProduction::setLimit(double limit) {
	this->limit = limit;
}

void EMTripletPairProduction::initRate(std::string filename) {
	std::ifstream infile(filename.c_str());

	if (!infile.good())
		throw std::runtime_error(
				"EMTripletPairProduction: could not open file " + filename);

	// clear previously loaded interaction rates
	tabElectronEnergy.clear();
	tabInteractionRate.clear();

	while (infile.good()) {
		if (infile.peek() != '#') {
			double a, b;
			infile >> a >> b;
			if (infile) {
				tabElectronEnergy.push_back(pow(10, a) * eV);
				tabInteractionRate.push_back(b / Mpc);
			}
		}
		infile.ignore(std::numeric_limits < std::streamsize > ::max(), '\n');
	}
	infile.close();
}

void EMTripletPairProduction::initCumulativeRate(std::string filename) {
	std::ifstream infile(filename.c_str());

	if (!infile.good())
		throw std::runtime_error(
				"EMTripletPairProduction: could not open file " + filename);

	// clear previously loaded interaction rates
	tabE.clear();
  tabs.clear();
	tabCumulativeRate.clear();

	while (infile.good()) {
		if (infile.peek() != '#') {
			double a, b, c;
			infile >> a >> b >> c;
			if (infile) {
				tabE.push_back(pow(10, a) * eV);
        tabs.push_back(pow(10,b) * eV * eV);
				tabCumulativeRate.push_back(c / Mpc);
			}
		}
		infile.ignore(std::numeric_limits < std::streamsize > ::max(), '\n');
	}
	infile.close();
}

void EMTripletPairProduction::performInteraction(Candidate *candidate) const {

  //approximation based on A. Mastichiadis et al.,
  //Astroph. Journ. 300:178-189 (1986), eq. 30.
  //This approx is valid only for   alpha >=100
  //where alpha = p0*eps*costheta - E0*eps;
  //for our purposes, me << E0 --> p0~ E0 -->
  //alpha = E0*eps*(costheta - 1) >= 100;

  int id = candidate->current.getId();
  double z = candidate->getRedshift();
  double E = candidate->current.getEnergy() * (1+z);
  double Epp = 0.;

  if (haveElectrons){
    // interpolate between tabulated electron energies to get corresponding cdf
    size_t i = std::upper_bound(tabE.begin(), tabE.end(), E) - tabE.begin() - 500; 
    double a = (E - tabE[i]) / (tabE[i + 500] - tabE[i]);

    std::vector<double> cdf(500);
    for (size_t j = 0; j < 500; j++)
      cdf[j] = tabCumulativeRate[i+j] + a * (tabCumulativeRate[i+500+j] - tabCumulativeRate[i+j]);

    // draw random value between 0. and maximum of corresponding cdf
    // choose bin of eps where cdf(eps) = cdf_rand -> eps_rand
    Random &random = Random::instance();
    size_t j = random.randBin(cdf); // draw random bin
    double binWidth = (tabs[i+j+1] - tabs[i+j]);
    double eps = (tabs[i+j] + random.rand() * binWidth)/4./E; // draw random uniform background photon energy in bin, note that one needs to convert tablulated eps back to s_kin and calculate right eps with real electron energy
    eps /= (1 + z); // dN/dE(Ep,Ee,z) = (1+z)^4 * dN/dE(Ep*(1+z),Ee*(1+z),0)  TODO: check if scaling needed 

    Epp = 5.7e-1 * pow(eps/eV, -0.56) * pow(E/eV, 0.44) * eV;
    
    candidate->addSecondary(11, Epp / (1+z));
    candidate->addSecondary(-11, Epp / (1+z));
    if (std::isfinite(Epp) == false)
      std::cout << "EMTripletPairProduction: " <<  Epp/eV << " " << E/eV << " " << tabE[0]/eV << " " << tabE[tabE.size() -1]/eV << " " << eps/eV << " " << i << " " << std::endl;
    if (Epp /eV > E/eV)
      std::cout << "EMTriplet: " << Epp / eV << " " << E / eV << std::endl; 
  }
  candidate->current.setEnergy((E - 2.*Epp)/(1+z));
}

void EMTripletPairProduction::process(Candidate *candidate) const {
	double step = candidate->getCurrentStep();
	double z = candidate->getRedshift();
	
  // check if electron / positron
  int id = candidate->current.getId();
  if (id != 11 && id != -11)
    return; 

  // instead of scaling the background photon energies, scale the electron energy
  double E = (1 + z) * candidate->current.getEnergy();

  // check if in tabulated energy range
  if (E < tabElectronEnergy.front() or (E > tabElectronEnergy.back()))
//  if (E < tabE.front() or (E > tabE.back()))
    return;

  // find interaction with minimum random distance
  Random &random = Random::instance();
  double randDistance = std::numeric_limits<double>::max();

  // comological scaling of interaction distance (comoving)
  double scaling = pow(1 + z, 2) * photonFieldScaling(photonField, z);
  double rate = scaling * interpolate(E, tabElectronEnergy, tabInteractionRate);
  randDistance = -log(random.rand()) / rate;

  // check if interaction does not happen
  if (step < randDistance) {
    candidate->limitNextStep(limit / rate);
    return;
  }

  // interact 
  performInteraction(candidate);
}

} // namespace crpropa
