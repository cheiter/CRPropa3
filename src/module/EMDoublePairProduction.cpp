#include "crpropa/module/EMDoublePairProduction.h"
#include "crpropa/Units.h"
#include "crpropa/ParticleID.h"
#include "crpropa/ParticleMass.h"
#include "crpropa/Random.h"

#include <fstream>
#include <limits>
#include <stdexcept>

namespace crpropa {

EMDoublePairProduction::EMDoublePairProduction(PhotonField photonField,
		bool haveElectrons, double limit) {
	setPhotonField(photonField);
	this->haveElectrons = haveElectrons;
	this->limit = limit;
}

void EMDoublePairProduction::setPhotonField(PhotonField photonField) {
	this->photonField = photonField;
	switch (photonField) {
	case CMB:
		setDescription("EMDoublePairProduction: CMB");
		initRate(getDataPath("EMDoublePairProduction_CMB.txt"));
    initEleCaStuff(getDataPath("cdf_table_EleCa_CMB.txt"));
//    initEleCaStuff(getDataPath("cdf_table_EleCa_All.txt"));
		break;
	case IRB:  // default: Kneiske '04 IRB model
	case IRB_Kneiske04:
		setDescription("EMDoublePairProduction: IRB (Kneiske 2004)");
		initRate(getDataPath("EMDoublePairProduction_IRB_Kneiske04.txt"));
    initEleCaStuff(getDataPath("cdf_table_EleCa_IRB.txt"));
//    initEleCaStuff(getDataPath("cdf_table_EleCa_All.txt"));
		break;
	case IRB_Stecker05:
		setDescription("EMDoublePairProduction: IRB (Stecker 2005)");
		initRate(getDataPath("EMDoublePairProduction_IRB_Stecker05.txt"));
    initEleCaStuff(getDataPath("cdf_table_EleCa_IRB.txt"));
//    initEleCaStuff(getDataPath("cdf_table_EleCa_All.txt"));
		break;
	case IRB_Franceschini08:
		setDescription("EMDoublePairProduction: IRB (Franceschini 2008)");
		initRate(getDataPath("EMDoublePairProduction_IRB_Franceschini08.txt"));
    initEleCaStuff(getDataPath("cdf_table_EleCa_IRB.txt"));
//    initEleCaStuff(getDataPath("cdf_table_EleCa_All.txt"));
		break;
	case IRB_Finke10:
		setDescription("EMDoublePairProduction: IRB (Finke 2010)");
		initRate(getDataPath("EMDoublePairProduction_IRB_Finke10.txt"));
    initEleCaStuff(getDataPath("cdf_table_EleCa_IRB.txt"));
//    initEleCaStuff(getDataPath("cdf_table_EleCa_All.txt"));
		break;
	case IRB_Dominguez11:
		setDescription("EMDoublePairProduction: IRB (Dominguez 2011)");
		initRate(getDataPath("EMDoublePairProduction_IRB_Dominguez11.txt"));
    initEleCaStuff(getDataPath("cdf_table_EleCa_IRB.txt"));
//    initEleCaStuff(getDataPath("cdf_table_EleCa_All.txt"));
		break;
	case IRB_Gilmore12:
		setDescription("EMDoublePairProduction: IRB (Gilmore 2012)");
		initRate(getDataPath("EMDoublePairProduction_IRB_Gilmore12.txt"));
    initEleCaStuff(getDataPath("cdf_table_EleCa_IRB.txt"));
//    initEleCaStuff(getDataPath("cdf_table_EleCa_All.txt"));
		break;
	case URB_Protheroe96:
		setDescription("EMDoublePairProduction: URB (Protheroe 1996)");
		initRate(getDataPath("EMDoublePairProduction_URB_Protheroe96.txt"));
    initEleCaStuff(getDataPath("cdf_table_EleCa_URB.txt"));
//    initEleCaStuff(getDataPath("cdf_table_EleCa_All.txt"));
		break;
	default:
		throw std::runtime_error(
				"EMDoublePairProduction: unknown photon background");
	}
}

void EMDoublePairProduction::setHaveElectrons(bool haveElectrons) {
	this->haveElectrons = haveElectrons;
}

void EMDoublePairProduction::setLimit(double limit) {
	this->limit = limit;
}

void EMDoublePairProduction::initRate(std::string filename) {
	std::ifstream infile(filename.c_str());

	if (!infile.good())
		throw std::runtime_error(
				"EMDoublePairProduction: could not open file " + filename);

	// clear previously loaded interaction rates
	tabPhotonEnergy.clear();
	tabInteractionRate.clear();

	while (infile.good()) {
		if (infile.peek() != '#') {
			double a, b;
			infile >> a >> b;
			if (infile) {
				tabPhotonEnergy.push_back(pow(10, a) * eV);
				tabInteractionRate.push_back(b / Mpc);
			}
		}
		infile.ignore(std::numeric_limits < std::streamsize > ::max(), '\n');
	}
	infile.close();
}

void EMDoublePairProduction::initEleCaStuff(std::string filename) {
	std::ifstream infile(filename.c_str());

	if (!infile.good())
		throw std::runtime_error(
				"EMPairProduction: could not open file " + filename);

	// clear previously loaded interaction rates
	tabEps.clear();
  tabCDF.clear();

	while (infile.good()) {
		if (infile.peek() != '#') {
			double a, b;
			infile >> a >> b;
			if (infile) {
				tabEps.push_back(a*eV);
        tabCDF.push_back(b);
			}
		}
		infile.ignore(std::numeric_limits < std::streamsize > ::max(), '\n');
	}
	infile.close();
}

void EMDoublePairProduction::performInteraction(Candidate *candidate) const {
  double E = candidate->current.getEnergy();
  double mec2 = mass_electron * c_squared;
  
	//EleCa Method:
  Random &random = Random::instance();
  double eps = 0.;
  double epsMin = 16. * mec2 * mec2 / 4. / E; // Minimum neccessary eps to have sufficient value of Mandelstam s for interaction process
  std::vector<double>::const_iterator it;
  it = std::lower_bound(tabEps.begin(), tabEps.end(), epsMin);
  size_t iE;
  if (it == tabEps.begin())
    iE = 0;
  else if (it == tabEps.end())
    iE = tabEps.size() - 1;
  else
    iE = it - tabEps.begin();
  double h = random.rand() * (1-tabCDF[iE]) + tabCDF[iE];
  it = std::upper_bound(tabCDF.begin(), tabCDF.end(), h);
  if (it == tabCDF.begin())
    eps = tabEps.front();
  else if (it == tabCDF.end())
    eps = tabEps.back();
  else
    eps =  tabEps[it - tabCDF.begin()];
  double binWidth = (tabEps[it-tabCDF.begin()+1] - tabEps[it-tabCDF.begin()]);
//  eps += random.rand() * binWidth; // draw random eps uniformly distributed in bin

//    Random &random = Random::instance();
//    size_t j = random.randBin(tabCDF); // draw random bin
//    double binWidth = (tabEps[j+1] - tabEps[j]);
//    double eps = tabEps[j] + random.rand() * binWidth; // draw random eps uniformly distributed in bin
//    double Eps = tabEps[j] + random.rand() * binWidth; // draw random s uniformly distributed in bin
  if (eps < epsMin)  //TODO: Abbruchbedingung interaction kann nciht stattfinden mit diesem eps, vor oder nach scaling + ist das ok oder muss anderes eps gewählt werden ??
    return;

	if (haveElectrons){
		double Ee = (E-2.*mass_electron*c_squared)/2.; // Use assumption of Lee 96 (i.e., all the energy goes equaly shared between only 1 couple of e+e- but take mass of second e+e- pair into account. In DPPpaper has been shown that this approximation is valid within -1.5%
		Vector3d pos = randomPositionInPropagationStep(candidate);
		candidate->addSecondary(11, Ee, pos);
		candidate->addSecondary(-11, Ee, pos);
	}
	candidate->setActive(false);
}

void EMDoublePairProduction::process(Candidate *candidate) const {
	double step = candidate->getCurrentStep();
	double z = candidate->getRedshift();

	// check if photon
	int id = candidate->current.getId();
	if (id != 22)
		return;

	// instead of scaling the background photon energies, scale the photon energy
	double E = (1 + z) * candidate->current.getEnergy();

	// check if in tabulated energy range
	if (E < tabPhotonEnergy.front() or (E > tabPhotonEnergy.back()))
		return;

	// find interaction with minimum random distance
	Random &random = Random::instance();
	double randDistance = std::numeric_limits<double>::max();

	// comological scaling of interaction distance (comoving)
	double scaling = pow(1 + z, 3) * photonFieldScaling(photonField, z);
	double rate = scaling * interpolate(E, tabPhotonEnergy, tabInteractionRate);
	randDistance = -log(random.rand()) / rate;

	candidate->limitNextStep(limit / rate);
	// check if interaction does not happen
	if (step < randDistance) {
		return;
	}

	// interact
	performInteraction(candidate);
}

} // namespace crpropa
