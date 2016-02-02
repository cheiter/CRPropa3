#include "crpropa/module/EMInverseComptonScattering.h"
#include "crpropa/Units.h"
#include "crpropa/ParticleID.h"
#include "crpropa/ParticleMass.h"
#include "crpropa/Random.h"

#include <fstream>
#include <limits>
#include <stdexcept>

namespace crpropa {

EMInverseComptonScattering::EMInverseComptonScattering(PhotonField photonField,
		bool havePhotons, double limit) {
	setPhotonField(photonField);
	this->havePhotons = havePhotons;
	this->limit = limit;
  out.open("/home/home1/institut_3a/heiter/Desktop/Energy_Secondary_Electrons_Photons_Directly_After_Interaction/data/CRPropa_ICS_electron.txt");
}

void EMInverseComptonScattering::setPhotonField(PhotonField photonField) {
	this->photonField = photonField;
	switch (photonField) {
	case CMB:
		setDescription("EMInverseComptonScattering: CMB");
		initRate(getDataPath("EMInverseComptonScattering_CMB.txt"));
		initCumulativeRate(getDataPath("EMInverseComptonScattering_CDF_CMB.txt"));
    initEleCaStuff(getDataPath("cdf_table_EleCa_CMB.txt"));
		break;
	case IRB:  // default: Kneiske '04 IRB model
	case IRB_Kneiske04:
		setDescription("EMInverseComptonScattering: IRB (Kneiske 2004)");
		initRate(getDataPath("EMInverseComptonScattering_IRB_Kneiske04.txt"));
		initCumulativeRate(getDataPath("EMInverseComptonScattering_CDF_IRB_Kneiske04.txt"));
    initEleCaStuff(getDataPath("cdf_table_EleCa_IRB.txt"));
		break;
	case IRB_Stecker05:
		setDescription("EMInverseComptonScattering: IRB (Stecker 2005)");
		initRate(getDataPath("EMInverseComptonScattering_IRB_Stecker05.txt"));
		initCumulativeRate(getDataPath("EMInverseComptonScattering_CDF_IRB_Stecker05.txt"));
    initEleCaStuff(getDataPath("cdf_table_EleCa_IRB.txt"));
		break;
	case IRB_Franceschini08:
		setDescription("EMInverseComptonScattering: IRB (Franceschini 2008)");
		initRate(getDataPath("EMInverseComptonScattering_IRB_Franceschini08.txt"));
		initCumulativeRate(getDataPath("EMInverseComptonScattering_CDF_IRB_Franceschini08.txt"));
    initEleCaStuff(getDataPath("cdf_table_EleCa_IRB.txt"));
		break;
	case IRB_Finke10:
		setDescription("EMInverseComptonScattering: IRB (Finke 2010)");
		initRate(getDataPath("EMInverseComptonScattering_IRB_Finke10.txt"));
		initCumulativeRate(getDataPath("EMInverseComptonScattering_CDF_IRB_Finke10.txt"));
    initEleCaStuff(getDataPath("cdf_table_EleCa_IRB.txt"));
		break;
	case IRB_Dominguez11:
		setDescription("EMInverseComptonScattering: IRB (Dominguez 2011)");
		initRate(getDataPath("EMInverseComptonScattering_IRB_Dominguez11.txt"));
		initCumulativeRate(getDataPath("EMInverseComptonScattering_CDF_IRB_Dominguez11.txt"));
    initEleCaStuff(getDataPath("cdf_table_EleCa_IRB.txt"));
		break;
	case IRB_Gilmore12:
		setDescription("EMInverseComptonScattering: IRB (Gilmore 2012)");
		initRate(getDataPath("EMInverseComptonScattering_IRB_Gilmore12.txt"));
		initCumulativeRate(getDataPath("EMInverseComptonScattering_CDF_IRB_Gilmore12.txt"));
    initEleCaStuff(getDataPath("cdf_table_EleCa_IRB.txt"));
		break;
	case URB_Protheroe96:
		setDescription("EMInverseComptonScattering: URB (Protheroe 1996)");
		initRate(getDataPath("EMInverseComptonScattering_URB_Protheroe96.txt"));
		initCumulativeRate(getDataPath("EMInverseComptonScattering_CDF_URB_Protheroe96.txt"));
    initEleCaStuff(getDataPath("cdf_table_EleCa_URB.txt"));
		break;
	default:
		throw std::runtime_error(
				"EMInverseComptonScattering: unknown photon background");
	}
}

void EMInverseComptonScattering::setHavePhotons(bool havePhotons) {
	this->havePhotons = havePhotons;
}

void EMInverseComptonScattering::setLimit(double limit) {
	this->limit = limit;
}

void EMInverseComptonScattering::initRate(std::string filename) {
	std::ifstream infile(filename.c_str());

	if (!infile.good())
		throw std::runtime_error(
				"EMInverseComptonScattering: could not open file " + filename);

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

void EMInverseComptonScattering::initCumulativeRate(std::string filename) {
	std::ifstream infile(filename.c_str());

	if (!infile.good())
		throw std::runtime_error(
				"EMInverseComptonScattering: could not open file " + filename);

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

void EMInverseComptonScattering::initEleCaStuff(std::string filename) {
	std::ifstream infile(filename.c_str());

	if (!infile.good())
		throw std::runtime_error(
				"EMInverseComptonScattering: could not open file " + filename);

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
///Differential cross-section for inverse Compton scattering. from Lee, eq. 23 for x = Ee'/Ee
double dSigmadE_ICSx(double x, double beta) {
  double q = ((1 - beta) / beta) * (1 - 1./x);
  double A = x + 1./x;
  return ((1 + beta) / beta) * (A + 2 * q + q * q);
}

/// Hold an data array to interpolate the energy distribution on 
class ICSSecondariesEnergyDistribution
{
	private:
		double *_data;
		size_t _Ns;
		size_t _Nrer;
		double _s_min;
		double _s_max;
		double _dls;

	public:
		ICSSecondariesEnergyDistribution(double s_min = 1.01 * mass_electron*c_squared * mass_electron*c_squared, double s_max =1e23,
				size_t Ns = 1000, size_t Nrer = 1000 )
		{
			// TODO: this boundary is just an estimate
			const double l = 1.01;
			if (s_min < l * mass_electron*c_squared * mass_electron*c_squared)
			{
				std::cerr << "Warning: Minimum COM Energy in ICS Interpolation s = " << s_min << " < " << l << " m_e**2 selected. Setting to s_min = " << l << " m_e**2.\n" ;
				s_min = l * mass_electron*c_squared * mass_electron*c_squared;
			}
			_Ns = Ns;
			_Nrer = Nrer;
			_s_min =s_min;
			_s_max = s_max;
			_data = new double[Ns*Nrer];
			_dls = (log(s_max) - log(s_min)) / (Ns);
      double ElectronMass = mass_electron*c_squared;
			for (size_t i = 0; i < Ns; i++)
			{
				const double s = s_min * exp(i*_dls);
        double beta = (s - ElectronMass * ElectronMass) / (s + ElectronMass * ElectronMass);
        double x0 = log((1.-beta) / (1.+beta));
        double dx = -log((1. - beta)/(1.+beta)) / (Nrer);
				_data[i * Nrer] = dSigmadE_ICSx(exp(x0), beta); 
				for (size_t j = 1; j < Nrer; j++)
				{
					double x = exp(x0 + j*dx); 
					_data[i * Nrer + j] =	dSigmadE_ICSx(x, beta) + _data[i * Nrer + j - 1];
				}
			}
		}

		// returns pointer to the the integrated distribution for a given s
		double* getDistribution(double s)
		{
			size_t idx = (log(s / _s_min)) / _dls;
			double *s0 = &_data[idx * _Nrer];
			return s0;
		}

		//samples the integrated distribution and returns Eer(Ee, s)
		double sample(double Ee, double s)
		{
      double ElectronMass = mass_electron*c_squared;
			double *s0 = getDistribution(s); 
      Random &random = Random::instance();
			double rnd = random.rand() *s0[_Nrer-1];
			for (size_t i=0; i < _Nrer; i++)
			{
				if (rnd < s0[i])
				{
					double beta = (s - ElectronMass * ElectronMass) / (s + ElectronMass * ElectronMass);
					double x0 = log((1-beta) / (1+beta));
					double dx =  - log((1-beta) / (1+beta)) / (_Nrer );
					return exp(x0 + i*dx) * Ee; 
				}
			}
			throw std::runtime_error("Grave logic error in sampling ICSSecondariesEnergyDistribution!");	
		}
};


// Helper function for actual Monte Carlo sampling to avoid code-duplication
double __extractICSSecondaries(double Ee, double s)
{
	static ICSSecondariesEnergyDistribution interpolation;
	return interpolation.sample(Ee, s);
}


void EMInverseComptonScattering::performInteraction(Candidate *candidate) const {

  int id = candidate->current.getId();
  double z = candidate->getRedshift();
  double E = candidate->current.getEnergy();
  double Epost = 0.;

  //    // interpolate between tabulated electron energies to get corresponding cdf
  //    size_t i = std::upper_bound(tabE.begin(), tabE.end(), E) - tabE.begin() - 500; 
  //    double a = (E - tabE[i]) / (tabE[i + 500] - tabE[i]);
  //
  //    std::vector<double> cdf(500);
  //    for (size_t j = 0; j < 500; j++)
  //      cdf[j] = tabCumulativeRate[i+j] + a * (tabCumulativeRate[i+500+j] - tabCumulativeRate[i+j]);
  //
  //    // draw random value between 0. and maximum of corresponding cdf
  //    // choose bin of s where cdf(s) = cdf_rand -> s_rand
  //    Random &random = Random::instance();
  //    size_t j = random.randBin(cdf); // draw random bin
  //    double binWidth = (tabs[i+j+1] - tabs[i+j]);
  //    double s_kin = tabs[i+j] + random.rand() * binWidth; // draw random s uniformly distributed in bin
  //    double s = s_kin + (mass_electron*c_squared)*(mass_electron*c_squared);
  //    s /= (1 + z) * (1 + z); // dN/dE(Ep,Ee,z) = (1+z)^4 * dN/dE(Ep*(1+z),Ee*(1+z),0) TODO: check if scaling needed

  //EleCa method: 
  // draw random value between 0. and maximum of corresponding cdf
  // choose bin of s where cdf(s) = cdf_rand -> s_rand
  double mec2 = mass_electron * c_squared;
  Random &random = Random::instance();
  double eps = 0.;
  double epsMin = 1. * mec2 * mec2 / 4. / E; // Minimum neccessary eps to have sufficient value of Mandelstam s for interaction process
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
  eps += random.rand() * binWidth; // draw random eps uniformly distributed in bin

  //    Random &random = Random::instance();
  //    size_t j = random.randBin(tabCDF); // draw random bin
  //    double binWidth = (tabEps[j+1] - tabEps[j]);
  //    double eps = tabEps[j] + random.rand() * binWidth; // draw random s uniformly distributed in bin
  if (eps < epsMin)  //TODO: Abbruchbedingung interaction kann nciht stattfinden mit diesem eps, vor oder nach scaling + ist das ok oder muss anderes eps gewählt werden ??
    return;
  eps *= (1.+z);
  double s = 4.*E*eps + (mass_electron*c_squared)*(mass_electron*c_squared);

//  s = 4.*E*1e-3 * eV;
  Epost = __extractICSSecondaries(E,s); 
  if (id ==11)
    out << Epost / eV << "\n";
  if (havePhotons)
    candidate->addSecondary(22, (E-Epost));
  //    if (std::isfinite(Epost) == false)
  //      std::cout << "EMInverseComptonScattering: " << Epost/eV << " " << E/eV << " " << tabE[0]/eV << " " << tabE[tabE.size() -1]/eV << " " << s << " " << i << " " << std::endl;
  candidate->current.setEnergy(Epost); 
}

void EMInverseComptonScattering::process(Candidate *candidate) const {
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
    return;

  // find interaction with minimum random distance
  Random &random = Random::instance();
  double randDistance = std::numeric_limits<double>::max();

  // cosmological scaling of interaction distance (comoving)
  double scaling = pow(1 + z, 2) * photonFieldScaling(photonField, z);
  double rate = scaling * interpolate(E, tabElectronEnergy, tabInteractionRate);
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
