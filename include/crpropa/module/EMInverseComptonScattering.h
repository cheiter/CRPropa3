#ifndef CRPROPA_EMINVERSECOMPTONSCATTERING_H
#define CRPROPA_EMINVERSECOMPTONSCATTERING_H

#include "crpropa/Module.h"
#include "crpropa/PhotonBackground.h"
#include <fstream>

namespace crpropa {

/**
 @class EMInverseComptonScattering
 @brief Inverse Compton scattering of electrons with background photons.

 This module simulates inverse Compton scattering in electron background photon interactions.\n
 Several photon fields can be selected.\n
 By default, the module limits the step size to 10% of the energy loss length of the particle.
 */
class EMInverseComptonScattering: public Module {
private:
	PhotonField photonField;

	std::vector<double> tabInteractionRate; /*< tabulated interaction rate in [1/m] */
	std::vector<double> tabElectronEnergy; /*< tabulated electron energy in [J] */
	std::vector<double> tabCumulativeRate; /*< tabulated cumulative interaction rate in [1/m] */
	std::vector<double> tabE; /*< tabulated electron energy in [J], 111 steps from 10^12 - 10^23 eV other stepsize than in tabElectronEnergy*/
	std::vector<double> tabs; /*< tabulated Mandelstam s in [J**2], 500 steps*/
	double limit; ///< fraction of energy loss length to limit the next step
	bool havePhotons;

public:
	EMInverseComptonScattering(PhotonField photonField = CMB, bool havePhotons =
			false, double limit = 0.1);

	void setPhotonField(PhotonField photonField);
	void setHavePhotons(bool havePhotons);
	void setLimit(double limit);

	void initRate(std::string filename);
	void initCumulativeRate(std::string filename);
	void process(Candidate *candidate) const;
	void performInteraction(Candidate *candidate) const;

};

} // namespace crpropa

#endif // CRPROPA_EMINVERSECOMPTONSCATTERING_H
