#include "crpropa/module/SaveEM.h"
#include "crpropa/Units.h"

#include <iostream>
#include <sstream>
#include <stdexcept>

using namespace std;

namespace crpropa {

SaveEM::SaveEM(const string &filename) :
		filename(filename), output(filename.c_str()) {
	output << "#ID\tE\tD\tpID\tpE\tiID\tiE\n";
	output << "#\n";
	output << "# ID          Id of particle (photon, electron, positron)\n";
	output << "# E           Energy [EeV]\n";
	output << "# D           Comoving trajectory length [Mpc]\n";
	output << "# pID         Id of parent particle\n";
	output << "# pE          Energy [EeV] of parent particle\n";
	output << "# iID         Id of source particle\n";
	output << "# iE          Energy [EeV] of source particle\n";
	output << "#\n";
	Ethreshold = 0.1 * EeV;
}

SaveEM::SaveEM(const string &filename, double Ethr) :
		filename(filename), output(filename.c_str()) {
	output << "#ID\tE\tD\tpID\tpE\tiID\tiE\n";
	output << "#\n";
	output << "# ID          Id of particle (photon, electron, positron)\n";
	output << "# E           Energy [EeV]\n";
	output << "# D           Comoving trajectory length [Mpc]\n";
	output << "# pID         Id of parent particle\n";
	output << "# pE          Energy [EeV] of parent particle\n";
	output << "# iID         Id of source particle\n";
	output << "# iE          Energy [EeV] of source particle\n";
	output << "#\n";
	Ethreshold = Ethr;
}

SaveEM::~SaveEM() {
}

void SaveEM::process(Candidate *candidate) const {
	int pid = candidate->current.getId();
	if ((pid != 22) and (abs(pid) != 11))
		return;
	if (candidate->isActive()==false || candidate->current.getEnergy() > Ethreshold)
		return;

	char buffer[1024];
	size_t p = 0;

	p += sprintf(buffer + p, "%4i\t", pid);
	p += sprintf(buffer + p, "%g\t", candidate->current.getEnergy() / EeV);
	p += sprintf(buffer + p, "%8.4f\t", candidate->created.getPosition().getR() / Mpc);

	p += sprintf(buffer + p, "%10i\t", candidate->created.getId());
	p += sprintf(buffer + p, "%8.4f\t", candidate->created.getEnergy() / EeV);

	p += sprintf(buffer + p, "%10i\t", candidate->source.getId());
	p += sprintf(buffer + p, "%8.4f\n", candidate->source.getEnergy() / EeV);


#pragma omp critical
	{
		output.write(buffer, p);
	}

	candidate->setActive(false);
}

void SaveEM::endRun() {
	output.flush();
}

string SaveEM::getDescription() const {
	std::stringstream s;
	s << "SaveEM: Output file = " << filename;
	return s.str();
}

} // namespace crpropa
