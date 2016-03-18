#ifndef CRPROPA_SAVE_EM_H
#define CRPROPA_SAVE_EM_H

#include "crpropa/Module.h"
#include "crpropa/Units.h"

#include <memory>
#include <fstream>

namespace crpropa {

class SaveEM: public Module {
private:
	std::string filename;
	mutable std::ofstream output;
	mutable double Ethreshold;
public:
	SaveEM(const std::string &filename);
	SaveEM(const std::string &filename, const double Ethr);
	~SaveEM();
	void process(Candidate *candidate) const;
	std::string getDescription() const;
	void endRun();
};

} // namespace crpropa

#endif // CRPROPA_SAVE_EM_H
