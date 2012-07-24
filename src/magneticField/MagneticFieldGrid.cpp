#include "mpc/magneticField/MagneticFieldGrid.h"
#include "mpc/Random.h"
#include <fstream>

namespace mpc {

MagneticFieldGrid::MagneticFieldGrid(VectorFieldGrid *grid) {
	setGrid(grid);
}

void MagneticFieldGrid::setGrid(VectorFieldGrid *grid) {
	this->grid = grid;
}

VectorFieldGrid *MagneticFieldGrid::getGrid() {
	return grid;
}

Vector3d MagneticFieldGrid::getField(const Vector3d &pos) const {
	return grid->interpolate(pos);
}

void ModulatedMagneticFieldGrid::setGrid(VectorFieldGrid *g) {
	grid = g;
}

VectorFieldGrid *ModulatedMagneticFieldGrid::getGrid() {
	return grid;
}

void ModulatedMagneticFieldGrid::setModulationGrid(ScalarFieldGrid *g) {
	modGrid = g;
}

ScalarFieldGrid *ModulatedMagneticFieldGrid::getModulationGrid() {
	return modGrid;
}

Vector3d ModulatedMagneticFieldGrid::getField(const Vector3d &pos) const {
	float m = modGrid->interpolate(pos);
	Vector3d b = grid->interpolate(pos);
	return b * m;
}

Vector3f meanFieldStrength(VectorFieldGrid *m) {
	size_t Nx = m->getNx();
	size_t Ny = m->getNy();
	size_t Nz = m->getNz();
	Vector3f mean(0.);
	for (int ix = 0; ix < Nx; ix++)
		for (int iy = 0; iy < Ny; iy++)
			for (int iz = 0; iz < Nz; iz++)
				mean += m->get(ix, iy, iz);
	return mean / Nx / Ny / Nz;
}

double rmsFieldStrength(VectorFieldGrid *m) {
	size_t Nx = m->getNx();
	size_t Ny = m->getNy();
	size_t Nz = m->getNz();
	double sumB2 = 0;
	for (int ix = 0; ix < Nx; ix++)
		for (int iy = 0; iy < Ny; iy++)
			for (int iz = 0; iz < Nz; iz++)
				sumB2 += m->get(ix, iy, iz).getMag2();
	return sqrt(sumB2 / Nx / Ny / Nz);
}

void scale(VectorFieldGrid *m, double a) {
	for (int ix = 0; ix < m->getNx(); ix++)
		for (int iy = 0; iy < m->getNy(); iy++)
			for (int iz = 0; iz < m->getNz(); iz++)
				m->get(ix, iy, iz) *= a;
}

#ifdef MPC_HAVE_FFTW3F
#include "fftw3.h"

void initTurbulence(VectorFieldGrid *m, double Brms, double lMin, double lMax,
		double spectralIndex, int seed) {
	size_t Nx = m->getNx();
	size_t Ny = m->getNy();
	size_t Nz = m->getNz();
	if ((Nx != Ny) or (Ny != Nz))
	throw std::runtime_error("turbulentField: only cubic grid supported");

	double spacing = m->getSpacing();
	if (lMin < 2 * spacing)
	throw std::runtime_error("turbulentField: lMin < 2 * spacing");
	if (lMin >= lMax)
	throw std::runtime_error("turbulentField: lMin >= lMax");
	if (lMax > Nx * spacing / 2)
	throw std::runtime_error("turbulentField: lMax > size / 2");

	size_t n = Nx; // size of array
	size_t n2 = (size_t) floor(n / 2) + 1;// size array in z-direction in configuration space

	// arrays to hold the complex vector components of the B(k)-field
	fftwf_complex *Bkx, *Bky, *Bkz;
	Bkx = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * n * n * n2);
	Bky = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * n * n * n2);
	Bkz = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * n * n * n2);

	Random random;
	if (seed != 0)
	random.seed(seed);// use given seed

	// calculate the n possible discrete wave numbers
	double K[n];
	for (int i = 0; i < n; i++)
	K[i] = (double) i / n - i / (n / 2);

	// construct the field in configuration space
	int i;
	double k, theta, phase, cosPhase, sinPhase;
	double kMin = spacing / lMax;
	double kMax = spacing / lMin;
	Vector3f b;// real b-field vector
	Vector3f ek, e1, e2;// orthogonal base
	Vector3f n0(1, 1, 1);// arbitrary vector to construct orthogonal base

	for (size_t ix = 0; ix < n; ix++) {
		for (size_t iy = 0; iy < n; iy++) {
			for (size_t iz = 0; iz < n2; iz++) {

				i = ix * n * n2 + iy * n2 + iz;
				ek.setXYZ(K[ix], K[iy], K[iz]);
				k = ek.getMag();

				// wave outside of turbulent range -> B(k) = 0
				if ((k < kMin) || (k > kMax)) {
					Bkx[i][0] = 0;
					Bkx[i][1] = 0;
					Bky[i][0] = 0;
					Bky[i][1] = 0;
					Bkz[i][0] = 0;
					Bkz[i][1] = 0;
					continue;
				}

				// construct an orthogonal base ek, e1, e2
				if (ek.getAngleTo(n0) < 1e-3) {
					// ek parallel to (1,1,1)
					e1.setXYZ(-1., 1., 0);
					e2.setXYZ(1., 1., -2.);
				} else {
					// ek not parallel to (1,1,1)
					e1 = n0.cross(ek);
					e2 = ek.cross(e1);
				}
				e1 /= e1.getMag();
				e2 /= e2.getMag();

				// random orientation perpendicular to k
				theta = 2 * M_PI * random.rand();
				b = e1 * cos(theta) + e2 * sin(theta);

				// standard normal distributed amplitude weighted with k^alpha/2
				b *= random.randNorm() * pow(k, spectralIndex / 2.);

				// uniform random phase
				phase = 2 * M_PI * random.rand();
				cosPhase = cos(phase);// real part
				sinPhase = sin(phase);// imaginary part

				Bkx[i][0] = b.x * cosPhase;
				Bkx[i][1] = b.x * sinPhase;
				Bky[i][0] = b.y * cosPhase;
				Bky[i][1] = b.y * sinPhase;
				Bkz[i][0] = b.z * cosPhase;
				Bkz[i][1] = b.z * sinPhase;
			}
		}
	}

	// in-place, complex to real, inverse Fourier transformation on each component
	// note that the last elements of B(x) are unused now
	float *Bx = (float*) Bkx;
	fftwf_plan plan_x = fftwf_plan_dft_c2r_3d(n, n, n, Bkx, Bx, FFTW_ESTIMATE);
	fftwf_execute(plan_x);
	fftwf_destroy_plan(plan_x);

	float *By = (float*) Bky;
	fftwf_plan plan_y = fftwf_plan_dft_c2r_3d(n, n, n, Bky, By, FFTW_ESTIMATE);
	fftwf_execute(plan_y);
	fftwf_destroy_plan(plan_y);

	float *Bz = (float*) Bkz;
	fftwf_plan plan_z = fftwf_plan_dft_c2r_3d(n, n, n, Bkz, Bz, FFTW_ESTIMATE);
	fftwf_execute(plan_z);
	fftwf_destroy_plan(plan_z);

	// save to grid
	for (size_t ix = 0; ix < n; ix++) {
		for (size_t iy = 0; iy < n; iy++) {
			for (size_t iz = 0; iz < n; iz++) {
				i = ix * n * 2 * n2 + iy * 2 * n2 + iz;
				Vector3f &b = m->get(ix, iy, iz);
				b.x = Bx[i];
				b.y = By[i];
				b.z = Bz[i];
			}
		}
	}

	fftwf_free(Bkx);
	fftwf_free(Bky);
	fftwf_free(Bkz);

	scale(m, Brms / rmsFieldStrength(m)); // normalize to Brms
}

double turbulentCorrelationLength(double lMin, double lMax, double spectralIndex) {
	double r = lMin / lMax;
	double a = -spectralIndex - 2;
	return lMax / 2 * (a - 1) / a * (1 - pow(r, a)) / (1 - pow(r, a - 1));
}
#endif // MPC_HAVE_FFTW3F

void load(VectorFieldGrid *m, std::string filename, double c) {
	std::ifstream fin(filename.c_str(), std::ios::binary);
	if (!fin)
		throw std::runtime_error("MagneticFieldGrid: File not found");
	for (int ix = 0; ix < m->getNx(); ix++) {
		for (int iy = 0; iy < m->getNy(); iy++) {
			for (int iz = 0; iz < m->getNz(); iz++) {
				Vector3f &b = m->get(ix, iy, iz);
				fin.read((char*) &(b.x), sizeof(float));
				fin.read((char*) &(b.y), sizeof(float));
				fin.read((char*) &(b.z), sizeof(float));
				b *= c;
			}
		}
	}
	fin.close();
}

void dump(VectorFieldGrid *m, std::string filename, double c) {
	std::ofstream fout(filename.c_str(), std::ios::binary);
	if (!fout)
		throw std::runtime_error("MagneticFieldGrid: Could not open file");
	for (int ix = 0; ix < m->getNx(); ix++) {
		for (int iy = 0; iy < m->getNy(); iy++) {
			for (int iz = 0; iz < m->getNz(); iz++) {
				Vector3f b = m->get(ix, iy, iz) * c;
				fout.write((char*) &(b.x), sizeof(float));
				fout.write((char*) &(b.y), sizeof(float));
				fout.write((char*) &(b.z), sizeof(float));
			}
		}
	}
	fout.close();
}

void loadTxt(VectorFieldGrid *m, std::string filename, double unit) {
	std::ifstream fin(filename.c_str());
	if (!fin)
		throw std::runtime_error("MagneticFieldGrid: file not found");

	// skip header lines
	while (fin.peek() == '#')
		fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

	for (int ix = 0; ix < m->getNx(); ix++) {
		for (int iy = 0; iy < m->getNy(); iy++) {
			for (int iz = 0; iz < m->getNz(); iz++) {
				Vector3f &b = m->get(ix, iy, iz);
				fin >> b.x >> b.y >> b.z;
				b *= unit;
				if (fin.eof())
					throw std::runtime_error(
							"MagneticFieldGrid: file too short");
			}
		}
	}
	fin.close();
}

void dumpTxt(VectorFieldGrid *m, std::string filename, double unit) {
	std::ofstream fout(filename.c_str());
	if (!fout)
		throw std::runtime_error("MagneticFieldGrid: could not open file");
	for (int ix = 0; ix < m->getNx(); ix++) {
		for (int iy = 0; iy < m->getNy(); iy++) {
			for (int iz = 0; iz < m->getNz(); iz++) {
				Vector3f b = m->get(ix, iy, iz) * unit;
				fout << b << "\n";
			}
		}
	}
	fout.close();
}

} // namespace mpc
