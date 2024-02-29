#include "PlanarGeometry.h"
#include <cmath>

// Restrict vector to parallel epipede spanned by pbc. 
// Maximum coefficient: +-0.5 in any direction after restriction
v2d restrict(const v2d vector, const std::pair<v2d, v2d> &pbc){
	const v2d coeff = represent(vector, pbc);

	double b1_coefficient = coeff.real();
	double b2_coefficient = coeff.imag();

	double rescaled_b1 = b1_coefficient - floor(b1_coefficient + 0.5);
	double rescaled_b2 = b2_coefficient - floor(b2_coefficient + 0.5);

	return pbc.first*rescaled_b1 + pbc.second*rescaled_b2;
}

// Does a1->b1 intersect segment a2->b2
bool intersectsSegment(v2d a1, v2d b1, v2d a2, v2d b2) {

	v2d dir1 = b1 - a1;
	v2d a2_ = (a2 - a1) / dir1;
	v2d b2_ = (b2 - a1) / dir1;


	// Both a2 and b2 on the same side of the a1 b1line
	if (a2_.imag()*b2_.imag() > 0) {
		return false;
	}

	// signs of imaginary parts is opposite, so lets take the absolute for interpolation
	double ia_abs = std::abs(a2_.imag());
	double ib_abs = std::abs(b2_.imag());
	double ra = a2_.real();
	double rb = b2_.real();
	/*std::cerr << dir1.real() << "\t" << dir1.imag() << std::endl;
	std::cerr << a2_.real() << "\t" << a2_.imag() << std::endl;
	std::cerr << b2_.real() << "\t" << b2_.imag() << std::endl;
	
	std::cerr << std::endl;
	std::cerr << ra << "\t" << ia_abs << std::endl;
	std::cerr << rb << "\t" << ib_abs << std::endl;*/
	
	// Calculate position of a2_->b2_ intersecting new x-axis
	double cutposx = (ra*ib_abs + rb*ia_abs) / (ia_abs + ib_abs);
	//std::cerr << "cutposx="<< cutposx << std::endl;

	return 0.0 <= cutposx && cutposx <= 1.0;
}

// Represent vector in terms of basis given by base
v2d represent(const v2d vector, const std::pair<v2d, v2d> &base) {
	double x = vector.real();
	double y = vector.imag();

	const v2d &base1 = base.first;
	const v2d &base2 = base.second;
	double a = base1.real();
	double c = base1.imag();
	double b = base2.real();
	double d = base2.imag();

	double b1_coefficient = (d*x - b*y) / (a*d - b*c);
	double b2_coefficient = (-c*x + a*y) / (a*d - b*c);
	return v2d(b1_coefficient, b2_coefficient);
}
