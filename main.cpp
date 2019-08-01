#include "pp.h"

int main()
{
	PP *pp = new PP;

	// number:
	// 0 - quark-diquark density
	// 1 - three-quarks density
	// 2 - mixed-density
	// 3 - three-quarks-one-gluon-body density
	// 6 - triangular density
	// 7 - quark-diquark trianfular density

	int number = 6;
	double kappa = 0.5;
	
	pp->performCollisions(number, kappa);

	delete pp;

    return 0;
}
