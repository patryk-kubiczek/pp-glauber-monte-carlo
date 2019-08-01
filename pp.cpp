#include "pp.h"

using namespace std;

inline double event_sigma(double N_c, double N_g, double b) // per b
{
    return (N_c > N_g * N_g ? 2 * M_PI * b : 2 * M_PI * b * (1 - pow(1 - N_c / (N_g * N_g), N_g * N_g)));
}

inline double gaussian(double x, double y, double r)
{
    return 1 / (2 * M_PI * r * r) * exp(- (x * x + y * y) / (2 * r * r));
}

inline double gaussian_tube(double x, double y, double theta, double phi, double r_s, double r_l)
{
    double det = r_s * r_s * (r_s * r_s + r_l * r_l + (r_s * r_s - r_l * r_l) * cos(2 * theta)) / 2.;
    double A = 0.25 * (3 * r_s * r_s + r_l * r_l + (r_s * r_s - r_l * r_l) * (cos(2 * theta) + 2 * cos(2 * phi) * sin(theta) * sin(theta))) / det;
    double B = (r_s * r_s - r_l * r_l) * sin(theta) * sin(theta) * cos(phi) * sin(phi) / det;
    double C = 0.25 * (3 * r_s * r_s + r_l * r_l + (r_s * r_s - r_l * r_l) * (cos(2 * theta) - 2 * cos(2 * phi) * sin(theta) * sin(theta))) / det;

    return 1. / (2 * M_PI * sqrt(det)) * exp(- 1 * (A * x * x + 2 * B * x * y + C * y * y) / 2.);
}

inline double three_quarks_thickness(double x, double y, double (&R)[9], double r_q, double N_g, double kappa)
{
    double R1 = sqrt(R[0] * R[0] + R[1] * R[1] + R[2] * R[2]);
    double R2 = sqrt(R[3] * R[3] + R[4] * R[4] + R[5] * R[5]);
    double R3 = sqrt(R[6] * R[6] + R[7] * R[7] + R[8] * R[8]);
    double theta1 = (R1 == 0 ? 0 : acos(R[2] / R1));
    double theta2 = (R2 == 0 ? 0 : acos(R[5] / R2));
    double theta3 = (R3 == 0 ? 0 : acos(R[8] / R3));
    double phi1 = (R[0] == 0 ? M_PI / 2. : atan2(R[1], R[0]));
    double phi2 = (R[3] == 0 ? M_PI / 2. : atan2(R[4], R[3]));
    double phi3 = (R[6] == 0 ? M_PI / 2. : atan2(R[7], R[6]));

    return (1 - kappa) * (N_g / 3.) * (gaussian(x - R[0], y - R[1], r_q)
                                     + gaussian(x - R[3], y - R[4], r_q)
                                     + gaussian(x - R[6], y - R[7], r_q))
                + kappa * (N_g / 3.) * (gaussian_tube(x - R[0] / 2., y - R[1] / 2., theta1, phi1, r_q, R1 / 2.)
                                      + gaussian_tube(x - R[3] / 2., y - R[4] / 2., theta2, phi2, r_q, R2 / 2.)
                                      + gaussian_tube(x - R[6] / 2., y - R[7] / 2., theta3, phi3, r_q, R3 / 2.));
//    + kappa * N_g * gaussian(x, y, sqrt(0.5 * 0.5 - r_q * r_q));
}

inline double three_quarks_one_gluon_body_thickness(double x, double y, double (&R)[9], double r_q, double N_g, double kappa)
{
    return (1 - kappa) * (N_g / 3.) * (gaussian(x - R[0], y - R[1], r_q)
                                     + gaussian(x - R[3], y - R[4], r_q)
                                     + gaussian(x - R[6], y - R[7], r_q))
           + kappa * N_g * gaussian(x, y, sqrt(0.5 * 0.5 - r_q * r_q));
}

inline double quark_diquark_thickness(double x, double y, double theta, double phi, double r_q, double d, double N_g, double kappa)
{
    return (1 - kappa) * (N_g / 2.) * (gaussian(x - d / 2. * cos(phi) * sin(theta), y - d / 2. * sin(phi) * sin(theta), r_q)
                                     + gaussian(x + d / 2. * cos(phi) * sin(theta), y + d / 2. * sin(phi) * sin(theta), r_q))
          + (kappa * N_g / 2.) * (gaussian_tube(x - d / 4. * cos(phi) * sin(theta),  y - d / 4. * sin(phi) * sin(theta), theta, phi, r_q, d / 4.) 
								+ gaussian_tube(x + d / 4. * cos(phi) * sin(theta),  y + d / 4. * sin(phi) * sin(theta), theta, phi, r_q, d / 4.));
}

struct conf { double b; double N_g; double sigma; double r_q; double kappa; double d;
              double (&RA)[9]; double (&RB)[9];
              double thetaA; double thetaB; double phiA; double phiB; };


// Density functions

// 0

inline double quark_diquark_density(double x, double y, void * conf)
{
    struct conf * c = (struct conf *)conf;
    return c->sigma * quark_diquark_thickness(x + c->b / 2., y, c->thetaA, c->phiA, c->r_q, c->d, c->N_g, c->kappa)
                    * quark_diquark_thickness(x - c->b / 2., y, c->thetaB, c->phiB, c->r_q, c->d, c->N_g, c->kappa);
}

// 1

inline double three_quarks_density(double x, double y, void * conf)
{
    struct conf * c = (struct conf *)conf;
    return c->sigma * three_quarks_thickness(x + c->b / 2., y, c->RA, c->r_q, c->N_g, c->kappa)
                    * three_quarks_thickness(x - c->b / 2., y, c->RB, c->r_q, c->N_g, c->kappa);
}

// 2

inline double mixed_density(double x, double y, void * conf)
{
    struct conf * c = (struct conf *)conf;
    return c->sigma * three_quarks_thickness(x + c->b / 2., y, c->RA, c->r_q, c->N_g, c->kappa)
                    * quark_diquark_thickness(x - c->b / 2., y, c->thetaB, c->phiB, c->r_q, c->d, c->N_g, c->kappa);
}

// 3

inline double three_quarks_one_gluon_body_density(double x, double y, void * conf)
{
    struct conf * c = (struct conf *)conf;
    return c->sigma * three_quarks_one_gluon_body_thickness(x + c->b / 2., y, c->RA, c->r_q, c->N_g, c->kappa)
                    * three_quarks_one_gluon_body_thickness(x - c->b / 2., y, c->RB, c->r_q, c->N_g, c->kappa);
}

// 6
inline double triangular_density(double x, double y, void * conf)
{
    return three_quarks_density(x, y, conf);
}

// 7
inline double qd_tr_density(double x, double y, void * conf)
{
    struct conf * c = (struct conf *)conf;
    return c->sigma * three_quarks_thickness(x + c->b / 2., y, c->RA, c->r_q, c->N_g, c->kappa)
                    * quark_diquark_thickness(x - c->b / 2., y, c->thetaB, c->phiB, c->r_q, c->d, c->N_g, c->kappa);
}



// Integrands

struct density { double (*density_function)(double, double, void *); double x_0; double y_0; conf c; };

double n_coll(double x[], size_t dim, void * density)
{
    struct density * d = (struct density *)density;
    return d->density_function(x[0] - d->x_0, x[1] - d->y_0, &d->c);
}

double mean_x(double x[], size_t dim, void * density)
{
    struct density * d = (struct density *)density;
    return x[0] * d->density_function(x[0] - d->x_0, x[1] - d->y_0, &d->c);
}

double mean_y(double x[], size_t dim, void * density)
{
    struct density * d = (struct density *)density;
    return x[1] * d->density_function(x[0] - d->x_0, x[1] - d->y_0, &d->c);
}

double mean_r2(double x[], size_t dim, void * density)
{
    struct density * d = (struct density *)density;
    return (x[0] * x[0] + x[1] * x[1]) * d->density_function(x[0] - d->x_0, x[1] - d->y_0, &d->c);
}

double mean_r2_cos2(double x[], size_t dim, void * density)
{
    struct density * d = (struct density *)density;
    return (x[0] * x[0] + x[1] * x[1]) * cos(2 * atan2(x[1], x[0])) * d->density_function(x[0] - d->x_0, x[1] - d->y_0, &d->c);
}

double mean_r2_sin2(double x[], size_t dim, void * density)
{
    struct density * d = (struct density *)density;
    return (x[0] * x[0] + x[1] * x[1]) * sin(2 * atan2(x[1], x[0])) * d->density_function(x[0] - d->x_0, x[1] - d->y_0, &d->c);
}

double mean_r3(double x[], size_t dim, void * density)
{
    struct density * d = (struct density *)density;
    return (x[0] * x[0] + x[1] * x[1]) * sqrt(x[0] * x[0] + x[1] * x[1]) * d->density_function(x[0] - d->x_0, x[1] - d->y_0, &d->c);
}

double mean_r3_cos3(double x[], size_t dim, void * density)
{
    struct density * d = (struct density *)density;
    return (x[0] * x[0] + x[1] * x[1]) * sqrt(x[0] * x[0] + x[1] * x[1]) * cos(3 * atan2(x[1], x[0])) * d->density_function(x[0] - d->x_0, x[1] - d->y_0, &d->c);
}

double mean_r3_sin3(double x[], size_t dim, void * density)
{
    struct density * d = (struct density *)density;
    return (x[0] * x[0] + x[1] * x[1]) * sqrt(x[0] * x[0] + x[1] * x[1]) * sin(3 * atan2(x[1], x[0])) * d->density_function(x[0] - d->x_0, x[1] - d->y_0, &d->c);
}

// Experimental parameters
const double exp_sigma_pp = 60 * 0.1;
const double exp_N = 30;
const double exp_dNdY = 5.8;
const double v2_over_eps = 0.3;
const double K_0 = 0.7;
const double R = 0.5;


// Functions used to rotate the proton density

inline double new_x(double x, double y, double theta, double phi, double psi)
{
    return x * (cos(theta) * cos(phi) * cos(psi) - sin(phi) * sin(psi)) + y * (-cos(psi) * sin(phi) - cos(theta) * cos(phi) * sin(psi));
}

inline double new_y(double x, double y, double theta, double phi, double psi)
{
    return x * (cos(theta) * cos(psi) * sin(phi) + cos(phi) * sin(psi)) + y * (cos(phi) * cos(psi) - cos(theta) * sin(phi) * sin(psi));
}

inline double new_z(double x, double y, double theta, double phi, double psi)
{
    return x * (-cos(psi) * sin(theta)) + y * (sin(theta) * sin(psi));
}

// Generate event

void generate_event(density * dens, gsl_rng * r)
{
    if(dens->density_function == &three_quarks_density || dens->density_function == &three_quarks_one_gluon_body_density)
    {
        for(int i = 0; i < 9; i++)
        {
            dens->c.RA[i] = gsl_ran_gaussian(r, sqrt(R * R - dens->c.r_q * dens->c.r_q));
            dens->c.RB[i] = gsl_ran_gaussian(r, sqrt(R * R - dens->c.r_q * dens->c.r_q));
        }

        // Shifting quarks positions so their center of mass is in the origin
        double mean_RA[3] = {0,0,0};
        double mean_RB[3] = {0,0,0};
        for(int i = 0; i < 3; i++)
        {
            mean_RA[i] = (dens->c.RA[i] + dens->c.RA[i + 3] + dens->c.RA[i + 6]) / 3.;
            mean_RB[i] = (dens->c.RB[i] + dens->c.RB[i + 3] + dens->c.RB[i + 6]) / 3.;
        }
        for(int i = 0; i < 9; i++)
        {
            dens->c.RA[i] = dens->c.RA[i] - mean_RA[i % 3];
            dens->c.RB[i] = dens->c.RB[i] - mean_RB[i % 3];
        }
    }
    if(dens->density_function == &quark_diquark_density)
    {
        dens->c.phiA = gsl_ran_flat(r, 0, 2 * M_PI);
        dens->c.phiB = gsl_ran_flat(r, 0, 2 * M_PI);
        double cos_thetaA = gsl_ran_flat(r, -1, 1);
        double cos_thetaB = gsl_ran_flat(r, -1, 1);
        dens->c.thetaA = acos(cos_thetaA);
        dens->c.thetaB = acos(cos_thetaB);
    }
    if(dens->density_function == &mixed_density)
    {
        for(int i = 0; i < 9; i++)
        {
            dens->c.RA[i] = gsl_ran_gaussian(r, sqrt(R * R - dens->c.r_q * dens->c.r_q));
        }

        // Shifting quarks positions so their center of mass is in the origin
        double mean_RA[3] = {0,0,0};
        for(int i = 0; i < 3; i++)
        {
            mean_RA[i] = (dens->c.RA[i] + dens->c.RA[i + 3] + dens->c.RA[i + 6]) / 3.;
        }
        for(int i = 0; i < 9; i++)
        {
            dens->c.RA[i] = dens->c.RA[i] - mean_RA[i % 3];
        }
        dens->c.phiB = gsl_ran_flat(r, 0, 2 * M_PI);
        double cos_thetaB = gsl_ran_flat(r, -1, 1);
        dens->c.thetaB = acos(cos_thetaB);
    }
	if(dens->density_function == &triangular_density)
    {
        double thetaA = acos(gsl_ran_flat(r, -1, 1));
		double phiA = gsl_ran_flat(r, 0, 2 * M_PI);
        double psiA = gsl_ran_flat(r, 0, 2 * M_PI);
		double thetaB = acos(gsl_ran_flat(r, -1, 1));
		double phiB = gsl_ran_flat(r, 0, 2 * M_PI);
        double psiB = gsl_ran_flat(r, 0, 2 * M_PI);

        dens->c.RA[0] = new_x(0, dens->c.d / 2., thetaA, phiA, psiA);
		dens->c.RA[1] = new_y(0, dens->c.d / 2., thetaA, phiA, psiA);
		dens->c.RA[2] = new_z(0, dens->c.d / 2., thetaA, phiA, psiA);

		dens->c.RA[3] = new_x(-sqrt(3) * dens->c.d / 4., -dens->c.d / 4., thetaA, phiA, psiA);
		dens->c.RA[4] = new_y(-sqrt(3) * dens->c.d / 4., -dens->c.d / 4., thetaA, phiA, psiA);
		dens->c.RA[5] = new_z(-sqrt(3) * dens->c.d / 4., -dens->c.d / 4., thetaA, phiA, psiA);

		dens->c.RA[6] = new_x(sqrt(3) * dens->c.d / 4., -dens->c.d / 4., thetaA, phiA, psiA);
		dens->c.RA[7] = new_y(sqrt(3) * dens->c.d / 4., -dens->c.d / 4., thetaA, phiA, psiA);
		dens->c.RA[8] = new_z(sqrt(3) * dens->c.d / 4., -dens->c.d / 4., thetaA, phiA, psiA);

        dens->c.RB[0] = new_x(0, dens->c.d / 2., thetaB, phiB, psiB);
        dens->c.RB[1] = new_y(0, dens->c.d / 2., thetaB, phiB, psiB);
        dens->c.RB[2] = new_z(0, dens->c.d / 2., thetaB, phiB, psiB);

        dens->c.RB[3] = new_x(-sqrt(3) * dens->c.d / 4., -dens->c.d / 4., thetaB, phiB, psiB);
        dens->c.RB[4] = new_y(-sqrt(3) * dens->c.d / 4., -dens->c.d / 4., thetaB, phiB, psiB);
        dens->c.RB[5] = new_z(-sqrt(3) * dens->c.d / 4., -dens->c.d / 4., thetaB, phiB, psiB);

        dens->c.RB[6] = new_x(sqrt(3) * dens->c.d / 4., -dens->c.d / 4., thetaB, phiB, psiB);
        dens->c.RB[7] = new_y(sqrt(3) * dens->c.d / 4., -dens->c.d / 4., thetaB, phiB, psiB);
        dens->c.RB[8] = new_z(sqrt(3) * dens->c.d / 4., -dens->c.d / 4., thetaB, phiB, psiB);
    }
	if(dens->density_function == &qd_tr_density)
    {
        double thetaA = acos(gsl_ran_flat(r, -1, 1));
		double phiA = gsl_ran_flat(r, 0, 2 * M_PI);
		double thetaB = acos(gsl_ran_flat(r, -1, 1));
		double phiB = gsl_ran_flat(r, 0, 2 * M_PI);
		
        dens->c.RA[0] = cos(thetaA) * cos(phiA) * 0 - cos(thetaA) * sin(phiA) * (dens->c.d / 2.);
		dens->c.RA[1] = sin(phiA) * 0 + cos(phiA) * (dens->c.d / 2.);
		dens->c.RA[2] = -sin(thetaA) * cos(phiA) * 0 + sin(thetaA) * sin(phiA) * (dens->c.d / 2.);

		dens->c.RA[3] = cos(thetaA) * cos(phiA) * (-sqrt(3) * dens->c.d / 4.) - cos(thetaA) * sin(phiA) * (-dens->c.d / 4.);
		dens->c.RA[4] = sin(phiA) * (-sqrt(3) * dens->c.d / 4.) + cos(phiA) * (-dens->c.d / 4.);
		dens->c.RA[5] = -sin(thetaA) * cos(phiA) * (-sqrt(3) * dens->c.d / 4.) + sin(thetaA) * sin(phiA) * (-dens->c.d / 4.);

		dens->c.RA[6] = cos(thetaA) * cos(phiA) * (sqrt(3) * dens->c.d / 4.) - cos(thetaA) * sin(phiA) * (-dens->c.d / 4.);
		dens->c.RA[7] = sin(phiA) * (sqrt(3) * dens->c.d / 4.) + cos(phiA) * (-dens->c.d / 4.);
		dens->c.RA[8] = -sin(thetaA) * cos(phiA) * (sqrt(3) * dens->c.d / 4.) + sin(thetaA) * sin(phiA) * (-dens->c.d / 4.);  

		dens->c.thetaB = thetaB;
		dens->c.phiB = phiB;

    }

}


void PP::performCollisions(int number, double Kappa)
{

    // Mixed types of collisions (number 5)
    double p = 0.2;  // probability of quark-diquark
    double random = 0;

    // Proton parameters
    double r_q = R / 2.;                                                                // !!!!!!!!!!!
    double d = 3 * R;                                                                   // !!!!!!!!!!!
    double kappa = Kappa;
    double N_g = 11.56;                                                                  // !!!!!!!!!!!
    double sigma = 4.3 * 0.1; //(in [mb] * 0.1)

    // Event parameters
    double b = 0;
    double RA[9] = {-0.5 + 0.5,0,0, 0.5 + 0.5,0,0, 0.5,sqrt(3)/2.,0}; // (x,y,z) of 1, 2 and 3 quark
    double RB[9] = {-0.5 + 0.5,0,0, 0.5 + 0.5,0,0, 0.5,sqrt(3)/2.,0}; // (x,y,z) of 1, 2 and 3 quark
    double thetaA = M_PI / 2.;
    double thetaB = M_PI / 2.;
    double phiA = M_PI / 2.;
    double phiB = M_PI / 2.;

    // Multiplicity coefficients
    double mmean_N_c = 6.9;
	double alpha = 1;
    double gamma = 1;
	alpha = exp_N / mmean_N_c;
    gamma = exp_dNdY / mmean_N_c;


    // Records creation - they keep track of parameter changes
    conf c = {b, N_g, sigma, r_q, kappa, d,
              RA, RB,
              thetaA, thetaB, phiA, phiB};
    density dens = {(number == 0 ? &quark_diquark_density : (
                    number == 1 ? &three_quarks_density : (
                    number == 2 ? &mixed_density : (
	            number == 3 ? &three_quarks_one_gluon_body_density : (
    		    number == 6 ? &triangular_density : &qd_tr_density ))))), 
	            0, 0, c};



    // Numerics parameters
    const double cutoff = 3 * R;                                                        // !!!!!!!!!!!
    double xl[2] = {-cutoff, -cutoff};
    double xu[2] = {cutoff, cutoff};
    size_t calls = 1000;                                                                // !!!!!!!!!!!


    // Monte Carlo integrands
    gsl_monte_function Mean_x = {&mean_x, 2, &dens};
    gsl_monte_function Mean_y = {&mean_y, 2, &dens};
    gsl_monte_function N_coll = {&n_coll, 2, &dens};
    gsl_monte_function Mean_r2 = {&mean_r2, 2, &dens};
    gsl_monte_function Mean_r2_cos2 = {&mean_r2_cos2, 2, &dens};
    gsl_monte_function Mean_r2_sin2 = {&mean_r2_sin2, 2, &dens};
    gsl_monte_function Mean_r3 = {&mean_r3, 2, &dens};
    gsl_monte_function Mean_r3_cos3 = {&mean_r3_cos3, 2, &dens};
    gsl_monte_function Mean_r3_sin3 = {&mean_r3_sin3, 2, &dens};

    // Random number generator
    const gsl_rng_type *T;
    gsl_rng *r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, time(0));

    // Monte Carlo integration
    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(2);
    double err;
    double err1;
    double N_c, x, y, r2, r3, r2_sin2, r2_cos2, r3_sin3, r3_cos3, eps_2, eps_3, S, K, v_2, v_3, ev_sigma;


    // Integration parameters
    const int n_samples_prerun = 5000;                                                   // !!!!!!!!!!!
    const int n_samples = 300000;                                                              // !!!!!!!!!!!

    // Averaged parameters
    double sigma_pp = 0;
    double mean_N_c = 0;
    double mean_eps_2 = 0;
    double mean_sq_eps_2 = 0;
    double mean_eps_3 = 0;
    double mean_sq_eps_3 = 0;
    double mean_r2 = 0;

    // Setting up the files

    ostringstream filename0;
    filename0 << "pp_results_" << number << "_" << kappa << ".dat";

    ofstream file0;
    file0.open(filename0.str().c_str());
    file0 << "PP MONTE CARLO" << endl;
    file0 << "samples prerun: " << n_samples_prerun << endl;
    file0 << "samples: " << n_samples << endl;
    file0 << "Model: " << number << endl;
    file0 << "r_q: " << r_q << endl;
    file0 << "d: " << d << endl;
    file0 << "kappa: " << kappa << endl;
    file0 << "sigma_gg: " << sigma << endl;

    ostringstream filename1;
    filename1 << "pp_data_" << number << "_" << kappa << ".dat";

    ofstream file1;
    file1.open(filename1.str().c_str());

    file1 << "b \t N_coll \t N \t sigma \t eps_2 \t eps_3 \t r^2" << endl;


    // RUN 1 - obtaining the correct sigma_pp

    double B = 0;
    const double B_max = 6 * R;
    double dB = B_max / (n_samples_prerun - 1.);
    //cout <<  dB << endl;

    double N_max = 5;
    double N_min = 15;

    while(true) // !!!!
    {
        sigma_pp = 0;

        dens.c.N_g = (N_max + N_min) / 2.;

        B = 0;
        for(int i = 0; i < n_samples_prerun; i++)
        {

            dens.c.b = B;
            //cout << "b: " << dens.c.b << endl;

            if(number == 5)
            {
                random = gsl_ran_flat(r, 0, 1);
                if(random < (1 - p) * (1 - p))
                {
                    dens.density_function = &three_quarks_density;
                }
                else
                {
                    if(random < 2 * (1 - p) * p + (1 - p) * (1 - p))
                        dens.density_function = &mixed_density;
                    else
                        dens.density_function = &quark_diquark_density;
                }
            }
			if(number == 8)
            {
                random = gsl_ran_flat(r, 0, 1);
                if(random < (1 - p) * (1 - p))
                {
                    dens.density_function = &triangular_density;
                }
                else
                {
                    if(random < 2 * (1 - p) * p + (1 - p) * (1 - p))
                        dens.density_function = &qd_tr_density;
                    else
                        dens.density_function = &quark_diquark_density;
                }
            }
            generate_event(&dens, r);

            // Number of binary collisions
            N_c = 0;
            gsl_monte_vegas_init(s);
            gsl_monte_vegas_integrate (&N_coll, xl, xu, 2, calls / 5, r, s, &N_c, &err);
            do
            {
                gsl_monte_vegas_integrate (&N_coll, xl, xu, 2, calls, r, s, &N_c, &err);
            }
            while(fabs(gsl_monte_vegas_chisq(s) - 1.0) > 0.5 || N_c > 10 * N_g * N_g);
            //cout << "N_coll: " << N_c << endl;


            // Event probability
            ev_sigma = event_sigma(N_c, dens.c.N_g, B);
            sigma_pp += ev_sigma;
            //cout << "Event sigma: " << ev_sigma << endl;

            // Increase B
            B = B + dB;
        }

        sigma_pp = dB * sigma_pp;
        cout << "Sigma_pp: " << sigma_pp << endl;
        cout << "N_g: " << dens.c.N_g << endl;

        if(abs(sigma_pp - exp_sigma_pp) < 0.05)
            break;
        else
        {
            if(sigma_pp - exp_sigma_pp < 0)
                N_max = dens.c.N_g;
            else
                N_min = dens.c.N_g;
        }
    }

    // RUN 2 - estimating mean number of collisions

    if(true) {

        sigma_pp = 0;
        mean_N_c = 0;

        B = 0;

        for (int i = 0; i < n_samples_prerun; i++) {

            dens.c.b = B;
            //cout << "b: " << dens.c.b << endl;

            if (number == 5) {
                random = gsl_ran_flat(r, 0, 1);
                if (random < (1 - p) * (1 - p)) {
                    dens.density_function = &three_quarks_density;
                }
                else {
                    if (random < 2 * (1 - p) * p + (1 - p) * (1 - p))
                        dens.density_function = &mixed_density;
                    else
                        dens.density_function = &quark_diquark_density;
                }
            }
            if (number == 8) {
                random = gsl_ran_flat(r, 0, 1);
                if (random < (1 - p) * (1 - p)) {
                    dens.density_function = &triangular_density;
                }
                else {
                    if (random < 2 * (1 - p) * p + (1 - p) * (1 - p))
                        dens.density_function = &qd_tr_density;
                    else
                        dens.density_function = &quark_diquark_density;
                }
            }
            generate_event(&dens, r);

            // Number of binary collisions
            N_c = 0;
            gsl_monte_vegas_init(s);
            gsl_monte_vegas_integrate(&N_coll, xl, xu, 2, calls / 5, r, s, &N_c, &err);
            do {
                gsl_monte_vegas_integrate(&N_coll, xl, xu, 2, calls, r, s, &N_c, &err);
            }
            while (fabs(gsl_monte_vegas_chisq(s) - 1.0) > 0.5);
            //cout << "N_coll: " << N_c << endl;

            // Event probability
            ev_sigma = event_sigma(N_c, dens.c.N_g, B);
            sigma_pp += ev_sigma;
            //cout << "Event sigma: " << ev_sigma << endl;

            // Calculate mean number of collisions
            mean_N_c += N_c * ev_sigma;

            // Increase B
            B = B + dB;

        }
        sigma_pp = sigma_pp * dB;
        cout << "Sigma_pp: " << sigma_pp << endl;
        mean_N_c = mean_N_c * dB / sigma_pp;
        cout << "Mean N_coll: " << mean_N_c << endl;

        // Determine alpha and gamma
        alpha = exp_N / mean_N_c;
        gamma = exp_dNdY / mean_N_c;
        cout << "alpha: " << alpha << endl;
        cout << "gamma: " << gamma << endl;
    }

    // RUN 3 - calculations of eccentricity and triangularity

    sigma_pp = 0;
    mean_N_c = 0;
    mean_eps_2 = 0;
    mean_sq_eps_2 = 0;
    mean_eps_3 = 0;
    mean_sq_eps_3 = 0;
    mean_r2 = 0;

    B = 0;
    dB = B_max / (n_samples - 1.);

    for(int i = 0; i < n_samples; i++)
    {

        dens.c.b = B;
        //cout << "b: " << dens.c.b << endl;

        file1 << B << "\t";

        dens.x_0 = 0;
        dens.y_0 = 0;

        if(number == 5)
        {
            random = gsl_ran_flat(r, 0, 1);
            if(random < (1 - p) * (1 - p))
            {
                dens.density_function = &three_quarks_density;
            }
            else
            {
                if(random < 2 * (1 - p) * p + (1 - p) * (1 - p))
                    dens.density_function = &mixed_density;
                else
                    dens.density_function = &quark_diquark_density;
            }
        }
		if(number == 8)
        {
            random = gsl_ran_flat(r, 0, 1);
            if(random < (1 - p) * (1 - p))
            {
                dens.density_function = &triangular_density;
            }
            else
            {
                if(random < 2 * (1 - p) * p + (1 - p) * (1 - p))
                    dens.density_function = &qd_tr_density;
                else
                    dens.density_function = &quark_diquark_density;
            }
        }
        generate_event(&dens, r);

        // Number of binary collisions
        N_c = 0;
        gsl_monte_vegas_init(s);
        gsl_monte_vegas_integrate (&N_coll, xl, xu, 2, calls / 5, r, s, &N_c, &err1);
        do
        {
            gsl_monte_vegas_integrate (&N_coll, xl, xu, 2, calls, r, s, &N_c, &err1);
        }
        while(fabs(gsl_monte_vegas_chisq(s) - 1.0) > 0.5 || N_c > 10 * N_g * N_g);
        //cout << "N_coll: " << N_c << endl;
        //cout << "N_coll_err: " << err1 << endl;
        //cout << "N: " << alpha * N_c << endl;

        file1 << N_c << "\t";
        file1 << N_c * alpha << "\t";


        // Event probability
        ev_sigma = event_sigma(N_c, dens.c.N_g, B);
        sigma_pp += ev_sigma;
        //cout << "Event sigma: " << ev_sigma << endl;

        file1 << ev_sigma << "\t";

        // Calculate mean number of collisions
        mean_N_c += N_c * ev_sigma;


        // Mean x
        x = 0;
        gsl_monte_vegas_init(s);
        gsl_monte_vegas_integrate (&Mean_x, xl, xu, 2, calls / 5, r, s, &x, &err);
        do
        {
            gsl_monte_vegas_integrate (&Mean_x, xl, xu, 2, calls, r, s, &x, &err);
        }
        while(fabs(gsl_monte_vegas_chisq(s) - 1.0) > 0.5);
        x = x / N_c;
        //cout << "Mean x: " << x << endl;


        // Mean y
        y = 0;
        gsl_monte_vegas_init(s);
        gsl_monte_vegas_integrate (&Mean_y, xl, xu, 2, calls / 5, r, s, &y, &err);
        do
        {
            gsl_monte_vegas_integrate (&Mean_y, xl, xu, 2, calls, r, s, &y, &err);
        }
        while(fabs(gsl_monte_vegas_chisq(s) - 1.0) > 0.5);
        y = y / N_c;
        // cout << "Mean y: " << y << endl;

        // Shifting density so its center is in the origin
        dens.x_0 = -x;
        dens.y_0 = -y;
        // cout << "Shifting..." << endl;

        // Mean x once more - should be 0
        x = 0;
        gsl_monte_vegas_init(s);
        gsl_monte_vegas_integrate (&Mean_x, xl, xu, 2, calls / 5, r, s, &x, &err);
        do
        {
            gsl_monte_vegas_integrate (&Mean_x, xl, xu, 2, calls, r, s, &x, &err);
        }
        while(fabs(gsl_monte_vegas_chisq(s) - 1.0) > 0.5);
        x = x / N_c;
        //cout << "Mean x: " << x << endl;


        // Eccentricity
        r2 = 0;
        gsl_monte_vegas_init(s);
        gsl_monte_vegas_integrate (&Mean_r2, xl, xu, 2, calls / 5, r, s, &r2, &err);
        do
        {
            gsl_monte_vegas_integrate (&Mean_r2, xl, xu, 2, calls, r, s, &r2, &err);
        }
        while(fabs(gsl_monte_vegas_chisq(s) - 1.0) > 0.5);
        r2_cos2 = 0;
        gsl_monte_vegas_init(s);
        gsl_monte_vegas_integrate (&Mean_r2_cos2, xl, xu, 2, calls / 5, r, s, &r2_cos2, &err);
        do
        {
            gsl_monte_vegas_integrate (&Mean_r2_cos2, xl, xu, 2, calls, r, s, &r2_cos2, &err);
        }
        while(fabs(gsl_monte_vegas_chisq(s) - 1.0) > 0.5);
        r2_sin2 = 0;
        gsl_monte_vegas_init(s);
        gsl_monte_vegas_integrate (&Mean_r2_sin2, xl, xu, 2, calls / 5, r, s, &r2_sin2, &err);
        do
        {
            gsl_monte_vegas_integrate (&Mean_r2_sin2, xl, xu, 2, calls, r, s, &r2_sin2, &err);
        }
        while(fabs(gsl_monte_vegas_chisq(s) - 1.0) > 0.5);
        eps_2 = sqrt(r2_cos2 * r2_cos2 + r2_sin2 * r2_sin2) / r2;
        //cout << "Eccentricity: " << eps_2 << endl;
        // cout << "Second harmonic plane angle: " << (atan2(r2_sin2, r2_cos2) + M_PI) / 2. / M_PI << " PI" << endl;

        file1 << eps_2 << "\t";

        // Triangularity
        r3 = 0;
        gsl_monte_vegas_init(s);
        gsl_monte_vegas_integrate (&Mean_r3, xl, xu, 2, calls / 5, r, s, &r3, &err);
        do
        {
            gsl_monte_vegas_integrate (&Mean_r3, xl, xu, 2, calls, r, s, &r3, &err);
        }
        while(fabs(gsl_monte_vegas_chisq(s) - 1.0) > 0.5);
        r3_cos3 = 0;
        gsl_monte_vegas_init(s);
        gsl_monte_vegas_integrate (&Mean_r3_cos3, xl, xu, 2, calls / 5, r, s, &r3_cos3, &err);
        do
        {
            gsl_monte_vegas_integrate (&Mean_r3_cos3, xl, xu, 2, calls, r, s, &r3_cos3, &err);
        }
        while(fabs(gsl_monte_vegas_chisq(s) - 1.0) > 0.5);
        r3_sin3 = 0;
        gsl_monte_vegas_init(s);
        gsl_monte_vegas_integrate (&Mean_r3_sin3, xl, xu, 2, calls / 5, r, s, &r3_sin3, &err);
        do
        {
            gsl_monte_vegas_integrate (&Mean_r3_sin3, xl, xu, 2, calls, r, s, &r3_sin3, &err);
        }
        while(fabs(gsl_monte_vegas_chisq(s) - 1.0) > 0.5);
        eps_3 = sqrt(r3_cos3 * r3_cos3 + r3_sin3 * r3_sin3) / r3;
        //cout << "Triangularity: " << eps_3 << endl;
        // cout << "Third harmonic plane angle: " << (atan2(r3_sin3, r3_cos3) + M_PI) / 3. /M_PI << " PI" << endl;

        file1 << eps_3 << "\t";
        file1 << r2 << endl;

        // Calculate mean eccentricity, triangularity and r^2
        mean_eps_2 += eps_2 * ev_sigma;
        mean_sq_eps_2 += eps_2 * eps_2 * ev_sigma;
        mean_eps_3 += eps_3 * ev_sigma;
        mean_sq_eps_3 += eps_3 * eps_3 * ev_sigma;
        mean_r2 += r2 * ev_sigma;

        // Estimate v_2 and v_3
        S = 2 * M_PI * r2 * sqrt(1 - eps_2);
        K = sqrt(3) * S / (sigma * gamma * N_c);
        v_2 = v2_over_eps * eps_2 / (1 + K / K_0);
        // cout << "Elliptic flow: " << v_2 << endl;
        v_3 = v2_over_eps * eps_3 / (1 + K / K_0);
        // cout << "Triangular flow: " << v_3 << endl;

        // Increase B
        B = B + dB;
    }

    sigma_pp *= dB;
    cout << "Sigma_pp: " << sigma_pp << endl;
    mean_N_c  = mean_N_c * dB / sigma_pp;
    cout << "Mean N_coll: " << mean_N_c << endl;

    mean_eps_2  = mean_eps_2 * dB / sigma_pp;
    cout << "Mean eccentricity: " << mean_eps_2 << endl;
    mean_sq_eps_2 = mean_sq_eps_2 * dB / sigma_pp;
    cout << "RMS eccentricity: " << sqrt(mean_sq_eps_2) << endl;
    mean_eps_3  = mean_eps_3 * dB / sigma_pp;
    cout << "Mean triangularity: " << mean_eps_3 << endl;
    mean_sq_eps_3 = mean_sq_eps_3 * dB / sigma_pp;
    cout << "RMS triangularity: " << sqrt(mean_sq_eps_3) << endl;
    mean_r2  = mean_r2 * dB / sigma_pp;
    cout << "Square root of mean r^2: " << sqrt(mean_r2) << endl;

    gsl_monte_vegas_free(s);

    file1.close();

    file0 << "N_g: " << dens.c.N_g << endl;
    file0 << "Sigma_pp: " << sigma_pp << endl;
    file0 << "Mean N_coll: " << mean_N_c << endl;
    file0 << "Alpha: " << alpha << endl;
    file0 << "Gamma: " << gamma << endl;
    file0 << "Mean eccentricity: " << mean_eps_2 << endl;
    file0 << "RMS eccentricity: " << sqrt(mean_sq_eps_2) << endl;
    file0 << "Mean triangularity: " << mean_eps_3 << endl;
    file0 << "RMS triangularity: " << sqrt(mean_sq_eps_3) << endl;
    file0 << "Square root of mean r^2: " << sqrt(mean_r2) << endl;

    file0.close();
}
