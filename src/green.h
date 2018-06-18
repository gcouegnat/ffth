#ifndef GREEN_H_GS3VYBQ0
#define GREEN_H_GS3VYBQ0

/* #include <Eigen/Core> */

void compute_green(const double q0, const double q1, const double q2,
                   const double lmbda, const double mu,
                   double &G00, double &G01, double &G02, double &G03, double &G04,
                   double &G05, double &G11, double &G12, double &G13,
                   double &G14, double &G15, double &G22, double &G23,
                   double &G24, double &G25, double &G33, double &G34,
                   double &G35, double &G44, double &G45, double &G55)

// Eigen::Matrix<double, 6, 6> &G)

{

    /* const    double k0 = 2.0*M_PI*q0; */
    /* const    double k1 = 2.0*M_PI*q1; */
    /* const    double k2 = 2.0*M_PI*q2; */

    const double k0 = sin(M_PI * q0) * cos(M_PI*q1) * cos(M_PI*q2); 
    const double k1 = cos(M_PI * q0) * sin(M_PI*q1) * cos(M_PI*q2); 
    const double k2 = cos(M_PI * q0) * cos(M_PI*q1) * sin(M_PI*q2); 

    const double k0k0 = k0 * k0;
    const double k1k1 = k1 * k1;
    const double k2k2 = k2 * k2;

    const double k0k0k0 = k0 * k0k0;
    const double k1k1k1 = k1 * k1k1;
    const double k2k2k2 = k2 * k2k2;

    const double k0k0k0k0 = k0k0 * k0k0;
    const double k1k1k1k1 = k1k1 * k1k1;
    const double k2k2k2k2 = k2k2 * k2k2;

    const double k = (k0k0 + k1k1 + k2k2);
    const double ksq = k * k;

    const double nu = 0.5 * lmbda / (lmbda + mu);
    const double tmp1 = 1.0 / (mu * (-nu + 1) * ksq);
    const double tmp2 = 1.0 / (mu * k);

    /* g00 = -0.5 * k0k0k0k0 * tmp1 + k0k0 * tmp2; */
    /* g01 = -0.5 * k0k0 * k1k1 * tmp1; */
    /* g02 = -0.5 * k0k0 * k2k2 * tmp1; */
    /* g03 = -k0k0k0 * k1 * tmp1 + k0 * k1 * tmp2; */
    /* g04 = -k0k0 * k1 * k2 * tmp1; */
    /* g05 = -k0k0k0 * k2 * tmp1 + k0 * k2 * tmp2; */

    /* g11 = -0.5 * k1k1k1k1 * tmp1 + k1k1 * tmp2; */
    /* g12 = -0.5 * k1k1 * k2k2 * tmp1; */
    /* g13 = -k0 * k1k1k1 * tmp1 + k0 * k1 * tmp2; */
    /* g14 = -k1k1k1 * k2 * tmp1 + k1 * k2 * tmp2; */
    /* g15 = -k0 * k1k1 * k2 * tmp1; */

    /* g22 = -0.5 * k2k2k2k2 * tmp1 + k2k2 * tmp2; */
    /* g23 = -k0 * k1 * k2k2 * tmp1; */
    /* g24 = -k1 * k2k2k2 * tmp1 + 0.5 * k1 * k2 * tmp2; */
    /* g25 = -k0 * k2k2k2 * tmp1 + 0.5 * k0 * k2 * tmp2; */

    /* g33 = -2.0 * k0k0 * k1k1 * tmp1 + (k0k0 + k1k1) * tmp2; */
    /* g34 = -2.0 * k0 * k1k1 * k2 * tmp1 + k0 * k2 * tmp2; */
    /* g35 = -2.0 * k0k0 * k1 * k2 * tmp1 + k1 * k2 * tmp2; */

    /* g44 = -2.0 * k1k1 * k2k2 * tmp1 + (k1k1 + k2k2) * tmp2; */
    /* g45 = -2.0 * k0 * k1 * k2k2 * tmp1 + k0 * k1 * tmp2; */

    /* g55 = -2.0 * k0k0 * k2k2 * tmp1 + (k0k0 + k2k2) * tmp2; */

    const double g00 = -0.5 * k0k0k0k0 * tmp1 + k0k0 * tmp2;
    const double g01 = -0.5 * k0k0 * k1k1 * tmp1;
    const double g02 = -0.5 * k0k0 * k2k2 * tmp1;
    const double g03 = -0.5 * k0k0k0 * k1 * tmp1 + 0.5 * k0 * k1 * tmp2;
    const double g04 = -0.5 * k0k0 * k1 * k2 * tmp1;
    const double g05 = -0.5 * k0k0k0 * k2 * tmp1 + 0.5 * k0 * k2 * tmp2;
    const double g11 = -0.5 * k1k1k1k1 * tmp1 + k1k1 * tmp2;
    const double g12 = -0.5 * k1k1 * k2k2 * tmp1;
    const double g13 = -0.5 * k0 * k1k1k1 * tmp1 + 0.5 * k0 * k1 * tmp2;
    const double g14 = -0.5 * k1k1k1 * k2 * tmp1 + 0.5 * k1 * k2 * tmp2;
    const double g15 = -0.5 * k0 * k1k1 * k2 * tmp1;
    const double g22 = -0.5 * k2k2k2k2 * tmp1 + k2k2 * tmp2;
    const double g23 = -0.5 * k0 * k1 * k2k2 * tmp1;
    const double g24 = -0.5 * k1 * k2k2k2 * tmp1 + 0.5 * k1 * k2 * tmp2;
    const double g25 = -0.5 * k0 * k2k2k2 * tmp1 + 0.5 * k0 * k2 * tmp2;
    const double g33 = -0.5 * k0k0 * k1k1 * tmp1 + 0.25 * (k0k0 + k1k1) * tmp2;
    const double g34 = -0.5 * k0 * k1k1 * k2 * tmp1 + 0.25 * k0 * k2 * tmp2;
    const double g35 = -0.5 * k0k0 * k1 * k2 * tmp1 + 0.25 * k1 * k2 * tmp2;
    const double g44 = -0.5 * k1k1 * k2k2 * tmp1 + 0.25 * (k1k1 + k2k2) * tmp2;
    const double g45 = -0.5 * k0 * k1 * k2k2 * tmp1 + 0.25 * k0 * k1 * tmp2;
    const double g55 = -0.5 * k0k0 * k2k2 * tmp1 + 0.25 * (k0k0 + k2k2) * tmp2;


    /* const double g00 = -0.5 * k0k0k0k0 * tmp1 + k0k0 * tmp2; */
    /* const double g01 = -0.5 * k0k0 * k1k1 * tmp1; */
    /* const double g02 = -0.5 * k0k0 * k2k2 * tmp1; */
    /* const double g03 = -0.5 * k0k0k0 * k1 * tmp1 + 0.5 * k0 * k1 * tmp2; */
    /* const double g04 = -0.5 * k0k0 * k1 * k2 * tmp1; */
    /* const double g05 = -0.5 * k0k0k0 * k2 * tmp1 + 0.5 * k0 * k2 * tmp2; */
    /* const double g11 = -0.5 * k1k1k1k1 * tmp1 + k1k1 * tmp2; */
    /* const double g12 = -0.5 * k1k1 * k2k2 * tmp1; */
    /* const double g13 = -0.5 * k0 * k1k1k1 * tmp1 + 0.5 * k0 * k1 * tmp2; */
    /* const double g14 = -0.5 * k1k1k1 * k2 * tmp1 + 0.5 * k1 * k2 * tmp2; */
    /* const double g15 = -0.5 * k0 * k1k1 * k2 * tmp1; */
    /* const double g22 = -0.5 * k2k2k2k2 * tmp1 + k2k2 * tmp2; */
    /* const double g23 = -0.5 * k0 * k1 * k2k2 * tmp1; */
    /* const double g24 = -0.5 * k1 * k2k2k2 * tmp1 + 0.5 * k1 * k2 * tmp2; */
    /* const double g25 = -0.5 * k0 * k2k2k2 * tmp1 + 0.5 * k0 * k2 * tmp2; */
    /* const double g33 = -0.5 * k0k0 * k1k1 * tmp1 + 0.25 * (k0k0 + k1k1) * tmp2;
     */
    /* const double g34 = -0.5 * k0 * k1k1 * k2 * tmp1 + 0.25 * k0 * k2 * tmp2; */
    /* const double g35 = -0.5 * k0k0 * k1 * k2 * tmp1 + 0.25 * k1 * k2 * tmp2; */
    /* const double g44 = -0.5 * k1k1 * k2k2 * tmp1 + 0.25 * (k1k1 + k2k2) * tmp2;
     */
    /* const double g45 = -0.5 * k0 * k1 * k2k2 * tmp1 + 0.25 * k0 * k1 * tmp2; */
    /* const double g55 = -0.5 * k0k0 * k2k2 * tmp1 + 0.25 * (k0k0 + k2k2) * tmp2;
     */

    /* G(0, 0) = g00; */
    /* G(1, 1) = g11; */
    /* G(2, 2) = g22; */

    /* G(0, 1) = G(1, 0) = g01; */
    /* G(0, 2) = G(2, 0) = g02; */
    /* G(1, 2) = G(2, 1) = g12; */

    /* G(0, 3) = G(3, 0) = 2 * g03; */
    /* G(0, 4) = G(4, 0) = 2 * g04; */
    /* G(0, 5) = G(5, 0) = 2 * g05; */
    /* G(1, 3) = G(3, 1) = 2 * g13; */
    /* G(1, 4) = G(4, 1) = 2 * g14; */
    /* G(1, 5) = G(5, 1) = 2 * g15; */
    /* G(2, 3) = G(3, 2) = 2 * g23; */
    /* G(2, 4) = G(4, 2) = 2 * g24; */
    /* G(2, 5) = G(5, 2) = 2 * g25; */

    /* G(3, 3) = 4 * g33; */
    /* G(4, 4) = 4 * g44; */
    /* G(5, 5) = 4 * g55; */
    /* G(3, 4) = G(4, 3) = 4 * g34; */
    /* G(3, 5) = G(5, 3) = 4 * g35; */
    /* G(4, 5) = G(5, 4) = 4 * g45; */



    G00 = g00;
    G01 = g01;
    G02 = g02;
    G03 = 2*g03;
    G04 = 2*g04;
    G05 = 2*g05;

    G11 = g11;
    G12 = g12;
    G13 = 2*g13;
    G14 = 2*g14;
    G15 = 2*g15;
    
    G22 = g22;
    G23 = 2*g23;
    G24 = 2*g24;
    G25 = 2*g25;
    
    G33 = 4*g33;
    G34 = 4*g34;
    G35 = 4*g35;
    
    G44 = 4*g44;
    G45 = 4*g45;
    
    G55 = 4*g55;

}

#endif /* end of include guard: GREEN_H_GS3VYBQ0 */
