#ifndef UTILS_H_NR803MPH
#define UTILS_H_NR803MPH

#include <time.h>
#include <sys/time.h>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>


double get_wall_time(){
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

inline double abs2(const std::complex<double> &z) {
  return z.real() * z.real() + z.imag() * z.imag();
}


void eigenvalues(Material& mat, double& eigmin, double& eigmax) {


    typedef Eigen::Matrix<double,6,6> Elasticity;

    Elasticity C;
    int k=0;
    for(int i=0; i<6; ++i) {
        for(int j=0; j<6; ++j) {
            C(i,j) = mat.props[k++];
        }
    }

    /* std::cout << C << std::endl; */

    Eigen::SelfAdjointEigenSolver<Elasticity> es;
    es.compute(C);

    Eigen::Matrix<double,6,1> ev = es.eigenvalues();

    /* std::cout << ev << std::endl; */
    eigmin = ev[0];
    eigmax = ev[5];


}


Eigen::Matrix<double,6,6> mat_to_C(Material& mat) {
    Eigen::Matrix<double,6,6> C;
    int k=0;
    for(int i=0; i<6; ++i) {
        for(int j=0; j<6; ++j) {
            C(i,j) = mat.props[k++];
        }
    }
    return C;
}

Eigen::Matrix<double,6,1> mat_to_a(Material& mat) {
    Eigen::Matrix<double,6,1> a;
    for(int i=0; i<6; ++i) {
            a[i] = mat.props[i+36];
    }
    return a;
}


Eigen::Matrix<double,6,6>
rotate_matrix(Eigen::Matrix<double,6,6> C,
              double alpha, double beta, double gamma) {

    typedef Eigen::Matrix<double,6,6> Mat6;
    typedef Eigen::Matrix<double,3,3> Mat3;
    typedef Eigen::Matrix<double,3,1> Vec3;

    Mat3 Rx, Ry, Rz;
    Rz << cos(alpha), -sin(alpha), 0.0, sin(alpha), cos(alpha), 0.0, 0.0, 0.0, 1.0;
    Ry  << cos(beta), 0, sin(beta), 0, 1, 0, -sin(beta), 0, cos(beta);
    Rx << 1,0,0, 0,cos(gamma), -sin(gamma), 0,sin(gamma),  cos(gamma);

    Mat3 Q;
    Q = Rz*Ry*Rx;

    double R11,R12,R13;
    double R21,R22,R23;
    double R31,R32,R33;

    R11 = Q(0,0);
    R12 = Q(0,1);
    R13 = Q(0,2);
    R21 = Q(1,0);
    R22 = Q(1,1);
    R23 = Q(1,2);
    R31 = Q(2,0);
    R32 = Q(2,1);
    R33 = Q(2,2);

    const double sqrt2 = sqrt(2.0);

    Mat6 R;
    R  << R11*R11, R12*R12, R13*R13, sqrt2*R11*R12, sqrt2*R12*R13, sqrt2*R11*R13,
    R21*R21, R22*R22, R23*R23, sqrt2*R21*R22, sqrt2*R22*R23, sqrt2*R21*R23,
    R31*R31, R32*R32, R33*R33, sqrt2*R31*R32, sqrt2*R32*R33, sqrt2*R31*R33,
    sqrt2*R11*R21, sqrt2*R12*R22, sqrt2*R13*R23, R11*R22+R12*R21, R12*R23+R13*R22, R11*R23+R13*R21,
    sqrt2*R21*R31, sqrt2*R22*R32, sqrt2*R23*R33, R21*R32+R22*R31, R22*R33+R23*R32, R21*R33+R23*R31,
    sqrt2*R11*R31, sqrt2*R12*R32, sqrt2*R13*R33, R11*R32+R12*R31, R12*R33+R13*R32, R11*R33+R13*R31;

    // Add sqrt2
    for(int i = 3; i < 6; ++i) {
        for (int j=3; j < 6; ++j) {
            C(i,j) *= 2.0;
        }
    }

    for (int i = 0; i< 3; ++i) {
        for(int j = 3; j < 6; ++j) {
            C(i,j) *= sqrt2;
            C(j,i) *= sqrt2;
        }
    }

    Mat6 Crot;
    Crot = R.transpose() * C * R;

    for(int i = 3; i < 6; ++i) {
        for (int j=3; j < 6; ++j) {
            Crot(i,j) /= 2.0;
        }
    }

    for (int i = 0; i< 3; ++i) {
        for(int j = 3; j < 6; ++j) {
            Crot(i,j) /= sqrt2;
            Crot(j,i) /= sqrt2;
        }
    }

    return Crot;
}

Eigen::Matrix<double,6,1>
rotate_matrix(Eigen::Matrix<double,6,1> a,
              double alpha, double beta, double gamma) {

    typedef Eigen::Matrix<double,6,6> Mat6; 
    typedef Eigen::Matrix<double,3,3> Mat3;
    typedef Eigen::Matrix<double,6,1> Vec6; 

    Mat3 Rx, Ry, Rz;
    Rz << cos(alpha), -sin(alpha), 0.0, sin(alpha), cos(alpha), 0.0, 0.0, 0.0, 1.0;
    Ry  << cos(beta), 0, sin(beta), 0, 1, 0, -sin(beta), 0, cos(beta);
    Rx << 1,0,0, 0,cos(gamma), -sin(gamma), 0,sin(gamma),  cos(gamma);

    Mat3 Q;
    Q = Rz*Ry*Rx;

    double R11,R12,R13;
    double R21,R22,R23;
    double R31,R32,R33;

    R11 = Q(0,0);
    R12 = Q(0,1);
    R13 = Q(0,2);
    R21 = Q(1,0);
    R22 = Q(1,1);
    R23 = Q(1,2);
    R31 = Q(2,0);
    R32 = Q(2,1);
    R33 = Q(2,2);

    const double sqrt2 = sqrt(2.0);

    Mat6 R;
    R  << R11*R11, R12*R12, R13*R13, sqrt2*R11*R12, sqrt2*R12*R13, sqrt2*R11*R13,
    R21*R21, R22*R22, R23*R23, sqrt2*R21*R22, sqrt2*R22*R23, sqrt2*R21*R23,
    R31*R31, R32*R32, R33*R33, sqrt2*R31*R32, sqrt2*R32*R33, sqrt2*R31*R33,
    sqrt2*R11*R21, sqrt2*R12*R22, sqrt2*R13*R23, R11*R22+R12*R21, R12*R23+R13*R22, R11*R23+R13*R21,
    sqrt2*R21*R31, sqrt2*R22*R32, sqrt2*R23*R33, R21*R32+R22*R31, R22*R33+R23*R32, R21*R33+R23*R31,
    sqrt2*R11*R31, sqrt2*R12*R32, sqrt2*R13*R33, R11*R32+R12*R31, R12*R33+R13*R32, R11*R33+R13*R31;

    for(int i = 3; i < 6; ++i) {
            a[i] /= sqrt2;
    }

    Vec6 arot;
    arot = R.transpose() * a;
    for(int i = 3; i < 6; ++i) {
        a[i] *= sqrt2;
    }
    
    return arot;
}


#endif /* end of include guard: UTILS_H_NR803MPH */

