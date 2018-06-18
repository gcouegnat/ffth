#include "fft.h"
#include "green.h"
#include "image.h"
#include "material.h"
#include "utils.h"
#include "argparser.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include <omp.h>


static double average(fft_real_view& f);
static void init_freq(std::vector<double> &ksi1, std::vector<double> &ksi2, std::vector<double> &ksi3);
static void write_binary_data( std::ofstream& output, fft_real_view& data);
static void write_raw( std::ofstream& output, fft_real_view& data);


int main(int argc, const char *argv[]) {

    typedef unsigned int uint;

    Argparser parser;

    /* parser.add_option<uint>("size", ""); */
    parser.add_option<uint>("nx", "", 1);
    parser.add_option<uint>("ny", "", 1);
    parser.add_option<uint>("nz", "", 1);
    parser.add_option<std::string>("input", "image");
    parser.add_option<std::string>("output", "output", "output.vti");
    parser.add_option<std::string>("phase", "phase");
    parser.add_option<std::string>("material", "material");
    parser.add_option<std::string>("loading", "loading");
    parser.add_option<double>("temperature", "temperature", 0);
    parser.add_option<double>("rel", "", 1e-6);
    parser.add_option<double>("div", "", 1e-4);

    parser.parse(argc, argv);
    parser.print_values();


    /* const uint N = 128; */
    /* const uint n = parser.get<int>("size"); */
    /* std::cout << n << std::endl; */
    const uint nx = parser.get<int>("nx");
    const uint ny = parser.get<int>("ny");
    const uint nz = parser.get<int>("nz");

    // Image
    const std::string imfilename = parser.get<std::string>("input");
    std::cout << "Loading image " << imfilename << std::endl;
    Image<unsigned char> IM(nx, ny, nz);
    IM.load(imfilename);

    // Materials
    const std::string matfilename = parser.get<std::string>("material");
    std::cout << "Loading materials from " << matfilename << std::endl;
    std::ifstream fp(matfilename);

    unsigned int nmat;
    fp >> nmat;
    std::cout << "Reading " << nmat << " materials" << std::endl;
    
    std::vector<Material> mat(nmat);

    for (int k = 0; k < nmat; ++k) {
        for (int i = 0; i < 42; ++i) {
            fp >> mat[k].props[i];
        }
    }
    fp.close();


    double gammamin, gammamax;
    eigenvalues(mat[0], gammamin, gammamax);

    double eigmin, eigmax;
    eigenvalues(mat[1], eigmin, eigmax);
    gammamin = std::min(gammamin, eigmin);
    gammamax = std::max(gammamax, eigmax);

    double gamma = 0.5*(gammamin+gammamax);

    // Reference material
    /* const double lmbda0 = 213327.464789; */
    /* const double mu0 = 40633.8028169; */
    /* const double lmbda0 = 152878.521127; */
    /* const double mu0 = 29119.7183099; */
    const double lmbda0 = 0.0;
    const double mu0 = 0.5*gamma;

    /* std::cout << "lambda0 " << lmbda0 << '\n'; */
    /* std::cout << "mu0     " << mu0    << '\n'; */

    // Phases
    std::cout << "Loading phases" << std::endl;
    int nphases;
    const std::string phfilename = parser.get<std::string>("phase");
    std::ifstream fp2(phfilename);

    fp2 >> nphases;
    std::cout << nphases << std::endl;

    std::vector<Phase> phases(nphases);
    for(int i=0; i < nphases; ++i) {
        fp2 >> phases[i].matid;
        for(int k=0; k < 3; ++k) {
            fp2 >> phases[i].angles[k];
        }
        /* std::cout << i << ", " << phases[i].matid << std::endl; */
    }

    // Loading
    //

    double EMACRO[6];
    const std::string loadfilename = parser.get<std::string>("loading");
    std::ifstream fp3(loadfilename);
    for(int i=0; i < 6; ++i) {
        fp3 >> EMACRO[i];
    }
    fp3.close();

    const double T = parser.get<double>("temperature");

    // Tolerance
    const double reltol = parser.get<double>("rel");
    const double divtol = parser.get<double>("div");
    const int nitermax = 2500;


    std::cout << "Creating FFT buffers" << std::endl;
    const int nthreads = omp_get_max_threads();

    std::cout << "Using " << nthreads << " threads" << std::endl;

    fftw_init_threads();
    fftw_plan_with_nthreads(nthreads);


    fft E11(nx, ny, nz), E22(nx, ny, nz), E33(nx, ny, nz), E12(nx, ny, nz),
        E23(nx, ny, nz), E31(nx, ny, nz);
    fft_real_view eto11(E11), eto22(E22), eto33(E33), eto12(E12), eto23(E23),
                  eto31(E31);
    fft_complex_view etobar11(E11), etobar22(E22), etobar33(E33), etobar12(E12),
                     etobar23(E23), etobar31(E31);

    fft S11(nx, ny, nz), S22(nx, ny, nz), S33(nx, ny, nz), S12(nx, ny, nz),
        S23(nx, ny, nz), S31(nx, ny, nz);
    fft_real_view sig11(S11), sig22(S22), sig33(S33), sig12(S12), sig23(S23),
                  sig31(S31);
    fft_complex_view sigbar11(S11), sigbar22(S22), sigbar33(S33), sigbar12(S12),
                     sigbar23(S23), sigbar31(S31);

    std::cout << "Init FFT frequencies" << std::endl;
    std::vector<double> ksi1(nx);
    std::vector<double> ksi2(ny);
    std::vector<double> ksi3(nz / 2 + 1);
    init_freq(ksi1, ksi2, ksi3);

    for (uint i = 0; i < nx; ++i) {
        for (uint j = 0; j < ny; ++j) {
            for (uint k = 0; k < nz; ++k) {
                eto11(i, j, k) = EMACRO[0];
                eto22(i, j, k) = EMACRO[1];
                eto33(i, j, k) = EMACRO[2];
                eto12(i, j, k) = EMACRO[3];
                eto23(i, j, k) = EMACRO[4];
                eto31(i, j, k) = EMACRO[5];
            }
        }
    }

    double tinit = get_wall_time();

    double fft_cum_time = 0.0;
    double green_cum_time = 0.0;
    double behavior_cum_time = 0.0;
    double conv_cum_time = 0.0;

    double norm_old = 1.0;
    double rel=1.0;
    double diverr=1.0;

    bool converged = false;

    for (int niter = 0; niter < nitermax; ++niter) {


        // --- Apply behavior -------------------------------------------------------
        double t0 = get_wall_time();

        /* #pragma omp parallel for */
        for (uint i = 0; i < nx; ++i) {
            for (uint j = 0; j < ny; ++j) {
                for (uint k = 0; k < nz; ++k) {

                    const long iph = IM(i, j, k);
                    
                    const int matid = phases[iph].matid;
                    const double alpha = phases[iph].angles[0];
                    const double beta  = phases[iph].angles[1];
                    const double gamma = phases[iph].angles[2];

                    /* if (iph > 0) */
                    /*     std::cout << matid << alpha << beta << gamma << std::endl; */

                    /* const Material &m = mat[matid]; */
                    Eigen::Matrix<double,6,6> C    = mat_to_C(mat[matid]);
                    Eigen::Matrix<double,6,6> Crot;

                    Eigen::Matrix<double,6,1> a = mat_to_a(mat[matid]);
                    Eigen::Matrix<double,6,1> arot;

                    if ((alpha*alpha + beta*beta + gamma*gamma) < 1e-6) {
                        Crot = C;
                        arot = a;
                    } else {
                        Crot = rotate_matrix(C, alpha, beta, gamma);
                        arot = rotate_matrix(a, alpha, beta, gamma);
                    }

                    Eigen::Matrix<double,6,1> e, s;
                    e <<  eto11(i, j, k), eto22(i, j, k), eto33(i, j, k), eto12(i, j, k), eto23(i, j, k), eto31(i, j, k);

                    s = Crot * (e - T*arot);

                    sig11(i,j,k) = s[0];
                    sig22(i,j,k) = s[1];
                    sig33(i,j,k) = s[2];
                    sig12(i,j,k) = s[3];
                    sig23(i,j,k) = s[4];
                    sig31(i,j,k) = s[5];

                    /* const double e1 = eto11(i, j, k); */
                    /* const double e2 = eto22(i, j, k); */
                    /* const double e3 = eto33(i, j, k); */
                    /* const double e4 = eto12(i, j, k); */
                    /* const double e5 = eto23(i, j, k); */
                    /* const double e6 = eto31(i, j, k); */

                    /* sig11(i, j, k) = Crot(0,0) * (e1 - arot[0]*T) + Crot(0,1) * (e2 - arot[1]*T)  +  Crot(0,2) * (e3 - arot[2]*T)  +  Crot(0,3) * (e4  - arot[3]*T) +  Crot(0,4) * (e5  - arot[4]*T) +  Crot(0,5) * (e6 - arot[5]*T) ; */
                    /* sig22(i, j, k) = Crot(1,0) * (e1 - arot[0]*T) + Crot(1,1) * (e2 - arot[1]*T)  +  Crot(1,2) * (e3 - arot[2]*T)  +  Crot(1,3) * (e4  - arot[3]*T) +  Crot(1,4) * (e5  - arot[4]*T) +  Crot(1,5) * (e6 - arot[5]*T) ; */
                    /* sig33(i, j, k) = Crot(2,0) * (e1 - arot[0]*T) + Crot(2,1) * (e2 - arot[1]*T)  +  Crot(2,2) * (e3 - arot[2]*T)  +  Crot(2,3) * (e4  - arot[3]*T) +  Crot(2,4) * (e5  - arot[4]*T) +  Crot(2,5) * (e6 - arot[5]*T) ; */
                    /* sig12(i, j, k) = Crot(3,0) * (e1 - arot[0]*T) + Crot(3,1) * (e2 - arot[1]*T)  +  Crot(3,2) * (e3 - arot[2]*T)  +  Crot(3,3) * (e4  - arot[3]*T) +  Crot(3,4) * (e5  - arot[4]*T) +  Crot(3,5) * (e6 - arot[5]*T) ; */
                    /* sig23(i, j, k) = Crot(4,0) * (e1 - arot[0]*T) + Crot(4,1) * (e2 - arot[1]*T)  +  Crot(4,2) * (e3 - arot[2]*T)  +  Crot(4,3) * (e4  - arot[3]*T) +  Crot(4,4) * (e5  - arot[4]*T) +  Crot(4,5) * (e6 - arot[5]*T) ; */
                    /* sig31(i, j, k) = Crot(5,0) * (e1 - arot[0]*T) + Crot(5,1) * (e2 - arot[1]*T)  +  Crot(5,2) * (e3 - arot[2]*T)  +  Crot(5,3) * (e4  - arot[3]*T) +  Crot(5,4) * (e5  - arot[4]*T) +  Crot(5,5) * (e6 - arot[5]*T) ; */

                    /* sig11(i, j, k) = m.props[0] * e1 + m.props[1] * e2 + m.props[2] * e3 + m.props[3] * e4 + m.props[4] * e5 + m.props[5] * e6; */
                    /* sig22(i, j, k) = m.props[6] * e1 + m.props[7] * e2 + m.props[8] * e3 + m.props[9] * e4 + m.props[10] * e5 + m.props[11] * e6; */
                    /* sig33(i, j, k) = m.props[12] * e1 + m.props[13] * e2 + m.props[14] * e3 + m.props[15] * e4 + m.props[16] * e5 + m.props[17] * e6; */
                    /* sig12(i, j, k) = m.props[18] * e1 + m.props[19] * e2 + m.props[20] * e3 + m.props[21] * e4 + m.props[22] * e5 + m.props[23] * e6; */
                    /* sig23(i, j, k) = m.props[24] * e1 + m.props[25] * e2 + m.props[26] * e3 + m.props[27] * e4 + m.props[28] * e5 + m.props[29] * e6; */
                    /* sig31(i, j, k) = m.props[30] * e1 + m.props[31] * e2 + m.props[32] * e3 + m.props[33] * e4 + m.props[34] * e5 + m.props[35] * e6; */
                }
            }
        }
        behavior_cum_time += get_wall_time() - t0;

        // --- FFT forward -------------------------------------------------------
        //
        /* std::cout << "FFT Forward" << std::endl; */
        double t1 = get_wall_time();
        E11.forward_fft();
        E22.forward_fft();
        E33.forward_fft();
        E12.forward_fft();
        E23.forward_fft();
        E31.forward_fft();
        S11.forward_fft();
        S22.forward_fft();
        S33.forward_fft();
        S12.forward_fft();
        S23.forward_fft();
        S31.forward_fft();
        fft_cum_time += get_wall_time() - t1;


        // --- Check convergence -------------------------------------------------------

        double t4 = get_wall_time();
        // norm of stress
        double norm_sig = 1.0 / (nx * ny * nz) * std::sqrt(abs2(sigbar11(0, 0, 0)) + abs2(sigbar22(0, 0, 0)) + abs2(sigbar33(0, 0, 0)) + abs2(sigbar12(0, 0, 0)) + abs2(sigbar23(0, 0, 0)) + abs2(sigbar31(0, 0, 0)));
        rel = std::abs(norm_sig - norm_old) / norm_old;

        // divergence of stress
        double tmp = 0.0;
        #pragma omp parallel for reduction(+:tmp)
        for (uint i = 0; i < nx; ++i) {
            for (uint j = 0; j < ny; ++j) {
                for (uint k = 0; k < nz / 2 + 1; ++k) {
                    /* const double k0 = 2.0*M_PI*ksi1[i]; */
                    /* const double k1 = 2.0*M_PI*ksi2[j]; */
                    /* const double k2 = 2.0*M_PI*ksi3[k]; */
                    const double cc0 = cos(M_PI * ksi1[i]);
                    const double ss0 = sin(M_PI * ksi1[i]);
                    const double cc1 = cos(M_PI * ksi2[j]);
                    const double ss1 = sin(M_PI * ksi2[j]);
                    const double cc2 = cos(M_PI * ksi3[k]);
                    const double ss2 = sin(M_PI * ksi3[k]);
                    const double k0 = ss0 * cc1 * cc2;
                    const double k1 = cc0 * ss1 * cc2;
                    const double k2 = cc0 * cc1 * ss2;
                    const std::complex<double> s1 = sigbar11(i, j, k);
                    const std::complex<double> s2 = sigbar22(i, j, k);
                    const std::complex<double> s3 = sigbar33(i, j, k);
                    const std::complex<double> s4 = sigbar12(i, j, k);
                    const std::complex<double> s5 = sigbar23(i, j, k);
                    const std::complex<double> s6 = sigbar31(i, j, k);
                    const std::complex<double> div1 = k0 * s1 + k1 * s4 + k2 * s6;
                    const std::complex<double> div2 = k0 * s4 + k1 * s2 + k2 * s5;
                    const std::complex<double> div3 = k0 * s6 + k1 * s5 + k2 * s3;
                    tmp += abs2(div1) + abs2(div2) + abs2(div3);
                }
            }
        }

        diverr = std::sqrt(tmp)/(nx*ny*nz*norm_sig);
        if (niter > 0) {
            if ((rel < reltol) || (diverr < divtol)) {
                converged = true;
            }
        }
        norm_old = norm_sig;
        conv_cum_time += get_wall_time() - t4;


        // --- Apply Green operator -------------------------------------------------------

        double t2 = get_wall_time();

        #pragma omp parallel for
        for (uint i = 0; i < nx; ++i) {
            for (uint j = 0; j < ny; ++j) {
                for (uint k = 0; k < nz / 2 + 1; ++k) {

                    if (i == 0 && j == 0 && k == 0) {
                        etobar11(0, 0, 0) = nx * ny * nz * EMACRO[0];
                        etobar22(0, 0, 0) = nx * ny * nz * EMACRO[1];
                        etobar33(0, 0, 0) = nx * ny * nz * EMACRO[2];
                        etobar12(0, 0, 0) = nx * ny * nz * EMACRO[3];
                        etobar23(0, 0, 0) = nx * ny * nz * EMACRO[4];
                        etobar31(0, 0, 0) = nx * ny * nz * EMACRO[5];
                    } else {

                        double g00, g01, g02, g03, g04, g05, g11, g12, g13, g14, g15, g22,
                               g23, g24, g25, g33, g34, g35, g44, g45, g55;

                        compute_green(ksi1[i], ksi2[j], ksi3[k], lmbda0, mu0,
                                      g00, g01, g02, g03, g04, g05,
                                      g11, g12, g13, g14, g15,
                                      g22, g23, g24, g25,
                                      g33, g34, g35,
                                      g44, g45,
                                      g55);

                        const std::complex<double> s1 = sigbar11(i, j, k);
                        const std::complex<double> s2 = sigbar22(i, j, k);
                        const std::complex<double> s3 = sigbar33(i, j, k);
                        const std::complex<double> s4 = sigbar12(i, j, k);
                        const std::complex<double> s5 = sigbar23(i, j, k);
                        const std::complex<double> s6 = sigbar31(i, j, k);

                        etobar11(i, j, k) -= (g00 * s1 + g01 * s2 + g02 * s3 + g03 * s4 + g04 * s5 + g05 * s6);
                        etobar22(i, j, k) -= (g01 * s1 + g11 * s2 + g12 * s3 + g13 * s4 + g14 * s5 + g15 * s6);
                        etobar33(i, j, k) -= (g02 * s1 + g12 * s2 + g22 * s3 + g23 * s4 + g24 * s5 + g25 * s6);
                        etobar12(i, j, k) -= (g03 * s1 + g13 * s2 + g23 * s3 + g33 * s4 + g34 * s5 + g35 * s6);
                        etobar23(i, j, k) -= (g04 * s1 + g14 * s2 + g24 * s3 + g34 * s4 + g44 * s5 + g45 * s6);
                        etobar31(i, j, k) -= (g05 * s1 + g15 * s2 + g25 * s3 + g35 * s4 + g45 * s5 + g55 * s6);
                    }
                }
            }
        }
        green_cum_time += get_wall_time() - t2;

        // --- FFT backward -------------------------------------------------------
        double t3 = get_wall_time();
        E11.backward_fft();
        E22.backward_fft();
        E33.backward_fft();
        E12.backward_fft();
        E23.backward_fft();
        E31.backward_fft();
        S11.backward_fft();
        S22.backward_fft();
        S33.backward_fft();
        S12.backward_fft();
        S23.backward_fft();
        S31.backward_fft();
        fft_cum_time += get_wall_time() - t3;


        // --- Info -------------------------------------------------------

        if (!(niter % 5)) {
            std::cout << std::scientific <<  niter << "  " << rel << "   " << diverr << " | "
                      << average(sig11) << "  " << average(sig22) << "  " << average(sig33) << "  "
                      << average(sig12) << "  " << average(sig23) << "  " << average(sig31) << '\n';
        }

        if (converged) {
            std::cout << "Converged in " << niter << " iterations (diverr = " << diverr << ")\n";
            break;
        }

    }


    double tfinal = get_wall_time() - tinit;
    std::cout << std::fixed;
    std::cout << "Total time:       " << tfinal << " s\n";
    std::cout << "   FFT time:      " << fft_cum_time      << " s  (" << fft_cum_time / tfinal       << ") \n";
    std::cout << "   Behavior time: " << behavior_cum_time << " s  (" << behavior_cum_time / tfinal  << ") \n";
    std::cout << "   Green time:    " << green_cum_time    << " s  (" << green_cum_time / tfinal     << ") \n";
    std::cout << "   Conv. time:    " << conv_cum_time    << " s  (" << conv_cum_time / tfinal     << ") \n";



    // --- Output -----
    //
    //

    std::cout << std::scientific ;
    std::cout << "<E11>  <E22>  <E33>  <E12>  <E23>  <E13>" << std::endl;
    std::cout << average(eto11) << "  " << average(eto22) << "  " << average(eto33) << "  " << average(eto12) << "  " << average(eto23) << "  " << average(eto31) << "\n";
    
    std::cout << "<S11>       <S22>  <S33>  <S12>  <S23>  <S13>" << std::endl;
    std::cout << average(sig11) << "  " << average(sig22) << "  " << average(sig33) << "  " << average(sig12) << "  " << average(sig23) << "  " << average(sig31) << "\n";
    

    std::ofstream res;
    res.open("res.out");
    res << std::scientific << average(sig11) << "  " << average(sig22) << "  " << average(sig33) << "  " << average(sig12) << "  " << average(sig23) << "  " << average(sig31) << "\n";
    res.close();

    std::ofstream output;
    std::string filename = parser.get<std::string>("output");

    output.open(filename.c_str(), std::ios::binary);
    output << "<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\">" << std::endl;
    output << "  <ImageData WholeExtent=\"0 " << nx << " 0 " << ny << " 0 " << nz << "\" Origin=\"0 0 0\" Spacing=\"1 1 1\">" << std::endl;
    output << "    <Piece Extent=\"0 " << nx << " 0 " << ny << " 0 " << nz << "\">" << std::endl;

    output << "      <CellData>" << std::endl;

    int offset = 0;

    output << "        <DataArray Name=\"phase\" type=\"Int32\" format=\"appended\" offset=\"" << offset << "\" />" << std::endl;
    offset += nx*ny*nz*sizeof(int) + 4;

    output << "        <DataArray Name=\"E11\" type=\"Float64\" format=\"appended\" offset=\"" << offset << "\" />" << std::endl;
    offset += nx*ny*nz*sizeof(double) + 4;
    output << "        <DataArray Name=\"E22\" type=\"Float64\" format=\"appended\" offset=\"" << offset << "\" />" << std::endl;
    offset += nx*ny*nz*sizeof(double) + 4;
    output << "        <DataArray Name=\"E33\" type=\"Float64\" format=\"appended\" offset=\"" << offset << "\" />" << std::endl;
    offset += nx*ny*nz*sizeof(double) + 4;
    output << "        <DataArray Name=\"E12\" type=\"Float64\" format=\"appended\" offset=\"" << offset << "\" />" << std::endl;
    offset += nx*ny*nz*sizeof(double) + 4;
    output << "        <DataArray Name=\"E23\" type=\"Float64\" format=\"appended\" offset=\"" << offset << "\" />" << std::endl;
    offset += nx*ny*nz*sizeof(double) + 4;
    output << "        <DataArray Name=\"E31\" type=\"Float64\" format=\"appended\" offset=\"" << offset << "\" />" << std::endl;
    offset += nx*ny*nz*sizeof(double) + 4;

    output << "        <DataArray Name=\"S11\" type=\"Float64\" format=\"appended\" offset=\"" << offset << "\" />" << std::endl;
    offset += nx*ny*nz*sizeof(double) + 4;
    output << "        <DataArray Name=\"S22\" type=\"Float64\" format=\"appended\" offset=\"" << offset << "\" />" << std::endl;
    offset += nx*ny*nz*sizeof(double) + 4;
    output << "        <DataArray Name=\"S33\" type=\"Float64\" format=\"appended\" offset=\"" << offset << "\" />" << std::endl;
    offset += nx*ny*nz*sizeof(double) + 4;
    output << "        <DataArray Name=\"S12\" type=\"Float64\" format=\"appended\" offset=\"" << offset << "\" />" << std::endl;
    offset += nx*ny*nz*sizeof(double) + 4;
    output << "        <DataArray Name=\"S23\" type=\"Float64\" format=\"appended\" offset=\"" << offset << "\" />" << std::endl;
    offset += nx*ny*nz*sizeof(double) + 4;
    output << "        <DataArray Name=\"S31\" type=\"Float64\" format=\"appended\" offset=\"" << offset << "\" />" << std::endl;
    offset += nx*ny*nz*sizeof(double) + 4;

    output << "      </CellData>" << std::endl;


    output << "    </Piece>" << std::endl;
    output << "  </ImageData>" << std::endl;


    output << "  <AppendedData encoding=\"raw\">_";

    unsigned int data_size;
    data_size = nx * ny * nz * sizeof(int);

    output.write(reinterpret_cast<char*>(&data_size), sizeof(unsigned int));
    for(uint k = 0; k < nz; ++k) {
        for (uint j = 0; j < ny; ++j) {
            for (uint i = 0; i < nx; ++i) {
                int tmp = IM(i,j,k);
                output.write(reinterpret_cast<char *>(&tmp), sizeof(int));
            }
        }
    }

    write_binary_data(output, eto11);
    write_binary_data(output, eto22);
    write_binary_data(output, eto33);
    write_binary_data(output, eto12);
    write_binary_data(output, eto23);
    write_binary_data(output, eto31);

    write_binary_data(output, sig11);
    write_binary_data(output, sig22);
    write_binary_data(output, sig33);
    write_binary_data(output, sig12);
    write_binary_data(output, sig23);
    write_binary_data(output, sig31);

    output << "</AppendedData>\n";


    output << "</VTKFile>" << std::endl;
    output.close();

#define WRITE_RAW(s)  filename = #s".raw"; output.open(filename.c_str(), std::ios::binary); write_raw(output, s); output.close();
   WRITE_RAW(sig11);
   WRITE_RAW(sig22);
   WRITE_RAW(sig33);
   WRITE_RAW(sig12);
   WRITE_RAW(sig23);
   WRITE_RAW(sig31);
#undef WRITE_RAW

    return 0;
}


void write_binary_data( std::ofstream& output, fft_real_view& data) {

    typedef unsigned int uint;
    int nx = data.size(0);
    int ny = data.size(1);
    int nz = data.size(2);

    unsigned int  data_size = nx * ny * nz * sizeof(double);
    output.write(reinterpret_cast<char*>(&data_size), sizeof(unsigned int));

    for(uint k = 0; k < nz; ++k) {
        for (uint j = 0; j < ny; ++j) {
            for (uint i = 0; i < nx; ++i) {
                double tmp = data(i,j,k);
                output.write(reinterpret_cast<char *>(&tmp), sizeof(double));
            }
        }
    }
}


void write_raw( std::ofstream& output, fft_real_view& data) {

    typedef unsigned int uint;
    int nx = data.size(0);
    int ny = data.size(1);
    int nz = data.size(2);

    for(uint k = 0; k < nz; ++k) {
        for (uint j = 0; j < ny; ++j) {
            for (uint i = 0; i < nx; ++i) {
                double tmp = data(i,j,k);
                output.write(reinterpret_cast<char *>(&tmp), sizeof(double));
            }
        }
    }

}


double average(fft_real_view &f) {
    int size = f.size(0) * f.size(1) * f.size(2);
    double ave = 0.0;
    for (int i = 0; i < f.size(0); ++i) {
        for (int j = 0; j < f.size(1); ++j) {
            for (int k = 0; k < f.size(2); ++k) {
                ave += f(i, j, k);
            }
        }
    }
    return ave / size;
}

void init_freq(std::vector<double> &ksi1, std::vector<double> &ksi2,
               std::vector<double> &ksi3) {

    typedef unsigned int uint;

    uint n1 = ksi1.size();
    uint n2 = ksi2.size();
    uint n3 = ksi3.size();

    for (uint i = 0; i < n1 / 2; ++i) {
        ksi1[i] = 1.0 * i / n1;
    }
    for (uint i = n1 / 2; i > 0; --i) {
        ksi1[n1 - i] = -1.0 * i / n1;
    }

    for (uint i = 0; i < n2 / 2; ++i) {
        ksi2[i] = 1.0 * i / n2;
    }
    for (uint i = n2 / 2; i > 0; --i) {
        ksi2[n2 - i] = -1.0 * i / n2;
    }

    for (uint i = 0; i < n3; ++i) {
        ksi3[i] = 0.5 * i / (n3-1);
    }
}





