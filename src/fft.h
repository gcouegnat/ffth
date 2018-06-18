#ifndef FFTIMAGE_H_NPTANXPJ
#define FFTIMAGE_H_NPTANXPJ

#include <complex>
#include <fftw3.h>
#include <fstream>
#include <iostream>

#ifdef WITH_MKL
#include<mkl.h>
#endif


class fft {

    public:
        fft(int ni, int nj, int nk);
        ~fft();

        int size(int i) { return isize[i]; }
        double* real_buffer() { return buffer; }
        std::complex<double>* complex_buffer() { return reinterpret_cast<std::complex<double>*>(buffer); }

        void forward_fft();
        void backward_fft();

    private:
        int buffer_size;
        double *buffer;
        int isize[3];
        
        fftw_plan plan_r2c;
        fftw_plan plan_c2r;
};

class fft_real_view  {

    public:
        fft_real_view(fft& in);
        inline double  operator() (int i, int j, int k) const;
        inline double& operator() (int i, int j, int k);

        int size(int i) { return isize[i]; }

    private:
        double *buffer;
        int isize[3];
        int osize[3];

};



class fft_complex_view {
    public:
        fft_complex_view(fft& in);
        inline std::complex<double>  operator() (int i, int j, int k) const;
        inline std::complex<double>& operator() (int i, int j, int k);

    private:
        std::complex<double> *buffer;
        int isize[3];
        int osize[3];
};






// ----------------------------------

fft::fft(int ni, int nj, int nk) {

    buffer_size = ni * nj * 2 * (nk/2 + 1);
    /* buffer = fftw_alloc_real(buffer_size); */
    buffer = (double *) fftw_malloc(sizeof(double) * buffer_size);
    
    isize[0] = ni;
    isize[1] = nj;
    isize[2] = nk;

    plan_r2c = fftw_plan_dft_r2c_3d(ni, nj, nk,
            buffer, reinterpret_cast<fftw_complex*>(buffer),
            FFTW_ESTIMATE);
    plan_c2r = fftw_plan_dft_c2r_3d(ni, nj, nk,
            reinterpret_cast<fftw_complex*>(buffer), buffer,
            FFTW_ESTIMATE);
}

fft::~fft() {
    fftw_destroy_plan(plan_r2c);
    fftw_destroy_plan(plan_c2r);
    fftw_free(buffer);
}

void fft::forward_fft() {
    fftw_execute(plan_r2c);
}

void fft::backward_fft() {
    fftw_execute(plan_c2r);
    // normalize output
    const int num = isize[0]*isize[1]*isize[2];
    double inv = 1.0 / num;
    
    for(int i = 0; i < buffer_size; ++i) {
        buffer[i] *= inv;
    }

}

// ----------------------------------

fft_real_view::fft_real_view(fft& in) {

    buffer = in.real_buffer();

    isize[0] = in.size(0);
    isize[1] = in.size(1);
    isize[2] = in.size(2);
    
    osize[0] = in.size(0);
    osize[1] = in.size(1);
    osize[2] = 2*(in.size(2)/2+1);
}

double fft_real_view::operator() (int i, int j, int k) const {
    return buffer[i * osize[1] * osize[2] + j * osize[2] + k];
}

double& fft_real_view::operator() (int i, int j, int k) {
    return buffer[i * osize[1] * osize[2] + j * osize[2] + k];
}

// ----------------------------------

fft_complex_view::fft_complex_view(fft& in) {
    buffer = in.complex_buffer();
    isize[0] = in.size(0);
    isize[1] = in.size(1);
    isize[2] = in.size(2);
    osize[0] = in.size(0);
    osize[1] = in.size(1);
    osize[2] = in.size(2)/2+1;
}

std::complex<double> fft_complex_view::operator() (int i, int j, int k) const {
    return buffer[i * osize[1] * osize[2] + j * osize[2] + k];
}

std::complex<double>& fft_complex_view::operator() (int i, int j, int k) {
    return buffer[i * osize[1] * osize[2] + j * osize[2] + k];
}

#endif /* end of include guard: FFTIMAGE_H_NPTANXPJ */
