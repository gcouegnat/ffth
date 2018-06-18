CXXFLAGS=-O3 -fopenmp

ffth: src/main.cpp
	${CXX} ${CXXFLAGS} $< -o $@ -I${EIGEN_INC} -I${FFTW_INC} -L${FFTW_LIB} -lfftw3_omp -lfftw3

clean:
	rm ffth
