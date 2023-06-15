#include "getspec.h"

#include <complex>
#include <cmath> 
#include <random>
#include <iostream>


#define ONE_OVER_SQRT_8 3.535533905932737e-01


class Matrix {
    private:
        floating* data;
        int dimension;
    public:
        Matrix(int dim) {
            dimension = dim;
            data = new floating[dim*dim];
            for (int i = 0; i < dim; i++) {
                for (int j = 0; j < dim; j++) {
                    data[dim*i + j] = 0.;
                }
            }
        };

        Matrix(const floating* const init, int dim) {
            dimension = dim;
            data = new floating[dim*dim];
            for (int i = 0; i < dim; i++) {
                for (int j = 0; j < dim; j++) {
                    data[dim*i + j] = init[dim*i + j];
                }
            }
        };

        Matrix(const Matrix& init) {
            dimension = init.dimension;
            data = new floating[dimension*dimension];
            for (int i = 0; i < dimension; i++) {
                for (int j = 0; j < dimension; j++) {
                    data[dimension*i + j] = init.data[dimension*i + j];
                }
            }
        };

        Matrix(Matrix&& init) {
            dimension = init.dimension;
            data = init.data;
            init.data = nullptr;
            init.dimension = 0;
        };

        ~Matrix() {
            delete [] data;
        };

        Matrix& operator=(const Matrix& other) {
            if (this == &other)
                return *this;

            if (dimension != other.dimension)           // resource in *this cannot be reused
            {
                delete [] data;
                data = new floating[other.dimension];   // allocate resource, if throws, do nothing
                dimension = other.dimension;
            } 
            for (int i = 0; i < dimension; i++) {
                for (int j = 0; j < dimension; j++) {
                    data[dimension*i + j] = other.data[dimension*i + j];
                }
            }
            
            return *this;
        };

        Matrix& operator=(Matrix&& other) {
            if (this == &other)
                return *this;

            delete [] data;
            data = other.data;
            other.data = nullptr;
            dimension = other.dimension;
            other.dimension = 0;
            
            return *this;
        };

        Matrix& operator+=(const Matrix& rhs) {
            if (rhs.dimension != dimension) {
                throw("DIMENSION INCOMPATIBLE ADD");
                return *this;
            }
            for (int i = 0; i < dimension; i++) {
                for (int j = 0; j < dimension; j++) {
                    data[dimension*i + j] += rhs.data[dimension*i + j];
                }
            }
            return *this;
        };

        friend Matrix operator+(Matrix lhs, const Matrix& rhs) {
            lhs += rhs;
            return lhs;
        };

        Matrix& operator*=(const Matrix& rhs) {
            if (rhs.dimension != dimension) {
                throw("DIMENSION INCOMPATIBLE MULT");
                return *this;
            }
            floating* res = new floating[dimension*dimension];
            for (int i = 0; i < dimension; i++) {
                for (int j = 0; j < dimension; j++) {
                    res[dimension*i + j] = 0.;
                    for (int k = 0; k < dimension; k++) {
                        res[dimension*i + j] += data[dimension*i + k] * rhs.data[dimension*k + j];
                    }
                }
            };
            for (int i = 0; i < dimension; i++) {
                for (int j = 0; j < dimension; j++) {
                    data[dimension*i + j] = res[dimension*i + j];
                }
            };
            delete [] res;
            return *this;
        }

        friend Matrix operator*(Matrix lhs, const Matrix& rhs) {
            lhs *= rhs;
            return lhs;
        };

        Matrix& operator*=(const floating rhs) {
            for (int i = 0; i < dimension; i++) {
                for (int j = 0; j < dimension; j++) {
                    data[dimension*i + j] *= rhs;
                }
            };
            return *this;
        };

        friend Matrix operator*(Matrix lhs, const floating rhs) {
            lhs *= rhs;
            return lhs;
        };

        friend Matrix operator*(const floating lhs, Matrix rhs) {
            rhs *= lhs;
            return rhs;
        };

        floating trace() {
            floating res = 0.;
            for (int i = 0; i < dimension; i++) {
                res += data[dimension*i + i];
            };
            return res;
        };
};


// Hard-coded FFT on 2^n dots

void fourier_transform(cfloating* data, cfloating* twist, int n) {
    unsigned int xxor = 1u << n-1;
    unsigned int step = 1u;
    
    for (int it = 0; it < n; it++) {
        int base = 0;
        for (int i = 1; i <= step; i++) {
            for (int ind = 0; ind < xxor; ind++) {
                cfloating tmp = data[base+ind];
                data[base+ind] += data[base+ind+xxor];
                data[base+ind+xxor] = tmp - data[base+ind+xxor];
                data[base+ind+xxor] *= twist[ind*step];
            }
            base += 2*xxor;
        }
        xxor = xxor >> 1;
        step = step << 1;
    }
}


// Some memory managment functions for python

cfloating* reserve_complex_array(int size) {
    cfloating* ptr = new cfloating [size];
    return ptr;
};

floating* reserve_array(int size) {
    floating* ptr = new floating [size];
    return ptr;
};

cfloating* reserve_twists(int Npow) {
    cfloating* twist = new cfloating[1u << (Npow - 1)];
    twist[0] = 1.; 
    for (int j = 1; j < (1u << (Npow-1)); j++) {
        twist[j] = std::exp(cfloating(0., -M_PI/(1u << (Npow-1))*j));
    }
    return twist;
}; 

// internally shifts fft, so zero frequency now at 2^(Npow-1) place
unsigned int* reserve_bitrev(int Npow) {
    unsigned int* bitrev = new unsigned int [1u << Npow];
    for (int j = 0; j < (1u << Npow); j++) {
        bitrev[j] = 0;
        unsigned int tmpui = j;
        for(int k = Npow-1; k >= 0; k--){
         bitrev[j] |= (tmpui & 1) << k;
         tmpui>>=1;
        }
        bitrev[j] = bitrev[j] ^ (1u << (Npow-1));
    }
    return bitrev;
};

void free_complex_array(cfloating* ptr) {
    delete [] ptr;
};

void free_array(floating* ptr) {
    delete [] ptr;
};

void free_uint_array(unsigned int* ptr) {
    delete [] ptr;
};


// Main function

bool sdeeval(floating* specx, floating* specy, cfloating* Ex, cfloating* Ey, cfloating* tmpEx, cfloating* tmpEy,  // arrays for SDE samples and output spectras 
             cfloating* twist, unsigned int* bitrev, int Npow, int skip, int Nav, int tauDt, floating Dt,            // some necessary variables
             floating alpha, floating kappa, floating gamma, floating gamma_d, floating gamma_a,                  // parameters
             floating gamma_p, floating beta, floating mu, floating C_sp, floating N_rel) {

    // TODO:
    // + return boolean value: true if success, false if there was NaN
    // + Adaptive Dt step in order to avoid negative Qp/Qm
    // * Test different approaches to RND sampling


    // define useful constants
    floating M = N_rel/(1.0 - N_rel);
    int L = (1u << Npow);

    // variables for SDE solver run   
    floating Qp = 0.;
    floating Qm = 0.;
    floating phi = 0.;
    floating psi = 0.;
    floating G = 1.;
    floating d = 0.;
    floating Qp_prev = 0.4892;
    floating Qm_prev = 0.4892;
    floating phi_prev = 0;
    floating psi_prev = 0.;
    floating G_prev = 1.0110;
    floating d_prev = 0.;

    // Initialization based on stability
    // TODO - refactor this!
    floating Q = (-gamma_a + kappa*(mu-1.f+2.f*M*C_sp) + std::sqrt(4.f*(2.f*C_sp-1.f)*kappa*mu*(gamma_a+kappa) + (gamma_a+kappa*(1.f+2.f*C_sp*M+mu))*(gamma_a+kappa*(1.f+2.f*C_sp*M+mu))) )/(4.f*(gamma_a+kappa));

    G_prev = mu/(1.f + 2.f*Q);
    Qp_prev = Q;
    Qm_prev = Q;

    floating Lmat[3*3] = {2.f*kappa*(G_prev-1.f), -8.f*Q*gamma_p, 4.f*kappa*(C_sp+Q),
                          gamma_p/2.f/Q, 2.f*gamma_a, alpha*kappa,
                          -G_prev*gamma, 0.f, -gamma_d-2.f*gamma*Q};

    floating charpoly[4] = {1, 0, 0, 0};
    floating Mat[3*3] = {0, 0, 0, 
                       0, 0, 0, 
                       0, 0, 0};
    floating Idenmat[9] = {1,0,0,0,1,0,0,0,1};
    Matrix LLmat = Matrix(Lmat, 3);
    Matrix MMmat = Matrix(Mat, 3);
    for (int k = 1; k < 4; k++) {
        MMmat = LLmat*MMmat + charpoly[k-1]*Matrix(Idenmat,3);
        charpoly[k] = - 1./k * ((LLmat*MMmat).trace());
    };

    if (charpoly[1] < 0 || charpoly[3] < 0 || charpoly[1]*charpoly[2] - charpoly[3] < 0) {
        floating Lmat2[9] = {2.f*kappa*(G_prev-1.f), 8.f*Q*gamma_p, 4.f*kappa*(C_sp+Q),
                             -gamma_p/2.f/Q, -2.f*gamma_a, alpha*kappa,
                             -G_prev*gamma, 0.f, -gamma_d-2.f*gamma*Q};
        Matrix LLmat2 = Matrix(Lmat2, 3);
        MMmat = Matrix(Mat, 3);
        floating charpoly2[4] = {1, 0, 0, 0};
        for (int k = 1; k < 4; k++) {
            MMmat = LLmat2*MMmat + charpoly2[k-1]*Matrix(Idenmat,3);
            charpoly2[k] = - 1./k * ((LLmat2*MMmat).trace());
        };
        if (charpoly2[1] < 0 || charpoly2[3] < 0 || charpoly2[1]*charpoly2[2] - charpoly2[3] < 0) {
            std::cout << "UNSTABLE SYSTEM\n";
        } else {
            psi_prev = M_PI_2;
        };
    };

    // clearing spectra arrays
    for (int i = 0; i < L; i++) {
        specx[i] = 0.;
        specy[i] = 0.;
    };

    // initializing random generator
    std::random_device rd {};
    std::mt19937 gen {rd()};

    // some cache
    floating sqrtDt = std::sqrt(Dt);
    floating sqrtkappa = std::sqrt(kappa);
    constexpr cfloating M_I = cfloating(0., 1.);

    // some memory allocations
    floating U;
    floating V;
    floating ksi_plus;
    floating ksi_minus;
    floating ksi_phi;
    floating ksi_psi;
    floating C_plus;
    floating C_minus;
    floating sqrtQp;
    floating sqrtQm;
    floating Fp;
    floating Fm;
    floating Fphi;
    floating Fpsi;
    floating sqrtQpQm;

    cfloating Epsq2;
    cfloating Emsq2;


    int start = 0;
    for (int av = 0; av < Nav; av++) {
        for (int i = start; i < L; i++) {
            for (int j = 0; j < tauDt; j++) {
                // Evaluating SDE

                // obtaining rnd variables
                U = (floating)gen() / gen.max();
                V = (floating)gen() / gen.max();
                ksi_plus = std::sqrt(-2.*std::log(U));
                ksi_minus = ksi_plus * std::sin(2*M_PI*V);
                ksi_plus *= std::cos(2*M_PI*V);
                U = (floating)gen() / gen.max();
                V = (floating)gen() / gen.max();
                ksi_phi = std::sqrt(-2.*std::log(U));
                ksi_psi = ksi_phi * std::sin(2*M_PI*V);
                ksi_phi *= std::cos(2*M_PI*V);

                // cache for forces
                C_plus = std::sqrt(C_sp*(G_prev+d_prev+M));
                C_minus = std::sqrt(C_sp*(G_prev-d_prev+M));
                sqrtQp = std::sqrt(Qp_prev);
                sqrtQm = std::sqrt(Qm_prev);
                // obtaining noise forces
                Fp = 2.*sqrtkappa*sqrtQp*C_plus*ksi_plus;
                Fm = 2.*sqrtkappa*sqrtQm*C_minus*ksi_minus;
                Fphi = ONE_OVER_SQRT_8*C_plus/sqrtQp*(ksi_phi+ksi_psi) + 
                                ONE_OVER_SQRT_8*C_minus/sqrtQm*(ksi_phi-ksi_psi);
                Fpsi = ONE_OVER_SQRT_8*C_plus/sqrtQp*(ksi_phi+ksi_psi) - 
                                ONE_OVER_SQRT_8*C_minus/sqrtQm*(ksi_phi-ksi_psi);

                // performing step of Euler-Maruyama method
                sqrtQpQm = sqrtQp*sqrtQm;
                Qp  = Qp_prev + 2.*( kappa*(G_prev+d_prev-1)*Qp_prev + kappa*C_sp*(G_prev+d_prev+M) - 
                    sqrtQpQm*(gamma_a*std::cos(2.*psi_prev-2*beta) + gamma_p*std::sin(2.*psi_prev) ) )*Dt + 
                    Fp*sqrtDt;
                Qm  = Qm_prev + 2.*( kappa*(G_prev-d_prev-1.)*Qm_prev + kappa*C_sp*(G_prev-d_prev+M) - 
                    sqrtQpQm*(gamma_a*std::cos(2.*psi_prev-2*beta) - gamma_p*std::sin(2.*psi_prev) ) )*Dt + 
                    Fm*sqrtDt;
                phi = phi_prev + ( (G_prev-1)*alpha*kappa - 
                    0.5/sqrtQpQm*( (Qp_prev+Qm_prev)*gamma_p*std::cos(2.*psi_prev) + 
                                   (Qp_prev-Qm_prev)*gamma_a*std::sin(2.*psi_prev-2*beta) ))*Dt + 
                    Fphi*sqrtDt;
                psi = psi_prev + (alpha*kappa*d_prev + 
                    0.5/sqrtQpQm*( (Qp_prev-Qm_prev)*gamma_p*std::cos(2.*psi_prev) +
                                   (Qp_prev+Qm_prev)*gamma_a*std::sin(2.*psi_prev-2*beta) ))*Dt + 
                    Fpsi*sqrtDt;
                G   = G_prev + gamma*Dt*( (mu-G_prev) - G_prev*(Qp_prev + Qm_prev) - d_prev*(Qp_prev - Qm_prev) );
                d   = d_prev + Dt*(-gamma_d*d_prev -gamma*G_prev*(Qp_prev - Qm_prev) - gamma*d_prev*(Qp_prev + Qm_prev) );
                Qp_prev = Qp;
                Qm_prev = Qm;
                phi_prev = phi;
                psi_prev = psi;
                G_prev = G;
                d_prev = d;
            };

            if (std::isnan(Qp) || std::isnan(Qm) || std::isnan(phi) || std::isnan(psi) || std::isnan(G) || std::isnan(d)) {
                std::cout << "NaN encountered. Aborting computation; av = " << av << ", i = " << i <<"\n";
                return false;
            }

            // evaluating fields in sample point
            Epsq2 = M_SQRT1_2 * std::sqrt(Qp) * std::exp(cfloating(0., phi + psi));
            Emsq2 = M_SQRT1_2 * std::sqrt(Qm) * std::exp(cfloating(0., phi - psi));
            Ex[i] = Epsq2 + Emsq2;
            Ey[i] = M_I * (Emsq2 - Epsq2);
        };

        // FFT part
        // saving necessary previous values
        for (int i = skip; i < L; i++) {
            tmpEx[i-skip] = Ex[i];
            tmpEy[i-skip] = Ey[i];
        };
        // obtain bit-reversed spectra
        fourier_transform(Ex, twist, Npow);
        fourier_transform(Ey, twist, Npow);
        // adding to existing spectra
        for (int i = 0; i < L; i++) {
            specx[bitrev[i]] += std::norm(Ex[i]);
            specy[bitrev[i]] += std::norm(Ey[i]);
        };
        // change Ei <-> tmpEi 
        std::swap<cfloating*>(Ex, tmpEx);
        std::swap<cfloating*>(Ey, tmpEy);
        start = L - skip;
    };

    // averaging spectra
    for (int i = 0; i < L; i++) {
            specx[i] /= Nav;
            specy[i] /= Nav;
            //outSPEC << specx[i] << " " << specy[i] << "\n";
    };

    return true;
};