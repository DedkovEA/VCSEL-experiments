#include <cmath>
#include <random>
#include <fstream>
#include <chrono>
#include <iostream>


#define INTRINSIC_TYPE double
#define ONE_OVER_SQRT_8 3.535533905932737e-01


constexpr float pi = 3.14159265358979323846;

typedef INTRINSIC_TYPE floating;


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


void matrix_mult_to_second(const floating* const left, floating* const right, int n) {
    floating* res = new floating[n*n];
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            res[n*i + j] = 0.;
            for (int k = 0; k < n; k++) {
                res[n*i + j] += left[n*i + k] * right[n*k + j];
            }
        }
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            right[n*i + j] = res[n*i + j];
        }
    }
    delete [] res;
};

void sdeeval() {
    floating alpha = 3.;
    floating kappa = 80.;
    floating gamma = 1.;
    floating gamma_d = 1000.;
    floating gamma_a =  2.5;
    floating gamma_p = 2*3.14159265*9;
    floating beta = 0.;
    floating mu = 6.;

    floating C_sp = 5e-4;
    floating N_th = 6.25e6;
    floating N_tr = 5.935e6;
    floating M = N_tr/(N_th - N_tr);

    //mu = (N_th*3.179f - N_tr)/(N_th - N_tr);

    floating Dt = 1e-6;
    floating tau = 1e-4;
    floating T = 50;

    int tauDt = std::floor(tau/Dt);
    if (tauDt < 1) { tauDt = 1; };
    tau = tauDt * Dt;  // fix tau in order to make it right value
    int N = std::ceil(T/Dt);  // total number of evaluating dots
    int L = N / tauDt;  // number of exporting dots

    
    floating Qp = 0.4892;
    floating Qm = 0.4892;
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

    std::cout << charpoly[0] << " " << charpoly[1] << " " << charpoly[2] << "\n" << charpoly[3] << "\n";

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
            psi_prev = pi/2;
        };
    };

    floating* Qp_out = new floating[L];
    floating* Qm_out = new floating[L];
    floating* phi_out = new floating[L];
    floating* psi_out = new floating[L];
    floating* G_out = new floating[L];
    floating* d_out = new floating[L];

    // initializing random generator
    std::normal_distribution<double> ndistr(0.0, 1.0);
    std::random_device rd {};
    std::mt19937 gen {rd()};

    // some cache
    floating sqrtDt = std::sqrt(Dt);
    floating sqrtkappa = std::sqrt(kappa);


    // just for testing purpose open output file
    std::ofstream out("SDEsolution.txt", std::ios::out);

    for (int i = 0; i < L; i++) {
        for (int j = 0; j < tauDt; j++) {
            // Evaluating SDE

            // obtaining rnd variables

            // OLD WAY
            // floating ksi_plus = ndistr(gen);
            // floating ksi_minus = ndistr(gen);
            // floating ksi_phi = ndistr(gen);
            // floating ksi_psi = ndistr(gen);

            // MUCH FASTER
            floating U = (float)gen() / gen.max();
            floating V = (float)gen() / gen.max();
            floating ksi_plus = std::sqrt(-2.*std::log(U));
            floating ksi_minus = ksi_plus * std::sin(2*pi*V);
            ksi_plus *= std::cos(2*pi*V);
            U = (float)gen() / gen.max();
            V = (float)gen() / gen.max();
            floating ksi_phi = std::sqrt(-2.*std::log(U));
            floating ksi_psi = ksi_phi * std::sin(2*pi*V);
            ksi_phi *= std::cos(2*pi*V);

            // cache for forces
            floating C_plus = std::sqrt(C_sp*(G_prev+d_prev+M));
            floating C_minus = std::sqrt(C_sp*(G_prev-d_prev+M));
            floating sqrtQp = std::sqrt(Qp_prev);
            floating sqrtQm = std::sqrt(Qm_prev);
            // obtaining noise forces
            floating Fp = 2.*sqrtkappa*sqrtQp*C_plus*ksi_plus;
            floating Fm = 2.*sqrtkappa*sqrtQm*C_minus*ksi_minus;
            floating Fphi = ONE_OVER_SQRT_8*C_plus/sqrtQp*(ksi_phi+ksi_psi) + 
                            ONE_OVER_SQRT_8*C_minus/sqrtQm*(ksi_phi-ksi_psi);
            floating Fpsi = ONE_OVER_SQRT_8*C_plus/sqrtQp*(ksi_phi+ksi_psi) - 
                            ONE_OVER_SQRT_8*C_minus/sqrtQm*(ksi_phi-ksi_psi);

            // performing step of Euler-Maruyama method
            floating sqrtQpQm = sqrtQp*sqrtQm;
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

        // Saving sample point
        Qp_out[i] = Qp;
        Qm_out[i] = Qm;
        phi_out[i] = phi;
        psi_out[i] = psi;
        G_out[i] = G;
        d_out[i] = d;
        out << Qp << " " << Qm << " " << phi << " " << psi << " " << G << " " << d << "\n";
    };
};

int main() {
    std::chrono::time_point start = std::chrono::steady_clock::now();
    sdeeval();
    std::chrono::time_point end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";
    return 0;
}