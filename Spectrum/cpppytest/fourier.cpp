#include <cmath>
#include <complex>
#include <iostream>

typedef double floating;
typedef std::complex<floating> cfloating;

unsigned int pow(int base, int power) {
    unsigned int pw = (unsigned int) power;
    int res = 1;
    while (pw != 0) {
        if (pw % 2)
            res *= base;
        base *= base;
        pw = pw >> 1;
    };
    return res;
};

void fourier_transform(cfloating* data, int n) {
    unsigned int xxor = 1u << n-1;
    unsigned int step = 1u;
    cfloating* twist = new cfloating[xxor];
    unsigned int* bitrev = new unsigned int [1u << n];

    twist[0] = 1; 
    for (int j = 1; j < xxor; j++) {
        twist[j] = std::exp(cfloating(0., -M_PI/xxor*j));
    }

    for (int j = 0; j < (1u << n); j++) {
        bitrev[j] = 0;
        unsigned int tmpui = j;
        for(int k = n-1; k >= 0; k--){
         bitrev[j] |= (tmpui & 1) << k;
         tmpui>>=1;
      }
    }

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

    for (int i = 0; i < (1u << n); i++) {
        if (bitrev[i] > i) {
            std::swap(data[i], data[bitrev[i]]);
        }
    }

    delete [] twist;
    delete [] bitrev;
}

int main() {
    cfloating dat[16] = {1., 2., 3., 4., 5., 6., 7., 8., 9, 10, 11, 12, 13, 14, 15, 16};
    fourier_transform(dat, 4);
    for (unsigned int i = 0; i < 16; i++) {
        std::cout << dat[i] << " ";
    }
    std::cout << '\n';
    return 0;
}