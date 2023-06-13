#include <random>
#include <iostream>
#include <fstream>
#include <chrono>
#include <cmath>

constexpr float pi = 3.14159265358979323846;

int main() 
{
    std::ofstream out("output.txt", std::ios::out);
    std::normal_distribution<double> distr(0.0, 1.0);
    std::uniform_real_distribution<double> udistr(-2., 11.);
    std::random_device rd {};
    std::mt19937 gen {rd()};

    std::chrono::time_point start = std::chrono::steady_clock::now();
    double tmp = 0;
    std::cout << (gen.min()) << " " << (gen.max());
    for (int i = 0; i < 500000; i++) 
    {
        // out << distr(gen) << '\n';
        // tmp = udistr(gen);
        

        float U = (float)gen() / gen.max();
        float V = (float)gen() / gen.max();
        float n1 = std::sqrt(-2.*std::log(U));
        float n2 = n1 * std::sin(2*pi*V);
        n1 *= std::cos(2*pi*V);
        if (n1 > 5.) {
            out << n1 << "\n";
        };
        if (n2 > 5.) {
            out << n2 << "\n";
        };
    }
    std::chrono::time_point end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";
}