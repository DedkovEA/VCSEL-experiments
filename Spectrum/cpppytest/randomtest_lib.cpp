#include <random>
#include <iostream>
#include <fstream>
#include <chrono>

int main() 
{
    std::ofstream out("output_lib.txt", std::ios::out);
    std::normal_distribution<double> distr(0.0, 1.0);
    std::random_device rd {};
    std::mt19937 gen {rd()};

    std::chrono::time_point start = std::chrono::steady_clock::now();
    double tmp = 0;
    std::cout << (gen.min()) << " " << (gen.max());
    for (int i = 0; i < 100000000; i++) 
    {
        tmp = distr(gen);
        if (tmp > 5.) {
            out << tmp << "\n";
        };
    }
    std::chrono::time_point end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";
}