#include <iostream>

#pragma once

#include <cmath>
#include <vector>
#include <complex>
#include <random>

using std::vector;
using std::complex;

class FFT {
    static vector<complex<double>>
    slicing(const std::vector<std::complex<double>> &vec, size_t X, size_t Y, size_t stride) {
        vector<complex<double>> result;
        int i = X;
        while (result.size() < Y) {
            result.push_back(vec[i]);
            i = i + stride;
        }
        return result;
    }

public:

    void fft(vector<complex<double>> &x) {
        const size_t N = x.size();
        if (N <= 1)
            return;

        else if (N % 2 == 0)    //Radix-2
        {
            vector<complex<double>> even = slicing(x, 0, N / 2, 2);
            vector<complex<double>> odd = slicing(x, 1, N / 2, 2);

            fft(even);
            fft(odd);

            for (size_t k = 0; k < N / 2; ++k) {
                complex<double> t = std::polar<double>(1.0, -2 * M_PI * k / N) * odd[k];
                x[k] = even[k] + t;
                x[k + N / 2] = even[k] - t;
            }
        } else if (N % 3 == 0)    //Radix-3
        {
            vector<complex<double>> p0 = slicing(x, 0, N / 3, 3);
            vector<complex<double>> p1 = slicing(x, 1, N / 3, 3);
            vector<complex<double>> p2 = slicing(x, 2, N / 3, 3);

            fft(p0);
            fft(p1);
            fft(p2);

            for (int i = 0; i < N; i++) {
                complex<double> temp = p0[i % ((int) N / 3)];
                temp += (p1[i % ((int) N / 3)] * std::polar<double>(1.0, -2 * M_PI * i / N));
                temp += (p2[i % ((int) N / 3)] * std::polar<double>(1.0, -4 * M_PI * i / N));
                x[i] = temp;
            }
        } else if (N % 5 == 0)    //Radix-5
        {
            vector<complex<double>> p0 = slicing(x, 0, N / 5, 5);
            vector<complex<double>> p1 = slicing(x, 1, N / 5, 5);
            vector<complex<double>> p2 = slicing(x, 2, N / 5, 5);
            vector<complex<double>> p3 = slicing(x, 3, N / 5, 5);
            vector<complex<double>> p4 = slicing(x, 4, N / 5, 5);

            fft(p0);
            fft(p1);
            fft(p2);
            fft(p3);
            fft(p4);
            for (int i = 0; i < N; i++) {
                complex<double> temp = p0[i % (int) N / 5];
                temp += (p1[i % (int) N / 5] * std::polar<double>(1.0, -2 * M_PI * i / N));
                temp += (p2[i % (int) N / 5] * std::polar<double>(1.0, -4 * M_PI * i / N));
                temp += (p3[i % (int) N / 5] * std::polar<double>(1.0, -6 * M_PI * i / N));
                temp += (p4[i % (int) N / 5] * std::polar<double>(1.0, -8 * M_PI * i / N));
                x[i] = temp;
            }
        }
    }

    void ifft(vector<complex<double>> &x) {
        for (auto &xi: x)
            xi = std::conj(xi);

        fft(x);

        for (auto &xi: x)
            xi = std::conj(xi) / static_cast<double>(x.size());;
    }
};


int main() {
    srand(time(NULL));

    int powerOfTwo = 1 << (rand() % 11 + 1);
    int powerOfThree = pow(3, rand() % 5 + 1);
    int powerOfFive = pow(5, rand() % 7 + 1);
    int variables[] = {powerOfTwo, powerOfThree, powerOfFive};
    int randomIndex = rand() % 3;
    int randomSize = variables[randomIndex];

    vector<complex<double>> arr(randomSize, 0);
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(-1.0, 1.0);
    for (size_t i = 0; i < arr.size(); ++i) {
        double real = distribution(generator);
        double imag = distribution(generator);
        arr[i] = {real, imag};
    }

    auto inp = arr;

    std::cout << "Before:" << std::endl;
    for (auto i: arr)
        std::cout << i << std::endl;

    std::cout << "After:" << std::endl;
    FFT f;
    f.fft(arr);
    f.ifft(arr);
    for (auto i: arr)
        std::cout << i << std::endl;


    std::cout << "Diff:" << std::endl;
    vector<complex<double>> diff(arr.size());
    for (size_t i = 0; i < arr.size(); ++i) {
        diff[i] = arr[i] - inp[i];
        std::cout << diff[i] << std::endl;
    }

    double norm = 0.0;
    for (auto i: diff)
        norm += std::norm(i);

    std::cout << "==============================================" << std::endl;
    std::cout << "Random size : " << randomSize << std::endl;
    std::cout << "Total Diff: " << std::sqrt(norm) << std::endl;

    return 0;
}
