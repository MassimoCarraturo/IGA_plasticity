// qi_test.cpp — test harness for the C++ wrapper.
//
// Build:
//   g++ -O2 -std=c++17 -fopenmp -I../include -I../src \
//       qi_test.cpp ../build/libqi.a -fopenmp -lm -o qi_test_cpp
//
// or link against the shared library:
//   g++ -O2 -std=c++17 -fopenmp -I../include -I../src \
//       qi_test.cpp -L../build -lqi -lm -o qi_test_cpp

#include "qi.hpp"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

static double halton(int i, int b)
{
    double f = 1.0, r = 0.0;
    while (i > 0) { f /= b; r += f * (i % b); i /= b; }
    return r;
}

int main()
{
    std::cout << "libqi C++ wrapper -- " << qi::version() << "\n";

    // ----------------- Test 1: polynomial reproduction -----------------
    {
        const int n = 2000;
        std::vector<double> x(n), y(n), f(n);
        for (int i = 0; i < n; ++i) {
            x[i] = halton(i+1, 2);
            y[i] = halton(i+1, 3);
            f[i] = 0.7 - 1.3*x[i] + 0.5*y[i] + 2.1*x[i]*x[i]
                   - 1.7*x[i]*y[i] + 0.9*y[i]*y[i];
        }
        qi::HSpace hs({2, 2}, {{8, 8}});
        auto c = qi::compute_coeffs(hs, x, y, f);

        const int G = 50;
        std::vector<double> gx(G*G), gy(G*G);
        for (int i = 0; i < G; ++i)
            for (int j = 0; j < G; ++j) {
                gx[i*G+j] = double(i)/(G-1);
                gy[i*G+j] = double(j)/(G-1);
            }
        auto fhat = qi::eval_tp(hs, 0, c, gx, gy);
        double err = 0.0;
        for (int i = 0; i < G*G; ++i) {
            double exact = 0.7 - 1.3*gx[i] + 0.5*gy[i] + 2.1*gx[i]*gx[i]
                            - 1.7*gx[i]*gy[i] + 0.9*gy[i]*gy[i];
            err = std::max(err, std::abs(fhat[i] - exact));
        }
        std::cout << "  polynomial reproduction max error = "
                  << std::scientific << std::setprecision(3) << err << "\n";
        if (err >= 1e-9) { std::cerr << "FAIL\n"; return 1; }
    }

    // ----------------- Test 2: smooth + diagnostics --------------------
    {
        const int n = 4000;
        std::vector<double> x(n), y(n), f(n);
        for (int i = 0; i < n; ++i) {
            x[i] = halton(i+1, 2);
            y[i] = halton(i+1, 3);
            double t1 = (std::tanh(9*y[i] - 9*x[i]) + 1)/9;
            double dx = 10*x[i] - 6, dy = 10*y[i] + 7;
            double t2 = std::exp(-(dx*dx + dy*dy))/1.5;
            f[i] = t1 + t2;
        }
        qi::HSpace hs({3, 3}, {{16, 16}});
        qi::Diagnostics diag;
        auto c = qi::compute_coeffs(hs, x, y, f, qi::Options{}, &diag);

        std::cout << "  smooth: ndof=" << hs.ndof()
                  << "  unique-degrees={";
        std::vector<bool> seen(8, false);
        for (int d : diag.deg_used)
            if (d >= 0 && d < 8) seen[d] = true;
        for (int d = 0; d < 8; ++d) if (seen[d]) std::cout << d << ' ';
        std::cout << "}\n";
        if ((int)c.size() != hs.ndof()) {
            std::cerr << "FAIL: bad coeff size\n"; return 1;
        }
    }

    std::cout << "ALL C++ TESTS PASSED\n";
    return 0;
}
