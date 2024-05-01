#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iomanip>


class Poly1d {
private:
    std::vector<double> coef;

public:
    Poly1d(std::vector<double> c) : coef(c) {}

    double operator()(double x) {
        double result = 0.0;
        for (int i = 0; i < coef.size(); i++) {
            result += coef[i] * std::pow(x, i);
        }
        return result;
    }

    Poly1d deriv() {
        std::vector<double> deriv_coef(coef.size() - 1);
        for (int i = 1; i < coef.size(); i++) {
            deriv_coef[i - 1] = coef[i] * i;
        }
        return Poly1d(deriv_coef);
    }
};

std::vector<double> calculate_a0(std::vector<double> coef, std::vector<double> x, std::vector<double> y) {
    Poly1d P(coef);
    Poly1d dP = P.deriv();
    int n = x.size();
    std::vector<double> a0(n);
    for (int i = 0; i < n - 1; i++) {
        a0[i] = (-dP(x[i + 1]) * std::pow(x[i], 2) * (x[i + 1]) * (x[i + 1] - x[i]) + y[i + 1] * std::pow(x[i], 2) * (3 * x[i + 1] - x[i])) / std::pow(x[i + 1] - x[i], 3) + (y[i] * std::pow(x[i + 1], 2) * (x[i + 1] - 3 * x[i]) - dP(x[i]) * x[i] * std::pow(x[i + 1], 2) * (x[i + 1] - x[i])) / std::pow(x[i + 1] - x[i], 3);
    }
    return a0;
}

std::vector<double> calculate_a1(std::vector<double> coef, std::vector<double> x, std::vector<double> y) {
    Poly1d P(coef);
    Poly1d dP = P.deriv();
    int n = x.size();
    std::vector<double> a1(n);
    for (int i = 0; i < n - 1; i++) {
        a1[i] = (dP(x[i + 1]) * x[i] * (x[i + 1] - x[i]) * (2 * x[i + 1] + x[i]) - 6 * (y[i + 1] - y[i]) * x[i] * x[i + 1]) / std::pow(x[i + 1] - x[i], 3) + (dP(x[i]) * x[i + 1] * (x[i + 1] - x[i]) * (x[i + 1] + 2 * x[i])) / std::pow(x[i + 1] - x[i], 3);
    }
    return a1;
}

std::vector<double> calculate_a2(std::vector<double> coef, std::vector<double> x, std::vector<double> y) {
    Poly1d P(coef);
    Poly1d dP = P.deriv();
    int n = x.size();
    std::vector<double> a2(n);
    for (int i = 0; i < n - 1; i++) {
        a2[i] = (-dP(x[i + 1]) * (x[i + 1] - x[i]) * (x[i + 1] + 2 * x[i]) + 3 * (y[i + 1] - y[i]) * (x[i + 1] + x[i])) / std::pow(x[i + 1] - x[i], 3) - (dP(x[i]) * (x[i + 1] - x[i]) * (2 * x[i + 1] + x[i])) / std::pow(x[i + 1] - x[i], 3);
    }
    return a2;
}

std::vector<double> calculate_a3(std::vector<double> coef, std::vector<double> x, std::vector<double> y) {
    Poly1d P(coef);
    Poly1d dP = P.deriv();
    int n = x.size();
    std::vector<double> a3(n);
    for (int i = 0; i < n - 1; i++) {
        a3[i] = (dP(x[i + 1]) * (x[i + 1] - x[i]) - 2 * (y[i + 1] - y[i]) + dP(x[i]) * (x[i + 1] - x[i])) / std::pow(x[i + 1] - x[i], 3);
    }
    return a3;
}
double calculate_spline(int i, const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c, const std::vector<double>& d, double x) {
    return a[i] + b[i] * x + (c[i]) * std::pow(x, 2) + (d[i]) * std::pow(x, 3);
}

double extrapolate(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c, const std::vector<double>& d, const std::vector<double>& x, double point) {
    int n = a.size();
    for (int i = 0; i < n; i++) {
        if (point <= x[i]) {
            return calculate_spline(i, a, b, c, d, point);
        }
    }
    return calculate_spline(n - 1, a, b, c, d, point - x[n - 1]);
}

int main() {
    std::vector<double> x;
    std::vector<double> y;
    int N;
    std::cin >> N;
    for (int i = 0; i < N; i++) {
        double x_val, y_val;
        std::cin >> x_val >> y_val;
        x.push_back(x_val);
        y.push_back(y_val);
    }
    int n;
    std::cin >> n;
    std::vector<double> coefs(n);
    for (int i = 0; i < n; i++) {
        std::cin >> coefs[i];
    }

    double value;
    std::cin >> value;
    std::vector<double> a0, a1, a2, a3;
    a0 = calculate_a0(coefs,x,y);
    a1 = calculate_a1(coefs,x,y);
    a2 = calculate_a2(coefs,x,y);
    a3 = calculate_a3(coefs,x,y);
    for (int i = 0; i < a1.size() - 1; i++) {
        std::cout << std::setprecision(5) << std::scientific << a3[i] << " " << a2[i] << " " << a1[i] << " " << a0[i] << std::endl;
    }
    double solution = extrapolate(a0, a1, a2, a3, x, value);
    std::cout << std::setprecision(5) << std::scientific << value << " " << solution << std::endl;
    return 0;
}