#include <iostream>
#include <vector>
#include <cmath>
#include <functional>

std::vector<double> Hermite(std::vector<double> x, std::vector<double> y, std::vector<double> der1, std::vector<double> der2, std::vector<int> der1_coord, std::vector<int> der2_coord) {
    int n = x.size() + der1.size() + der2.size();
    std::vector<double> a(n);
    std::vector<double> H;
    std::vector<double> derH;
    std::vector<double> der2H;
    double x0;
    double h = a[0];
    for (int i = 1; i < n; i++) {
        h += a[i] * pow(x0, i);
    }
    std::function<double(double)> H = [&](double x0) { return h; };
    std::function<double(double)> dH = [&](double x0) { return h.diff(x0); };
    std::function<double(double)> d2H = [&](double x0) { return dH.diff(x0); };
    std::vector<double> eq;
    return eq;
}

int main() {
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> der1;
    std::vector<double> der2;
    int N;
    std::cin >> N;
    for (int i = 0; i < N; i++) {
        double x0;
        std::cin >> x0;
        x.push_back(x0);
    }
    for (int i = 0; i < N; i++) {
        double y0;
        std::cin >> y0;
        y.push_back(y0);
    }
    std::vector<int> der1_coord;
    std::vector<int> der2_coord;
    int N1;
    std::cin >> N1;
    for (int i = 0; i < N1; i++) {
        int point;
        std::cin >> point;
        double der;
        std::cin >> der;
        der1_coord.push_back(point - 1);
        der1.push_back(der);
    }
    int N2;
    std::cin >> N2;
    for (int i = 0; i < N2; i++) {
        int point;
        std::cin >> point;
        double der;
        std::cin >> der;
        der2_coord.push_back(point - 1);
        der2.push_back(der);
    }
    std::vector<double> result = Hermite(x, y, der1, der2, der1_coord, der2_coord);
    for (int i = 0; i < result.size(); i++) {
        std::cout << result[i] << std::endl;
    }
    return 0;
}