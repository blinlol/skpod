#include "solution.hpp"

#include <chrono>
#include <cmath>
#include <iostream>
#include <numbers>
#include <stdexcept>
#include <utility>

Solution::Solution(double Lx, double Ly, double Lz, double T, int N, int K) {
    if (Lx <= 0 || Ly <= 0 || Lz <= 0 || T <= 0) {
        throw std::invalid_argument("Lx, Ly, Lz and T must be positive");
    }
    if (N < 2) {
        throw std::invalid_argument("N must be at least 2");
    }
    if (K < 1) {
        throw std::invalid_argument("K must be at least 1");
    }

    Lx_ = Lx;
    Ly_ = Ly;
    Lz_ = Lz;
    T_ = T;
    N_ = N;
    K_ = K;
    hx_ = Lx / N;
    hy_ = Ly / N;
    hz_ = Lz / N;
    tau_ = T / K;
    seconds_ = 0;
}

double Solution::frequency(double Lx, double Ly, double Lz) {
    double pi = std::numbers::pi;
    return pi * std::sqrt(1.0 / (Lx * Lx) + 1.0 / (Ly * Ly) + 4.0 / (Lz * Lz));
}

double Solution::analytical(double x, double y, double z, double t,
                            double Lx, double Ly, double Lz) {
    double pi = std::numbers::pi;
    double at = frequency(Lx, Ly, Lz);
    return std::sin(pi * x / Lx) * std::sin(pi * y / Ly) * std::sin(2.0 * pi * z / Lz)
           * std::cos(at * t + 2.0 * pi);
}

double Solution::phi(double x, double y, double z, double Lx, double Ly, double Lz) {
    return analytical(x, y, z, 0.0, Lx, Ly, Lz);
}

double Solution::laplacian(double u_im, double u, double u_ip,
                           double u_jm, double u_jp,
                           double u_km, double u_kp,
                           double hx, double hy, double hz) {
    double dx = (u_im - 2.0 * u + u_ip) / (hx * hx);
    double dy = (u_jm - 2.0 * u + u_jp) / (hy * hy);
    double dz = (u_km - 2.0 * u + u_kp) / (hz * hz);
    return dx + dy + dz;
}

double Solution::layer1(double u0, double lap_phi, double tau) {
    return u0 + tau * tau / 2.0 * lap_phi;
}

double Solution::nextLayer(double u_prev, double u_curr, double lap, double tau) {
    return 2.0 * u_curr - u_prev + tau * tau * lap;
}

int Solution::index(int i, int j, int k) const {
    int nodes = N_ + 1;
    return (i * nodes + j) * nodes + k;
}

double Solution::pointLaplacian(const std::vector<double>& u, int i, int j, int k) const {
    int k_left = k - 1;
    int k_right = k + 1;
    if (k == 0) {
        k_left = N_ - 1;
    }
    if (k == N_ - 1) {
        k_right = 0;
    }

    double center = u[index(i, j, k)];
    return laplacian(u[index(i - 1, j, k)], center, u[index(i + 1, j, k)],
                     u[index(i, j - 1, k)], u[index(i, j + 1, k)],
                     u[index(i, j, k_left)], u[index(i, j, k_right)],
                     hx_, hy_, hz_);
}

void Solution::applyBoundaries(std::vector<double>& u) const {
    int nodes = N_ + 1;
    for (int j = 0; j < nodes; ++j) {
        for (int k = 0; k < nodes; ++k) {
            u[index(0, j, k)] = 0.0;
            u[index(N_, j, k)] = 0.0;
        }
    }
    for (int i = 0; i < nodes; ++i) {
        for (int k = 0; k < nodes; ++k) {
            u[index(i, 0, k)] = 0.0;
            u[index(i, N_, k)] = 0.0;
        }
    }
    for (int i = 0; i < nodes; ++i) {
        for (int j = 0; j < nodes; ++j) {
            u[index(i, j, N_)] = u[index(i, j, 0)];
        }
    }
}

void Solution::fillLayer0(std::vector<double>& u) const {
    int nodes = N_ + 1;
    for (int i = 0; i < nodes; ++i) {
        double x = i * hx_;
        for (int j = 0; j < nodes; ++j) {
            double y = j * hy_;
            for (int k = 0; k < nodes; ++k) {
                double z = k * hz_;
                u[index(i, j, k)] = phi(x, y, z, Lx_, Ly_, Lz_);
            }
        }
    }
    applyBoundaries(u);
}

void Solution::fillLayer1(const std::vector<double>& u0, std::vector<double>& u1) const {
    for (int i = 1; i < N_; ++i) {
        for (int j = 1; j < N_; ++j) {
            for (int k = 0; k < N_; ++k) {
                double value = u0[index(i, j, k)];
                double lap = pointLaplacian(u0, i, j, k);
                u1[index(i, j, k)] = layer1(value, lap, tau_);
            }
        }
    }
    applyBoundaries(u1);
}

void Solution::step(const std::vector<double>& prev,
                    const std::vector<double>& curr,
                    std::vector<double>& next) const {
    for (int i = 1; i < N_; ++i) {
        for (int j = 1; j < N_; ++j) {
            for (int k = 0; k < N_; ++k) {
                int id = index(i, j, k);
                double lap = pointLaplacian(curr, i, j, k);
                next[id] = nextLayer(prev[id], curr[id], lap, tau_);
            }
        }
    }
    applyBoundaries(next);
}

void Solution::run() {
    double inv_h2 = 1.0 / (hx_ * hx_) + 1.0 / (hy_ * hy_) + 1.0 / (hz_ * hz_);
    double cfl_limit = 1.0 / std::sqrt(inv_h2);
    if (tau_ > cfl_limit) {
        std::cerr << "warning: scheme is unstable, tau " << tau_
                  << " exceeds CFL limit " << cfl_limit << '\n';
    }

    int nodes = N_ + 1;
    int cells = nodes * nodes * nodes;
    std::vector<double> prev(cells, 0.0);
    std::vector<double> curr(cells, 0.0);
    std::vector<double> next(cells, 0.0);

    auto start = std::chrono::steady_clock::now();
    fillLayer0(prev);
    fillLayer1(prev, curr);
    for (int n = 1; n < K_; ++n) {
        step(prev, curr, next);
        std::swap(prev, curr);
        std::swap(curr, next);
    }
    auto finish = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    seconds_ = elapsed.count();
    result_ = curr;
}

double Solution::maxError() const {
    if (result_.empty()) {
        throw std::logic_error("run() must be called before maxError()");
    }

    double delta = 0.0;
    int nodes = N_ + 1;
    for (int i = 0; i < nodes; ++i) {
        double x = i * hx_;
        for (int j = 0; j < nodes; ++j) {
            double y = j * hy_;
            for (int k = 0; k < nodes; ++k) {
                double z = k * hz_;
                double exact = analytical(x, y, z, T_, Lx_, Ly_, Lz_);
                double diff = std::abs(result_[index(i, j, k)] - exact);
                if (diff > delta) {
                    delta = diff;
                }
            }
        }
    }
    return delta;
}

double Solution::seconds() const {
    return seconds_;
}
