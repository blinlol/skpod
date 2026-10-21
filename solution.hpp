#pragma once

#include <vector>

class Solution {
public:
    Solution(double Lx, double Ly, double Lz, double T, int N, int K);

    static double frequency(double Lx, double Ly, double Lz);
    static double analytical(double x, double y, double z, double t,
                             double Lx, double Ly, double Lz);
    static double phi(double x, double y, double z,
                      double Lx, double Ly, double Lz);
    static double laplacian(double u_im, double u, double u_ip,
                            double u_jm, double u_jp,
                            double u_km, double u_kp,
                            double hx, double hy, double hz);
    static double layer1(double u0, double lap_phi, double tau);
    static double nextLayer(double u_prev, double u_curr, double lap, double tau);

    void run();
    double maxError() const;
    double seconds() const;

private:
    int index(int i, int j, int k) const;
    double pointLaplacian(const std::vector<double>& u, int i, int j, int k) const;
    void applyBoundaries(std::vector<double>& u) const;
    void fillLayer0(std::vector<double>& u) const;
    void fillLayer1(const std::vector<double>& u0, std::vector<double>& u1) const;
    void step(const std::vector<double>& prev,
              const std::vector<double>& curr,
              std::vector<double>& next) const;

    double Lx_;
    double Ly_;
    double Lz_;
    double T_;
    int N_;
    int K_;
    double hx_;
    double hy_;
    double hz_;
    double tau_;
    double seconds_;
    std::vector<double> result_;
};
