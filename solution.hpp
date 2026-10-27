#pragma once

#include <mpi.h>

#include <vector>

class Solution {
public:
    Solution(double Lx, double Ly, double Lz, double T, int N, int K);
    ~Solution();

    Solution(const Solution&) = delete;
    Solution& operator=(const Solution&) = delete;

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
    int rank() const;

private:
    static void splitRange(int length, int parts, int part, int& begin, int& count);

    int index(int li, int lj, int lk) const;
    int wrapZ(int k) const;
    int localI(int i) const;
    int localJ(int j) const;
    int localK(int k) const;
    int faceCount(int axis) const;
    double pointLaplacian(const std::vector<double>& u, int li, int lj, int lk) const;
    void applyBoundaries(std::vector<double>& u) const;
    void fillLayer0(std::vector<double>& u) const;
    void fillLayer1(const std::vector<double>& u0, std::vector<double>& u1) const;
    void step(const std::vector<double>& prev,
              const std::vector<double>& curr,
              std::vector<double>& next) const;
    void packPlane(const std::vector<double>& u, int axis, int localFixed,
                   std::vector<double>& buf) const;
    void unpackPlane(std::vector<double>& u, int axis, int localFixed,
                     const std::vector<double>& buf) const;
    void exchangeAxis(std::vector<double>& u, int axis, int minus, int plus);
    void copyPeriodicZ(std::vector<double>& u) const;
    void exchange(std::vector<double>& u);

    double Lx_ = 0;
    double Ly_ = 0;
    double Lz_ = 0;
    double T_ = 0;
    int N_ = 0;
    int K_ = 0;
    double hx_ = 0;
    double hy_ = 0;
    double hz_ = 0;
    double tau_ = 0;
    double seconds_ = 0;
    int rank_ = 0;
    int dims_[3] = {1, 1, 1};
    int coords_[3] = {0, 0, 0};
    int i0_ = 0;
    int nx_ = 0;
    int j0_ = 0;
    int ny_ = 0;
    int k0_ = 0;
    int nz_ = 0;
    int i_halo_left_ = 0;
    int i_halo_right_ = 0;
    int j_halo_left_ = 0;
    int j_halo_right_ = 0;
    int k_halo_left_ = 1;
    int k_halo_right_ = 1;
    int nxl_ = 0;
    int nyl_ = 0;
    int nzl_ = 0;
    int xm_ = MPI_PROC_NULL;
    int xp_ = MPI_PROC_NULL;
    int ym_ = MPI_PROC_NULL;
    int yp_ = MPI_PROC_NULL;
    int zm_ = MPI_PROC_NULL;
    int zp_ = MPI_PROC_NULL;
    MPI_Comm cart_ = MPI_COMM_NULL;
    std::vector<double> result_;
    std::vector<double> send_;
    std::vector<double> recv_;
};
