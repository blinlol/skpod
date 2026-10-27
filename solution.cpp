#include "solution.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <utility>

constexpr double kPi = 3.14159265358979323846;

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

    int size = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int dims[3] = {0, 0, 0};
    MPI_Dims_create(size, 3, dims);
    if (dims[0] > N + 1 || dims[1] > N + 1 || dims[2] > N) {
        throw std::invalid_argument("process grid is finer than the mesh");
    }

    int periods[3] = {0, 0, 1};
    MPI_Comm cart = MPI_COMM_NULL;
    if (MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, 0, &cart) != MPI_SUCCESS) {
        throw std::runtime_error("MPI_Cart_create failed");
    }
    cart_ = cart;
    dims_[0] = dims[0];
    dims_[1] = dims[1];
    dims_[2] = dims[2];
    MPI_Comm_rank(cart_, &rank_);
    MPI_Cart_coords(cart_, rank_, 3, coords_);
    MPI_Cart_shift(cart_, 0, 1, &xm_, &xp_);
    MPI_Cart_shift(cart_, 1, 1, &ym_, &yp_);
    MPI_Cart_shift(cart_, 2, 1, &zm_, &zp_);

    splitRange(N + 1, dims_[0], coords_[0], i0_, nx_);
    splitRange(N + 1, dims_[1], coords_[1], j0_, ny_);
    splitRange(N, dims_[2], coords_[2], k0_, nz_);

    i_halo_left_ = i0_ > 0 ? 1 : 0;
    i_halo_right_ = i0_ + nx_ < N + 1 ? 1 : 0;
    j_halo_left_ = j0_ > 0 ? 1 : 0;
    j_halo_right_ = j0_ + ny_ < N + 1 ? 1 : 0;
    k_halo_left_ = 1;
    k_halo_right_ = 1;
    nxl_ = nx_ + i_halo_left_ + i_halo_right_;
    nyl_ = ny_ + j_halo_left_ + j_halo_right_;
    nzl_ = nz_ + k_halo_left_ + k_halo_right_;
}

Solution::~Solution() {
    if (cart_ != MPI_COMM_NULL) {
        MPI_Comm_free(&cart_);
    }
}

double Solution::frequency(double Lx, double Ly, double Lz) {
    return kPi * std::sqrt(1.0 / (Lx * Lx) + 1.0 / (Ly * Ly) + 4.0 / (Lz * Lz));
}

double Solution::analytical(double x, double y, double z, double t,
                            double Lx, double Ly, double Lz) {
    double at = frequency(Lx, Ly, Lz);
    return std::sin(kPi * x / Lx) * std::sin(kPi * y / Ly) * std::sin(2.0 * kPi * z / Lz)
           * std::cos(at * t + 2.0 * kPi);
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

void Solution::splitRange(int length, int parts, int part, int& begin, int& count) {
    int q = length / parts;
    int r = length % parts;
    count = q + (part < r ? 1 : 0);
    begin = part * q + std::min(part, r);
}

int Solution::index(int li, int lj, int lk) const {
    return (li * nyl_ + lj) * nzl_ + lk;
}

int Solution::wrapZ(int k) const {
    if (k < 0) {
        return k + N_;
    }
    if (k >= N_) {
        return k - N_;
    }
    return k;
}

int Solution::localI(int i) const {
    return i - i0_ + i_halo_left_;
}

int Solution::localJ(int j) const {
    return j - j0_ + j_halo_left_;
}

int Solution::localK(int k) const {
    return k - k0_ + k_halo_left_;
}

int Solution::faceCount(int axis) const {
    if (axis == 0) {
        return ny_ * nz_;
    }
    if (axis == 1) {
        return nx_ * nz_;
    }
    return nx_ * ny_;
}

double Solution::pointLaplacian(const std::vector<double>& u, int li, int lj, int lk) const {
    double center = u[index(li, lj, lk)];
    return laplacian(u[index(li - 1, lj, lk)], center, u[index(li + 1, lj, lk)],
                     u[index(li, lj - 1, lk)], u[index(li, lj + 1, lk)],
                     u[index(li, lj, lk - 1)], u[index(li, lj, lk + 1)],
                     hx_, hy_, hz_);
}

void Solution::applyBoundaries(std::vector<double>& u) const {
    if (i0_ == 0 || i0_ + nx_ == N_ + 1) {
#pragma omp parallel for schedule(static)
        for (int j = 0; j < ny_; ++j) {
            int lj = localJ(j0_ + j);
            for (int k = 0; k < nz_; ++k) {
                int lk = localK(k0_ + k);
                if (i0_ == 0) {
                    u[index(localI(0), lj, lk)] = 0.0;
                }
                if (i0_ + nx_ == N_ + 1) {
                    u[index(localI(N_), lj, lk)] = 0.0;
                }
            }
        }
    }
    if (j0_ == 0 || j0_ + ny_ == N_ + 1) {
#pragma omp parallel for schedule(static)
        for (int i = 0; i < nx_; ++i) {
            int li = localI(i0_ + i);
            for (int k = 0; k < nz_; ++k) {
                int lk = localK(k0_ + k);
                if (j0_ == 0) {
                    u[index(li, localJ(0), lk)] = 0.0;
                }
                if (j0_ + ny_ == N_ + 1) {
                    u[index(li, localJ(N_), lk)] = 0.0;
                }
            }
        }
    }
}

void Solution::fillLayer0(std::vector<double>& u) const {
#pragma omp parallel for schedule(static)
    for (int li = 0; li < nxl_; ++li) {
        int gi = i0_ + li - i_halo_left_;
        double x = gi * hx_;
        for (int lj = 0; lj < nyl_; ++lj) {
            int gj = j0_ + lj - j_halo_left_;
            double y = gj * hy_;
            for (int lk = 0; lk < nzl_; ++lk) {
                int gk = wrapZ(k0_ + lk - k_halo_left_);
                double z = gk * hz_;
                u[index(li, lj, lk)] = phi(x, y, z, Lx_, Ly_, Lz_);
            }
        }
    }
    applyBoundaries(u);
}

void Solution::fillLayer1(const std::vector<double>& u0, std::vector<double>& u1) const {
    int iBegin = std::max(i0_, 1);
    int iEnd = std::min(i0_ + nx_, N_);
    int jBegin = std::max(j0_, 1);
    int jEnd = std::min(j0_ + ny_, N_);
#pragma omp parallel for schedule(static)
    for (int gi = iBegin; gi < iEnd; ++gi) {
        int li = localI(gi);
        for (int gj = jBegin; gj < jEnd; ++gj) {
            int lj = localJ(gj);
            for (int gk = k0_; gk < k0_ + nz_; ++gk) {
                int lk = localK(gk);
                int id = index(li, lj, lk);
                double lap = pointLaplacian(u0, li, lj, lk);
                u1[id] = layer1(u0[id], lap, tau_);
            }
        }
    }
    applyBoundaries(u1);
}

void Solution::step(const std::vector<double>& prev,
                    const std::vector<double>& curr,
                    std::vector<double>& next) const {
    int iBegin = std::max(i0_, 1);
    int iEnd = std::min(i0_ + nx_, N_);
    int jBegin = std::max(j0_, 1);
    int jEnd = std::min(j0_ + ny_, N_);
#pragma omp parallel for schedule(static)
    for (int gi = iBegin; gi < iEnd; ++gi) {
        int li = localI(gi);
        for (int gj = jBegin; gj < jEnd; ++gj) {
            int lj = localJ(gj);
            for (int gk = k0_; gk < k0_ + nz_; ++gk) {
                int lk = localK(gk);
                int id = index(li, lj, lk);
                double lap = pointLaplacian(curr, li, lj, lk);
                next[id] = nextLayer(prev[id], curr[id], lap, tau_);
            }
        }
    }
    applyBoundaries(next);
}

void Solution::packPlane(const std::vector<double>& u, int axis, int localFixed,
                         std::vector<double>& buf) const {
    int p = 0;
    if (axis == 0) {
        for (int j = 0; j < ny_; ++j) {
            int lj = localJ(j0_ + j);
            for (int k = 0; k < nz_; ++k) {
                buf[p++] = u[index(localFixed, lj, localK(k0_ + k))];
            }
        }
    } else if (axis == 1) {
        for (int i = 0; i < nx_; ++i) {
            int li = localI(i0_ + i);
            for (int k = 0; k < nz_; ++k) {
                buf[p++] = u[index(li, localFixed, localK(k0_ + k))];
            }
        }
    } else {
        for (int i = 0; i < nx_; ++i) {
            int li = localI(i0_ + i);
            for (int j = 0; j < ny_; ++j) {
                buf[p++] = u[index(li, localJ(j0_ + j), localFixed)];
            }
        }
    }
}

void Solution::unpackPlane(std::vector<double>& u, int axis, int localFixed,
                           const std::vector<double>& buf) const {
    int p = 0;
    if (axis == 0) {
        for (int j = 0; j < ny_; ++j) {
            int lj = localJ(j0_ + j);
            for (int k = 0; k < nz_; ++k) {
                u[index(localFixed, lj, localK(k0_ + k))] = buf[p++];
            }
        }
    } else if (axis == 1) {
        for (int i = 0; i < nx_; ++i) {
            int li = localI(i0_ + i);
            for (int k = 0; k < nz_; ++k) {
                u[index(li, localFixed, localK(k0_ + k))] = buf[p++];
            }
        }
    } else {
        for (int i = 0; i < nx_; ++i) {
            int li = localI(i0_ + i);
            for (int j = 0; j < ny_; ++j) {
                u[index(li, localJ(j0_ + j), localFixed)] = buf[p++];
            }
        }
    }
}

void Solution::exchangeAxis(std::vector<double>& u, int axis, int minus, int plus) {
    if (minus == MPI_PROC_NULL && plus == MPI_PROC_NULL) {
        return;
    }

    int count = faceCount(axis);
    int ownedFirst = 0;
    int ownedLast = 0;
    int haloMinus = 0;
    int haloPlus = 0;
    if (axis == 0) {
        ownedFirst = localI(i0_);
        ownedLast = localI(i0_ + nx_ - 1);
        haloMinus = localI(i0_ - 1);
        haloPlus = localI(i0_ + nx_);
    } else if (axis == 1) {
        ownedFirst = localJ(j0_);
        ownedLast = localJ(j0_ + ny_ - 1);
        haloMinus = localJ(j0_ - 1);
        haloPlus = localJ(j0_ + ny_);
    } else {
        ownedFirst = localK(k0_);
        ownedLast = localK(k0_ + nz_ - 1);
        haloMinus = localK(k0_ - 1);
        haloPlus = localK(k0_ + nz_);
    }

    int sendPlus = plus == MPI_PROC_NULL ? 0 : count;
    int recvMinus = minus == MPI_PROC_NULL ? 0 : count;
    if (sendPlus > 0) {
        packPlane(u, axis, ownedLast, send_);
    }
    MPI_Sendrecv(send_.data(), sendPlus, MPI_DOUBLE, plus, 0,
                 recv_.data(), recvMinus, MPI_DOUBLE, minus, 0,
                 cart_, MPI_STATUS_IGNORE);
    if (recvMinus > 0) {
        unpackPlane(u, axis, haloMinus, recv_);
    }

    int sendMinus = minus == MPI_PROC_NULL ? 0 : count;
    int recvPlus = plus == MPI_PROC_NULL ? 0 : count;
    if (sendMinus > 0) {
        packPlane(u, axis, ownedFirst, send_);
    }
    MPI_Sendrecv(send_.data(), sendMinus, MPI_DOUBLE, minus, 0,
                 recv_.data(), recvPlus, MPI_DOUBLE, plus, 0,
                 cart_, MPI_STATUS_IGNORE);
    if (recvPlus > 0) {
        unpackPlane(u, axis, haloPlus, recv_);
    }
}

void Solution::copyPeriodicZ(std::vector<double>& u) const {
    int srcLeft = localK(k0_ + nz_ - 1);
    int srcRight = localK(k0_);
    int dstLeft = localK(k0_ - 1);
    int dstRight = localK(k0_ + nz_);
    for (int i = 0; i < nx_; ++i) {
        int li = localI(i0_ + i);
        for (int j = 0; j < ny_; ++j) {
            int lj = localJ(j0_ + j);
            u[index(li, lj, dstLeft)] = u[index(li, lj, srcLeft)];
            u[index(li, lj, dstRight)] = u[index(li, lj, srcRight)];
        }
    }
}

void Solution::exchange(std::vector<double>& u) {
    exchangeAxis(u, 0, xm_, xp_);
    exchangeAxis(u, 1, ym_, yp_);
    if (dims_[2] == 1) {
        copyPeriodicZ(u);
    } else {
        exchangeAxis(u, 2, zm_, zp_);
    }
}

void Solution::run() {
    double inv_h2 = 1.0 / (hx_ * hx_) + 1.0 / (hy_ * hy_) + 1.0 / (hz_ * hz_);
    double cfl_limit = 1.0 / std::sqrt(inv_h2);
    if (rank_ == 0 && tau_ > cfl_limit) {
        std::cerr << "warning: scheme is unstable, tau " << tau_
                  << " exceeds CFL limit " << cfl_limit << '\n';
    }

    std::size_t cells = static_cast<std::size_t>(nxl_) * static_cast<std::size_t>(nyl_) * static_cast<std::size_t>(nzl_);
    std::vector<double> prev(cells, 0.0);
    std::vector<double> curr(cells, 0.0);
    std::vector<double> next(cells, 0.0);
    int face = std::max(faceCount(0), std::max(faceCount(1), faceCount(2)));
    send_.assign(static_cast<std::size_t>(face), 0.0);
    recv_.assign(static_cast<std::size_t>(face), 0.0);

    MPI_Barrier(cart_);
    double started = MPI_Wtime();
    fillLayer0(prev);
    fillLayer1(prev, curr);
    for (int n = 1; n < K_; ++n) {
        exchange(curr);
        step(prev, curr, next);
        std::swap(prev, curr);
        std::swap(curr, next);
    }
    MPI_Barrier(cart_);
    seconds_ = MPI_Wtime() - started;
    result_ = std::move(curr);
}

double Solution::maxError() const {
    if (result_.empty()) {
        throw std::logic_error("run() must be called before maxError()");
    }

    double delta = 0.0;
#pragma omp parallel for schedule(static) reduction(max:delta)
    for (int gi = i0_; gi < i0_ + nx_; ++gi) {
        double x = gi * hx_;
        for (int gj = j0_; gj < j0_ + ny_; ++gj) {
            double y = gj * hy_;
            for (int gk = k0_; gk < k0_ + nz_; ++gk) {
                double z = gk * hz_;
                double exact = analytical(x, y, z, T_, Lx_, Ly_, Lz_);
                int id = index(localI(gi), localJ(gj), localK(gk));
                double diff = std::abs(result_[id] - exact);
                if (diff > delta) {
                    delta = diff;
                }
            }
        }
    }

    double global = 0.0;
    MPI_Allreduce(&delta, &global, 1, MPI_DOUBLE, MPI_MAX, cart_);
    return global;
}

double Solution::seconds() const {
    return seconds_;
}

int Solution::rank() const {
    return rank_;
}
