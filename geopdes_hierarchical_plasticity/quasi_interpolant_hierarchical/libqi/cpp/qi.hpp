/*
 * qi.hpp — Header-only modern C++ wrapper around libqi.
 *
 * #include "qi.hpp" and link against libqi.{so,dylib,a}.
 *
 * Example:
 *
 *   #include "qi.hpp"
 *   #include <vector>
 *
 *   qi::HSpace hs({3, 3}, {{16, 16}});      // degree (3,3), one level 16x16
 *
 *   std::vector<double> x(N), y(N), f(N);
 *   // ... fill x, y, f ...
 *
 *   auto coeffs = qi::compute_coeffs(hs, x, y, f);   // std::vector<double>
 *
 *   std::vector<double> gx(M), gy(M);
 *   auto fhat = qi::eval_tp(hs, 0, coeffs, gx, gy);
 */

#ifndef QI_HPP
#define QI_HPP

#include "qi.h"

#include <stdexcept>
#include <string>
#include <vector>
#include <cstddef>
#include <cstdint>
#include <utility>

namespace qi {

inline std::string version() { return std::string(qi_version()); }

// ------------------------------------------------------------------------- //
class Error : public std::runtime_error {
public:
    explicit Error(qi_status_t st)
        : std::runtime_error(make_msg(st)), status(st) {}
    qi_status_t status;
private:
    static std::string make_msg(qi_status_t st) {
        switch (st) {
            case QI_OK:                    return "OK";
            case QI_ERR_INVALID_ARG:       return "QI_ERR_INVALID_ARG";
            case QI_ERR_OUT_OF_MEMORY:     return "QI_ERR_OUT_OF_MEMORY";
            case QI_ERR_NOT_ENOUGH_DATA:   return "QI_ERR_NOT_ENOUGH_DATA";
            case QI_ERR_NO_ADMISSIBLE_DEG: return "QI_ERR_NO_ADMISSIBLE_DEG";
            case QI_ERR_SINGULAR_MATRIX:   return "QI_ERR_SINGULAR_MATRIX";
        }
        return "qi_unknown_error_" + std::to_string((int)st);
    }
};

inline void check(qi_status_t st) { if (st != QI_OK) throw Error(st); }

// ------------------------------------------------------------------------- //
// Hierarchical-space description with RAII storage of the marshalled arrays.
// ------------------------------------------------------------------------- //
class HSpace {
public:
    using Pair = std::pair<int, int>;

    HSpace(Pair degree, std::vector<Pair> nel,
           std::vector<std::vector<int>> active_1based = {})
        : degree_x_(degree.first), degree_y_(degree.second),
          nlevels_(static_cast<int>(nel.size()))
    {
        if (nlevels_ < 1)
            throw std::invalid_argument("HSpace: at least one level needed");

        nel_x_.reserve(nlevels_);
        nel_y_.reserve(nlevels_);
        ndof_per_level_.reserve(nlevels_);
        for (auto& p : nel) {
            nel_x_.push_back(p.first);
            nel_y_.push_back(p.second);
            ndof_per_level_.push_back((p.first + degree_x_) *
                                      (p.second + degree_y_));
        }

        active_offset_.push_back(0);
        if (active_1based.empty()) {
            // every basis function active on every level
            for (int lv = 0; lv < nlevels_; ++lv) {
                int n = ndof_per_level_[lv];
                int base = active_offset_.back();
                for (int i = 0; i < n; ++i) active_indices_.push_back(i);
                active_offset_.push_back(base + n);
            }
        } else {
            if ((int)active_1based.size() != nlevels_)
                throw std::invalid_argument("HSpace: active.size() != nlevels");
            for (auto& a : active_1based) {
                int base = active_offset_.back();
                for (int idx : a) active_indices_.push_back(idx - 1);
                active_offset_.push_back(base + (int)a.size());
            }
        }
        ndof_ = active_offset_.back();
    }

    /* Per-level non-uniform breaks (optional). */
    void set_breaks(const std::vector<std::pair<std::vector<double>,
                                                std::vector<double>>>& breaks)
    {
        if ((int)breaks.size() != nlevels_)
            throw std::invalid_argument("set_breaks: size != nlevels");
        breaks_x_data_.clear(); breaks_y_data_.clear();
        breaks_x_offset_.clear(); breaks_y_offset_.clear();
        breaks_x_offset_.push_back(0); breaks_y_offset_.push_back(0);
        for (auto& b : breaks) {
            breaks_x_data_.insert(breaks_x_data_.end(),
                                  b.first.begin(), b.first.end());
            breaks_y_data_.insert(breaks_y_data_.end(),
                                  b.second.begin(), b.second.end());
            breaks_x_offset_.push_back((int)breaks_x_data_.size());
            breaks_y_offset_.push_back((int)breaks_y_data_.size());
        }
    }

    int ndof()      const { return ndof_; }
    int degree_x()  const { return degree_x_; }
    int degree_y()  const { return degree_y_; }
    int nlevels()   const { return nlevels_; }

    /* Build a qi_hspace_t pointing at our owned buffers. */
    qi_hspace_t to_c() const {
        qi_hspace_t H;
        H.degree_x         = degree_x_;
        H.degree_y         = degree_y_;
        H.nlevels          = nlevels_;
        H.nel_x            = nel_x_.data();
        H.nel_y            = nel_y_.data();
        H.breaks_x_data    = breaks_x_data_.empty() ? nullptr : breaks_x_data_.data();
        H.breaks_x_offset  = breaks_x_offset_.empty() ? nullptr : breaks_x_offset_.data();
        H.breaks_y_data    = breaks_y_data_.empty() ? nullptr : breaks_y_data_.data();
        H.breaks_y_offset  = breaks_y_offset_.empty() ? nullptr : breaks_y_offset_.data();
        H.ndof             = ndof_;
        H.ndof_per_level   = ndof_per_level_.data();
        H.active_offset    = active_offset_.data();
        H.active_indices   = active_indices_.data();
        return H;
    }

private:
    int degree_x_, degree_y_, nlevels_, ndof_ = 0;
    std::vector<int> nel_x_, nel_y_, ndof_per_level_;
    std::vector<int> active_offset_, active_indices_;
    std::vector<double> breaks_x_data_, breaks_y_data_;
    std::vector<int>    breaks_x_offset_, breaks_y_offset_;
};

// ------------------------------------------------------------------------- //
// Algorithm options (thin wrapper over qi_options_t).
// ------------------------------------------------------------------------- //
struct Options {
    int    max_degree      = 0;     // 0 = use d
    double sigma_threshold = 0.0;   // 0 = library default (1e2)
    int    n_threads       = 0;     // 0 = OpenMP default
    int    verbose         = 0;

    qi_options_t to_c() const {
        qi_options_t o;
        qi_options_defaults(&o);
        if (max_degree)      o.max_degree      = max_degree;
        if (sigma_threshold) o.sigma_threshold = sigma_threshold;
        if (n_threads)       o.n_threads       = n_threads;
        o.verbose = verbose;
        return o;
    }
};

// ------------------------------------------------------------------------- //
// Diagnostics returned by compute_coeffs (one entry per basis function).
// ------------------------------------------------------------------------- //
struct Diagnostics {
    std::vector<int>    deg_used;
    std::vector<int>    nlocal;
    std::vector<double> cond_estimate;
    std::vector<int>    level;
};

// ------------------------------------------------------------------------- //
// compute_coeffs — main entry point.  Coefficients are returned row-major
// (size = ndof * ncomp).  ncomp is inferred from f.size() / x.size().
// ------------------------------------------------------------------------- //
inline std::vector<double>
compute_coeffs(const HSpace& hspace,
               const std::vector<double>& x,
               const std::vector<double>& y,
               const std::vector<double>& f,
               const Options& opt = Options{},
               Diagnostics* diag_out = nullptr)
{
    if (x.size() != y.size())
        throw std::invalid_argument("compute_coeffs: x.size() != y.size()");
    if (f.size() % x.size() != 0)
        throw std::invalid_argument("compute_coeffs: f.size() not multiple of n");
    int npts  = (int)x.size();
    int ncomp = (int)(f.size() / npts);

    qi_hspace_t  H = hspace.to_c();
    qi_options_t O = opt.to_c();
    qi_data_t    D{ npts, ncomp, x.data(), y.data(), f.data() };

    std::vector<double> coeffs((std::size_t)H.ndof * ncomp);

    qi_diag_t diag{};
    Diagnostics tmp;
    if (diag_out) {
        tmp.deg_used.assign(H.ndof, 0);
        tmp.nlocal.assign(H.ndof, 0);
        tmp.cond_estimate.assign(H.ndof, 0.0);
        tmp.level.assign(H.ndof, 0);
        diag.deg_used      = tmp.deg_used.data();
        diag.nlocal        = tmp.nlocal.data();
        diag.cond_estimate = tmp.cond_estimate.data();
        diag.level         = tmp.level.data();
    }

    qi_status_t st = qi_compute_coeffs(&H, &D, &O, coeffs.data(),
                                       diag_out ? &diag : nullptr);
    check(st);

    if (diag_out) *diag_out = std::move(tmp);
    return coeffs;
}

// ------------------------------------------------------------------------- //
// eval_tp — evaluate a single-level tensor product.  out resized as needed.
// ------------------------------------------------------------------------- //
inline void
eval_tp(const HSpace& hspace, int level,
        const std::vector<double>& coeffs_lev, int ncomp,
        const std::vector<double>& x,
        const std::vector<double>& y,
        std::vector<double>& out)
{
    if (x.size() != y.size())
        throw std::invalid_argument("eval_tp: x.size() != y.size()");
    qi_hspace_t H = hspace.to_c();
    out.assign(x.size() * (std::size_t)ncomp, 0.0);
    qi_status_t st = qi_eval_tp(&H, level, coeffs_lev.data(),
                                ncomp, (int)x.size(),
                                x.data(), y.data(), out.data());
    check(st);
}

// Convenience scalar overload
inline std::vector<double>
eval_tp(const HSpace& hspace, int level,
        const std::vector<double>& coeffs_lev,
        const std::vector<double>& x,
        const std::vector<double>& y)
{
    std::vector<double> out;
    eval_tp(hspace, level, coeffs_lev, 1, x, y, out);
    return out;
}

// ------------------------------------------------------------------------- //
// to_finest — convert to fine-grid tensor-product coefficients.
// ------------------------------------------------------------------------- //
inline std::vector<double>
to_finest(const HSpace& hspace, const std::vector<double>& coeffs, int ncomp)
{
    qi_hspace_t H = hspace.to_c();
    int M = hspace.nlevels();
    int nelx = H.nel_x[M-1];
    int nely = H.nel_y[M-1];
    int nfin = (nelx + H.degree_x) * (nely + H.degree_y);
    std::vector<double> out((std::size_t)nfin * ncomp);
    check(qi_to_finest(&H, coeffs.data(), ncomp, out.data()));
    return out;
}

} // namespace qi

#endif // QI_HPP
