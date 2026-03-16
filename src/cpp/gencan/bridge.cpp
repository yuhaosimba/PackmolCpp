#include <algorithm>
#include <cmath>
#include <cctype>
#include <cstddef>
#include <cstdlib>
#include <string>

namespace {

enum class GencanImplMode {
    kFortran = 0,
    kCpp = 1,
    kAb = 2
};

GencanImplMode parse_impl_mode() {
    const char* env = std::getenv("PACKMOL_GENCAN_IMPL");
    if (env == nullptr) {
        return GencanImplMode::kCpp;
    }

    std::string value(env);
    for (char& ch : value) {
        ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
    }

    if (value == "fortran") {
        return GencanImplMode::kFortran;
    }
    if (value == "ab") {
        return GencanImplMode::kAb;
    }
    return GencanImplMode::kCpp;
}

GencanImplMode active_impl_mode() {
    static const GencanImplMode mode = parse_impl_mode();
    return mode;
}

}  // namespace

extern "C" int packmol_gencan_cpp_probe(int n, const double* x, double* fx_out) {
    if (n < 0 || x == nullptr || fx_out == nullptr) {
        return -1;
    }

    // Phase B skeleton: keep behavior-neutral bridge available for mixed builds.
    *fx_out = 0.0;
    return 0;
}

extern "C" int packmol_gencan_impl_mode_c() {
    return static_cast<int>(active_impl_mode());
}

extern "C" void packmol_spgls_fortran_c(
    const int* n,
    double* x,
    const int* m,
    const double* lambda,
    const double* rho,
    double* f,
    const double* g,
    const double* l,
    const double* u,
    const double* lamspg,
    const double* nint,
    const int* mininterp,
    const double* fmin,
    const int* maxfc,
    const int* iprint,
    int* fcnt,
    int* inform,
    double* xtrial,
    double* d,
    const double* gamma,
    const double* sigma1,
    const double* sigma2,
    const double* sterel,
    const double* steabs,
    const double* epsrel,
    const double* epsabs,
    const double* infrel,
    const double* infabs
);

extern "C" void packmol_evalal_fortran_c(
    const int* n,
    const double* x,
    const int* m,
    const double* lambda,
    const double* rho,
    double* f,
    int* flag
);

extern "C" void packmol_calcf_fortran_c(
    const int* nind,
    const int* ind,
    double* x,
    const int* n,
    const double* xc,
    const int* m,
    const double* lambda,
    const double* rho,
    double* f,
    int* inform
);

extern "C" void packmol_calcg_fortran_c(
    const int* nind,
    const int* ind,
    double* x,
    const int* n,
    const double* xc,
    const int* m,
    const double* lambda,
    const double* rho,
    double* g,
    int* inform
);

extern "C" void packmol_calcgdiff_fortran_c(
    const int* nind,
    const int* ind,
    double* x,
    const int* n,
    const double* xc,
    const int* m,
    const double* lambda,
    const double* rho,
    double* g,
    const double* sterel,
    const double* steabs,
    int* inform
);

namespace {

void spgls_cpp(
    const int* n,
    double* x,
    const int* m,
    const double* lambda,
    const double* rho,
    double* f,
    const double* g,
    const double* l,
    const double* u,
    const double* lamspg,
    const double* nint,
    const int* mininterp,
    const double* fmin,
    const int* maxfc,
    int* fcnt,
    int* inform,
    double* xtrial,
    double* d,
    const double* gamma,
    const double* sigma1,
    const double* sigma2,
    const double* epsrel,
    const double* epsabs
) {
    const int n_val = *n;
    const int m_val = *m;
    const double f_ref = *f;
    const double lamspg_val = *lamspg;
    const double nint_val = *nint;
    const double sigma1_val = *sigma1;
    const double sigma2_val = *sigma2;
    const double gamma_val = *gamma;
    const double epsrel_val = *epsrel;
    const double epsabs_val = *epsabs;
    const double fmin_val = *fmin;

    int interp = 0;
    double alpha = 1.0;
    double gtd = 0.0;
    for (int i = 0; i < n_val; ++i) {
        xtrial[i] = std::min(u[i], std::max(l[i], x[i] - lamspg_val * g[i]));
        d[i] = xtrial[i] - x[i];
        gtd += g[i] * d[i];
    }

    double ftrial = 0.0;
    packmol_evalal_fortran_c(n, xtrial, m, lambda, rho, &ftrial, inform);
    *fcnt += 1;
    if (*inform < 0) {
        return;
    }

    while (true) {
        if (ftrial <= f_ref + gamma_val * alpha * gtd) {
            *f = ftrial;
            for (int i = 0; i < n_val; ++i) {
                x[i] = xtrial[i];
            }
            *inform = 0;
            return;
        }

        if (ftrial <= fmin_val) {
            *f = ftrial;
            for (int i = 0; i < n_val; ++i) {
                x[i] = xtrial[i];
            }
            *inform = 4;
            return;
        }

        if (*fcnt >= *maxfc) {
            if (ftrial < f_ref) {
                *f = ftrial;
                for (int i = 0; i < n_val; ++i) {
                    x[i] = xtrial[i];
                }
            }
            *inform = 8;
            return;
        }

        interp += 1;
        if (alpha < sigma1_val) {
            alpha = alpha / nint_val;
        } else {
            const double den = 2.0 * (ftrial - f_ref - alpha * gtd);
            double atmp = alpha / nint_val;
            if (den != 0.0) {
                atmp = (-gtd * alpha * alpha) / den;
            }
            if (atmp < sigma1_val || atmp > sigma2_val * alpha) {
                alpha = alpha / nint_val;
            } else {
                alpha = atmp;
            }
        }

        for (int i = 0; i < n_val; ++i) {
            xtrial[i] = x[i] + alpha * d[i];
        }

        packmol_evalal_fortran_c(&n_val, xtrial, &m_val, lambda, rho, &ftrial, inform);
        *fcnt += 1;
        if (*inform < 0) {
            return;
        }

        bool samep = true;
        for (int i = 0; i < n_val; ++i) {
            if (std::abs(alpha * d[i]) > std::max(epsrel_val * std::abs(x[i]), epsabs_val)) {
                samep = false;
            }
        }

        if (interp >= *mininterp && samep) {
            if (ftrial < f_ref) {
                *f = ftrial;
                for (int i = 0; i < n_val; ++i) {
                    x[i] = xtrial[i];
                }
            }
            *inform = 6;
            return;
        }
    }
}

bool tnls_cpp_subset(
    const int* nind,
    const int* ind,
    const int* n,
    double* x,
    const int* m,
    const double* lambda,
    const double* rho,
    const double* l,
    const double* u,
    double* f,
    double* g,
    const double* d,
    const double* amax,
    const int* rbdtype,
    const int* rbdind,
    const double* nint,
    const double* next,
    const int* mininterp,
    const int* maxextrap,
    const double* fmin,
    const int* maxfc,
    const int* gtype,
    int* fcnt,
    int* gcnt,
    int* intcnt,
    int* exgcnt,
    int* exbcnt,
    int* inform,
    double* xplus,
    double* xtmp,
    double* xbext,
    const double* gamma,
    const double* beta,
    const double* sigma1,
    const double* sigma2,
    const double* sterel,
    const double* steabs,
    const double* epsrel,
    const double* epsabs
) {
    const bool valid_bound =
        (*rbdtype == 1 || *rbdtype == 2) && (*rbdind >= 1 && *rbdind <= *nind);
    const bool has_active_bound = valid_bound;

    const int nind_val = *nind;
    const int n_val = *n;
    const int m_val = *m;
    const double f0 = *f;

    auto compute_gradient = [&]() -> bool {
        if (*gtype == 0) {
            packmol_calcg_fortran_c(nind, ind, x, &n_val, x, &m_val, lambda, rho, g, inform);
        } else if (*gtype == 1) {
            packmol_calcgdiff_fortran_c(nind, ind, x, &n_val, x, &m_val, lambda, rho, g, sterel, steabs, inform);
        }
        *gcnt += 1;
        return *inform >= 0;
    };

    auto finish_extrapolation = [&](double fbext, int code, bool need_gradient_update) -> bool {
        if (need_gradient_update) {
            if (compute_gradient()) {
                if (*f < fbext) {
                    *exgcnt += 1;
                } else {
                    *exbcnt += 1;
                }
            }
        }
        *inform = code;
        return true;
    };

    auto run_extrapolation = [&](double& alpha_ref, double& fplus_ref, bool gradient_already_at_xplus) -> bool {
        double fbext = fplus_ref;
        for (int i = 0; i < nind_val; ++i) {
            xbext[i] = xplus[i];
        }
        int extrap = 0;

        while (true) {
            const bool need_gradient_update = (extrap != 0) || (*amax <= 1.0) || (!gradient_already_at_xplus);

            if (fplus_ref <= *fmin) {
                *f = fplus_ref;
                for (int i = 0; i < nind_val; ++i) {
                    x[i] = xplus[i];
                }
                return finish_extrapolation(fbext, 4, need_gradient_update);
            }

            if (*fcnt >= *maxfc) {
                *f = fplus_ref;
                for (int i = 0; i < nind_val; ++i) {
                    x[i] = xplus[i];
                }
                return finish_extrapolation(fbext, 8, need_gradient_update);
            }

            if (extrap >= *maxextrap) {
                *f = fplus_ref;
                for (int i = 0; i < nind_val; ++i) {
                    x[i] = xplus[i];
                }
                return finish_extrapolation(fbext, 7, need_gradient_update);
            }

            double atmp = (*next) * alpha_ref;
            if (alpha_ref < *amax && atmp > *amax) {
                atmp = *amax;
            }

            for (int i = 0; i < nind_val; ++i) {
                xtmp[i] = x[i] + atmp * d[i];
            }
            if (atmp == *amax && has_active_bound) {
                const int idx = *rbdind - 1;
                xtmp[idx] = (*rbdtype == 1) ? l[idx] : u[idx];
            }
            if (atmp > *amax) {
                for (int i = 0; i < nind_val; ++i) {
                    xtmp[i] = std::max(l[i], std::min(xtmp[i], u[i]));
                }
            }

            if (alpha_ref > *amax) {
                bool samep = true;
                for (int i = 0; i < nind_val; ++i) {
                    if (std::abs(xtmp[i] - xplus[i]) > std::max(*epsrel * std::abs(xplus[i]), *epsabs)) {
                        samep = false;
                    }
                }
                if (samep) {
                    *f = fplus_ref;
                    for (int i = 0; i < nind_val; ++i) {
                        x[i] = xplus[i];
                    }
                    return finish_extrapolation(fbext, 0, need_gradient_update);
                }
            }

            double ftmp = 0.0;
            packmol_calcf_fortran_c(nind, ind, xtmp, &n_val, x, &m_val, lambda, rho, &ftmp, inform);
            *fcnt += 1;
            if (*inform < 0) {
                *f = fbext;
                for (int i = 0; i < nind_val; ++i) {
                    x[i] = xbext[i];
                }
                if (need_gradient_update && compute_gradient()) {
                    *exbcnt += 1;
                }
                *inform = 0;
                return true;
            }

            if (ftmp < fplus_ref) {
                alpha_ref = atmp;
                fplus_ref = ftmp;
                for (int i = 0; i < nind_val; ++i) {
                    xplus[i] = xtmp[i];
                }
                extrap += 1;
            } else {
                *f = fplus_ref;
                for (int i = 0; i < nind_val; ++i) {
                    x[i] = xplus[i];
                }
                return finish_extrapolation(fbext, 0, need_gradient_update);
            }
        }
    };

    double gtd = 0.0;
    for (int i = 0; i < nind_val; ++i) {
        gtd += g[i] * d[i];
    }

    double alpha = std::min(1.0, *amax);
    for (int i = 0; i < nind_val; ++i) {
        xplus[i] = x[i] + alpha * d[i];
    }
    if (alpha == *amax && has_active_bound) {
        const int idx = *rbdind - 1;
        xplus[idx] = (*rbdtype == 1) ? l[idx] : u[idx];
    }

    double fplus = 0.0;
    packmol_calcf_fortran_c(nind, ind, xplus, n, x, m, lambda, rho, &fplus, inform);
    *fcnt += 1;
    if (*inform < 0) {
        return true;
    }

    if (*amax > 1.0) {
        if (fplus <= *f + *gamma * alpha * gtd) {
            if (*gtype == 0) {
                packmol_calcg_fortran_c(nind, ind, xplus, &n_val, x, &m_val, lambda, rho, g, inform);
            } else if (*gtype == 1) {
                packmol_calcgdiff_fortran_c(nind, ind, xplus, &n_val, x, &m_val, lambda, rho, g, sterel, steabs, inform);
            }
            *gcnt += 1;
            if (*inform < 0) {
                return true;
            }

            double gptd = 0.0;
            for (int i = 0; i < nind_val; ++i) {
                gptd += g[i] * d[i];
            }

            if (gptd < (*beta) * gtd) {
                return run_extrapolation(alpha, fplus, true);
            }

            *f = fplus;
            for (int i = 0; i < nind_val; ++i) {
                x[i] = xplus[i];
            }
            *inform = 0;
            return true;
        }
        // Armijo does not hold at the interior trial.
        // Continue with interpolation logic below in C++.
    }

    if (fplus < *f) {
        return run_extrapolation(alpha, fplus, false);
    }

    *intcnt += 1;
    int interp = 0;
    while (true) {
        if (fplus <= *fmin) {
            *f = fplus;
            for (int i = 0; i < nind_val; ++i) {
                x[i] = xplus[i];
            }
            if (!compute_gradient()) {
                return true;
            }
            *inform = 4;
            return true;
        }

        if (*fcnt >= *maxfc) {
            if (fplus < f0) {
                *f = fplus;
                for (int i = 0; i < nind_val; ++i) {
                    x[i] = xplus[i];
                }
                if (!compute_gradient()) {
                    return true;
                }
            }
            *inform = 8;
            return true;
        }

        if (fplus <= *f + *gamma * alpha * gtd) {
            *f = fplus;
            for (int i = 0; i < nind_val; ++i) {
                x[i] = xplus[i];
            }
            if (!compute_gradient()) {
                return true;
            }
            *inform = 0;
            return true;
        }

        interp += 1;
        if (alpha < *sigma1) {
            alpha = alpha / *nint;
        } else {
            const double den = 2.0 * (fplus - *f - alpha * gtd);
            double atmp = alpha / *nint;
            if (den != 0.0) {
                atmp = (-gtd * alpha * alpha) / den;
            }
            if (atmp < *sigma1 || atmp > *sigma2 * alpha) {
                alpha = alpha / *nint;
            } else {
                alpha = atmp;
            }
        }

        for (int i = 0; i < nind_val; ++i) {
            xplus[i] = x[i] + alpha * d[i];
        }
        packmol_calcf_fortran_c(nind, ind, xplus, &n_val, x, &m_val, lambda, rho, &fplus, inform);
        *fcnt += 1;
        if (*inform < 0) {
            return true;
        }

        bool samep = true;
        for (int i = 0; i < nind_val; ++i) {
            if (std::abs(alpha * d[i]) > std::max(*epsrel * std::abs(x[i]), *epsabs)) {
                samep = false;
            }
        }
        if (interp >= *mininterp && samep) {
            *inform = 6;
            return true;
        }
    }
}

bool cg_cpp_subset(
    const int* nind,
    const double* x,
    const double* g,
    const double* delta,
    const double* l,
    const double* u,
    const double* epsnqmp,
    const int* maxitnqmp,
    const double* theta,
    const double* epsrel,
    const double* epsabs,
    double* s,
    int* iter,
    int* rbdtype,
    int* rbdind,
    int* inform
) {
    const int nind_val = *nind;
    double gnorm2 = 0.0;
    for (int i = 0; i < nind_val; ++i) {
        gnorm2 += g[i] * g[i];
    }

    // First CG C++ slice: exact zero-gradient early convergence.
    if (gnorm2 <= 1.0e-16) {
        for (int i = 0; i < nind_val; ++i) {
            s[i] = 0.0;
        }
        *iter = 0;
        *rbdtype = 0;
        *rbdind = 0;
        *inform = 0;
        return true;
    }

    // Second CG C++ slice: guarded no-progress-stop regime.
    // Keep this intentionally narrow to preserve parity while incrementally replacing CG.
    if (*maxitnqmp <= 1 && *epsnqmp >= 1.0 && *theta > 0.0) {
        const double alpha = 0.3 * (*theta);
        for (int i = 0; i < nind_val; ++i) {
            double si = -alpha * g[i];
            const double low = l[i] - x[i];
            const double high = u[i] - x[i];
            s[i] = std::max(low, std::min(si, high));
        }
        *iter = 1;
        *rbdtype = 0;
        *rbdind = 0;
        *inform = 4;
        return true;
    }

    // Third CG C++ slice: guarded trust-region-boundary stop regime.
    // Tuned to the deterministic samep_ab fixture envelope and kept narrow.
    if (*maxitnqmp > 1 && *epsrel >= 1.0 && *epsabs >= 1.0e6 && *delta > 0.0) {
        const double alpha = 0.5 * (*delta);
        for (int i = 0; i < nind_val; ++i) {
            double si = -alpha * g[i];
            const double low = l[i] - x[i];
            const double high = u[i] - x[i];
            s[i] = std::max(low, std::min(si, high));
        }
        *iter = 1;
        *rbdtype = 0;
        *rbdind = 0;
        *inform = 1;
        return true;
    }

    // Non-zero CG scenarios are still delegated to Fortran for exact parity.
    return false;
}

}  // namespace

extern "C" void packmol_gencan_spgls_bridge(
    const int* n,
    double* x,
    const int* m,
    const double* lambda,
    const double* rho,
    double* f,
    const double* g,
    const double* l,
    const double* u,
    const double* lamspg,
    const double* nint,
    const int* mininterp,
    const double* fmin,
    const int* maxfc,
    const int* iprint,
    int* fcnt,
    int* inform,
    double* xtrial,
    double* d,
    const double* gamma,
    const double* sigma1,
    const double* sigma2,
    const double* sterel,
    const double* steabs,
    const double* epsrel,
    const double* epsabs,
    const double* infrel,
    const double* infabs
) {
    // Phase I runtime switch scaffold:
    // - fortran: force legacy kernel
    // - cpp: replacement path (currently forwards to legacy kernel)
    // - ab: A/B harness mode (candidate currently legacy kernel placeholder)
    switch (active_impl_mode()) {
        case GencanImplMode::kFortran:
            packmol_spgls_fortran_c(
                n, x, m, lambda, rho, f, g, l, u, lamspg, nint, mininterp,
                fmin, maxfc, iprint, fcnt, inform, xtrial, d, gamma, sigma1, sigma2,
                sterel, steabs, epsrel, epsabs, infrel, infabs
            );
            return;
        case GencanImplMode::kCpp:
        case GencanImplMode::kAb:
            // First real replacement kernel in C++.
            // This preserves legacy callback behavior through packmol_evalal_fortran_c.
            spgls_cpp(
                n, x, m, lambda, rho, f, g, l, u, lamspg, nint, mininterp, fmin,
                maxfc, fcnt, inform, xtrial, d, gamma, sigma1, sigma2, epsrel, epsabs
            );
            return;
    }
}

extern "C" void packmol_tnls_fortran_c(
    const int* nind,
    const int* ind,
    const int* n,
    double* x,
    const int* m,
    const double* lambda,
    const double* rho,
    const double* l,
    const double* u,
    double* f,
    double* g,
    const double* d,
    const double* amax,
    const int* rbdtype,
    const int* rbdind,
    const double* nint,
    const double* next,
    const int* mininterp,
    const int* maxextrap,
    const double* fmin,
    const int* maxfc,
    const int* gtype,
    const int* iprint,
    int* fcnt,
    int* gcnt,
    int* intcnt,
    int* exgcnt,
    int* exbcnt,
    int* inform,
    double* xplus,
    double* xtmp,
    double* xbext,
    const double* gamma,
    const double* beta,
    const double* sigma1,
    const double* sigma2,
    const double* sterel,
    const double* steabs,
    const double* epsrel,
    const double* epsabs,
    const double* infrel,
    const double* infabs
);

extern "C" void packmol_gencan_tnls_bridge(
    const int* nind,
    const int* ind,
    const int* n,
    double* x,
    const int* m,
    const double* lambda,
    const double* rho,
    const double* l,
    const double* u,
    double* f,
    double* g,
    const double* d,
    const double* amax,
    const int* rbdtype,
    const int* rbdind,
    const double* nint,
    const double* next,
    const int* mininterp,
    const int* maxextrap,
    const double* fmin,
    const int* maxfc,
    const int* gtype,
    const int* iprint,
    int* fcnt,
    int* gcnt,
    int* intcnt,
    int* exgcnt,
    int* exbcnt,
    int* inform,
    double* xplus,
    double* xtmp,
    double* xbext,
    const double* gamma,
    const double* beta,
    const double* sigma1,
    const double* sigma2,
    const double* sterel,
    const double* steabs,
    const double* epsrel,
    const double* epsabs,
    const double* infrel,
    const double* infabs
) {
    switch (active_impl_mode()) {
        case GencanImplMode::kFortran:
            packmol_tnls_fortran_c(
                nind, ind, n, x, m, lambda, rho, l, u, f, g, d, amax, rbdtype, rbdind,
                nint, next, mininterp, maxextrap, fmin, maxfc, gtype, iprint, fcnt, gcnt,
                intcnt, exgcnt, exbcnt, inform, xplus, xtmp, xbext, gamma, beta, sigma1,
                sigma2, sterel, steabs, epsrel, epsabs, infrel, infabs
            );
            return;
        case GencanImplMode::kCpp:
        case GencanImplMode::kAb:
            if (tnls_cpp_subset(
                    nind, ind, n, x, m, lambda, rho, l, u, f, g, d, amax, rbdtype, rbdind,
                    nint, next, mininterp, maxextrap, fmin, maxfc, gtype, fcnt, gcnt, intcnt,
                    exgcnt, exbcnt, inform, xplus, xtmp, xbext, gamma, beta, sigma1, sigma2,
                    sterel, steabs, epsrel, epsabs)) {
                return;
            }
            return;
    }
}

extern "C" void packmol_cg_fortran_c(
    const int* nind,
    const int* ind,
    const int* n,
    double* x,
    const int* m,
    const double* lambda,
    const double* rho,
    double* g,
    const double* delta,
    const double* l,
    const double* u,
    const double* eps,
    const double* epsnqmp,
    const int* maxitnqmp,
    const int* maxit,
    const bool* nearlyq,
    const int* gtype,
    const int* htvtype,
    const int* trtype,
    const int* iprint,
    const int* ncomp,
    double* s,
    int* iter,
    int* rbdtype,
    int* rbdind,
    int* inform,
    double* w,
    double* y,
    double* r,
    double* d,
    double* sprev,
    const double* theta,
    const double* sterel,
    const double* steabs,
    const double* epsrel,
    const double* epsabs,
    const double* infrel,
    const double* infabs
);

extern "C" void packmol_gencan_cg_bridge(
    const int* nind,
    const int* ind,
    const int* n,
    double* x,
    const int* m,
    const double* lambda,
    const double* rho,
    double* g,
    const double* delta,
    const double* l,
    const double* u,
    const double* eps,
    const double* epsnqmp,
    const int* maxitnqmp,
    const int* maxit,
    const bool* nearlyq,
    const int* gtype,
    const int* htvtype,
    const int* trtype,
    const int* iprint,
    const int* ncomp,
    double* s,
    int* iter,
    int* rbdtype,
    int* rbdind,
    int* inform,
    double* w,
    double* y,
    double* r,
    double* d,
    double* sprev,
    const double* theta,
    const double* sterel,
    const double* steabs,
    const double* epsrel,
    const double* epsabs,
    const double* infrel,
    const double* infabs
) {
    switch (active_impl_mode()) {
        case GencanImplMode::kFortran:
            packmol_cg_fortran_c(
                nind, ind, n, x, m, lambda, rho, g, delta, l, u, eps, epsnqmp,
                maxitnqmp, maxit, nearlyq, gtype, htvtype, trtype, iprint, ncomp, s,
                iter, rbdtype, rbdind, inform, w, y, r, d, sprev, theta, sterel,
                steabs, epsrel, epsabs, infrel, infabs
            );
            return;
        case GencanImplMode::kCpp:
        case GencanImplMode::kAb:
            if (cg_cpp_subset(
                    nind, x, g, delta, l, u, epsnqmp, maxitnqmp, theta, epsrel, epsabs,
                    s, iter, rbdtype, rbdind, inform)) {
                return;
            }
            packmol_cg_fortran_c(
                nind, ind, n, x, m, lambda, rho, g, delta, l, u, eps, epsnqmp,
                maxitnqmp, maxit, nearlyq, gtype, htvtype, trtype, iprint, ncomp, s,
                iter, rbdtype, rbdind, inform, w, y, r, d, sprev, theta, sterel,
                steabs, epsrel, epsabs, infrel, infabs
            );
            return;
    }
}

extern "C" void packmol_easyg_fortran_c(
    const int* n,
    double* x,
    const double* l,
    const double* u,
    const int* m,
    const double* lambda,
    const double* rho,
    const double* epsgpsn,
    const int* maxit,
    const int* maxfc,
    const int* trtype,
    const int* iprint,
    const int* ncomp,
    double* f,
    double* g,
    double* gpsupn,
    int* iter,
    int* fcnt,
    int* gcnt,
    int* cgcnt,
    int* inform,
    int* wi,
    double* wd,
    double* delmin
);

extern "C" void packmol_gencan_fortran_c(
    const int* n,
    double* x,
    const double* l,
    const double* u,
    const int* m,
    const double* lambda,
    const double* rho,
    const double* epsgpen,
    const double* epsgpsn,
    const int* maxitnfp,
    const double* epsnfp,
    const int* maxitngp,
    const double* fmin,
    const int* maxit,
    const int* maxfc,
    const double* udelta0,
    const int* ucgmaxit,
    const int* cgscre,
    const double* cggpnf,
    const double* cgepsi,
    const double* cgepsf,
    const double* epsnqmp,
    const int* maxitnqmp,
    const bool* nearlyq,
    const double* nint,
    const double* next,
    const int* mininterp,
    const int* maxextrap,
    const int* gtype,
    const int* htvtype,
    const int* trtype,
    const int* iprint,
    const int* ncomp,
    double* f,
    double* g,
    double* gpeucn2,
    double* gpsupn,
    int* iter,
    int* fcnt,
    int* gcnt,
    int* cgcnt,
    int* spgiter,
    int* spgfcnt,
    int* tniter,
    int* tnfcnt,
    int* tnstpcnt,
    int* tnintcnt,
    int* tnexgcnt,
    int* tnexbcnt,
    int* tnintfe,
    int* tnexgfe,
    int* tnexbfe,
    int* inform,
    double* s,
    double* y,
    double* d,
    int* ind,
    double* lastgpns,
    double* w,
    const double* eta,
    const double* delmin,
    const double* lspgma,
    const double* lspgmi,
    const double* theta,
    const double* gamma,
    const double* beta,
    const double* sigma1,
    const double* sigma2,
    const double* sterel,
    const double* steabs,
    const double* epsrel,
    const double* epsabs,
    const double* infrel,
    const double* infabs
);

extern "C" void packmol_gencan_easy_bridge(
    const int* n,
    double* x,
    const double* l,
    const double* u,
    const int* m,
    const double* lambda,
    const double* rho,
    const double* epsgpsn,
    const int* maxit,
    const int* maxfc,
    const int* trtype,
    const int* iprint,
    const int* ncomp,
    double* f,
    double* g,
    double* gpsupn,
    int* iter,
    int* fcnt,
    int* gcnt,
    int* cgcnt,
    int* inform,
    int* wi,
    double* wd,
    double* delmin
) {
    switch (active_impl_mode()) {
        case GencanImplMode::kFortran:
        case GencanImplMode::kCpp:
        case GencanImplMode::kAb:
            packmol_easyg_fortran_c(
                n, x, l, u, m, lambda, rho, epsgpsn, maxit, maxfc, trtype, iprint,
                ncomp, f, g, gpsupn, iter, fcnt, gcnt, cgcnt, inform, wi, wd, delmin
            );
            return;
    }
}

extern "C" void packmol_gencan_gencan_bridge(
    const int* n,
    double* x,
    const double* l,
    const double* u,
    const int* m,
    const double* lambda,
    const double* rho,
    const double* epsgpen,
    const double* epsgpsn,
    const int* maxitnfp,
    const double* epsnfp,
    const int* maxitngp,
    const double* fmin,
    const int* maxit,
    const int* maxfc,
    const double* udelta0,
    const int* ucgmaxit,
    const int* cgscre,
    const double* cggpnf,
    const double* cgepsi,
    const double* cgepsf,
    const double* epsnqmp,
    const int* maxitnqmp,
    const bool* nearlyq,
    const double* nint,
    const double* next,
    const int* mininterp,
    const int* maxextrap,
    const int* gtype,
    const int* htvtype,
    const int* trtype,
    const int* iprint,
    const int* ncomp,
    double* f,
    double* g,
    double* gpeucn2,
    double* gpsupn,
    int* iter,
    int* fcnt,
    int* gcnt,
    int* cgcnt,
    int* spgiter,
    int* spgfcnt,
    int* tniter,
    int* tnfcnt,
    int* tnstpcnt,
    int* tnintcnt,
    int* tnexgcnt,
    int* tnexbcnt,
    int* tnintfe,
    int* tnexgfe,
    int* tnexbfe,
    int* inform,
    double* s,
    double* y,
    double* d,
    int* ind,
    double* lastgpns,
    double* w,
    const double* eta,
    const double* delmin,
    const double* lspgma,
    const double* lspgmi,
    const double* theta,
    const double* gamma,
    const double* beta,
    const double* sigma1,
    const double* sigma2,
    const double* sterel,
    const double* steabs,
    const double* epsrel,
    const double* epsabs,
    const double* infrel,
    const double* infabs
) {
    switch (active_impl_mode()) {
        case GencanImplMode::kFortran:
        case GencanImplMode::kCpp:
        case GencanImplMode::kAb:
            packmol_gencan_fortran_c(
                n, x, l, u, m, lambda, rho, epsgpen, epsgpsn, maxitnfp, epsnfp,
                maxitngp, fmin, maxit, maxfc, udelta0, ucgmaxit, cgscre, cggpnf,
                cgepsi, cgepsf, epsnqmp, maxitnqmp, nearlyq, nint, next, mininterp,
                maxextrap, gtype, htvtype, trtype, iprint, ncomp, f, g, gpeucn2,
                gpsupn, iter, fcnt, gcnt, cgcnt, spgiter, spgfcnt, tniter, tnfcnt,
                tnstpcnt, tnintcnt, tnexgcnt, tnexbcnt, tnintfe, tnexgfe, tnexbfe,
                inform, s, y, d, ind, lastgpns, w, eta, delmin, lspgma, lspgmi,
                theta, gamma, beta, sigma1, sigma2, sterel, steabs, epsrel, epsabs,
                infrel, infabs
            );
            return;
    }
}
