'use strict'; // -*- js2 -*-

/** @module lib/rpoly
 *
 * @file
 * Implementation of the Jenkins--Traub algorithm for polynomial factorization.
 *
 * This is a (hopefully faithful) translation of the Fortran implementation:
 *
 *    http://www.netlib.org/toms/493
 *
 * The algorithm has been written in idiomatic Javascript---especially when
 * working with arrays---and includes minor simplifications made possible by
 * eliminating the shared variables used in 493.FOR.
 */


const {abs, max, min, log, exp, floor, ceil, sqrt, sin, cos, PI} = Math;

// Number.MIN_VALUE is the minimum *subnormal* positive integer.
const SMALNO = 2.2250738585072014e-308;
const INFIN = Number.MAX_VALUE;
const ETA = Number.EPSILON;
const ARE = ETA;
const MRE = ETA;
const LO = SMALNO / ETA;
const LOG2 = log(2);
const COSR = cos(94*PI/180); // ~ -0.069756474
const SINR = sin(94*PI/180); // ~ 0.99756405
const SQRT2INV = sqrt(0.5); // ~ 0.70710678

/**
 * A polynomial is represented as a list of real coefficients.
 *
 * The zeroth coefficient is the highest-order term.
 *
 * @typedef {number[]} Polynomial
 */

/**
 * A complex number is represented as an array `[Re, Im]`.
 *
 * @typedef {number[]} Complex
 */

/**
 * Compute the roots of the polynomial `P`.
 *
 * The return value is an Array of roots.  If this Array has fewer roots than
 * the degree of `P`, then the algorithm has failed to find all roots.
 *
 * @param {Polynomial} P - The polynomial to factor.
 * @return {Complex[]} The roots found.
 */
function rpoly(P) {
    // Make a copy of the coefficients
    P = P.slice();
    // Translate so the leading coefficient is nonzero
    while(P[0] === 0) P.shift();

    let nn = P.length;
    let n = nn - 1;
    let zeros = [];

    // Remove the zeros at the origin if any
    while(P[n] === 0) {
        zeros.push([0, 0]);
        P.pop();
        nn -= 1;
        n  -= 1;
    }

    let xx = SQRT2INV;
    let yy = -xx;

    while(n >= 1) {// Main loop
        if(n == 1) {
            zeros.push([-P[1]/P[0], 0]);
            break;
        } else if(n == 2) {
            zeros.push(...quad(...P));
            break;
        }

        // Find the largest and smallest moduli of coefficients
        let PT = P.map(abs);
        const maxc = max(...PT);
        const minc = min(...PT.filter(x => x !== 0));

        // Scale if there are large or very small coefficients.  Computes a
        // scale factor to multiply the coefficients of the polynomial.  The
        // scaling is done to avoid overflow and to avoid undetected
        // underflow interfering with the convergence criterion.
        const sc = LO/minc;
        if((sc > 1 && INFIN/sc >= maxc) || (sc <= 1 && maxc >= 10)) {
            if(sc === 0) sc = SMALNO;
            // Round towards zero
            let l = log(sc)/LOG2 + 0.5;
            l = (l > 0) ? floor(l) : ceil(l);
            if(l != 0) {
                const factor = 2**l;
                for(let i = 0; i < P.length; ++i) {
                    P[i] *= factor;
                    PT[i] *= factor;
                }
            }
        }

        // Compute lower bound on moduli of zeros.
        PT[n] *= -1;

        // Compute upper estimate of bound
        let x = exp((log(-PT[n])-log(PT[0]))/n);
        if(PT[n-1] !== 0)
            // If Newton step at the origin is better, use it.
            x = min(x, -PT[n]/PT[n-1]);

        // Chop the interval (0,X) until ff <= 0
        let xm = x * 0.1;
        while(PT.reduce((a, x) => a*xm + x, 0) > 0) {
            x = xm;
            xm *= 0.1;
        }

        // Do Newton iteration until X converges to two decimal places
        let dx = x;
        while(abs(dx/x) > .005) {
            let [ff, df] = PT.reduce(
                ([ff, df], xx) => [ff*x + xx, df*x + ff], [0, 0]);
            dx = ff/df;
            x -= dx;
        }
        const bnd = x;

        // Compute the derivative as the initial K polynomial and do 5 steps
        // with no shift
        let K = P.map((x, i) => (n-i) * x / n);
        K.pop();

        const aa = P[n];
        const bb = P[n-1];
        let zerok = K[n-1] === 0;
        for(let jj = 0; jj < 5; ++jj) {
            const cc = K[n-1];
            if(zerok) {
                // Use unscaled form of recurrence
                K.unshift(0);
                K.pop();
                zerok = K[n-1] === 0;
            } else {
                // Use scaled form of recurrence if value of K at 0 is nonzero
                const t = -aa/cc;
                K = P.map((x, i) => t * (i > 0 ? K[i-1] : 0) + x);
                K.pop();
                zerok = abs(K[n-1]) <= abs(bb)*ETA*10;
            }
        }

        // Loop to select the quadratic corresponding to each new shift
        let cnt;
        for(cnt = 1; cnt <= 20; ++cnt) {
            // Quadratic corresponds to a double shift to a non-real point and its
            // complex conjugate.  The point has modulus BND and amplitude rotated
            // by 94 degrees from the previous shift
            [xx, yy] = [COSR*xx - SINR*yy, SINR*xx + COSR*yy];
            const u = -2*bnd*xx;
            const v = bnd;

            const [newZeros, QP] = fxshfr(20*cnt, P, K, u, v);
            if(newZeros.length !== 0) {
                // The second stage jumps directly to one of the third stage
                // iterations and returns here if successful.  Deflate the
                // polynomial, store the zero or zeros and return to the main
                // algorithm.
                zeros.push(...newZeros);
                nn -= newZeros.length;
                n   = nn - 1;
                P   = QP.slice(0, nn);
                break;
            }
        }
        if(cnt > 20)
            break; // Return with failure if no convergence with 20 shifts
    }

    return zeros;
}

/**
 * Compute up to `l2` fixed shift K-polynomials, testing for convergence in the
 * linear or quadratic case.  Initiates one of the variable shift iterations and
 * return with the zeros found.
 *
 * @param {integer} l2 - The number of K-polynomials to compute.
 * @param {Polynomial} P - The polynomial.
 * @param {Polynomial} K - The K-polynomial.
 * @param {number} u - The starting linear coefficient.
 * @param {number} v - The starting constant coefficient.
 * @return {Complex[]} The zeros found.
 *
 * @private
 */
function fxshfr(l2, P, K, u, v) {
    const n = K.length;
    K = K.slice();
    let betav = 0.25;
    let betas = 0.25;
    // Evaluate polynomial by synthetic division
    const [a, b, QP] = quadsd(u, v, P);
    let [type, QK, scalars] = calcsc(u, v, a, b, K);

    let ovv, oss, otv, ots;
    for(let j = 0; j < l2; ++j) {
        // Calculate next K polynomial and estimate V
        K = nextk(type, K, QK, QP, a, b, scalars);
        [type, QK, scalars] = calcsc(u, v, a, b, K);
        let [ui, vi] = newest(type, K, P, a, b, u, v, scalars);
        const vv = vi;
        // Estimate S
        const ss = K[n-1] !== 0 ? -P[n]/K[n-1] : 0;
        let tv = 1, ts = 1;

        if(j > 0 && type !== 3) {
            // Compute relative measures of convergence of S and V sequences
            if(vv !== 0) tv = abs((vv-ovv)/vv);
            if(ss !== 0) ts = abs((ss-oss)/ss);
            // If decreasing, multiply two most recent convergence measures
            const tvv = tv < otv ? tv*otv : 1;
            const tss = ts < ots ? ts*ots : 1;
            // Compare with convergence criteria
            const vpass = tvv < betav;
            const spass = tss < betas;

            if(spass || vpass) {
                // At least one sequence has passed the convergence test.
                let s = ss;
                // Choose iteration according to the fastest converging
                // sequence.
                let vtry = false;
                let stry = false;
                // 1: S; 2: V
                let which = spass && (!vpass || tss < tvv) ? 1 : 2;

                while(true) {
                    if(which == 1) {
                        const ret = realit(s, P, K);
                        if(ret.success)
                            return [[[ret.s, 0]], ret.QP];
                        // Linear iteration has failed.  Flag that it has been tried
                        // and decrease the convergence criterion
                        stry = true;
                        betas *= 0.25;
                        if(ret.iflag) {
                            s = ret.s;
                            K = ret.K;
                            // If linear iteration signals an almost double real
                            // zero attempt quadratic interation
                            ui = -(s+s);
                            vi = s*s;
                            which = 2;
                            continue;
                        }
                    } else if(which == 2) {
                        const [zeros, QP] = quadit(ui, vi, P, K);
                        if(zeros) return [zeros, QP];
                        // Quadratic iteration has failed.  Flag that it has been
                        // tried and decrease the convergence criterion.
                        vtry = true;
                        betav *= 0.25;
                        if(spass && !stry) {
                            which = 1;
                            continue;
                        }
                    }

                    if(vpass && !vtry) which = 2;
                    else break;
                }
            }
        }

        ovv = vv;
        oss = ss;
        otv = tv;
        ots = ts;
    }

    return [[], QP];
}

/**
 * Variable-shift K-polynomial iteration for a quadratic factor converges
 * only if the zeros are equimodular or nearly so.
 *
 * @param {number} u - The starting linear coefficient.
 * @param {number} v - The starting constant coefficient.
 * @param {Polynomial} P - The polynomial.
 * @param {Polynomial} K - The K-polynomial.
 * @return {Array} An array `[roots, QP]`, where `roots` contains the two zeros
 *   found and `QP` is the quotient of `P` by `u*x+v`.  If unsuccessful, return
 *   `[null, null]`.
 *
 * @private
 */
function quadit(u, v, P, K) {
    K = K.slice();
    const n = K.length;
    let tried = false;
    let j = 0;

    let relstp, omp;
    while(true) {
        // Main loop
        const roots = quad(1, u, v);
        const [[szr, szi], [lzr, lzi]] = roots;
        // Return if roots of the quadratic are real and not close to multiple
        // or nearly equal and of opposite sign
        if(abs(abs(szr)-abs(lzr)) > 0.01*abs(lzr))
            return [null, null];
        // Evaluate polynomial by quadratic synthetic division
        let [a, b, QP] = quadsd(u, v, P);
        const mp = abs(a-szr*b) + abs(szi*b);
        // Compute a rigorous bound on the rounding error in evaluting P
        const zm = sqrt(abs(v));
        let ee = 2*abs(QP[0]);
        for(let i = 1; i < n; ++i)
            ee = ee*zm + abs(QP[i]);
        const t = -szr*b;
        ee = ee*zm + abs(a+t);
        ee = (5*MRE + 4*ARE) * ee
            - (5*MRE + 2*ARE) * (abs(a+t)+abs(b)*zm)
            + 2*ARE*abs(t);
        // Iteration has converged sufficiently if the polynomial value is less
        // than 20 times this bound
        if(mp <= 20*ee)
            return [roots, QP];

        // Stop iteration after 20 steps
        if(++j > 20) return [null, null];

        if(j >= 2 && relstp <= .01 && mp >= omp && !tried) {
            // A cluster appears to be stalling the convergence.  Five fixed
            // shift steps are taken with a U,V close to the cluster.
            if(relstp < ETA) relstp = ETA;
            relstp = sqrt(relstp);
            u -= u*relstp;
            v += v*relstp;
            [a, b, QP] = quadsd(u, v, P);
            for(let i = 0; i < 5; ++i) {
                const [type, QK, scalars] = calcsc(u, v, a, b, K);
                K = nextk(type, K, QK, QP, a, b, scalars);
            }
            tried = true;
            j = 0;
        }

        // Calculate next K polynomial and new U and V
        omp = mp;
        let [type, QK, scalars] = calcsc(u, v, a, b, K);
        K = nextk(type, K, QK, QP, a, b, scalars);
        [type, QK, scalars] = calcsc(u, v, a, b, K);
        const [ui, vi] = newest(type, K, P, a, b, u, v, scalars);
        // If vi is zero the iteration is not converging
        if(vi === 0) return [null, null];
        relstp = abs((vi-v)/vi);
        u = ui;
        v = vi;
    }
}

/**
 * Variable-shift H-polynomial iteration for a real zero.
 *
 * @param {number} s - The starting point.
 * @param {Polynomial} P - The polynomial.
 * @param {Polynomial} K - The K-polynomial.
 * @return {Object} If iteration does not converge, return `{}`.  If a cluster
 *   of zeros near the real axis is encountered, return `{s, K, iflag: true}` to
 *   indicate that quadratic iteration should be tried, where `s` determines the
 *   value of `u` and `v` as in {@link fxshfr} and `K` is the current
 *   K-polynomial.  If iteration converges, return `{success: true, s, QP}`,
 *   where `s` is the zero and `QP` is the quotient polynomial.
 *
 * @private
 */
function realit(s, P, K) {
    K = K.slice();
    const n = K.length, nn = P.length;
    let j = 0;
    let t, omp;
    let QP = new Array(P.length);
    let QK = new Array(K.length);
    while(true) {
        // Main loop
        const pv = P.reduce((a, x, i) => (QP[i] = a*s + x), 0);
        const mp = abs(pv);
        // Compute a rigorous bound on the error in evaluating P
        const ms = abs(s);
        let ee = (MRE/(ARE+MRE))*abs(QP[0]);
        for(let i = 1; i < nn; ++i)
            ee = ee*ms + abs(QP[i]);

        // Iteration has converged sufficiently if the polynomial value is less
        // than 20 times this bound
        if(mp <= 20*((ARE+MRE)*ee-MRE*mp))
            return {success: true, s, QP};
        // Stop iteration after 10 steps
        if(++j > 10) return {};
        if(j >= 2 && abs(t) <= 0.001*abs(s-t) && mp > omp)
            // A cluster of zeros near the real axis has been encountered.
            // Return with IFLAG set to initiate a quadratic iteration.
            return {s, K, iflag: true};

        // Return if the polynomial value has increased significantly
        omp = mp;
        // Compute T, the next polynomial, and the new iterate
        let kv = K.reduce((a, x, i) => (QK[i] = a*s + x), 0);
        if(abs(kv) > abs(K[n-1])*10*ETA) {
            // Use the scaled form of the recurrence if the value of K at S is
            // nonzero
            t = -pv/kv;
            K[0] = QP[0];
            for(let i = 1; i < n; ++i)
                K[i] = t*QK[i-1] + QP[i];
        } else
            // Use unscaled form
            K = [0, ...QK.slice(0, n-1)];
        kv = K.reduce((a, x) => a*s + x, 0);
        t = abs(kv) > abs(K[n-1])*10*ETA ? -pv/kv : 0;
        s += t;
    }
}

/**
 * Calculate the scalar quantities used to compute the next K polynomial and new
 * estimates of the quadratic coefficients.
 *
 * @param {number} u - The linear coefficient.
 * @param {number} v - The constant coefficient.
 * @param {number} a - The linear coefficient of the remainder.
 * @param {number} b - The constant coefficient of the remainder.
 * @param {Polynomial} K - The current K-polynomial.
 * @return {Array} Return `[type, QK, scalars]`, where `type` is described
 *   below, `QK` is the quotient of `K` by `u*x + v`, and `scalars` contains the
 *   computed scalar quantities `{c, d, f, g, h, a1, a3, a7}`.
 *
 * @private
 */
function calcsc(u, v, a, b, K) {
    // Synthetic division of K by the quadratic 1,U,V
    const n = K.length;
    const [c, d, QK] = quadsd(u, v, K);
    if(abs(c) <= abs(K[n-1])*100*ETA && abs(d) <= abs(K[n-2])*100*ETA)
        // type 3 indicates the quadratic is almost a factor of K
        return [3, QK, {c, d}];
    let type;
    let e, f, g, h, a1, a3, a7;
    if(abs(d) >= abs(c)) {
        e = a/d;
        f = c/d;
        g = u*b;
        h = v*b;
        a3 = (a+g)*e + h*(b/d);
        a1 = b*f - a;
        a7 = (f+u)*a + h;
        // type 2 indicates that all formulas are divided by D
        type = 2;
    } else {
        e = a/c;
        f = d/c;
        g = u*e;
        h = v*b;
        a3 = a*e + (h/c+g)*b;
        a1 = b - a*(d/c);
        a7 = a + g*d + h*f;
        // type 1 indicates that all formulas are divided by C
        type = 1;
    }
    return [type, QK, {c, d, f, g, h, a1, a3, a7}];
}

/**
 * Compute the next K polynomial using scalars computed in {@link calcsc}.
 *
 * Note that `K` and `scalars` are modified in-place.
 *
 * @param type - The type from {@link calcsc}.
 * @param {Polynomial} K - The current K-polynomial.
 * @param {Polynomial} QK - The quotient K-polynomial.
 * @param {Polynomial} QP - The quotient polynomial.
 * @param {number} a - The linear coefficient of the remainder.
 * @param {number} b - The constant coefficient of the remainder.
 * @param {Object} scalars - The scalars from {@link calcsc}.
 * @return {Polynomial} The new K-polynomial.
 *
 * @private
 */
function nextk(type, K, QK, QP, a, b, scalars) {
    const n = K.length;
    if(type == 3) {
        // Use unscaled form of the recurrence if type is 3
        return [0, 0, ...QK.slice(n-2)];
    }
    let {a1, a3, a7} = scalars;
    if(abs(a1) > abs(type === 1 ? b : a) * ETA * 10) {
        // Use scaled form of the recurrence
        a7 /= a1; scalars.a7 = a7;
        a3 /= a1; scalars.a3 = a3;
        K[0] = QP[0];
        K[1] = QP[1] - a7*QP[0];
        for(let i = 2; i < n; ++i)
            K[i] = a3*QK[i-2] - a7*QP[i-1] + QP[i];
        return K;
    }
    // If A1 is nearly zero then use a special form of the recurrence
    K[0] = 0;
    K[1] = -a7*QP[1];
    for(let i = 2; i < n; ++i)
        K[i] = a3*QK[i-2] - a7*QP[i-1];
    return K;
}

/**
 * Compute new estimates of the quadratic coefficients using the scalars
 * computed in {@link calcsc}.
 *
 * @param type - The type from {@link calcsc}.
 * @param {Polynomial} K - The current K-polynomial.
 * @param {Polynomial} P - The polynomial.
 * @param {number} a - The linear coefficient of the remainder.
 * @param {number} b - The constant coefficient of the remainder.
 * @param {number} u - The linear coefficient.
 * @param {number} v - The constant coefficient.
 * @param {Object} scalars - The scalars from {@link calcsc}.
 * @return {number[]} An array `[u, v]` containing the new quadratic
 *   coefficients.
 *
 * @private
 */
function newest(type, K, P, a, b, u, v, scalars) {
    if(type === 3)
        // If type==3 the quadratic is zeroed
        return [0, 0];
    // Use formulas appropriate to setting of type.
    const {a1, a3, a7, c, d, f, g, h} = scalars;
    let a4, a5;
    if(type === 2) {
        a4 = (a+g)*f + h;
        a5 = (f+u)*c + v*d;
    } else {
        a4 = a + u*b + h*f;
        a5 = c + (u+v*f)*d;
    }
    // Evaluate new quadratic coefficients.
    const n = K.length, nn = P.length;
    const b1 = -K[n-1]/P[nn-1];
    const b2 = -(K[n-2]+b1*P[n-1])/P[nn-1];
    const c1 = v*b2*a1;
    const c2 = b1*a7;
    const c3 = b1*b1*a3;
    const c4 = c1 - c2 - c3;
    const temp = a5 + b1*a4 - c4;
    if(temp === 0)
        return [0, 0];
    const uu = u - (u*(c3+c2)+v*(b1*a1+b2*a7))/temp;
    const vv = v*(1.+c4/temp);
    return [uu, vv];
}

/**
 * Synthetic division of `P` by the quadratic `x^2 + u*x + v`.
 *
 * @param {number} u - The linear coefficient.
 * @param {number} v - The constant coefficient.
 * @param {Polynomial} P - The polynomial to divide.
 * @return {Array} An array `[a, b, Q]` containing the remainder and the
 *   quotient polynomial.  The polynomial `Q` has the same number of entries as
 *   `P`; the last two entries are `a` and `b` again.
 *
 * @private
 */
function quadsd(u, v, P) {
    let a = 0, b = 0;
    let Q = P.map(x => {
        [a, b] = [x - u*a - v*b, a];
        return a;
    });
    return [a, b, Q];
}

/**
 * Calculate the zeros of the quadratic a*x^2+b1*x+c.
 *
 * The quadratic formula, modified to avoid overflow, is used to find the larger
 * zero if the zeros are real and both zeros are complex.  The smaller real zero
 * is found directly from the product of the zeros c/a.
 *
 * Returns `[[sr, si], [lr, li]]`, where `[sr, si]` is the smaller root and
 * `[lr, li]` is the larger.
 *
 * @param {number} a  - The quadratic coefficient.
 * @param {number} b1 - The linear coefficient.
 * @param {number} c  - The constant coefficient.
 * @return {Complex[]} The roots.
 *
 * @private
 */
function quad(a, b1, c) {
    if(a === 0)
        return [[b1 === 0 ? 0 : -c/b1, 0], [0, 0]];
    if(c === 0)
        return [[0, 0], [-b1/a, 0]];
    // Compute discriminant avoiding overflow
    const b = b1/2;
    let d, e;
    if(abs(b) >= abs(c)) {
        e = 1 - (a/b) * (c/b);
        d = sqrt(abs(e)) * abs(b);
    } else {
        e = c < 0 ? -a : a;
        e = b * (b/abs(c)) - e;
        d = sqrt(abs(e)) * sqrt(abs(c));
    }
    if(e >= 0) {
        // Real zeros
        if(b >= 0) d = -d;
        const lr = (-b+d)/a;
        return [[lr !== 0 ? (c/lr)/a : 0, 0], [lr, 0]];
    }
    // Complex conjugate zeros
    const sr = -b/a, si = abs(d/a);
    return [[sr, si], [sr, -si]];
}

export default rpoly;
