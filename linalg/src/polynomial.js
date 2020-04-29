'use strict';

/** @module polynomial
 *
 * @file
 * A class for polynomial manipulations.
 */

const { abs, sqrt, cbrt, atan2, cos, min, max } = Math;
const π = Math.PI;

import Complex from './complex.js';
import { range } from "./util.js";

// Convenience
const C = (a, b=0) => new Complex(a, b);

// Memoized Legendre polynomials
let legendres;

// A cube root of unity
const ζ = new Complex(-1/2, sqrt(3)/2);


/**
 * @summary
 * Sort an Array of Roots according to the ordering described in {@link
 * Polynomial#factor}.
 *
 * @desc
 * This assumes that complex roots come in conjugate pairs.
 *
 * @param {number} roots - The roots to sort.
 * @param {number} [translate=0] - Translate the roots by this real number.
 * @return {Root[]} The sorted roots.
 *
 * @private
 */
function sortRoots(roots, translate) {
    // Take one of each conjugate pair with positive real part
    return roots
        .map(([x, m]) => [C(x), m])
        .filter(([x, ]) => x.Im >= 0)
        .sort(([x,], [y,]) => x.Re === y.Re ? x.Im - y.Im : x.Re - y.Re)
        .flatMap(([x, m]) => x.Im === 0
                 ? [[x.Re+translate, m]]
                 : [[x.clone().add(translate), m], [x.clone().conj().add(translate), m]]);
}


/**
 * @summary
 * Calculate the zeros of a quadratic polynomial.
 *
 * @desc
 * This function finds the roots of the polynomial `x^2 + b x + c` using the
 * quadratic formula.
 *
 * @param {number} b - The linear coefficient.
 * @param {number} c - The constant coefficient.
 * @param {number} [ε=1e-10] - The roots will be considered equal if the
 *   discriminant is smaller than this.
 * @return {Root[]} The roots.  See {@link Polynomial#factor} for the ordering.
 *
 * @private
 */
function quadratic(b, c, ε=1e-10) {
    // Handle roots at zero
    if(abs(c) < ε) {
        if(abs(b) < ε)
            return [[0, 2]];
        return b < 0 ? [[0, 1], [-b, 1]] : [[-b, 1], [0, 1]];
    }
    let Δ = b*b-4*c;
    if(abs(Δ) <= ε)
        return [[-b/2, 2]]; // Double real root
    if(Δ > 0) {
        // Distinct real roots
        let D = sqrt(Δ);
        return [[-b/2 - D/2, 1], [-b/2 + D/2, 1]];
    }
    // Complex roots
    let D = sqrt(-Δ);
    return [[C(-b/2, D/2), 1], [C(-b/2, -D/2), 1]];
}


/**
 * @summary
 * Find the roots of a cubic polynomial.
 *
 * @desc
 * This function finds the roots of the polynomial `x^3 + b x^2 + c x + d` using
 * [Cardano's formula]{@link
 * https://www.encyclopediaofmath.org/index.php/Cardano_formula}.
 *
 * @example
 * cardano(-4, 5, -2); // [[1, 2], [2, 1]] (approximately)
 *
 * @param {number} b - The quadratic coefficient.
 * @param {number} c - The linear coefficient.
 * @param {number} d - The constant coefficient.
 * @param {number} [ε=1e-10] - The roots will be considered equal if the
 *   discriminant is smaller than this.
 * @return {Root[]} The roots.  See {@link Polynomial#factor} for the ordering.
 *
 * @private
 */
function cardano(b, c, d, ε=1e-10) {
    // Change of variables x --> x-b/3
    const [p, q] = [-1/3*b*b+c, 2/27*b*b*b - 1/3*b*c + d];
    // Now we're solving x^3 + px + q

    // Handle roots at zero
    if(abs(q) <= ε) {
        if(abs(p) <= ε)
            return [[-b/3, 3]];  // Triple root
        let s = sqrt(abs(p));
        if(p < 0)
            return [[-s-b/3, 1], [-b/3, 1], [s-b/3, 1]];
        return [[-b/3, 1], [C(-b/3, -s), 1], [C(-b/3, s), 1]];
    }

    // Discriminant
    const Δ = -27*q*q - 4*p*p*p;
    if(abs(Δ) <= ε) {
        // Simple root and double root
        const cr = cbrt(-q/2);
        return cr < 0
            ? [[2*cr-b/3, 1], [ -cr-b/3, 2]]
            : [[ -cr-b/3, 2], [2*cr-b/3, 1]];
    }

    if(Δ > 0) {
        // Three distinct real roots: 2*Re(cube roots of -q/2 + i sqrt(Δ/108))
        const D = sqrt(Δ/108);
        const mod = sqrt(cbrt(q*q/4 + Δ/108)) * 2;
        const arg = atan2(D, -q/2);
        return [cos(arg/3), cos(arg/3 + 2*π/3), cos(arg/3 - 2*π/3)]
            .sort((x, y) => x-y).map(x => [mod*x-b/3, 1]);
    }

    // Simple real root and conjugate complex roots
    const D = sqrt(-Δ/108);
    const α = cbrt(-q/2 + D), β = cbrt(-q/2 - D), r = α + β - b/3;
    let ret = [[C(α + β), 1]];
    const z = ζ.clone().scale(α).add(ζ.clone().conj().scale(β)).sub(b/3);
    if(z.Im < 0) z.conj();
    return (r <= z.Re ? [r, z, z.clone().conj()] : [z, z.clone().conj(), r])
        .map(x => [x, 1]);
}


/**
 * @summary
 * Find the roots of a quartic polynomial.
 *
 * @desc
 * This function finds the roots of the polynomial `x^4 + b x^3 + c x^2 + d x + e` using
 * [Descarte's method]{@link
 * https://en.wikipedia.org/wiki/Quartic_function#Descartes'_solution}.
 *
 * @example
 * descartes(-4, -13, 28, 60).factor();  // [[-2, 2], [3, 1], [5, 1]]
 *
 * @param {number} b - The cubic coefficient.
 * @param {number} c - The quadratic coefficient.
 * @param {number} d - The linear coefficient.
 * @param {number} e - The constant coefficient.
 * @param {number} [ε=1e-10] - The roots will be considered equal if a suitable
 *   discriminant is smaller than this.
 * @return {Root[]} The roots.  See {@link Polynomial#factor} for the ordering.
 *
 * @private
 */
function descartes(b, c, d, e, ε=1e-10) {
    // Change of variables x --> x-b/4
    const p = c - 3/8*b*b;
    const q = 1/8*b*b*b - 1/2*b*c + d;
    const r = -3/256*b*b*b*b + 1/16*c*b*b - 1/4*b*d + e;
    // Now we're solving x^4 + px^2 + qx + r

    const insertRoot = (roots, a, m) => {
        const i = roots.findIndex(
            ([x,]) => x instanceof Complex ? x.Re >= a : x >= a);
        if(i === -1)
            roots.push([a, m]);
        else
            roots.splice(i, 0, [a, m]);
        return roots.map(
            ([x, m]) => [x instanceof Complex ? x.sub(b/4) : x - b/4, m]);
    };

    // Handle roots at zero
    if(abs(r) <= ε) {
        if(abs(q) <= ε) {
            if(abs(p) <= ε)
                return [[-b/4, 4]];
            let s = sqrt(abs(p));
            if(p < 0)
                return [[-s-b/4, 1], [-b/4, 2], [s-b/4, 1]];
            return [[-b/4, 2], [C(-b/4, s), 1], [C(-b/4, -s), 1]];
        }
        return insertRoot(cardano(0, p, q, ε), 0, 1);
    }

    // Here r is nonzero.
    if(abs(q) <= ε) {
        // Biquadratic x^4 + px^2 + r.  Solve y^2 + py + r and take square
        // roots, handling multiplicities.
        let roots = quadratic(p, r, ε);
        if(roots.length === 1) {
            // Two double roots
            const r = roots[0][0];
            const s = sqrt(abs(r));
            return r >= 0
                ? [[  -b/4  -s,  2], [  -b/4 + s,  2]]
                : [[C(-b/4,  s), 2], [C(-b/4, -s), 2]];
        }
        const [[r1, ], [r2, ]] = roots;
        if(r1 instanceof Complex) {
            // Simple complex roots.
            let z = r1.pow(1/2);  // This has positive real and imaginary part.
            let w = z.clone().scale(-1);
            return [w.clone().conj(), w, z, z.clone().conj()]
                .map(x => [x.sub(b/4), 1]);
        }
        let s1 = sqrt(abs(r1));
        let s2 = sqrt(abs(r2));
        if(r2 < 0) // Simple complex roots
            return [C(-b/4, s2), C(-b/4, -s2), C(-b/4, s1), C(-b/4, -s1)]
                .map(x => [x, 1]);
        if(r1 < 0) // Two simple complex roots and two simple real roots
            return [-s2-b/4, C(-b/4, s1), C(-b/4, -s1), s2-b/4]
                .map(x => [x, 1]);
        // Four simple real roots
        return [-s2-b/4, -s1-b/4, s1-b/4, s2-b/4].map(x => [x, 1]);
    }

    // Here q and r are both nonzero.
    // Discriminant:
    const Δ = 16*r*p*p*p*p - 4*q*q*p*p*p - 128*r*r*p*p
          + 144*r*q*q*p - 27*q*q*q*q + 256*r*r*r;
    if(abs(Δ) <= ε) {
        // There is a multiple root.  It must be real because otherwise the
        // polynomial would be a square of a quadratic, which is biquadratic and
        // was handled above.  It does not have multiplicity 4 because then it
        // would be zero, which was also handled above.

        // Say the double root is at a, so we can factor as
        //             (x - a)^2 (x^2 + bx + c),
        // where b = 2a because there is no cubic term.
        const Δ0 = p*p + 12*r;
        if(abs(Δ0) <= ε) {
            // Triple real root a = sqrt(-p/6); the other is b = -3a
            let a = sqrt(max(-p/6, ε));
            if(q < 0) a *= -1;
            return a < 0
                ? [[   a-b/4, 3], [-3*a-b/4, 1]]
                : [[-3*a-b/4, 1], [   a-b/4, 3]];
        }

        // Here a is a root with exact multiplicity 2, and the other roots are
        // simple (otherwise the polynomial would be biquadratic).  Hence (x-a)
        // is the gcd of the polynomial and its derivative, so we can use
        // Euclid's algorithm to compute a.  The result is this:
        const a = -q*Δ0/(2*p*p*p - 8*r*p + 9*q*q);
        // If you write the denominator in terms of a and c then you find
        // that it can only be zero when Δ0 = 0 or q = 0.
        const c = p + 3*a*a;
        return insertRoot(quadratic(2*a, c, ε), a, 2);
    }

    // The discriminant is nonzero, so all roots are simple.
    // The resolvent cubic always has a positive real root by the intermediate
    // value theorem.  This retrieves that number:
    let [uu, ] = cardano(2*p, p*p-4*r, -q*q, ε).filter(
        ([x,]) => !(x instanceof Complex)).pop();
    uu = max(uu, ε);  // Paranoia
    const u = sqrt(uu);
    const s = -u;
    const t = (p + uu + q/u)/2;
    const v = t - q/u;
    // We have factored the polynomial into (x^2 + sx + t) (x^2 + ux + v).
    let roots1 = quadratic(s, t, ε);
    let roots2 = quadratic(u, v, ε);
    return sortRoots(roots1.concat(roots2), -b/4);
}


/**
 * @summary
 * Data type holding roots of a polynomial with multiplicity.
 *
 * @desc
 * Roots are counted with multiplicity: `[x, m]` indicates a root of a
 * polynomial at the (Complex or real) number `x` with multiplicity `m`.  A Root
 * can also be a simple number or Complex number, which is taken to have
 * multiplicity equal to 1.
 *
 * @typedef {(number|Complex|Array)} Root
 */

/**
 * @summary
 * Class representing a polynomial.
 *
 * @desc
 * A polynomial is stored as an Array of coefficients `[an, ..., a1, a0]` which
 * represents the polynomial `an x^n + ... + a1 x + a0` in the variable `x`.
 * Polynomials can be added and scaled like Vectors, but can also be multiplied
 * and evaluated at a number.
 *
 * Do not use the constructor directly; instead use the convenience routine
 * {@link Polynomial.create}.
 *
 * @example {@lang javascript}
 * Polynomial.create(1, 2, 1).toString(1); // "x^2 + 2.0 x + 1.0"
 *
 * @extends Array
 */
class Polynomial extends Array {
    /**
     * @summary
     * Create a Polynomial with the given coefficients.
     *
     * @desc
     * This is like `Polynomial.of`, except it strips all leading zeros from the
     * coefficients.
     *
     * @example {@lang javascript}
     * Polynomial.create(1, 2, 1).toString(1); // "x^2 + 2.0 x + 1.0"
     * Polynomial.create(1).toString(1);       // "1.0"
     * Polynomial.create(0, 1, 2).toString(1); // "x + 2.0"
     *
     * @param {...number} coeffs - The coefficients of the resulting
     *   Polynomial.
     * @return {Polynomial} The polynomial with the given coefficients.
     */
    static create(...coeffs) {
        while(coeffs[0] === 0)
            coeffs.shift();
        return Polynomial.from(coeffs);
    }

    /**
     * @summary
     * Create a monic Polynomial with the given roots.
     *
     * @desc
     * Given roots `λ1, λ2, ..., λn`, this returns the polynomial that factors
     * as `(x-λ1)(x-λ2)...(x-λn)`.  Roots are assumed to be provided in
     * complex-conjugate pairs; the complex part of the product is dropped.
     *
     * @example {@lang javascript}
     * Polynomial.fromRoots(1, -1).toString(1);    // x^2 - 1.0
     * Polynomial.fromRoots([1, 3]).toString(1);   // x^3 - 3.0 x^2 + 3.0 x - 1.0
     * Polynomial.fromRoots(1, Complex.i, Complex.i.conj()).toString(1);
     *    // x^3 - 1.0 x^2 + 1.0 x - 1.0
     *
     * @param {...Root} roots - The roots of the resulting polynomial.
     * @return {Polynomial} The monic polynomial with the given roots.
     */
    static fromRoots(...roots) {
        let P = [C(1)];
        let n = 1;
        roots = roots.flatMap(
            r => r instanceof Complex ? [r]
                : r instanceof Array ? Array.from(range(r[1]), () => C(r[0]))
                : [r]);
        while(roots.length > 0) {
            let c = roots.pop();
            P.forEach((a, i) => a.mult(c).scale(-1).add((i < n-1 ? P[i+1] : 0)));
            P.unshift(C(1));
            n++;
        }
        return Polynomial.from(P, x => x.Re);
    }

    /**
     * @summary
     * Compute the Legendre polynomial of degree `n`.
     *
     * @desc
     * See [the Wikipedia article]{@link
     * https://en.wikipedia.org/wiki/Legendre_polynomials}
     * for the definition and many wonderful properties of Legendre
     * polynomials.
     *
     * @example {@lang javascript}
     * Polynomial.legendre(5).toString(); // 7.8750 x^5 - 8.7500 x^3 + 1.8750 x
     *
     * @param {integer} n - Compute the Legendre polynomial of this degree.
     * @return {Polynomial} The `n`th Legendre polynomial.
     */
    static legendre(n) {
        if(legendres[n]) return legendres[n].clone();
        let p = Polynomial.legendre(n-1).clone().scale(2*n-1);
        p.push(0);
        let q = Polynomial.legendre(n-2).clone().scale(n-1);
        legendres[n] = p.sub(q).scale(1/n);
        return legendres[n].clone();
    }

    /**
     * @summary
     * The degree of the polynomial.
     *
     * @desc
     * This is one less than the number of coefficients, or `-Infinity` for the
     * zero polynomial.
     *
     * @example {@lang javascript}
     * Polynomial.create(1, 0, 0).deg;  // 2
     * Polynomia.lcreate().deg;         // -Infinity
     *
     * @type {integer}
     */
    get deg() {
        if(this.length > 0)
            return this.length - 1;
        return -Infinity;
    }

    /**
     * @summary
     * Create a new Polynomial with the same coefficients.
     *
     * @return {Polynomial} The new polynomial.
     */
    clone() {
        return this.slice();
    }

    /**
     * @summary
     * Return a string representation of the polynomial.
     *
     * @example {@lang javascript}
     * Polynomial.create(1, 2, 1).toString(1);       // "x^2 + 2.0 x + 1.0"
     * Polynomial.create(1, 0, -2).toString(1, 'z'); // "z^2 - 2.0"
     *
     * @param {integer} [precision=4] - The number of decimal places to include.
     * @param {string} [variable='x'] - The dummy variable.
     * @return {string} A string representation of the polynomial.
     */
    toString(precision=4, variable='x') {
        if(this.isZero()) return (0).toFixed(precision);
        let ret = '';
        for(let i = 0; i <= this.deg; ++i) {
            let n = this.deg - i;
            let c = this[i];
            if(c === 0) continue;
            if(i > 0) {
                if(c < 0) {
                    ret += ' - ';
                    c *= -1;
                } else
                    ret += ' + ';
            } else {
                if(c === -1) {
                    ret += '-';
                    c = 1;
                }
            }
            if(c !== 1 || n == 0) {
                ret += c.toFixed(precision);
                if(n > 0) ret += ' ';
            }
            if(n > 1) ret += `${variable}^${n}`;
            else if(n === 1) ret += variable;
        }
        return ret;
    }

    /**
     * @summary
     * Test if a polynomial is zero.
     *
     * @desc
     * This means that the degree is `-Infinity`, i.e. that the coefficient list
     * is empty.
     *
     * @example {@lang javascript}
     * Polynomial.create(1, 0, 0).isZero();  // false
     * Polynomial.create(0, 0, 0).isZero();  // true
     *
     * @return {boolean} True if the polynomial is zero.
     */
    isZero() {
        return this.length === 0;
    }

    /**
     * @summary
     * Test if this Polynomial is equal to `P`.
     *
     * @desc
     * Two polynomials are equal if they have the same degree, and all
     * coefficients are equal.
     *
     * @example {@lang javascript}
     * let p = Polynomial.create(1, 0.01, -0.01, 0);
     * let q = Polynomial.create(1, 0, 0, 0);
     * p.equals(q);                          // false
     * p.equals(q, 0.05);                    // true
     * q.equals(Polynomial.create(1, 0, 0)); // false
     *
     * @param {Polynomial} P - The polynomial to compare.
     * @param {number} [ε=0] - Coefficients will test as equal if they are
     *   within `ε` of each other.  This is provided in order to account for
     *   rounding errors.
     * @return {boolean} True if the polynomials are equal.
     */
    equals(P, ε=0) {
        if(this.deg !== P.deg)
            return false;
        if(ε === 0)
            return this.every((v, i) => v === P[i]);
        return this.every((v, i) => abs(v - P[i]) <= ε);
    }

    /**
     * @summary
     * Add a polynomial.
     *
     * @desc
     * This creates a new Polynomial equal to the sum of `this` and `P`.
     * Leading zero coefficients are stripped.
     *
     * @example {@lang javascript}
     * let p = Polynomial.create(1, 2), q = Polynomial.create(1, 3, 4);
     * p.add(q).toString(1);      // "x^2 + 4.0 x + 6.0"
     * p.add(q, 2).toString(1);   // "2.0 x^2 + 7.0 x + 10.0"
     * let q1 = Polynomial.create(-1, 0, 0);
     * q.add(q1).toString(1);     // "3.0 x + 4.0"
     *
     * @param {Polynomial} P - The polynomial to add.
     * @param {number} [factor=1] - Add `factor` times `P` instead of just
     *   adding `P`.
     * @return {Polynomial} A new polynomial equal to the sum.
     */
    add(P, factor=1) {
        let ret = this.clone();
        while(ret.deg < P.deg)
            ret.unshift(0);
        let i0 = ret.deg - P.deg;
        for(let i = 0; i < P.length; ++i)
            ret[i + i0] += factor*P[i];
        while(ret[0] === 0)
            ret.shift();
        return ret;
    }

    /**
     * @summary
     * Subtract a polynomial.
     *
     * @desc
     * This creates a new Polynomial equal to the difference of `this` and
     * `P`.  This is an alias for `this.add(P, -1)`.
     *
     * @example {@lang javascript}
     * let p = Polynomial.create(1, 2), q = Polynomial.create(1, 3, 4);
     * p.sub(q).toString(1);   // "-x^2 - 2.0 x - 2.0"
     * let q1 = Polynomial.create(1, 0, 0);
     * q.sub(q1).toString(1);  // "3.0 x + 4.0"
     *
     * @param {Polynomial} P - The polynomial to subtract.
     * @return {Polynomial} A new polynomial equal to the difference.
     */
    sub(P) {
        return this.add(P, -1);
    }

    /**
     * @summary
     * Multiply by another polynomial.
     *
     * @desc
     * This creates a new Polynomial equal to the product of `this` and `P`.
     * The degree of the resulting polynomial is the sum of the degrees.
     *
     * @example {@lang javascript}
     * let p = Polynomial.create(1, 2), q = Polynomial.create(1, 3, 4);
     * p.mult(q).toString(1);   // "x^3 + 5.0 x^2 + 10.0 x + 8.0"
     *
     * @param {Polynomial} P - The polynomial to multiply.
     * @return {Polynomial} A new polynomial equal to the product.
     */
    mult(P) {
        if(this.isZero() || P.isZero()) return Polynomial.create();
        return Polynomial.from(
            range(this.deg + P.deg + 1),
            i => {
                let acc = 0;
                for(let j = max(0, i-P.deg); j <= min(i, this.deg); ++j)
                    acc += this[j] * P[i-j];
                return acc;
            });
    }

    /**
     * @summary
     * Synthetic division by another polynomial.
     *
     * @desc
     * This returns the quotient and remainder obtained by dividing `this` by
     * `P`.  The quotient has degree equal to `this.deg - P.deg`, and the
     * remainder has degree less than `P.deg`.
     *
     * @example {@lang javascript}
     * let p = Polynomial.create(6, 5, 0, -7);
     * let q = Polynomial.create(3, -2, -1);
     * let [Q, R] = p.div(q);
     * Q.toString(1);  // 2.0 x + 3.0
     * R.toString(1);  // 8.0 x - 4.0
     *
     * @param {Polynomial} P - The polynomial to divide.
     * @param {number} [ε=0] - Leading coefficients of the remainder are
     *   considered to be zero if they are smaller than this.  This is provided
     *   in order to account for rounding errors.
     * @return {Polynomial[]} An array `[Q, R]`, where `Q` is the quotient of
     *   `this` by `P` and `R` is the remainder.
     * @throws Will throw an error if `P` has larger degree or is the zero
     *   polynomial.
     */
    div(P, ε=0) {
        if(P.deg > this.deg)
            throw new Error("Tried to divide by a polynomial of larger degree");
        if(P.isZero())
            throw new Error("Tried to divide by the zero polynomial");
        let ret = this.clone();
        let normalizer = P[0];
        for(let i = 0; i <= this.deg - P.deg; ++i) {
            ret[i] /= normalizer;
            let coef = ret[i];
            if(coef === 0) continue;
            for(let j = 1; j <= P.deg; ++j)
                ret[i+j] -= P[j] * coef;
        }
        let quot = ret.slice(0, this.deg - P.deg + 1);
        let rem = ret.slice(this.deg - P.deg + 1, ret.length);
        while(Math.abs(rem[0]) <= ε) rem.shift();
        return [quot, rem];
    }

    /**
     * @summary
     * Multiply in-place by a scalar.
     *
     * @desc
     * This multiplies all coefficients of `this` by the number `c`;
     *
     * @example {@lang javascript}
     * let p = Polynomial.create(1, 2, 3);
     * p.mult(2);
     * p.toString(1);   // "2.0 x^2 + 4.0 x + 6.0"
     *
     * @param {number} c - The scalar to multiply.
     * @return {Polynomial} `this`
     */
    scale(c) {
        if(c === 1) return this;
        for(let i = 0; i < this.length; ++i)
            this[i] *= c;
        return this;
    }

    /**
     * @summary
     * Scale in-place to be monic.
     *
     * @desc
     * This scales by the reciprocal of the highest-order coefficient.
     *
     * @example {@lang javascript}
     * let p = Polynomial.create(2, 4, 6);
     * p.monic();
     * p.toString(1);   // "x^2 + 2.0 x + 3.0"
     *
     * @return {Polynomial} `this`
     */
    monic() {
        return this.scale(1/this[0]);
    }

    /**
     * @summary
     * Raise to an integer power.
     *
     * @desc
     * This returns a new polynomial eaqual to the `n`th power of `this`.
     *
     * @example {@lang javascript}
     * Polynomial.create(1, 1).pow(4).toString(1);
     *    // x^4 + 4.0 x^3 + 6.0 x^2 + 4.0 x + 1.0
     *
     * @param {integer} n - The power to raise.
     * @return {Polynomial} A new Polynomial equal to the `n`th power of `this`.
     */
    pow(n) {
        let ret = Polynomial.create(1);
        let base = this;
        while(true) {
            if(n & 1)
                ret = ret.mult(base);
            n >>= 1;
            if(n === 0) break;
            base = base.mult(base);
        }
        return ret;
    }

    /**
     * @summary
     * Evaluate the polynomial at a number.
     *
     * @desc
     * This substitutes the dummy variable for the number `z`, which may be real
     * or Complex, and returns the result.
     *
     * @example {@lang javascript}
     * let p = Polynomial.create(1, 1, 1, 1);
     * p.evaluate(2);                              // 2**3 + 2**2 + 2 + 1
     * p.evaluate(new Complex(1, 1)).toString(1);  // 0.0 + 5.0 i
     *
     * @param {number|Complex} z - The number to evaluate.
     * @return {number|Complex} The result of evaluating `this` at `z`.
     */
    eval(z) {
        if(z instanceof Complex)
            return this.reduce((a, x) => a.mult(z).add(x), C(0, 0));
        return this.reduce((a, x) => a*z + x, 0);
    }

    /**
     * @summary
     * Compose with another polynomial.
     *
     * @desc
     * This substitutes the dummy variable for the polynomial `P`, and returns a
     * new Polynomial.  The degree of the resulting polynomial is the product of
     * the degrees.
     *
     * @example {@lang javascript}
     * let p = Polynomial.create(1, 1, 1, 1);
     * let q = Polynomial.create(1, 1, 1);
     * p.compose(q).toString(0);  // x^6 + 3 x^5 + 7 x^4 + 9 x^3 + 10 x^2 + 6 x + 4
     * q.compose(p).toString(0);  // x^6 + 2 x^5 + 3 x^4 + 5 x^3 + 4 x^2 + 3 x + 3
     *
     * @param {Polynomial} P - The polynomial to substitute.
     * @return {Polynomial} A new Polynomial equal to the result of composing
     *   `this` with `P`.
     */
    compose(P) {
        return this.reduce((a, x) => {
            a = a.mult(P);
            if(a.length == 0) {
                if(x !== 0) a.push(x);
            }
            else a[a.deg] += x;
            return a;
        }, Polynomial.create());
    }

    /**
     * @summary
     * Compute the `n`th formal derivative of `this`.
     *
     * @desc
     * The first derivative is the polynomial of degree one less, where the
     * coefficient of `x^i` is multiplied by `i`.  Higher derivatives are
     * computed recursively.
     *
     * @example {@lang javascript}
     * let p = Polynomial.create(1, 1, 1, 1);
     * p.derivative().toString(1);   // "3.0 x^2 + 2.0 x + 1.0"
     * p.derivative(2).toString(1);  // "6.0 x + 2.0"
     *
     * @param {integer} [n=1] - The number of derivatives to take.
     * @return {Polynomial} A new Polynomial obtained by taking `n` derivatives
     *   of `this`.
     */
    derivative(n=1) {
        if(n === 0) return this.clone();
        if(this.deg < n)
            return Polynomial.create();
        if(n === 1) {
            let d = this.deg;
            let ret = this.map((x, i) => (d-i) * x);
            ret.pop();
            return ret;
        }
        return this.derivative(n-1).derivative();
    }

    /**
     * @summary
     * Find all (real and complex) roots of the polynomial.
     *
     * @desc
     * This uses the quadratic formula in degree 2,
     * [Cardano's formula]{@link
     * https://www.encyclopediaofmath.org/index.php/Cardano_formula}
     * in degree 3, and
     * [Descarte's method]{@link
     * https://en.wikipedia.org/wiki/Quartic_function#Descartes'_solution}
     * in degree 4.  It is not implemented for higher degrees.
     *
     * Closed-form solutions like the quadratic formula do not exist in degrees
     * 5 and above: see the [Abel&ndash;Ruffini theorem]{@link
     * https://en.wikipedia.org/wiki/Abel-Ruffini_theorem}.  Of course, there
     * are numerical methods for finding roots of polynomials&mdash;the best of
     * which use linear algebraic methods to directly find the eigenvalues of a
     * matrix with a given characteristic polynomial&mdash;but they are beyond
     * the scope of this library.
     *
     * The return value is ordered as follows.
     *  * Roots with smaller real part come first.
     *  * Real roots come before complex roots with the same real part.
     *  * Complex roots come in adjacent conjugate pairs; the one with positive
     *    imaginary part comes first.
     *  * Pairs of complex roots with smaller imaginary part (in absolute value)
     *    come first.
     *
     * @example {@lang javascript}
     * Polynomial.fromRoots([-2, 2], 3, 5).factor();  // [[-2, 2], [3, 1], [5, 1]]
     *
     * @param {number} [ε=1e-10] - If certain discriminant quantities are
     *   smaller than this, then multiple roots have been found.
     * @return {Root[]} The roots found.
     * @throws Will throw an error if the degree is not 1, 2, 3, or 4.
     */
    factor(ε=1e-10) {
        if(this.isZero())
            throw new Error("The zero polynomial does not have discrete roots");
        const [a, b, c, d, e] = this;
        if(this.deg === 1)
            return [[-b/a, 1]];
        if(this.deg === 2)
            return quadratic(b/a, c/a, ε);
        if(this.deg === 3)
            return cardano(b/a, c/a, d/a, ε);
        if(this.deg === 4)
            return descartes(b/a, c/a, d/a, e/a, ε);
        throw new Error("Root finding is not implemented for degrees greater than 4");
    }
}

legendres = [Polynomial.create(1), Polynomial.create(1, 0)];

export default Polynomial;
