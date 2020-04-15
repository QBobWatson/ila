/** @module lib/polynomial
 *
 * @file
 * A class for polynomial manipulations.
 */

'use strict'; // -*- js2 -*-

import Complex from './complex.js';

// Convenience
const C = (a, b=0) => a instanceof Complex ? a : new Complex(a, b);

// Memoized Legendre polynomials
let legendres;

/**
 * Roots are counted with multiplicity: `[x, m]` indicates a root of a
 * polynomial at the (Complex or real) number `x` with multiplicity `m`.  A Root
 * can also be a simple number or Complex number, which is taken to have
 * multiplicity equal to 1.
 *
 * @typedef {(number|Complex|Array)} Root
 */

/**
 * Class representing a polynomial.
 *
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
     * Create a Polynomial with the given coefficients.
     *
     * Unlike `new Polynomial()`, this works when `coeffs` has length
     * one.  This strips all leading zeros from the coefficients.
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
        if(coeffs.length === 1)
            return new Polynomial(1).fill(coeffs[0]);
        return new Polynomial(...coeffs);
    }

    /**
     * Create a monic Polynomial with the given roots.
     *
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
     * @param {...Root} roots - The roots of the resulting Polynomial.
     * @return {Polynomial} The monic polynomial with the given roots.
     */
    static fromRoots(...roots) {
        let P = [C(1)];
        let n = 1;
        roots = roots.flatMap(r => {
            if(r instanceof Complex) return [r];
            if(r instanceof Array) {
                let [c, m] = r;
                return new Array(m).fill(0).map(() => C(c));
            }
            return [r];
        });
        while(roots.length > 0) {
            let c = roots.pop();
            P.forEach((a, i) => a.mult(c).scale(-1).add((i < n-1 ? P[i+1] : 0)));
            P.unshift(C(1));
            n++;
        }
        return Polynomial.create(...P.map(x => x.Re));
    }

    /**
     * Compute the Legendre polynomial of degree `n`.
     *
     * See [the Wikipedia article]{@link https://en.wikipedia.org/wiki/Legendre_polynomials}.
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
     * The degree of the polynomial.
     *
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
     * Create a new Polynomial with the same coefficients.
     *
     * @return {Polynomial} The new polynomial.
     */
    clone() {
        return this.slice();
    }

    /**
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
     * Decide if a polynomial is zero.
     *
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
     * Check if this Polynomial is equal to `P`.
     *
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
        return this.every((v, i) => Math.abs(v - P[i]) < ε);
    }

    /**
     * Add a polynomial.
     *
     * This creates a new Polynomial equal to the sum of `this` and `P`.
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
     * Subtract a polynomial.
     *
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
     * Multiply by another polynomial.
     *
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
        let ret = new Polynomial(this.deg + P.deg + 1).fill(0);
        for(let i = 0; i < ret.length; ++i) {
            for(let j = Math.max(0, i-P.deg); j <= Math.min(i, this.deg); ++j)
                ret[i] += this[j] * P[i-j];
        }
        return ret;
    }

    /**
     * Synthetic division by another polynomial.
     *
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
     * @return {Polynomial[]} An array `[Q, R]`, where `Q` is the quotient of
     *   `this` by `P` and `R` is the remainder.
     * @throws Will throw an error if `P` has larger degree or is the zero
     *   polynomial.
     */
    div(P) {
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
        while(rem[0] === 0) rem.shift();
        return [quot, rem];
    }

    /**
     * Multiply in-place by a scalar.
     *
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
        for(let i = 0; i < this.length; ++i)
            this[i] *= c;
        return this;
    }

    /**
     * Raise to an integer power.
     *
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
            if(!n) break;
            base = base.mult(base);
        }
        return ret;
    }

    /**
     * Evaluate the polynomial at a number.
     *
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
     * Compose with another polynomial.
     *
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
     * Compute the `n`th formal derivative of `this`.
     *
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
}

legendres = [Polynomial.create(1), Polynomial.create(1, 0)];

export default Polynomial;
