/** @module lib/vector
 *
 * @file
 * Implements a Vector class used for doing vector arithmetic.
 */

'use strict';

import Matrix from "./matrix.js";

var Vector;

export function initVector() {
    if(Vector) return;
    Vector =

/**
 * Class representing a vector.
 *
 * A vector is a sequence of numbers of a fixed length, called the "length" of
 * the vector.
 *
 * Note that `Vector` uses the `Array` constructor, so `new Vector(3)` will
 * return an "empty" vector of length 3.  Instead use the convenience routine
 * {@link Vector.create}.
 *
 * @example
 * Vector.create(1, 2, 3).toString(1);  // "[1.0 2.0 3.0]"
 *
 * @class Vector
 * @extends Array
 */
class Vector extends Array {
    /**
     * Create a Vector with the given entries.
     *
     * Unlike `new Vector()`, this works when `entries` has length one.
     *
     * @example {@lang javascript}
     * Vector.create(1).toString(1);    // "[1.0]"
     * Vector.create(1, 2).toString(1); // "[1.0 2.0]"
     *
     * @param {...number} entries - The entries of the resulting Vector.
     * @return {Vector} The vector with the given entries.
     */
    static create(...entries) {
        if(entries.length === 1)
            return Vector.constant(1, entries[0]);
        return new Vector(...entries);
    }

    /**
     * Create a Vector with all entries equal to `c`.
     *
     * @example {@lang javascript}
     * Vector.constant(3, 1).toString(1);  // "[1.0 1.0 1.0]"
     *
     * @param {integer} n - The size of the resulting Vector.
     * @param {number} c - The value of the entries.
     * @return {Vector} The vector `[c, c, ..., c]`.
     */
    static constant(n, c) {
        let ret = new Vector(n);
        ret.fill(c);
        return ret;
    }

    /**
     * Create a Vector with all entries equal to zero.
     *
     * @example {@lang javascript}
     * Vector.zero(3).toString(1);  // "[0.0 0.0 0.0]"
     *
     * @param {integer} n - The size of the resulting Vector.
     * @return {Vector} The zero vector of size `n`.
     */
    static zero(n) {
        return Vector.constant(n, 0);
    }

    /**
     * Create a unit coordinate vector.
     *
     * This is the Vector with all entries equal to 0, except the `i`th equal to
     * 1.
     *
     * @example {@lang javascript}
     * Vector.e(1, 3).toString(1); // "[0.0 1.0 0.0]"
     *
     * @param {integer} i - The nonzero entry.
     * @param {integer} n - The size of the resulting Vector.
     * @param {number} [λ=1] - The nonzero entry is set to this.
     * @return {Vector} The `i`th unit coordinate vector in `R^n`.
     */
    static e(i, n, λ=1) {
        let ret = Vector.zero(n);
        ret[i] = λ;
        return ret;
    }

    /**
     * Check if an iterable of vectors is linearly independent.
     *
     * This means that the vector equation
     *   `c1 v1 + c2 v2 + ... + cn vn = 0`
     * has only the solution `c1 = c2 = ... = cn = 0`.
     *
     * @example {@lang javascript}
     * Vector.isLinearlyIndependent(
     *     [Vector.create(1, 0, 0),
     *      Vector.create(0, 1, 0),
     *      Vector.create(0, 0, 1)]);  // true
     * Vector.isLinearlyIndependent(
     *     [Vector.create(1, 0, 0),
     *      Vector.create(0, 1, 0),
     *      Vector.create(1, 1, 0)]);  // false
     *
     * @param {Vector[]} vecs - The vectors to check.
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting.
     * @return {boolean} True if the vectors are linearly independent.
     * @throws Will throw an error if the vectors do not have the same length.
     */
    static isLinearlyIndependent(vecs, ε=1e-10) {
        let M = new Matrix(...vecs);
        // Tall matrices never have linearly independent rows.
        if(M.m > M.n) return false;
        return M.rank(ε) == vecs.length;
    }

    /**
     * Check if an iterable of vectors is linearly dependent.
     *
     * This is an alias for `!Vector.isLinearlyIndependent(vecs)`.
     *
     * @example {@lang javascript}
     * Vector.isLinearlyDependent(
     *     [Vector.create(1, 0, 0),
     *      Vector.create(0, 1, 0),
     *      Vector.create(0, 0, 1)]);  // false
     * Vector.isLinearlyDependent(
     *     [Vector.create(1, 0, 0),
     *      Vector.create(0, 1, 0),
     *      Vector.create(1, 1, 0)]);  // true
     *
     * @param {Vector[]} vecs - The vectors to check.
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting.
     * @return {boolean} True if the vectors are linearly dependent.
     * @throws Will throw an error if the vectors do not have the same length.
     */
    static isLinearlyDependent(vecs, ε=1e-10) {
        return !Vector.isLinearlyIndependent(vecs, ε);
    }

    /**
     * Return a linearly independent subset.
     *
     * This returns an Array containing a maximal linearly independent subset of
     * vectors from `vecs`.
     *
     * @example {@lang javascript}
     * let v1 = Vector.create(1, 0, 0);
     * let v2 = Vector.create(0, 1, 0);
     * let v3 = Vector.create(1, 1, 0);
     * Vector.linearlyIndependentSubset([v1, v2, v3]);  // [v1, v2]
     *
     * @param {Vector[]} vecs - The vectors to use.
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting.
     * @return {Vector[]} A linearly independent subset of `vecs`.
     * @throws Will throw an error if the vectors do not have the same length.
     */
    static linearlyIndependentSubset(vecs, ε=1e-10) {
        let M = new Matrix(...vecs).transpose;
        let pivots = M.pivots(ε);
        return pivots.map(([, j]) => vecs[j]);
    }

    /**
     * Compute a linear combination of vectors.
     *
     * The linear combination of the vectors `v1, v2, ..., vn` with coefficients
     * `c1, c2, ..., cn` is the vector `c1 v1 + c2 v2 + ... + cn vn`.
     *
     * @example {@lang javascript}
     * let v1 = Vector.create(1, 0, 0);
     * let v2 = Vector.create(0, 1, 0);
     * let v3 = Vector.create(1, 1, 0);
     * Vector.linearCombination([1, 2, 3], [v1, v2, v3]).toString(1);
     *    // "[4.0 5.0 0.0]"
     *
     * @param {number[]} coeffs - A non-empty array of coefficients.
     * @param {Vector[]} vecs - A non-empty array of vectors, of the same length
     *   as `coeffs`.
     * @return The sum of the vectors scaled by the coefficients.
     * @throws Will throw an error if the vectors do not have the same length,
     *   or if `coeffs` is empty.
     */
    static linearCombination(coeffs, vecs) {
        return coeffs.reduce((a, c, i) => a.add(vecs[i].clone().scale(c)),
                             Vector.zero(vecs[0].length));
    }


    /**
     * The square of the geometric length of the vector.
     *
     * @example {@lang javascript}
     * Vector.create(3, 4).sizesq;   // 25
     *
     * @type {number}
     */
    get sizesq() {
        return this.dot(this);
    }

    /**
     * The geometric length of the vector.
     *
     * This is the square root of `this.sizesq`.
     *
     * @example {@lang javascript}
     * Vector.create(3, 4).size;   // 5
     *
     * @type {number}
     */
    get size() {
        return Math.sqrt(this.sizesq);
    }

    /**
     * Check if this vector is equal to `other`.
     *
     * Two vectors are equal if they have the same number of entries, and all
     * entries are equal.
     *
     * @example {@lang javascript}
     * let v = Vector.create(0.01, -0.01, 0);
     * let w = Vector.zero(3);
     * v.equals(w);              // false
     * v.equals(w, 0.05);        // true
     * w.equals(Vector.zero(2)); // false
     *
     * @param {Vector} other - The vector to compare.
     * @param {number} [ε=0] - Entries will test as equal if they are within `ε`
     *   of each other.  This is provided in order to account for rounding
     *   errors.
     * @return {boolean} True if the vectors are equal.
     */
    equals(other, ε=0) {
        if(this.length !== other.length)
            return false;
        if(ε === 0)
            return this.every((v, i) => v === other[i]);
        return this.every((v, i) => Math.abs(v - other[i]) < ε);
    }

    /**
     * Create a new Vector with the same entries.
     *
     * @return {Vector} The new vector.
     */
    clone() {
        return this.slice();
    }

    /**
     * Return a string representation of the vector.
     *
     * @example {@lang javascript}
     * Vector.create(1, 2, 3).toString(2); // "[1.00 2.00 3.00]"
     *
     * @param {integer} [precision=4] - The number of decimal places to include.
     * @return {string} A string representation of the vector.
     */
    toString(precision=4) {
        // this.map() returns a new Vector
        const strings = [...this.map(v => v.toFixed(precision))];
        return "[" + strings.join(' ') + "]";
    }

    /**
     * Decide if a vector is zero.
     *
     * This is functionally equivalent to
     *    `this.equals(Vector.zero(this.length), ε)`.
     *
     * @param {number} [ε=0] - Entries smaller than this in absolute value will
     *   be considered zero.
     * @return {boolean} True if the vector has all zero entries.
     */
    isZero(ε=0) {
        return this.every(x => Math.abs(x) <= ε);
    }

    /**
     * Scale the vector by the reciprocal of its length.
     *
     * This modifies the vector in-place to have length 1.
     *
     * @example {@lang javascript}
     * let v = Vector.create(3, 4);
     * v.normalize();
     * v.toString(2); // "[0.60 0.80]"
     *
     * @return {Vector} `this`
     * @throws Will throw an error if `this` is the zero vector.
     */
    normalize() {
        const s = 1/this.size;
        if(!isFinite(s))
            throw new Error("Tried to normalize the zero vector");
        return this.scale(s);
    }

    /**
     * Add a Vector in-place.
     *
     * This modifies the vector in-place by adding the entries of `other`.
     *
     * @example {@lang javascript}
     * let v = Vector.create(1, 2), w = Vector.create(3, 4);
     * v.add(w);
     * v.toString(1);  // "[4.0 6.0]"
     *
     * @param {Vector} other - The vector to add.
     * @param {number} [factor=1] - Add `factor` times `other` instead of just
     *   adding `other`.
     * @param {integer} [start=0] - Only add the entries `start...this.length`.
     *   Provided for optimizations when the entries of `other` before `start`
     *   are known to be zero.
     * @return {Vector} `this`
     * @throws Will throw an error if the vectors have different lengths.
     */
    add(other, factor=1, start=0) {
        if(this.length !== other.length)
            throw new Error(
                'Tried to add vectors of different lengths');
        for(let i = start; i < this.length; ++i)
            this[i] += factor*other[i];
        return this;
    }

    /**
     * Subtract a Vector in-place.
     *
     * This modifies the vector in-place by subtracting the entries of `other`.
     *
     * @example {@lang javascript}
     * let v = Vector.create(1, 2), w = Vector.create(3, 4);
     * v.sub(w);
     * v.toString(1);  // "[-2.0 -2.0]"
     *
     * @param {Vector} other - The vector to subtract.
     * @param {integer} [start=0] - Only subtract the entries
     *   `start...this.length`.  Provided for optimizations when the entries of
     *   `other` before `start` are known to be zero.
     * @return {Vector} `this`
     * @throws Will throw an error if the vectors have different lengths.
     */
    sub(other, start=0) {
        return this.add(other, -1, start);
    }

    /**
     * Multiply a Vector by a scalar in-place.
     *
     * This modifies the vector in-place by multiplying all entries by `c`.
     *
     * @example {@lang javascript}
     * let v = Vector.create(1, 2);
     * v.scale(2);
     * v.toString(1);  // "[2.0 4.0]"
     *
     * @param {number} c - The scaling factor.
     * @param {integer} [start=0] - Only scale the entries
     *   `start...this.length`.  Provided for optimizations when the entries
     *   before `start` are known to be zero.
     * @return {Vector} `this`
     */
    scale(c, start=0) {
        for(let i = start; i < this.length; ++i)
            this[i] *= c;
        return this;
    }

    /**
     * Compute the dot product with another vector.
     *
     * This is the sum of the pairwise products of the entries of `this` and
     * `other`.
     *
     * @example {@lang javascript}
     * let v = Vector.create(1, 2), w = Vector.create(3, 4);
     * v.dot(w);  // 1*3 + 2*4
     *
     * @param {Vector} other - The vector to dot.
     * @return {number} The dot product.
     * @throws Will throw an error if the vectors have different lengths.
     */
    dot(other) {
        if(this.length !== other.length)
            throw new Error(
                'Tried to take the dot product of vectors of different lengths');
        return this.reduce((a, v, i) => a + v * other[i], 0);
    }
};

}

initVector();

export {Vector as default};
