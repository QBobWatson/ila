/*
 * Copyright (c) 2020 Joseph Rabinoff
 * All rights reserved
 *
 * This file is part of linalg.js.
 *
 * linalg.js is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * linalg.js is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with linalg.js.  If not, see <https://www.gnu.org/licenses/>.
 */

/** @module vector
 *
 * @file
 * Implements a Vector class used for doing vector arithmetic.
 */

import { range } from "./util";
import Matrix from "./matrix3";


/**
 * @summary
 * Class representing a vector.
 *
 * @desc
 * A vector is a sequence of numbers of a fixed length, called the "size" of the
 * vector.
 *
 * @example
 * new Vector(1, 2, 3).toString(1);  // "[1.0 2.0 3.0]"
 *
 * @extends Array
 */
class Vector implements Iterable<number> {
    /**
     * @summary
     * The entries of the vector.
     *
     * @example {@lang javascript}
     * let v = new Vector(1, 2, 3);
     * v[0];  // 1
     * v[1];  // 2
     * v[2];  // 3
     */
    [index: number]: number;

    /**
     * @summary
     * The size of the vector.
     *
     * @desc
     * This is the number of entries in the vector.
     */
    size: number;

    /**
     * @summary
     * Create a Vector with the given entries using the given map function.
     *
     * @desc
     * This is a convenience method meant to mirror `Array.from`.
     *
     * @example {@lang javascript}
     * Vector.from([1,2,3], k => k*2).toString(0);  // "[2 4 6]"
     *
     * @param entries - The entries to pass through the map function.
     * @param mapfn - The map function to pass the entries through (by default
     *   the identity function)
     * @return The new Vector.
     */
    static from<T>(entries: Iterable<T> | ArrayLike<T>,
                   mapfn?: (x: T, k: number) => number): Vector {
        mapfn = mapfn || ((x, _) => (x as any) as number);
        return new Vector(...Array.from(entries, mapfn));
    }

    /**
     * @summary
     * Create a Vector with `n` entries equal to `c`.
     *
     * @example {@lang javascript}
     * Vector.constant(3, 1).toString(1);  // "[1.0 1.0 1.0]"
     *
     * @param n - The size of the resulting Vector, an integer.
     * @param c - The value of the entries.
     * @return The vector `[c, c, ..., c]`.
     */
    static constant(n: number, c: number): Vector {
        return Vector.from(range(n), () => c);
    }

    /**
     * @summary
     * Create a Vector with `n` entries equal to zero.
     *
     * @example {@lang javascript}
     * Vector.zero(3).toString(1);  // "[0.0 0.0 0.0]"
     *
     * @param n - The size of the resulting Vector, an integer.
     * @return The zero vector of size `n`.
     */
    static zero(n: number): Vector {
        return Vector.constant(n, 0);
    }

    /**
     * @summary
     * Create a unit coordinate vector.
     *
     * @desc
     * This is the Vector with all entries equal to 0, except the `i`th equal to
     * 1.
     *
     * @example {@lang javascript}
     * Vector.e(1, 3).toString(1); // "[0.0 1.0 0.0]"
     *
     * @param i - The nonzero entry, an integer.
     * @param n - The size of the resulting Vector, an integer.
     * @param [λ=1] - The nonzero entry is set to this.
     * @return The `i`th unit coordinate vector in `R^n`.
     */
    static e(i: number, n: number, λ=1): Vector {
        return Vector.from(range(n), j => j === i ? λ : 0);
    }

    /**
     * @summary
     * Check if an iterable of vectors is linearly independent.
     *
     * @desc
     * This means that the vector equation
     *   `c1 v1 + c2 v2 + ... + cn vn = 0`
     * has only the solution `c1 = c2 = ... = cn = 0`.  Equivalently, the matrix
     * with columns `v1, v2, ..., vn` has full column rank.
     *
     * @example {@lang javascript}
     * Vector.isLinearlyIndependent(
     *     [new Vector(1, 0, 0),
     *      new Vector(0, 1, 0),
     *      new Vector(0, 0, 1)]);  // true
     * Vector.isLinearlyIndependent(
     *     [new Vector(1, 0, 0),
     *      new Vector(0, 1, 0),
     *      new Vector(1, 1, 0)]);  // false
     *
     * @param vecs - The vectors to check.
     * @param [ε=1e-10] - Entries smaller than this value are taken to be zero
     *   for the purposes of pivoting.
     * @return True if the vectors are linearly independent.
     * @throws Error if the vectors do not have the same size.
     */
    static isLinearlyIndependent(vecs: Vector[], ε: number=1e-10): boolean {
        let M = Matrix.create(...vecs);
        // Tall matrices never have linearly independent rows.
        if(M.m > M.n) return false;
        return M.rank(M.PLU(ε)) == vecs.length;
    }

    /**
     * @summary
     * Check if an iterable of vectors is linearly dependent.
     *
     * @desc
     * This is an alias for `!Vector.isLinearlyIndependent(vecs)`.
     *
     * @example {@lang javascript}
     * Vector.isLinearlyDependent(
     *     [new Vector(1, 0, 0),
     *      new Vector(0, 1, 0),
     *      new Vector(0, 0, 1)]);  // false
     * Vector.isLinearlyDependent(
     *     [new Vector(1, 0, 0),
     *      new Vector(0, 1, 0),
     *      new Vector(1, 1, 0)]);  // true
     *
     * @param vecs - The vectors to check.
     * @param [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting.
     * @return True if the vectors are linearly dependent.
     * @throws Error if the vectors do not have the same size.
     */
    static isLinearlyDependent(vecs: Vector[], ε: number=1e-10): boolean {
        return !Vector.isLinearlyIndependent(vecs, ε);
    }

    /**
     * @summary
     * Return a linearly independent subset.
     *
     * @desc
     * This returns an Array containing a maximal linearly independent subset of
     * vectors from `vecs`.
     *
     * @example {@lang javascript}
     * let v1 = new Vector(1, 0, 0);
     * let v2 = new Vector(0, 1, 0);
     * let v3 = new Vector(1, 1, 0);
     * Vector.linearlyIndependentSubset([v1, v2, v3]);  // [v1, v2]
     *
     * @param {Vector[]} vecs - The vectors to use.
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting.
     * @return {Vector[]} A linearly independent subset of `vecs`.
     * @throws Will throw an error if the vectors do not have the same size.
     */
    static linearlyIndependentSubset(vecs: Vector[],
                                     ε: number=1e-10): Vector[] {
        let M = Matrix.create(...vecs).transpose;
        return Array.from(M.pivots(M.PLU(ε)), ([,j]) => vecs[j]);
    }

    /**
     * @summary
     * Compute a linear combination of vectors.
     *
     * @desc
     * The linear combination of the vectors `v1, v2, ..., vn` with coefficients
     * `c1, c2, ..., cn` is the vector `c1 v1 + c2 v2 + ... + cn vn`.
     *
     * @example {@lang javascript}
     * let v1 = new Vector(1, 0, 0);
     * let v2 = new Vector(0, 1, 0);
     * let v3 = new Vector(1, 1, 0);
     * Vector.linearCombination([1, 2, 3], [v1, v2, v3]).toString(1);
     *    // "[4.0 5.0 0.0]"
     *
     * @param coeffs - A non-empty array of coefficients.
     * @param vecs - A non-empty array of vectors, of the same size as `coeffs`.
     * @return The sum of the vectors scaled by the coefficients.
     * @throws Error if the vectors do not have the same size,
     *   or if `coeffs` is empty.
     */
    static linearCombination(coeffs: number[], vecs: Vector[]): Vector {
        return coeffs.reduce(
            (a, c, i) => a.add(vecs[i].clone().scale(c)),
            Vector.zero(vecs[0].size));
    }

    /**
     * @summary
     * The square of the geometric length of the vector.
     *
     * @desc
     * This is the sum of the squares of the entries, which is the dot product
     * of `this` with itself.
     *
     * @example {@lang javascript}
     * new Vector(3, 4).lensq;   // 25
     */
    get lensq(): number {
        return this.dot(this);
    }

    /**
     * @summary
     * The geometric length of the vector.
     *
     * @desc
     * This is the square root of `this.lensq`.
     *
     * @example {@lang javascript}
     * new Vector(3, 4).len;   // 5
     */
    get len(): number {
        return Math.sqrt(this.lensq);
    }


    /**
     * @summary
     * Create a Vector with the given entries.
     *
     * @example {@lang javascript}
     * new Vector(1).toString(1);    // "[1.0]"
     * new Vector(1, 2).toString(1); // "[1.0 2.0]"
     *
     * @param entries - The entries of the resulting Vector.
     */
    constructor(...entries: number[]) {
        entries.forEach((x, i) => this[i] = x);
        this.size = entries.length;
    }

    /**
     * @summary
     * Iterate over the entries.
     *
     * @example {@lang javascript}
     * let v = new Vector(1, 2, 1).toString(1);
     * [...v]; // [1, 2, 1]
     */
    [Symbol.iterator](): IterableIterator<number> {
        let self = this;
        let f = function*() {
            for (let i = 0; i < self.size; ++i)
                yield self[i];
        };
        return f();
    }

    /**
     * @summary
     * Check if this vector is equal to `other`.
     *
     * @desc
     * Two vectors are equal if they have the same number of entries, and all
     * entries are equal.
     *
     * @example {@lang javascript}
     * let v = new Vector(0.01, -0.01, 0);
     * let w = Vector.zero(3);
     * v.equals(w);              // false
     * v.equals(w, 0.05);        // true
     * w.equals(Vector.zero(2)); // false
     *
     * @param other - The vector to compare.
     * @param [ε=0] - Entries will test as equal if they are within `ε`
     *   of each other.  This is provided in order to account for rounding
     *   errors.
     * @return True if the vectors are equal.
     */
    equals(other: Vector, ε: number=0): boolean {
        if(this.size !== other.size)
            return false;
        if(ε === 0) {
            for(let i = 0; i < this.size; ++i)
                if(this[i] != other[i])
                    return false;
            return true;
        }
        for(let i = 0; i < this.size; ++i)
            if(Math.abs(this[i] - other[i]) > ε)
                return false;
        return true;
    }

    /**
     * @summary
     * Create a new Vector with the same entries.
     *
     * @return The new vector.
     */
    clone(): Vector {
        return new Vector(...this);
    }

    /**
     * @summary
     * Return a string representation of the vector.
     *
     * @example {@lang javascript}
     * new Vector(1, 2, 3).toString(2); // "[1.00 2.00 3.00]"
     *
     * @param [precision=4] - The number of decimal places to include.
     * @return A string representation of the vector.
     */
    toString(precision: number=4): string {
        const strings = Array.from(this, v => v.toFixed(precision));
        return "[" + strings.join(' ') + "]";
    }

    /**
     * @summary
     * Decide if a vector is zero.
     *
     * @desc
     * This is functionally equivalent to
     *    `this.equals(Vector.zero(this.size), ε)`.
     *
     * @param [ε=0] - Entries smaller than this in absolute value will
     *   be considered zero.
     * @return True if the vector has all zero entries.
     */
    isZero(ε: number=0): boolean {
        for(const x of this) {
            if(Math.abs(x) > ε)
                return false;
        }
        return true;
    }

    /**
     * @summary
     * Scale the vector by the reciprocal of its length.
     *
     * @desc
     * This modifies the vector in-place to have [length]{@link Vector#len} 1.
     *
     * @example {@lang javascript}
     * let v = new Vector(3, 4);
     * v.normalize();
     * v.toString(2); // "[0.60 0.80]"
     *
     * @return `this`
     * @throws Error if `this` is the zero vector.
     */
    normalize(): this {
        const s = 1/this.len;
        if(!isFinite(s))
            throw new Error("Tried to normalize the zero vector");
        return this.scale(s);
    }

    /**
     * @summary
     * Add a Vector in-place.
     *
     * @desc
     * This modifies the vector in-place by adding the entries of `other`.
     *
     * @example {@lang javascript}
     * let v = new Vector(1, 2), w = new Vector(3, 4);
     * v.add(w);
     * v.toString(1);  // "[4.0 6.0]"
     *
     * @param other - The vector to add.
     * @param [factor=1] - Add `factor` times `other` instead of just
     *   adding `other`.
     * @param [start=0] - Only add the entries `start...this.size`.
     *   Provided for optimizations when the entries of `other` before `start`
     *   are known to be zero.
     * @return `this`
     * @throws Error if the vectors have different sizes.
     */
    add(other: Vector, factor: number=1, start: number=0): this {
        if(this.size !== other.size)
            throw new Error(
                'Tried to add vectors of different sizes');
        for(let i = start; i < this.size; ++i)
            this[i] += factor*other[i];
        return this;
    }

    /**
     * @summary
     * Subtract a Vector in-place.
     *
     * @desc
     * This modifies the vector in-place by subtracting the entries of `other`.
     *
     * @example {@lang javascript}
     * let v = new Vector(1, 2), w = new Vector(3, 4);
     * v.sub(w);
     * v.toString(1);  // "[-2.0 -2.0]"
     *
     * @param other - The vector to subtract.
     * @param [start=0] - Only subtract the entries `start...this.size`.
     *   Provided for optimizations when the entries of `other` before `start`
     *   are known to be zero.
     * @return `this`
     * @throws Error if the vectors have different sizes.
     */
    sub(other: Vector, start: number=0): this {
        return this.add(other, -1, start);
    }

    /**
     * @summary
     * Multiply a Vector by a scalar in-place.
     *
     * @desc
     * This modifies the vector in-place by multiplying all entries by `c`.
     *
     * @example {@lang javascript}
     * let v = new Vector(1, 2);
     * v.scale(2);
     * v.toString(1);  // "[2.0 4.0]"
     *
     * @param c - The scaling factor.
     * @param [start=0] - Only scale the entries
     *   `start...this.size`.  Provided for optimizations when the entries
     *   before `start` are known to be zero.
     * @return `this`
     */
    scale(c: number, start: number=0): this {
        for(let i = start; i < this.size; ++i)
            this[i] *= c;
        return this;
    }

    /**
     * @summary
     * Compute the dot product with another vector.
     *
     * @desc
     * This is the sum of the pairwise products of the entries of `this` and
     * `other`.
     *
     * @example {@lang javascript}
     * let v = new Vector(1, 2), w = new Vector(3, 4);
     * v.dot(w);  // 1*3 + 2*4
     *
     * @param other - The vector to dot.
     * @return The dot product.
     * @throws Error if the vectors have different sizes.
     */
    dot(other: Vector): number {
        if(this.size !== other.size)
            throw new Error(
                'Tried to take the dot product of vectors of different sizes');
        let a = 0;
        for(let i = 0; i < this.size; ++i)
            a += this[i] * other[i];
        return a;
    }
};


export default Vector;
