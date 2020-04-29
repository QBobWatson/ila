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

'use strict';

/** @module subspace
 *
 * @file
 * Implements a Subspace class to represent subspaces of `R^n`.
 */

// TODO: principal angles?

import Matrix from "./matrix.js";
import { range } from "./util.js";


class Subspace {
    /**
     * @summary
     * The subspace `R^n`.
     *
     * @desc
     * This creates the subspace of `R^n` with
     * [basis]{@link Subspace#basis} equal to the `n`x`n` identity
     * matrix.
     *
     * @example {@lang javascript}
     * Subspace.Rn(3).toString();  // "The full subspace R^3"
     *
     * @param {integer} n - The dimension of the space.
     * @return {Subspace} The subspace `R^n`.
     */
    static Rn(n) {
        return new Subspace(Matrix.identity(n), {n: n, isBasis: true, isON: true});
    }

    /**
     * @summary
     * The zero subspace of `R^n`.
     *
     * @desc
     * This creates the subspace of `R^n` with no generators.
     *
     * @example {@lang javascript}
     * Subspace.zero(3).toString();  // "The zero subspace of R^3"
     *
     * @param {integer} n - The ambient dimension.
     * @return {Subspace} The subspace `{0}`.
     */
    static zero(n) {
        return new Subspace([], {n: n, isBasis: true, isON: true});
    }

    /**
     * @summary
     * Class representing a Subspace.
     *
     * @desc
     * A subspace of `R^n` is a collection of vectors `V` satisfying:
     *  * The zero vector is contained in `V`.
     *  * The sum of two vectors `V` is also in `V`.
     *  * All scalar multiples of a vector in `V` is also in `V`.
     *
     * Subspaces often arise as a span of some number of vectors, or as the null
     * space of a matrix.  Any subspace can be expressed as a span, which is the
     * same as the column space of a matrix.  This class represents a subspace
     * as the column space of the matrix {@link Subspace#basis}.
     *
     * This constructor takes a spanning set for the Subspace.  (To produce a
     * null space, use {@link Matrix#nullSpace}).  It constructs a
     * [basis]{@link Subspace#basis} from a maximal linearly independent subset
     * of the generators.
     *
     * @example {@lang javascript}
     * let A = Matrix.create([ 0, -3, -6,  4,  9],
     *                       [-1, -2, -1,  3,  1],
     *                       [-2, -3,  0,  3, -1],
     *                       [ 1,  4,  5, -9, -7]);
     * let V = new Subspace(A);
     * V.toString(1);
     *   // "Subspace of R^4 of dimension 3 with basis
     *   //  [ 0.0] [-3.0] [ 4.0]
     *   //  [-1.0] [-2.0] [ 3.0]
     *   //  [-2.0] [-3.0] [ 3.0]
     *   //  [ 1.0] [ 4.0] [-9.0]"
     *
     * @param {(Vector[]|Matrix)} generators - A spanning set for the subspace.
     *   If it is a matrix, use the columns of the matrix as the spanning set.
     * @param {Object} [hints={}] - Precomputed properties of the subspace.
     * @param {integer} [hints.n=generators.m] - The dimension of the ambient
     *   R^n.  Only necessary to provide if `generators` is empty.
     * @param {boolean} [hints.isBasis=false] - The generators are linearly
     *   independent.
     * @param {boolean} [hints.isON=false] - The generators are orthonormal.
     * @param {number} [hints.ε=1e-10] - Use this value in {@link
     *   Matrix#colBasis} when computing a basis.
     * @throws Will throw an error if the generators do not all have the same
     *   length, or if it can't determine the ambient space R^n.
     * @see Matrix#colBasis
     */
    constructor(generators, {n, isBasis=false, isON=false, ε=1e-10}={}) {
        if(n !== undefined)
            /**
             * @summary
             * The dimension of the ambient `R^n`.
             *
             * @desc
             * All vectors in this subspace have `n` entries.
             *
             * @example {@lang javascript}
             * new Subspace([[1, 2, 3], [4, 5, 6]]).n;  // 3
             *
             * @type {integer}
             */
            this.n = n;
        else if(generators instanceof Matrix)
            this.n = generators.m;
        else if(generators.length > 0)
            this.n = generators[0].length;
        else
            throw new Error("Cannot determine the ambient R^n from zero generators");

        generators = generators instanceof Matrix
            ? generators : Matrix.create(...generators).transpose;

        /**
         * @summary
         * The columns of this matrix form a basis for this subspace.
         *
         * @desc
         * A nonzero subspace has infinitely many bases.  This is the basis
         * created in the constructor from the passed generators.
         *
         * @example {@lang javascript}
         * let A = Matrix.create([ 0, -3, -6,  4,  9],
         *                       [-1, -2, -1,  3,  1],
         *                       [-2, -3,  0,  3, -1],
         *                       [ 1,  4,  5, -9, -7]);
         * let V = new Subspace(A);
         * V.basis.toString(1);
         *   // "[ 0.0 -3.0  4.0]
         *   //  [-1.0 -2.0  3.0]
         *   //  [-2.0 -3.0  3.0]
         *   //  [ 1.0  4.0 -9.0]"
         *
         * @type {Matrix}
         */
        this.basis = isBasis
            ? generators
            : Matrix.from(generators.colBasis(ε)).transpose;

        this._ONbasis = isON ? generators : undefined;
    }

    /**
     * @summary
     * The dimension of the subspace.
     *
     * @desc
     * The dimension of a subspace is by definition the number of vectors in any
     * basis of the subspace.  Hence this is an alias for `this.basis.n`.
     *
     * @example {@lang javascript}
     * let A = Matrix.create([ 0, -3, -6,  4,  9],
     *                       [-1, -2, -1,  3,  1],
     *                       [-2, -3,  0,  3, -1],
     *                       [ 1,  4,  5, -9, -7]);
     * let V = new Subspace(A);
     * V.dim;  // 3
     *
     * @type {integer}
     */
    get dim() {
        return this.basis.n;
    }

    /**
     * @summary
     * Return a string representation of the subspace.
     *
     * @example {@lang javascript}
     * let A = Matrix.create([ 0, -3, -6,  4,  9],
     *                       [-1, -2, -1,  3,  1],
     *                       [-2, -3,  0,  3, -1],
     *                       [ 1,  4,  5, -9, -7]);
     * let V = new Subspace(A);
     * V.toString(1);
     *   // "Subspace of R^4 of dimension 3 with basis
     *   //  [ 0.0] [-3.0] [ 4.0]
     *   //  [-1.0] [-2.0] [ 3.0]
     *   //  [-2.0] [-3.0] [ 3.0]
     *   //  [ 1.0] [ 4.0] [-9.0]"
     *
     * @param {integer} [precision=4] - The number of decimal places to include.
     * @return {string} A string representation of the subspace.
     */
    toString(precision=4) {
        if(this.dim === this.n)
            return `The full subspace R^${this.n}`;
        if(this.dim === 0)
            return `The zero subspace of R^${this.n}`;
        let ret = `Subspace of R^${this.n} of dimension ${this.dim} with basis\n`;
        let rowStrings = Array.from(
            this.basis.rows(), row => Array.from(row, x => x.toFixed(precision)));
        let colLengths = Array.from(
            range(this.dim), i => Math.max(...rowStrings.map(row => row[i].length)));
        return ret + rowStrings.map(
            row => row.map((x, i) => `[${x.padStart(colLengths[i])}]`).join(' '))
            .join('\n');
     }

    /**
     * @summary
     * Test whether the subspace is all of `R^n`.
     *
     * @desc
     * The only `n`-dimensional subspace of `R^n` is `R^n` itself, so this is a
     * shortcut for `this.dim === this.n`.
     *
     * @example {@lang javascript}
     * Subspace.Rn(3).isMaximal();  // true
     *
     * @example {@lang javascript}
     * let A = Matrix.create([ 0, -3, -6,  4,  9],
     *                       [-1, -2, -1,  3,  1],
     *                       [-2, -3,  0,  3, -1],
     *                       [ 1,  4,  5, -9, -7]);
     * new Subspace(A).isMaximal();  // false
     *
     * @return {boolean} True if this subspace is equal to `R^n`.
     */
    isMaximal() {
        return this.dim === this.n;
    }

    /**
     * @summary
     * Test whether the subspace only contains the zero vector.
     *
     * @desc
     * The only zero-dimensional subspace of `R^n` is `{0}`, so this is a
     * shortcut for `this.dim === 0`.
     *
     * @example {@lang javascript}
     * Subspace.zero(3).isZero();  // true
     *
     * @example {@lang javascript}
     * let A = Matrix.create([ 0, -3, -6,  4,  9],
     *                       [-1, -2, -1,  3,  1],
     *                       [-2, -3,  0,  3, -1],
     *                       [ 1,  4,  5, -9, -7]);
     * new Subspace(A).isZero();  // false
     *
     * @return {boolean} True if this subspace is equal to `{0}`.
     */
    isZero() {
        return this.dim === 0;
    }

    /**
     * @summary
     * Return the sum of two subspaces.
     *
     * @desc
     * The sum is the subspace generated by the bases of `this` and `other`.  It
     * is the smallest subspace containing both `this` and `other`.
     *
     * @example {@lang javascript}
     * let V = new Subspace([[ 0, -1, -2, 1], [-3, -2, -3,  4]]);
     * let W = new Subspace([[-6, -1,  0, 5], [ 4,  3,  3, -9]]);
     * V.add(W).toString(1);
     *   // "Subspace of R^4 of dimension 3 spanned by
     *   //  [ 0.0] [-3.0] [ 4.0]
     *   //  [-1.0] [-2.0] [ 3.0]
     *   //  [-2.0] [-3.0] [ 3.0]
     *   //  [ 1.0] [ 4.0] [-9.0]"
     *
     * @param {Subspace} other - The subspace to add.
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting.
     * @return {Subspace} The sum of `this` and `other`.
     * @throws Will throw an error of `this.n != other.n`.
     */
    add(other, ε=1e-10) {
        if(this.n != other.n)
            throw new Error("Tried to add subspaces of different R^n");
        return new Subspace(
            [...this.basis.cols(), ...other.basis.cols()],
            {n: this.n, ε});
    }

    /**
     * @summary
     * Take the intersection of two subspaces.
     *
     * @desc
     * This is the subspace consisting of all vectors contained both in `this`
     * and in `other`.  It is computed by taking the orthogonal complement of
     * the sum of the orthogonal complements.
     *
     * @example {@lang javascript}
     * let V = new Subspace([[ 0, -1, -2, 1], [-3, -2, -3,  4]]);
     * let W = new Subspace([[-6, -1,  0, 5], [ 4,  3,  3, -9]]);
     * V.intersect(W).toString();
     *   // "Subspace of R^4 of dimension 1 spanned by
     *   //  [ 1.0000]
     *   //  [ 0.1667]
     *   //  [ 0.0000]
     *   //  [-0.8333]"
     *
     * @param {Subspace} other - The subspace to intersect.
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting.
     * @return {Subspace} The largest subspace contained in `this` and `other`.
     * @throws Will throw and error of `this.n != other.n`.
     */
    intersect(other, ε=1e-10) {
        if(this.n != other.n)
            throw new Error("Tried to add subspaces of different R^n");
        return this.perp(ε).add(other.perp(ε)).perp(ε);
    }

    /**
     * @summary
     * Compute the projection matrix onto `this`.
     *
     * @desc
     * The projection matrix is the `n`x`n` matrix `P` such that `Pv` is the
     * projection of `v` onto this subspace.  This method computes `P` using the
     * formula `A(A^TA)^(-1)A^T`, where `A` is the matrix `this.basis`, unless
     * an orthonormal basis has been computed, in which case it returns `QQ^T`,
     * where the columns of `Q` form an orthonormal basis.
     *
     * @example {@lang javascript}
     * let V = new Subspace(Matrix.create([ 0, -3, -6,  4,  9],
     *                                    [-1, -2, -1,  3,  1],
     *                                    [-2, -3,  0,  3, -1],
     *                                    [ 1,  4,  5, -9, -7]));
     * let P = V.projectionMatrix();
     * P.toString();
     *   // "[1.0000  0.0000  0.0000  0.0000]
     *   //  [0.0000  0.1667  0.3333 -0.1667]
     *   //  [0.0000  0.3333  0.8667  0.0667]
     *   //  [0.0000 -0.1667  0.0667  0.9667]"
     * P.mult(P).equals(P, 1e-10);      // true
     * P.colSpace().equals(V);          // true
     * P.nullSpace().equals(V.perp());  // true
     *
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting.
     * @return {Matrix} The `n`x`n` matrix for projection onto `this`.
     * @see Subspace#basis
     * @see Subspace#ONbasis
     */
    projectionMatrix(ε=1e-10) {
        if(this.isZero())    this._projectionMatrix = Matrix.zero(this.n);
        if(this.isMaximal()) this._projectionMatrix = Matrix.identity(this.n);
        if(this._projectionMatrix) return this._projectionMatrix;
        if(this._ONbasis) {
            let Q = this._ONbasis;
            this._projectionMatrix = Q.mult(Q.transpose);
        } else {
            let A = this.basis;
            this._projectionMatrix = A.mult(A.normal.inverse(ε)).mult(A.transpose);
        }
        return this._projectionMatrix;
    }

    /**
     * @summary
     * Compute the orthogonal projection of a vector onto `this`.
     *
     * @desc
     * The orthogonal projection `v1` of `v` is the closest vector to `v` in the
     * subspace.  It is defined by the property that `v-v1` is orthogonal to
     * `V`.
     *
     * If the [projection matrix]{@link Subspace#projectionMatrix} `P` has been
     * cached, `v1` is computed as `Pv`.  Otherwise `v1` is computed using
     * `basis.projectColSpace()`.  If you want to compute many orthogonal
     * projections, run `this.projectionMatrix()` first.
     *
     * @example {@lang javascript}
     * let V = new Subspace(Matrix.create([ 0, -3, -6,  4,  9],
     *                                    [-1, -2, -1,  3,  1],
     *                                    [-2, -3,  0,  3, -1],
     *                                    [ 1,  4,  5, -9, -7]));
     *
     * let u = Vector.create(1, 2, 3, 4);
     * let u1 = V.project(u);
     * u1.toString();                 // "[1.0000 0.6667 3.5333 3.7333]"
     * V.contains(u1);                // true
     * V.isOrthogonalTo(u.sub(u1));   // true
     *
     * let v = Vector.create(0, -1, -2, 1); // in V
     * V.project(v).toString(0);      // "[0 -1 -2 1]"
     *
     * let w = Vector.create(0, 5, -2, 1); // in Vperp
     * V.project(w).toString(0);      // "[0 0 0 0]"
     *
     * @param {Vector} v - The vector to project.
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting.
     * @return {Vector} The orthogonal projection.
     * @throws Throws an error if `v.length != this.n`.
     * @see Subspace#orthoDecomp
     * @see Subspace#projectionMatrix
     * @see Matrix#projectColSpace
     */
    project(v, ε=1e-10) {
        if(this._projectionMatrix)
            return this._projectionMatrix.apply(v);
        return this.basis.projectColSpace(v, ε);
    }

    /**
     * @summary
     * Compute the orthogonal decomposition of a vector with respect to `this`.
     *
     * @desc
     * This returns the unique pair of vectors `[v1, v2]` such that `v1` is in
     * `this`, `v2` is orthogonal to `this`, and `v1 + v2 = v`.  The vector `v1`
     * is the [orthogonal projection]{@link Subspace#project} of `v`, and `v2`
     * is just `v-v1`.
     *
     * @example {@lang javascript}
     * let V = new Subspace(Matrix.create([ 0, -3, -6,  4,  9],
     *                                    [-1, -2, -1,  3,  1],
     *                                    [-2, -3,  0,  3, -1],
     *                                    [ 1,  4,  5, -9, -7]));
     *
     * let u = Vector.create(1, 2, 3, 4);
     * let [u1, u2] = V.orthoDecomp(u);
     * u1.toString();                 // "[1.0000 0.6667 3.5333 3.7333]"
     * u2.toString();                 // "[0.0000 1.3333 -0.5333 0.2667]"
     * u1.clone().add(u2).equals(u);  // true
     * V.contains(u1);                // true
     * V.isOrthogonalTo(u2);          // true
     *
     * let v = Vector.create(0, -1, -2, 1); // in V
     * let [v1, v2] = V.orthoDecomp(v);
     * v1.toString(0);                // "[0 -1 -2 1]"
     * v2.toString(0);                // "[0 0 0 0]"
     *
     * let w = Vector.create(0, 5, -2, 1); // in Vperp
     * let [w1, w2] = V.orthoDecomp(w);
     * w1.toString(0);                // "[0 0 0 0]"
     * w2.toString(0);                // "[0 5 -2 1]"
     *
     * @param {Vector} v - The vector to decompose.
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting.
     * @return {Vector[]} Returns `[v1, v2]` as described above.
     * @throws Throws an error if `v.length != this.n`.
     * @see Subspace#project
     */
    orthoDecomp(v, ε=1e-10) {
        let v1 = this.project(v, ε);
        return [v1, v.clone().sub(v1)];
    }

    /**
     * @summary
     * This is an alias for `this.orthoDecomp(v, ε)[1]`.
     *
     * @desc
     * The complement `v2` of `v` with respect to `this` is the shortest vector
     * from `this` to `v`.  It is defined by the property that `v-v2` is the
     * orthogonal projection of `v` onto `this`.  Equivalently, `v2` is the
     * orthogonal projection of `v` onto the orthogonal complement.
     *
     * @example {@lang javascript}
     * let V = new Subspace(Matrix.create([ 0, -3, -6,  4,  9],
     *                                    [-1, -2, -1,  3,  1],
     *                                    [-2, -3,  0,  3, -1],
     *                                    [ 1,  4,  5, -9, -7]));
     *
     * let u = Vector.create(1, 2, 3, 4);
     * let u2 = V.complement(u);
     * u2.toString();                 // "[0.0000 1.3333 -0.5333 0.2667]"
     * V.isOrthogonalTo(u2);          // true
     * V.contains(u2.sub(u));         // true
     *
     * let v = Vector.create(0, -1, -2, 1); // in V
     * let v2 = V.complement(v);
     * v2.toString(0);                // "[0 0 0 0]"
     *
     * let w = Vector.create(0, 5, -2, 1); // in Vperp
     * let w2 = V.complement(w);
     * w2.toString(0);                // "[0 5 -2 1]"
     *
     * @param {Vector} v - The vector to project.
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting.
     * @return {Vector} The orthogonal projection onto the orthogonal
     *   complement.
     * @throws Throws an error if `v.length != this.n`.
     * @see Subspace#project
     * @see Subspace#orthoDecomp
     * @see Subspace#perp
     */
    complement(v, ε=1e-10) {
        return this.orthoDecomp(v, ε)[1];
    }

    /**
     * @summary
     * Compute the distance of `v` from `this`.
     *
     * @desc
     * This is an alias for `this.complement(v, ε).size`.
     *
     * @param {Vector} v - The vector to measure.
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting.
     * @return {number} The distance to `this`.
     * @throws Throws an error if `v.length != this.n`.
     */
    distanceTo(v, ε=1e-10) {
        return this.complement(v, ε).size;
    }

    /**
     * @summary
     * Test if this subspace contains a vector.
     *
     * @desc
     * To say that `this` contains `v` means that `v` is a linear combination of
     * the [basis vectors]{@link Subspace#basis}.
     *
     * If `v.size <= ε` then this method returns true.  Otherwise this method
     * measures the [distance]{@link Subspace#distanceTo} from `v.normalize()`
     * to `this`, and returns `true` if the distance is at most `ε`.
     *
     * @example {@lang javascript}
     * let V = new Subspace(Matrix.create([ 0, -3, -6,  4,  9],
     *                                    [-1, -2, -1,  3,  1],
     *                                    [-2, -3,  0,  3, -1],
     *                                    [ 1,  4,  5, -9, -7]));
     * let v = Vector.create(-9, -4, -5, 10);  // Sum of the first three columns
     * V.contains(v);    // true
     * let w = Vector.create(-9, -4, -5, 9);
     * V.contains(w);    // false
     * V.distanceTo(w);  // 0.1825  (approximately)
     *
     * @param {Vector} v - The vector to test.
     * @param {number} [ε=1e-10] - Vectors shorter than this value are taken to
     *   be zero.  Also used for pivoting if no projection matrix has been
     *   cached.
     * @return {boolean} True if the distance from `v` to `this` is at most `ε`.
     * @throws Will throw an error if `v.length != this.n`.
     * @see Subspace#distanceTo
     * @see Subspace#project
     */
    contains(v, ε=1e-10) {
        if(this.n != v.length)
            throw new Error("Vector has the wrong number of entries");
        if(this.isMaximal())
            return true; // Nothing to check!
        if(this.isZero())
            return v.sizesq <= ε*ε;
        return v.sizesq <= ε*ε || this.distanceTo(v.clone().normalize(), ε) <= ε;
    }

    /**
     * @summary
     * Test if a vector is orthogonal to `this`.
     *
     * @desc
     * To say that `v` is orthogonal to `this` means that the dot product of `v`
     * with the [basis vectors]{@link Subspace#basis} is equal to zero.
     *
     * If `v.size <= ε` then this method returns true.  Otherwise this method
     * measures the size of the [projection]{@link Subspace#project} of
     * `v.normalize()` onto `this`, and returns `true` if the size is at most
     * `ε`.  This is (mathematically if not numerically) equivalent to
     * `this.perp().contains(v, ε)`.
     *
     * @example {@lang javascript}
     * let V = new Subspace(Matrix.create([ 0, -3, -6,  4,  9],
     *                                    [-1, -2, -1,  3,  1],
     *                                    [-2, -3,  0,  3, -1],
     *                                    [ 1,  4,  5, -9, -7]));
     *
     * let v = Vector.create(0, 5, -2, 1);
     * V.isOrthogonalTo(v);                  // true
     * V.basis.transpose.apply(v).isZero();  // true
     * let w = Vector.create(0, 5, -2, 2);
     * V.isOrthogonalTo(w);                  // false
     * V.project(w).size;                    // 0.9831  (approximately)
     *
     * @param {Vector} v - The vector to test.
     * @param {number} [ε=1e-10] - Vectors shorter than this value are taken to
     *   be zero.  Also used for pivoting if no projection matrix has been
     *   cached.
     * @return {boolean} True if the distance from `v` to the orthogonal
     *   complement of `this` is at most `ε`.
     * @throws Will throw an error if `v.length != this.n`.
     * @see Subspace#project
     * @see Subspace#distanceTo
     * @see Subspace#perp
     */
    isOrthogonalTo(v, ε=1e-10) {
        if(this.n != v.length)
            throw new Error("Vector has the wrong number of entries");
        if(this.isMaximal())
            return v.isZero(ε); // Nothing to check!
        if(this.isZero())
            return true;
        return v.sizesq <= ε*ε || this.project(v.clone().normalize(), ε).sizesq <= ε*ε;
    }

    /**
     * @summary
     * Test if this subspace is contained in another subspace.
     *
     * @desc
     * This subspace `V` is contained in another subspace `W` when all of the
     * [basis vectors]{@link Subspace#basis} of `V` are [contained in]{@link
     * Subspace#contains} `W`.
     *
     * @example {@lang javascript}
     * let V = new Subspace([[ 0, -3, -6,  4,  9],
     *                       [-1, -2, -1,  3,  1],
     *                       [-2, -3,  0,  3, -1],
     *                       [ 1,  4,  5, -9, -7]]);
     * let W = new Subspace([[ 1,  0, -3,  0,  5],
     *                       [ 0,  0,  0,  1,  0]]);
     * W.isSubspaceOf(V);       // true
     * W.perp().isSubspaceOf(V);  // false
     *
     * @param {Subspace} other - The subspace to test.
     * @param {number} [ε=1e-10] - Parameter passed to {@link Subspace#contains}.
     * @return {boolean} True if `this` is contained in `other`.
     * @throws Will throw an error if `this.n != other.n`.
     * @see Subspace#basis
     * @see Subspace#contains
     */
    isSubspaceOf(other, ε=1e-10) {
        if(this.n != other.n)
            throw new Error("Tried to test containment of subspaces in different R^n");
        if(this.dim > other.dim)
            return false;
        for(const col of this.basis.cols()) {
            if(!other.contains(col, ε))
                return false;
        }
        return true;
    }

    /**
     * @summary
     * Test if this Subspace is equal to `other`.
     *
     * @desc
     * Two subspaces are equal if and only if they have the same dimension and
     * one contains the other.
     *
     * @example {@lang javascript}
     * let V = new Subspace([[ 0, -3, -6,  4,  9],
     *                       [-1, -2, -1,  3,  1],
     *                       [-2, -3,  0,  3, -1],
     *                       [ 1,  4,  5, -9, -7]]);
     * let W = new Subspace([[-1, -5, -7,  7, 10],
     *                       [-3, -5, -1,  6,  0],
     *                       [-1,  1,  5, -6, -8]]);
     * let U = new Subspace([[ 1,  0,  0,  0,  0],
     *                       [ 0,  1,  0,  0,  0],
     *                       [ 0,  0,  1,  0,  0]]);
     * V.equals(W);  // true
     * V.equals(U);  // false
     *
     * @param {Subspace} other - The subspace to compare.
     * @param {number} [ε=1e-10] - Parameter passed to {@link Subspace#contains}.
     * @return {boolean} True if the subspaces are equal.
     * @throws Will throw an error if `this.n != other.n`.
     * @see Subspace#dim
     * @see Subspace#isSubspaceOf
     */
    equals(other, ε=1e-10) {
        if(this.n != other.n)
            throw new Error("Tried to test equality of subspaces in different R^n");
        if(this.dim !== other.dim)
            return false;
        return this.isSubspaceOf(other, ε);
    }

    /**
     * @summary
     * Return the orthogonal complement of `this`.
     *
     * @desc
     * The orthogonal complement of a subspace `V` is the set of all vectors `w`
     * such that `w` is orthogonal to every vector in `V`.  Its dimension is `n`
     * minus the dimension of `V`.
     *
     * If `V` is the column space of a matrix `A`, then its orthogonal
     * complement is the left null space of `A`.
     *
     * @example {@lang javascript}
     * let V = new Subspace([[ 0, -3, -6,  4,  9],
     *                       [-1, -2, -1,  3,  1],
     *                       [-2, -3,  0,  3, -1],
     *                       [ 1,  4,  5, -9, -7]]);
     * V.toString(1);
     *   // "Subspace of R^5 of dimension 3 with basis
     *   //  [ 0.0] [-1.0] [-2.0]
     *   //  [-3.0] [-2.0] [-3.0]
     *   //  [-6.0] [-1.0] [ 0.0]
     *   //  [ 4.0] [ 3.0] [ 3.0]
     *   //  [ 9.0] [ 1.0] [-1.0]"
     * V.perp().toString(1);
     *   // "Subspace of R^5 of dimension 2 with basis
     *   //  [ 0.0] [ 1.0]
     *   //  [-0.2] [-0.6]
     *   //  [ 1.0] [ 0.0]
     *   //  [ 0.0] [-0.0]
     *   //  [ 0.6] [-0.2]"
     * Array.from(V.perp().basis.cols()).every(col => V.isOrthogonalTo(col));
     *   // true
     *
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting.
     * @return {Subspace} The orthogonal complement.
     * @see Subspace#isOrthogonalTo
     * @see Matrix#leftNullSpace
     */
    perp(ε=1e-10) {
        if(this.isZero())
            return Subspace.Rn(this.n);
        if(this.isMaximal())
            return Subspace.zero(this.n);
        return this.basis.leftNullSpace(ε);
    }

    /**
     * @summary
     * Compute an orthonormal basis for the subspace.
     *
     * @desc
     * This is a shortcut for `this.basis.QR(ε).Q`.
     *
     * @example {@lang javascript}
     * let V = new Subspace([[ 3,  1, -1,  3],
     *                       [-5,  1,  5, -7],
     *                       [ 1,  1, -2,  8]]);
     * let Q = V.ONbasis();
     * Q.toString();
     *   // "[ 0.6708  0.2236 -0.6708]
     *   //  [ 0.2236  0.6708  0.2236]
     *   //  [-0.2236  0.6708  0.2236]
     *   //  [ 0.6708 -0.2236  0.6708]"
     * Q.colSpace().equals(V);  // true
     * Q.transpose.mult(Q).toString(1);
     *   // "[1.0 0.0 0.0]
     *   //  [0.0 1.0 0.0]
     *   //  [0.0 0.0 1.0]"
     *
     * @param {number} [ε=1e-10] - Vectors shorter than this value are taken to
     *   be zero.
     * @return {Matrix} A matrix whose columns form an orthonormal basis for the
     *   Subspace.
     * @see Matrix#QR
     */
    ONbasis(ε=1e-10) {
        if(this._ONbasis) return this._ONbasis;
        this._ONbasis = this.basis.QR(ε).Q;
        return this._ONbasis;
    }
}

export default Subspace;
