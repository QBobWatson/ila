'use strict';

/** @module matrix
 *
 * @file
 * Implements a Matrix class containing algorithms from basic linear algebra.
 */

import Vector from "./vector.js";
import Complex from "./complex.js";
import Subspace from "./subspace.js";
import Polynomial from "./polynomial.js";

import { range } from "./util.js";

// TODO:
//  * PCA?

/**
 * @summary
 * Class representing a matrix.
 *
 * @desc
 * A matrix is a 2-dimensional array of numbers.  The first dimension indexes
 * the row, written horizontally; the second indexes the column.
 *
 * The rows of a `Matrix` are {@link Vector} instances.  All rows must have the
 * same length.
 *
 * Use the static method {@link Matrix.create} instead of the constructor.
 *
 * @example {@lang javascript}
 * Matrix.create([1, 2], [3, 4], [5, 6]).toString(0);
 *   // "[1 2]
 *   //  [3 4]
 *   //  [5 6]"
 *
 * @extends Array
 */
class Matrix extends Array {
    /**
     * @summary
     * Create a Matrix.
     *
     * @desc
     * The arguments are {@link Vector} instances or Arrays of numbers, used as
     * the rows of the matrix.  All rows must have the same length.
     *
     * @example {@lang javascript}
     * Matrix.create([1, 2], [3, 4], [5, 6]).toString(0);
     *   // "[1 2]
     *   //  [3 4]
     *   //  [5 6]"
     *
     * @param {...(Array<number>|Vector)} rows - The rows of the matrix.  Arrays
     *   will be promoted to Vector instances.  All rows must have the same
     *   length.
     * @return {Matrix} The new matrix.
     * @throws Will throw an error if the rows do not have the same length.
     */
    static create(...rows) {
        if(rows.length === 0)
            return new Matrix();
        let ret = Matrix.from(
            rows, row => row instanceof Vector ? row : Vector.from(row));
        let n = ret[0].length;
        if(ret.some(row => row.length !== n))
            throw new Error("Matrix rows must have the same length.");
        return ret;
    }

    /**
     * @summary
     * Create an `n`x`n` identity Matrix.
     *
     * @desc
     * This is the square matrix with ones on the diagonal and zeros elsewhere.
     *
     * @example {@lang javascript}
     * Matrix.identity(3).toString(0);
     *   // "[1 0 0]
     *   //  [0 1 0]
     *   //  [0 0 1]"
     *
     * @param {integer} n - The resulting Matrix will have this many rows and columns.
     * @param {number} [λ=1] - The diagonal entries will be equal to `λ`.
     * @return {Matrix} The `n`x`n` (scaled) identity matrix.
     */
    static identity(n, λ=1) {
        return Matrix.from(range(n), i => Vector.e(i, n, λ));
    }

    /**
     * @summary
     * Create an `m`x`n` zero Matrix.
     *
     * @example {@lang javascript}
     * Matrix.zero(2, 3).toString(0);
     *   // "[0 0 0]
     *   //  [0 0 0]"
     *
     * @param {integer} m - The resulting Matrix will have this many rows.
     * @param {integer} [n=m] - The resulting Matrix will have this many columns.
     * @return {Matrix} The `m`x`n` zero matrix.
     */
    static zero(m, n=m) {
        return Matrix.from(range(m), i => Vector.zero(n));
    }

    /**
     * @summary
     * Create a diagonal matrix with specified diagonal entries.
     *
     * @desc
     * A diagonal matrix has nonzero entries only along the main diagonal.  This
     * method creates a `d`x`d` diagonal matrix with diagonal entries equal to
     * `entries`, where `d = entries.length`.  If `m` is specified, then it
     * creates an `m`x`d` matrix, and if `m` and `n` are both specified, then it
     * creates an `m`x`n` matrix.
     *
     * @example {@lang javascript}
     * Matrix.diagonal([1, 2, 3]).toString(0);
     *  // "[1 0 0]
     *  //  [0 2 0]
     *  //  [0 0 3]"
     * Matrix.diagonal([1, 2, 3], 4).toString(0);
     *  // "[1 0 0]
     *  //  [0 2 0]
     *  //  [0 0 3]
     *  //  [0 0 0]"
     * Matrix.diagonal([1, 2, 3], 3, 4).toString(0);
     *  // "[1 0 0 0]
     *  //  [0 2 0 0]
     *  //  [0 0 3 0]"
     *
     * @param {number[]} entries - The diagonal entries of the resulting Matrix.
     * @param {integer} [m=entries.length] - The resulting Matrix will have this
     *   many rows.
     * @param {integer} [n=entries.length] - The resulting Matrix will have this
     *   many columns.
     * @return {Matrix} The `m`x`n` diagonal matrix with `entries` along the
     *   main diagonal.
     */
    static diagonal(entries, m=entries.length, n=entries.length) {
        let ret = Matrix.zero(m, n);
        for(let i = 0; i < Math.min(m, n, entries.length); ++i)
            ret[i][i] = entries[i];
        return ret;
    }

    /**
     * @summary
     * Create a permutation Matrix.
     *
     * @desc
     * This is the square matrix whose `i`th row is the `vals[i]`th unit
     * coordinate vector `Vector.e(vals[i], n)`.
     *
     * @example {@lang javascript}
     * Matrix.permutation([1, 0]).toString(0);
     *   // "[0 1]
     *   //  [1 0]"
     *
     * @example {@lang javascript}
     * Matrix.permutation([2, 0, 1]).toString(0);
     *   // "[0 0 1]
     *   //  [1 0 0]
     *   //  [0 1 0]"
     *
     * @param {number[]} vals - If `n` numbers are given, these should be a
     *   permutation of the set `{0,1,...,n-1}`.
     * @return {Matrix} The permutation matrix `M` with the property that the
     *   `i`th entry of `M.apply(v)` is the `vals[i]`th entry of `v`.
     */
    static permutation(vals) {
        let n = vals.length;
        return Matrix.from(range(n), i => Vector.e(vals[i], n));
    }

    /**
     * @private
     *
     * @type {object}
     * @property {Matrix} transpose - Transpose matrix.
     * @property {integer} rank - The rank of the matrix.
     * @property {PLUData} PLU - PLU factorization; computed in `PLU()`.
     * @property {Matrix} rref - Reduced row echelon form.
     * @property {Matrix} E - Matrix such that `E*this = rref`.
     * @property {Matrix} adjugate - Transpose matrix of cofactors.
     * @property {Matrix} normal - The normal matrix `A^TA`.
     * @property {Vector[]} nullBasis - Basis for the null space.
     * @property {Subspace} nullSpace - The null space.
     * @property {Vector[]} colBasis - Basis for the column space.
     * @property {Subspace} nullSpace - The column space.
     * @property {Vector[]} rowBasis - Basis for the row space.
     * @property {Subspace} rowSpace - The row space.
     * @property {Vector[]} leftNullBasis - Basis for the left null space.
     * @property {Subspace} leftNullSpace - The left null space.
     * @property {QRData} QR - QR factorization; computed in `QR()`.
     * @property {Polynomial} charpoly - Characteristic polynomial.
     * @property {Root[]} eigenvalues - Eigenvalues.
     * @property {Map.<number, Subspace>} eigenspaces - Real eigenspaces.
     * @property {Map.<Complex, Array.<Complex[]>>} cplxEigenspaces - Complex
     *   eigenspaces.
     * @property {Matrix} pinv - The pseudo-inverse.
     * @property {boolean} isSymmetric - Whether the matrix is symmetric.
     * @property {boolean} isUpperTri - Whether the matrix is upper-triangular.
     * @property {boolean} isUpperUni - Whether the matrix is upper-unitriangular.
     * @property {boolean} isLowerTri - Whether the matrix is lower-triangular.
     * @property {boolean} isLowerUni - Whether the matrix is lower-unitriangular.
     * @property {boolean} isEchelon - Whether the matrix is in row-echelon form.
     * @property {boolean} hasONCols - Whether the matrix has orthonormal columns.
     */
    get _cache() {
        if(!this.__cache) this.__cache = {};
        return this.__cache;
    }

    /**
     * @summary
     * The number of rows of the matrix.
     *
     * @example {@lang javascript}
     * Matrix.create([1, 2], [3, 4], [5, 6]).m;  // 3
     *
     * @type {integer}
     */
    get m() { return this.length; }

    /**
     * @summary
     * The number of columns of the matrix.
     *
     * @example {@lang javascript}
     * Matrix.create([1, 2], [3, 4], [5, 6]).n;  // 2
     *
     * @type {integer}
     */
    get n() { return this.length === 0 ? 0 : this[0].length; }

    /**
     * @summary
     * The transpose Matrix.
     *
     * @desc
     * The rows of the transpose are the columns of the matrix.
     *
     * @example {@lang javascript}
     * Matrix.create([1, 2], [3, 4], [5, 6]).transpose.toString(0);
     *   // "[1 3 5]
     *   //  [2 4 6]"
     *
     * @type {Matrix}
     */
    get transpose() {
        if(this._cache.transpose) return this._cache.transpose;
        this._cache.transpose = Matrix.from(this.cols());
        return this._cache.transpose;
    }

    /**
     * @summary
     * The normal matrix `A^TA`.
     *
     * @desc
     * The normal matrix is a symmetric matrix with the same null space.  It is
     * used in least-squares computations.
     *
     * @example {@lang javascript}
     * Matrix.create([1, 2], [3, 4], [5, 6]).normal.toString(0);
     *   // "[35 44]
     *   //  [44 56]"
     *
     * @type {Matrix}
     */
    get normal() {
        if(this._cache.normal) return this._cache.normal;
        this._cache.normal = this.transpose.mult(this);
        this._cache.normal.hint({isSymmetric: true});
        return this._cache.normal;
    }

    /**
     * @summary
     * The sum of the diagonal entries of the matrix.
     *
     * @example {@lang javascript}
     * Matrix.create([1, 2],
     *               [3, 4],
     *               [5, 6]).trace;  // 1 + 4
     *
     * @type {number}
     */
    get trace() {
        let acc = 0;
        for(const d of this.diag()) acc += d;
        return acc;
    }

    _fadeev_leverrier() {
        let n = this.n;
        let ret = Polynomial.create(1);
        let AM = Matrix.zero(n), adjugate;
        let c = 1;
        for(let k = 1; k <= n; ++k) {
            for(let i = 0; i < n; ++i)
                AM[i][i] += c;
            if(k == n) adjugate = AM;
            AM = this.mult(AM);
            c = -AM.trace/k;
            ret.push(c);
        }
        // This is det(λI - A); multiply by (-1)^n now
        if(n % 2 === 1) ret.scale(-1);
        if(n % 2 === 0) adjugate.scale(-1);
        this._cache.charpoly = ret;
        this._cache.adjugate = adjugate;
    }

    /**
     * @summary
     * The characteristic polynomial of the matrix.
     *
     * @desc
     * The characteristic polynomial of an `n`x`n` matrix `A` is the determinant
     * of `(A - λI_n)`, where λ is an indeterminate.  This is a polynomial of
     * degree `n` and leading coefficient `(-1)^n`:
     * > `(-1)^n λ^n + c1 λ^(n-1) + c2 λ^(n-2) + ... + c(n-1) λ + cn`
     *
     * This only makes sense for square matrices.  Throws an error if the matrix
     * is not square.
     *
     * This method implements the
     * [Fadeev&ndash;LeVerrier algorithm]{@link
     * https://en.wikipedia.org/wiki/Faddeev-LeVerrier_algorithm}.
     * This is not the most efficient algorithm: it runs in approximately
     * `O(n^4)` time.  However, it has the advantage that it computes the
     * [adjugate matrix]{@link Matrix#adjugate} at the same time.
     *
     * @example {@lang javascript}
     * Matrix.create([1, 6,4],
     *               [2,-1,3],
     *               [5, 0,1]).charpoly.toString(0);  // "-x^3 + x^2 + 33 x + 97"
     *
     * @type {Polynomial}
     * @throws Will throw an error if the matrix is not square.
     */
    get charpoly() {
        if(!this.isSquare())
            throw new Error("Tried to compute the characteristic polynomial of a non-square matrix");
        if(this._cache.charpoly) return this._cache.charpoly;
        this._fadeev_leverrier();
        return this._cache.charpoly;
    }

    /**
     * @summary
     * The adjugate matrix of the matrix.
     *
     * @desc
     * The adjugate matrix of `A` is the matrix `B` whose `(i, j)` entry is the
     * `(j, i)` cofactor of `A`.  It has the property that
     * > `AB = BA = det(A) I`
     *
     * In particular, if `A` is invertible then `B = det(A) A^(-1)`.
     *
     * This method uses the
     * [Fadeev&ndash;LeVerrier algorithm]{@link
     * https://en.wikipedia.org/wiki/Faddeev-LeVerrier_algorithm},
     * which computes the [characteristic polynomial]{@link Matrix#charpoly} at
     * the same time.  It runs in about `O(n^4)` time, which is much better than
     * computing `n^2` determinants.
     *
     * This only makes sense for square matrices.  Throws an error if the matrix
     * is not square.
     *
     * @example {@lang javascript}
     * let A = Matrix.create([1, 6,4],
     *                       [2,-1,3],
     *                       [5, 0,1]);
     * let B = A.adjugate;
     * B.toString(0);
     *   // "[-1  -6  22]
     *   //  [13 -19   5]
     *   //  [ 5  30 -13]"
     * A.mult(B).toString(0);
     *   // "[97  0  0]
     *   //  [ 0 97  0]
     *   //  [ 0  0 97]"
     * B.mult(A).toString(0);
     *   // "[97  0  0]
     *   //  [ 0 97  0]
     *   //  [ 0  0 97]"
     * A.det();  // 97
     *
     * @type {Matrix}
     * @throws Will throw an error if the matrix is not square.
     * @see Matrix#charpoly
     */
    get adjugate() {
        if(!this.isSquare())
            throw new Error("Tried to compute the adjugate of a non-square matrix");
        this._fadeev_leverrier(); // Computes the adjugate
        return this._cache.adjugate;
    }


    /**
     * @summary
     * Invalidate cached computations.
     *
     * @desc
     * Call this after modifying any matrix entries.
     *
     * @return {undefined}
     */
    invalidate() {
        if(this.__cache) delete this.__cache;
    }

    /**
     * @summary
     * Hint precomputed or known properties of the matrix.
     *
     * @desc
     * This is provided so that the methods in this class can select the best
     * algorithm depending on special properties of the matrix.
     *
     * @param {Object} hints - Hinted properties.
     * @param {boolean} hints.isSymmetric - Whether the matrix is symmetric.
     * @param {boolean} hints.isUpperTri - Whether the matrix is upper-triangular.
     * @param {boolean} hints.isUpperUni - Whether the matrix is upper-unitriangular.
     * @param {boolean} hints.isLowerTri - Whether the matrix is lower-triangular.
     * @param {boolean} hints.isLowerUni - Whether the matrix is lower-unitriangular.
     * @param {boolean} hints.isEchelon - Whether the matrix is in row-echelon form.
     * @param {boolean} hints.hasONCols - Whether the matrix has orthonormal columns.
     * @param {Root[]} hints.eigenvalues - The eigenvalues of the matrix.
     * @return {undefined}
     */
    hint({isSymmetric,
          isUpperTri,
          isUpperUni,
          isLowerTri,
          isLowerUni,
          isEchelon,
          hasONCols,
          eigenvalues}) {
        if(isSymmetric !== undefined)
            this._cache.isSymmetric = isSymmetric;
        if(hasONCols !== undefined) {
            this._cache.hasONCols = hasONCols;
            this._cache.rank = this.n;
        }
        if(isUpperTri !== undefined)
            this._cache.isUpperTri = isUpperTri;
        if(isUpperUni !== undefined) {
            this._cache.isUpperUni = isUpperUni;
            if(isUpperUni) {
                this._cache.isUpperTri = true;
                this._cache.isEchelon = true;
            }
        }
        if(isLowerTri !== undefined)
            this._cache.isLowerTri = isLowerTri;
        if(isLowerUni !== undefined) {
            this._cache.isLowerUni = isLowerUni;
            if(isLowerUni)
                this._cache.isLowerTri = true;
        }
        if(isEchelon !== undefined) {
            this._cache.isEchelon = isEchelon;
            if(isEchelon)
                this._cache.isUpperTri = true;
        }
        if(eigenvalues)
            this._cache.eigenvalues = eigenvalues;
    }

    /**
     * @summary
     * Insert a submatrix into the matrix at position `(i,j)`.
     *
     * @desc
     * This overwrites the entries of `this` with the relevant entries of `M`.
     *
     * @example {@lang javascript}
     * Matrix.zero(4, 5).insertSubmatrix(
     *     1, 2, Matrix.create([1, 1],
     *                         [1, 1])).toString(0);
     *   // "[0 0 0 0 0]
     *   //  [0 0 1 1 0]
     *   //  [0 0 1 1 0]
     *   //  [0 0 0 0 0]"
     *
     * @param {integer} i - The row to insert.
     * @param {integer} j - The column to insert.
     * @param {Matrix} M - The submatrix.
     * @return {Matrix} `this`, after modification.
     */
    insertSubmatrix(i, j, M) {
        for(let ii = 0; ii < M.m; ++ii)
            this[ii+i].splice(j, M.n, ...M[ii]);
        return this;
    }

    /**
     * @summary
     * Test if this matrix is equal to `other`.
     *
     * @desc
     * Two matrices are equal if they have the same dimensions, and all
     * entries are equal.
     *
     * @example {@lang javascript}
     * let A = Matrix.zero(2, 3);
     * let B = Matrix.create([0, 0,  0.01],
     *                       [0, 0, -0.01]);
     * A.equals(B);                 // false
     * A.equals(B, 0.05);           // true
     * A.equals(Matrix.zero(3, 2)); // false
     *
     * @param {Matrix} other - The matrix to compare.
     * @param {number} [ε=0] - Entries will test as equal if they are within `ε`
     *   of each other.  This is provided in order to account for rounding
     *   errors.
     * @return {boolean} True if the matrices are equal.
     */
    equals(other, ε=0) {
        if(this.m !== other.m || this.n !== other.n)
            return false;
        return this.every((v, i) => v.equals(other[i], ε));
    }

    /**
     * @summary
     * Create a new Matrix with the same entries.
     *
     * @return {Matrix} The new matrix.
     */
    clone() {
        return Matrix.from(this, row => row.clone());
    }

    /**
     * @summary
     * Return a string representation of the matrix.
     *
     * @example {@lang javascript}
     * Matrix.create([1, 2], [3, 4], [5, 6]).toString(1);
     *   // "[1.0 2.0]
     *   //  [3.0 4.0]
     *   //  [5.0 6.0]"
     *
     * @param {integer} [precision=4] - The number of decimal places to include.
     * @return {string} A string representation of the matrix.
     */
    toString(precision=4) {
        let strings = Array.from(
            this, row => Array.from(row, v => v.toFixed(precision)));
        let colLens = Array.from(
            range(this.n), j => Math.max(...strings.map(row => row[j].length)));
        return strings.map(
            row => '[' + row.map(
                (val, j) => val.padStart(colLens[j], ' ')).join(' ') + ']')
            .join('\n');
    }

    /**
     * @summary
     * Return the `i`th row.
     *
     * @example {@lang javascript}
     * Matrix.create([1, 2], [3, 4], [5, 6]).row(1).toString(0);  // "[3 4]"
     *
     * @param {integer} i - The row to return.
     * @return {Vector} The `i`th row of `this`.
     */
    row(i) {
        return this[i];
    }

    /**
     * @summary
     * Return an iterable over the rows.
     *
     * @example {@lang javascript}
     * let A = Matrix.create([1, 2], [3, 4], [5, 6]);
     * Array.from(A.rows(), row => row.toString(0));
     *   // ["[1 2]", "[3 4]", "[5 6]"]
     *
     * @return {Iterable.<Vector>} An iterable over the rows.
     */
    rows() {
        return this[Symbol.iterator]();
    }

    /**
     * @summary
     * Return the `j`th column.
     *
     * @example {@lang javascript}
     * Matrix.create([1, 2], [3, 4], [5, 6]).col(1).toString(0);  // "[2 4 6]"
     *
     * @param {integer} j - The column to return.
     * @return {Vector} The `j`th column of `this`.
     */
    col(j) {
        return Vector.from(this, row => row[j]);
    }

    /**
     * @summary
     * Return an iterable over the columns.
     *
     * @example {@lang javascript}
     * let A = Matrix.create([1, 2], [3, 4], [5, 6]);
     * Array.from(A.cols(), col => col.toString(0));
     *   // ["[1 3 5]", "[2 4 6]"]
     *
     * @return {Iterable.<Vector>} An iterable over the columns.
     */
    cols() {
        let self = this;
        return (function*() {
            for(let j = 0; j < self.n; ++j)
                yield self.col(j);
        })();
    }

    /**
     * @summary
     * Return an iterable over the diagonal entriese.
     *
     * @example {@lang javascript}
     * let A = Matrix.create([1, 2], [3, 4], [5, 6]);
     * Array.from(A.diag());  // [1, 4]
     *
     * @return {Iterable.<number>} An iterable over the diagonal entries.
     */
    diag() {
        let self = this;
        return (function*() {
            for(let j = 0; j < Math.min(self.m, self.n); ++j)
                yield self[j][j];
        })();
    }

    /**
     * @summary
     * Return the `(i, j)` minor of the matrix.
     *
     * @desc
     * The `(i, j)` minor of a matrix is the matrix obtained by deleting the
     * `i`th row and the `j`th column.
     *
     * @example {@lang javascript}
     * let A = Matrix.create([1, 2, 3],
     *                       [4, 5, 6],
     *                       [7, 8, 9]);
     * A.minor(0, 1).toString(0);
     *  // "[4 6]
     *  //  [7 9]"
     *
     * @param {integer} i - The row to delete.
     * @param {integer} j - The column to delete.
     * @return {Matrix} The `(i, j)` minor.
     * @see Matrix#cofactor
     */
    minor(i, j) {
        let ret = this.clone();
        ret.splice(i, 1);
        ret.forEach(row => row.splice(j, 1));
        return ret;
    }

    /**
     * @summary
     * Return the `(i, j)` cofactor of the matrix.
     *
     * @desc
     * The `(i, j)` cofactor of a matrix is `(-1)^(i+j)` times the determinant
     * of the `(i, j)` minor.  This only makes sense for square matrices.
     *
     * @example {@lang javascript}
     * let A = Matrix.create([1, 2, 3],
     *                       [4, 5, 6],
     *                       [7, 8, 9]);
     * A.cofactor(0, 1);  // 6
     *
     * @param {integer} i - The row.
     * @param {integer} j - The column.
     * @return {number} The `(i, j)` cofactor.
     * @throws Will throw an error if the matrix is not square.
     * @see Matrix#minor
     * @see Matrix#adjugate
     */
    cofactor(i, j) {
        return this.minor(i, j).det() * ((i + j) % 2 === 0 ? 1 : -1);
    }

    /**
     * @summary
     * A pivot position is recorded as a pair `[row, column]`.
     *
     * @typedef {number[]} Pivot
     */

    /**
     * @summary
     * A list of leading entries of each row.
     *
     * @desc
     * Zero rows do not have a leading entry, and are not included.
     *
     * @example {@lang javascript}
     * Matrix.create([0, 0, 0],
     *               [0, 1, 2],
     *               [0, 0, 0],
     *               [0, 0, 2]).leadingEntries();  // [[1, 1], [3, 2]]
     *
     * @param {number} [ε=0] - Entries smaller than this value are taken
     *   to be zero.
     * @return {Pivot[]}
     */
    leadingEntries(ε=0) {
        let entries = [];
        for(let [i, row] of this.entries()) {
            for(let j = 0; j < this.n; ++j) {
                if(Math.abs(row[j]) > ε) {
                    entries.push([i, j]);
                    break;
                }
            }
        }
        return entries;
    }


    // Properties the matrix can have

    /**
     * @summary
     * Test whether the matrix is square.
     *
     * @desc
     * A square matrix has the same number of rows as columns.
     *
     * @example {@lang javascript}
     * Matrix.create([1, 2], [3, 4]).isSquare();          // true
     * Matrix.create([1, 2], [3, 4], [5, 6]).isSquare();  // false
     *
     * @return {boolean} True if the matrix is square.
     */
    isSquare() {
        return this.m == this.n;
    }

    /**
     * @summary
     * Test whether the matrix is zero.
     *
     * @desc
     * All entries of the zero matrix are equal to zero.
     *
     * @example {@lang javascript}
     * Matrix.create([0, 0], [0, 0]).isZero();        // true
     * Matrix.create([0, 0], [0, 0.01]).isZero();     // false
     * Matrix.create([0, 0], [0, 0.01]).isZero(0.02); // true
     *
     * @param {number} [ε=0] - Entries smaller than this value are assumed
     *   to be zero.
     * @return {boolean} True if the matrix is zero.
     */
    isZero(ε=0) {
        return this.every(row => row.isZero(ε));
    }

    /**
     * @summary
     * Test whether the matrix is upper-triangular.
     *
     * @desc
     * A matrix is upper-triangular if all entries below the main diagonal are
     * equal to zero.  Equivalently, `this.transpose` is lower-triangular.
     *
     * @example {@lang javascript}
     * Matrix.create([1, 1, 1],
     *               [0, 1, 1],
     *               [0, 0, 2]).isUpperTri(); // true
     * Matrix.create([3, 1],
     *               [0, 1],
     *               [0, 0]).isUpperTri();    // true
     * Matrix.create([1, 1, 1],
     *               [0, 1, 1],
     *               [0, 2, 1]).isUpperTri(); // false
     *
     * @param {number} [ε=0] - Entries smaller than this value are assumed
     *   to be zero.
     * @return {boolean} True if the matrix is upper-triangular.
     * @see Matrix#isLowerTri
     */
    isUpperTri(ε=0) {
        if(this._cache.isUpperTri !== undefined)
            return this._cache.isUpperTri;
        for(let i = 1; i < this.m; ++i) {
            for(let j = 0; j < i; ++j) {
                if(Math.abs(this[i][j]) > ε) {
                    this._cache.isUpperTri = false;
                    return false;
                }
            }
        }
        this._cache.isUpperTri = true;
        return true;
    }

    /**
     * @summary
     * Test whether the matrix is upper-unitriangular.
     *
     * @desc
     * An upper-unitriangular matrix is an upper-triangular matrix with ones on
     * the diagonal.  Equivalently, `this.transpose` is lower-unitriangular.
     *
     * @example {@lang javascript}
     * Matrix.create([1, 2, 1],
     *               [0, 1, 1],
     *               [0, 0, 1]).isUpperUni(); // true
     * Matrix.create([1, 3],
     *               [0, 1],
     *               [0, 0]).isUpperUni();    // true
     * Matrix.create([1, 1, 1],
     *               [0, 1, 1],
     *               [0, 2, 1]).isUpperUni(); // false
     * Matrix.create([1, 1, 1],
     *               [0, 1, 1],
     *               [0, 0, 2]).isUpperUni(); // false
     *
     * @param {number} [ε=0] - Entries smaller than this value are assumed
     *   to be zero.
     * @return {boolean} True if the matrix is upper-unipotent.
     * @see Matrix#isLowerUni
     */
    isUpperUni(ε=0) {
        for(let d of this.diag()) {
            if(Math.abs(d - 1) > ε)
                return false;
        }
        return this.isUpperTri(ε);
    }

    /**
     * @summary
     * Test whether the matrix is lower-triangular.
     *
     * @desc
     * A matrix is lower-triangular if all entries above the main diagonal are
     * equal to zero.  Equivalently, `this.transpose` is upper-triangular.
     *
     * @example {@lang javascript}
     * Matrix.create([1, 0, 0],
     *               [1, 1, 0],
     *               [1, 1, 2]).isLowerTri(); // true
     * Matrix.create([3, 0],
     *               [1, 1],
     *               [2, 2]).isLowerTri();    // true
     * Matrix.create([1, 0, 0],
     *               [1, 1, 2],
     *               [1, 1, 1]).isLowerTri(); // false
     *
     * @param {number} [ε=0] - Entries smaller than this value are assumed
     *   to be zero.
     * @return {boolean} True if the matrix is lower-triangular.
     */
    isLowerTri(ε=0) {
        if(this._cache.isLowerTri !== undefined)
            return this._cache.isLowerTri;
        for(let i = 0; i < this.m; ++i) {
            for(let j = i+1; j < this.n; ++j) {
                if(Math.abs(this[i][j]) > ε) {
                    this._cache.isLowerTri = false;
                    return false;
                }
            }
        }
        this._cache.isLowerTri = true;
        return true;
    }

    /**
     * @summary
     * Test whether the matrix is lower-unitriangular.
     *
     * @desc
     * A lower-unitriangular matrix is a lower-triangular matrix with ones on
     * the diagonal.  Equivalently, `this.transpose` is upper-unitriangular.
     *
     * @example {@lang javascript}
     * Matrix.create([1, 0, 0],
     *               [2, 1, 0],
     *               [1, 1, 1]).isLowerUni(); // true
     * Matrix.create([1, 0],
     *               [2, 1],
     *               [3, 3]).isLowerUni();    // true
     * Matrix.create([1, 0, 0],
     *               [1, 1, 2],
     *               [1, 1, 1]).isLowerUni(); // false
     * Matrix.create([1, 0, 0],
     *               [1, 1, 0],
     *               [1, 1, 2]).isLowerUni(); // false
     *
     * @param {number} [ε=0] - Entries smaller than this value are assumed
     *   to be zero.
     * @return {boolean} True if the matrix is lower-unipotent.
     * @see Matrix#isUpperUni
     */
    isLowerUni(ε=0) {
        for(let d of this.diag()) {
            if(Math.abs(d - 1) > ε)
                return false;
        }
        return this.isLowerTri(ε);
    }

    /**
     * @summary
     * Test whether the matrix is diagonal.
     *
     * @desc
     * A matrix is diagonal if the only nonzero entries of the matrix are on the
     * diagonal.  Equivalently, the matrix is both upper- and lower-triangular.
     *
     * @example {@lang javascript}
     * Matrix.create([1, 0, 0],
     *               [0, 2, 0],
     *               [0, 0, 3]).isDiagonal(); // true
     * Matrix.create([1, 0],
     *               [0, 1],
     *               [0, 0]).isDiagonal();    // true
     * Matrix.create([0, 0, 0],
     *               [0, 0, 0],
     *               [0, 0, 0]).isDiagonal(); // true
     * Matrix.create([1, 0, 0],
     *               [0, 1, 2],
     *               [0, 0, 1]).isDiagonal(); // false
     *
     * @param {number} [ε=0] - Entries smaller than this value are assumed
     *   to be zero.
     * @return {boolean} True if the matrix is diagonal.
     */
    isDiagonal(ε=0) {
        return this.isLowerTri(ε) && this.isUpperTri(ε);
    }

    /**
     * @summary
     * Test whether the matrix is in row-echelon form.
     *
     * @desc
     * A matrix is in row-echelon form if the first nonzero entry of each row is
     * to the right of the first nonzero entry of the previous row, which
     * implies that all entries below a pivot are zero and all zero rows are at
     * the bottom.
     *
     * @example {@lang javascript}
     * Matrix.create([1,  0,  2],
     *               [0,  1, -1]).isEchelon();         // true
     * Matrix.create([0, 1, 8, 0]).isEchelon();        // true
     * Matrix.create([1, 17,  0],
     *               [0,  0,  1]).isEchelon();         // true
     * Matrix.zero(2, 3).isEchelon();                  // true
     * Matrix.create([2, 1],
     *               [0, 1]).isEchelon();              // true
     * Matrix.create([2,  7, 1, 4],
     *               [0,  0, 2, 1],
     *               [0,  0, 0, 3]).isEchelon();       // true
     * Matrix.create([1, 17, 0],
     *               [0,  1, 1]).isEchelon();          // true
     * Matrix.create([2,  1, 3],
     *               [0,  0, 0]).isEchelon();          // true
     * Matrix.create([2,  7, 1, 4],
     *               [0,  0, 2, 1],
     *               [0,  0, 1, 3]).isEchelon();       // false
     * Matrix.create([0, 17, 0],
     *               [0,  2, 1]).isEchelon();          // false
     * Matrix.create([2,  1],
     *               [2,  1]).isEchelon();             // false
     * Matrix.create([0,1,0,0]).transpose.isEchelon(); // false
     *
     * @param {number} [ε=0] - Entries smaller than this value are assumed
     *   to be zero.
     * @return {boolean} True if the matrix is in row echelon form.
     */
    isEchelon(ε=0) {
        if(this._cache.isEchelon) return this._cache.isEchelon;
        let i0 = -1;
        let j0 = -1;
        for(let [i, j] of this.leadingEntries(ε)) {
            if(i !== i0 + 1 || j <= j0) {
                this._cache.isEchelon = false;
                return false;
            }
            i0 = i;
            j0 = j;
        }
        this._cache.isEchelon = true;
        return true;
    }

    /**
     * @summary
     * Test whether the matrix is in reduced row-echelon form.
     *
     * @desc
     * A matrix is in reduced row-echelon form if it is in row-echelon form, and
     * in addition, all pivots are equal to one, and all entries above a pivot
     * are equal to zero.
     *
     * @example {@lang javascript}
     * Matrix.create([1,  0,  2],
     *               [0,  1, -1]).isRREF();         // true
     * Matrix.create([0, 1, 8, 0]).isRREF();        // true
     * Matrix.create([1, 17,  0],
     *               [0,  0,  1]).isRREF();         // true
     * Matrix.zero(2, 3).isRREF();                  // true
     * Matrix.create([2, 1],
     *               [0, 1]).isRREF();              // false
     * Matrix.create([2,  7, 1, 4],
     *               [0,  0, 2, 1],
     *               [0,  0, 0, 3]).isRREF();       // false
     * Matrix.create([1, 17, 0],
     *               [0,  1, 1]).isRREF();          // false
     * Matrix.create([2,  1, 3],
     *               [0,  0, 0]).isRREF();          // false
     * Matrix.create([2,  7, 1, 4],
     *               [0,  0, 2, 1],
     *               [0,  0, 1, 3]).isRREF();       // false
     * Matrix.create([0, 17, 0],
     *               [0,  2, 1]).isRREF();          // false
     * Matrix.create([2,  1],
     *               [2,  1]).isRREF();             // false
     * Matrix.create([0,1,0,0]).transpose.isRREF(); // false
     *
     * @param {number} [ε=0] - Entries smaller than this value are assumed
     *   to be zero.
     * @return {boolean} True if the matrix is in reduced row-echelon form.
     */
    isRREF(ε=0) {
        let i0 = -1;
        let j0 = -1;
        for(let [i, j] of this.leadingEntries(ε)) {
            if(i !== i0 + 1 || j <= j0)
                return false;
            if(this[i][j] != 1)
                return false;
            for(let k = 0; k < i; ++k)
                if(this[k][j] !== 0)
                    return false;
            i0 = i;
            j0 = j;
        }
        return true;
    }

    /**
     * @summary
     * Test whether the matrix has full row rank.
     *
     * @desc
     * A matrix has full row rank if the following equivalent conditions hold:
     *  * There is a pivot in every row.
     *  * The number of pivots equals `m`.
     *  * The column space has dimension `m`.
     *  * The matrix equation `Ax=b` has a solution for every choice of `b`.
     *
     * Running this method involves computing the rank if it has not been
     * computed already, unless `m > n`, in which case the matrix cannot have
     * full row rank.
     *
     * @example {@lang javascript}
     * Matrix.create([2,  7, 1, 4],
     *               [0,  0, 2, 1],
     *               [0,  0, 0, 3]).isFullRowRank();  // true
     * Matrix.create([2,  7, 1, 4],
     *               [0,  0, 2, 1],
     *               [0,  0, 0, 0]).isFullRowRank();  // false
     *
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting and projecting.
     * @return {boolean} True if the matrix has full row rank.
     * @see Matrix#rank
     */
    isFullRowRank(ε=1e-10) {
        if(this.m > this.n) return false; // Don't need to row reduce
        return this.rank(ε) == this.m;
    }

    /**
     * @summary
     * Test whether the matrix has full column rank.
     *
     * @desc
     * A matrix has full column rank if the following equivalent conditions hold:
     *  * There is a pivot in every column.
     *  * The number of pivots equals `n`.
     *  * The null space is zero.
     *  * The matrix equation `Ax=b` has at most one solution for every choice
     *    of `b`.
     *
     * Running this method involves computing the rank if it has not been
     * computed already, unless `m < n`, in which case the matrix cannot have
     * full column rank.
     *
     * @example {@lang javascript}
     * Matrix.create([2,  7, 1],
     *               [0,  1, 2],
     *               [0,  0, 7]).isFullColRank();  // true
     * Matrix.create([2,  7, 1],
     *               [0,  0, 2],
     *               [0,  0, 0]).isFullColRank();  // false
     *
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting and projecting.
     * @return {boolean} True if the matrix has full column rank.
     * @see Matrix#rank
     */
    isFullColRank(ε=1e-10) {
        if(this.m < this.n) return false; // Don't need to row reduce
        return this.rank(ε) == this.n;
    }

    /**
     * @summary
     * Test whether the matrix is invertible.
     *
     * @desc
     * A matrix is invertible if the following equivalent conditions hold:
     *  * The matrix is square and has the maximum number of pivots.
     *  * The matrix has full row rank and full column rank.
     *  * There is another matrix (namely, the inverse) such that the product
     *    with `this` on both sides is equal to the identity.
     *  * The number zero is not an eigenvalue of `this`.
     *  * The determinant of `this` is nonzero.
     *
     * Running this method involves computing the rank if it has not been
     * computed already.
     *
     * @example {@lang javascript}
     * Matrix.create([2,  7, 1],
     *               [0,  1, 2],
     *               [0,  0, 7]).isInvertible();     // true
     * Matrix.create([2,  7, 1],
     *               [0,  0, 2],
     *               [0,  0, 0]).isInvertible();     // false
     * Matrix.create([2,  7, 1, 4],
     *               [0,  0, 2, 1],
     *               [0,  0, 0, 3]).isInvertible();  // false
     *
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting and projecting.
     * @return {boolean} True if the matrix is invertible.
     * @see Matrix#inverse
     */
    isInvertible(ε=1e-10) {
        return this.isFullRowRank(ε) && this.isFullColRank(ε);
    }

    /**
     * @summary
     * Alias for `!this.isInvertible()`.
     *
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting and projecting.
     * @return {boolean} True if the matrix is not invertible.
     * @see Matrix#isInvertible
     */
    isSingular(ε=1e-10) {
        return !this.isInvertible(ε);
    }

    /**
     * @summary
     * Test if the matrix is has orthonormal columns.
     *
     * @desc
     * This means that `A^TA` is the identity matrix.
     *
     * @example {@lang javascript}
     * Matrix.create([1, 0],
     *               [0, 1],
     *               [0, 0]).hasONCols(); // true
     *
     * @param {number} [ε=1e-10] - Numbers smaller than this value are taken
     *   to be zero.
     * @return {boolean} True if the matrix has orthonormal columns.
     */
    hasONCols(ε=1e-10) {
        if(this._cache.hasONCols) return this._cache.hasONCols;
        this._cache.hasONCols =
            this.normal.equals(Matrix.identity(this.n), ε);
        return this._cache.hasONCols;
    }

    /**
     * @summary
     * Test if the matrix is orthogonal.
     *
     * @desc
     * A matrix is orthogonal if it is square and has orthonormal columns.
     * Equivalently `A` is square and `A^TA` is the identity matrix.
     *
     * @example {@lang javascript}
     * Matrix.create([ 3/5, 4/5],
     *               [-4/5, 3/5]).isOrthogonal();   // true
     * Matrix.create([ 3/5, 1],
     *               [-4/5, 0]).isOrthogonal();     // false
     * Matrix.create([ 3, 4],
     *               [-4, 3]).isOrthogonal();       // false
     * Matrix.create([1, 0],
     *               [0, 1],
     *               [0, 0]).isOrthogonal();        // false
     *
     * @param {number} [ε=1e-10] - Numbers smaller than this value are taken
     *   to be zero.
     * @return {boolean} True if the matrix is orthogonal.
     */
    isOrthogonal(ε=1e-10) {
        return this.isSquare() && this.hasONCols(ε);
    }

    /**
     * @summary
     * Test if the matrix is symmetric.
     *
     * @desc
     * A matrix is symmetric if it is equal to its transpose.  This method is a
     * shortcut for `this.equals(this.transpose, ε)`.
     *
     * @example {@lang javascript}
     * Matrix.create([1, 2, 3],
     *               [2, 4, 5],
     *               [3, 5, 6]).isSymmetric();  // true
     * Matrix.create([1, 2, 3],
     *               [4, 5, 6],
     *               [7, 8, 9]).isSymmetric();  // false
     * Matrix.create([1, 2, 3],
     *               [2, 4, 5]).isSymmetric();  // false
     *
     * @param {number} [ε=1e-10] - Entries will test as equal if they are within
     *   `ε` of each other.
     * @return {boolean} True if the matrix is symmetric.
     */
    isSymmetric(ε=1e-10) {
        if(this._cache.isSymmetric) return this._cache.isSymmetric;
        this._cache.isSymmetric = this.equals(this.transpose, ε);
        return this._cache.isSymmetric;
    }

    /**
     * @summary
     * Test if a symmetric matrix is positive-definite.
     *
     * @desc
     * A symmetric matrix is positive-definite if and only if `x^T A x > 0`
     * whenever `x` is a nonzero vector.  Equivalently, all eigenvalues of `A`
     * are positive, or all diagonal entries of `D` in the LDLT decomposition
     * are positive.
     *
     * This method works by computing the LDLT decomposition.
     *
     * @example {@lang javascript}
     * let A = Matrix.create([3, 2, 3],
     *                       [2, 7, 4],
     *                       [3, 4, 8]);
     * A.isPosDef();  // true
     *
     * @example {@lang javascript}
     * let A = Matrix.create([1, 2, 3],
     *                       [2, 5, 4],
     *                       [3, 4, 2]);
     * A.isPosDef();  // false
     * // Here is a witness to the indefiniteness of A:
     * let v = Vector.create(-7, 2, 1);
     * v.dot(A.apply(v));  // -11
     *
     * @param {number} [ε=1e-10] - Entries are considered to be zero if they are
     *   smaller than this.
     * @return {boolean} True if the matrix is positive-definite.
     * @throws Will throw an error if the matrix is not symmetric.
     * @see Matrix#LDLT
     */
    isPosDef(ε=1e-10) {
        if(!this.isSymmetric())
            throw new Error("Tried to test positive-definiteness of a non-symmetric matrix.");
        let LDLT = this.LDLT();
        if(!LDLT) return false;
        return LDLT.D.every(x => x > 0);
    }

    /**
     * @summary
     * Test if the matrix is diagonalizable.
     *
     * @desc
     * A matrix is diagonalizable if it is square and there exists an
     * invertible matrix `C` such that `CAC^(-1)` is diagonal.  Equivalently,
     * the matrix admits `n` linearly independent eigenvectors (the columns of
     * `C`).
     *
     * @example {@lang javascript}
     * Matrix.create([11/13, 22/39,  2/39],
     *               [-4/13, 83/39,  4/39],
     *               [-1/13, 11/39, 40/39]).isDiagonalizable();  // true
     * Matrix.create([    1,   1/2,     0],
     *               [-4/13, 83/39,  4/39],
     *               [ 5/13,  7/78, 34/39]).isDiagonalizable();  // false
     *
     * @param {number} [ε=1e-10] - Rounding factor.
     * @return {boolean} True if the matrix is diagonalizable.
     * @throws Will throw an error if the matrix is not square.
     * @see Matrix#diagonalize
     */
    isDiagonalizable(ε=1e-10) {
        if(this.isSymmetric(ε))
            return true;
        return !!this.diagonalize({ε});
    }


    // Matrix arithmetic

    /**
     * @summary
     * Add a Matrix in-place.
     *
     * @desc
     * This modifies the matrix in-place by adding the entries of `other`.
     *
     * @example {@lang javascript}
     * let A = Matrix.create([1, 2], [3, 4]);
     * let B = Matrix.create([2, 1], [4, 3]);
     * A.add(B);
     * A.toString(0);
     *   // "[3 3]
     *   //  [7 7]"
     *
     * @param {Matrix} other - The matrix to add.
     * @param {number} [factor=1] - Add `factor` times `other` instead of just
     *   adding `other`.
     * @return {Matrix} `this`
     * @throws Will throw an error if the matrices have different sizes.
     */
    add(other, factor=1) {
        if(this.m !== other.m || this.n !== other.n)
            throw new Error('Tried to add matrices of different sizes');
        this.forEach((row, i) => row.add(other[i], factor));
        return this;
    }

    /**
     * @summary
     * Subtract a Matrix in-place.
     *
     * @desc
     * This modifies the matrix in-place by subtracting the entries of `other`.
     *
     * @example {@lang javascript}
     * let A = Matrix.create([1, 2], [3, 4]);
     * let B = Matrix.create([2, 1], [4, 3]);
     * A.sub(B);
     * A.toString(0);
     *   // "[-1 1]
     *   //  [-1 1]"
     *
     * @param {Matrix} other - The matrix to subtract.
     * @return {Matrix} `this`
     * @throws Will throw an error if the matrices have different sizes.
     */
    sub(other) {
        return this.add(other, -1);
    }

    /**
     * @summary
     * Scale a Matrix in-place.
     *
     * @desc
     * This modifies the matrix in-place by multiplying all entries by `c`.
     *
     * @example {@lang javascript}
     * let A = Matrix.create([1, 2], [3, 4]);
     * A.scale(2);
     * A.toString(0);
     *   // "[2 4]
     *   //  [6 8]"
     *
     * @param {number} c - The number to multiply.
     * @return {Matrix} `this`
     */
    scale(c) {
        this.forEach(row => row.scale(c));
        return this;
    }

    /**
     * @summary
     * Multiply by another matrix.
     *
     * @desc
     * This creates a new matrix equal to the product of `this` and `other`.
     * The `(i, j)` entry of the product is the dot product of the `i`th row of
     * `this` and the `j`th column of `other`.  This only makes sense when the
     * number of rows of `other` equals the number of columns of `this`.
     *
     * This method runs in about `O(n^3)` time for an `n`x`n` matrix.
     * Amazingly, this is known not to be optimal: see [sub-cubic algorithms]{@link
     * https://en.wikipedia.org/wiki/Matrix_multiplication_algorithm#Sub-cubic_algorithms}
     * on Wikipedia.
     *
     * @example {@lang javascript}
     * let A = Matrix.create([1, 2], [3, 4], [5, 6]);
     * let B = Matrix.create([1, 2, 3], [4, 5, 6]);
     * A.mult(B).toString(0);
     *   // "[ 9 12 15]
     *   //  [19 26 33]
     *   //  [29 40 51]"
     *
     * @param {Matrix} other - The matrix to multiply.
     * @return {Matrix} The matrix product.
     * @throws Will throw an error if the matrices have incompatible
     *   dimensions.
     */
    mult(other) {
        if(other.m !== this.n)
            throw new Error('Cannot multiply matrices of incompatible dimensions');
        return Matrix.from(this,
            row => Vector.from(range(other.n), i => row.reduce(
                (a, v, j) => a + v * other[j][i], 0)));
    }

    /**
     * @summary
     * Multiply a matrix times a vector.
     *
     * @desc
     * This creates a new vector equal to the product of `this` and `v`.  The
     * `i`th entry of the product is the dot product of the `i`th row of `this`
     * with `v`.  This only makes sense when the number of entries of `v` equals
     * the number of columns of `this`.
     *
     * This is an optimized special case of {@link Matrix#mult}.
     *
     * @example {@lang javascript}
     * let A = Matrix.create([1, 2, 3], [4, 5, 6]);
     * B.apply(Vector.create(1, -1, 1)).toString(0); // "[2 5]"
     *
     * @param {Vector|number[]} v - The vector to multiply.  An array of numbers
     *   is promoted to a vector.
     * @return {Vector} The matrix-vector product.
     * @throws Will throw an error if the matrix and vector have incompatible
     *   dimensions.
     */
    apply(v) {
        if(v.length !== this.n)
            throw new Error('Cannot multiply matrix and vector of incompatible dimensions');
        return Vector.from(this, row => row.dot(v));
    }


    /**
     * @summary
     * Compute the inverse matrix.
     *
     * @desc
     * The inverse matrix is the unique matrix `A^(-1)` such that `A A^(-1)` and
     * `A^(-1) A` are both the identity matrix.
     *
     * This only makes sense for square matrices.  Throws an error if the matrix
     * is not square or is not invertible.
     *
     * This method runs Gauss&ndash;Jordan elimination, keeping track of the
     * row operations involved to produce the inverse matrix.  It runs in
     * about `O(n^3)` time, which is not optimal.
     *
     * @example {@lang javascript}
     * let A = Matrix.create([0,  1, 2],
     *                       [1,  0, 3],
     *                       [4, -3, 8]);
     * let Ainv = A.inverse();
     * Ainv.toString(1);
     *   // "[-4.5  7.0 -1.5]
     *   //  [-2.0  4.0 -1.0]
     *   //  [ 1.5 -2.0  0.5]"
     * Ainv.mult(A).toString(1);
     *   // "[1.0 0.0 0.0]
     *   //  [0.0 1.0 0.0]
     *   //  [0.0 0.0 1.0]"
     * A.mult(Ainv).toString(1);
     *   // "[1.0 0.0 0.0]
     *   //  [0.0 1.0 0.0]
     *   //  [0.0 0.0 1.0]"
     *
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting.
     * @return {Matrix} The inverse matrix.
     * @throws Will throw an error if the matrix is not square or not
     *   invertible.
     * @see Matrix#rowOps
     */
    inverse(ε=1e-10) {
        if(!this.isSquare())
            throw new Error("Tried to invert a non-square matrix");
        let E = this.rowOps(ε);
        if(!this.isInvertible(ε))
            throw new Error("Tried to invert a singular matrix");
        E.hint({inverse: this});
        return E;
    }


    // Row operations

    /**
     * @summary
     * Scale a row by a constant.
     *
     * @desc
     * This modifies the matrix in-place by multiplying all entries of row `i`
     * by the scalar `c`.  This is one of the three fundamental row operations.
     *
     * @example {@lang javascript}
     * let A = Matrix.create([1, 2, 3],
     *                       [4, 5, 6],
     *                       [7, 8, 9]);
     * A.rowScale(1, 2);
     * A.toString(0);
     *   // "[1  2  3]
     *   //  [8 10 12]
     *   //  [7  8  9]"
     *
     * @param {integer} i - The row to scale.
     * @param {number} c - The scaling factor.
     * @param {integer} [start=0] - Only scale the entries `start...this.n`.
     *   Provided for optimizations when the entries before `start` are known to
     *   be zero.
     * @return {Matrix} `this`
     */
    rowScale(i, c, start=0) {
        this[i].scale(c, start);
        return this;
    }

    /**
     * @summary
     * Add a constant times one row to another row.
     *
     * @desc
     * This modifies the matrix in-place by adding row `i2` times the scalar `c`
     * to row `i1`.  This is one of the three fundamental row operations.
     *
     * @example {@lang javascript}
     * let A = Matrix.create([1, 2, 3],
     *                       [4, 5, 6],
     *                       [7, 8, 9]);
     * A.rowReplace(1, 2, -1);
     * A.toString(0);
     *   // "[ 1  2  3]
     *   //  [-3 -3 -3]
     *   //  [ 7  8  9]"
     *
     * @param {integer} i1 - The row to replace.
     * @param {integer} i2 - The row to add.
     * @param {number} c - The scaling factor.
     * @param {integer} [start=0] - Only add the entries `start...this.n`.
     *   Provided for optimizations when the entries before `start` are known to
     *   be zero.
     * @return {Matrix} `this`
     */
    rowReplace(i1, i2, c, start=0) {
        this[i1].add(this[i2], c, start);
        return this;
    }

    /**
     * @summary
     * Swap two rows.
     *
     * @desc
     * This exchanges rows `i1` and `i2`.  This is one of the three fundamental
     * row operations.
     *
     * @example {@lang javascript}
     * let A = Matrix.create([1, 2, 3],
     *                       [4, 5, 6],
     *                       [7, 8, 9]);
     * A.rowSwap(0, 1);
     * A.toString(0);
     *   // "[4 5 6]
     *   //  [1 2 3]
     *   //  [7 8 9]"
     *
     * @param {integer} i1 - The first row.
     * @param {integer} i2 - The second row.
     * @return {Matrix} `this`
     */
    rowSwap(i1, i2) {
        [this[i1], this[i2]] = [this[i2], this[i1]];
        return this;
    }


    // Gaussian elimination

    /**
     * @summary
     * The core data computed by Gaussian elimination.
     *
     * @desc
     * When performing Gaussian elimination on a computer, it is easy to keep
     * track of the row operations performed, and various other data.
     *
     * @typedef PLUData
     * @type {Object}
     * @property {number[]} P - A permutation of the numbers `1...m-1`.
     *   by {@link Matrix.permutation}.
     * @property {Matrix} L - An `m`x`m` lower-triangular matrix with ones on
     *   the diagonal.
     * @property {Matrix} U - An `m`x`n` matrix in row echelon form.
     * @property {Matrix} E - An `m`x`m` invertible matrix equal to `L^(-1)P`,
     *   so `EA = U`.
     * @property {Pivot[]} pivots - An array of pivot positions.
     * @property {number} det - The determinant of the matrix (only set for
     *   square matrices).
     *
     * @see Matrix#PLU
     */

    /**
     * @summary
     * Compute a `PA = LU` factorization.
     *
     * @desc
     * This computes an `m`x`m` permutation matrix `P`, an `m`x`m`
     * lower-triangular matrix `L` with ones on the diagonal, and a matrix `U`
     * in row echelon form, such that `PA = LU`.  Along the way, it computes the
     * sign of the permutation `P`, the pivot positions of the matrix, and an
     * invertible `m`x`m` matrix `E` such that `EA = U`.  (The matrix `E` is the
     * product of the elementary matrices corresponding to the row operations
     * performed on `A`).
     *
     * This is the core method that implements Gaussian elimination.  It uses
     * maximal partial pivoting for numerical stability.  This algorithm
     * requires about `2n^3/3` operations for an `n`x`n` matrix, which seems to
     * be standard.
     *
     * The permutation matrix `P` is returned as a list of `n` numbers defining
     * the permutation.  Use {@link Matrix.permutation} to turn it into a
     * matrix.
     *
     * @example {@lang javascript}
     * let A = Matrix.create([ 0, -3, -6,  4,  9],
     *                       [-1, -2, -1,  3,  1],
     *                       [-2, -3,  0,  3, -1],
     *                       [ 1,  4,  5, -9, -7]);
     * let {P, L, U, E, pivots} = A.PLU();
     * P;  // [2, 0, 3, 1]
     * // let PM = Matrix.permutation(P);
     * PM.toString(0);
     *   // "[0 0 1 0]
     *   //  [1 0 0 0]
     *   //  [0 0 0 1]
     *   //  [0 1 0 0]"
     * L.toString();
     *   // "[ 1.0000  0.0000  0.0000 0.0000]
     *   //  [ 0.0000  1.0000  0.0000 0.0000]
     *   //  [-0.5000 -0.8333  1.0000 0.0000]
     *   //  [ 0.5000  0.1667 -0.2000 1.0000]"
     * L.isLowerUni();  // true
     * U.toString();
     *   // "[-2.0000 -3.0000  0.0000  3.0000 -1.0000]
     *   //  [ 0.0000 -3.0000 -6.0000  4.0000  9.0000]
     *   //  [ 0.0000  0.0000  0.0000 -4.1667  0.0000]
     *   //  [ 0.0000  0.0000  0.0000  0.0000  0.0000]"
     * U.isEchelon();  // true
     * U.leadingEntries();  // [[0, 0], [1, 1], [2, 3]], the same as pivots
     * E.mult(A).toString();
     *   // "[-2.0000 -3.0000  0.0000  3.0000 -1.0000]
     *   //  [ 0.0000 -3.0000 -6.0000  4.0000  9.0000]
     *   //  [ 0.0000  0.0000  0.0000 -4.1667  0.0000]
     *   //  [ 0.0000  0.0000  0.0000  0.0000  0.0000]"
     * PM.mult(A).toString(2);
     *   // "[-2.00 -3.00  0.00  3.00 -1.00]
     *   //  [ 0.00 -3.00 -6.00  4.00  9.00]
     *   //  [ 1.00  4.00  5.00 -9.00 -7.00]
     *   //  [-1.00 -2.00 -1.00  3.00  1.00]"
     * L.mult(U).toString(2);
     *   // "[-2.00 -3.00  0.00  3.00 -1.00]
     *   //  [ 0.00 -3.00 -6.00  4.00  9.00]
     *   //  [ 1.00  4.00  5.00 -9.00 -7.00]
     *   //  [-1.00 -2.00 -1.00  3.00  1.00]"
     *
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting.
     * @return {PLUData} The `PA=LU` factorization, along with other data
     *   computed at the same time.
     * @see Matrix#isEchelon
     * @see Matrix#isLowerUni
     * @see Matrix.permutation
     * @see Matrix#rref
     * @see Matrix#det
     */
    PLU(ε=1e-10) {
        if(this._cache.PLU) return this._cache.PLU;

        let P = Array.from(range(this.m));
        let L = Matrix.identity(this.m);
        let E = Matrix.identity(this.m);
        let U = this.clone();
        let {m, n} = this;
        let pivots = [];
        let signP = 1;
        let det = 1;

        for(let curRow = 0, curCol = 0; curRow < m && curCol < n; ++curCol) {
            // Find maximal pivot
            let pivot = U[curRow][curCol], row = curRow;
            for(let i = curRow+1; i < m; ++i) {
                if(Math.abs(U[i][curCol]) > Math.abs(pivot)) {
                    pivot = U[i][curCol];
                    row = i;
                }
            }

            if(Math.abs(pivot) > ε) {
                // curCol is a pivot column
                if(row != curRow) {
                    // Row swap
                    [P[row], P[curRow]] = [P[curRow], P[row]];
                    signP *= -1;
                    U.rowSwap(row, curRow);
                    E.rowSwap(row, curRow);
                    for(let j = 0; j < curRow; ++j)
                        [L[row][j], L[curRow][j]]
                            = [L[curRow][j], L[row][j]];
                }
                // Eliminate
                for(let i = curRow+1; i < m; ++i) {
                    let l = U[i][curCol] / pivot;
                    L[i][curRow] = l;
                    U[i][curCol] = 0;
                    U.rowReplace(i, curRow, -l, curCol+1);
                    E.rowReplace(i, curRow, -l);
                }
                pivots.push([curRow, curCol]);
                if(m === n) det *= pivot;
                curRow++;

            } else {
                // Clear the column so U is really upper-triangular
                for(let i = curRow; i < m; ++i)
                    U[i][curCol] = 0;
            }
        }

        L.hint({isLowerUni: true});
        U.hint({isEchelon: true});
        this._cache.PLU = {P, L, U, E, pivots};
        if(m === n)
            this._cache.PLU.det = pivots.length === n ? det * signP : 0;
        return this._cache.PLU;
    }

    /**
     * @summary
     * Compute the pivot positions of the matrix.
     *
     * @desc
     * The pivot positions are the leading entries of a row-echelon form of the
     * matrix.
     *
     * This is a shortcut for `this.PLU(ε).pivots`.
     *
     * @example {@lang javascript}
     * let A = Matrix.create([ 0, -3, -6,  4,  9],
     *                       [-1, -2, -1,  3,  1],
     *                       [-2, -3,  0,  3, -1],
     *                       [ 1,  4,  5, -9, -7]);
     * A.PLU().U.toString();
     *   // "[-2.0000 -3.0000  0.0000  3.0000 -1.0000]
     *   //  [ 0.0000 -3.0000 -6.0000  4.0000  9.0000]
     *   //  [ 0.0000  0.0000  0.0000 -4.1667  0.0000]
     *   //  [ 0.0000  0.0000  0.0000  0.0000  0.0000]"
     * A.pivots();  // [[0, 0], [1, 1], [2, 3]]
     *
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting.
     * @return {Pivot[]} The pivot positions.
     * @see Matrix#PLU
     * @see Matrix#leadingEntries
     */
    pivots(ε=1e-10) {
        return this.PLU(ε).pivots;
    }

    /**
     * @summary
     * Compute the rank of the matrix.
     *
     * @desc
     * The rank of the matrix is the dimension of the column space.  It is equal
     * to the number of pivots.
     *
     * The rank is computed in {@link Matrix#PLU} and {@link Matrix#QR}, and is
     * cached.  If the rank has not been computed yet, then
     * `this.pivots(ε).length` is returned.
     *
     * @example {@lang javascript}
     * let A = Matrix.create([ 0, -3, -6,  4,  9],
     *                       [-1, -2, -1,  3,  1],
     *                       [-2, -3,  0,  3, -1],
     *                       [ 1,  4,  5, -9, -7]);
     * A.PLU().U.toString();
     *   // "[-2.0000 -3.0000  0.0000  3.0000 -1.0000]
     *   //  [ 0.0000 -3.0000 -6.0000  4.0000  9.0000]
     *   //  [ 0.0000  0.0000  0.0000 -4.1667  0.0000]
     *   //  [ 0.0000  0.0000  0.0000  0.0000  0.0000]"
     * A.pivots();  // [[0, 0], [1, 1], [2, 3]]
     * A.rank();    // 3
     *
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting.
     * @return {integer} The rank of the matrix.
     * @see Matrix#PLU
     * @see Matrix#QR
     */
    rank(ε=1e-10) {
        if(this._cache.rank === undefined) // Hasn't been computed yet
            this._cache.rank = this.pivots(ε).length;
        return this._cache.rank;
    }

    /**
     * @summary
     * Compute the nullity of the matrix.
     *
     * @desc
     * The nullity of a matrix is the dimension of its null space.  This is
     * equal to the number of free variables, which is the number of columns
     * without pivots.
     *
     * According to the Rank-Nullity Theorem, the rank plus the nullity equals
     * `n`, the number of columns.  Therefore, this method simply returns
     * `this.n - this.rank(ε)`.
     *
     * @example {@lang javascript}
     * let A = Matrix.create([ 0, -3, -6,  4,  9],
     *                       [-1, -2, -1,  3,  1],
     *                       [-2, -3,  0,  3, -1],
     *                       [ 1,  4,  5, -9, -7]);
     * A.nullity(); // 2
     * A.rank();    // 3
     *
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting.
     * @return {integer} The nullity of the matrix.
     * @see Matrix#rank
     */
    nullity(ε=1e-10) {
        return this.n - this.rank(ε);
    }

    /**
     * @summary
     * Compute the matrix determinant.
     *
     * @desc
     * The determinant is computed as a side-effect of {@link Matrix#PLU}: the
     * determinant of `A` equals the determinant of `P` times the determinant of
     * `U`.  Since `P` is a permutation matrix and `U` is upper-triangular, both
     * quantities are easy to compute.
     *
     * The `PLU` decomposition does not use exact arithmetic: pivots smaller
     * than `ε` are considered to be zero.  Also, computing the echelon form may
     * involve dividing by large pivots.  Therefore, `this.det()` may fail to be
     * an integer even when `this` has all integer entries.  On the other hand,
     * the constant coefficient of the [characteristic polynomial]{@link
     * Matrix#charpoly} is also equal to the determinant, is computed using
     * exact arithmetic (but in `O(n^4)` time), and will always produce an
     * integer determinant for matrices with integer entries.
     *
     * This only makes sense for square matrices.  Throws an error if the matrix
     * is not square.
     *
     * @example {@lang javascript}
     * Matrix.create([1, 6,4],
     *               [2,-1,3],
     *               [5, 0,1]).det();  // 97
     *
     * @example {@lang javascript}
     * let A = Matrix.create([3, 4],
     *                       [5, 6]);
     * A.det();             // -2.0000000000000018
     * A.charpoly.eval(0);  // -2
     *
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting.
     * @return {number} The matrix determinant.
     * @throws Will throw an error if the matrix is not square.
     * @see Matrix#PLU
     */
    det(ε) {
        if(!this.isSquare())
            throw new Error("Tried to compute the determinant of a non-square matrix");
        return this.PLU(ε).det;
    }

    /**
     * @summary
     * Compute the reduced row-echelon form of the matrix.
     *
     * @desc
     * The reduced row-echelon form of the matrix is obtained from a row-echelon
     * form by scaling so all pivots are equal to 1, and performing row
     * replacements so that all entries above a pivot are equal to zero.  This
     * is called Gauss&ndash;Jordan elimination.
     *
     * This method runs `this.PLU(ε)` to first compute the row-echelon form
     * `U`.  Running `rref()` after `PLU()` takes about 50% more time.
     *
     * @example {@lang javascript}
     * let A = Matrix.create([ 0, -3, -6,  4,  9],
     *                       [-1, -2, -1,  3,  1],
     *                       [-2, -3,  0,  3, -1],
     *                       [ 1,  4,  5, -9, -7]);
     * let R = A.rref();
     * R.toString(1);
     *   // "[1.0 0.0 -3.0 0.0  5.0]
     *   //  [0.0 1.0  2.0 0.0 -3.0]
     *   //  [0.0 0.0  0.0 1.0  0.0]
     *   //  [0.0 0.0  0.0 0.0  0.0]"
     * R.isRREF();  // true
     *
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting.
     * @return {Matrix} The reduced row-echelon form of the matrix.
     * @see Matrix#PLU
     * @see Matrix#isRREF
     */
    rref(ε=1e-10) {
        if(this._cache.rref) return this._cache.rref;
        let {U, E, pivots} = this.PLU(ε);
        let rowOps = E.clone();
        let rref = U.clone();
        // Start with the last pivot
        for(let k = pivots.length-1; k >= 0; --k) {
            let [row, col] = pivots[k];
            let pivot = rref[row][col];
            for(let i = 0; i < row; ++i) {
                rref.rowReplace(i, row, -rref[i][col]/pivot, col+1);
                rowOps.rowReplace(i, row, -rref[i][col]/pivot);
                rref[i][col] = 0;
            }
            rref.rowScale(row, 1/pivot, col+1);
            rowOps.rowScale(row, 1/pivot);
            rref[row][col] = 1;
        }
        this._cache.rref = rref;
        this._cache.rowOps = rowOps;
        return this._cache.rref;
    }

    /**
     * @summary
     * Return the row operations used in Gauss&ndash;Jordan elimination.
     *
     * @desc
     * This returns an invertible `m`x`m` matrix `E` such that `EA` is the
     * reduced row-echelon form of `A`.  This matrix is the product of the
     * elementary matrices corresponding to the row operations performed on `A`
     * to reduce it to reduced row-echelon form.  It is computed as a side
     * effect of running `this.rref(ε)`.
     *
     * If `A` is invertible then `E` is the inverse of `A`.
     *
     * @example {@lang javascript}
     * let A = Matrix.create([ 0, -3, -6,  4,  9],
     *                       [-1, -2, -1,  3,  1],
     *                       [-2, -3,  0,  3, -1],
     *                       [ 1,  4,  5, -9, -7]);
     * let E = A.rowOps();
     * E.toString(2);
     *   // "[ 0.60 0.00 -0.44  0.12]
     *   //  [-0.60 0.00 -0.16 -0.32]
     *   //  [-0.20 0.00 -0.12 -0.24]
     *   //  [ 0.00 1.00 -0.40  0.20]"
     * E.mult(A).toString(1);
     *   // "[1.0 0.0 -3.0 0.0  5.0]
     *   //  [0.0 1.0  2.0 0.0 -3.0]
     *   //  [0.0 0.0  0.0 1.0  0.0]
     *   //  [0.0 0.0  0.0 0.0  0.0]"
     * A.rref().toString(1);
     *   // "[1.0 0.0 -3.0 0.0  5.0]
     *   //  [0.0 1.0  2.0 0.0 -3.0]
     *   //  [0.0 0.0  0.0 1.0  0.0]
     *   //  [0.0 0.0  0.0 0.0  0.0]"
     *
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting.
     * @return {Matrix} The matrix recording the row operations.
     * @see Matrix#rref
     */
    rowOps(ε=1e-10) {
        if(this._cache.rowOps) return this._cache.rowOps;
        this.rref(ε);
        return this._cache.rowOps;
    }

    /**
     * @summary
     * Solve `Ax=b` using reverse-substitution.
     *
     * @desc
     * This method assumes `A` is in row-echelon form.  Solving using
     * reverse-substitution is very efficient: it requires about `m(m+1)/2`
     * operations.
     *
     * @private
     */
    _solveReverseSubst(b, ε) {
        let x = Vector.zero(this.n);
        let pivots = this.leadingEntries(ε);
        let r = pivots.length;
        // Check if a solution exists
        for(let i = r; i < this.m; ++i) {
            if(Math.abs(b[i]) > ε)
                return null;
        }
        for(let p = r-1; p >= 0; --p) {
            let [row, col] = pivots[p];
            x[col] = b[row];
            for(let pp = p+1; pp < r; ++pp) {
                let [, col1] = pivots[pp];
                x[col] -= this[row][col1] * x[col1];
            }
            x[col] /= this[row][col];
        }
        return x;
    }

    /**
     * @summary
     * Solve `Ax=b` using forward-substitution.
     *
     * @desc
     * This method assumes `A` is lower-unitriangular.  Solving using
     * forward-substitution is very efficient: it requires about `m(m+1)/2`
     * operations.
     *
     * @private
     */
    _solveForwardSubst(b) {
        let x = b.clone();
        for(let i = 1; i < this.m; ++i) {
            for(let ii = 0; ii < i; ++ii)
                x[i] -= this[i][ii] * x[ii];
        }
        return x;
    }

    /**
     * @summary
     * Find some solution `x` of the equation `Ax=b`.
     *
     * @desc
     * This method tries to find the most efficient algorithm based on
     * properties of the matrix:
     *  * If the matrix is in echelon form, solve using reverse-substitution.
     *    This requires about `m(m+1)/2` operations.
     *  * If the matrix is lower-unitriangular, solve using
     *    forward-substitution.  This requires about `m(m+1)/2` operations.
     *  * If the matrix is symmetric, try an LDLT decomposition.  If it exists,
     *    solve using forward- and reverse-substitution.  This requires about
     *    `m(m+1)` operations once the `LDLT` decomposition has been computed.
     *  * If the matrix has no special properties, compute a PLU decomposition
     *    and then solve using forward- and reverse-substitution.  This requires
     *    about `m(m+1)` operations once the `PA=LU` decomposition has been
     *    computed.
     *
     * @example {@lang javascript}
     * let A = Matrix.create([ 0, -3, -6,  4,  9],
     *                       [-1, -2, -1,  3,  1],
     *                       [-2, -3,  0,  3, -1],
     *                       [ 1,  4,  5, -9, -7]);
     * let b = Vector.create(-3, 7, 10, -15);
     * let x = A.solve(b);
     * x.toString(0);           // "[-8 5 0 3 0]"
     * A.apply(x).toString(0);  // "[-3 7 10 15]"
     * let b1 = Vector.create(0, 1, 0, 0);
     * A.solve(b1);             // null
     *
     * @param {Vector} b - The right-hand side of the equation `Ax=b`.
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting.
     * @return {?Vector} Returns a solution `x`, or `null` if no solution
     *   exists.
     * @throws Will throw an error if `b.length != this.m`.
     * @see Matrix#PLU
     * @see Matrix#LDLT
     * @see Matrix#isEchelon
     * @see Matrix#isLowerUni
     */
    solve(b, ε=1e-10) {
        if(b.length != this.m)
            throw new Error("Incompatible dimensions of matrix and vector");

        if(this.isEchelon(ε))
            return this._solveReverseSubst(b, ε);
        if(this.isLowerUni(ε))
            return this._solveForwardSubst(b);
        if(this.isSymmetric(ε)) {
            let LDLT = this.LDLT(ε);
            if(LDLT) {
                let {L, D} = LDLT;
                let y = L.solve(b);
                for(let i = 0; i < this.n; ++i)
                    y[i] /= D[i];
                return L.transpose.solve(y);
            }
        }

        let {P, L, U} = this.PLU(ε);
        let r = this.rank(ε);
        // Solve LUx = PAx = Pb
        let Pb = Vector.from(range(this.m), i => b[P[i]]);
        let y = L.solve(Pb);
        return U.solve(y, ε);
    }

    /**
     * @summary
     * Find the shortest solution `x` of the equation `Ax=b`.
     *
     * @desc
     * This is obtained by finding some solution using {@link Matrix#solve},
     * then projecting onto the row space using {@link Matrix#projectRowSpace}.
     *
     * This takes about twice as long as {@link Matrix#solve} once `PA=LU`
     * decompositions for `this` and `this.transpose.normal` have been
     * computed.
     *
     * @example {@lang javascript}
     * let A = Matrix.create([ 0, -3, -6,  4,  9],
     *                       [-1, -2, -1,  3,  1],
     *                       [-2, -3,  0,  3, -1],
     *                       [ 1,  4,  5, -9, -7]);
     * let b = Vector.create(-3, 7, 10, -15);
     * let x = A.solveShortest(b);
     * x.toString();            // "[-0.1429 0.1429 0.7143 3.0000 -1.1429]"
     * A.apply(x).toString(0);  // "[-3 7 10 15]"
     * // No other solution `x` has smaller size
     * x.size.toFixed(4);       // "9.8995"
     * let b1 = Vector.create(0, 1, 0, 0);
     * A.solveShortest(b1);     // null
     *
     * @param {Vector} b - The right-hand side of the equation `Ax=b`.
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting.
     * @return {?Vector} Returns the shortest solution `x`, or `null` if no
     *   solution exists.
     * @throws Will throw an error if `b.length != this.m`.
     * @see Matrix#solve
     * @see Matrix#projectRowSpace
     */
    solveShortest(b, ε=1e-10) {
        let x = this.solve(b, ε);
        if(x === null) return null;
        return this.projectRowSpace(x, ε);
    }

    /**
     * @summary
     * Find some least-squares solution `x` of the equation `Ax=b`.
     *
     * @desc
     * A least-squares solution is a vector `x` such that `Ax` is as close to
     * `b` as possible.  Equivalently, it is a vector `x` such that `Ax-b` is
     * orthogonal to the columns of `A`.  This method finds one such solution;
     * the others are obtained by adding vectors in the null space.
     *
     * A least-squares solution is obtained by solving the normal equation
     * `A^TAx = A^Tb`, which is always consistent.  This takes about `O(n^2)`
     * time once a `PA=LU` decomposition of `this.normal` has been computed.
     *
     * If `Ax=b` is consistent, then a least-squares solution is the same as an
     * actual solution of the equation `Ax=b`.
     *
     * @example {@lang javascript}
     * let A = Matrix.create([ 0, -3, -6,  4,  9],
     *                       [-1, -2, -1,  3,  1],
     *                       [-2, -3,  0,  3, -1],
     *                       [ 1,  4,  5, -9, -7]);
     * let b = Vector.create(0, 1, 0, 0);
     * let x = A.solveLeastSquares(b);
     * x.toString();           // "[-0.1667 0.0000 0.0000 0.0000 0.0000]"
     * A.apply(x).toString();  // "[0.0000 0.1667 0.3333 -0.1667]"
     * let Axmb = A.apply(x).sub(b);
     * // No other value of `x` makes `Axmb` shorter than this.
     * Axmb.size.toFixed(4);  // "0.9129"
     * // Axmb is orthogonal to the columns of A
     * A.transpose.apply(Axmb).isZero(1e-10);  // true
     *
     * // If `Ax=b` is consistent, then a least-squares solution is an ordinary
     * // solution.
     * let b1 = Vector.create(-3, 7, 10, -15);
     * let x1 = A.solveLeastSquares(b1);
     * x1.toString(0); // "[-8 5 0 3 0]"
     * let x2 = A.solve(b1);
     * x2.toString(0); // "[-8 5 0 3 0]"
     *
     * @param {Vector} b - The right-hand side of the equation `Ax=b`.
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting.
     * @return {Vector} A least-squares solution.
     * @throws Will throw an error if `b.length != this.m`.
     */
    solveLeastSquares(b, ε=1e-10) {
        return this.normal.solve(this.transpose.apply(b), ε);
    }

    /**
     * @summary
     * Find the shortest least-squares solution `x` to the equation `Ax=b`.
     *
     * @desc
     * This is obtained by finding some solution using {@link
     * Matrix#solveLeastSquares}, then projecting onto the row space using
     * {@link Matrix#projectRowSpace}.  It takes about twice as long as {@link
     * Matrix#solve} once `PA=LU` decompositions for `this` and
     * `this.transpose.normal` have been computed.
     *
     * If the [pseudo-inverse]{@link Matrix#pseudoInverse} `A^+` has been
     * computed, then this method multiplies by `A^+`.  If you need to run this
     * method many times, compute `A^+` first.
     *
     * @example {@lang javascript}
     * let A = Matrix.create([ 0, -3, -6,  4,  9],
     *                       [-1, -2, -1,  3,  1],
     *                       [-2, -3,  0,  3, -1],
     *                       [ 1,  4,  5, -9, -7]);
     * let b = Vector.create(0, 1, 0, 0);
     * let x = A.solveLeastSquaresShortest(b);
     * x.toString();           // "[-0.0476 -0.0714 0.0000 0.0000 -0.0238]"
     * A.apply(x).toString();  // "[0.0000 0.1667 0.3333 -0.1667]"
     * // No other least-squares solution `x` has smaller size
     * x.size.toFixed(4);      // "0.0891"
     *
     * @param {Vector} b - The right-hand side of the equation.
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting.
     * @return {Vector} The shortest least-squares solution.
     * @throws Will throw an error if `b.length != this.m`.
     * @see Matrix#solveLeastSquares
     * @see Matrix#projectRowSpace
     * @see Matrix#pseudoInverse
     */
    solveLeastSquaresShortest(b, ε=1e-10) {
        if(this._cache.pinv)
            return this._cache.pinv.apply(b);
        return this.projectRowSpace(this.solveLeastSquares(b, ε));
    }

    /**
     * @summary
     * Compute the pseudo-inverse matrix.
     *
     * @desc
     * The pseudo-inverse `A^+` of a matrix `A` has the following geometric
     * description: if `v` is a vector in `R^m`, then `A^+v` is obtained by
     * projecting onto `Col(A)` to obtain a vector `v1`, then solving `Ax=v1` to
     * obtain a solution `x`, then projecting `x` onto the row space of `A`.
     * The resulting vector is simply the shortest least-squares solution of
     * `Ax=v`: that is, `A^+v = A.solveLeastSquaresShortest(v)`.  Hence to
     * compute the columns of `A^+`, one has to find the shortest least-squares
     * solutions of `Ax=ei` for `i=1,...,m` (the unit coordinate vectors).
     *
     * This method runs {@link Matrix#solveLeastSquaresShortest} for each of the
     * `m` unit coordinate vectors in `R^m`.  In practice, the pseudo-inverse is
     * usually computed by finding the singular value decomposition `A = U Σ
     * V^T`; then `A^+ = U Σ^+ V^T`, where `Σ^+` is the same as `Σ` with the
     * nonzero entries inverted.  Since this matrix library does not include a
     * high-quality SVD algorithm, this method is implemented in the naïve way.
     *
     * @example {@lang javascript}
     * let A = Matrix.create([ 0, -3, -6,  4,  9],
     *                       [-1, -2, -1,  3,  1],
     *                       [-2, -3,  0,  3, -1],
     *                       [ 1,  4,  5, -9, -7]);
     * let Ap = A.pseudoInverse();
     * Ap.toString();
     *   // "[-0.0857 -0.0476 -0.1752 -0.1124]
     *   //  [-0.1714 -0.0714 -0.2743 -0.1914]
     *   //  [-0.0857  0.0000 -0.0229 -0.0457]
     *   //  [-0.2000  0.0000 -0.1200 -0.2400]
     *   //  [ 0.0857 -0.0238 -0.0533  0.0124]"
     * A.mult(Ap).equals(A.colSpace().projectionMatrix(), 1e-10);  // true
     * Ap.mult(A).equals(A.rowSpace().projectionMatrix(), 1e-10);  // true
     *
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting.
     * @return {Matrix} The pseudo-inverse matrix.
     * @see Matrix#solveLeastSquaresShortest
     */
    pseudoInverse(ε=1e-10) {
        if(this._cache.pinv) return this._cache.pinv;
        this._cache.pinv = Matrix.from(
            range(this.m),
            i => this.solveLeastSquaresShortest(
                Vector.e(i, this.m), ε)).transpose;
        return this._cache.pinv;
    }

    /**
     * @summary
     * Project a vector onto the column space of the matrix.
     *
     * @desc
     * If `x` is a least-squares solution of `Ax=b`, then `Ax` is by definition
     * the projection of `b` onto the column space of `A`.  Therefore, this
     * method is a shortcut for `this.apply(this.solveLeastSquares(b, ε))`.
     *
     * @example {@lang javascript}
     * let A = Matrix.create([ 0, -3, -6,  4,  9],
     *                       [-1, -2, -1,  3,  1],
     *                       [-2, -3,  0,  3, -1],
     *                       [ 1,  4,  5, -9, -7]);
     * let b = Vector.create(0, 1, 0, 0);
     * let x = A.projectColSpace(b);
     * x.toString();  // "[0.0000 0.1667 0.3333 -0.1667]"
     * // b-x is orthogonal to the columns of A
     * A.transpose.apply(b.sub(x)).isZero(1e-10);  // true
     *
     * @param {Vector} b - A vector of length `m`.
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting.
     * @return {Vector} The projection of `b` onto the column space.
     * @throws Will throw an error if `b.length != this.m`.
     * @see Matrix#solveLeastSquares
     */
    projectColSpace(b, ε=1e-10) {
        return this.apply(this.solveLeastSquares(b, ε));
    }

    /**
     * @summary
     * Project a vector onto the row space of the matrix.
     *
     * @desc
     * This is an alias for `this.transpose.projectColSpace(b, ε)`.
     *
     * @param {Vector} b - A vector of length `n`.
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting.
     * @return {Vector} The projection of `b` onto the row space.
     * @throws Will throw an error if `b.length != this.n`.
     * @see Matrix#solveLeastSquares
     * @see Matrix#projectColSpace
     */
    projectRowSpace(b, ε=1e-10) {
        return this.transpose.projectColSpace(b, ε);
    }

    // The four fundamental subspaces

    /**
     * @summary
     * Return a basis for the null space of the matrix.
     *
     * @desc
     * The null space is the subspace of `R^n` consisting of solutions of the
     * matrix equation `Ax=0`.  This method computes a basis for the null space
     * by finding the parametric vector form of the general solution of `Ax=0`
     * using the reduced row echelon form of the matrix.  These are all vectors
     * of length `n`.
     *
     * If the null space is zero (i.e., the matrix has full column rank), then
     * the empty list is returned.
     *
     * This method takes negligible time once the reduced row-echelon form has
     * been computed.
     *
     * @example {@lang javascript}
     * let A = Matrix.create([ 0, -3, -6,  4,  9],
     *                       [-1, -2, -1,  3,  1],
     *                       [-2, -3,  0,  3, -1],
     *                       [ 1,  4,  5, -9, -7]);
     * let R = A.rref();
     * R.toString(0);
     *   // "[1 0 -3 0  5]
     *   //  [0 1  2 0 -3]
     *   //  [0 0  0 1  0]
     *   //  [0 0  0 0  0]"
     * A.nullBasis().map(vec => vec.toString(0));
     *   // ["[3 -2 1 0 0]", "[-5 3 0 0 1]"]
     * A.isFullColRank();  // false
     *
     * @example {@lang javascript}
     * let A = Matrix.create([1, 0],
     *                       [0, 1],
     *                       [0, 0]);
     * A.nullBasis();      // []
     * A.isFullColRank();  // true
     *
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting.
     * @return {Vector[]} A basis for the null space of the matrix.
     * @see Matrix#rref
     * @see Matrix#isFullColRank
     */
    nullBasis(ε=1e-10) {
        if(this._cache.nullBasis) return this._cache.nullBasis;
        let rref = this.rref(ε), pivots = this.pivots(ε).slice(), previous = [];
        let basis = [];
        for(let j = 0; j < this.n; ++j) {
            if(pivots.length && pivots[0][1] === j) {
                // Pivot column
                previous.push(pivots.shift());
                continue;
            }
            // Free column
            let v = Vector.zero(this.n);
            for(let [row, col] of previous)
                v[col] = -rref[row][j];
            v[j] = 1;
            basis.push(v);
        }
        this._cache.nullBasis = basis;
        return this._cache.nullBasis;
    }

    /**
     * @summary
     * Return the null space of the matrix.
     *
     * @desc
     * This is essentially a shortcut for `new Subspace(this.nullBasis(ε))`.
     *
     * @example {@lang javascript}
     * let A = Matrix.create([ 0, -3, -6,  4,  9],
     *                       [-1, -2, -1,  3,  1],
     *                       [-2, -3,  0,  3, -1],
     *                       [ 1,  4,  5, -9, -7]);
     * A.nullBasis().map(vec => vec.toString(0));
     *   // ["[3 -2 1 0 0]", "[-5 3 0 0 1]"]
     * A.nullSpace().toString(0);
     *   // "Subspace of R^5 of dimension 2 with basis
     *   //  [ 3] [-5]
     *   //  [-2] [ 3]
     *   //  [ 1] [ 0]
     *   //  [ 0] [ 0]
     *   //  [ 0] [ 1]"
     *
     * @example {@lang javascript}
     * let A = Matrix.create([1, 0],
     *                       [0, 1],
     *                       [0, 0]);
     * A.nullBasis(); // []
     * A.nullSpace().toString();  // "The zero subspace of R^2"
     *
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting.
     * @return {Subspace} The null space of the matrix.
     * @see Matrix#nullBasis
     */
    nullSpace(ε=1e-10) {
        if(this._cache.nullSpace) return this._cache.nullSpace;
        this._cache.nullSpace = new Subspace(
            this.nullBasis(ε), {n: this.n, isBasis: true});
        return this._cache.nullSpace;
    }

    /**
     * @summary
     * Return a basis for the column space of the matrix.
     *
     * @desc
     * The column space of a matrix is the subspace of `R^m` consisting of all
     * linear combinations of the columns.  This is the same as the set of all
     * vectors of the form `Ax` for `x` in `R^n`.
     *
     * This method computes a basis for the column space by finding the pivots
     * and returning the pivot columns.  These are all vectors of length `m`.
     *
     * This method takes negligible time once a `PA=LU` decomposition has been
     * computed.
     *
     * @example {@lang javascript}
     * let A = Matrix.create([ 0, -3, -6,  4,  9],
     *                       [-1, -2, -1,  3,  1],
     *                       [-2, -3,  0,  3, -1],
     *                       [ 1,  4,  5, -9, -7]);
     * A.pivots();  // [[0, 0], [1, 1], [2, 3]]
     * A.colBasis().map(vec => vec.toString(0));
     *   // ["[0 -1 -2 1]", "[-3 -2 -3 4]", "[4 3 3 -9]"]
     *
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting.
     * @return {Vector[]} A basis for the column space of the matrix.
     * @see Matrix#pivots
     */
    colBasis(ε=1e-10) {
        if(this._cache.colBasis) return this._cache.colBasis;
        this._cache.colBasis = Array.from(this.pivots(ε), ([, j]) => this.col(j));
        return this._cache.colBasis;
    }

    /**
     * @summary
     * Return the column space of the matrix.
     *
     * @example {@lang javascript}
     * let A = Matrix.create([ 0, -3, -6,  4,  9],
     *                       [-1, -2, -1,  3,  1],
     *                       [-2, -3,  0,  3, -1],
     *                       [ 1,  4,  5, -9, -7]);
     * A.colBasis().map(vec => vec.toString(0));
     *   // ["[0 -1 -2 1]", "[-3 -2 -3 4]", "[4 3 3 -9]"]
     * A.colSpace().toString(0);
     *   // "Subspace of R^4 of dimension 3 with basis
     *   //  [ 0] [-3] [ 4]
     *   //  [-1] [-2] [ 3]
     *   //  [-2] [-3] [ 3]
     *   //  [ 1] [ 4] [-9]"
     *
     * @desc
     * This is essentially a shortcut for `new Subspace(Array.from(this.cols()))`.
     *
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting.
     * @return {Subspace} The column space of the matrix.
     * @see Matrix#colBasis
     */
    colSpace(ε=1e-10) {
        if(this._cache.colSpace) return this._cache.colSpace;
        this._cache.colSpace = new Subspace(this, {n: this.m, ε});
        return this._cache.colSpace;
    }

    /**
     * @summary
     * Return a basis for the row space of the matrix.
     *
     * @desc
     * The row space of a matrix is the subspace of `R^n` consisting of all
     * linear combinations of the rows.  This is the same as the column space of
     * the transpose.
     *
     * This method computes a basis for the row space by returning the nonzero
     * rows of a row echelon-form of the matrix.  These are all vectors of
     * length `n`.
     *
     * This method takes negligible time once a `PA=LU` decomposition has been
     * computed.
     *
     * @example {@lang javascript}
     * let A = Matrix.create([ 0, -3, -6,  4,  9],
     *                       [-1, -2, -1,  3,  1],
     *                       [-2, -3,  0,  3, -1],
     *                       [ 1,  4,  5, -9, -7]);
     * A.PLU().U.toString();
     *   // "[-2.0000 -3.0000  0.0000  3.0000 -1.0000]
     *   //  [ 0.0000 -3.0000 -6.0000  4.0000  9.0000]
     *   //  [ 0.0000  0.0000  0.0000 -4.1667  0.0000]
     *   //  [ 0.0000  0.0000  0.0000  0.0000  0.0000]"
     * A.rowBasis().map(vec => vec.toString());
     *   // ["[-2.0000 -3.0000 0.0000 3.0000 -1.0000]",
     *   //  "[0.0000 -3.0000 -6.0000 4.0000 9.0000]",
     *   //  "[0.0000 0.0000 0.0000 -4.1667 0.0000]"]
     *
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting.
     * @return {Vector[]} A basis for the row space of the matrix.
     * @see Matrix#PLU
     */
    rowBasis(ε=1e-10) {
        if(this._cache.rowBasis) return this._cache.rowBasis;
        let {pivots, U} = this.PLU(ε);
        this._cache.rowBasis = pivots.map(([i,]) => U[i].clone());
        return this._cache.rowBasis;
    }

    /**
     * @summary
     * Return the row space of the matrix.
     *
     * @desc
     * This is essentially a shortcut for `new Subspace(Array.from(this.rows()))`,
     * except that the returned subspace will come equipped with the cached
     * output of `this.rowBasis()` if it has been computed.
     *
     * @example {@lang javascript}
     * let A = Matrix.create([ 0, -3, -6,  4,  9],
     *                       [-1, -2, -1,  3,  1],
     *                       [-2, -3,  0,  3, -1],
     *                       [ 1,  4,  5, -9, -7]);
     * // Note that the Subspace computed a basis by finding the pivots of the
     * // matrix whose columns are the generators it was provided.
     * A.rowSpace().toString(0);
     *   // "Subspace of R^5 of dimension 3 with basis
     *   //  [ 0] [ 1] [-2]
     *   //  [-3] [-2] [-3]
     *   //  [-6] [-1] [ 0]
     *   //  [ 4] [ 3] [ 3]
     *   //  [ 9] [ 1] [-1]"
     *
     * @example {@lang javascript}
     * let A = Matrix.create([ 0, -3, -6,  4,  9],
     *                       [-1, -2, -1,  3,  1],
     *                       [-2, -3,  0,  3, -1],
     *                       [ 1,  4,  5, -9, -7]);
     * // Precompute a basis of the row space
     * A.rowBasis();
     * A.rowSpace().toString();
     *   // "Subspace of R^5 of dimension 3 with basis
     *   //  [-2.0000] [ 0.0000] [ 0.0000]
     *   //  [-3.0000] [-3.0000] [ 0.0000]
     *   //  [ 0.0000] [-6.0000] [ 0.0000]
     *   //  [ 3.0000] [ 4.0000] [-4.1667]
     *   //  [-1.0000] [ 9.0000] [ 0.0000]"
     *
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting.
     * @return {Subspace} The row space of the matrix.
     * @see Matrix#rowBasis
     */
    rowSpace(ε=1e-10) {
        if(this._cache.rowSpace) return this._cache.rowSpace;
        if(this._cache.rowBasis)
            this._cache.rowSpace = new Subspace(this._cache.rowBasis,
                                                {n: this.n, isBasis: true, ε});
        else
            this._cache.rowSpace = new Subspace(this.transpose, {n: this.n, ε});
        return this._cache.rowSpace;
    }

    /**
     * @summary
     * Return a basis for the left null space of the matrix.
     *
     * @desc
     * The left null space of a matrix is the null space of the transpose.  This
     * is a subspace of `R^m`.
     *
     * This method computes a basis for the left null space by returning the
     * last `m-r` rows of `L^(-1) P`, where `P` and `L` come from the `PA = LU`
     * decomposition and `r` is the rank.  These are all vectors of length `m`.
     *
     * This method takes negligible time once a `PA=LU` decomposition has been
     * computed.
     *
     * @example {@lang javascript}
     * let A = Matrix.create([ 0, -3, -6,  4,  9],
     *                       [-1, -2, -1,  3,  1],
     *                       [-2, -3,  0,  3, -1],
     *                       [ 1,  4,  5, -9, -7]);
     * A.PLU().E.toString();
     *   // "[0.0000 0.0000  1.0000 0.0000]
     *   //  [1.0000 0.0000  0.0000 0.0000]
     *   //  [0.8333 0.0000  0.5000 1.0000]
     *   //  [0.0000 1.0000 -0.4000 0.2000]"
     * A.rank();  // 3
     * A.leftNullBasis().map(vec => vec.toString());
     *   // ["[0.0000 1.0000 -0.4000 0.2000]"]
     *
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting.
     * @return {Vector[]} A basis for the left null space of the matrix.
     * @see Matrix#PLU
     */
    leftNullBasis(ε=1e-10) {
        if(this._cache.leftNullBasis) return this._cache.leftNullBasis;
        let {E} = this.PLU(ε);
        let r = this.rank(ε);
        let ret = new Array(this.m - r);
        for(let i = r; i < this.m; ++i)
            ret[i - r] = E[i].clone();
        this._cache.leftNullBasis = ret;
        return this._cache.leftNullBasis;
    }

    /**
     * @summary
     * Return the left null space of the matrix.
     *
     * @example {@lang javascript}
     * let A = Matrix.create([ 0, -3, -6,  4,  9],
     *                       [-1, -2, -1,  3,  1],
     *                       [-2, -3,  0,  3, -1],
     *                       [ 1,  4,  5, -9, -7]);
     * A.leftNullBasis().map(vec => vec.toString());
     *   // ["[0.0000 1.0000 -0.4000 0.2000]"]
     * A.leftNullSpace().toString(1);
     *   // "Subspace of R^4 of dimension 1 with basis
     *   //  [ 0.0]
     *   //  [ 1.0]
     *   //  [-0.4]
     *   //  [ 0.2]"
     *
     * @desc
     * This is essentially a shortcut for
     * `new Subspace(this.leftNullBasis(ε))`.
     *
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting.
     * @return {Subspace} The left null space of the matrix.
     */
    leftNullSpace(ε=1e-10) {
        if(this._cache.leftNullSpace) return this._cache.leftNullSpace;
        this._cache.leftNullSpace = new Subspace(
            this.leftNullBasis(ε), {n: this.m, isBasis: true});
        return this._cache.leftNullSpace;
    }

    // Orthogonality

    /**
     * @summary
     * The data computed in the Gram&ndash;Schmidt algorithm.
     *
     * @desc
     * This primarily consists of an `m`x`n` matrix `Q` with orthogonal
     * columns and an upper-triangular `n`x`n` matrix `R` such that `A = QR`.
     * Also included is a list of which columns of `Q` are zero.
     *
     * @typedef QRData
     * @type {object}
     * @property {Matrix} Q - An `m`x`n` matrix with orthogonal columns.
     * @property {Matrix} R - An `n`x`n` upper-triangular matrix.
     * @property {number[]} LD - A list of the zero columns of `Q`.
     */

    /**
     * @summary
     * Compute a QR decomposition.
     *
     * @desc
     * This computes a matrix `Q` and an upper-triangular matrix `R` such that
     * `A = QR`.  If `A` has full column rank then `Q` has orthonormal columns
     * and `R` is invertible, and the decomposition is unique.  Otherwise `Q`
     * may have zero columns and `R` may have zero rows.  In any case the
     * nonzero columns of `Q` form an orthonormal basis for the column space of
     * `A`.
     *
     * This method implements the modified Gram&ndash;Schmidt algorithm, which
     * is a simple tweak of the usual Gram&ndash;Schmidt algorithm with better
     * numerical properties.  It runs in about `O(n^3)` time for an `n`x`n`
     * matrix.  (The [Householder algorithm]{@link
     * https://en.wikipedia.org/wiki/Householder_transformation} is faster and
     * is the one generally used in practice.)
     *
     * As a side-effect, this method also computes the rank of the matrix.
     *
     * @example {@lang javascript}
     * let A = Matrix.create([ 3, -5,  1],
     *                       [ 1,  1,  1],
     *                       [-1,  5, -2],
     *                       [ 3, -7,  8]);
     * let {Q, R} = A.QR();
     * Q.toSring();
     *   // "[ 0.6708  0.2236 -0.6708]
     *   //  [ 0.2236  0.6708  0.2236]
     *   //  [-0.2236  0.6708  0.2236]
     *   //  [ 0.6708 -0.2236  0.6708]"
     * // Q has orthogonal columns
     * Q.transpose.mult(Q).toString(1);
     *   // "[1.0 0.0 0.0]
     *   //  [0.0 1.0 0.0]
     *   //  [0.0 0.0 1.0]"
     * Q.colSpace().equals(A.colSpace());  // true
     * R.toString();
     *   // "[4.4721 -8.9443  6.7082]
     *   //  [0.0000  4.4721 -2.2361]
     *   //  [0.0000  0.0000  4.4721]"
     * R.isUpperTri();                     // true
     * Q.mult(R).toString(1);
     *   // "[ 3.0 -5.0  1.0]
     *   //  [ 1.0  1.0  1.0]
     *   //  [-1.0  5.0 -2.0]
     *   //  [ 3.0 -7.0  8.0]"
     *
     * @example {@lang javascript}
     * // This matrix has rank 3
     * let A = Matrix.create([ 0, -3, -6,  4,  9],
     *                       [-1, -2, -1,  3,  1],
     *                       [-2, -3,  0,  3, -1],
     *                       [ 1,  4,  5, -9, -7]);
     * let {Q, R, LD} = A.QR();
     * Q.toString();
     *   // "[ 0.0000 -0.8018 0.0000 -0.5976 0.0000]
     *   //  [-0.4082  0.0000 0.0000  0.0000 0.0000]
     *   //  [-0.8165  0.2673 0.0000 -0.3586 0.0000]
     *   //  [ 0.4082  0.5345 0.0000 -0.7171 0.0000]"
     * // Columns 3 and 5 are zero: that means the third column was in the span
     * // of the first 2, and the fifth was in the span of the first 4.
     * LD;  // [2, 4]
     * // The nonzero columns form an orthonormal basis of the column space of A.
     * Q.colSpace().equals(A.colSpace());  // true
     * R.toString();
     *   // "[2.4495 4.8990 2.4495 -7.3485  -2.4495]
     *   //  [0.0000 3.7417 7.4833 -7.2161 -11.2250]
     *   //  [0.0000 0.0000 0.0000  0.0000   0.0000]
     *   //  [0.0000 0.0000 0.0000  2.9881   0.0000]
     *   //  [0.0000 0.0000 0.0000  0.0000   0.0000]"
     * Q.mult(R).toString(1);
     *   // "[ 0.0 -3.0 -6.0  4.0  9.0]
     *   //  [-1.0 -2.0 -1.0  3.0  1.0]
     *   //  [-2.0 -3.0  0.0  3.0 -1.0]
     *   //  [ 1.0  4.0  5.0 -9.0 -7.0]"
     *
     * @param {number} [ε=1e-10] - Vectors smaller than this value are taken
     *   to be zero.
     * @return {QRData} The QR factorization of the matrix.
     * @see Matrix#isFullColRank
     */
    QR(ε=1e-10) {
        if(this._cache.QR) return this._cache.QR;
        let {m, n} = this;
        let ui = new Array(n), LD=[];
        let vi = Array.from(this.cols());
        let R = Matrix.zero(n);
        for(let j = 0; j < n; ++j) {
            let u = vi[j].clone();
            // Compute ui[j] and the jth column of R at the same time
            for(let jj = 0; jj < j; ++jj) {
                // Modified Gram--Schmidt: ui[jj].dot(u) instead of
                // ui[jj].dot(vi[j]).
                let factor = ui[jj].dot(u);
                u.sub(ui[jj].clone().scale(factor));
                R[jj][j] = factor;
            }
            let l = u.size;
            if(l > ε) {
                R[j][j] = l;
                ui[j] = u.scale(1/l);
            } else {
                R[j][j] = 0;
                ui[j] = Vector.zero(m);
                LD.push(j);
            }
        }
        let Q = Matrix.from(ui).transpose;
        if(LD.length === 0)
            Q.hint({hasONCols: true});
        this._cache.QR = {Q, R, LD};
        this._cache.rank = n - LD.length;
        return this._cache.QR;
    }

    // Eigenstuff

    /**
     * @summary
     * Compute the (real and complex) eigenvalues of the matrix.
     *
     * @desc
     * The eigenvalues are the numbers `λ` such that there exists a nonzero
     * vector `v` with `Av = λv`.  They are the roots of the characteristic
     * polynomial.
     *
     * This method factors the characteristic polynomial and returns the roots.
     * This is not how eigenvalues are generally computed in practice.
     *
     * [Factorization]{@link Polynomial#factor} of polynomials is only
     * implemented in degrees at most 4; for matrices larger than 4x4, compute
     * the eigenvalues some other way and use {@link Matrix#hint}.
     *
     * @example {@lang javascript}
     * Matrix.create([1, 1], [1, 1]).eigenvalues();  // [[0, 1], [2, 1]]
     * Matrix.create([1, 1], [0, 1]).eigenvalues();  // [[1, 2]]
     * Matrix.create([1, 1], [-1,1]).eigenvalues();
     *   // [[new Complex(1, 1), 1], [new Complex(1, -1), 1]]
     *
     * @param {number} [ε=1e-10] - Rounding factor to determine multiplicities,
     *   as in {@link Polynomial#factor}.
     * @return {Root[]} The eigenvalues with algebraic multiplicity.  They are
     *   returned in the order specified in {@link Polynomial#factor}.
     * @throws Will throw an error if the matrix is not square, or if it is
     *   larger than 4x4 and the eigenvalues have not been
     *   [hinted]{@link Matrix#hint}.
     * @see Matrix#charpoly
     * @see Matrix#hint
     * @see Polynomial#factor
     */
    eigenvalues(ε=1e-10) {
        if(this._cache.eigenvalues) return this._cache.eigenvalues;
        if(!this.isSquare())
            throw new Error("Tried to compute the eigenvalues of a non-square matrix");
        this._cache.eigenvalues = this.charpoly.factor(ε);
        return this._cache.eigenvalues;
    }

    /**
     * @summary
     * Compute the `λ`-eigenspace of the matrix.
     *
     * @desc
     * The `λ`-eigenspace is the null space of `A - λI`.
     *
     * This method works for any size matrix if you know an eigenvalue.  For
     * complex eigenvalues, it returns a basis for the eigenspace represented as
     * an Array of pairs of vectors `[v_r, v_i]`, where `v_r + i v_i` is the
     * eigenvector.
     *
     * @example {@lang javascript}
     * let A = Matrix.create([1, 1], [1, 1]);
     * A.eigenspace(0).toString(1);
     *   // "Subspace of R^2 of dimension 1 with basis
     *   //  [-1.0]
     *   //  [ 1.0]"
     * A.eigenspace(2).toString(1);
     *   // "Subspace of R^2 of dimension 1 with basis
     *   //  [1.0]
     *   //  [1.0]"
     *
     * @example {@lang javascript}
     * let A = Matrix.create([2, 0], [0, 2]);
     * A.eigenspace(2).toString();  // "The full subspace R^2"
     *
     * @example {@lang javascript}
     * let A = Matrix.create([1, 1], [0, 1]);
     * A.eigenspace(1).toString(1);
     *   // "Subspace of R^2 of dimension 1 with basis
     *   //  [1.0]
     *   //  [0.0]"
     *
     * @example {@lang javascript}
     * let A = Matrix.create([24, -53/2], [20, -22]);
     * // The (1+i)-eigenspace is spanned by (1.15 + 0.05i, 1).
     * A.eigenspace(new Complex(1, 1));
     *   // [Vector.create(1.15, 1), Vector.create(0.05, 0)]
     * A.eigenspace(new Complex(1, -1));
     *   // [Vector.create(1.15, 1), Vector.create(-0.05, 0)]
     *
     * @param {number} λ - The eigenvalue.
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting and rounding.
     * @return {(Subspace|Array.<Vector[]>)} The eigenspace.
     * @throws Will throw an error if the matrix is not square, or if `λ` is not
     *   an eigenvalue of the matrix.
     */
    eigenspace(λ, ε=1e-10) {
        if(λ instanceof Complex) {
            if(Math.abs(λ.Im) > ε)
                return this._complexEigenspace(λ, ε);
            λ = λ.Re;
        }
        return this._realEigenspace(λ, ε);
    }

    _realEigenspace(λ, ε=1e-10) {
        if(!this._cache.eigenspaces) this._cache.eigenspaces = new Map();
        let closest = Infinity, best;
        // Find best matching eigenvalue
        for(let [λ1, V] of this._cache.eigenspaces.entries()) {
            let c = Math.abs(λ1 - λ);
            if(c < closest && c <= ε) {
                closest = c;
                best = V;
            }
        }
        if(best) return best;
        // Compute the eigenspace
        let V = this.clone().sub(Matrix.identity(this.n, λ)).nullSpace(ε);
        if(V.dim == 0)
            throw new Error("λ is not an eigenvalue of this matrix");
        this._cache.eigenspaces.set(λ, V);
        return V;
    }

    _complexEigenspace(λ, ε=1e-10) {
        if(!this._cache.cplxEigenspaces) this._cache.cplxEigenspaces = new Map();
        let closest = Infinity, best;
        // Find best matching eigenvalue
        for(let [λ1, V] of this._cache.cplxEigenspaces.entries()) {
            let c = λ1.clone().sub(λ).sizesq;
            if(c < closest && c <= ε*ε) {
                closest = c;
                best = V;
            }
        }
        if(best) return best;

        // The row reduction algorithm in PLU() won't work for complex
        // matrices.  We implement a simplified version here.
        let {m, n} = this;
        let pivots = [];
        let U = Array.from(this, row => Array.from(row, x => new Complex(x)));
        for(let i = 0; i < n; ++i) U[i][i].sub(λ);

        for(let curRow = 0, curCol = 0; curRow < m && curCol < n; ++curCol) {
            // Find maximal pivot
            let pivot = U[curRow][curCol], row = curRow;
            for(let i = curRow+1; i < m; ++i) {
                if(U[i][curCol].sizesq > pivot.sizesq) {
                    pivot = U[i][curCol];
                    row = i;
                }
            }
            if(pivot.sizesq > ε*ε) {
                // curCol is a pivot column
                [U[row], U[curRow]] = [U[curRow], U[row]];
                // Eliminate
                for(let i = curRow+1; i < m; ++i) {
                    let l = U[i][curCol].clone().div(pivot).mult(-1);
                    U[i][curCol].Re = U[i][curCol].Im = 0;
                    for(let j = curCol+1; j < n; ++j)
                        U[i][j].add(U[curRow][j].clone().mult(l));
                }
                pivots.push([curRow, curCol]);
                curRow++;

            } else {
                // Clear the column so U is really upper-triangular
                for(let i = curRow; i < m; ++i)
                    U[i][curCol].Re = U[i][curCol].Im = 0;
            }
        }

        if(pivots.length == n)
            throw new Error("λ is not an eigenvalue of this matrix");

        // Transform into rref
        for(let k = pivots.length-1; k >= 0; --k) {
            let [row, col] = pivots[k];
            let pivot = U[row][col];
            for(let i = 0; i < row; ++i) {
                let l = U[i][col].clone().div(pivot).mult(-1);
                for(let j = col+1; j < n; ++j)
                    U[i][j].add(U[row][j].clone().mult(l));
                U[i][col].Re = U[i][col].Im = 0;
            }
            let l = pivot.recip();
            for(let j = col+1; j < n; ++j)
                U[row][j].mult(l);
            U[row][col].Re = 1;
            U[row][col].Im = 0;
        }

        // Now extract the null basisa
        let basis = [], previous = [];
        for(let j = 0; j < n; ++j) {
            if(pivots.length && pivots[0][1] === j) {
                // Pivot column
                previous.push(pivots.shift());
                continue;
            }
            // Free column
            let Re_v = Vector.zero(n), Im_v = Vector.zero(n);
            for(let [row, col] of previous) {
                Re_v[col] = -U[row][j].Re;
                Im_v[col] = -U[row][j].Im;
            }
            Re_v[j] = 1;
            basis.push([Re_v, Im_v]);
        }

        this._cache.cplxEigenspaces.set(λ.clone(), basis);
        return basis;
    }

    /**
     * @summary
     * Diagonalizability data: `this = CDC^(-1)`.
     *
     * @typedef Diagonalization
     * @type {object}
     * @property {Matrix} C - An `n`x`n` invertible matrix.
     * @property {Matrix} D - An `n`x`n` diagonal matrix, or a block diagonal
     *   matrix in the case of block diagonalization.
     */

    /**
     * @summary
     * Diagonalize the matrix.
     *
     * @desc
     * For usual diagonalization (`opts.block` is falsy), it returns an
     * invertible matrix `C` and a diagonal matrix `D` such that `A = CDC^(-1)`.
     * This is possible exactly when `A` has `n` linearly independent real
     * eigenvectors, which form the columns of `C`; the eigenvalues are the
     * diagonal entries of `D`.
     *
     * For block diagonalization (`opts.block` is truthy), it returns an
     * invertible matrix `C` and a matrix `D` with diagonal blocks consisting of
     * numbers and rotation-scaling matrices, such that `A = CDC^(-1)`.  This is
     * possible exactly when `A` has `n` linearly independent (real and complex)
     * eigenvectors.  If all eigenvalues are real, then diagonalization is the
     * same as block diagonalization.
     *
     * If `opts.ortho` is truthy, compute orthonormal eigenbases for all real
     * eigenvalues.  Then if (and only if) `A` is symmetric, the matrix `C` will
     * be orthogonal, thus resulting in an orthogonal decomposition `A = CDC^T`.
     *
     * This is only implemented for matrices up to 4x4, unless the eigenvalues
     * have been [hinted]{@link Matrix#hint}.
     *
     * @example {@lang javascript}
     * let A = Matrix.create([11/13, 22/39,  2/39],
     *                       [-4/13, 83/39,  4/39],
     *                       [-1/13, 11/39, 40/39]);
     * A.eigenvalues();  // [[1, 2], [2, 1]]  (approximately)
     * let {C, D} = A.diagonalize();
     * D.toString(1);
     *   // "[1.0 0.0 0.0]
     *   //  [0.0 1.0 0.0]
     *   //  [0.0 0.0 2.0]"
     * C.toString();
     *   // "[3.6667 0.3333 2.0000]
     *   //  [1.0000 0.0000 4.0000]
     *   //  [0.0000 1.0000 1.0000]"
     * C.mult(D).mult(C.inverse()).scale(39).toString(1);
     *   // "[ 33.0 22.0  2.0]
     *   //  [-12.0 83.0  4.0]
     *   //  [ -3.0 11.0 40.0]"
     *
     * @example {@lang javascript}
     * let A = Matrix.create([    1,   1/2,     0],
     *                       [-4/13, 83/39,  4/39],
     *                       [ 5/13,  7/78, 34/39]);
     * A.eigenvalues();  // [[1, 2], [2, 1]]  (approximately)
     * A.diagonalize();  // null
     *
     * @example {@lang javascript}
     * let A = Matrix.create([33/29, -23/29,   9/29],
     *                       [22/29,  33/29, -23/29],
     *                       [19/29,  14/29,  50/29]);
     * A.eigenvalues();
     *   // [[new Complex(1, 1), 1], [new Complex(1, -1), 1], [2, 1]]
     *   // (approximately)
     * A.diagonalize();  // null
     * let {C, D} = A.diagonalize({block: true});
     * D.toString(1);
     *   // "[ 1.0 1.0 0.0]
     *   //  [-1.0 1.0 0.0]
     *   //  [ 0.0 0.0 2.0]"
     * C.toString();
     *   // "[-1.4000 0.2000  0.6667]
     *   //  [ 0.4000 1.8000 -0.3333]
     *   //  [ 1.0000 0.0000  1.0000]"
     * C.mult(D).mult(C.inverse()).scale(29).toString(1);
     *   // "[33.0 -23.0   9.0]
     *   //  [22.0  33.0 -23.0]
     *   //  [19.0  14.0  50.0]"
     *
     * @example {@lang javascript}
     * let A = Matrix.create([-1,  5,  4],
     *                       [ 5, -2,  3],
     *                       [ 4,  3, -9]);
     * A.isSymmetric();  // true
     * let {C, D} = A.diagonalize({ortho: true});
     * C.toString();
     *   // "[-0.3105  0.6337 0.7085]
     *   //  [-0.1442 -0.7681 0.6239]
     *   //  [ 0.9396  0.0916 0.3299]"
     * C.isOrthogonal();  // true
     * C.mult(D).mult(C.transpose).equals(A, 1e-10);  // true
     *
     * @param {Object} [opts={}] - Options.
     * @param {boolean} [opts.block=false] - Perform block diagonalization.
     * @param {boolean} [opts.ortho=false] - Use orthonormal bases for real
     *   eigenspaces.
     * @param {number} [opts.ε=1e-10] - Rounding factor.
     * @return {?Diagonalization} The diagonalization, or `null` if the matrix
     *   is not (block) diagonalizable.
     * @throws Will throw an error if the matrix is not square or if the matrix is
     *   larger than 4x4 and the eigenvalues have not been
     *   [hinted]{@link Matrix#hint}.
     * @see Matrix#eigenvalues
     * @see Matrix#hint
     * @see Matrix#eigenspace
     */
    diagonalize({block=false, ortho=false, ε=1e-10}={}) {
        let eigenbasis = [];
        let D = Matrix.zero(this.n);
        let i = 0;
        // Only use one of a conjugate pair of eigenvalues
        for(let [λ,m] of this.eigenvalues(ε).filter(
            ([λ,]) => !(λ instanceof Complex) || λ.Im >= 0)) {
            if(λ instanceof Complex) {
                if(!block) return null;
                let B = this.eigenspace(λ, ε);
                if(B.length < m) // Impossible for matrices <= 3x3
                    return null;
                for(let j = 0; j < m; ++j, i += 2) {
                    // The columns are the real and complex parts of the eigenvectors
                    eigenbasis.push(...B[j]);
                    D.insertSubmatrix(i, i, Matrix.create([λ.Re, λ.Im], [-λ.Im, λ.Re]));
                }
            } else {
                let V = this.eigenspace(λ, ε);
                if(V.dim < m)
                    return null;
                let B = ortho ? V.ONbasis(ε) : V.basis;
                for(let j = 0; j < m; ++j, ++i) {
                    D[i][i] = λ;
                    eigenbasis[i] = B.col(j);
                }
            }
        }
        let C = Matrix.from(eigenbasis).transpose;
        return {C, D};
    }

    /**
     * @summary
     * Singular value decomposition data.
     *
     * @typedef SVDData
     * @type {object}
     * @property {Matrix} U - An `m`x`m` orthogonal matrix.
     * @property {Matrix} V - An `n`x`n` orthogonal matrix.
     * @property {number[]} Σ - The singular values.
     */

    /**
     * @summary
     * Singular Value decomposition.
     *
     * @desc
     * This computes an orthogonal `m`x`m` matrix `U`, an orthogonal `n`x`n`
     * matrix `V`, and a diagonal `m`x`n` matrix `Σ`, such that `A = U Σ V^T`.
     * The nonzero diagonal entries of `Σ` are the singular values of the
     * matrix, in descending order.  The first `r = rank(A)` columns of `V` form
     * an orthonormal basis for the row space of `A`, and the last `n-r` columns
     * form an orthonormal basis for the null space.  The first `r` columns of
     * `U` form an orthonormal basis for the column space of `A`, and the last
     * `m-r` columns form an orthonormal basis for the left null space.  If `vi`
     * is the `i`th column of `V`, with `i <= r`, then `Avi = σi ui`, where
     * `σi` is the `i`th singular value and `ui` is the `i`th column of `U`.
     *
     * Only the nonzero diagonal entries of the matrix `Σ` are returned.  Use
     * {@link Matrix.diagonal} to turn it into a matrix, as in
     * `Matrix.diagonal(Σ, m, n)`.
     *
     * This method will fail if `min(m, n) > 4`, as the eigenvalue computations
     * have not been implemented for matrices larger that 4x4.  This method
     * implements the naïve schoolbook algorithm for computing the SVD, which is
     * not numerically accurate and is rarely used in practice.
     *
     * @example {@lang javascript}
     * let A = Matrix.create([ 0, -3, -6,  4,  9],
     *                       [-1, -2, -1,  3,  1],
     *                       [-2, -3,  0,  3, -1],
     *                       [ 1,  4,  5, -9, -7]);
     * let {U, V, Σ} = A.SVD();
     * U.toString();
     *   // "[-0.6428  0.5630 0.5194 -0.0000]
     *   //  [-0.1951 -0.3367 0.1235  0.9129]
     *   //  [-0.1224 -0.6974 0.6045 -0.3651]
     *   //  [ 0.7306  0.2887 0.5912  0.1826]"
     * U.isOrthogonal();  // true
     * V.toString();
     *   // "[ 0.0660  0.3409 -0.4064  0.0000  0.8452]
     *   //  [ 0.3162  0.3765 -0.6874 -0.1690 -0.5071]
     *   //  [ 0.4344 -0.2697 -0.1557  0.8452  0.0000]
     *   //  [-0.5694 -0.5819 -0.5806  0.0000 -0.0000]
     *   //  [-0.6186  0.5750  0.0304  0.5071 -0.1690]"
     * V.isOrthogonal();  // true
     * // The number of singular values is the rank.
     * Σ;  //  [17.73589205992891, 5.925426481232183, 1.824130985994107]
     * U.mult(Matrix.diagonal(Σ, 4, 5)).mult(V.transpose).equals(A, 1e-10);  // true
     * let vi = Array.from(V.cols());
     * new Subspace(vi.slice(0, 3)).equals(A.rowSpace());    // true
     * new Subspace(vi.slice(3)).equals(A.nullSpace());      // true
     * let ui = Array.from(U.cols());
     * new Subspace(ui.slice(0, 3)).equals(A.colSpace());    // true
     * new Subspace(ui.slice(3)).equals(A.leftNullSpace());  // true
     *
     * @param {number} [ε=1e-10] - Rounding factor used when computing
     *   eigenvalues and eigenvectors of the normal matrix.
     * @return {SVDData} The singular value decomposition.
     * @throws Will throw an error if the `min(m, n) > 4`, or if the eigenvalue
     *   / eigenspace computations otherwise fail numerically.
     */
    SVD(ε=1e-10) {
        if(this._cache.SVDData) return this._cache.SVDData;
        // If the matrix is wide, then the normal matrix of the transpose is
        // smaller, so compute the SVD of the transpose instead.
        let {m, n} = this;
        if(m < n) {
            let {U, V, Σ} = this.transpose.SVD(ε);
            this._cache.SVDData = {U: V, V: U, Σ};
            return this._cache.SVDData;
        }
        let ATA = this.normal;
        let eigenvals = ATA.eigenvalues(ε);
        let Σ = [], ui = [], vi = [];
        for(let i = eigenvals.length - 1; i >= 0; --i) {
            let [λ, m] = eigenvals[i];
            if(λ instanceof Complex)
                throw new Error("Eigenvalue computation failed for normal matrix");
            λ = Math.max(λ, 0);  // Paranoia
            let σ = Math.sqrt(λ);
            let Q = ATA.eigenspace(λ, ε).ONbasis(ε);
            if(Q.n < m)
                throw new Error("Eigenspace computation failed for normal matrix");
            for(let v of Q.cols()) {
                vi.push(v);
                if(λ > ε) {
                    // Singular value
                    Σ.push(σ);
                    ui.push(this.apply(v).scale(1/σ));
                }
            }
        }
        // Now we have computed V and the singular values, as well as the first
        // `r` columns of U, which form a basis for the column space of A.  It
        // remains to compute an orthonormal basis of the left null space.
        ui.push(...this.leftNullSpace(ε).ONbasis(ε).cols());
        let U = Matrix.from(ui).transpose;
        let V = Matrix.from(vi).transpose;
        U.hint({hasONCols: true});
        V.hint({hasONCols: true});
        this._cache.SVDData = { U, V, Σ };
        return this._cache.SVDData;
    }

    /**
     * @summary
     * LDLT decomposition data.
     *
     * @typedef LDLTData
     * @type {object}
     * @property {Matrix} L - An `n`x`n` lower-unitriangular matrix.
     * @property {Matrix} D - The diagonal entries.
     */

    /**
     * @summary
     * Compute the LDLT decomposition of a symmetric matrix.
     *
     * @desc
     * For a symmetric matrix `A`, this computes a lower-unitriangular matrix
     * `L` and a diagonal matrix `D` with nonzero entries such that `A = L D
     * L^T`, if possible.  This should be seen as a variant of Gaussian
     * elimination that takes advantage of the symmetry of `A`.  It uses about
     * `n^3/3` operations, which is about twice as fast as Gaussian
     * elimination.
     *
     * The LDLT decomposition exists for symmetric, positive-definite matrices
     * (e.g. the normal matrix of a matrix with full column rank).  It exists
     * for some indefinite matrices, and does not exist for singular matrices.
     *
     * Only the diagonal entries of the matrix `D` are returned.  Use
     * {@link Matrix.diagonal} to turn it into a matrix, as in
     * `Matrix.diagonal(D)`.
     *
     * @example {@lang javascript}
     * let A = Matrix.create([3, 2, 3],
     *                       [2, 7, 4],
     *                       [3, 4, 8]);
     * A.isPosDef();  // true
     * let {L, D} = A.LDLT();
     * L.toString();
     *   // "[1.0000 0.0000 0.0000]
     *   //  [0.6667 1.0000 0.0000]
     *   //  [1.0000 0.3529 1.0000]"
     * let DM = Matrix.diagonal(D);
     * DM.toString();
     *   // "[3.0000 0.0000 0.0000]
     *   //  [0.0000 5.6667 0.0000]
     *   //  [0.0000 0.0000 4.2941]"
     * L.mult(DM).mult(L.transpose).toString();
     *   // "[3.0000 2.0000 3.0000]
     *   //  [2.0000 7.0000 4.0000]
     *   //  [3.0000 4.0000 8.0000]"
     *
     * @example {@lang javascript}
     * let A = Matrix.create([1, 2, 3],
     *                       [2, 5, 4],
     *                       [3, 4, 2]);
     * A.isPosDef();  // false
     * // Not positive-definite, but the LDLT decomposition still exists
     * let {L, D} = A.LDLT();
     * L.toString(0);
     *   // "[1  0 0]
     *   //  [2  1 0]
     *   //  [3 -2 1]"
     * let DM = Matrix.diagonal(D);
     * DM.toString(0);
     *   // "[1 0   0]
     *   //  [0 1   0]
     *   //  [0 0 -11]"
     * L.mult(DM).mult(L.transpose).equals(A);  // true
     *
     * @example {@lang javascript}
     * let A = Matrix.create([ 27, -27,  18],
     *                       [-27,  45,  12],
     *                       [ 18,  12,  62]);
     * A.isSingular();  // true
     * A.LDLT();        // null
     *
     * @param {number} [ε=1e-10] - Entries smaller than this are considered to
     *   be zero.
     * @return {?LDLTData} The LDLT decomposition, or `null` if it does not exist.
     * @throws Will throw an error if the matrix is not symmetric.
     * @see Matrix#cholesky
     */
    LDLT(ε=1e-10) {
        if(this._cache.LDLT !== undefined) return this._cache.LDLT;
        if(!this.isSymmetric(ε))
            throw new Error("Tried to compute an LDLT decomposition of a non-symmetric matrix");
        let n = this.n;
        let D = new Array(n);
        let L = Matrix.identity(n);
        for(let j = 0; j < n; ++j) {
            D[j] = this[j][j];
            for(let k = 0; k < j; ++k)
                D[j] -= L[j][k]*L[j][k]*D[k];
            if(Math.abs(D[j]) <= ε) {
                this._cache.LDLT = null;
                return null;
            }
            for(let i = j+1; i < n; ++i) {
                let l = this[i][j];
                for(let k = 0; k < j; ++k)
                    l -= L[i][k]*L[j][k]*D[k];
                L[i][j] = l / D[j];
            }
        }
        L.hint({isLowerUni: true});
        L.transpose.hint({isUpperUni: true});
        this._cache.LDLT = {L, D};
        return this._cache.LDLT;
    }

    /**
     * @summary
     * Compute the [Cholesky decomposition]{@link https://en.wikipedia.org/wiki/Cholesky_decomposition} of a positive-definite symmetric matrix.
     *
     * @desc
     * For a positive-definite symmetric matrix `A`, this computes a
     * lower-triangular matrix `L` with positive diagonal entries such that `A =
     * L L^T`.  These matrices are obtained from an [LDLT decomposition]{@link
     * Matrix#LDLT} by multiplying the columns of `L` by the square roots of the
     * corresponding entries of `D`:
     * > `A = L D L^T = (L D^(1/2)) (L D^(1/2))^T`
     *
     * @example {@lang javascript}
     * let A = Matrix.create([3, 2, 3],
     *                       [2, 7, 4],
     *                       [3, 4, 8]);
     * A.isPosDef();  // true
     * let L = A.cholesky();
     * L.toString();
     *   // "[1.7321 0.0000 0.0000]
     *   //  [1.1547 2.3805 0.0000]
     *   //  [1.7321 0.8402 2.0722]"
     * L.mult(L.transpose).equals(A, 1e-10);  // true
     *
     * @example {@lang javascript}
     * let A = Matrix.create([1, 2, 3],
     *                       [2, 5, 4],
     *                       [3, 4, 2]);
     * A.isPosDef();  // false
     * A.cholesky();  // null
     *
     * @param {number} [ε=1e-10] - Entries smaller than this are considered to
     *   be zero.
     * @return {?Matrix} The matrix `L` such that `A = L L^T`, or `null` if the
     *   matrix is not positive-definite.
     * @throws Will throw an error if the matrix is not symmetric.
     * @see Matrix#LDLT
     */
    cholesky(ε=1e-10) {
        let LDLT = this.LDLT(ε);
        if(!LDLT) return null;
        let {L, D} = LDLT;
        L = L.clone();
        for(let [i, d] of D.entries()) {
            if(d < 0) return null;
            d = Math.sqrt(d);
            for(let j = i; j < this.n; ++j)
                L[j][i] *= d;
        }
        return L;
    }
};


export default Matrix;
