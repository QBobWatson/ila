/** @module lib/matrix
 *
 * @file
 * Implements a Matrix class containing most algorithms from basic linear
 * algebra.
 */

'use strict';

import Vector from "./vector.js";
import Complex from "./complex.js";
import Subspace from "./subspace.js";

//import { quadratic, cardano } from "./linalg.js";


/**
 * Class representing a matrix.
 *
 * A matrix is a 2 dimensional array of numbers.  The first dimension indexes
 * the row, written horizontally; the second indexes the column.
 *
 * The rows of a `Matrix` are `Vector` instances.  All rows must have the same
 * length.
 *
 * @extends Array
 */
class Matrix extends Array {
    /**
     * Create an `n`x`n` identity Matrix.
     *
     * @example
     * Matrix.identity(3).equals(new Matrix([1,0,0],[0,1,0],[0,0,1])) == true;
     *
     * @param {integer} n - The resulting Matrix will have this many rows and columns.
     * @param {number} [λ=1] - The diagonal entries will be equal to `λ`.
     * @return {Matrix}
     */
    static identity(n, λ=1) {
        let rows = new Array(n);
        for(let i = 0; i < n; ++i)
            rows[i] = Vector.e(i, n, λ);
        return new Matrix(...rows);
    }

    /**
     * Create an `m`x`n` zero Matrix.
     *
     * @example
     * Matrix.zero(2, 3).equals(new Matrix([0,0,0],[0,0,0])) == true;
     *
     * @param {integer} m - The resulting Matrix will have this many rows.
     * @param {integer} [n=m] - The resulting Matrix will have this many columns.
     * @return {Matrix}
     */
    static zero(m, n=m) {
        let rows = new Array(m);
        for(let i = 0; i < m; ++i)
            rows[i] = Vector.zero(n);
        return new Matrix(...rows);
    }

    /**
     * Create a permutation Matrix.
     *
     * This is the matrix whose `i`th row is the `j`th unit coordinate
     * vector.
     *
     * @example
     * Matrix.permutation(1, 0).equals(new Matrix([0,1],[1,0])) == true;
     *
     * @param {number[]} vals - If `n` numbers are given, these should be a
     *   permutation of the set {0,1,...,n-1}, although no error is thrown if
     *   each argument is between 0 and n-1.
     * @return {Matrix}
     */
    static permutation(vals) {
        let ret = Matrix.zero(vals.length);
        for(let i = 0; i < vals.length; ++i)
            ret[i][vals[i]] = 1;
        return ret;
    }

    /**
     * Create a Matrix.
     *
     * @param {...(Array<number>|Vector)} rows - The rows of the matrix.  Arrays
     *   will be promoted to Vector instances.  All rows must have the same
     *   length.
     * @throws Will throw an error if the rows do not have the same length.
     */
    constructor(...rows) {
        if(rows.length === 1 && typeof rows[0] === "number") {
            // Need to support Array(n) semantics so this.map() works
            super(rows[0]);
            return;
        }
        if(rows.length === 0) {
            super();
            return;
        }
        super(...rows.map(
            row => row instanceof Vector ? row : Vector.create(...row)));
        let n = this[0].length;
        if(this.some(row => row.length !== n))
            throw new Error("Matrix rows must have the same length.");
    }

    /**
     * Invalidate cached computations.
     *
     * Call this after modifying any matrix entries.  It is not called
     * automatically for performance reasons.
     *
     * @return {undefined}
     */
    invalidate() {
        if(this.__cache) delete this.__cache;
    }

    /**
     * @private
     *
     * @type {object}
     * @property {Matrix} transpose - Transpose matrix.
     * @property {integer} rank - The rank of the matrix.
     * @property {PLUData} PLU - PLU factorization; computed in PLU().
     * @property {Matrix} rref - Reduced row echelon form.
     * @property {Matrix} E - Matrix such that `E*this = rref`.
     * @property {Vector[]} nullBasis - Basis for the null space.
     * @property {Subspace} nullSpace - The null space.
     * @property {Vector[]} colBasis - Basis for the column space.
     * @property {Subspace} nullSpace - The column space.
     * @property {Vector[]} rowBasis - Basis for the row space.
     * @property {Subspace} rowSpace - The row space.
     * @property {Vector[]} leftNullBasis - Basis for the left null space.
     * @property {Subspace} leftNullSpace - The left null space.
     * @property {QRData} QR - QR factorization; computed in QR().
     * @property {number} det - The determinant.
     * @property {number[]} charpoly - Characteristic polynomial.
     * @property {Root[]} eigenvalues - Eigenvalues.
     * @property {Map.<number, Subspace>} eigenspaces - Real eigenspaces.
     * @property {Map.<Complex, Array.<Complex[]>>} cplxEigenspaces - Complex
     *   eigenspaces.
     * @property {Matrix} normal - The normal matrix A^TA.
     */
    get _cache() {
        if(!this.__cache) this.__cache = {};
        return this.__cache;
    }

    /**
     * Insert a submatrix into the matrix at position `(i,j)`.
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
     * The number of rows of the matrix.
     *
     * @type {integer}
     */
    get m() { return this.length; }

    /**
     * The number of columns of the matrix.
     *
     * @type {integer}
     */
    get n() { return this.length === 0 ? 0 : this[0].length; }

    /**
     * Whether the matrix is square.
     *
     * @return {boolean}
     */
    isSquare() {
        return this.m == this.n;
    }

    /**
     * Check if this matrix is equal to `other`.
     *
     * Two matrices are equal if they have the same dimensions, and all
     * entries are equal.
     *
     * @param {Matrix} other - The matrix to compare.
     * @param {number} [ε=0] - Entries will test as equal if they are within `ε`
     *   of each other.  This is provided in order to account for rounding
     *   errors.
     * @return {boolean}
     */
    equals(other, ε=0) {
        if(this.m !== other.m || this.n !== other.n)
            return false;
        return this.every((v, i) => v.equals(other[i], ε));
    }

    /**
     * Create a new Matrix with the same entries.
     *
     * @return {Matrix}
     */
    clone() {
        // this.map() calls the Matrix() constructor
        return this.map(row => row.clone());
    }

    /**
     * Return a string representation of the matrix.
     *
     * @param {integer} [precision=4] - The number of decimal places to include.
     * @return {string}
     */
    toString(precision=4) {
        // this.map returns a Matrix
        let strings = new Array(this.m), i = 0, j;
        for(let row of this) {
            let rowStrs = new Array(row.length);
            j = 0;
            for(let v of row)
                rowStrs[j++] = v.toFixed(precision);
            strings[i++] = rowStrs;
        }
        let colLens = new Array(this.n);
        for(j = 0; j < this.n; ++j)
            colLens[j] = Math.max(...strings.map(row => row[j].length));
        let ret = new Array(strings.length);
        i = 0;
        for(let row of strings)
            ret[i++] = row.map((val, j) => val.padStart(colLens[j], ' ')).join(' ');
        return ret.join('\n');
    }

    /**
     * Return the `i`th row.
     *
     * @param {integer} i - The row to return.
     * @return {Vector} The `i`th row of `this`.
     */
    row(i) {
        return this[i];
    }

    /**
     * Return an iterable over the rows.
     *
     * @return {Object}
     */
    rows() {
        return this[Symbol.iterator]();
    }

    /**
     * Return the `j`th column.
     *
     * @param {integer} j - The column to return.
     * @return {Vector} The `j`th column of `this`.
     */
    col(j) {
        let v = Vector.zero(this.m), i = 0;
        for(let row of this)
            v[i++] = row[j];
        return v;
    }

    /**
     * Return an iterable over the columns.
     *
     * @return {Object}
     */
    cols() {
        let self = this;
        return (function*() {
            for(let j = 0; j < self.n; ++j)
                yield self.col(j);
        })();
    }

    /**
     * Add a Matrix in-place.
     *
     * This modifies the matrix in-place by adding the entries of `other`.
     *
     * @param {Matrix} other - The matrix to add.
     * @param {number} [factor=1] - Add `factor` times `other` instead of just
     *   adding `other`.
     * @return {Matrix} `this`
     * @throws Will throw an error if the matrices have different sizes.
     */
    add(other, factor=1) {
        if(this.m !== other.m || this.n !== other.n)
            throw new Error(
                'Tried to add matrices of different sizes');
        this.forEach((row, i) => row.add(other[i], factor));
        return this;
    }

    /**
     * Subtract a Matrix in-place.
     *
     * This modifies the matrix in-place by subtracting the entries of `other`.
     *
     * @param {Matrix} other - The matrix to add.
     * @return {Matrix} `this`
     * @throws Will throw an error if the matrices have different sizes.
     */
    sub(other) {
        return this.add(other, -1);
    }

    /**
     * Scale a Matrix in-place.
     *
     * This modifies the matrix in-place by multiplying all entries by `c`.
     *
     * @param {number} c - The number to multiply.
     * @return {Matrix} `this`
     */
    scale(c) {
        for(let row of this)
            row.scale(c);
        return this;
    }

    /**
     * The transpose Matrix.
     *
     * The rows of the transpose are the columns of the matrix.
     *
     * @type {Matrix}
     */
    get transpose() {
        if(this._cache.transpose) return this._cache.transpose;
        this._cache.transpose = new Matrix(this.n);
        for(let j = 0; j < this.n; ++j)
            this._cache.transpose[j] = this.col(j);
        return this._cache.transpose;
    }

    /**
     * The normal matrix `A^TA`.
     *
     * @type {Matrix}
     */
    get normal() {
        if(this._cache.normal) return this._cache.normal;
        this._cache.normal = this.transpose.mult(this);
        return this._cache.normal;
    }

    /**
     * Multiply by another matrix.
     *
     * This creates a new matrix equal to the product of `this` and `other`.
     * This only makes sense when the number of rows of `other` equals the
     * number of columns of `this`.
     *
     * @param {Matrix} other - The matrix to multiply.
     * @return {Matrix} The matrix product.
     * @throws Will throw an error if the matrices have incompatible
     * dimensions.
     */
    mult(other) {
        if(other.m !== this.n)
            throw new Error('Cannot multiply matrices of incompatible dimensions');
        return this.map(
            row => {
                let newRow = Vector.zero(other.n);
                for(let i = 0; i < other.n; ++i) {
                    for(let j = 0; j < this.n; ++j)
                        newRow[i] += row[j] * other[j][i];
                }
                return newRow;
            });
    }

    /**
     * Multiply a matrix times a vector.
     *
     * This creates a new vector equal to the product of `this` and `v`.  This
     * only makes sense when the number of entries of `v` equals the number of
     * columns of `this`.
     *
     * This is an optimized special case of {@link Matrix#mult}.
     *
     * @param {Vector} v - The vector to multiply.
     * @return {Vector} The matrix-vector product.
     * @throws Will throw an error if the matrix and vector have incompatible
     *   dimensions.
     */
    apply(v) {
        if(v.length !== this.n)
            throw new Error('Cannot multiply matrix and vector of incompatible dimensions');
        let ret = Vector.zero(this.m), i = 0;
        for(let row of this)
            ret[i++] = row.dot(v);
        return ret;
    }

    /**
     * Whether the matrix is zero.
     *
     * @param {number} [ε=0] - Entries smaller than this value are assumed
     *   to be zero.
     * @return {boolean}
     */
    isZero(ε=0) {
        return this.every(row => row.isZero(ε));
    }

    /**
     * Whether the matrix is upper-triangular.
     *
     * This means that all entries below the main diagonal are equal to zero.
     *
     * @param {number} [ε=0] - Entries smaller than this value are assumed
     *   to be zero.
     * @return {boolean}
     */
    isUpperTri(ε=0) {
        for(let i = 1; i < this.m; ++i) {
            for(let j = 0; j < i; ++j) {
                if(Math.abs(this[i][j]) > ε)
                    return false;
            }
        }
        return true;
    }

    /**
     * Whether the matrix is upper-triangular with ones on the diagonal.
     *
     * @param {number} [ε=0] - Entries smaller than this value are assumed
     *   to be zero.
     * @return {boolean}
     */
    isUpperUnip(ε=0) {
        for(let i = 1; i < Math.min(this.m, this.n); ++i) {
            if(Math.abs(this[i][i] - 1) > ε)
                return false;
        }
        return this.isUpperTri(ε);
    }

    /**
     * Whether the matrix is lower-triangular.
     *
     * This means that all entries above the main diagonal are equal to zero.
     *
     * @param {number} [ε=0] - Entries smaller than this value are assumed
     *   to be zero.
     * @return {boolean}
     */
    isLowerTri(ε=0) {
        for(let i = 0; i < this.m; ++i) {
            for(let j = i+1; j < this.n; ++j) {
                if(Math.abs(this[i][j]) > ε)
                    return false;
            }
        }
        return true;
    }

    /**
     * Whether the matrix is lower-triangular with ones on the diagonal.
     *
     * @param {number} [ε=0] - Entries smaller than this value are assumed
     *   to be zero.
     * @return {boolean}
     */
    isLowerUnip(ε=0) {
        for(let i = 1; i < Math.min(this.m, this.n); ++i) {
            if(Math.abs(this[i][i] - 1) > ε)
                return false;
        }
        return this.isLowerTri(ε);
    }

    /**
     * Whether the matrix is diagonal.
     *
     * This means the matrix is both upper- and lower-triangular.
     *
     * @param {number} [ε=0] - Entries smaller than this value are assumed
     *   to be zero.
     * @return {boolean}
     */
    isDiagonal(ε=0) {
        return this.isLowerTri(ε) && this.isUpperTri(ε);
    }


    /**
     * A pivot position is recorded as a pair `[row, column]`.
     *
     * @typedef {number[]} Pivot
     */

    /**
     * A list of leading entries of each row.
     *
     * Zero rows do not have a leading entry, and are not included.
     *
     * @param {number} [ε=0] - Entries smaller than this value are assumed
     *   to be zero.
     * @return {Pivot[]}
     */
    leadingEntries(ε=0) {
        let entries = [];
        for(let i = 0; i < this.m; ++i) {
            for(let j = 0; j < this.n; ++j) {
                if(Math.abs(this[i][j]) > ε) {
                    entries.push([i, j]);
                    break;
                }
            }
        }
        return entries;
    }

    /**
     * Whether the matrix is in row-echelon form.
     *
     * This means that the first nonzero entry of each row is to the right of
     * the first nonzero entry of the previous row, which implies that all
     * entries below a pivot are zero and all zero rows are at the bottom.
     *
     * @param {number} [ε=0] - Entries smaller than this value are assumed
     *   to be zero.
     * @return {boolean}
     */
    isEchelon(ε=0) {
        let i0 = -1;
        let j0 = -1;
        for(let [i, j] of this.leadingEntries(ε)) {
            if(i !== i0 + 1 || j <= j0)
                return false;
            i0 = i;
            j0 = j;
        }
        return true;
    }

    /**
     * Whether the matrix is in reduced row-echelon form.
     *
     * This means that it is in row-echelon form, and in addition, all pivots
     * are equal to one, and all entries above a pivot are equal to zero.
     *
     * @param {number} [ε=0] - Entries smaller than this value are assumed
     *   to be zero.
     * @return {boolean}
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
     * The sum of the diagonal entries of the matrix.
     *
     * @type {number}
     */
    get trace() {
        let acc = 0;
        for(let i = 0; i < Math.min(this.m, this.n); ++i)
            acc += this[i][i];
        return acc;
    }

    // Row operations

    /**
     * Scale a row by a constant.
     *
     * This modifies the matrix in-place by multiplying all entries of row `i`
     * by the scalar `c`.
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
     * Add a constant times one row to another row.
     *
     * This modifies the matrix in-place by adding all entries of row `i2` times
     * the scalar `c` to the entries of row `i1`.
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
     * Swap two rows.
     *
     * This exchanges rows `i1` and `i2`.
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
     * The core data computed by Gaussian elimination.
     *
     * @typedef PLUData
     * @type {object}
     * @property {number[]} P - A permutation of the numbers `1...m-1`.
     * @property {number} detP - The sign of the permutation P.
     * @property {Matrix} L - An `m`x`m` lower-triangular matrix with ones on
     *   the diagonal.
     * @property {Matrix} U - An `m`x`n` matrix in row echelon form.
     * @property {Matrix} E - An `m`x`m` invertible matrix equal to `L^(-1)P`,
     *   so `EA = U`.
     * @property {Pivot[]} pivots - An array of pivot positions.
     */

    /**
     * PA = LU factorization.
     *
     * This computes an `m`x`m` permutation matrix `P`, an `m`x`m`
     * lower-triangular matrix `L` with ones on the diagonal, and a matrix `U`
     * in row echelon form, such that `PA = LU`.
     *
     * This is the core method that implements Gaussian elimination.  It uses
     * maximal partial pivoting.  Many properties of the matrix are derived from
     * the output of this method, which is cached.
     *
     * The permutation matrix `P` is returned as a list of `n` numbers defining
     * the permutation.  Use {@link permutation} to turn it into a matrix.
     *
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting.
     * @return {PLUData}
     */
    PLU(ε=1e-10) {
        if(this._cache.PLU) return this._cache.PLU;

        let P = new Array(this.m);
        for(let i = 0; i < this.m; ++i) P[i] = i;
        let L = Matrix.identity(this.m);
        let E = Matrix.identity(this.m);
        let U = this.clone();
        let {m, n} = this;
        let pivots = [];
        let detP = 1;

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
                    detP *= -1;
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
                curRow++;

            } else {
                // Clear the column so U is really upper-triangular
                for(let i = curRow; i < m; ++i)
                    U[i][curCol] = 0;
            }
        }

        this._cache.PLU = {P, L, U, E, pivots, detP};
        return this._cache.PLU;
    }

    /**
     * Find some solution `x` to the equation `Ax=b`.
     *
     * This is the solution obtained by forward- and back-substituting in the
     * `PA=LU` decomposition, taking the free variables equal to zero.  All
     * other solutions can be obtained by adding elements of the null space.
     *
     * @param {Vector} b - The right-hand side of the equation.
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting.
     * @return {(Vector|null)} Returns `null` if no solution exists.
     * @throws Will throw an error if `b.length != this.m`.
     */
    solve(b, ε=1e-10) {
        if(b.length != this.m)
            throw new Error("Incompatible dimensions of matrix and vector");
        let {P, L, U, pivots} = this.PLU(ε);
        let r = this.rank(ε);
        let Pb = Vector.zero(this.m);
        // Solve LUx = PAx = Pb
        for(let i = 0; i < this.m; ++i) Pb[i] = b[P[i]];
        // Solve Ly = Pb by forward-substitution
        let y = Pb.clone();
        for(let i = 1; i < this.m; ++i) {
            for(let ii = 0; ii < i; ++ii)
                y[i] -= L[i][ii] * y[ii];
        }
        // Check if a solution exists
        for(let i = r; i < this.m; ++i) {
            if(Math.abs(y[i]) > ε)
                return null;
        }
        // Solve Ux = y by reverse-substitution
        let x = Vector.zero(this.n);
        for(let p = r-1; p >= 0; --p) {
            let [row, col] = pivots[p];
            x[col] = y[row];
            for(let pp = p+1; pp < r; ++pp) {
                let [, col1] = pivots[pp];
                x[col] -= U[row][col1] * x[col1];
            }
            x[col] /= U[row][col];
        }
        return x;
    }

    /**
     * Find the shortest solution `x` to the equation `Ax=b`.
     *
     * This is obtained by finding some solution using {@link solve}, then
     * projecting onto the row space.
     *
     * @param {Vector} b - The right-hand side of the equation.
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting and projecting.
     * @return {(Vector|null)} Returns `null` if no solution exists.
     * @throws Will throw an error if `b.length != this.m`.
     */
    solveShortest(b, ε=1e-10) {
        let x = this.solve(b, ε);
        if(x === null) return null;
        return this.rowSpace().project(x, ε);
    }

    /**
     * Find some least-squares solution `x` to the equation `Ax=b`.
     *
     * This is obtained by solving the normal equation `A^TAx = A^Tb`.
     *
     * Use {@link pseudoinverse} if you have to run this routine many times.
     *
     * @param {Vector} b - The right-hand side of the equation.
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting and projecting.
     * @return {Vector} A least-squares solution.
     * @throws Will throw an error if `b.length != this.m`.
     */
    solveLeastSquares(b, ε=1e-10) {
        let ATA = this.normal;
        let x = ATA.solve(this.transpose.apply(b), ε);
        return x;
    }

    /**
     * Find the shortest least-squares solution `x` to the equation `Ax=b`.
     *
     * This is obtained by finding some solution using
     * {@link solveLeastSquares}, then projecting onto the row space.
     *
     * Use {@link pseudoinverse} if you have to run this routine many times.
     *
     * @param {Vector} b - The right-hand side of the equation.
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting and projecting.
     * @return {Vector} The shortest least-squares solution.
     * @throws Will throw an error if `b.length != this.m`.
     */
    solveLeastSquaresShortest(b, ε=1e-10) {
        let ATA = this.normal;
        let x = ATA.solve(this.transpose.apply(b), ε);
        return this.rowSpace().project(x, ε);
    }

    /**
     * Compute the pivot positions of the matrix.
     *
     * This requires Gaussian elimination to compute.
     *
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting.
     * @return {Pivot[]}
     */
    pivots(ε=1e-10) {
        let {pivots} = this.PLU(ε);
        return pivots;
    }

    /**
     * Compute the rank of the matrix.
     *
     * This is computed as a side-effect of other computations.
     *
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting and projecting.
     * @return {integer}
     */
    rank(ε=1e-10) {
        if(this._cache.rank === undefined) // Hasn't been computed yet
            this._cache.rank = this.pivots(ε).length;
        return this._cache.rank;
    }

    /**
     * Compute the nullity of the matrix.
     *
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting and projecting.
     * @return {integer}
     */
    nullity(ε=1e-10) {
        return this.n - this.rank(ε);
    }

    /**
     * The reduced row echelon form of the matrix.
     *
     * This completes the Gauss--Jordan elimination on the row echelon form of
     * the matrix, by scaling so all pivots are equal to 1 and performing row
     * replacements so that all entries above a pivot are equal to zero.
     *
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting.
     * @return {Matrix}
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
     * Return the row operations used in Gauss--Jordan elimination.
     *
     * This returns an invertible `m`x`m` matrix `E` such that `EA = rref`.
     *
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting.
     * @return {Matrix}
     */
    rowOps(ε=1e-10) {
        if(this._cache.rowOps) return this._cache.rowOps;
        this.rref(ε);
        return this._cache.rowOps;
    }

    // Bases

    /**
     * Basis for the null space of the matrix.
     *
     * This is computed using the reduced row echelon form and parametric vector
     * form.
     *
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting.
     * @return {Vector[]}
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
     * The null space of the matrix.
     *
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting.
     * @return {Subspace}
     */
    nullSpace(ε=1e-10) {
        if(this._cache.nullSpace) return this._cache.nullSpace;
        this._cache.nullSpace = new Subspace(
            this.nullBasis(ε), {n: this.n, isBasis: true});
        return this._cache.nullSpace;
    }

    /**
     * Basis for the column space of the matrix.
     *
     * This is composed of the pivot columns.
     *
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting.
     * @return {Vector[]}
     */
    colBasis(ε=1e-10) {
        if(this._cache.colBasis) return this._cache.colBasis;
        let {pivots} = this.PLU(ε);
        this._cache.colBasis = pivots.map(([, j]) => this.col(j));
        return this._cache.colBasis;
    }

    /**
     * The column space of the matrix.
     *
     * @return {Subspace}
     */
    colSpace() {
        if(this._cache.colSpace) return this._cache.colSpace;
        this._cache.colSpace = new Subspace(this, {n: this.m});
        return this._cache.colSpace;
    }

    /**
     * Basis for the row space of the matrix.
     *
     * This is composed of the nonzero rows of a REF.
     *
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting.
     * @return {Vector[]}
     */
    rowBasis(ε=1e-10) {
        if(this._cache.rowBasis) return this._cache.rowBasis;
        let {pivots, U} = this.PLU(ε);
        this._cache.rowBasis = pivots.map(([i,]) => U[i].clone());
        return this._cache.rowBasis;
    }

    /**
     * The row space of the matrix.
     *
     * @return {Subspace}
     */
    rowSpace() {
        if(this._cache.rowSpace) return this._cache.rowSpace;
        this._cache.rowSpace = new Subspace(this.transpose, {n: this.n});
        return this._cache.rowSpace;
    }

    /**
     * Basis for the left null space of the matrix.
     *
     * This is composed of the last `m-r` rows of `L^-1 P`, where `P` and `L`
     * come from the `PA = LU` decomposition and `r` is the rank.
     *
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting.
     * @return {Vector[]}
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
     * The left null space of the matrix.
     *
     * @return {Subspace}
     */
    leftNullSpace(ε=1e-10) {
        if(this._cache.leftNullSpace) return this._cache.leftNullSpace;
        this._cache.leftNullSpace = new Subspace(
            this.leftNullBasis(ε), {n: this.m, isBasis: true});
        return this._cache.leftNullSpace;
    }

    /**
     * Test whether the matrix has full row rank.
     *
     * This means that there is a pivot in every row.
     *
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting and projecting.
     * @return {boolean}
     */
    isFullRowRank(ε=1e-10) {
        if(this.m > this.n) return false; // Don't need to row reduce
        return this.rank(ε) == this.m;
    }

    /**
     * Test whether the matrix has full column rank.
     *
     * This means that there is a pivot in every column.
     *
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting and projecting.
     * @return {boolean}
     */
    isFullColRank(ε=1e-10) {
        if(this.m < this.n) return false; // Don't need to row reduce
        return this.rank(ε) == this.n;
    }

    /**
     * Test whether the matrix is invertible.
     *
     * This means that the matrix is square and has the maximum number of
     * pivots.  Equivalently, there is another matrix (namely `this.inverse()`)
     * such that the product with `this` on both sides is equal to the identity.
     *
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting and projecting.
     * @return {boolean}
     */
    isInvertible(ε=1e-10) {
        return this.isFullRowRank(ε) && this.isFullColRank(ε);
    }

    /**
     * Alias for `!this.isInvertible()`.
     *
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting and projecting.
     * @return {boolean}
     */
    isSingular(ε=1e-10) {
        return !this.isInvertible(ε);
    }

    /**
     * Compute the inverse matrix.
     *
     * This only makes sense for square matrices.  Throws an error if the matrix
     * is not square or is not invertible.
     *
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting.
     * @return {Matrix}
     * @throws Will throw an error if the matrix is not square or not
     *   invertible.
     */
    inverse(ε=1e-10) {
        if(!this.isSquare())
            throw new Error("Tried to invert a non-square matrix");
        let E = this.rowOps(ε);
        if(!this.isInvertible(ε))
            throw new Error("Tried to invert a singular matrix");
        return E;
    }

    /**
     * Compute the matrix determinant.
     *
     * This only makes sense for square matrices.  Throws an error if the matrix
     * is not square.
     *
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting.
     * @return {number}
     * @throws Will throw an error if the matrix is not square.
     */
    det(ε=1e-10) {
        if(this._cache.det !== undefined) return this._cache.det;
        if(!this.isSquare())
            throw new Error("Tried to take the determinant of a non-square matrix");
        let {U, detP} = this.PLU(ε);
        let det = 1;
        for(let i = 0; i < this.n; ++i)
            det *= U[i][i];
        this._cache.det = det * detP;
        return this._cache.det;
    }

    // Orthogonality

    /**
     * Check if the matrix is orthogonal.
     *
     * That means that the matrix is square and has orthonormal columns.
     *
     * @param {number} [ε=1e-10] - Numbers smaller than this value are taken
     *   to be zero.
     * @return {QRData}
     */
    isOrthogonal(ε=1e-10) {
        return this.isSquare() &&
            this.transpose.mult(this).equals(Matrix.identity(this.n), ε);
    }

    /**
     * The data computed Gram--Schmidt.
     *
     * @typedef QRData
     * @type {object}
     * @property {Matrix} Q - An `m`x`n` matrix with orthogonal columns.
     * @property {Matrix} R - An `n`x`n` upper-triangular matrix.
     * @property {number[]} LD - A list of the zero columns of `Q`.
     */

    /**
     * QR decomposition.
     *
     * This computes a matrix `Q` and an upper-triangular matrix `R` such that
     * `this = QR`.  If `this` has full column rank then `Q` has orthonormal
     * columns and `R` is invertible.  Otherwise `Q` may have zero columns and
     * `R` may have zero rows; the nonzero columns of `Q` are orthonormal.
     *
     * This implements the Gram--Schmidt algorithm, which is not a very
     * efficient or numerically accurate.
     *
     * @param {number} [ε=1e-10] - Vectors smaller than this value are taken
     *   to be zero.
     * @return {QRData}
     */
    QR(ε=1e-10) {
        if(this._cache.QR) return this._cache.QR;
        let {m, n} = this;
        let ui = new Array(n), LD=[];
        let vi = [...this.cols()];
        let R = Matrix.zero(n);
        for(let j = 0; j < n; ++j) {
            let u = vi[j].clone();
            // Compute ui[j] and the jth column of R at the same time
            for(let jj = 0; jj < j; ++jj) {
                let factor = ui[jj].dot(vi[j]);
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
        let Q = new Matrix(...ui).transpose;
        this._cache.QR = {Q, R, LD};
        if(this._cache.rank === undefined)
            this._cache.rank = n - LD.length;
        return this._cache.QR;
    }

    // Eigenstuff

    /**
     * Compute the characteristic polynomial.
     *
     * The characteristic polynomial of an `n`x`n` matrix `A` is the determinant
     * of `(A - λI_n)`, where λ is an indeterminate.  This is a polynomial of
     * degree `n` and leading coefficient `(-1)^n`.  The polynomial
     *      (-1)^n λ^n + c1 λ^(n-1) + c2 λ^(n-2) + ... + c(n-1) λ + cn
     * is represented as an array of numbers `[c1, c2, c3, ..., cn]`.
     *
     * This only makes sense for square matrices.  Throws an error if the matrix
     * is not square.
     *
     * This algorithm uses a recursive formula for the characteristic polynomial
     * in terms of traces of powers of the matrix that can be found in:
     *      R. R. Silva, J. Math. Phys. 39, 6206-6213 (1998)
     * I make no claim that this is the most efficient algorithm -- but it is
     * easy to implement.
     *
     * @return {number[]} The coefficients of λ^{n-1}, λ^{n-2}, ..., 1.
     * @throws Will throw an error if the matrix is not square.
     */
    charpoly() {
        if(this._cache.charpoly) return this._cache.charpoly;
        if(!this.isSquare())
            throw new Error("Tried to compute the characteristic polynomial of a non-square matrix");
        let n = this.n;
        let ret = new Array(n);
        let traces = new Array(n);
        let power = this;
        for(let i = 0; i < n; ++i) {
            traces[i] = power.trace;
            power = power.mult(this);
            ret[i] = traces[i];
            for(let j = 0; j < i; ++j)
                ret[i] += ret[j] * traces[i-j-1];
            ret[i] = -ret[i]/(i+1);
        }
        // This is det(λI - A); multiply by (-1)^n now
        if(n % 2 == 0)
            this._cache.charpoly = ret;
        else
            this._cache.charpoly = ret.map(x => -x);
        return this._cache.charpoly;
    }

    /**
     * Compute the (real and complex) eigenvalues of the matrix.
     *
     * These are the roots of the characteristic polynomial.
     *
     * Only implemented for matrices of size 1x1, 2x2, and 3x3.
     *
     * @param {number} [ε=1e-10] - Eigenvalues will be considered equal if the
     *   discriminant of the characteristic polynomial is smaller than this.
     * @return {Root[]} The eigenvalues with algebraic multiplicity.
     * @throws Will throw an error if the matrix is not square, or if the matrix
     *   is bigger than 3x3.
     */
    eigenvalues(ε=1e-10) {
        if(this._cache.eigenvalues) return this._cache.eigenvalues;
        if(!this.isSquare())
            throw new Error("Tried to compute the eigenvalues of a non-square matrix");
        switch(this.n) {
        case 1:
            this._cache.eigenvalues = [[this[0][0], 1]];
            break;
        case 2:
            this._cache.eigenvalues = quadratic(...this.charpoly(), ε);
            break;
        case 3:
            let [b, c, d] = this.charpoly();
            // replace λ by -λ and re-sort
            this._cache.eigenvalues = cardano(b, -c, d, ε).map(
                ([x, m]) => x instanceof Complex ? [x.mult(-1), m] : [-x, m])
                .sort(([a,], [b,]) => {
                    if(typeof a === "number" && typeof b === "number")
                        return a - b;
                    if(typeof a === "number") return -1;
                    if(typeof b === "number") return 1;
                    if(Math.abs(a.Re - b.Re) < ε) return a.Im - b.Im;
                    return a.Re - b.Re;
                });
            break;
        default:
            throw new Error("Eigenvalue computations are only implemented for matrices up to 3x3");
        }
        return this._cache.eigenvalues;
    }

    /**
     * Cache the (real and complex) eigenvalues of the matrix.
     *
     * This is so that (block) diagonalization can be performed for matrices
     * larger than 3x3.
     *
     * @param {...Root} eigenvalues - The eigenvalues with algebraic multiplicity.
     * @return {undefined}
     */
    hintEigenvalues(...eigenvalues) {
        this._cache.eigenvalues = eigenvalues;
    }

    /**
     * Compute the `λ`-eigenspace of the matrix.
     *
     * This works for any size matrix if you know an eigenvalue.  For complex
     * eigenvalues, it returns a basis for the eigenspace represented as an
     * Array of pairs of vectors `[v_Re, v_Im]`, where `v_Re + iv_Im` is the
     * eigenvector.
     *
     * @param {number} λ - The eigenvalue.
     * @param {number} [ε=1e-10] - Rounding factor.
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
            if(c < closest && c < ε) {
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
            if(c < closest && c < ε*ε) {
                closest = c;
                best = V;
            }
        }
        if(best) return best;

        // The row reduction algorithm in PLU() won't work for complex
        // matrices.  We implement a simplified version here.
        let {m, n} = this;
        let pivots = [];
        let U = [...this].map(row => [...row].map(x => new Complex(x)));
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
     * Diagonalizability data: `this = CDC^(-1)`.
     *
     * @typedef Diagonalization
     * @type {object}
     * @property {Matrix} C - An `n`x`n` invertible matrix.
     * @property {Matrix} D - An `n`x`n` diagonal matrix, or a block diagonal
     *   matrix in the case of block diagonalization.
     */

    /**
     * Diagonalize the matrix.
     *
     * This is only implemented for matrices up to 3x3, unless the eigenvalues
     * have been hinted with {@link hintEigenvalues}.
     *
     * For usual diagonalization, it returns an invertible matrix `C` and a
     * diagonal matrix `D` such that `this = CDC^(-1)`.  For block
     * diagonalization, it returns an invertible matrix `C` and a matrix `D`
     * with diagonal blocks consisting of numbers and rotation-scaling matrices,
     * such that `this = CDC^(-1)`.  If all eigenvalues are real, then
     * diagonalization is the same as block diagonalization.
     *
     * @param {Object} [opts={}] - Options.
     * @param {boolean} [opts.block=false] - Perform block diagonalization.
     * @param {boolean} [opts.ortho=false] - Use orthonormal bases for real
     *   eigenspaces.
     * @param {number} [opts.ε=1e-10] - Rounding factor.
     * @return {?Diagonalization} The diagonalization, or `null` if the matrix
     *   is not (block) diagonalizable.
     * @throws Will throw an error if the matrix is not square or if the matrix is
     *   larger than 3x3 and the eigenvalues have not been hinted.
     */
    diagonalize({block=false, ortho=false, ε=1e-10}={}) {
        let eigenbasis = new Array(this.n);
        let D = Matrix.zero(this.n);
        let i = 0;
        let seen = [];
        for(let [λ,m] of this.eigenvalues(ε)) {
            if(λ instanceof Complex) {
                if(!block) return null;
                if(seen.find(z => z.equals(λ)))
                    continue; // Only use one of a conjugate pair of eigenvalues
                let B = this.eigenspace(λ, ε);
                if(B.length < m) // Impossible for matrices <= 3x3
                    return null;
                for(let j = 0; j < m; ++j, i += 2) {
                    // The columns are the real and complex parts of the eigenvectors
                    eigenbasis[i  ] = B[j][0];
                    eigenbasis[i+1] = B[j][1];
                    D.insertSubmatrix(i, i, new Matrix([λ.Re, λ.Im], [-λ.Im, λ.Re]));
                }
                seen.push(λ.clone().conj());
            } else {
                let V = this.eigenspace(λ, ε);
                if(V.dim < m)
                    return null;
                let B = ortho ? V.ONbasis(ε) : V.basis(ε);
                for(let j = 0; j < m; ++j, ++i) {
                    D[i][i] = λ;
                    eigenbasis[i] = B.col(j);
                }
            }
        }
        let C = new Matrix(...eigenbasis).transpose;
        return {C, D};
    }

    /**
     * Test if the matrix is diagonalizable.
     *
     * This is only implemented for matrices up to 3x3, unless the eigenvalues
     * have been hinted with {@link hintEigenvalues}.
     *
     * @param {number} [ε=1e-10] - Rounding factor.
     * @return {boolean}
     * @throws Will throw an error if the matrix is not square or if the matrix is
     *   larger than 3x3 and the eigenvalues have not been hinted.
     */
    isDiagonalizable(ε=1e-10) {
        return !!this.diagonalize({ε});
    }

    /**
     * Singular value decomposition data.
     *
     * @typedef SVDData
     * @type {object}
     * @property {Matrix} U - An `m`x`m` orthogonal matrix.
     * @property {Matrix} V - An `n`x`n` orthogonal matrix.
     * @property {number[]} Σ - The singular values.
     */


    /**
     * Singular Value decomposition.
     *
     * This computes an orthogonal `m`x`m` matrix `U`, an orthogonal `n`x`n`
     * matrix `V`, and a diagonal `m`x`n` matrix `Σ`, such that
     * `this = U Σ V^T`.  The diagonal entries of `Σ` are the singular values of
     * the matrix.
     *
     * Only the nonzero diagonal entries of the matrix `Σ` are returned.  Use
     * {@link diagonal} to turn it into a matrix.
     *
     * @param {number} [ε=1e-10] - Rounding factor.
     * @return {SVDData}
     */
    SVD(ε=1e-10) {
        if(this._cache.SVDData) return this._cache.SVDData;
        return this._cache.SVDData;
    }
};


export default Matrix;
