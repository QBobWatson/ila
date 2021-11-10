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

/** @module matrix
 *
 * @file
 * Implements a Matrix class containing algorithms from basic linear algebra.
 */


import Vector from "./vector";
import Subspace from "./subspace";
import Polynomial from "./polynomial";
import { Root, MultRoot } from  "./polynomial";
import Complex from "./complex";

import { range } from "./util";


class MatrixError extends Error {
    constructor(message: string) {
        super(message);
        this.name = "MatrixError";
    }
}


type Pivot = [number, number];

export interface PLUData {
    /**
     * A permutation of the numbers `1...m-1` by {@link Matrix.permutation}.
     */
    P: number[];

    /**
     * An `m`x`m` lower-triangular matrix with ones on the diagonal.
     */
    L: SquareMatrix;

    /**
     * An `m`x`n` matrix in row echelon form.
     */
    U: Matrix;

    /**
     * An `m`x`m` invertible matrix equal to `L^(-1)P`, so `EA = U`.
     */
    E: SquareMatrix;

    /**
     * An array of pivot positions.
     */
    pivots: Pivot[];

    /**
     * The determinant of the matrix (only set for square matrices).
     */
    det?: number;
}

export interface JordanData extends PLUData {
    /**
     * An `m`x`m` invertible matrix which gives the reduced row echelon form of
     * `A` when left-multiplied by `A`.  When `A` is invertible, this is the
     * inverse of `A`.
     */
    rowOps: SquareMatrix;

    /**
     * The reduced row echelon form of `A`.
     */
    rref: Matrix;
}

export interface LDLTData {
    /**
     * An `n`x`n` lower-unitriangular matrix.
     */
    L: SquareMatrix;

    /**
     * The diagonal entries.
     */
    D: number[];
}

export interface QRData {
    /** An `m`x`n` matrix with orthogonal columns. */
    Q: Matrix;

    /** An `n`x`n` upper-triangular matrix. */
    R: SquareMatrix;

    /** A list of the zero columns of `Q` / zero rows of `R`. */
    LD: number[];
};

export interface DiagonalizationData {
    /** An `n`x`n` invertible matrix. */
    C: SquareMatrix;

    /** An `n`x`n` (block-)diagonal matrix. */
    D: SquareMatrix;
}

export interface SVDData {
    /** An `m`x`m` orthogonal matrix. */
    U: SquareMatrix;

    /** An `n`x`n` orthogonal matrix. */
    V: SquareMatrix;

    /** The singular values. */
    Σ: number[];
}


class Matrix {

    /***********************************************************************/
    /* Properties                                                          */
    /***********************************************************************/

    [index: number]: Vector;
    m: number;
    n: number;


    /***********************************************************************/
    /* Constructors                                                        */
    /***********************************************************************/

    static from<T>(entries: Iterable<T> | ArrayLike<T>,
                   mapfn?: (x: T, k: number) => number[] | Vector)
    : Matrix | SquareMatrix {
        mapfn = mapfn || ((x, _) => (x as any));
        return Matrix.create(...Array.from(entries, mapfn));
    }

    static create(...rows: (number[] | Vector)[]): Matrix | SquareMatrix {
        const row = rows[0];
        if(rows.length === (row ?
            (row instanceof Vector ? row.size : row.length) : 0))
            return new SquareMatrix(...rows);
        return new Matrix(...rows);
    }

    protected constructor(...rows: (number[] | Vector)[]) {
        if(rows.length === 0)
            throw new MatrixError("A matrix must have at least one row.")
        let rows2 = rows.map(r => r instanceof Vector
            ? r : new Vector(...r)) as Vector[];
        this.m = rows2.length;
        this.n = rows2[0]!.size;
        if(rows2.some(v => v.size !== this.n))
            throw new MatrixError("Matrix rows must have the same length.");
        for(let i = 0; i < this.m; ++i)
            this[i] = rows2[i];
    }

    static identity(n: number, λ: number=1): SquareMatrix {
        return Matrix.from(range(n), i => Vector.e(i, n, λ)) as SquareMatrix;
    }

    static zero(m: number, n: number=m): Matrix | SquareMatrix {
        return Matrix.from(range(m), _ => Vector.zero(n));
    }

    static constant(c: number, m: number, n: number=m): Matrix | SquareMatrix {
        return Matrix.from(range(m), _ => Vector.constant(n, c));
    }

    static diagonal(entries: number[], m: number=entries.length,
                    n: number=entries.length): Matrix | SquareMatrix {
        let ret = Matrix.zero(m, n);
        for(let i = 0; i < Math.min(m, n, entries.length); ++i)
            ret[i][i] = entries[i];
        return ret;
    }

    static permutation(vals: number[]): SquareMatrix {
        let n = vals.length;
        return Matrix.from(range(n), i => Vector.e(vals[i], n)) as SquareMatrix;
    }

    clone(): this {
        return new (this.constructor as any)(
            ...Array.from(this.rows, r => r.clone()));
    }


    /***********************************************************************/
    /* Accessors                                                           */
    /***********************************************************************/

    get transpose(): this {
        return new (this.constructor as any)(...this.cols);
    }

    get normal(): SquareMatrix {
        // TODO: This is not so efficient...
        return this.transpose.mult(this) as SquareMatrix;
    }

    get trace(): number {
        let acc = 0;
        for(const d of this.diag) acc += d;
        return acc;
    }

    row(i: number): Vector {
        return this[i];
    }

    get rows(): Iterable<Vector> {
        let self = this;
        return (function*() {
            for(let i = 0; i < self.m; ++i)
                yield self.row(i);
        })();
    }

    col(j: number): Vector {
        return Vector.from(this.rows, row => row[j]);
    }

    get cols(): Iterable<Vector> {
        let self = this;
        return (function*() {
            for(let j = 0; j < self.n; ++j)
                yield self.col(j);
        })();
    }

    get diag(): Iterable<number> {
        let self = this;
        return (function*() {
            for(let j = 0; j < Math.min(self.m, self.n); ++j)
                yield self[j][j];
        })();
    }


    /***********************************************************************/
    /* Utility                                                             */
    /***********************************************************************/

    toString(precision: number=4): string {
        let strings = Array.from(
            this.rows, row => Array.from(row, v => v.toFixed(precision)));
        let colLens = Array.from(
            range(this.n), j => Math.max(...strings.map(row => row[j].length)));
        return strings.map(
            row => '[' + row.map(
                (val, j) => val.padStart(colLens[j], ' ')).join(' ') + ']')
            .join('\n');
    }

    insertSubmatrix(i: number, j: number, M: Matrix): this {
        for(let ii = 0; ii < M.m; ++ii) {
            for(let jj = 0; jj < M.n; ++jj)
                this[ii+i][jj+j] = M[ii][jj];
        }
        return this;
    }

    equals(other: Matrix, ε: number=0): boolean {
        if(this.m !== other.m || this.n !== other.n)
            return false;
        for(let i = 0; i < this.m; ++i) {
            if(!this[i].equals(other[i], ε))
                return false;
        }
        return true;
    }

    leadingEntries(ε: number=0): Pivot[] {
        let entries: Pivot[] = [];
        for(let i = 0; i < this.m; ++i) {
            const row = this[i];
            for(let j = 0; j < this.n; ++j) {
                if(Math.abs(row[j]) > ε) {
                    entries.push([i, j]);
                    break;
                }
            }
        }
        return entries;
    }


    /***********************************************************************/
    /* Traits                                                              */
    /***********************************************************************/

    isSquare(): this is SquareMatrix {
        return this.m === this.n;
    }

    isZero(ε: number=0): boolean {
        for(let i = 0; i < this.m; ++i) {
            if(!this[i].isZero(ε))
                return false;
        }
        return true;
    }

    isUpperTri(ε: number=0): boolean {
        for(let i = 1; i < this.m; ++i) {
            for(let j = 0; j < i; ++j) {
                if(Math.abs(this[i][j]) > ε)
                    return false;
            }
        }
        return true;
    }

    isUpperUni(ε: number=0): boolean {
        for(let d of this.diag) {
            if(Math.abs(d - 1) > ε)
                return false;
        }
        return this.isUpperTri(ε);
    }

    isLowerTri(ε: number=0): boolean {
        for(let i = 0; i < this.m; ++i) {
            for(let j = i+1; j < this.n; ++j) {
                if(Math.abs(this[i][j]) > ε)
                    return false;
            }
        }
        return true;
    }

    isLowerUni(ε: number=0): boolean {
        for(let d of this.diag) {
            if(Math.abs(d - 1) > ε)
                return false;
        }
        return this.isLowerTri(ε);
    }

    isDiagonal(ε: number=0): boolean {
        return this.isLowerTri(ε) && this.isUpperTri(ε);
    }

    isREF(ε: number=0): boolean {
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

    isRREF(ε: number=0): boolean {
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

    hasONCols(ε: number=1e-10): boolean {
        return this.normal.equals(Matrix.identity(this.n), ε);
    }


    /***********************************************************************/
    /* Arithmetic                                                          */
    /***********************************************************************/

    add(other: Matrix, factor: number=1): this {
        if(this.m !== other.m || this.n !== other.n)
            throw new MatrixError('Tried to add matrices of different sizes');
        for(let i = 0; i < this.m; ++i)
            this[i].add(other[i], factor);
        return this;
    }

    sub(other: Matrix): this {
        return this.add(other, -1);
    }

    scale(c: number): this {
        for(let i = 0; i < this.m; ++i)
            this[i].scale(c);
        return this;
    }

    mult(other: Matrix): Matrix | SquareMatrix {
        if(other.m !== this.n)
            throw new MatrixError(
                'Cannot multiply matrices of incompatible dimensions');
        return Matrix.from(
            this.rows,
            row => Vector.from(
                range(other.n), i => [...row].reduce(
                    (a, v, j) => a + v * other[j][i], 0)));
    }

    apply(v: Vector): Vector {
        if(v.size !== this.n)
            throw new MatrixError(
                'Cannot multiply matrix and vector of incompatible dimensions');
        return Vector.from(this.rows, row => row.dot(v));
    }


    /***********************************************************************/
    /* Row Operations                                                      */
    /***********************************************************************/

    rowScale(i: number, c: number, start: number=0): this {
        this[i].scale(c, start);
        return this;
    }

    rowReplace(i1: number, i2: number, c: number, start: number=0): this {
        this[i1].add(this[i2], c, start);
        return this;
    }

    rowSwap(i1: number, i2: number): this {
        [this[i1], this[i2]] = [this[i2], this[i1]];
        return this;
    }


    /***********************************************************************/
    /* Gauss--Jordan Elimination                                           */
    /***********************************************************************/

    PLU(ε: number=1e-10): PLUData {
        let P = Array.from(range(this.m));
        let L = Matrix.identity(this.m);
        let E = Matrix.identity(this.m);
        let U = this.clone();
        let {m, n} = this;
        let pivots: Pivot[] = [];
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

        const ret: PLUData = {P, L, U, E, pivots};
        if(m === n)
            ret.det = pivots.length === n ? det * signP : 0;
        return ret;
    }

    jordanSubst(PLU: PLUData=this.PLU()): JordanData {
        let {U, E, pivots} = PLU;
        let rowOps = E.clone() as SquareMatrix;
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

        return {...PLU, rref, rowOps};
    }

    nullBasis(JD: JordanData=this.jordanSubst()): Vector[] {
        let { rref, pivots } = JD;
        pivots = pivots.slice();
        let previous: Pivot[] = [];
        let basis = [];
        for(let j = 0; j < this.n; ++j) {
            if(pivots.length && pivots[0][1] === j) {
                // Pivot column
                previous.push(pivots.shift()!);
                continue;
            }
            // Free column
            let v = Vector.zero(this.n);
            for(let [row, col] of previous)
                v[col] = -rref[row][j];
            v[j] = 1;
            basis.push(v);
        }
        return basis;
    }

    nullSpace(JD: JordanData=this.jordanSubst(), ε: number=1e-10): Subspace {
        return Subspace.spanOf(this.nullBasis(JD), {n: this.n, ε});
    }

    colBasis(PLU: PLUData=this.PLU()): Vector[] {
        return Array.from(PLU.pivots, ([, j]) => this.col(j));
    }

    colSpace(PLU: PLUData=this.PLU(), ε: number=1e-10): Subspace {
        return Subspace.spanOf(this.colBasis(PLU), {n: this.m, ε});
    }

    rowBasis(PLU: PLUData=this.PLU()): Vector[] {
        const { pivots, U } = PLU;
        return pivots.map(([i,]) => U[i].clone());
    }

    rowSpace(PLU: PLUData=this.PLU(), ε: number=1e-10): Subspace {
        return Subspace.spanOf(this.rowBasis(PLU), {n: this.n, ε});
    }

    leftNullBasis(PLU: PLUData=this.PLU()): Vector[] {
        const { E, pivots } = PLU;
        const r = pivots.length;
        let ret: Vector[] = new Array(this.m - r);
        for(let i = r; i < this.m; ++i)
            ret[i - r] = E[i].clone();
        return ret;
    }

    leftNullSpace(PLU: PLUData=this.PLU(), ε: number=1e-10): Subspace {
        return Subspace.spanOf(this.leftNullBasis(PLU), {n: this.m, ε});
    }

    solveReverseSubst(b: Vector, pivots?: Pivot[],
                      ε: number=1e-10): Vector | null {
        let x = Vector.zero(this.n);
        pivots = pivots || this.leadingEntries(ε);
        let r = pivots.length;
        // Check if a solution exists
        for(let i = r; i < this.m; ++i) {
            if(Math.abs(b[i]) > ε)
                return null;
        }
        for(let p = r - 1; p >= 0; --p) {
            let [row, col] = pivots[p];
            x[col] = b[row];
            for(let pp = p + 1; pp < r; ++pp) {
                let [, col1] = pivots[pp];
                x[col] -= this[row][col1] * x[col1];
            }
            x[col] /= this[row][col];
        }
        return x;
    }

    solve(b: Vector, hints?: {
        PLU?: PLUData,
        ε?: number}): Vector | null {

        if(b.size != this.m)
            throw new MatrixError("Incompatible dimensions of matrix and vector");

        hints = {ε: 1e-10, ...hints};
        const { ε, PLU } = hints;
        const { P, L, U, pivots } = PLU || this.PLU(ε);
        // Solve LUx = PAx = Pb
        const Pb = Vector.from(range(this.m), i => b[P[i]]);
        const y = L.solveForwardSubst(Pb);
        return U.solveReverseSubst(y, pivots, ε);
    }

    solveLeastSquares(b: Vector, hints?: {
        QR?: QRData,
        ε?: number}): Vector {

        hints = {ε: 1e-10, ...hints};
        const {ε, QR} = hints;

        if(QR) {
            // This is essentially SquareMatrix.solveReverseSubst1, except that
            // we ignore zero rows and columns.
            const { Q, R, LD } = QR;
            const QTb = Q.transpose.apply(b);
            let x = QTb.clone();
            for(let p = R.n - 1; p >= 0; --p) {
                if(p in LD) continue;
                for(let pp = p + 1; pp < R.n; ++pp) {
                    x[p] -= R[p][pp] * x[pp];
                }
                x[p] /= R[p][p];
            }
            return x;
        }

        return this.normal.solve(this.transpose.apply(b), {ε})!;
    }

    projectColSpace(b: Vector, hints?: {
        QR?: QRData, ε?: number}): Vector {
        return this.apply(this.solveLeastSquares(b, hints));
    }

    projectRowSpace(b: Vector, hints?: {
        transposeQR?: QRData, ε?: number}): Vector {
        if(!hints) hints = {};
        return this.transpose.projectColSpace(
            b, {QR: hints.transposeQR, ε: hints.ε});
    }

    solveShortest(b: Vector, hints?: {
        PLU?: PLUData, transposeQR?: QRData, ε?: number}): Vector | null {
        let x = this.solve(b, hints);
        if(x === null) return null;
        if(!hints) hints = {}
        return this.projectRowSpace(x, hints);
    }

    solveLeastSquaresShortest(b: Vector, hints?: {
        QR?: QRData, transposeQR?: QRData, ε?: number}): Vector {
        return this.projectRowSpace(this.solveLeastSquares(b, hints), hints);
    }

    pseudoInverse(hints?: {QR?: QRData, transposeQR?:
                           QRData, ε?: number}): Matrix | SquareMatrix {
        let { QR, transposeQR, ε } = {ε: 1e-10, ...hints};
        if(!QR) QR = this.QR(ε);
        if(!transposeQR) transposeQR = this.transpose.QR(ε);
        return Matrix.from(
            range(this.m),
            i => this.solveLeastSquaresShortest(
                Vector.e(i, this.m), { QR, transposeQR, ε })).transpose;
    }


    /***********************************************************************/
    /* Convenience                                                         */
    /***********************************************************************/

    pivots(PLU: PLUData=this.PLU()): Pivot[] {
        return PLU.pivots;
    }

    rank(PLU?: PLUData): number;
    rank(QR?: QRData): number;
    rank(x?: PLUData | QRData): number {
        if(x === undefined)
            x = this.PLU();
        if("P" in x)
            return x.pivots.length;
        return this.n - x.LD.length;
    }

    nullity(PLU: PLUData=this.PLU()): number {
        return this.n - PLU.pivots.length;
    }

    hasFullRowRank(PLU?: PLUData): boolean;
    hasFullRowRank(QR?: QRData): boolean;
    hasFullRowRank(x?: PLUData | QRData): boolean {
        if(x === undefined)
            x = this.PLU();
        let r = "P" in x ? this.rank(x) : this.rank(x);
        return r == this.m;
    }

    hasFullColRank(PLU?: PLUData): boolean;
    hasFullColRank(QR?: QRData): boolean;
    hasFullColRank(x?: PLUData | QRData): boolean {
        if(x === undefined)
            x = this.PLU();
        let r = "P" in x ? this.rank(x) : this.rank(x);
        return r == this.n;
    }

    rref(JD: JordanData=this.jordanSubst()): Matrix {
        return JD.rref;
    }


    /***********************************************************************/
    /* Orthogonality                                                       */
    /***********************************************************************/

    QR(ε: number=1e-10): QRData {
        let { m, n } = this;
        let ui: Vector[] = new Array(n), LD: number[]=[];
        let vi = Array.from(this.cols);
        let R = Matrix.zero(n) as SquareMatrix;
        for(let j = 0; j < n; ++j) {
            let u = vi[j].clone();
            // Compute ui[j] and the jth column of R at the same time
            for(let jj = 0; jj < j; ++jj) {
                // Modified Gram--Schmidt: ui[jj].dot(u) instead of
                // ui[jj].dot(vi[j]).
                let factor = ui[jj].dot(u);
                u.add(ui[jj], -factor);
                R[jj][j] = factor;
            }
            let l = u.len;
            if(l > ε) {
                R[j][j] = l;
                ui[j] = u.scale(1/l);
            } else {
                R[j][j] = 0;
                ui[j] = Vector.zero(m);
                LD.push(j);
            }
        }
        let Q = Matrix.create(...ui).transpose;
        return { Q, R, LD };
    }


    /***********************************************************************/
    /* SVD                                                                 */
    /***********************************************************************/

    SVD(hints?: {ε?: number, singularValues?: Root[]}): SVDData {
        // If the matrix is wide, then the normal matrix of the transpose is
        // smaller, so compute the SVD of the transpose instead.
        let { m, n } = this;
        if(m < n) {
            let { U, V, Σ } = this.transpose.SVD(hints);
            return { U: V, V: U, Σ };
        }
        let { ε, singularValues } = { ε: 1e-10, ...hints };
        let ATA = this.normal;
        let eigenvals: MultRoot[];
        if(singularValues) {
            eigenvals = singularValues.map<MultRoot>(
                x => Array.isArray(x) ? x : [x, 1]).map(
                    ([λ, m]) => {
                        if(typeof λ == "number") {
                            if(λ <= ε)
                                throw new MatrixError(
                                    "Singular values must be positive");
                            return [λ*λ, m];
                        }
                        throw new MatrixError(
                            "Singular values must be real")
                    });
            (eigenvals as [number, number][]).sort(
                ([λ1, _m1], [λ2, _m2]) => λ1-λ2);
            let sumams = eigenvals.reduce((acc, [_, m]) => acc + m, 0);
            if(sumams < n)
                eigenvals.push([0, n - sumams]);
        }
        else {
            eigenvals = ATA.eigenvalues(ε);
            let sumams = eigenvals.reduce((acc, [_, m]) => acc + m, 0);
            if(sumams < n)
                throw new MatrixError(
                    "Eigenvalue computation failed for normal matrix");
        }
        let Σ: number[] = [], ui: Vector[] = [], vi: Vector[] = [];
        for(let i = eigenvals.length - 1; i >= 0; --i) {
            let [λ, m] = eigenvals[i];
            if(λ instanceof Complex || λ < -ε)
                throw new MatrixError(
                    "Eigenvalue computation failed for normal matrix");
            let σ = Math.sqrt(λ);
            let B = ATA.eigenspaceBasis(λ, ε) as Vector[];
            if(B.length < m)
                throw new MatrixError(
                    "Eigenspace computation failed for normal matrix");
            let Q = B.length > 1
                ? [...Matrix.from(B).transpose.QR(ε).Q.cols]
                : [B[0].normalize()];
            for(let v of Q) {
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
        let r = ui.length;
        if(r < m) {
            let B = this.leftNullBasis(this.PLU(ε));
            if(B.length == 1)
                ui.push(B[0].normalize());
            else
                ui.push(...Matrix.from(B).transpose.QR(ε).Q.cols);
            if(ui.length < m)
                throw new MatrixError(
                    "Left null space computation failed");
        }
        let U = Matrix.from(ui).transpose as SquareMatrix;
        let V = Matrix.from(vi).transpose as SquareMatrix;
        return { U, V, Σ };
    }
}


class SquareMatrix extends Matrix {
    [x: string]: any;

    constructor(...rows: (number[] | Vector)[]) {
        super(...rows);
        if(this.m !== this.n)
            throw new MatrixError(
                "Tried to create a non-square SquareMatrix.");
    }


    /***********************************************************************/
    /* Traits                                                              */
    /***********************************************************************/

    isOrthogonal(ε: number=1e-10): boolean {
        return this.hasONCols(ε);
    }

    isSymmetric(ε: number=1e-10): boolean {
        return this.equals(this.transpose, ε);
    }


    /***********************************************************************/
    /* Solving Ax=b                                                        */
    /***********************************************************************/

    solveForwardSubst(b: Vector): Vector {
        let x = b.clone();
        for(let i = 1; i < this.m; ++i) {
            for(let ii = 0; ii < i; ++ii)
                x[i] -= this[i][ii] * x[ii];
        }
        return x;
    }

    solveReverseSubst1(b: Vector): Vector {
        let x = b.clone();
        for(let p = this.n - 1; p >= 0; --p) {
            for(let pp = p + 1; pp < this.n; ++pp) {
                x[p] -= this[p][pp] * x[pp];
            }
            x[p] /= this[p][p];
        }
        return x;
    }

    solve(b: Vector, hints?: {
        PLU?: PLUData,
        LDLT?: LDLTData,
        ε?: number}): Vector | null {

        if(hints?.LDLT) {
            let { L, D } = hints.LDLT;
            let y = L.solveForwardSubst(b);
            for(let i = 0; i < this.n; ++i)
                y[i] /= D[i];
            return L.transpose.solveReverseSubst1(y);
        }

        return super.solve(b, hints);
    }


    inverse(JD: JordanData=this.jordanSubst()): SquareMatrix | null {
        if(this.rank(JD) < this.n)
            return null;
        return JD.rowOps;
    }

    isInvertible(JD: JordanData=this.jordanSubst()): boolean {
        return !!this.inverse(JD);
    }

    isSingular(JD: JordanData=this.jordanSubst()): boolean {
        return !this.isInvertible(JD);
    }


    /***********************************************************************/
    /* Determinants and Eigenvalues                                        */
    /***********************************************************************/

    det(PLU: PLUData=this.PLU()): number {
        return PLU.det!;
    }

    fadeevLeverrier(): {charpoly: Polynomial, adjugate: SquareMatrix} {
        const n = this.n;
        let ret = [1];
        let AM = Matrix.zero(n) as SquareMatrix;
        let adjugate: SquareMatrix;
        let c = 1;
        for(let k = 1; k <= n; ++k) {
            for(let i = 0; i < n; ++i)
                AM[i][i] += c;
            if(k == n) adjugate = AM;
            AM = this.mult(AM) as SquareMatrix;
            c = -AM.trace/k;
            ret.push(c);
        }
        // This is det(λI - A); multiply by (-1)^n now
        if(n % 2 === 1) {
            for(let i = 0; i < ret.length; ++i)
                ret[i] *= -1;
        }
        if(n % 2 === 0) adjugate!.scale(-1);
        return {charpoly: new Polynomial(...ret), adjugate: adjugate!};
    }

    get charpoly(): Polynomial {
        return this.fadeevLeverrier().charpoly;
    }

    get adjugate(): SquareMatrix {
        return this.fadeevLeverrier().adjugate;
    }

    minor(i: number, j: number): SquareMatrix {
        let rows = [...this.rows];
        rows.splice(i, 1)
        let ret = Array.from(
            rows,
            row => {
                let entries = [...row];
                entries.splice(j, 1);
                return Vector.from(entries);
            });
        return Matrix.from(ret) as SquareMatrix;
    }

    cofactor(i: number, j: number, ε: number=1e-10): number {
        let A = this.minor(i, j);
        return A.det(A.PLU(ε)) * ((i + j) % 2 === 0 ? 1 : -1);
    }

    eigenvalues(ε: number=1e-10): MultRoot[] {
        return this.charpoly.factor(ε);
    }

    eigenspaceBasis(λ: number | Complex,
                    ε: number=1e-10): Vector[] | [Vector, Vector][] {
        if(λ instanceof Complex) {
            if(Math.abs(λ.Im) > ε)
                return this._complexEigenspaceBasis(λ, ε);
            λ = λ.Re;
        }
        return this._realEigenspaceBasis(λ, ε);
    }

    _realEigenspaceBasis(λ: number, ε: number=1e-10): Vector[] {
        let AmlI = this.clone().sub(Matrix.identity(this.n, λ));
        let B = AmlI.nullBasis(AmlI.jordanSubst(AmlI.PLU(ε)));
        if(B.length == 0)
            throw new MatrixError("λ is not an eigenvalue of this matrix");
        return B;
    }

    _complexEigenspaceBasis(λ: Complex, ε=1e-10): [Vector, Vector][] {
        // The row reduction algorithm in PLU() won't work for complex
        // matrices.  We implement a simplified version here.
        let { m, n } = this;
        let pivots: Pivot[] = [];
        let U = Array.from(
            this.rows, row => Array.from(row, x => new Complex(x)));
        for(let i = 0; i < n; ++i) U[i][i].sub(λ);

        for(let curRow = 0, curCol = 0; curRow < m && curCol < n; ++curCol) {
            // Find maximal pivot
            let pivot = U[curRow][curCol], row = curRow;
            for(let i = curRow+1; i < m; ++i) {
                if(U[i][curCol].modsq > pivot.modsq) {
                    pivot = U[i][curCol];
                    row = i;
                }
            }
            if(pivot.modsq > ε*ε) {
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
            throw new MatrixError("λ is not an eigenvalue of this matrix");

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

        // Now extract the null basis
        let basis: [Vector, Vector][] = [], previous: Pivot[] = [];
        for(let j = 0; j < n; ++j) {
            if(pivots.length && pivots[0][1] === j) {
                // Pivot column
                previous.push(pivots.shift()!);
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

        return basis;
    }

    eigenspace(λ: number, ε: number=1e-10) {
        return Subspace.spanOf(this._realEigenspaceBasis(λ, ε));
    }

    diagonalize(hints?: {block?: boolean, ortho?:
                         boolean, ε?: number}): DiagonalizationData | null {
        let { block, ortho, ε }
            = {block: false, ortho: false, ε:1e-10, ...hints};
        let eigenbasis: Vector[] = [];
        let D = Matrix.zero(this.n) as SquareMatrix;
        let i = 0;
        // Only use one of a conjugate pair of eigenvalues
        for(let [λ,m] of this.eigenvalues(ε).filter(
            ([λ,]) => !(λ instanceof Complex) || λ.Im >= 0)) {
            if(λ instanceof Complex) {
                if(!block) return null;
                let B = this.eigenspaceBasis(λ, ε) as [Vector, Vector][];
                if(B.length < m) // Impossible for matrices <= 3x3
                    return null;
                for(let j = 0; j < m; ++j, i += 2) {
                    // The columns are the real and complex parts of the eigenvectors
                    eigenbasis.push(...B[j]);
                    D.insertSubmatrix(i, i, new Matrix([λ.Re, λ.Im], [-λ.Im, λ.Re]));
                }
            } else {
                let B = this.eigenspaceBasis(λ, ε) as Vector[];
                if(B.length < m)
                    return null;
                if(ortho) {
                    if(B.length == 1)
                        B[0].normalize();
                    else {
                        let M = Matrix.from(B).transpose;
                        B = [...M.QR().Q.cols];
                    }
                }
                for(let j = 0; j < m; ++j, ++i) {
                    D[i][i] = λ;
                    eigenbasis[i] = B[j];
                }
            }
        }
        let C = Matrix.from(eigenbasis).transpose as SquareMatrix;
        return { C, D };
    }

    isDiagonalizable(hints?: {block?: boolean, ortho?:
                              boolean, ε?: number}): boolean {
        return !!this.diagonalize(hints);
    }

    LDLT(ε: number=1e-10): LDLTData | null {
        let n = this.n;
        let D = new Array(n);
        let L = Matrix.identity(n);
        for(let j = 0; j < n; ++j) {
            D[j] = this[j][j];
            for(let k = 0; k < j; ++k)
                D[j] -= L[j][k]*L[j][k]*D[k];
            if(Math.abs(D[j]) <= ε)
                return null;
            for(let i = j+1; i < n; ++i) {
                let l = this[i][j];
                for(let k = 0; k < j; ++k)
                    l -= L[i][k]*L[j][k]*D[k];
                L[i][j] = l / D[j];
            }
        }
        return { L, D };
    }

    cholesky(LDLT: LDLTData | null=this.LDLT()): SquareMatrix | null {
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

    isPosDef(LDLT: LDLTData | null=this.LDLT()): boolean {
        if(!LDLT) return false;
        return LDLT.D.every(x => x > 0);
    }
}

export default Matrix;
export { SquareMatrix};
