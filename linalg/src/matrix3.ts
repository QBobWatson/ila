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

// TODO: MatrixError, and the like in the other modules

import Vector from "./vector";

import { range } from "./util";


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

    static create(...rows: (number[] | Vector)[]): Matrix | SquareMatrix {
        const row = rows[0];
        if(rows.length === (row ?
            (row instanceof Vector ? row.size : row.length) : 0))
            return new SquareMatrix(...rows);
        return new Matrix(...rows);
    }

    protected constructor(...rows: (number[] | Vector)[]) {
        if(rows.length === 0)
            throw new Error("A matrix must have at least one row.")
        let rows2 = rows.map(r => r instanceof Vector
            ? r : new Vector(...r)) as Vector[];
        this.m = rows2.length;
        this.n = rows2[0]!.size;
        if(rows2.some(v => v.size !== this.n))
            throw new Error("Matrix rows must have the same length.");
        for(let i = 0; i < this.m; ++i)
            this[i] = rows2[i];
    }

    static identity(n: number, λ: number=1): SquareMatrix {
        return new SquareMatrix(
            ...Array.from(range(n), i => Vector.e(i, n, λ)));
    }

    static zero(m: number, n: number=m): Matrix | SquareMatrix {
        return Matrix.create(
            ...Array.from(range(m), _ => Vector.zero(n)));
    }

    static constant(c: number, m: number, n: number=m): Matrix | SquareMatrix {
        return Matrix.create(
            ...Array.from(range(m), _ => Vector.constant(n, c)));
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
        return new SquareMatrix(
            ...Array.from(range(n), i => Vector.e(vals[i], n)));
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
        return new Vector(...Array.from(this.rows, row => row[j]));
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
            throw new Error('Tried to add matrices of different sizes');
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
            throw new Error('Cannot multiply matrices of incompatible dimensions');
        return Matrix.create(
            ...Array.from(
                this.rows,
                row => new Vector(...Array.from(
                    range(other.n), i => [...row].reduce(
                        (a, v, j) => a + v * other[j][i], 0)))));
    }

    apply(v: Vector): Vector {
        if(v.size !== this.n)
            throw new Error('Cannot multiply matrix and vector of incompatible dimensions');
        return new Vector(...Array.from(this.rows, row => row.dot(v)));
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

    colBasis(PLU: PLUData=this.PLU()): Vector[] {
        return Array.from(PLU.pivots, ([, j]) => this.col(j));
    }

    rowBasis(PLU: PLUData=this.PLU()): Vector[] {
        const { pivots, U } = PLU;
        return pivots.map(([i,]) => U[i].clone());
    }

    leftNullBasis(PLU: PLUData=this.PLU()): Vector[] {
        const { E, pivots } = PLU;
        const r = pivots.length;
        let ret: Vector[] = new Array(this.m - r);
        for(let i = r; i < this.m; ++i)
            ret[i - r] = E[i].clone();
        return ret;
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
            throw new Error("Incompatible dimensions of matrix and vector");

        hints = {ε: 1e-10, ...hints};
        const { ε, PLU } = hints;
        const { P, L, U, pivots } = PLU || this.PLU(ε);
        // Solve LUx = PAx = Pb
        const Pb = new Vector(
            ...Array.from(range(this.m), i => b[P[i]]));
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


    /***********************************************************************/
    /* Convenience                                                         */
    /***********************************************************************/

    pivots(PLU: PLUData=this.PLU()): Pivot[] {
        return PLU.pivots;
    }

    rank(PLU: PLUData=this.PLU()): number {
        return PLU.pivots.length;
    }

    nullity(PLU: PLUData=this.PLU()): number {
        return this.n - PLU.pivots.length;
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
}


class SquareMatrix extends Matrix {

    constructor(...rows: (number[] | Vector)[]) {
        super(...rows);
        if(this.m !== this.n)
            throw new Error("Tried to create a non-square SquareMatrix.");
    }

    isOrthogonal(ε: number=1e-10): boolean {
        return this.hasONCols(ε);
    }

    isSymmetric(ε: number=1e-10): boolean {
        return this.equals(this.transpose, ε);
    }

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
}

export default Matrix;
export { SquareMatrix};
