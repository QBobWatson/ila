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

/** @module subspace
 *
 * @file
 * Implements a Subspace class to represent subspaces of `R^n`.
 */

// TODO: principal angles?

import Matrix from "./matrix3";
import { SquareMatrix } from "./matrix3";
import Vector from "./vector";
import { range } from "./util";


class Subspace {
    readonly n: number;
    readonly Q: Matrix | null;
    readonly QT: Matrix | null;
    protected _projectionMatrix: SquareMatrix | undefined;

    static spanOf(generators: Matrix | (Vector | number[])[], hints?:
                  {n?: number, isON?: boolean, ε?: number}): Subspace {
        let { n, isON, ε } = hints || {};
        if(n === undefined) {
            if(generators instanceof Matrix)
                n = generators.m;
            else if(generators.length > 0) {
                if(generators[0] instanceof Vector)
                    n = generators[0].size;
                else
                    n = generators[0].length;
            }
            else
                throw new Error("Cannot determine the ambient R^n from zero generators");
        }
        if(generators instanceof Matrix)
            return new Subspace(generators, n, isON);
        if(generators.length > 0)
            return new Subspace(
                Matrix.create(...generators).transpose, n, isON, ε);
        return new Subspace(null, n);
    }

    static nullSpaceOf(A: Matrix, ε=1e-10): Subspace {
        let generators = A.nullBasis(A.jordanSubst(A.PLU(ε)));
        if(generators.length == 0)
            return new Subspace(null, A.n);
        let M = Matrix.from(generators).transpose;
        return new Subspace(M, A.n, undefined, ε);
    }

    static Rn(n: number): Subspace {
        return new Subspace(Matrix.identity(n), n, true);
    }

    static zero(n: number): Subspace {
        return new Subspace(null, n);
    }

    protected constructor(A: Matrix | null, n: number,
                          isON?: boolean, ε=1e-10) {
        this.n = n;
        if(A === null)
            this.Q = this.QT = null;
        else if(isON)  {
            this.Q = A;
            this.QT = A.transpose;
        } else {
            let { Q, LD } = A.QR(ε);
            let rows = Array.from(Q.cols).filter((_, i) => !LD.includes(i));
            if(rows.length == 0)
                this.Q = this.QT = null;
            else {
                this.QT = Matrix.create(...rows);
                this.Q = this.QT.transpose;
            }
        }
    }

    clone(): Subspace {
        let Q = this.isZero() ? null : this.Q!.clone();
        return new Subspace(Q, this.n, true);
    }

    get dim(): number {
        return this.Q === null ? 0 : this.Q.n;
    }

    get basis(): Matrix | null {
        return this.Q;
    }

    get projectionMatrix(): SquareMatrix {
        if(this._projectionMatrix)
            return this._projectionMatrix;
        if(this.isZero())
            this._projectionMatrix = Matrix.zero(this.n) as SquareMatrix;
        else if(this.isMaximal())
            this._projectionMatrix = Matrix.identity(this.n) as SquareMatrix;
        else
            this._projectionMatrix = this.Q!.mult(this.QT!) as SquareMatrix;
        return this._projectionMatrix;
    }

    toString(precision=4): string {
        if(this.dim === this.n)
            return `The full subspace R^${this.n}`;
        if(this.dim === 0)
            return `The zero subspace of R^${this.n}`;
        let ret = `Subspace of R^${this.n} of dimension ${this.dim} with basis\n`;
        let rowStrings = Array.from(
            this.Q!.rows, row => Array.from(row, x => x.toFixed(precision)));
        let colLengths = Array.from(
            range(this.dim), i => Math.max(...rowStrings.map(row => row[i].length)));
        return ret + rowStrings.map(
            row => row.map((x, i) => `[${x.padStart(colLengths[i])}]`).join(' '))
            .join('\n');
     }

    isMaximal(): boolean {
        return this.dim === this.n;
    }

    isZero(): boolean {
        return this.dim === 0;
    }

    add(other: Subspace, ε=1e-10): Subspace {
        if(this.n != other.n)
            throw new Error("Tried to add subspaces of different R^n");
        if(this.isZero())
            return other.clone();
        if(other.isZero())
            return this.clone();
        if(this.isMaximal() || other.isMaximal())
            return Subspace.Rn(this.n);
        return Subspace.spanOf([...this.QT!.rows, ...other.QT!.rows],
                               {n: this.n, ε});
    }

    intersect(other: Subspace, ε=1e-10): Subspace {
        if(this.n != other.n)
            throw new Error("Tried to add subspaces of different R^n");
        if(this.isZero() || other.isZero())
            return Subspace.zero(this.n);
        if(this.isMaximal())
            return other.clone();
        if(other.isMaximal())
            return this.clone();
        return this.perp(ε).add(other.perp(ε)).perp(ε);
    }

    project(v: Vector): Vector {
        return this.projectionMatrix.apply(v);
    }

    complement(v: Vector): Vector {
        return this.orthoDecomp(v)[1];
    }

    orthoDecomp(v: Vector): [Vector, Vector] {
        let v1 = this.project(v);
        return [v1, v.clone().sub(v1)];
    }

    distanceTo(v: Vector): number {
        return this.complement(v).len;
    }

    contains(v: Vector | number[], ε=1e-10): boolean {
        if(!(v instanceof Vector))
            v = new Vector(...v);
        if(this.n != v.size)
            throw new Error("Vector has the wrong number of entries");
        if(this.isMaximal())
            return true; // Nothing to check!
        if(this.isZero())
            return v.lensq <= ε*ε;
        return v.lensq <= ε*ε || this.distanceTo(v.clone().normalize()) <= ε;
    }

    isOrthogonalTo(v: Vector, ε=1e-10): boolean {
        if(this.n != v.size)
            throw new Error("Vector has the wrong number of entries");
        if(this.isMaximal())
            return v.isZero(ε); // Nothing to check!
        if(this.isZero())
            return true;
        return v.lensq <= ε*ε
            || this.project(v.clone().normalize()).lensq <= ε*ε;
    }

    isSubspaceOf(other: Subspace, ε=1e-10): boolean {
        if(this.n != other.n)
            throw new Error("Tried to test containment of subspaces in different R^n");
        if(this.dim > other.dim)
            return false;
        if(this.isZero() || other.isMaximal())
            return true;
        if(other.isZero())
            return false;
        for(const col of this.Q!.cols) {
            if(!other.contains(col, ε))
                return false;
        }
        return true;
    }

    equals(other: Subspace, ε=1e-10): boolean {
        if(this.n != other.n)
            throw new Error("Tried to test equality of subspaces in different R^n");
        if(this.dim !== other.dim)
            return false;
        return this.isSubspaceOf(other, ε);
    }

    perp(ε=1e-10) {
        if(this.isZero())
            return Subspace.Rn(this.n);
        if(this.isMaximal())
            return Subspace.zero(this.n);
        return Subspace.nullSpaceOf(this.QT!, ε);
    }
}


export default Subspace;
