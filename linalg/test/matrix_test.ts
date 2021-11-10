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

import { should, Assertion, expect } from 'chai';
should();

import './lib/resemble';

import Matrix from "../src/matrix3";
import { SquareMatrix, PLUData, QRData, JordanData } from "../src/matrix3";
import Vector from "../src/vector";
import Polynomial from '../src/polynomial';
import Complex from '../src/complex';

const mat = Matrix.create;
const smat = mat as (...rows: (number[] | Vector)[]) => SquareMatrix;
const vec = (...entries: number[]) => new Vector(...entries);
const poly = (...coeffs: number[]) => new Polynomial(...coeffs);
const C = (a: number | Complex, b=0) => new Complex(a, b);


type Factorization = 'PA=LU' | 'QR' | 'SVD' | 'LDLT';

Assertion.addMethod(
    'factorize', function(this: typeof Assertion, A: Matrix,
                          which: Factorization, ε: number=1e-10) {
        switch (which) {
            case 'PA=LU':
                let { P, L, U, E, pivots } = this._obj as PLUData;
                expect(L.isLowerUni(),
                       `L = \n${L.toString(2)}`
                       + `\nis not lower unipotent`).to.be.true;
                expect(U.isREF(),
                       `U =\n${U.toString(2)}\nis not in REF`).to.be.true;
                expect(U.leadingEntries()).to.eql(pivots);
                let PA = Matrix.permutation(P).mult(A);
                let LU = L.mult(U);
                expect(PA.equals(LU, ε),
                       `Matrix is not correctly factored:`
                       + ` PA =\n${PA.toString(2)}\nLU =\n${LU.toString(2)}`)
                        .to.be.true;
                let EA = E.mult(A);
                expect(EA.equals(U, ε),
                       `Matrix is not correctly factored:`
                       + ` EA =\n${EA.toString(2)}\nU =\n${U.toString(2)}`)
                        .to.be.true;
                break;

            case 'QR':
                let { Q, R, LD } = this._obj as QRData;
                expect(R.isUpperTri(),
                       `R is not upper-triangular: R=\n${R.toString(2)}`)
                        .to.be.true;
                let QTQ = Q.transpose.mult(Q);
                // Might not have orthonormal cols since some are zero
                expect(QTQ.isDiagonal(ε),
                       `Q does not have orthogonal columns:`
                       + ` Q=\n${Q.toString(2)}\nQTQ=\n${QTQ.toString(2)}`)
                        .to.be.true;
                for(let j = 0; j < QTQ.n; ++j) {
                    if(LD.includes(j)) {
                        expect(Q.col(j).isZero(ε),
                               `Column ${j} of Q should be zero:`
                               + ` LD=${LD}, Q=\n${Q.toString(2)}`)
                            .to.be.true;
                        expect(R.row(j).isZero(ε),
                               `Row ${j} of R should be zero:`
                               + ` LD=${LD}, R=\n${R.toString(2)}`)
                            .to.be.true;
                    } else {
                        expect(QTQ[j][j]).to.be.approximately(1, ε);
                        expect(R[j][j]).to.be.above(ε);
                    }
                }
                let QR = Q.mult(R);
                expect(QR.equals(A, ε),
                       `Matrix is not correctly factored:`
                       + ` QR =\n${QR.toString(2)}\nA =\n${A.toString(2)}`)
                        .to.be.true;
                break;
        }
    }
);

Assertion.addMethod(
    'solve', function(this: typeof Assertion, M: Matrix, b: Vector,
                      leastSquares: boolean=false, ε: number=1e-10) {
        let x = this._obj;
        expect(x, 'A solution was not found when one was expected')
            .not.to.be.null;
        let Mx = M.apply(x);
        if(!leastSquares)
            expect(Mx.equals(b, ε),
                   `Mx != b:\nx = ${x.toString(2)}\n`
                   + `Mx = ${Mx.toString(2)}\nb = ${b.toString(2)}`)
                .to.be.true;
        else {
            expect(M.transpose.apply(Mx.sub(b)).isZero(1e-10),
                   'Expected Mx-b to be orthogonal to Col(M)')
                .to.be.true;
        }
    }
);

Assertion.addMethod(
    'eigenvectors', function(this: typeof Assertion, M: SquareMatrix,
                             λ: number | Complex, m: number, ε: number=1e-10) {
        this._obj.length.should.equal(m);
        if(typeof λ === "number") {
            for(let v of this._obj) {
                let Mv = M.apply(v);
                let λv = v.clone().scale(λ);
                Mv.should.resemble(λv, ε);
            }
        } else {
            for(let [v_Re, v_Im] of this._obj) {
                let Mv_Re = M.apply(v_Re);
                let Mv_Im = M.apply(v_Im);
                let λv_Re = v_Re.clone().scale(λ.Re).sub(v_Im.clone().scale(λ.Im));
                let λv_Im = v_Re.clone().scale(λ.Im).add(v_Im.clone().scale(λ.Re));
                [Mv_Re, Mv_Im].should.resemble([λv_Re, λv_Im], ε);
            }
        }
    }
);

Assertion.addMethod(
    'diagonalize', function(this: typeof Assertion, A: SquareMatrix,
                            block=false, ε=1e-10) {
        expect(this._obj, `No diagonalization was found when one was expected`)
            .not.to.be.null;
        let { C, D } = this._obj as { C: SquareMatrix, D: SquareMatrix };
        if(!block) expect(D.isDiagonal()).to.be.true;
        let Cinv = C.inverse(C.jordanSubst(C.PLU(ε)));
        expect(Cinv).not.to.be.null;
        C.mult(D).equals(A.mult(C), ε).should.be.true;
    }
);

declare global {
    export namespace Chai {
        interface Assertion {
            factorize(A: Matrix, which: Factorization, ε?: number): void;
            solve(M: Matrix, b: Vector, leastSquares?: boolean, ε?: number): void;
            eigenvectors(M: SquareMatrix, λ: number | Complex,
                         m: number, ε?: number): void;
            diagonalize(A: SquareMatrix, block?: boolean, ε?: number): void;
        }
    }
}


describe('Matrix', () => {

    /***********************************************************************/
    /* Constructors                                                        */
    /***********************************************************************/

    describe('#create()', () => {
        it('should create a 3x2 matrix', () => {
            let A = mat([1, 2], [3, 4], [5, 6]);
            A.should.include({ m: 3, n: 2 });
            A.should.not.be.an.instanceOf(SquareMatrix);
        });
        it('should create a square matrix', () =>
            mat([1, 2, 3], [4, 5, 6], [7, 8, 9]).should
                .be.an.instanceOf(SquareMatrix));
    });
    describe('#identity()', () => {
        it('should return a 3x3 identity matrix', () => {
            let M = Matrix.identity(3);
            M.equals(mat([1, 0, 0], [0, 1, 0], [0, 0, 1])).should.be.true;
            M.should.be.an.instanceof(SquareMatrix);
        });
        it('should return a scaled 3x3 identity matrix', () =>
            Matrix.identity(3, 2).equals(
                mat([2, 0, 0], [0, 2, 0], [0, 0, 2])).should.be.true);
    });
    describe('#zero()', () => {
        it('should return the 2x3 zero matrix', () =>
            Matrix.zero(2, 3).equals(
                mat([0, 0, 0], [0, 0, 0])).should.be.true);
        it('should return the 2x2 zero matrix', () => {
            let M = Matrix.zero(2);
            M.equals(mat([0, 0], [0, 0])).should.be.true;
            M.should.be.an.instanceof(SquareMatrix);
        });
    });
    describe('#diagonal()', () => {
        it('should create a square diagonal matrix', () =>
            Matrix.diagonal([1, 2, 3]).equals(
                mat([1, 0, 0], [0, 2, 0], [0, 0, 3])).should.be.true);
        it('should create a tall diagonal matrix', () =>
            Matrix.diagonal([1, 2, 3], 4).equals(
                mat([1, 0, 0], [0, 2, 0], [0, 0, 3], [0, 0, 0])).should.be.true);
        it('should create a wide diagonal matrix', () =>
            Matrix.diagonal([1, 2, 3], 3, 4).equals(
                mat([1, 0, 0, 0], [0, 2, 0, 0], [0, 0, 3, 0])).should.be.true);
        it('should create a big diagonal matrix', () =>
            Matrix.diagonal([1, 2, 3], 4, 5).equals(
                mat([1, 0, 0, 0, 0], [0, 2, 0, 0, 0], [0, 0, 3, 0, 0], [0, 0, 0, 0, 0]))
                .should.be.true);
        it('should create a small diagonal matrix', () =>
            Matrix.diagonal([1, 2, 3, 4], 2, 3).equals(
                mat([1, 0, 0], [0, 2, 0])).should.be.true);
    });
    describe('#permutation()', () =>
        it('should create a 3x3 permutation matrix', () => {
            let A = Matrix.permutation([2, 0, 1]);
            A.equals(mat([0, 0, 1], [1, 0, 0], [0, 1, 0])).should.be.true;
            A.should.be.an.instanceof(SquareMatrix);
        }));
    describe('#clone()', () => {
        let M = mat([1, 2], [3, 4]), N = M.clone();
        it('should create distinct objects', () =>
            M.should.not.equal(N));
        it('should create identical objects', () =>
            M.equals(N).should.be.true);
        it('should create a SquareMatrix', () =>
            N.should.be.an.instanceOf(SquareMatrix));
        let M2 = mat([1, 2, 3], [4, 5, 6]), N2 = M2.clone();
        it('should create a non-square Matrix', () =>
            N2.should.not.be.an.instanceOf(SquareMatrix));
    });

    /***********************************************************************/
    /* Accessors                                                           */
    /***********************************************************************/

    describe('#transpose', () => {
        it('should construct the transpose', () => {
            let A = mat([1, 2, 3], [4, 5, 6]), B = A.transpose;
            B.equals(mat([1, 4], [2, 5], [3, 6])).should.be.true;
            B.should.not.be.an.instanceOf(SquareMatrix);
        });
    });
    describe('#normal', () =>
        it('should construct the normal matrix', () => {
            let A = mat([1, 2], [3, 4], [5, 6]), B = A.normal;
            B.equals(mat([35, 44], [44, 56])).should.be.true;
            B.should.be.an.instanceOf(SquareMatrix);
        }));
    describe('#trace', () => {
        it('should compute the sum of the diagonal entries', () =>
            mat([1, 2, 3], [4, 5, 6], [7, 8, 9]).trace.should.equal(1 + 5 + 9));
        it('should work for non-square matrices', () =>
            mat([1, 2, 3], [4, 5, 6]).trace.should.equal(1 + 5));
        it('should work for non-square matrices', () =>
            mat([1, 2], [3, 4], [5, 6]).trace.should.equal(1 + 4));
    });
    describe('#row() and #col()', () => {
        let M = mat([1, 2], [3, 4]);
        it('should return the second row', () =>
            M.row(1).should.eql(vec(3, 4)));
        it('should return the second column', () =>
            M.col(1).should.eql(vec(2, 4)));
    });
    describe('#rows() and #cols()', () => {
        let M = mat([1, 2], [3, 4]);
        it('should iterate over the rows', () =>
            Array.from(M.rows).should.eql([vec(1, 2), vec(3, 4)]));
        it('should iterate over the columns', () =>
            Array.from(M.cols).should.eql([vec(1, 3), vec(2, 4)]));
    });
    describe('#diag()', () => {
        it('should iterate over diagonal entries of square matrices', () =>
            Array.from(mat([1, 2, 3], [4, 5, 6], [7, 8, 9]).diag)
                .should.eql([1, 5, 9]));
        it('should iterate over diagonal entries of tall matrices', () =>
            Array.from(mat([1, 2], [4, 5], [7, 8]).diag)
                .should.eql([1, 5]));
        it('should iterate over diagonal entries of wide matrices', () =>
            Array.from(mat([1, 2, 3], [4, 5, 6]).diag)
                .should.eql([1, 5]));
    });

    /***********************************************************************/
    /* Utility                                                             */
    /***********************************************************************/

    describe('#toString()', () => {
        let M = mat([10, 2], [3, 4]);
        it('should have 4 decimal places by default', () =>
            M.toString().should.eql("[10.0000 2.0000]\n[ 3.0000 4.0000]"));
        it('can have other precision', () =>
            M.toString(2).should.eql("[10.00 2.00]\n[ 3.00 4.00]"));
    });
    describe('#insertSubmatrix()', () => {
        it('should insert a submatrix correctly', () => {
            Matrix.zero(4, 5).insertSubmatrix(1, 2, mat([1, 1], [1, 1]))
                .should.eql(mat([0, 0, 0, 0, 0],
                                [0, 0, 1, 1, 0],
                                [0, 0, 1, 1, 0],
                                [0, 0, 0, 0, 0]));
        });
    });
    describe('#equals()', () => {
        let M = mat([0, 1, 2], [3, 4, 5]), M1 = M.clone();
        let N = mat([0.01, 1.01, 2.01], [3, 4, 5]);
        let M3 = mat([0, 1, 2]);
        it('should compare two matrices as equal', () =>
            M.equals(M1).should.be.true);
        it('should compare different matrices as not equal', () =>
            M.equals(N).should.be.false);
        it('should compare similar matrices as equal with ε>0', () =>
            M.equals(N, 0.02).should.be.true);
        it('should compare matrices of different sizes as not equal', () =>
            M3.equals(M).should.be.false);
    });
    describe('#leadingEntries()', () => {
        it('should compute leading entries', () =>
            mat([1,0],[0,1]).leadingEntries().should.eql([[0,0], [1,1]]));
        it('should compute leading entries', () =>
            mat([0,0],[0,1]).leadingEntries().should.eql([[1,1]]));
        it('should compute leading entries', () =>
            mat([0,0,0],[0,1,0],[0,0,0],[0,0,1])
                .leadingEntries().should.eql([[1,1],[3,2]]));
    });

    /***********************************************************************/
    /* Traits                                                              */
    /***********************************************************************/

    describe('#isZero()', () => {
        it('should detect the zero matrix', () => {
            let M = Matrix.identity(3, 0.01);
            M.isZero().should.be.false;
            M.isZero(0.1).should.be.true;
        });
    });
    describe('#isUpperTri(), #isLowerTri(), #isDiagonal()', () => {
        let M = mat([1, 1, 1], [0, 1, 1], [0, 0, 1]);
        let N = mat([0, 0, 0], [0, 0, 0], [0, 1, 0]);
        it('should be upper-triangular', () =>
            M.isUpperTri().should.be.true);
        it('should not be upper-triangular', () =>
            N.isUpperTri().should.be.false);
        it('should be lower-triangular', () =>
            M.transpose.isLowerTri().should.be.true);
        it('should not be lower-triangular', () =>
            N.transpose.isLowerTri().should.be.false);
        it('should not be diagonal', () => {
            M.isDiagonal().should.be.false;
            N.isDiagonal().should.be.false;
            mat([1,0,0,0],[0,2,0,0],[0,0,3,0]).isDiagonal().should.be.true;
        });
    });
    describe('#isUpperUni() and #isLowerUni()', () => {
        let M = mat([1, 1, 1], [0, 1, 1], [0, 0, 1]);
        let N = mat([1, 1, 1], [0, 2, 1], [0, 0, 1]);
        it('should be upper-uniponent', () =>
            M.isUpperUni().should.be.true);
        it('should not be upper-unipotent', () =>
            N.isUpperUni().should.be.false);
        it('should be lower-unipotent', () =>
            M.transpose.isLowerUni().should.be.true);
        it('should not be lower-unipotent', () =>
            N.transpose.isLowerUni().should.be.false);
    });
    describe('#isREF() and #isRREF()', () => {
        let testMats1 = [
            mat([1,  0,  2],
                [0,  1, -1]),
            mat([0,  1,  8, 0]),
            mat([1, 17,  0],
                [0,  0,  1]),
            Matrix.zero(2, 3)
        ];
        it('should be in REF and RREF', () => {
            for(let M of testMats1) {
                M.isREF().should.be.true;
                M.isRREF().should.be.true;
            }
        });
        let testMats2 = [
            mat([2,  1],
                [0,  1]),
            mat([2,  7, 1, 4],
                [0,  0, 2, 1],
                [0,  0, 0, 3]),
            mat([1, 17, 0],
                [0,  1, 1]),
            mat([2,  1, 3],
                [0,  0, 0])
        ];
        it('should be in REF but not RREF', () => {
            for(let M of testMats2) {
                M.isREF().should.be.true;
                M.isRREF().should.be.false;
            }
        });
        let testMats3 = [
            mat([2,  7, 1, 4],
                [0,  0, 2, 1],
                [0,  0, 1, 3]),
            mat([0, 17, 0],
                [0,  2, 1]),
            mat([2,  1],
                [2,  1]),
            mat([0,1,0,0]).transpose
        ];
        it('should not be in REF or RREF', () => {
            for(let M of testMats3) {
                M.isREF().should.be.false;
                M.isRREF().should.be.false;
            }
        });
    });
    describe('#hasONCols()', () => {
        it('detects matrices with orthonormal columns', () =>
            mat([1,1],[1,-1],[1,0])
                .mult(Matrix.diagonal([1/Math.sqrt(3), 1/Math.sqrt(2)]))
                .hasONCols().should.be.true);
        it('detects matrices without orthonormal columns', () =>
            mat([1,1],[1,-1],[1,1]).scale(1/Math.sqrt(3)).hasONCols()
                .should.be.false);
    });

    /***********************************************************************/
    /* Arithmetic                                                          */
    /***********************************************************************/

    describe('#add() and #sub()', () => {
        it('should add componentwise', () =>
            mat([1, 2], [3, 4]).add(mat([2, 1], [4, 3])).equals(
                mat([3, 3], [7, 7])).should.be.true);
        it('should add with a factor', () =>
            mat([1, 2], [3, 4]).add(mat([2, 1], [4, 3]), 2).equals(
                mat([5, 4], [11, 10])).should.be.true);
        it('should subtract componentwise', () =>
            mat([1, 2], [3, 4]).sub(mat([2, 1], [4, 3])).equals(
                mat([-1, 1], [-1, 1])).should.be.true);
        it('should throw when matrices have different sizes', () => {
            let M = mat([1, 2, 3], [4, 5, 6]);
            M.add.bind(M, mat([1, 2], [3, 4])).should.throw(/different sizes/);
        });
    });
    describe('#scale()', () => {
        it('should scale entries', () =>
            mat([1, 2], [3, 4]).scale(2).equals(
                mat([2, 4], [6, 8])).should.be.true);
        it('should add with a factor', () =>
            mat([1, 2], [3, 4]).add(mat([2, 1], [4, 3]), 2).equals(
                mat([5, 4], [11, 10])).should.be.true);
        it('should subtract componentwise', () =>
            mat([1, 2], [3, 4]).sub(mat([2, 1], [4, 3])).equals(
                mat([-1, 1], [-1, 1])).should.be.true);
        it('should throw when matrices have different sizes', () => {
            let M = mat([1, 2, 3], [4, 5, 6]);
            M.add.bind(M, mat([1, 2], [3, 4])).should.throw(/different sizes/);
        });
    });
    describe('#mult()', () => {
        it('should compute the product', () => {
            let A = mat([1, 2], [3, 4], [5, 6]).mult(mat([1, 2, 3], [4, 5, 6]));
            A.equals(mat([9, 12, 15], [19, 26, 33], [29, 40, 51]))
                .should.be.true;
            A.should.be.an.instanceOf(SquareMatrix);
        });
        it('should throw when the matrices have incompatible dimensions', () => {
            let M = mat([1, 2], [3, 4], [5, 6]);
            M.mult.bind(M, M).should.throw(/incompatible dimensions/);
        });
    });
    describe('#apply()', () => {
        it('should compute the product', () =>
            mat([1, 2, 3], [4, 5, 6]).apply(vec(7, 8, 9))
                .equals(vec(50, 122)).should.be.true);
        it('should throw when the vector has incompatible length', () => {
            let M = mat([1, 2], [3, 4], [5, 6]);
            M.apply.bind(M, vec(1, 2, 3))
                .should.throw(/incompatible dimensions/);
        });
    });

    /***********************************************************************/
    /* Row Operations                                                      */
    /***********************************************************************/

    describe('#rowScale(), #rowReplace(), #rowSwap()', () => {
        let M = mat([1, 2, 3], [4, 5, 6], [7, 8, 9]);
        it('should scale one row', () =>
            M.rowScale(1, 2).equals(mat(
                [1, 2, 3], [8, 10, 12], [7, 8, 9])).should.be.true);
        it('should add 2 x one row to another', () =>
            M.rowReplace(0, 2, 2).equals(mat(
                [15, 18, 21], [8, 10, 12], [7, 8, 9])).should.be.true);
        it('should swap two rows', () =>
            M.rowSwap(1, 2).equals(mat(
                [15, 18, 21], [7, 8, 9], [8, 10, 12])).should.be.true);
    });

    /***********************************************************************/
    /* Gauss--Jordan Elimination                                           */
    /***********************************************************************/

    interface TestMats {
        M: Matrix,
        rref: Matrix,
        N: Vector[],
        C: Vector[],
        R: Vector[],
        PLU?: PLUData,
        JD?: JordanData,
    };

    describe('#PLU(), #jordanSubst(), #{null,col,row,leftNull}Basis()', () => {
        let testMats: TestMats[] = [
            {M: mat(
                [10,-7,0],
                [-3, 2,6],
                [ 5,-1,5]),
             rref: Matrix.identity(3),
             N: [],
             C: [vec(10,-3,  5),
                 vec(-7, 2, -1),
                 vec( 0, 6,  5)],
             R: [Vector.e(0, 3), Vector.e(1, 3), Vector.e(2, 3)]
            },
            {M: mat(
                [2, 1, 1, 0],
                [4, 3, 3, 1],
                [8, 7, 9, 5],
                [6, 7, 9, 8]),
             rref: Matrix.identity(4),
             N: [],
             C: [vec(2, 4, 8, 6),
                 vec(1, 3, 7, 7),
                 vec(1, 3, 9, 9),
                 vec(0, 1, 5, 8)],
             R: [Vector.e(0, 4), Vector.e(1, 4), Vector.e(2, 4), Vector.e(3, 4)]
            },
            {M: mat(
                [-1, 0, 1],
                [ 2, 1, 1],
                [-1, 2, 0]),
             rref: Matrix.identity(3),
             N: [],
             C: [vec(-1, 2, -1),
                 vec( 0, 1,  2),
                 vec( 1, 1,  0)],
             R: [Vector.e(0, 3), Vector.e(1, 3), Vector.e(2, 3)]
            },
            {M: mat(
                [ 2, -6,  6],
                [-4,  5, -7],
                [ 3,  5, -1],
                [-6,  4, -8],
                [ 8, -3,  9]),
             rref: mat(
                 [ 1, 0,  6/7],
                 [ 0, 1, -5/7],
                 [ 0, 0, 0],
                 [ 0, 0, 0],
                 [ 0, 0, 0]),
             N: [vec(-6/7, 5/7, 1)],
             C: [vec( 2, -4, 3, -6,  8),
                 vec(-6,  5, 5,  4, -3)],
             R: [vec(1, 0,  6/7),
                 vec(0, 1, -5/7)]
            },
            {M: mat(
                [ 0, -3, -6,  4,  9],
                [-1, -2, -1,  3,  1],
                [-2, -3,  0,  3, -1],
                [ 1,  4,  5, -9, -7]),
             rref: mat(
                 [1, 0, -3, 0,  5],
                 [0, 1,  2, 0, -3],
                 [0, 0,  0, 1,  0],
                 [0, 0,  0, 0,  0]),
             N: [vec(3, -2, 1, 0, 0),
                 vec(-5, 3, 0, 0, 1)],
             C: [vec( 0, -1, -2,  1),
                 vec(-3, -2, -3,  4),
                 vec( 4,  3,  3, -9)],
             R: [vec(1, 0, -3, 0,  5),
                 vec(0, 1,  2, 0, -3),
                 vec(0, 0,  0, 1,  0)]
            },
            {M: mat(
                [ 0,  3, -6,  6,  4, -5],
                [ 3, -7,  8, -5,  8,  9],
                [ 3, -9, 12, -9,  6, 15]),
             rref: mat(
                 [ 1, 0, -2, 3, 0, -24],
                 [ 0, 1, -2, 2, 0,  -7],
                 [ 0, 0,  0, 0, 1,   4]),
             N: [vec( 2,  2, 1, 0,  0, 0),
                 vec(-3, -2, 0, 1,  0, 0),
                 vec(24,  7, 0, 0, -4, 1)],
             C: [vec(0,  3,  3),
                 vec(3, -7, -9),
                 vec(4,  8,  6)],
             R: [vec( 1, 0, -2, 3, 0, -24),
                 vec( 0, 1, -2, 2, 0,  -7),
                 vec( 0, 0,  0, 0, 1,   4)]
            },
            {M: mat(
                [1, 2, 3, 4],
                [4, 5, 6, 7],
                [6, 7, 8, 9]),
             rref: mat(
                 [1, 0, -1, -2],
                 [0, 1,  2,  3],
                 [0, 0,  0,  0]),
             N: [vec(1, -2, 1, 0),
                 vec(2, -3, 0, 1)],
             C: [vec(1, 4, 6),
                 vec(2, 5, 7)],
             R: [vec(1, 0, -1, -2),
                 vec(0, 1,  2,  3)]
            },
            {M: mat(
                [1, 3, 5, 7],
                [3, 5, 7, 9],
                [5, 7, 9, 1]),
             rref: mat(
                 [1, 0, -1, 0],
                 [0, 1,  2, 0],
                 [0, 0,  0, 1]),
             N: [vec(1, -2, 1, 0)],
             C: [vec(1, 3, 5),
                 vec(3, 5, 7),
                 vec(7, 9, 1)],
             R: [vec(1, 0, -1, 0),
                 vec(0, 1,  2, 0),
                 vec(0, 0,  0, 1)]
            },
        ];

        it('should factorize PA=LU correctly', () => {
            for(let test of testMats) {
                let { M, C } = test;
                test.PLU = M.PLU();
                test.PLU.should.factorize(M, 'PA=LU');
                M.rank(test.PLU).should.equal(C.length);
            }
        });
        it('should compute the rref and row ops correctly', () => {
            for(let test of testMats) {
                let M = test.M;
                test.JD = M.jordanSubst(test.PLU!);
                const { rref, rowOps } = test.JD;
                rref.equals(test.rref, 1e-10).should.be.true;
                rowOps.mult(M).equals(test.rref, 1e-10).should.be.true;
            }
        });
        it('should compute a basis of the null space', () => {
            for(let { M, N, JD } of testMats)
                M.nullBasis(JD!).should.resemble(N);
        });
        it('should compute a basis of the column space', () => {
            for(let { M, C, PLU } of testMats)
                M.colBasis(PLU!).should.resemble(C);
        });
        it('should compute a basis of the row space', () => {
            for(let { M, N, PLU } of testMats) {
                let R = M.rowBasis(PLU!);
                R.length.should.equal(PLU!.pivots.length);
                let RR = Matrix.create(...R);
                // The row space is orthogonal to the null space
                N.every(v => RR.apply(v).isZero(1e-10)).should.be.true;
            }
        });
        it('should compute a basis of the left null space', () => {
            for(let { M, PLU } of testMats) {
                let B = M.leftNullBasis(PLU!);
                B.length.should.equal(M.m - PLU!.pivots.length);
                // The left null space is orthogonal to the column space
                let MT = M.transpose;
                B.every(v => MT.apply(v).isZero(1e-10)).should.be.true;
            }
        });
    });
    describe('#solveReverseSubst()', () => {
        // Test reverse-substitution in echelon form
        let testRevSub = [
            {M: mat([ 1, 0,  6/7],
                    [ 0, 1, -5/7],
                    [ 0, 0, 0],
                    [ 0, 0, 0],
                    [ 0, 0, 0]),
             b: vec(7, 3, 0, 0, 0),
             b2: vec(0, 0, 0, 0, 1)},
            {M: mat([1, 0, -3, 0,  5],
                    [0, 1,  2, 0, -3],
                    [0, 0,  0, 1,  0],
                    [0, 0,  0, 0,  0]),
             b: vec(7, 2, 1, 0),
             b2: vec(0, 0, 0, 1)}
        ];
        it('should solve Mx=b using reverse-substitution', () => {
            for(let {M, b, b2} of testRevSub) {
                expect(M.solveReverseSubst(b)).to.solve(M, b);
                expect(M.solveReverseSubst(b2)).to.be.null;
            }
        });
    });
    describe('#solve(), #solveLeastSquares()', () => {
        let testMats = [
            {M: mat([10,-7,0],
                    [-3, 2,6],
                    [ 5,-1,5]),
             b: vec(1,2,3)
            },
            {M: mat([2, 1, 1, 0],
                    [4, 3, 3, 1],
                    [8, 7, 9, 5],
                    [6, 7, 9, 8]),
             b: vec(1,2,3,4)
            },
            {M: mat([ 2, -6,  6],
                    [-4,  5, -7],
                    [ 3,  5, -1],
                    [-6,  4, -8],
                    [ 8, -3,  9]),
             b: vec(8, -15, 10, -22, 29)
            },
            {M: mat([ 0, -3, -6,  4,  9],
                    [-1, -2, -1,  3,  1],
                    [-2, -3,  0,  3, -1],
                    [ 1,  4,  5, -9, -7]),
             b: vec(37, 9, -1, -47)
            }
        ];
        it('should solve Mx=b', () => {
            for(let {M, b} of testMats) {
                expect(M.solve(b)).to.solve(M, b);
                /* let x = M.solveShortest(b); */
                /* x.should.solve(M, b); */
                /* M.nullSpace().isOrthogonalTo(x).should.be.true(); */
                M.solveLeastSquares(b).should.solve(M, b);
                M.solveLeastSquares(b, {QR: M.QR()}).should.solve(M, b);
                /* M.solveLeastSquaresShortest(b).equals(x, 1e-10).should.be.true(); */
                /* M.pseudoInverse(); */
                /* M.solveLeastSquaresShortest(b).equals(x, 1e-10).should.be.true(); */
            }
        });
        let testNoSoln = [
            {M: mat([ 2, -6,  6],
                    [-4,  5, -7],
                    [ 3,  5, -1],
                    [-6,  4, -8],
                    [ 8, -3,  9]),
             b: vec(0, -15, 10, -22, 29)
            },
            {M: mat([ 0, -3, -6,  4,  9],
                    [-1, -2, -1,  3,  1],
                    [-2, -3,  0,  3, -1],
                    [ 1,  4,  5, -9, -7]),
             b: vec(37, 0, -1, -47)
            }
        ];
        it('should have no solution when inconsistent', () => {
            for(let {M, b} of testNoSoln)
                expect(M.solve(b)).to.be.null;
        });
        it('should have a least-squares solution', () => {
            for(let {M, b} of testNoSoln) {
                M.solveLeastSquares(b).should.solve(M, b, true);
                M.solveLeastSquares(b, {QR: M.QR()}).should.solve(M, b, true);
                /* let x = M.solveLeastSquaresShortest(b); */
                /* x.should.solve(M, b, true); */
                /* M.nullSpace().isOrthogonalTo(x).should.be.true(); */
                /* M.pseudoInverse(); */
                /* M.solveLeastSquaresShortest(b).equals(x, 1e-10) */
                /*     .should.be.true(); */
            }
        });
        it('should throw for incompatible dimensions', () => {
            let M = mat([1,2],[3,4]);
            M.solve.bind(M, vec(1,2,3)).should.throw(/Incompatible/);
        });
    });

    /***********************************************************************/
    /* Orthogonality                                                       */
    /***********************************************************************/

    describe('#QR()', () => {
        let testMats = [
            mat([ 3, -5,  1],
                [ 1,  1,  1],
                [-1,  5, -2],
                [ 3, -7,  8]),
            mat([ 1,  2,  5],
                [-1,  1, -4],
                [-1,  4, -3],
                [ 1, -4,  7],
                [ 1,  2,  1]),
            mat([ 0,  3, -6,  6,  4, -5],
                [ 3, -7,  8, -5,  8,  9],
                [ 3, -9, 12, -9,  6, 15]).transpose,
            mat([-10,  13,  7, -11],
                [  2,   1, -5,   3],
                [ -6,   3, 13,  -3],
                [ 16, -16, -2,   5],
                [  2,   1, -5,  -7])
        ];
        it('should factorize matrices correctly', () => {
            for(let M of testMats) {
                let QR = M.QR();
                QR.should.factorize(M, 'QR', 1e-15);
                M.rank(QR).should.equal(M.rank());
            }
        });
        let testMats2 = [
            mat([ 0, -3, -6,  4,  9],
                [-1, -2, -1,  3,  1],
                [-2, -3,  0,  3, -1],
                [ 1,  4,  5, -9, -7]),
            mat([ 2, -6,  6],
                [-4,  5, -7],
                [ 3,  5, -1],
                [-6,  4, -8],
                [ 8, -3,  9])
        ];
        it('should factorize matrices with linearly dependent columns', () => {
            for(let M of testMats2) {
                let QR = M.QR();
                QR.should.factorize(M, 'QR');
                M.rank(QR).should.equal(M.rank());
            }
        });
    });
});


describe('SquareMatrix', () => {

    /***********************************************************************/
    /* Traits                                                              */
    /***********************************************************************/

    describe('#isOrthogonal()', () => {
        it('detects orthogonal matrices', () =>
            smat([1,1],[1,-1]).scale(1/Math.sqrt(2)).isOrthogonal()
                .should.be.true);
        it('detects non-orthogonal matrices', () =>
            smat([1/Math.sqrt(2), 1],[1/Math.sqrt(2), 0]).isOrthogonal()
                .should.be.false);
    });
    describe('#isSymmetric()', () => {
        it('detects symmetric matrices', () =>
            smat([1,2,3],[2,4,5],[3,5,6]).isSymmetric().should.be.true);
        it('detects non-symmetric square matrices', () =>
            smat([1,2,3],[4,5,6],[7,8,9]).isSymmetric().should.be.false);
    });

    /***********************************************************************/
    /* Solving Ax=b                                                        */
    /***********************************************************************/

    describe('#solveForwardSubst(), #solveReverseSubst1()', () => {
        it('should solve Ax=b using forward-substitution', () => {
            let M = smat([1, 0, 0, 0],
                         [2, 1, 0, 0],
                         [3, 2, 1, 0],
                         [4, 3, 2, 1]);
            let b = vec(7, 6, 5, 4);
            M.solveForwardSubst(b).should.solve(M, b);
        });
        it('should solve Ax=b using reverse-substitution', () => {
            let M = smat([2, 1, 3, 4],
                         [0, 5, 1, 1],
                         [0, 0, 1, 1],
                         [0, 0, 0, 2]);
            let b = vec(7, 6, 5, 4);
            M.solveReverseSubst1(b).should.solve(M, b);
        });
    });

    /***********************************************************************/
    /* Determinants and Eigenvalues                                        */
    /***********************************************************************/

    describe('#det', () => {
        it('should compute the determinant (1x1)', () =>
           Matrix.identity(1, 3).det().should.equal(3));
        it('should compute the determinant (2x2)', () =>
           smat([3, 4], [5, 6]).det().should.be.approximately(-2, 1e-10));
        it('should compute the determinant (3x3#1)', () =>
           smat([0,  1, 2],
                [1,  0, 3],
                [4, -3, 8]).det().should.be.approximately(-2, 1e-10));
        it('should compute the determinant (3x3#2)', () =>
           smat([ 1, -2, -1],
                [-1,  5,  6],
                [ 5, -4,  5]).det().should.equal(0));
        it('should compute the determinant (4x4)', () =>
           smat([ 1,  7,  4, 2],
                [ 3, 11,  9, 5],
                [-2, -3,  3, 3],
                [ 7,  8, -8, 9]).det().should.be.approximately(-1329, 1e-10));
    });
    describe('#inverse(), #adjugate()', () => {
        let testMats = [
            smat([3, 4],
                 [5, 6]),
            smat([0,  1, 2],
                 [1,  0, 3],
                 [4, -3, 8]),
            smat([ 1,  7,  4, 2],
                 [ 3, 11,  9, 5],
                 [-2, -3,  3, 3],
                 [ 7,  8, -8, 9])
        ];
        it('should compute the inverse', () => {
            for(let M of testMats) {
                let I = M.inverse();
                expect(I).not.to.be.null;
                M.mult(I!).equals(Matrix.identity(M.n), 1e-10).should.be.true;
            }
        });
        it('should compute the adjugate', () => {
            for(let M of testMats) {
                M.inverse()!.scale(M.det()).equals(M.adjugate, 1e-10)
                    .should.be.true;
            }
        });
        it('should return null for singular matrices', () => {
            let M = smat([ 1, -2, -1],
                         [-1,  5,  6],
                         [ 5, -4,  5]);
            expect(M.inverse()).to.be.null;
        });

    });
    describe('#charpoly', () => {
        it('should compute the characteristic polynomial of a 1x1 matrix', () => {
            smat([3]).charpoly.should.eql(poly(-1, 3));
        });
        it('should compute the characteristic polynomial of a 2x2 matrix', () => {
            smat([1,2],[3,4]).charpoly.should.eql(poly(1, -5, -2));
        });
        it('should compute the characteristic polynomial of a 3x3 matrix', () => {
            smat([1, 6,4],
                 [2,-1,3],
                 [5, 0,1]).charpoly.should.eql(poly(-1, 1, 33, 97));
        });
        it('should compute the characteristic polynomial of a 4x4 matrix', () => {
            smat([2, 1, 1, 0],
                 [4, 3, 3, 1],
                 [8, 7, 9, 5],
                 [6, 7, 9, 8]).charpoly.should.eql(poly(1, -22, 78, -50, 8));
        });

    });
    describe('#eigenvalues()', () => {
        it('should compute eigenvalues for 1x1 matrices', () =>
            smat([3]).eigenvalues().should.eql([[3, 1]]));
        it('should compute eigenvalues for 2x2 matrices', () => {
            smat([1,1],[1,1]).eigenvalues().should.eql([[0, 1], [2, 1]]);
            smat([1,1],[0,1]).eigenvalues().should.eql([[1, 2]]);
            smat([1,1],[-1,1]).eigenvalues()
                .should.resemble([[C(1, 1), 1], [C(1, -1), 1]]);
        });
        it('should compute eigenvalues for 3x3 matrices', () => {
            Matrix.identity(3, 3).eigenvalues().should.resemble([[3, 3]]);
            smat([11/13, 22/39,  2/39],
                 [-4/13, 83/39,  4/39],
                 [-1/13, 11/39, 40/39]).eigenvalues(1e-7)
                     .should.resemble([[1, 2], [2, 1]]);
            smat([    1,   1/2,     0],
                 [-4/13, 83/39,  4/39],
                 [ 5/13,  7/78, 34/39]).eigenvalues(1e-7)
                     .should.resemble([[1, 2], [2, 1]]);
            smat([23/13,  53/78, -10/39],
                 [-4/13, 122/39,   4/39],
                 [-4/13,  49/78,  43/39]).eigenvalues()
                     .should.resemble([1, 2, 3].map(x => [x, 1]));
            smat([43/13,   7/13, -88/13],
                 [ 6/13,  17/13, -28/13],
                 [24/13, -23/13,  57/13]).eigenvalues()
                     .should.resemble([1, C(4, 3), C(4, -3)].map(x => [x, 1]));
        });
        it('should compute eigenvalues for 4x4 matrices', () => {
            Matrix.identity(4, 3).eigenvalues().should.resemble([[3, 4]]);
            smat([-140,    325,   -730,   -964],
                 [  88,   -202,    457,    607],
                 [ 127, -585/2, 1317/2, 1741/2],
                 [ -45,  207/2, -465/2, -613/2]).eigenvalues()
                     .should.resemble([1, 2, 3, 4].map(x => [x, 1]));
            smat([-95,    220,   -490,   -634],
                 [ 22,    -48,    105,    123],
                 [ 73, -333/2,  741/2,  949/2],
                 [-33,  151/2, -337/2, -437/2]).eigenvalues()
                     .should.resemble([[1, 2], [3, 1], [4, 1]]);
            smat([-635,    1462,   -3298,   -4414],
                 [ 282,    -646,    1457,    1943],
                 [ 453, -2081/2,  4693/2,  6269/2],
                 [-153,   703/2, -1585/2, -2117/2]).eigenvalues()
                     .should.resemble([[1, 3], [4, 1]]);
            smat([-1229,  5651/2, -12725/2, -16955/2],
                 [  568, -2605/2,   5865/2,   7799/2],
                 [  871,   -2000,     4503,     5994],
                 [ -285,  1309/2,  -2947/2,  -3923/2]).eigenvalues()
                     .should.resemble([[1, 2], [4, 2]]);
            smat([-76, 345/2, -771/2, -983/2],
                 [ 10,   -18,     39,     33],
                 [ 60,  -134,    299,    377],
                 [-30,    68,   -152,   -196]).eigenvalues()
                     .should.resemble([C(1,1), C(1,-1), 3, 4].map(x => [x, 1]));
            smat([-1272,    2921,   -6584,   -8783],
                 [  582, -2665/2,  6007/2,  7997/2],
                 [  892,   -2046,    4611,    6145],
                 [ -290,  1331/2, -2999/2, -3997/2]).eigenvalues()
                     .should.resemble([[C(1,1), 1], [C(1,-1), 1], [3, 2]]);
            smat([ 3440,  -7907,    17832,    23865],
                 [ -722, 3329/2,  -7515/2, -10105/2],
                 [-1294, 5955/2, -13435/2, -18013/2],
                 [  232,   -534,     1205,     1617]).eigenvalues()
                     .should.resemble([C(1, 1), C(1, -1)].map(x => [x, 2]));
            smat([1996, -9177/2, 20695/2, 27703/2],
                 [ -26,      65,    -150,    -226],
                 [-276,     638,   -1441,   -1947],
                 [ -90,     206,    -464,    -616]).eigenvalues()
                     .should.resemble([C(1, 1), C(1, -1), C(1, 2), C(1, -2)]
                                          .map(x => [x, 1]));
        });
    });
    describe('#eigenspaceBasis()', () => {
        it('should work for 1x1 matrices', () =>
            smat([3]).eigenspaceBasis(3).should.have.eigenvectors(smat([3]), 3, 1));
        it('should work for 2x2 matrices with distinct eigenvalues', () => {
            let M = smat([1, 1], [1, 1]);
            M.eigenspaceBasis(0).should.have.eigenvectors(M, 0, 1);
            M.eigenspaceBasis(2).should.have.eigenvectors(M, 2, 1);
            M.eigenspaceBasis.bind(M, 2.001).should.throw(/not an eigenvalue/);
        });
        it('should work for 2x2 matrices with a complex eigenvalue', () => {
            let M = smat([-15/7, 52/7], [-40/7, 57/7]);
            let λ = C(3, 4);
            let B = M.eigenspaceBasis(λ);
            B.should.have.length(1);
            B.should.have.eigenvectors(M, λ, 1);
            λ.conj();
            B = M.eigenspaceBasis(λ);
            B.should.have.length(1);
            B.should.have.eigenvectors(M, λ, 1);
        });
        it('should work for 2x2 matrices with a repeated eigenvalue', () => {
            let M = smat([1,1],[0,1]);
            let B = M.eigenspaceBasis(1);
            B.should.have.eigenvectors(M, 1, 1);
            M = smat([2,0],[0,2]);
            B = M.eigenspaceBasis(2);
            B.should.have.eigenvectors(M, 2, 2);
        });
        it('should work for 3x3 matrices with distinct eigenvalues', () => {
            let M = smat([23/13,  53/78, -10/39],
                         [-4/13, 122/39,   4/39],
                         [-4/13,  49/78,  43/39]);
            M.eigenspaceBasis(1).should.have.eigenvectors(M, 1, 1);
            M.eigenspaceBasis(2).should.have.eigenvectors(M, 2, 1);
            M.eigenspaceBasis(3).should.have.eigenvectors(M, 3, 1);
        });
        it('should work for 3x3 matrices with repeated eigenvalues', () => {
            let M = smat([11/13, 22/39,  2/39],
                         [-4/13, 83/39,  4/39],
                         [-1/13, 11/39, 40/39]);
            M.eigenspaceBasis(1).should.have.eigenvectors(M, 1, 2);
            M.eigenspaceBasis(2).should.have.eigenvectors(M, 2, 1);
            M = smat([    1,   1/2,     0],
                     [-4/13, 83/39,  4/39],
                     [ 5/13,  7/78, 34/39]);
            M.eigenspaceBasis(1).should.have.eigenvectors(M, 1, 1);
            M.eigenspaceBasis(2).should.have.eigenvectors(M, 2, 1);
            M = Matrix.identity(3, 3);
            M.eigenspaceBasis(3).should.have.eigenvectors(M, 3, 3);
        });
        it('should work for 3x3 matrices with a complex eigenvalue', () => {
            let M = smat([33, -23,   9],
                         [22,  33, -23],
                         [19,  14,  50]).scale(1/29);
            M.eigenspaceBasis(2).should.have.eigenvectors(M, 2, 1);
            let λ = C(1, 1);
            M.eigenspaceBasis(λ).should.have.eigenvectors(M, λ, 1);
            λ.conj();
            M.eigenspaceBasis(λ).should.have.eigenvectors(M, λ, 1);
        });
        it('should work for 4x4 matrices with distinct complex eigenvalues', () => {
            let M = smat([-1226,   230,  1166, -989],
                         [ 1530,  -192,  -786,  932],
                         [ 8938, -1856, -6401, 6438],
                         [12222, -2471, -9114, 8883]).scale(1/133);
            let λ = C(3, 4);
            M.eigenspaceBasis(λ).should.have.eigenvectors(M, λ, 1);
            λ.conj();
            M.eigenspaceBasis(λ).should.have.eigenvectors(M, λ, 1);
            λ = C(1, 1);
            M.eigenspaceBasis(λ).should.have.eigenvectors(M, λ, 1);
            λ.conj();
            M.eigenspaceBasis(λ).should.have.eigenvectors(M, λ, 1);
        });
        it('should work for 4x4 matrices with repeated complex eigenvalues', () => {
            let M = smat([ -213,   396,   208, -160],
                         [-2616,  1059,  1056, -976],
                         [10256, -2036, -7221, 7212],
                         [10968, -2528, -8088, 7971]).scale(1/133);
            let λ = C(3, 4);
            M.eigenspaceBasis(λ).should.have.eigenvectors(M, λ, 2);
            λ.conj();
            M.eigenspaceBasis(λ).should.have.eigenvectors(M, λ, 2);

            M = smat([ -126,   371,   168, -119],
                     [-2442,  1009,   976, -894],
                     [10604, -2136, -7381, 7376],
                     [11229, -2603, -8208, 8094]).scale(1/133);
            λ = C(3, 4);
            M.eigenspaceBasis(λ).should.have.eigenvectors(M, λ, 1);
            λ.conj();
            M.eigenspaceBasis(λ).should.have.eigenvectors(M, λ, 1);
        });
        it('should work for 4x4 matrices with real and complex eigenvalues', () => {
            let M = smat([  276,    17,   240, -113],
                         [ 1419,   138, -1056, 1109],
                         [ 9654, -2022, -6773, 6806],
                         [10827, -2249, -8280, 8088]).scale(1/133);
            let λ = C(3, 4);
            M.eigenspaceBasis(λ).should.have.eigenvectors(M, λ, 1);
            λ.conj();
            M.eigenspaceBasis(λ).should.have.eigenvectors(M, λ, 1);
            M.eigenspaceBasis(3).should.have.eigenvectors(M, 3, 1);
            M.eigenspaceBasis(4).should.have.eigenvectors(M, 4, 1);

            M = smat([  363,    -8,   200,  -72],
                     [ 2463,  -162, -1536, 1601],
                     [ 9480, -1972, -6693, 6724],
                     [10827, -2249, -8280, 8088]).scale(1/133);
            λ = C(3, 4);
            M.eigenspaceBasis(λ).should.have.eigenvectors(M, λ, 1);
            λ.conj();
            M.eigenspaceBasis(λ).should.have.eigenvectors(M, λ, 1);
            M.eigenspaceBasis(3).should.have.eigenvectors(M, 3, 1);

            M.eigenspaceBasis.bind(M, C(3, 5)).should.throw(/not an eigenvalue/);
        });
    });
    describe('#diagonalize(), #isDiagonalizable()', () => {
        it('should work for 1x1 matrices', () => {
            let M = smat([3]);
            expect(M.diagonalize()).to.diagonalize(M);
            expect(M.diagonalize({block: true})).to.diagonalize(M);
        });
        it('should diagonalize 2x2 matrices', () => {
            let testMats = [
                smat([1, 1], [1, 1]),
                smat([2, 0], [0, 2])
            ];
            for(let M of testMats) {
                expect(M.diagonalize()).to.diagonalize(M);
                expect(M.diagonalize({block: true})).to.diagonalize(M);
            }
            expect(smat([1,1],[0,1]).diagonalize()).to.be.null;
        });
        it('should diagonalize 3x3 matrices', () => {
            let testMats = [
                smat([23/13,  53/78, -10/39],
                     [-4/13, 122/39,   4/39],
                     [-4/13,  49/78,  43/39]),
                smat([11/13, 22/39,  2/39],
                     [-4/13, 83/39,  4/39],
                     [-1/13, 11/39, 40/39]),
                Matrix.identity(3, 3)
            ];
            for(let M of testMats) {
                expect(M.diagonalize()).to.diagonalize(M);
                expect(M.diagonalize({block: true})).to.diagonalize(M);
            }
            expect(smat([    1,   1/2,     0],
                        [-4/13, 83/39,  4/39],
                        [ 5/13,  7/78, 34/39]).diagonalize()).to.be.null;
        });
        it('should diagonalize 4x4 matrices', () => {
            let testMats = [
                smat([1482, -247, -741,  703],
                     [ 963,  214, -828,  729],
                     [1572, -360, -842, 1016],
                     [ 105,  -21, -273,  476]).scale(1/133),
                smat([ 374, 132, 212, -186],
                     [ -18,  76, -12,   38],
                     [-132, -40, -32,   92],
                     [ 198, 116, 132,  -26]).scale(1/56),
                smat([344, 128, 192, -160],
                     [ 72,  88,  48,  -40],
                     [-72, -32,   8,   40],
                     [288, 128, 192, -104]).scale(1/56),
                Matrix.identity(4, 3),
            ];
            for(let M of testMats) {
                expect(M.diagonalize()).to.diagonalize(M);
                expect(M.diagonalize({block: true})).to.diagonalize(M);
            }
        });
        it('should block-diagonalize 2x2 matrices with a complex eigenvalue', () => {
            let testMats = [
                smat([2, -1], [2, 0]),
                smat([-Math.sqrt(3)+1, -2], [1, -Math.sqrt(3)-1])
            ];
            for(let M of testMats) {
                expect(M.diagonalize()).to.be.null;
                let CD = M.diagonalize({block: true});
                expect(CD).to.diagonalize(M, true);
                let D = CD!.D;
                D[0][0].should.be.approximately(D[1][1], 1e-10);
                D[0][1].should.be.approximately(-D[1][0], 1e-10);
            }
        });
        it('should block-diagonalize 3x3 matrices with a complex eigenvalue', () => {
            let M = smat([33, -23,   9],
                         [22,  33, -23],
                         [19,  14,  50]).scale(1/29);
            expect(M.diagonalize()).to.be.null;
            let CD = M.diagonalize({block: true});
            expect(CD).to.diagonalize(M, true);
            let D = CD!.D;
            D[0][2].should.equal(0);
            D[1][2].should.equal(0);
            D[2][0].should.equal(0);
            D[2][1].should.equal(0);
            D[0][0].should.be.approximately(D[1][1], 1e-10);
            D[0][1].should.be.approximately(-D[1][0], 1e-10);
        });
        it('should block-diagonalize 4x4 matrices'
            + 'with repeated complex eigenvalues', () => {
                let M = smat([ -213,   396,   208, -160],
                             [-2616,  1059,  1056, -976],
                             [10256, -2036, -7221, 7212],
                             [10968, -2528, -8088, 7971]).scale(1/133);
                expect(M.diagonalize()).to.be.null;
                expect(M.diagonalize({block: true, ε: 1e-8})).to.diagonalize(M, true);

                M = smat([ -126,   371,   168, -119],
                         [-2442,  1009,   976, -894],
                         [10604, -2136, -7381, 7376],
                         [11229, -2603, -8208, 8094]).scale(1/133);
                expect(M.diagonalize()).to.be.null;
                expect(M.diagonalize({block: true, ε: 1e-8})).to.be.null;
            });
        it('should block-diagonalize 4x4 matrices'
           + ' with real and complex eigenvalues', () => {
               let M = smat([  276,    17,   240, -113],
                            [ 1419,   138, -1056, 1109],
                            [ 9654, -2022, -6773, 6806],
                            [10827, -2249, -8280, 8088]).scale(1/133);
               expect(M.diagonalize()).to.be.null;
               expect(M.diagonalize({block: true})).to.diagonalize(M, true);

               M = smat([ 1360,  480,  1504, -1048],
                        [-2288, -824, -2608,  1856],
                        [ 1476,  600,  1712, -1212],
                        [ 2288,  880,  2608, -1800]).scale(1/56);
               expect(M.diagonalize()).to.be.null;
               expect(M.diagonalize({block: true})).to.diagonalize(M, true);
        });
        it('should orthogonally diagonalize symmetric matrices', () => {
            let testMats = [
                smat([3]),
                smat([1,2], [2,1]),
                smat([4,3], [3,4]),
                smat([1,2,3],
                     [2,4,7],
                     [3,7,9]),
                smat([-1, 5, 4],
                     [ 5,-2, 3],
                     [ 4, 3,-9]),
                smat([1.50199300338977629,
                      0.44200755514697981,
                      0.23372066474849376],
                     [0.44200755514698000,
                      1.38919004346224608,
                      0.20579230968403700],
                     [0.23372066474849394,
                      0.20579230968403707,
                      1.10881695314797761])
            ];
            for(let M of testMats) {
                let CD = M.diagonalize({ortho: true});
                expect(CD).to.diagonalize(M);
                CD!.C.isOrthogonal().should.be.true;
            }
            let M = testMats.pop()!;
            expect(M.eigenvalues()).to.resemble([[1, 2], [2, 1]]);
            M.eigenspaceBasis(1).length.should.equal(2);
        });
    });
});

/* // Test LDLT solving for positive-definite symmetric matrices */
/* let testLDLT = [ */
/*     {M: mat([10,-7,0], */
/*             [-3, 2,6], */
/*             [ 5,-1,5]), */
/*      b: vec(1, 2, 3)}, */
/*     {M: mat([2, 1, 1, 0], */
/*             [4, 3, 3, 1], */
/*             [8, 7, 9, 5], */
/*             [6, 7, 9, 8]), */
/*      b: vec(1, 2, 3, 4)}, */
/*     {M: mat([ 0,  3, -6,  6,  4, -5], */
/*             [ 3, -7,  8, -5,  8,  9], */
/*             [ 3, -9, 12, -9,  6, 15]).transpose, */
/*      b: vec(1, 2, 3)}, */
/*     {M: mat([1, 3, 5, 7], */
/*             [3, 5, 7, 9], */
/*             [5, 7, 9, 1]).transpose, */
/*      b: vec(1, 2, 3)}, */
/*     {M: mat([1,2,3],[2,7,4],[3,4,8]), */
/*      b: vec(1, 2, 3)} */
/* ]; */
/* it('should solve Mx=b using LDLT', () => { */
/*     for(let {M, b} of testLDLT) { */
/*         if(!M.isSquare() || !M.isSymmetric()) */
/*             M = M.normal; */
/*         M.solve(b).should.solve(M, b); */
/*     } */
/* }); */
/* // Test forward-substitution in lower-unitriangular form */
