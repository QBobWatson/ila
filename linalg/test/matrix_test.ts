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

import Vector from "../src/vector";
import Matrix from "../src/matrix3";
import { SquareMatrix, PLUData, QRData, JordanData } from "../src/matrix3";

const mat = Matrix.create;
const vec = (...entries: number[]) => new Vector(...entries);


type Factorization = 'PA=LU' | 'QR' | 'SVD' | 'LDLT';

Assertion.addMethod(
    'factorize', function(this: typeof Assertion, A: Matrix,
                          which: Factorization, ε: number=1e-10) {
        switch (which) {
            case 'PA=LU':
                let { P, L, U, E, pivots } = this._obj as PLUData;
                expect(L.isLowerUni(), `L = \n${L.toString(2)}\nis not lower unipotent`).to.be.true;
                expect(U.isREF(), `U =\n${U.toString(2)}\nis not in REF`).to.be.true;
                expect(U.leadingEntries()).to.eql(pivots);
                let PA = Matrix.permutation(P).mult(A);
                let LU = L.mult(U);
                expect(PA.equals(LU, ε),
                       `Matrix is not correctly factored: PA =\n${PA.toString(2)}\nLU =\n${LU.toString(2)}`)
                        .to.be.true;
                let EA = E.mult(A);
                expect(EA.equals(U, ε),
                       `Matrix is not correctly factored: EA =\n${EA.toString(2)}\nU =\n${U.toString(2)}`)
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
                       `Q does not have orthogonal columns: Q=\n${Q.toString(2)}\nQTQ=\n${QTQ.toString(2)}`)
                        .to.be.true;
                for(let j = 0; j < QTQ.n; ++j) {
                    if(LD.includes(j)) {
                        expect(Q.col(j).isZero(ε),
                               `Column ${j} of Q should be zero: LD=${LD}, Q=\n${Q.toString(2)}`)
                            .to.be.true;
                        expect(R.row(j).isZero(ε),
                               `Row ${j} of R should be zero: LD=${LD}, R=\n${R.toString(2)}`)
                            .to.be.true;
                    } else {
                        expect(QTQ[j][j]).to.be.approximately(1, ε);
                        expect(R[j][j]).to.be.above(ε);
                    }
                }
                let QR = Q.mult(R);
                expect(QR.equals(A, ε),
                       `Matrix is not correctly factored: QR =\n${QR.toString(2)}\nA =\n${A.toString(2)}`)
                        .to.be.true;
                break;
        }
    });

Assertion.addMethod(
    'solve', function(this: typeof Assertion, M: Matrix, b: Vector,
                      leastSquares: boolean=false, ε: number=1e-10) {
        let x = this._obj;
        expect(x, 'A solution was not found when one was expected').not.to.be.null;
        let Mx = M.apply(x);
        if(!leastSquares)
            expect(Mx.equals(b, ε),
                   `Mx != b:\nx = ${x.toString(2)}\nMx = ${Mx.toString(2)}\nb = ${b.toString(2)}`)
                .to.be.true;
        else {
            expect(M.transpose.apply(Mx.sub(b)).isZero(1e-10),
                   'Expected Mx-b to be orthogonal to Col(M)')
                .to.be.true;
        }
    });

declare global {
    export namespace Chai {
        interface Assertion {
            factorize(A: Matrix, which: Factorization, ε?: number): void;
            solve(M: Matrix, b: Vector, leastSquares?: boolean, ε?: number): void;
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
                let M = test.M;
                test.PLU = M.PLU();
                test.PLU.should.factorize(M, 'PA=LU');
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
            for(let M of testMats)
                M.QR().should.factorize(M, 'QR', 1e-15);
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
            for(let M of testMats2)
                M.QR().should.factorize(M, 'QR');
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
/* it('should solve Ax=b using forward-substitution', () => { */
/*     let M = mat([1, 0, 0, 0], [2, 1, 0, 0], [3, 2, 1, 0], [4, 3, 2, 1]); */
/*     let b = vec(7, 6, 5, 4); */
/*     M.solve(b).should.solve(M, b); */
/* }); */
