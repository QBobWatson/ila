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

import Subspace from "../src/subspace";
import Matrix from "../src/matrix3";
import Vector from "../src/vector";

const mat = Matrix.create;
const vec = (...entries: number[]) => new Vector(...entries);
const spanOf = Subspace.spanOf;


Assertion.addMethod(
    'decompose', function(this: typeof Assertion, vec: Vector,
                          V: Subspace, ε=1e-10) {
        let [v1, v2] = this._obj;
        expect(v1.clone().add(v2).equals(vec, 1e-10)).to.be.true;
        expect(V.contains(v1, ε)).to.be.true;
        expect(V.QT!.apply(v2)).to.resemble(Vector.zero(V.QT!.m));
    });

declare global {
    export namespace Chai {
        interface Assertion {
            decompose(vec: Vector, v: Subspace, ε?: number): void;
        }
    }
}


describe('Subspace', () => {
    describe('#spanOf()', () => {
        it('should construct a subspace from an Array of number[]', () => {
            let ns = [[1,2,3], [4,5,6]];
            let V = spanOf(ns);
            V.dim.should.equal(2);
            V.n.should.equal(3);
            for(let v of ns)
                V.contains(v).should.be.true;
        });
        it('should construct a subspace from an Array of Vector', () => {
            let vs = [vec(2,1,0), vec(-4,-2,0)];
            let V = spanOf(vs);
            V.dim.should.equal(1);
            V.n.should.equal(3);
            for(let v of vs)
                V.contains(v).should.be.true;
        });
        it('should construct a subspace from a Matrix', () => {
            let M = mat(
                [ 2, -6,  6],
                [-4,  5, -7],
                [ 3,  5, -1],
                [-6,  4, -8],
                [ 8, -3,  9]);
            let V = spanOf(M);
            V.dim.should.equal(2);
            V.n.should.equal(5);
            for(let v of M.cols)
                V.contains(v).should.be.true;
        });
        it('should construct a subspace from a Matrix with orthogonal columns',
           () => {
               let M = mat([1,1],[1,-1],[1,0])
                   .mult(Matrix.diagonal([1/Math.sqrt(3), 1/Math.sqrt(2)]));
               let V = spanOf(M, {isON: true});
               V.dim.should.equal(2);
               V.n.should.equal(3);
               expect(V.basis).not.to.be.null;
               V.basis!.should.equal(M);
           });
        it('should construct a the zero subspace', () => {
            let V = spanOf([], {n: 3});
            V.dim.should.equal(0);
            V.n.should.equal(3);
        });
        it('should throw if it can\'t determine n', () =>
            expect(() => spanOf([])).to.throw(/ambient/));
    });
    describe('#nullSpaceOf()', () => {
        it('should construct the zero subspace correctly', () => {
            let M = mat(
                [2, 1, 1, 0],
                [4, 3, 3, 1],
                [8, 7, 9, 5],
                [6, 7, 9, 8]);
            let V = Subspace.nullSpaceOf(M);
            expect(V.basis).to.be.null;
            expect(V.dim).to.equal(0);
        });
        it('should construct a nonzero null space correctly', () => {
            let M = mat(
                [ 0, -3, -6,  4,  9],
                [-1, -2, -1,  3,  1],
                [-2, -3,  0,  3, -1],
                [ 1,  4,  5, -9, -7]);
            let V = Subspace.nullSpaceOf(M);
            expect(V.dim).to.equal(2);
        });
    });
    describe('#Rn(), #zero()', () => {
        it('should construct all of Rn', () => {
            let V = Subspace.Rn(3);
            V.n.should.equal(3);
            expect(V.basis).not.to.be.null;
            V.basis!.should.eql(Matrix.identity(3));
            V.dim.should.equal(3);
        });
        it('should construct the zero space', () => {
            let V = Subspace.zero(3);
            V.n.should.equal(3);
            expect(V.basis).to.be.null;
            V.dim.should.equal(0);
        });
    });
    describe('#clone()', () => {
        it('should clone a subspace', () => {
            let V = spanOf([[1,2,3], [4,5,6]]);
            let W = V.clone();
            W.n.should.equal(V.n);
            expect(V.basis).to.eql(W.basis);
            expect(V.basis).not.to.equal(W.basis);
        });
    });
    describe('#toString()', () => {
        it('should describe all of R^n', () =>
            Subspace.Rn(3).toString()
                .should.eql("The full subspace R^3"));
        it('should describe the zero subspace', () =>
            Subspace.zero(3).toString()
                .should.eql("The zero subspace of R^3"));
        let A = mat([ 0, -3, -6,  4,  9],
                    [-1, -2, -1,  3,  1],
                    [-2, -3,  0,  3, -1],
                    [ 1,  4,  5, -9, -7]);
        let V = spanOf(A);
        it('should have 4 decimal places by default', () =>
            V.toString().should.eql("Subspace of R^4 of dimension 3 with basis\n"
                + "[ 0.0000] [-0.8018] [-0.5976]\n"
                + "[-0.4082] [ 0.0000] [ 0.0000]\n"
                + "[-0.8165] [ 0.2673] [-0.3586]\n"
                + "[ 0.4082] [ 0.5345] [-0.7171]"));
        it('can have other precision', () =>
            V.toString(2).should.eql("Subspace of R^4 of dimension 3 with basis\n"
                + "[ 0.00] [-0.80] [-0.60]\n"
                + "[-0.41] [ 0.00] [ 0.00]\n"
                + "[-0.82] [ 0.27] [-0.36]\n"
                + "[ 0.41] [ 0.53] [-0.72]"));
    });
    describe('#dim, #isMaximal(), #isZero()', () => {
        it('should construct a full subspace', () => {
            let V = spanOf(mat([ 0,  3, -6,  6,  4, -5],
                               [ 3, -7,  8, -5,  8,  9],
                               [ 3, -9, 12, -9,  6, 15]));
            V.dim.should.equal(3);
            V.isMaximal().should.be.true;
            V.isZero().should.be.false;
        });
        it('should construct a proper subspace', () => {
            let V = spanOf(mat([ 0, -3, -6,  4,  9],
                               [-1, -2, -1,  3,  1],
                               [-2, -3,  0,  3, -1],
                               [ 1,  4,  5, -9, -7]));
            V.dim.should.equal(3);
            V.isMaximal().should.be.false;
            V.isZero().should.be.false;
        });
        it('should detect the zero space', () => {
            Subspace.zero(3).isZero().should.be.true;
            spanOf([[.01, .01, .01], [-.01, .01, .01]])
                .isZero().should.be.false;
            spanOf([[.01, .01, .01], [-.01, .01, .01]], {ε: .02})
                .isZero().should.be.true;
        });
    });
    describe('#projectionMatrix', () => {
        let A = mat([ 0, -3, -6,  4,  9],
                    [-1, -2, -1,  3,  1],
                    [-2, -3,  0,  3, -1],
                    [ 1,  4,  5, -9, -7]);
        it('should find a nontrivial projection matrix', () => {
            let V = spanOf(A);
            let P = V.projectionMatrix;
            P.mult(P).equals(P, 1e-10).should.be.true;
            P.isSymmetric().should.be.true;
            for(let v of P.cols)
                V.contains(v).should.be.true;
            P.rank().should.equal(V.dim);
            let Vperp = V.perp();
            for(let v of P.nullBasis())
                Vperp.contains(v).should.be.true;
        });
        it('should return the zero matrix for the zero subspace', () =>
            Subspace.zero(3).projectionMatrix.equals(Matrix.zero(3))
                .should.be.true);
        it('should return the identity matrix for the full subspace', () =>
            Subspace.Rn(3).projectionMatrix.equals(Matrix.identity(3))
                .should.be.true);
    });
    describe('#add(), #intersect()', () => {
        let V = spanOf([[1,1,1], [1,0,0]]);
        let W = spanOf([[1,-1,1], [0,1,0]]);
        it('should add the subspaces', () =>
            V.add(W).equals(Subspace.Rn(3)).should.be.true);
        it('should intersect the subspaces', () =>
            V.intersect(W).equals(spanOf([[1,1,1]])).should.be.true);
        let V1 = spanOf([[0,-1,-2,1], [-3,-2,-3,4]]);
        let W1 = spanOf([[-6,-1,0,5], [4,3,3,-9]]);
        it('should add the subspaces', () =>
            V1.add(W1).equals(spanOf([[0,-1,-2,1], [-3,-2,-3,4], [4,3,3,-9]]))
                .should.be.true);
        it('should intersect the subspaces', () =>
            V1.intersect(W1).equals(spanOf([[-6,-1,0,5]])).should.be.true);
        let V2 = spanOf([[0,1,2,3]]);
        let W2 = spanOf([[2,0,-1,1]]);
        it('should form a direct sum', () =>
            V2.add(W2).equals(spanOf([[0,1,2,3],[2,0,-1,1]])).should.be.true);
        it('should form an empty intersection', () =>
            V2.intersect(W2).equals(Subspace.zero(4)).should.be.true);
        it('should add the zero space correctly', () => {
            V.add(Subspace.zero(3)).equals(V).should.be.true;
            Subspace.zero(3).add(V).equals(V).should.be.true;
        });
        it('should add Rn correctly', () => {
            V.add(Subspace.Rn(3)).isMaximal().should.be.true;
            Subspace.Rn(3).add(V).isMaximal().should.be.true;
        });
        it('should intersect the zero space correctly', () => {
            V.intersect(Subspace.zero(3)).isZero().should.be.true;
            Subspace.zero(3).intersect(V).isZero().should.be.true;
        });
    });
    describe('#orthoDecomp()', () => {
        let V = spanOf(mat([ 0, -3, -6,  4,  9],
                           [-1, -2, -1,  3,  1],
                           [-2, -3,  0,  3, -1],
                           [ 1,  4,  5, -9, -7]));
        it('should find a nontrivial decomposition', () => {
            let testVecs = [
                vec( 1, 2, 3, 4),
                vec(-2, 3, 7, 8),
                vec(11, 0, 1, 4),
            ];
            for(let v of testVecs)
                V.orthoDecomp(v).should.decompose(v, V);
        });
        it('should decompose a vector in V as itself plus zero', () => {
            let v = vec(-3, -3, -5, 5);
            let [v1, v2] = V.orthoDecomp(v);
            v1.equals(v, 1e-10).should.be.true;
            v2.isZero(1e-10).should.be.true;
        });
        it('should decompose a vector in the complement as zero plus itself', () => {
            let v = vec(0, 5, -2, 1);
            let [v1, v2] = V.orthoDecomp(v);
            v2.equals(v, 1e-10).should.be.true;
            v1.isZero(1e-10).should.be.true;
        });
    });
    describe('#contains(), #isOrthogonalTo()', () => {
        it('should contain every vector when it is maximal', () => {
            let V = spanOf(mat([ 0,  3, -6,  6,  4, -5],
                               [ 3, -7,  8, -5,  8,  9],
                               [ 3, -9, 12, -9,  6, 15]));
            V.contains(vec( 1, 2, 3)).should.be.true;
            V.contains(vec(-7, 3, 11)).should.be.true;
            Subspace.Rn(4).contains(vec(1,2,3,4)).should.be.true;
        });
        it('should contain some vectors but not others when not maximal', () => {
            let M = mat([ 0, -3, -6,  4,  9],
                        [-1, -2, -1,  3,  1],
                        [-2, -3,  0,  3, -1],
                        [ 1,  4,  5, -9, -7]);
            let vecs = Array.from(M.cols);
            let V = spanOf(vecs);
            let Vperp = V.perp();
            let v1 = Vector.linearCombination([1, 2, 3, 4, 5], vecs);
            V.contains(v1).should.be.true;
            Vperp.isOrthogonalTo(v1).should.be.true;
            v1 = Vector.linearCombination([-2, 1, -4, 3, 5], vecs)
            V.contains(v1).should.be.true;
            Vperp.isOrthogonalTo(v1).should.be.true;
            V.contains(vec(0,0,0,1)).should.be.false;
            Vperp.isOrthogonalTo(vec(0,0,0,1)).should.be.false;
            Subspace.zero(4).contains(vec(1,2,3,4)).should.be.false;
            Subspace.zero(4).isOrthogonalTo(vec(1,2,3,4)).should.be.true;
            Subspace.zero(4).contains(Vector.zero(4)).should.be.true;
            Subspace.zero(4).isOrthogonalTo(Vector.zero(4)).should.be.true;
        });
        it('should throw if the vector has the wrong length', () => {
            let V = Subspace.Rn(3);
            V.contains.bind(V, vec(1,2)).should.throw(/number of entries/);
        });
    });
    describe('#isSubspaceOf()', () => {
        it('should detect inclusion with different generating sets', () => {
            let V = spanOf(mat([ 0, -3, -6,  4,  9],
                               [-1, -2, -1,  3,  1],
                               [-2, -3,  0,  3, -1],
                               [ 1,  4,  5, -9, -7]).transpose);
            let W = spanOf(mat([1, 0, -3, 0,  5],
                               [0, 0,  0, 1,  0]).transpose);
            W.isSubspaceOf(V).should.be.true;
            V.isSubspaceOf(W).should.be.false;
        });
    });
    describe('#equals()', () => {
        it('should detect equality with different bases', () => {
            let V = spanOf([[ 0, -3, -6,  4,  9],
                            [-1, -2, -1,  3,  1],
                            [-2, -3,  0,  3, -1],
                            [ 1,  4,  5, -9, -7]]);
            let W = spanOf([[-1, -5, -7,  7, 10],
                            [-3, -5, -1,  6,  0],
                            [-1,  1,  5, -6, -8]]);
            V.equals(W).should.be.true;
        });
    });
    describe('#perp()', () => {
        it('should compute the orthogonal complement', () => {
            let V = spanOf(mat([ 0, -3, -6,  4,  9],
                               [-1, -2, -1,  3,  1],
                               [-2, -3,  0,  3, -1],
                               [ 1,  4,  5, -9, -7]));
            let W = V.perp();
            (V.dim + W.dim).should.equal(4);
            expect(V.basis).not.to.be.null;
            for(let v of V.basis!.cols)
                W.isOrthogonalTo(v).should.be.true;
        });
        it('should compute the orthogonal complement of {0}', () =>
           Subspace.zero(3).perp().isMaximal().should.be.true);
        it('should compute the orthogonal complement of R^n', () =>
           Subspace.Rn(3).perp().isZero().should.be.true);
    });
});

