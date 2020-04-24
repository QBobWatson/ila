'use strict'; // -*- js2 -*-

import { Subspace, Vector, Matrix } from "../lib/linalg.js";

import should from 'should';

const vec = Vector.create;
const mat = Matrix.create;
const subspc = (...args) => new Subspace(...args);


should.use(function(should, Assertion) {
    Assertion.add('decompose', function(vec, V, ε=1e-10) {
        this.params = {
            operator: 'to orthogonally decompose',
            expected: vec
        };
        let [v1, v2] = this.obj;
        v1.clone().add(v2).equals(vec, 1e-10).should.be.true();
        V.contains(v1).should.be.true();
        for(let v of V.basis.cols())
            v.dot(v2).should.be.approximately(0, 1e-10);
    });
});


describe('Subspace', () => {
    describe('#constructor()', () => {
        it('should construct a subspace from an Array of vectors', () =>
           subspc([[1,2,3], [4,5,6]]).basis.equals(mat([1,4],[2,5],[3,6]))
                   .should.be.true());
        it('should construct a subspace from a Matrix', () => {
            let M = mat([1,4],[2,5],[3,6]);
            new Subspace(M).basis.equals(M).should.be.true();
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
        let V = subspc(A);
        it('should have 4 decimal places by default', () =>
           V.toString().should.eql("Subspace of R^4 of dimension 3 with basis\n"
                                   + "[ 0.0000] [-3.0000] [ 4.0000]\n"
                                   + "[-1.0000] [-2.0000] [ 3.0000]\n"
                                   + "[-2.0000] [-3.0000] [ 3.0000]\n"
                                   + "[ 1.0000] [ 4.0000] [-9.0000]"));
        it('can have other precision', () =>
           V.toString(2).should.eql("Subspace of R^4 of dimension 3 with basis\n"
                                    + "[ 0.00] [-3.00] [ 4.00]\n"
                                    + "[-1.00] [-2.00] [ 3.00]\n"
                                    + "[-2.00] [-3.00] [ 3.00]\n"
                                    + "[ 1.00] [ 4.00] [-9.00]"));
    });
    describe('#basis, #dim, #isMaximal(), #isZero()', () => {
        it('should find a full basis', () => {
            let V = subspc(mat([ 0,  3, -6,  6,  4, -5],
                               [ 3, -7,  8, -5,  8,  9],
                               [ 3, -9, 12, -9,  6, 15]));
            V.basis.equals(mat([0,  3, 4],
                               [3, -7, 8],
                               [3, -9, 6])).should.be.true();
            V.dim.should.equal(3);
            V.isMaximal().should.be.true();
        });
        it('should find a proper basis', () => {
            let V = subspc(mat([ 0, -3, -6,  4,  9],
                               [-1, -2, -1,  3,  1],
                               [-2, -3,  0,  3, -1],
                               [ 1,  4,  5, -9, -7]));
            V.basis.equals(mat([ 0, -3,  4],
                               [-1, -2,  3],
                               [-2, -3,  3],
                               [ 1,  4, -9])).should.be.true();
            V.dim.should.equal(3);
            V.isMaximal().should.be.false();
        });
        it('should detect the zero space', () => {
            Subspace.zero(3).isZero().should.be.true();
            subspc([[.01, .01, .01], [-.01, .01, .01]])
                .isZero().should.be.false();
            subspc([[.01, .01, .01], [-.01, .01, .01]], {ε: .02})
                .isZero().should.be.true();
        });
    });
    describe('#ONbasis()', () => {
        it('should find an orthonormal basis', () => {
            let V = subspc(mat([ 0, -3, -6,  4,  9],
                               [-1, -2, -1,  3,  1],
                               [-2, -3,  0,  3, -1],
                               [ 1,  4,  5, -9, -7]));
            let Q = V.ONbasis();
            Q.colSpace().equals(V).should.be.true();
            Q.transpose.mult(Q).equals(Matrix.identity(V.dim), 1e-10)
                .should.be.true();
        });
    });
    describe('#add(), #intersect()', () => {
        let V = subspc([[1,1,1], [1,0,0]]);
        let W = subspc([[1,-1,1], [0,1,0]]);
        it('should add the subspaces', () =>
           V.add(W).equals(Subspace.Rn(3)).should.be.true());
        it('should intersect the subspaces', () =>
           V.intersect(W).equals(subspc([[1,1,1]])).should.be.true());
        let V1 = subspc([[0,-1,-2,1], [-3,-2,-3,4]]);
        let W1 = subspc([[-6,-1,0,5], [4,3,3,-9]]);
        it('should add the subspaces', () =>
           V1.add(W1).equals(subspc([[0,-1,-2,1], [-3,-2,-3,4], [4,3,3,-9]]))
                   .should.be.true());
        it('should intersect the subspaces', () =>
           V1.intersect(W1).equals(subspc([[-6,-1,0,5]])).should.be.true());
    });
    describe('#projectionMatrix()', () => {
        let A = mat([ 0, -3, -6,  4,  9],
                    [-1, -2, -1,  3,  1],
                    [-2, -3,  0,  3, -1],
                    [ 1,  4,  5, -9, -7]);
        it('should find a nontrivial projection matrix using A(A^TA)^(-1)A^T', () => {
            let V = subspc(A);
            let P = V.projectionMatrix();
            P.mult(P).equals(P, 1e-10).should.be.true();
            P.colSpace().equals(V).should.be.true();
            P.nullSpace().equals(V.perp()).should.be.true();
        });
        it('should find a nontrivial projection matrix using QQ^T', () => {
            let V = subspc(A);
            let Q = V.ONbasis();
            let P = V.projectionMatrix();
            P.mult(P).equals(P, 1e-10).should.be.true();
            P.colSpace().equals(V).should.be.true();
            P.nullSpace().equals(V.perp()).should.be.true();
        });
        it('should return the zero matrix for the zero subspace', () =>
           Subspace.zero(3).projectionMatrix().equals(Matrix.zero(3))
                   .should.be.true());
        it('should return the identity matrix for the full subspace', () =>
           Subspace.Rn(3).projectionMatrix().equals(Matrix.identity(3))
                   .should.be.true());
    });
    describe('#orthoDecomp()', () => {
        let V = subspc(mat([ 0, -3, -6,  4,  9],
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
            v1.equals(v, 1e-10).should.be.true();
            v2.isZero(1e-10).should.be.true();
        });
        it('should decompose a vector in the complement as zero plus itself', () => {
            let v = vec(0, 5, -2, 1);
            let [v1, v2] = V.orthoDecomp(v);
            v2.equals(v, 1e-10).should.be.true();
            v1.isZero(1e-10).should.be.true();
        });
        // Now do everything using the projection matrix
        let P = V.projectionMatrix();
        it('should find a nontrivial decomposition with P', () => {
            let testVecs = [
                vec( 1, 2, 3, 4),
                vec(-2, 3, 7, 8),
                vec(11, 0, 1, 4),
            ];
            for(let v of testVecs)
                V.orthoDecomp(v).should.decompose(v, V);
        });
        it('should decompose a vector in V as itself plus zero with P', () => {
            let v = vec(-3, -3, -5, 5);
            let [v1, v2] = V.orthoDecomp(v);
            v1.equals(v, 1e-10).should.be.true();
            v2.isZero(1e-10).should.be.true();
        });
        it('should decompose a vector in the complement as zero plus itself with P',
           () => {
            let v = vec(0, 5, -2, 1);
            let [v1, v2] = V.orthoDecomp(v);
            v2.equals(v, 1e-10).should.be.true();
            v1.isZero(1e-10).should.be.true();
        });
    });
    describe('#contains()', () => {
        it('should contain every vector when it is maximal', () => {
            let V = subspc(mat([ 0,  3, -6,  6,  4, -5],
                               [ 3, -7,  8, -5,  8,  9],
                               [ 3, -9, 12, -9,  6, 15]));
            V.contains(vec( 1, 2, 3)).should.be.true();
            V.contains(vec(-7, 3, 11)).should.be.true();
            Subspace.Rn(4).contains(vec(1,2,3,4)).should.be.true();
        });
        it('should contain some vectors but not others when not maximal', () => {
            let M = mat([ 0, -3, -6,  4,  9],
                        [-1, -2, -1,  3,  1],
                        [-2, -3,  0,  3, -1],
                        [ 1,  4,  5, -9, -7]);
            let vecs = Array.from(M.cols());
            let V = subspc(vecs);
            V.contains(Vector.linearCombination([1, 2, 3, 4, 5], vecs))
                .should.be.true();
            V.contains(Vector.linearCombination([-2, 1, -4, 3, 5], vecs))
                .should.be.true();
            V.contains(vec(0,0,0,1)).should.be.false();
            Subspace.zero(4).contains(vec(1,2,3,4)).should.be.false();
            Subspace.zero(4).contains(Vector.zero(4)).should.be.true();
        });
        it('should throw if the vector has the wrong length', () => {
            let V = Subspace.Rn(3);
            V.contains.bind(V, vec(1,2)).should.throw(/number of entries/);
        });
    });
    describe('#isSubspaceOf()', () => {
        it('should detect inclusion with different generating sets', () => {
            let V = subspc(mat([ 0, -3, -6,  4,  9],
                               [-1, -2, -1,  3,  1],
                               [-2, -3,  0,  3, -1],
                               [ 1,  4,  5, -9, -7]).transpose);
            let W = subspc(mat([1, 0, -3, 0,  5],
                               [0, 0,  0, 1,  0]).transpose);
            W.isSubspaceOf(V).should.be.true();
            V.isSubspaceOf(W).should.be.false();
        });
    });
    describe('#equals()', () => {
        it('should detect equality with different bases', () => {
            let V = subspc([[ 0, -3, -6,  4,  9],
                            [-1, -2, -1,  3,  1],
                            [-2, -3,  0,  3, -1],
                            [ 1,  4,  5, -9, -7]]);
            let W = subspc([[-1, -5, -7,  7, 10],
                            [-3, -5, -1,  6,  0],
                            [-1,  1,  5, -6, -8]]);
            V.equals(W).should.be.true();
        });
    });
    describe('#perp()', () => {
        it('should compute the orthogonal complement', () => {
            let V = subspc(mat([ 0, -3, -6,  4,  9],
                               [-1, -2, -1,  3,  1],
                               [-2, -3,  0,  3, -1],
                               [ 1,  4,  5, -9, -7]));
            let W = V.perp();
            (V.dim + W.dim).should.equal(4);
            for(let v of V.basis.cols())
                W.isOrthogonalTo(v).should.be.true();
        });
        it('should compute the orthogonal complement of {0}', () =>
           Subspace.zero(3).perp().isMaximal().should.be.true());
        it('should compute the orthogonal complement of R^n', () =>
           Subspace.Rn(3).perp().isZero().should.be.true());
    });
});
