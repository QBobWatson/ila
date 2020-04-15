'use strict'; // -*- js2 -*-

import { cardano, quadratic } from '../lib/linalg.js';
import Complex from "../lib/complex.js";
import Vector from "../lib/vector.js";
import Matrix from "../lib/matrix.js";
import Subspace from "../lib/subspace.js";

import should from 'should';

const vec = Vector.create;
const mat = (...args) => new Matrix(...args);


should.use(function(should, Assertion) {
    // Custom assertion for comparing arrays of numbers using 'approximately'
    function approximate(obj, other, ε) {
        if(typeof obj === "number")
            obj.should.be.approximately(other, ε);
        else if(obj instanceof Matrix) {
            other.should.be.an.instanceOf(Matrix);
            obj.m.should.equal(other.m);
            obj.n.should.equal(other.n);
            approximate([...obj.rows()], [...other.rows()], ε);
        }
        else {
            obj.should.be.an.Array();
            other.should.be.an.Array();
            other.should.have.length(obj.length);
            for(let i = 0; i < obj.length; ++i)
                approximate(obj[i], other[i], ε);
        }
    }

    Assertion.add('resemble', function(other, ε=1e-10) {
        this.params = {
            operator: 'to approximately equal',
            expected: other
        };
        approximate(this.obj, other, ε);
    });

    Assertion.add('decompose', function(vec, V, ε=1e-10) {
        this.params = {
            operator: 'to orthogonally decompose',
            expected: vec
        };
        let [v1, v2] = this.obj;
        v1.clone().add(v2).equals(vec, 1e-10).should.be.true();
        V.contains(v1).should.be.true();
        for(let v of V.generators.cols())
            v.dot(v2).should.be.approximately(0, 1e-10);
    });

    Assertion.add('factorize', function(A, which, ε=1e-10) {
        this.params = {
            operator: 'factorizes',
            expected: A
        };
        if(which === 'PA=LU') {
            let {P, L, U, E} = this.obj;
            L.isLowerUnip().should.be.true('L is not lower unipotent');
            U.isEchelon().should.be.true('U is not in echelon form');
            for(let [i, j] of U.leadingEntries())
                Math.abs(U[i][j]).should.be.above(ε, 'U has small pivots');
            let PA = Matrix.permutation(P).mult(A);
            let LU = L.mult(U);
            PA.equals(LU, 1e-10).should.be.true(
                `Matrix is not correctly factored: PA =\n${PA.toString(2)}\nLU =\n${LU.toString(2)}`);
            let EA = E.mult(A);
            EA.equals(U, 1e-10).should.be.true(
                `Matrix is not correctly factored: EA =\n${EA.toString(2)}\nU =\n${U.toString(2)}`);
        } else if(which === 'QR') {
            let {Q, R, LD} = this.obj;
            R.isUpperTri().should.be.true();
            let QTQ = Q.transpose.mult(Q);
            QTQ.isDiagonal(1e-10).should.be.true();
            for(let j = 0; j < QTQ.n; ++j) {
                if(LD.includes(j)) {
                    QTQ[j][j].should.be.approximately(0, 1e-10);
                    R[j][j].should.be.approximately(0, 1e-10);
                } else
                    QTQ[j][j].should.be.approximately(1, 1e-10);
            }
            let QR = Q.mult(R);
            QR.equals(A, 1e-10).should.be.true(
                `Matrix is not correctly factored: QR =\n${QR.toString(2)}\nA =\n${A.toString(2)}`);
        }
    });

    Assertion.add('diagonalize', function(A, block=false, ε=1e-10) {
        this.params = {
            operator: 'diagonalizes',
            expected: A
        };
        let {C, D} = this.obj;
        if(!block) D.isDiagonal().should.be.true();
        C.isInvertible().should.be.true();
        C.mult(D).equals(A.mult(C), ε).should.be.true();
    });

    Assertion.add('cplxEigenvectors', function(M, λ, ε=1e-10) {
        this.params = {
            operator: 'be complex eigenvectors'
        };
        for(let [v_Re, v_Im] of this.obj) {
            let Mv_Re = M.apply(v_Re);
            let Mv_Im = M.apply(v_Im);
            let λv_Re = v_Re.clone().scale(λ.Re).sub(v_Im.clone().scale(λ.Im));
            let λv_Im = v_Re.clone().scale(λ.Im).add(v_Im.clone().scale(λ.Re));
            [Mv_Re, Mv_Im].should.resemble([λv_Re, λv_Im], ε);
        }
    });

    Assertion.add('solve', function(M, b, leastSquares=false, ε=1e-10) {
        this.params = {
            operator: 'solves',
        };
        let x = this.obj;
        let Mx = M.apply(x);
        if(!leastSquares)
            Mx.equals(b, ε).should.be.true(
                `x = ${x.toString(2)}\nMx = ${Mx.toString(2)}\nb = ${b.toString(2)}`);
        else {
            let bb = M.colSpace().project(b);
            Mx.equals(bb, ε).should.be.true();
        }
    });
});


describe('Vector', () => {
    describe('#constant()', () => {
        it('should create a constant vector', () =>
           Vector.constant(3, 1).should.eql(vec(1, 1, 1)));
        it('should create vectors of length one correctly', () =>
           Vector.constant(1, 3).should.resemble([3]));
    });
    describe('#zero()', () =>
             it('should create the zero vector', () =>
                Vector.zero(3).should.eql(vec(0, 0, 0))));
    describe('#e()', () => {
        it('should create a unit coordinate vector', () =>
           Vector.e(1, 3).should.eql(vec(0, 1, 0)));
        it('should create a scaled coordinate vector', () =>
           Vector.e(1, 3, 2).should.eql(vec(0, 2, 0)));
    });
    describe('#isLinearlyIndependent()', () => {
        it('should say that m > n vectors are linearly dependent', () =>
           Vector.isLinearlyIndependent(
               [vec(1,2), vec(3,4), vec(5,6)])
                   .should.be.false());
        it('should detect linear dependence', () =>
           Vector.isLinearlyIndependent(
               [vec(1,2,3), vec(4,5,6), vec(3,3,3)])
                   .should.be.false());
        it('should detect linear independence', () =>
           Vector.isLinearlyIndependent(
               [vec(1,2,3), vec(4,5,6), vec(3,3,4)])
                   .should.be.true());
    });
    describe('#linearlyIndependentSubset()', () => {
        it('should return a proper linearly independent subset', () =>
           Vector.linearlyIndependentSubset(
               [vec(2, -4, 3, -6, 8),
                vec(-6, 5, 5, 4, -3),
                vec(6, -7, -1, -8, 9)]).should.eql(
                    [vec(2, -4, 3, -6, 8),
                     vec(-6, 5, 5, 4, -3)]));
        it('should return a full linearly independent subset', () => {
            let vecs = [vec(10,-3,  5),
                        vec(-7, 2, -1),
                        vec( 0, 6,  5)];
            Vector.linearlyIndependentSubset(vecs).should.eql(vecs);
        });
    });
    describe('#linearCombination()', () =>
             it('should compute a linear combination', () =>
                Vector.linearCombination(
                    [2, -1, 3],
                    [vec(1, 3), vec(2, 4), vec(1, 1)])
                        .should.eql(vec(3, 5))));
    describe('#equals()', () =>  {
        let v = vec(1,2,3), w = vec(1,2,3), u = vec(1,2);
        let v1 = vec(1.01, 1.99, 3);
        it('should compare two vectors as equal', () =>
           v.equals(w).should.be.true());
        it('should compare different vectors as not equal', () =>
           v.equals(v1).should.be.false());
        it('should compare similar vectors as equal with ε>0', () =>
           v.equals(v1, 0.05).should.be.true());
        it('should compare similar vectors as not equal with ε small', () =>
           v.equals(v1, 0.005).should.be.false());
        it('should compare vectors of different lengths as not equal', () =>
           u.equals(v).should.be.false());
    });
    describe('#clone()', () => {
        let v = vec(1,2,3), w = v.clone();
        it('should create distinct objects', () =>
           v.should.not.equal(w));
        it('should create equal vectors', () =>
           v.should.eql(w));
    });
    describe('#toString()', () => {
        let v = vec(1,2,3);
        it('should have 4 decimal places by default', () =>
           v.toString().should.eql("[1.0000 2.0000 3.0000]"));
        it('can have other precision', () =>
           v.toString(2).should.eql("[1.00 2.00 3.00]"));
    });
    describe('#isZero()', () => {
        it('should detect the zero vector', () =>
           Vector.zero(3).isZero().should.be.true());
    });
    describe('#isZero()', () => {
        it('should detect nonzero vectors', () => {
            Vector.constant(3, 0.01).isZero().should.be.false();
            Vector.constant(3, 0.01).isZero(.1).should.be.true();
        });
    });
    describe('#sizesq()', () =>
             it('should have the correct squared size', () =>
                vec(1, 2, 3).sizesq.should.equal(1 + 4 + 9)));
    describe('#size()', () =>
             it('should have the correct size', () =>
                vec(3, 4).size.should.be.approximately(5, 1e-10)));
    describe('#normalize()', () => {
        it('should give the unit vector in the same direction', () =>
           vec(3, 4).normalize().should.resemble(vec(3/5, 4/5)));
        it('should throw when normalizing the zero vector', () => {
            let v = Vector.zero(2);
            v.normalize.bind(v).should.throw(/zero vector/);
        });
    });
    describe('#add() and #sub()', () => {
        it('should add componentwise', () =>
           vec(3, 4).add(vec(1, 2)).should.eql(vec(4, 6)));
        it('should add with a factor', () =>
           vec(3, 4).add(vec(1, 2), -1).should.eql(vec(2, 2)));
        it('should subtract componentwise', () =>
           vec(3, 4).sub(vec(1, 2)).should.eql(vec(2, 2)));
        it('should throw when vectors have different lengths', () => {
            let v = vec(3, 4);
            v.add.bind(v, vec(1, 2, 3)).should.throw(/different lengths/);
        });
    });
    describe('#scale()', () =>
             it('should scale componentwise', () =>
                vec(3, 4).scale(2).should.eql(vec(6, 8))));
    describe('#dot()', () => {
        it('should compute the dot product', () =>
           vec(3, 4).dot(vec(1, 2)).should.equal(11));
        it('should throw when vectors have different lengths', () => {
            let v = vec(3, 4);
            v.dot.bind(v, vec(1, 2, 3)).should.throw(/different lengths/);
        });
    });
});


describe('Complex', () => {
    const π = Math.PI;

    describe('#fromPolar()', () => {
        it('should create positive real numbers with θ=0', () =>
           Complex.fromPolar(1, 0).should.resemble(new Complex(1, 0)));
        it('should create negative real numbers with θ=π', () =>
           Complex.fromPolar(1, π).should.resemble(new Complex(-1, 0)));
        it('should create purely imaginary numbers with θ=π/2', () =>
           Complex.fromPolar(1, π/2).should.resemble(new Complex(0, 1)));
        it('should create a cube root of -8 with r=2 and θ=π/3', () => {
            let z = Complex.fromPolar(2, π/3), w = z.clone();
            z.should.resemble(new Complex(1, Math.sqrt(3)));
            z.mult(w).mult(w).should.resemble(new Complex(-8, 0));
        });
    });
    describe('#i()', () => {
        it('should be 0 + 1i', () =>
           Complex.i.should.eql(new Complex(0, 1)));
        it('should be a square root of -1', () =>
           Complex.i.mult(Complex.i).should.eql(new Complex(-1,0)));
    });
    describe('#equals()', () => {
        let z = new Complex(3, 0), w = new Complex(3.01, 0.01);
        it('should compare equal to a real number if Im==0', () =>
           z.equals(3).should.be.true());
        it('should compare not equal to a real number if Im!=0', () =>
           w.equals(3).should.be.false());
        it('should compare equal when ε>0', () =>
           w.equals(3, .05).should.be.true());
    });
    describe('#toString()', () => {
        let z = new Complex(3, 4);
        it('should have 4 decimal places by default', () =>
           z.toString().should.eql("3.0000 + 4.0000 i"));
        it('can have other precision', () =>
           z.toString(2).should.eql("3.00 + 4.00 i"));
        it('subtracts a negative imaginary part', () =>
           z.conj().toString(2).should.eql("3.00 - 4.00 i"));
    });
    describe('#Re(), #im(), #mod()', () => {
        let z = new Complex(3, 4);
        it('should have the correct real part', () =>
           z.Re.should.equal(3));
        it('should have the correct imaginary part', () =>
           z.Im.should.equal(4));
        it('should have the correct modulus', () =>
           z.mod.should.be.approximately(5, 1e-10));
        it('should set the real part', () => {
            z.Re = 1;
            z.should.eql(new Complex(1, 4));
        });
        it('should set the imaginary part', () => {
            z.Im = 1;
            z.should.eql(new Complex(1, 1));
        });
    });
    describe('#arg()', () => {
        it('should compute the modulus of a 6th root of unity', () =>
           new Complex( 1/2,  Math.sqrt(3)/2).arg
                   .should.be.approximately(   π/3, 1e-10));
        it('should compute the modulus of a 3rd root of unity', () =>
           new Complex(-1/2,  Math.sqrt(3)/2).arg
                   .should.be.approximately( 2*π/3, 1e-10));
        it('should compute the modulus of another 6th root of unity', () =>
           new Complex( 1/2, -Math.sqrt(3)/2).arg
                   .should.be.approximately(  -π/3, 1e-10));
        it('should compute the modulus of another 3rd root of unity', () =>
           new Complex(-1/2, -Math.sqrt(3)/2).arg
                   .should.be.approximately(-2*π/3, 1e-10));
    });
    describe('#conj()', () => {
        let z = new Complex(1, 2);
        it('should negate the imaginary part', () =>
           z.conj().should.eql(new Complex(1, -2)));
    });
    describe('#add(), #sub()', () => {
        it('should add componentwise', () =>
           new Complex(3, 4).add(new Complex(1, 2)).should.eql(new Complex(4, 6)));
        it('should add componentwise with factor', () =>
           new Complex(3, 4).add(new Complex(1, 2), -1).should.eql(new Complex(2, 2)));
        it('should subtract componentwise', () =>
           new Complex(3, 4).sub(new Complex(1, 2)).should.eql(new Complex(2, 2)));
        it('should add real numbers correctly', () =>
           new Complex(3, 4).add(2).should.eql(new Complex(5, 4)));
        it('should add real numbers with factor', () =>
           new Complex(3, 4).add(2, -1).should.eql(new Complex(1, 4)));
        it('should subtract real numbers correctly', () =>
           new Complex(3, 4).sub(2).should.eql(new Complex(1, 4)));
    });
    describe('#mult()', () => {
        it('should multiply to -5 + 10 i', () =>
           new Complex(1, 2).mult(new Complex(3, 4))
                   .should.eql(new Complex(-5, 10)));
        it('should scalar-multiply to 2 + 4 i', () =>
           new Complex(1, 2).mult(2).should.eql(new Complex(2, 4)));
    });
    describe('#recip()', () => {
        it('should compute the reciprocal', () =>
           new Complex(3, 4).recip()
                   .should.eql(new Complex(3/25, -4/25)));
        it('should throw when the number is zero', () => {
            let z = new Complex(0);
            z.recip.bind(z).should.throw(/divide by zero/);
        });
    });
    describe('#div()', () => {
        it('should compute the quotient', () =>
           new Complex(1, 2).div(new Complex(3, 4))
                   .should.resemble(new Complex(11/25, 2/25)));
        it('should compute the quotient by a scalar', () =>
           new Complex(1, 2).div(3).should.eql(new Complex(1/3, 2/3)));
        it('should throw when the denominator is zero', () => {
            let z = new Complex(1, 2);
            z.div.bind(z, new Complex(0)).should.throw(/divide by zero/);
        });
        it('should throw when the denominator is the scalar zero', () => {
            let z = new Complex(1, 2);
            z.div.bind(z, 0).should.throw(/divide by zero/);
        });
    });
});


describe('Root finders', function() {
    const ζ = new Complex(-1/2, Math.sqrt(3)/2);

    describe('quadratic()', function() {
        it('should find a double root', function() {
            // x^2 - 2x + 1 = (x-1)^2
            quadratic(-2,  1).should.resemble([[1, 2]]);
        });
        it('should find simple real roots', function() {
            // x^2 - 1 = (x-1) (x+1)
            quadratic( 0, -1).should.resemble([[-1, 1], [1, 1]]);
            // x^2 - 3x + 2 = (x-1) (x-2)
            quadratic(-3,  2).should.resemble([[1, 1], [2, 1]]);
        });
        it('should find complex roots', function() {
            quadratic( 0,  1).should.resemble(
                [[Complex.i, 1], [Complex.i.conj(), 1]]);
            quadratic( 1,  1).should.resemble(
                [[ζ, 1], [ζ.clone().conj(), 1]]);
        });
    });

    describe('cardano()', () => {
        it('should find a triple root', () => {
            // x^3 - 6x^2 + 12x - 8 = (x-2)^3
            cardano(-6, 12, -8).should.resemble([[2, 3]]);
        });
        it('should find double roots', () => {
            // x^3 - 4x^2 + 5x - 2 = (x-2) (x-1)^2
            cardano(-4,  5, -2).should.resemble([[1, 2], [2, 1]]);
            // x^3 + 4x^2 + 5x + 2 = (x+2) (x+1)^2
            cardano( 4,  5,  2).should.resemble([[-2, 1], [-1, 2]]);
        });
        it('should find simple real roots', () => {
            // x^3 + 2x^2 - x - 2 = (x+2) (x+1) (x-1)
            cardano( 2, -1, -2).should.resemble([[-2, 1], [-1, 1], [1, 1]]);
            // x^3 +  x^2 - 4x - 4 = (x+2) (x+1) (x-2)
            cardano( 1, -4, -4).should.resemble([[-2, 1], [-1, 1], [2, 1]]);
        });
        it('should find complex Roots', () => {
            // x^2 - 2x^2 + x - 2 = (x-2) (x^2+1)
            cardano(-2,  1, -2).should.resemble(
                [[2, 1], [Complex.i.conj(), 1], [Complex.i, 1]]);
            // x^2 - x^2 - x - 2 = (x-2) (x^2+x+1)
            cardano(-1, -1, -2).should.resemble(
                [[2, 1], [ζ.clone().conj(), 1], [ζ, 1]]);
        });

    });
});


describe('Matrix', () => {
    describe('#identity()', () => {
        it('should return a 3x3 identity matrix', () =>
           Matrix.identity(3).equals(
               mat([1,0,0],[0,1,0],[0,0,1])).should.be.true());
        it('should return a scaled 3x3 identity matrix', () =>
           Matrix.identity(3, 2).equals(
               mat([2,0,0],[0,2,0],[0,0,2])).should.be.true());
    });
    describe('#zero()', () => {
        it('should return the 2x3 zero matrix', () =>
           Matrix.zero(2, 3).equals(
               mat([0, 0, 0], [0, 0, 0])).should.be.true());
        it('should return the 2x2 zero matrix', () =>
           Matrix.zero(2).equals(
               mat([0, 0], [0, 0])).should.be.true());
    });
    describe('#permutation()', () =>
             it('should create a 3x3 permutation matrix', () =>
                Matrix.permutation([2, 0, 1]).equals(
                    mat([0,0,1],[1,0,0],[0,1,0])).should.be.true()));
    describe('#constructor()', () =>
             it('should create a 3x2 matrix', () =>
                mat([1, 2], [3, 4], [5, 6]).should
                        .have.properties({m: 3, n: 2})));
    describe('#equals()', () => {
        let M = mat([0, 1, 2], [3, 4, 5]), M1 = M.clone();
        let N = mat([0.01, 1.01, 2.01], [3, 4, 5]);
        let M3 = mat([0, 1, 2]);
        it('should compare two matrices as equal', () =>
           M.equals(M1).should.be.true());
        it('should compare different matrices as not equal', () =>
           M.equals(N).should.be.false());
        it('should compare similar matrices as equal with ε>0', () =>
           M.equals(N, 0.02).should.be.true());
        it('should compare matrices of different sizes as not equal', () =>
           M3.equals(M).should.be.false());
    });
    describe('#clone()', () => {
        let M = mat([1,2],[3,4]), N = M.clone();
        it('should create distinct objects', () =>
           M.should.not.equal(N));
        it('should create identical objects', () =>
           M.equals(N).should.be.true());
    });
    describe('#toString()', () => {
        let M = mat([10,2],[3,4]);
        it('should have 4 decimal places by default', () =>
           M.toString().should.eql("10.0000 2.0000\n 3.0000 4.0000"));
        it('can have other precision', () =>
           M.toString(2).should.eql("10.00 2.00\n 3.00 4.00"));
    });
    describe('#row() and #col()', () => {
        let M = mat([1,2],[3,4]);
        it('should return the second row', () =>
           M.row(1).should.eql(vec(3, 4)));
        it('should return the second column', () =>
           M.col(1).should.eql(vec(2, 4)));
    });
    describe('#rows and #cols', () => {
        let M = mat([1,2],[3,4]);
        it('should iterate over the rows', () =>
           [...M.rows()].should.eql([vec(1,2),vec(3,4)]));
        it('should iterate over the columns', () =>
           [...M.cols()].should.eql([vec(1,3),vec(2,4)]));
    });
    describe('#add() and #sub()', () => {
        it('should add componentwise', () =>
           mat([1,2],[3,4]).add(mat([2,1],[4,3])).equals(
               mat([3,3],[7,7])).should.be.true());
        it('should add with a factor', () =>
           mat([1,2],[3,4]).add(mat([2,1],[4,3]), 2).equals(
               mat([5,4],[11,10])).should.be.true());
        it('should subtract componentwise', () =>
           mat([1,2],[3,4]).sub(mat([2,1],[4,3])).equals(
               mat([-1,1],[-1,1])).should.be.true());
        it('should throw when matrices have different sizes', () => {
            let M = mat([1,2,3],[4,5,6]);
            M.add.bind(M, mat([1,2],[3,4])).should.throw(/different sizes/);
        });
    });
    describe('#scale()', () => {
        it('should scale entries', () =>
           mat([1,2],[3,4]).scale(2).equals(
               mat([2,4],[6,8])).should.be.true());
        it('should add with a factor', () =>
           mat([1,2],[3,4]).add(mat([2,1],[4,3]), 2).equals(
               mat([5,4],[11,10])).should.be.true());
        it('should subtract componentwise', () =>
           mat([1,2],[3,4]).sub(mat([2,1],[4,3])).equals(
               mat([-1,1],[-1,1])).should.be.true());
        it('should throw when matrices have different sizes', () => {
            let M = mat([1,2,3],[4,5,6]);
            M.add.bind(M, mat([1,2],[3,4])).should.throw(/different sizes/);
        });
    });
    describe('#transpose', () =>
             it('should construct the transpose', () =>
                mat([1,2,3], [4,5,6]).transpose.equals(
                    mat([1,4], [2,5], [3,6])).should.be.true()));
    describe('#mult()', () => {
        it('should compute the product', () =>
           mat([1,2],[3,4],[5,6]).mult(mat([1,2,3],[4,5,6]))
                   .equals(mat([9,12,15],[19,26,33],[29,40,51]))
                   .should.be.true());
        it('should throw when the matrices have incompatible dimensions', () => {
            let M = mat([1,2],[3,4],[5,6]);
            M.mult.bind(M, M).should.throw(/incompatible dimensions/);
        });
    });
    describe('#apply()', () => {
        it('should compute the product', () =>
           mat([1,2,3],[4,5,6]).apply(vec(7,8,9))
                   .equals(vec(50,122)).should.be.true());
        it('should throw when the vector has incompatible length', () => {
            let M = mat([1,2],[3,4],[5,6]);
            M.mult.bind(M, vec(1,2,3))
                .should.throw(/incompatible dimensions/);
        });
    });
    describe('#isZero()', () => {
        it('should detect the zero matrix', () => {
            let M = Matrix.identity(3, 0.01);
            M.isZero().should.be.false();
            M.isZero(0.1).should.be.true();
        });
    });
    describe('#isUpperTri(), #isLowerTri(), #isDiagonal()', () => {
        let M = mat([1, 1, 1], [0, 1, 1], [0, 0, 1]);
        let N = mat([0, 0, 0], [0, 0, 0], [0, 1, 0]);
        it('should be upper-triangular', () =>
           M.isUpperTri().should.be.true());
        it('should not be upper-triangular', () =>
           N.isUpperTri().should.be.false());
        it('should be lower-triangular', () =>
           M.transpose.isLowerTri().should.be.true());
        it('should not be lower-triangular', () =>
           N.transpose.isLowerTri().should.be.false());
        it('should not be diagonal', () => {
            M.isDiagonal().should.be.false();
            N.isDiagonal().should.be.false();
            Matrix.identity(3).isDiagonal().should.be.true();
        });
    });
    describe('#isUpperUnip() and #isLowerUnip()', () => {
        let M = mat([1, 1, 1], [0, 1, 1], [0, 0, 1]);
        let N = mat([1, 1, 1], [0, 2, 1], [0, 0, 1]);
        it('should be upper-uniponent', () =>
           M.isUpperUnip().should.be.true());
        it('should not be upper-unipotent', () =>
           N.isUpperUnip().should.be.false());
        it('should be lower-unipotent', () =>
           M.transpose.isLowerUnip().should.be.true());
        it('should not be lower-unipotent', () =>
           N.transpose.isLowerUnip().should.be.false());
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
    describe('#isEchelon() and #isRREF()', () => {
        let testMats1 = [
            mat([1,  0,  2],
                [0,  1, -1]),
            mat([0,  1,  8, 0]),
            mat([1, 17,  0],
                [0,  0,  1]),
            Matrix.zero(2, 3)
        ];
        it('should be in EF and RREF', () => {
            for(let M of testMats1) {
                M.isEchelon().should.be.true();
                M.isRREF().should.be.true();
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
        it('should be in EF but not RREF', () => {
            for(let M of testMats2) {
                M.isEchelon().should.be.true();
                M.isRREF().should.be.false();
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
        it('should not be in EF or RREF', () => {
            for(let M of testMats3) {
                M.isEchelon().should.be.false();
                M.isRREF().should.be.false();
            }
        });
    });
    describe('#trace', () => {
        it('should compute the sum of the diagonal entries', () =>
           mat([1,2,3],[4,5,6],[7,8,9]).trace.should.equal(1+5+9));
        it('should work for non-square matrices', () =>
           mat([1,2,3],[4,5,6]).trace.should.equal(1+5));
        it('should work for non-square matrices', () =>
           mat([1,2],[3,4],[5,6]).trace.should.equal(1+4));
    });
    describe('#rowScale(), #rowReplace(), #rowSwap()', () => {
        let M = mat([1,2,3],[4,5,6],[7,8,9]);
        it('should scale one row', () =>
           M.rowScale(1, 2).equals(mat(
               [1,2,3],[8,10,12],[7,8,9])).should.be.true());
        it('should add 2 x one row to another', () =>
           M.rowReplace(0, 2, 2).equals(mat(
               [15,18,21],[8,10,12],[7,8,9])).should.be.true());
        it('should swap two rows', () =>
           M.rowSwap(1, 2).equals(mat(
               [15,18,21],[7,8,9],[8,10,12])).should.be.true());
    });
    describe('#PLU(), #rref(), #*Basis(), #isFull*Rank()', () => {
        let testMats = [
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

        it('should factorize matrices correctly', () => {
            for(let {M} of testMats)
                M.PLU().should.factorize(M, 'PA=LU');
        });
        it('should compute the rref correctly', () => {
            for(let {M, rref} of testMats)
                M.rref().equals(rref, 1e-10).should.be.true();
        });
        it('should compute the row operations correctly', () => {
            for(let {M, rref} of testMats)
                M.rowOps().mult(M).equals(M.rref(), 1e-10)
                    .should.be.true();
        });
        it('should compute the null space correctly', () => {
            for(let {M, N} of testMats) {
                M.nullBasis().should.resemble(N);
                M.nullSpace().equals(new Subspace(N, {n: M.n})).should.be.true();
            }
        });
        it('should compute the column space correctly', () => {
            for(let {M, C} of testMats) {
                M.colBasis().should.resemble(C);
                M.colSpace().equals(new Subspace(C, {n: M.m})).should.be.true();
            }
        });
        it('should compute the row space correctly', () => {
            for(let {M, R} of testMats) {
                let RS = new Subspace(R);
                new Subspace(M.rowBasis()).equals(RS).should.be.true();
                M.rowSpace().equals(RS).should.be.true();
            }
        });
        it('should compute the left null space correctly', () => {
            for(let {M} of testMats) {
                let b = M.leftNullBasis();
                Vector.isLinearlyIndependent(b).should.be.true();
                b.should.have.length(M.m - M.rank());
                let z = Vector.zero(M.n);
                let MT = M.transpose;
                for(let v of b)
                    MT.apply(v).equals(z, 1e-10).should.be.true();
            }
        });
        it('should know whether the matrix has foll row and column rank', () => {
            for(let {M, N, C} of testMats) {
                M.isFullRowRank().should.equal(C.length == M.m);
                M.isFullColRank().should.equal(N.length == 0);
            }
        });
    });
    describe('#solve*()', () => {
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
                M.solve(b).should.solve(M, b);
                let x = M.solveShortest(b);
                x.should.solve(M, b);
                M.nullSpace().isOrthogonalTo(x).should.be.true();
                M.solveLeastSquares(b).should.solve(M, b);
                M.solveLeastSquaresShortest(b).equals(x, 1e-10).should.be.true();
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
        it('should have no solution', () => {
            for(let {M, b} of testNoSoln)
                should(M.solve(b)).be.null();
        });
        it('should have a least-squares solution', () => {
            for(let {M, b} of testNoSoln) {
                should(M.solve(b)).be.null();
                M.solveLeastSquares(b).should.solve(M, b, true);
                let x = M.solveLeastSquaresShortest(b);
                x.should.solve(M, b, true);
                M.nullSpace().isOrthogonalTo(x).should.be.true();
            }
        });
        it('should throw for incompatible dimensions', () => {
            let M = mat([1,2],[3,4]);
            M.solve.bind(M, vec(1,2,3)).should.throw(/Incompatible/);
        });
    });
    describe('#inverse()', () => {
        let testMats = [
            mat([3, 4],
                [5, 6]),
            mat([0,  1, 2],
                [1,  0, 3],
                [4, -3, 8]),
            mat([ 1,  7,  4, 2],
                [ 3, 11,  9, 5],
                [-2, -3,  3, 3],
                [ 7,  8, -8, 9])
        ];
        it('should compute the inverse', () => {
            for(let M of testMats) {
                M.mult(M.inverse()).equals(Matrix.identity(M.n), 1e-10)
                    .should.be.true();
                M.isInvertible().should.be.true();
            }
        });
        it('should throw for non-square matrices', () => {
            let M = mat([1, 2, 3], [2, 4, 5]);
            M.inverse.bind(M).should.throw(/non-square/);
            M.isSingular().should.be.true();
        });
        it('should throw for singular matrices', () => {
            let M = mat([ 1, -2, -1],
                        [-1,  5,  6],
                        [ 5, -4,  5]);
            M.inverse.bind(M).should.throw(/singular/);
            M.isSingular().should.be.true();
        });
    });
    describe('#det()', () => {
        it('should compute the determinant (1x1)', () =>
           Matrix.identity(1, 3).det().should.be.approximately(3, 1e-10));
        it('should compute the determinant (2x2)', () =>
           mat([3, 4], [5, 6]).det().should.be.approximately(-2, 1e-10));
        it('should compute the determinant (3x3#1)', () =>
           mat([0,  1, 2],
               [1,  0, 3],
               [4, -3, 8]).det().should.be.approximately(-2, 1e-10));
        it('should compute the determinant (3x3#2)', () =>
           mat([ 1, -2, -1],
               [-1,  5,  6],
               [ 5, -4,  5]).det().should.equal(0));
        it('should compute the determinant (4x4)', () =>
           mat([ 1,  7,  4, 2],
               [ 3, 11,  9, 5],
               [-2, -3,  3, 3],
               [ 7,  8, -8, 9]).det().should.be.approximately(-1329, 1e-10));
        it('should throw for non-square matrices', () => {
            let M = mat([1,2,3],[4,5,6]);
            M.det.bind(M).should.throw(/non-square/);
        });
    });
    describe('#isOrthogonal()', () => {
        it('detects orthogonal matrices', () =>
           mat([1,1],[1,-1]).scale(1/Math.sqrt(2)).isOrthogonal()
                   .should.be.true());
        it('detects non-orthogonal matrices', () =>
           mat([1/Math.sqrt(2), 1],[1/Math.sqrt(2), 0]).isOrthogonal()
                   .should.be.false());
        it('says non-square matrices are not orthogonal', () =>
           mat([1/Math.sqrt(2)],[1/Math.sqrt(2)]).isOrthogonal()
                   .should.be.false());
    });
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
                M.QR().should.factorize(M, 'QR');
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
    describe('#charpoly()', () => {
        it('should compute the characteristic polynomial of a 1x1 matrix', () => {
            mat([3]).charpoly().should.eql([3]);
        });
        it('should compute the characteristic polynomial of a 2x2 matrix', () => {
            mat([1,2],[3,4]).charpoly().should.eql([-5, -2]);
        });
        it('should compute the characteristic polynomial of a 3x3 matrix', () => {
            mat([1, 6,4],
                [2,-1,3],
                [5, 0,1]).charpoly().should.eql([1, 33, 97]);
        });
        it('should compute the characteristic polynomial of a 4x4 matrix', () => {
            mat([2, 1, 1, 0],
                [4, 3, 3, 1],
                [8, 7, 9, 5],
                [6, 7, 9, 8]).charpoly().should.eql([-22, 78, -50, 8]);
        });
        it('should throw for non-square matrices', () => {
            let M = mat([1,2,3],[4,5,6]);
            M.charpoly.bind(M).should.throw(/non-square/);
        });
    });
    describe('#eigenvalues()', () => {
        it('should compute eigenvalues for 1x1 matrices', () =>
           mat([3]).eigenvalues().should.eql([[3, 1]]));
        it('should compute eigenvalues for 2x2 matrices', () => {
            mat([1,1],[1,1]).eigenvalues()
                .should.eql([[0,1], [2,1]]);
            mat([1,1],[0,1]).eigenvalues().should.eql([[1,2]]);
            mat([1,1],[-1,1]).eigenvalues()
                .should.resemble([[new Complex(1, 1), 1], [new Complex(1, -1), 1]]);
        });
        it('should compute eigenvalues for 3x3 matrices', () => {
            Matrix.identity(3, 3).eigenvalues().should.resemble([[3, 3]]);
            mat([11/13, 22/39,  2/39],
                [-4/13, 83/39,  4/39],
                [-1/13, 11/39, 40/39]).eigenvalues()
                .should.resemble([[1, 2], [2, 1]]);
            mat([    1,   1/2,     0],
                [-4/13, 83/39,  4/39],
                [ 5/13,  7/78, 34/39]).eigenvalues()
                .should.resemble([[1, 2], [2, 1]]);
            mat([23/13,  53/78, -10/39],
                [-4/13, 122/39,   4/39],
                [-4/13,  49/78,  43/39]).eigenvalues()
                .should.resemble([[1, 1], [2, 1], [3, 1]]);
            mat([43/13,   7/13, -88/13],
                [ 6/13,  17/13, -28/13],
                [24/13, -23/13,  57/13]).eigenvalues()
                .should.resemble([[1, 1],
                                  [new Complex(4, -3), 1],
                                  [new Complex(4,  3), 1]]);
        });
    });
    describe('#eigenspace()', () => {
        it('should work for 1x1 matrices', () =>
           mat([3]).eigenspace(3).equals(Subspace.Rn(1)).should.be.true());
        it('should work for 2x2 matrices with distinct eigenvalues', () => {
            let M = mat([1, 1], [1, 1]);
            M.eigenspace(0).equals(new Subspace([[1,-1]])).should.be.true();
            M.eigenspace(2).equals(new Subspace([[1, 1]])).should.be.true();
            // test caching
            M.eigenspace(2.001, .1).equals(
                new Subspace([[1, 1]])).should.be.true();
            M.eigenspace.bind(M, 2.001).should.throw(/not an eigenvalue/);
        });
        it('should work for 2x2 matrices with a complex eigenvalue', () => {
            let M = mat([-15/7, 52/7], [-40/7, 57/7]);
            let λ = new Complex(3, 4);
            let B = M.eigenspace(λ);
            B.should.have.length(1);
            B.should.be.cplxEigenvectors(M, λ);
            λ.conj();
            B = M.eigenspace(λ);
            B.should.have.length(1);
            B.should.be.cplxEigenvectors(M, λ);
        });
        it('should work for 2x2 matrices with a repeated eigenvalue', () => {
            let M = mat([1,1],[0,1]);
            M.eigenspace(1).equals(new Subspace([[1, 0]])).should.be.true();
            M = mat([2,0],[0,2]);
            M.eigenspace(2).equals(Subspace.Rn(2)).should.be.true();
        });
        it('should work for 3x3 matrices with distinct eigenvalues', () => {
            let M = mat([23/13,  53/78, -10/39],
                        [-4/13, 122/39,   4/39],
                        [-4/13,  49/78,  43/39]);
            M.eigenspace(1).equals(new Subspace([[1, 0, 3]]))
                .should.be.true();
            M.eigenspace(2).equals(new Subspace([[7, 2, -1]]))
                .should.be.true();
            M.eigenspace(3).equals(new Subspace([[2, 4, 1]]))
                .should.be.true();
        });
        it('should work for 3x3 matrices with repeated eigenvalues', () => {
            let M = mat([11/13, 22/39,  2/39],
                        [-4/13, 83/39,  4/39],
                        [-1/13, 11/39, 40/39]);
            M.eigenspace(1).equals(new Subspace(
                [[1, 0, 3], [7, 2, -1]])).should.be.true();
            M.eigenspace(2).equals(new Subspace([[2, 4, 1]]))
                .should.be.true();
            M = mat([    1,   1/2,     0],
                    [-4/13, 83/39,  4/39],
                    [ 5/13,  7/78, 34/39]);
            M.eigenspace(1).equals(new Subspace([[1, 0, 3]]))
                .should.be.true();
            M.eigenspace(2).equals(new Subspace([[2, 4, 1]]))
                .should.be.true();
            Matrix.identity(3,3).eigenspace(3).equals(Subspace.Rn(3))
                .should.be.true();
        });
        it('should work for 3x3 matrices with a complex eigenvalue', () => {
            let M = mat([33, -23,   9],
                        [22,  33, -23],
                        [19,  14,  50]).scale(1/29);
            M.eigenspace(2).equals(new Subspace([[2, -1, 3]]))
                .should.be.true();
            let λ = new Complex(1, 1);
            let B = M.eigenspace(λ);
            B.should.have.length(1);
            B.should.be.cplxEigenvectors(M, λ);
            λ.conj();
            B = M.eigenspace(λ);
            B.should.have.length(1);
            B.should.be.cplxEigenvectors(M, λ);
        });
        it('should work for 4x4 matrices with distinct complex eigenvalues', () => {
            let M = mat([-1226,   230,  1166, -989],
                        [ 1530,  -192,  -786,  932],
                        [ 8938, -1856, -6401, 6438],
                        [12222, -2471, -9114, 8883]).scale(1/133);
            let λ = new Complex(3, 4);
            let B = M.eigenspace(λ);
            B.should.have.length(1);
            B.should.be.cplxEigenvectors(M, λ);
            λ.conj();
            B = M.eigenspace(λ);
            B.should.have.length(1);
            B.should.be.cplxEigenvectors(M, λ);
            λ = new Complex(1, 1);
            B = M.eigenspace(λ);
            B.should.have.length(1);
            B.should.be.cplxEigenvectors(M, λ);
            λ.conj();
            B = M.eigenspace(λ);
            B.should.have.length(1);
            B.should.be.cplxEigenvectors(M, λ);
        });
        it('should work for 4x4 matrices with repeated complex eigenvalues', () => {
            let M = mat([ -213,   396,   208, -160],
                        [-2616,  1059,  1056, -976],
                        [10256, -2036, -7221, 7212],
                        [10968, -2528, -8088, 7971]).scale(1/133);
            let λ = new Complex(3, 4);
            let B = M.eigenspace(λ);
            B.should.have.length(2);
            B.should.be.cplxEigenvectors(M, λ);
            λ.conj();
            B = M.eigenspace(λ);
            B.should.have.length(2);
            B.should.be.cplxEigenvectors(M, λ);

            M = mat([ -126,   371,   168, -119],
                    [-2442,  1009,   976, -894],
                    [10604, -2136, -7381, 7376],
                    [11229, -2603, -8208, 8094]).scale(1/133);
            λ = new Complex(3, 4);
            B = M.eigenspace(λ);
            B.should.have.length(1);
            B.should.be.cplxEigenvectors(M, λ);
            λ.conj();
            B = M.eigenspace(λ);
            B.should.have.length(1);
            B.should.be.cplxEigenvectors(M, λ);
        });
        it('should work for 4x4 matrices with real and complex eigenvalues', () => {
            let M = mat([  276,    17,   240, -113],
                        [ 1419,   138, -1056, 1109],
                        [ 9654, -2022, -6773, 6806],
                        [10827, -2249, -8280, 8088]).scale(1/133);
            let λ = new Complex(3, 4);
            let B = M.eigenspace(λ);
            B.should.have.length(1);
            B.should.be.cplxEigenvectors(M, λ);
            λ.conj();
            B = M.eigenspace(λ);
            B.should.have.length(1);
            B.should.be.cplxEigenvectors(M, λ);
            M.eigenspace(3).equals(new Subspace([[4,9,0,-3]]))
                .should.be.true();
            M.eigenspace(4).equals(new Subspace([[3,-3,2,-3]]))
                .should.be.true();

            M = mat([  363,    -8,   200,  -72],
                    [ 2463,  -162, -1536, 1601],
                    [ 9480, -1972, -6693, 6724],
                    [10827, -2249, -8280, 8088]).scale(1/133);
            λ = new Complex(3, 4);
            B = M.eigenspace(λ);
            B.should.have.length(1);
            B.should.be.cplxEigenvectors(M, λ);
            λ.conj();
            B = M.eigenspace(λ);
            B.should.have.length(1);
            B.should.be.cplxEigenvectors(M, λ);
            M.eigenspace(3).equals(new Subspace([[4,9,0,-3],]))
                .should.be.true();

            M.eigenspace.bind(M, new Complex(3, 5)).should.throw(/not an eigenvalue/);
        });
    });
    describe('#diagonalize()', () => {
        it('should work for 1x1 matrices', () => {
            let M = mat([3]);
            M.diagonalize().should.diagonalize(M);
            M.diagonalize({block: true}).should.diagonalize(M);
        });
        it('should diagonalize 2x2 matrices', () => {
            let testMats = [
                mat([1, 1], [1, 1]),
                mat([2, 0], [0, 2])
            ];
            for(let M of testMats) {
                M.diagonalize().should.diagonalize(M);
                M.diagonalize({block: true}).should.diagonalize(M);
            }
            mat([1,1],[0,1]).isDiagonalizable().should.be.false();
        });
        it('should diagonalize 3x3 matrices', () => {
            let testMats = [
                mat([23/13,  53/78, -10/39],
                    [-4/13, 122/39,   4/39],
                    [-4/13,  49/78,  43/39]),
                mat([11/13, 22/39,  2/39],
                    [-4/13, 83/39,  4/39],
                    [-1/13, 11/39, 40/39]),
                Matrix.identity(3, 3)
            ];
            for(let M of testMats) {
                M.diagonalize().should.diagonalize(M);
                M.diagonalize({block: true}).should.diagonalize(M);
            }
            mat([    1,   1/2,     0],
                [-4/13, 83/39,  4/39],
                [ 5/13,  7/78, 34/39]).isDiagonalizable().should.be.false();
        });
        it('should diagonalize 4x4 matrices with given eigenvalues', () => {
            let testMats = [{
                M: mat([1482, -247, -741,  703],
                       [ 963,  214, -828,  729],
                       [1572, -360, -842, 1016],
                       [ 105,  -21, -273,  476]).scale(1/133),
                ev: [[1,1],[2,1],[3,1],[4,1]]
            }, {
                M: mat([ 374, 132, 212, -186],
                       [ -18,  76, -12,   38],
                       [-132, -40, -32,   92],
                       [ 198, 116, 132,  -26]).scale(1/56),
                ev: [[1,2],[2,1],[3,1]]
            }, {
                M: mat([344, 128, 192, -160],
                       [ 72,  88,  48,  -40],
                       [-72, -32,   8,   40],
                       [288, 128, 192, -104]).scale(1/56),
                ev: [[1,3],[3,1]]
            }, {
                M: Matrix.identity(4, 3),
                ev: [[3,4]]
            }];
            for(let {M, ev} of testMats) {
                M.hintEigenvalues(...ev);
                M.diagonalize().should.diagonalize(M);
                M.diagonalize({block: true}).should.diagonalize(M);
            }
        });
        it('should block-diagonalize 2x2 matrices with a complex eigenvalue', () => {
            let testMats = [
                mat([2, -1], [2, 0]),
                mat([-Math.sqrt(3)+1, -2], [1, -Math.sqrt(3)-1])
            ];
            for(let M of testMats) {
                let {C, D} = M.diagonalize({block: true});
                ({C, D}).should.diagonalize(M, true);
                D[0][0].should.be.approximately(D[1][1], 1e-10);
                D[0][1].should.be.approximately(-D[1][0], 1e-10);
            }
        });
        it('should block-diagonalize 3x3 matrices with a complex eigenvalue', () => {
            let M = mat([33, -23,   9],
                        [22,  33, -23],
                        [19,  14,  50]).scale(1/29);
            let {C, D} = M.diagonalize({block: true});
            ({C, D}).should.diagonalize(M, true);
            D[0][1].should.equal(0);
            D[0][2].should.equal(0);
            D[1][0].should.equal(0);
            D[2][0].should.equal(0);
            D[1][1].should.be.approximately(D[2][2], 1e-10);
            D[1][2].should.be.approximately(-D[2][1], 1e-10);
        });
        it('should block-diagonalize 4x4 matrices with repeated complex eigenvalues', () => {
            let M = mat([ -213,   396,   208, -160],
                        [-2616,  1059,  1056, -976],
                        [10256, -2036, -7221, 7212],
                        [10968, -2528, -8088, 7971]).scale(1/133);
            let λ = new Complex(3, 4);
            M.hintEigenvalues([λ, 2], [λ.clone().conj(), 2]);
            M.diagonalize({block: true}).should.diagonalize(M, true);

            M = mat([ -126,   371,   168, -119],
                    [-2442,  1009,   976, -894],
                    [10604, -2136, -7381, 7376],
                    [11229, -2603, -8208, 8094]).scale(1/133);
            M.hintEigenvalues([λ, 2], [λ.clone().conj(), 2]);
            should(M.diagonalize({block: true})).be.null();
        });
        it('should block-diagonalize 4x4 matrices with real and complex eigenvalues', () => {
            let M = mat([  276,    17,   240, -113],
                        [ 1419,   138, -1056, 1109],
                        [ 9654, -2022, -6773, 6806],
                        [10827, -2249, -8280, 8088]).scale(1/133);
            let λ = new Complex(3, 4);
            M.hintEigenvalues([3, 1], [4, 1], [λ, 1], [λ.clone().conj(), 1]);
            M.diagonalize({block: true}).should.diagonalize(M, true);

            M = mat([ 1360,  480,  1504, -1048],
                    [-2288, -824, -2608,  1856],
                    [ 1476,  600,  1712, -1212],
                    [ 2288,  880,  2608, -1800]).scale(1/56);
            M.hintEigenvalues([1, 2], [λ, 1], [λ.clone().conj(), 1]);
            M.diagonalize({block: true}).should.diagonalize(M, true);
        });
        it('should orthogonally diagonalize symmetric matrices', () => {
            let testMats = [
                mat([3]),
                mat([1,2], [2,1]),
                mat([4,3], [3,4]),
                mat([1,2,3],
                    [2,4,7],
                    [3,7,9]),
                mat([-1, 5, 4],
                    [ 5,-2, 3],
                    [ 4, 3,-9]),
                mat([1.5019930033897762998393795730931589398,
                     0.44200755514697981035225502694854564642,
                     0.23372066474849376426747114556730413180],
                    [0.44200755514698000019419379151116544908,
                     1.3891900434622460804985772933506269552,
                     0.20579230968403700444145102672004495632],
                    [0.23372066474849394525106620858971434640,
                     0.20579230968403707541083109726046249083,
                     1.1088169531479776196620431335562141050])
            ];
            for(let M of testMats) {
                let {C, D} = M.diagonalize({ortho: true});
                ({C, D}).should.diagonalize(M);
                C.isOrthogonal().should.be.true();
            }
            let M = testMats.pop();
            M.eigenvalues().should.resemble([[1, 2], [2, 1]]);
            M.eigenspace(1).dim.should.equal(2);
        });
    });
});


describe('Subspace', () => {
    describe('#constructor()', () => {
        it('should construct a subspace from an Array of vectors', () =>
           new Subspace([[1,2,3], [4,5,6]]).generators
                   .should.eql(mat([1,4],[2,5],[3,6])));
        it('should construct a subspace from a Matrix', () => {
            let M = mat([1,4],[2,5],[3,6]);
            new Subspace(M).generators.should.eql(M);
        });
    });
    describe('#toString()', () => {
        let V = new Subspace([[10,2],[3,4], [-1,1]]);
        it('should have 4 decimal places by default', () =>
           V.toString().should.eql("Subspace of R^2 spanned by\n"
                                   + "10.0000 3.0000 -1.0000\n"
                                   + " 2.0000 4.0000  1.0000"));
        it('can have other precision', () =>
           V.toString(2).should.eql("Subspace of R^2 spanned by\n"
                                    + "10.00 3.00 -1.00\n"
                                    + " 2.00 4.00  1.00"));
    });
    describe('#basis(), #dim, #isMaximal(), #isZero()', () => {
        it('should find a full basis', () => {
            let V = new Subspace(mat([ 0,  3, -6,  6,  4, -5],
                                     [ 3, -7,  8, -5,  8,  9],
                                     [ 3, -9, 12, -9,  6, 15]));
            V.basis().should.eql(mat([0,  3, 4],
                                     [3, -7, 8],
                                     [3, -9, 6]));
            V.dim.should.equal(3);
            V.isMaximal().should.be.true();
        });
        it('should find a proper basis', () => {
            let V = new Subspace(mat([ 0, -3, -6,  4,  9],
                                     [-1, -2, -1,  3,  1],
                                     [-2, -3,  0,  3, -1],
                                     [ 1,  4,  5, -9, -7]));
            V.basis().should.eql(mat([ 0, -3,  4],
                                     [-1, -2,  3],
                                     [-2, -3,  3],
                                     [ 1,  4, -9]));
            V.dim.should.equal(3);
            V.isMaximal().should.be.false();
        });
        it('should detect the zero space', () => {
            Subspace.zero(3).isZero().should.be.true();
            new Subspace([[.01, .01, .01], [-.01, .01, .01]])
                .isZero().should.be.false();
            new Subspace([[.01, .01, .01], [-.01, .01, .01]])
                .isZero(.1).should.be.true();
        });
    });
    describe('#ONbasis()', () => {
        it('should find an orthonormal basis', () => {
            let V = new Subspace(mat([ 0, -3, -6,  4,  9],
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
        let V = new Subspace([[1,1,1], [1,0,0]]);
        let W = new Subspace([[1,-1,1], [0,1,0]]);
        it('should add the subspaces', () =>
           V.add(W).equals(Subspace.Rn(3)).should.be.true());
        it('should intersect the subspaces', () =>
           V.intersect(W).equals(new Subspace([[1,1,1]])).should.be.true());
    });
    describe('#projectionMatrix()', () => {
        it('should find a nontrivial projection matrix', () => {
            let V = new Subspace(mat([ 0, -3, -6,  4,  9],
                                     [-1, -2, -1,  3,  1],
                                     [-2, -3,  0,  3, -1],
                                     [ 1,  4,  5, -9, -7]));
            let P = V.projectionMatrix();
            P.mult(P).equals(P, 1e-10).should.be.true();
            P.colSpace().equals(V).should.be.true();
        });
        it('should return the zero matrix for the zero subspace', () =>
           Subspace.zero(3).projectionMatrix().equals(Matrix.zero(3), 1e-10)
                   .should.be.true());
        it('should return the identity matrix for the full subspace', () =>
           Subspace.Rn(3).projectionMatrix().equals(Matrix.identity(3), 1e-10)
                   .should.be.true());
    });
    describe('#orthoDecomp()', () => {
        let V = new Subspace(mat([ 0, -3, -6,  4,  9],
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
    });
    describe('#contains()', () => {
        it('should contain every vector when it is maximal', () => {
            let V = new Subspace(mat([ 0,  3, -6,  6,  4, -5],
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
            let vecs = [...M.cols()];
            let V = new Subspace(vecs);
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
            let V = new Subspace(mat([ 0, -3, -6,  4,  9],
                                     [-1, -2, -1,  3,  1],
                                     [-2, -3,  0,  3, -1],
                                     [ 1,  4,  5, -9, -7]).transpose);
            let W = new Subspace(mat([1, 0, -3, 0,  5],
                                     [0, 0,  0, 1,  0]).transpose);
            W.isSubspaceOf(V).should.be.true();
            V.isSubspaceOf(W).should.be.false();
        });
    });
    describe('#orthoComplement()', () => {
        it('should compute the orthogonal complement', () => {
            let V = new Subspace(mat([ 0, -3, -6,  4,  9],
                                     [-1, -2, -1,  3,  1],
                                     [-2, -3,  0,  3, -1],
                                     [ 1,  4,  5, -9, -7]));
            let W = V.orthoComplement();
            (V.dim + W.dim).should.equal(4);
            for(let v of V.generators.cols())
                W.isOrthogonalTo(v).should.be.true();
        });
        it('should compute the orthogonal complement of {0}', () =>
           Subspace.zero(3).orthoComplement().isMaximal().should.be.true());
        it('should compute the orthogonal complement of R^n', () =>
           Subspace.Rn(3).orthoComplement().isZero().should.be.true());
    });
});
