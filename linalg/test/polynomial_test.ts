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

import chai from 'chai';
chai.should();

import './lib/resemble';

import Polynomial from "../src/polynomial";
import Complex from "../src/complex";


const poly = (...coeffs: number[]) => new Polynomial(...coeffs);
const fromRoots = Polynomial.fromRoots;
const C = (a: number | Complex, b=0) => new Complex(a, b);
const ζ = new Complex(-1/2, Math.sqrt(3)/2);


describe('Polynomial', () => {
    describe('#create(), @@iterator()', () => {
        it("should create polynomials of degree 2", () =>
           Array.from(poly(1, 2, 3)).should.eql([1, 2, 3]));
        it("should create polynomials of degree 0", () =>
           Array.from(poly(1)).should.eql([1]));
        it("should create polynomials of degree -Infinity", () =>
           Array.from(poly()).should.eql([]));
        it("should strip leading zeros", () =>
           Array.from(poly(0, 0, 1, 2)).should.eql([1, 2]));
    });
    describe('#coeff()', () => {
        it("should get the correct coefficient", () =>
            poly(1, 2, 3).coeff(2).should.equal(1));
        it("should return zero for zero coefficients", () =>
            poly(1, 2, 3).coeff(4).should.equal(0));
    });
    describe('#fromRoots()', () => {
        it("should create polynomials with real roots", () => {
            fromRoots(1, -1).equals(poly(1, 0, -1)).should.be.true;
            fromRoots(1, 1, 1).equals(poly(1, -3, 3, -1)).should.be.true;
        });
        it("should create polynomials with complex roots", () => {
            fromRoots(1, C(0,1), C(0,-1)).equals(poly(1, -1, 1, -1))
                .should.be.true;
            fromRoots(
                -1, ζ, ζ.clone().conj(),
                ζ.clone().mult(-1),
                ζ.clone().conj().mult(-1))
                .equals(poly(1, 1, 1, 1, 1, 1), 1e-10).should.be.true;
        });
        it("should create polynomials with multiple roots", () => {
            fromRoots([1, 4]).equals(poly(1, -4, 6, -4, 1)).should.be.true;
            fromRoots([C(0,1), 2], [C(0,-1), 2], 1)
                .equals(poly(1, -1, 2, -2, 1, -1)).should.be.true;
        });
    });
    describe('#legendre()', () => {
        const legendres = [
            poly(1),
            poly(1, 0),
            poly(3, 0, -1).mult(1/2),
            poly(5, 0, -3, 0).mult(1/2),
            poly(35, 0, -30, 0, 3).mult(1/8),
            poly(63, 0, -70, 0, 15, 0).mult(1/8),
            poly(231, 0, -315, 0, 105, 0, -5).mult(1/16),
            poly(429, 0, -693, 0, 315, 0, -35, 0).mult(1/16),
            poly(6435, 0, -12012, 0, 6930, 0, -1260, 0, 35).mult(1/128),
            poly(12155, 0, -25740, 0, 18018, 0, -4620, 0, 315, 0).mult(1/128),
            poly(46189, 0, -109395, 0, 90090, 0, -30030, 0, 3465, 0, -63).mult(1/256)
        ];
        for(let n = 0; n < legendres.length; ++n) {
            it("should compute the Legendre polynomial of degree " + n, () =>
               Polynomial.legendre(n).equals(legendres[n]).should.be.true);
        }
    });
    describe('#deg', () => {
        it("should return the degree", () => {
            poly(1, 0, 0).deg.should.equal(2);
            poly(0).deg.should.equal(-Infinity);
        });
    });
    describe('#clone()', () => {
        let v = poly(1,2,3), w = v.clone();
        it('should create distinct objects', () =>
           v.should.not.equal(w));
        it('should create equal polynomials', () =>
           v.should.eql(w));
    });
    describe('#isZero()', () => {
        it("should know if a polynomial is nonzero", () =>
           poly(1, 0, 0).isZero().should.be.false);
        it("should know if a polynomial is zero", () =>
           poly(0, 0, 0).isZero().should.be.true);
    });
    describe('#equals()', () => {
        let p = poly(1, 0.01, -0.01, 0);
        let q = poly(1, 0, 0, 0);
        it("should know if polynomials are not equal", () =>
           p.equals(q).should.be.false);
        it("should know if polynomials are equal within a range", () =>
           p.equals(q, 0.05).should.be.true);
        it("should know if polynomials are equal", () =>
           q.equals(poly(1, 0, 0, 0)).should.be.true);
        it("should test polynomials of different degrees as not equal", () =>
           q.equals(poly(1, 0, 0)).should.be.false);
    });
    describe('#toString()', () => {
        it("should work with different precisions", () => {
            poly(3, 2, -1).toString().should.equal("3.0000 x^2 + 2.0000 x - 1.0000");
            poly(3, 2, -1).toString(2).should.equal("3.00 x^2 + 2.00 x - 1.00");
        });
        it("should work with different variables", () =>
           poly(3, 2, -1).toString(1, 'z').should.equal("3.0 z^2 + 2.0 z - 1.0"));
        it("should work for the zero polynomial", () =>
           poly().toString().should.equal("0.0000"));
        it("should work for polynomials of degree zero", () => {
            poly(1).toString(1).should.equal("1.0");
            poly(-1).toString(1).should.equal("-1.0");
            poly(2).toString(1).should.equal("2.0");
        });
        it("should work for polynomials of degree one", () => {
            poly(1, 1).toString(1).should.equal("x + 1.0");
            poly(-1, 0).toString(1).should.equal("-x");
            poly(2, -1).toString(1, 'z').should.equal("2.0 z - 1.0");
        });
        it("should work for polynomials of higher degree", () => {
            poly(1, 0, 1).toString(1).should.equal("x^2 + 1.0");
            poly(-1, 0, -1).toString(1).should.equal("-x^2 - 1.0");
            poly(1, 1, -1).toString(1).should.equal("x^2 + x - 1.0");
            poly(1, -1, 1).toString(1).should.equal("x^2 - x + 1.0");
            poly(2, -1, 1).toString(1).should.equal("2.0 x^2 - x + 1.0");
            poly(-2, 1, 1).toString(1).should.equal("-2.0 x^2 + x + 1.0");
        });
    });
    describe('#add() and sub()', () => {
        it('should add polynomials of the same degree', () =>
           poly(1,2,3).add(poly(4,5,6)).should.eql(poly(5,7,9)));
        it('should add polynomials different degrees', () => {
            poly(1,2,3).add(poly(5,6)).should.eql(poly(1,7,9));
            poly(2,3).add(poly(4,5,6)).should.eql(poly(4,7,9));
        });
        it('should add with a factor', () =>
           poly(1,2,3).add(poly(5,6), 2).should.eql(poly(1,12,15)));
        it('should cancel leading coefficients', () => {
            poly(1,2,3).add(poly(-1,0,0)).should.eql(poly(2,3));
            poly(1,2,3).add(poly(-1,-2,0)).should.eql(poly(3));
            poly(1,2,3).sub(poly(1,2,3)).should.eql(poly());
        });
    });
    describe('#mult()', () => {
        it('should multiply by the zero polynomial', () => {
            poly(1,2,3).mult(poly()).should.eql(poly());
            poly().mult(poly(1,2,3)).should.eql(poly());
        });
        it('should multiply polynomials of the same degree', () =>
           poly(1,2,3).mult(poly(4,5,6)).should.eql(poly(4, 13, 28, 27, 18)));
        it('should multiply polynomials of different degrees', () => {
            poly(1,2,3).mult(poly(4,5,6,7)).should.eql(poly(4, 13, 28, 34, 32, 21));
            poly(4,5,6,7).mult(poly(1,2,3)).should.eql(poly(4, 13, 28, 34, 32, 21));
            poly(1,2,3).mult(poly(2)).should.eql(poly(2,4,6));
            poly().mult(poly(2,3,4)).should.eql(poly());
        });
        it('should scale by a scalar', () => {
            poly(1,2,3).mult(2).should.eql(poly(2,4,6));
            poly(1,2,3).mult(0).should.eql(poly());
        });
    });
    describe('#pxPlusA', () => {
        it('should compute p*x + a', () => {
            poly(1,2,3).pxPlusA(4).should.eql(poly(1,2,3,4));
            poly(1,2,3).pxPlusA().should.eql(poly(1,2,3,0));
        });
    });
    describe('#div()', () => {
        it('should divide with remainder', () => {
            poly(1,-12,0,-42).div(poly(1,1,-3))
                .should.eql([poly(1,-13),poly(16,-81)]);
            poly(6,5,0,-7).div(poly(3,-2,-1))
                .should.eql([poly(2,3),poly(8,-4)]);
            poly(2,-5,3,7).div(poly(1,-2))
                .should.eql([poly(2,-1,1),poly(9)]);
            poly(4,-8,-1,5).div(poly(2,-1))
                .should.eql([poly(2,-3,-2), poly(3)]);
            poly(1,0,0,-1).div(poly(1,-1))
                .should.eql([poly(1,1,1), poly(0)]);
        });
        it('should throw when input has larger degree', () => {
            let p = poly(1,2);
            p.div.bind(p, poly(1,2,3)).should.throw(/larger degree/);
        });
        it('should throw when input is zero', () => {
            let p = poly(1,2);
            p.div.bind(p, poly()).should.throw(/zero polynomial/);
        });
    });
    describe('#monic()', () => {
        it('should scale by the leading coefficient', () =>
           poly(2,4,6).monic().should.eql(poly(1,2,3)));
    });
    describe('#pow()', () => {
        it('should compute powers correctly', () => {
            let p = poly(1, 1);
            p.pow(0).should.eql(poly(1));
            p.pow(1).should.eql(p);
            p.pow(2).should.eql(poly(1, 2, 1));
            p.pow(3).should.eql(poly(1, 3, 3, 1));
            p.pow(4).should.eql(poly(1, 4, 6, 4, 1));
            p.pow(5).should.eql(poly(1, 5, 10, 10, 5, 1));
        });
    });
    describe('#eval()', () => {
        it('should evaluate at real numbers', () => {
            let p = poly(1, 1, 1, 1);
            p.eval(1).should.equal(4);
            p.eval(-1).should.equal(0);
            p.eval(2).should.equal(15);
            poly().eval(1).should.equal(0);
        });
        it('should evaluate at complex numbers', () => {
            let p = poly(1, 1, 1, 1);
            poly().eval(C(0)).should.eql(C(0));
            p.eval(C(0,1)).should.eql(C(0));
            p.eval(C(1,1)).should.eql(C(0, 5));
        });
    });
    describe('#compose()', () => {
        let p = poly(1, 1, 1, 1);
        let q = poly(1, 1, 1);
        it('should compose correctly', () => {
            p.compose(q).should.eql(poly(1,3,7,9,10,6,4));
            q.compose(p).should.eql(poly(1,2,3,5,4,3,3));
        });
        it('should handle corner cases', () => {
            p.compose(poly(1)).should.eql(poly(4));
            p.compose(poly()).should.eql(poly(1));
            poly().compose(p).should.eql(poly());
        });
    });
    describe('#derivative()', () => {
        it('should compute derivatives', () => {
            poly(1,1,1,1).derivative(0).should.eql(poly(1,1,1,1));
            poly(1,1,1,1).derivative().should.eql(poly(3,2,1));
            poly(1,1,1,1).derivative(2).should.eql(poly(6,2));
            poly(1,1,1,1).derivative(3).should.eql(poly(6));
            poly(1,1,1,1).derivative(4).should.eql(poly(0));
        });
        it('should handle corner cases', () => {
            poly().derivative().should.eql(poly());
            poly(1).derivative(2).should.eql(poly());
        });
    });

    describe('#factor()', () => {
        describe('factoring linear polynomials', () => {
            it('should find one root', () =>
               poly(1, -1).factor().should.eql([[1, 1]]));
        });
        describe('factoring quadratics', () => {
            it('should find a double root', () => {
                // x^2 - 2x + 1 = (x-1)^2
                poly(1, -2, 1).factor().should.resemble([[1, 2]]);
                poly(1, 0, 0).factor().should.eql([[0, 2]]);
            });
            it('should find simple real roots', () => {
                // x^2 - 1 = (x-1) (x+1)
                poly(1, 0, -1).factor().should.resemble([[-1, 1], [1, 1]]);
                // x^2 - 3x + 2 = (x-1) (x-2)
                poly(1, -3, 2).factor().should.resemble([[1, 1], [2, 1]]);
                poly(1, -2, 0).factor().should.eql([[0, 1], [2, 1]]);
            });
            it('should find complex roots', function() {
                poly(1, 0, 1).factor().should.resemble(
                    [[C(0, 1), 1], [C(0, -1), 1]]);
                poly(1, 1, 1).factor().should.resemble(
                    [[ζ, 1], [ζ.clone().conj(), 1]]);
            });
        });
        describe('factoring cubics', () => {
            it('should find a triple root', () => {
                // x^3 - 6x^2 + 12x - 8 = (x-2)^3
                poly(1, -6, 12, -8).factor().should.resemble([[2, 3]]);
            });
            it('should find double roots', () => {
                // x^3 - 4x^2 + 5x - 2 = (x-2) (x-1)^2
                poly(1, -4,  5, -2).factor().should.resemble([[1, 2], [2, 1]]);
                // x^3 + 4x^2 + 5x + 2 = (x+2) (x+1)^2
                poly(1,  4,  5,  2).factor().should.resemble([[-2, 1], [-1, 2]]);
            });
            it('should find simple real roots', () => {
                // x^3 + 2x^2 - x - 2 = (x+2) (x+1) (x-1)
                poly(1, 2, -1, -2).factor().should
                    .resemble([-2, -1, 1].map(x => [x, 1]));
                // x^3 +  x^2 - 4x - 4 = (x+2) (x+1) (x-2)
                poly(1, 1, -4, -4).factor().should
                    .resemble([-2, -1, 2].map(x => [x, 1]));
                // x^3 - 4x
                poly(1, 0, -4, 0).factor().should
                    .resemble([-2, 0, 2].map(x => [x, 1]));
                // x^3 + 4x
                poly(1, 0, 4, 0).factor().should
                    .resemble([0, C(0,-2), C(0,2)].map(x => [x, 1]));
            });
            it('should find complex Roots', () => {
                // x^3 - 2x^2 + x - 2 = (x-2) (x^2+1)
                poly(1, -2,  1, -2).factor().should.resemble(
                    [Complex.i, Complex.i.conj(), 2].map(x => [x, 1]));
                // x^3 - x^2 - x - 2 = (x-2) (x^2+x+1)
                poly(1, -1, -1, -2).factor().should.resemble(
                    [ζ, ζ.clone().conj(), 2].map(x => [x, 1]));
            });

        });
        describe('factoring quartics', () => {
            let xm1 = poly(1, -1);
            it('should find roots at zero', () => {
                // y^4 for y=x-1
                poly(1, 0, 0, 0, 0).compose(xm1).factor()
                    .should.resemble([[1, 4]]);
                // y^4 - 4y^2 for y=x-1
                poly(1, 0, -4, 0, 0).compose(xm1).factor()
                    .should.resemble([[-1, 1], [1, 2], [3, 1]]);
                // y^4 + 4y^2 for y=x-1
                poly(1, 0, 4, 0 ,0).compose(xm1).factor()
                    .should.resemble([[1, 2], [C(1, 2), 1], [C(1, -2), 1]]);
                // y (y^3 - 3y + 2) = y (y-1)^2 (y+2) for y=x-1
                poly(1, 0, -3, 2, 0).compose(xm1).factor()
                    .should.resemble([[-1, 1], [1, 1], [2, 2]]);
                // y (y^3 - 2y + 4) = y (y+2) (y^2 - 2y + 2) for y=x-1
                poly(1, 0, -2, 4, 0).compose(xm1).factor()
                    .should.resemble([-1, 1, C(2,1), C(2,-1)].map(x => [x, 1]));
                // y (y^3 - 2y - 4) = y (y-2) (y^2 + 2y + 2) for y=x-1
                poly(1, 0, -2, -4, 0).compose(xm1).factor()
                    .should.resemble([C(0, 1), C(0, -1), 1, 3].map(x => [x, 1]));
            });
            it('should factor biquadratics', () => {
                // y^4 - 2 y^2 + 1 for y=x-1
                poly(1, 0, -2, 0, 1).compose(xm1).factor()
                    .should.resemble([[0, 2], [2, 2]]);
                // y^4 + 2 y^2 + 1 for y=x-1
                poly(1, 0, 2, 0, 1).compose(xm1).factor()
                    .should.resemble([[C(1, 1), 2], [C(1, -1), 2]]);
                // y^4 + 6 y^2 + 25 for y=x-1
                poly(1, 0, 6, 0, 25).compose(xm1).factor()
                    .should.resemble([C(-1, 2), C(-1, -2), C(1, 2), C(1, -2)]
                                     .map(x => [x.add(1), 1]));
                // y^4 + 5 y^2 + 4 = (y^2+1) (y^2+4) for y=x-1
                poly(1, 0, 5, 0, 4).compose(xm1).factor()
                    .should.resemble([C(0, 1), C(0, -1), C(0, 2), C(0, -2)]
                                     .map(x => [x.add(1), 1]));
                // y^4 - 3 y^2 - 4 = (y^2+1) (y^2-4) for y=x-1
                poly(1, 0, -3, 0, -4).compose(xm1).factor()
                    .should.resemble([-1, C(1, 1), C(1, -1), 3].map(x => [x, 1]));
                // y^4 - 5 y^2 + 4 = (y^2-1) (y^2-4) for y=x-1
                poly(1, 0, -5, 0, 4).compose(xm1).factor()
                    .should.resemble([-2, -1, 1, 2].map(x => [x+1, 1]));
            });
            it('should handle triple roots', () => {
                // (y-1)^3 (y+3) for y=x-1
                poly(1, 0, -6, 8, -3).compose(xm1).factor()
                    .should.resemble([[-2, 1], [2, 3]]);
                // (y+1)^3 (y-3) for y=x-1
                poly(1, 0, -6, -8, -3).compose(xm1).factor()
                    .should.resemble([[0, 3], [4, 1]]);
            });
            it('should handle double roots', () => {
                // (y+3)^2 (y-2) (y-4) for y=x-1
                poly(1, 0, -19, -6, 72).compose(xm1).factor()
                    .should.resemble([[-2, 2], [3, 1], [5, 1]]);
                // (y-1)^2 (y+4) (y-2) for y=x-1
                poly(1, 0, -11, 18, -8).compose(xm1).factor()
                    .should.resemble([[-3, 1], [2, 2], [3, 1]]);
                // (y-1)^2 (y^2 + 2y + 2) for y=x-1
                poly(1, 0, -1, -2, 2).compose(xm1).factor()
                    .should.resemble([[C(0, 1), 1], [C(0, -1), 1], [2, 2]]);
                // (y+1)^2 (y^2 - 2y + 2) for y=x-1
                poly(1, 0, -1, 2, 2).compose(xm1).factor()
                    .should.resemble([[0, 2], [C(2, 1), 1], [C(2, -1), 1]]);
            });
            it('should handle simple roots', () => {
                // (y-1) (y-4) (y+2) (y+3) for y=x-1
                poly(1, 0, -15, -10, 24).compose(xm1).factor()
                    .should.resemble([-3, -2, 1, 4].map(x => [x+1, 1]));
                // (y-1) (y-3) (y^2 + 4y + 5) for y=x-1
                poly(1, 0, -8, -8, 15).compose(xm1).factor()
                    .should.resemble([C(-1,1), C(-1,-1), 2, 4].map(x => [x, 1]));
                // (y+1) (y+3) (y^2 - 4y + 5) for y=x-1
                poly(1, 0, -8, 8, 15).compose(xm1).factor()
                    .should.resemble([-2, 0, C(3, 1), C(3, -1)].map(x => [x, 1]));
                // (y^2 - 2y + 2) (y^2 + 2y + 5)for y=x-1
                poly(1, 0, 3, -6, 10).compose(xm1).factor()
                    .should.resemble([C(-1, 2), C(-1, -2), C(1, 1), C(1, -1)]
                                     .map(x => [x.add(1), 1]));
            });
        });
    });
});
