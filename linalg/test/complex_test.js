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

'use strict';

import { Complex } from "../src/linalg.js";

import should from 'should';
import './lib/resemble.js';

const C = (a, b=0) => new Complex(a, b);


describe('Complex', () => {
    const π = Math.PI;

    describe('#fromPolar()', () => {
        it('should create positive real numbers with θ=0', () =>
           Complex.fromPolar(1, 0).should.resemble(C(1, 0)));
        it('should create negative real numbers with θ=π', () =>
           Complex.fromPolar(1, π).should.resemble(C(-1, 0)));
        it('should create purely imaginary numbers with θ=π/2', () =>
           Complex.fromPolar(1, π/2).should.resemble(C(0, 1)));
        it('should create a cube root of -8 with r=2 and θ=π/3', () => {
            let z = Complex.fromPolar(2, π/3), w = z.clone();
            z.should.resemble(C(1, Math.sqrt(3)));
            z.mult(w).mult(w).should.resemble(C(-8, 0));
        });
    });
    describe('#i()', () => {
        it('should be 0 + 1i', () =>
           Complex.i.should.eql(C(0, 1)));
        it('should be a square root of -1', () =>
           Complex.i.mult(Complex.i).should.eql(C(-1,0)));
    });
    describe('#equals()', () => {
        let z = C(3, 0), w = C(3.01, 0.01);
        it('should compare equal to a real number if Im==0', () =>
           z.equals(3).should.be.true());
        it('should compare not equal to a real number if Im!=0', () =>
           w.equals(3).should.be.false());
        it('should compare equal when ε>0', () =>
           w.equals(3, .05).should.be.true());
    });
    describe('#toString()', () => {
        let z = C(3, 4);
        it('should have 4 decimal places by default', () =>
           z.toString().should.eql("3.0000 + 4.0000 i"));
        it('can have other precision', () =>
           z.toString(2).should.eql("3.00 + 4.00 i"));
        it('subtracts a negative imaginary part', () =>
           z.conj().toString(2).should.eql("3.00 - 4.00 i"));
    });
    describe('#Re(), #im(), #mod()', () => {
        let z = C(3, 4);
        it('should have the correct real part', () =>
           z.Re.should.equal(3));
        it('should have the correct imaginary part', () =>
           z.Im.should.equal(4));
        it('should have the correct modulus', () =>
           z.mod.should.be.approximately(5, 1e-10));
        it('should set the real part', () => {
            z.Re = 1;
            z.should.eql(C(1, 4));
        });
        it('should set the imaginary part', () => {
            z.Im = 1;
            z.should.eql(C(1, 1));
        });
    });
    describe('#arg()', () => {
        it('should compute the modulus of a 6th root of unity', () =>
           C( 1/2,  Math.sqrt(3)/2).arg
                   .should.be.approximately(   π/3, 1e-10));
        it('should compute the modulus of a 3rd root of unity', () =>
           C(-1/2,  Math.sqrt(3)/2).arg
                   .should.be.approximately( 2*π/3, 1e-10));
        it('should compute the modulus of another 6th root of unity', () =>
           C( 1/2, -Math.sqrt(3)/2).arg
                   .should.be.approximately(  -π/3, 1e-10));
        it('should compute the modulus of another 3rd root of unity', () =>
           C(-1/2, -Math.sqrt(3)/2).arg
                   .should.be.approximately(-2*π/3, 1e-10));
    });
    describe('#conj()', () => {
        let z = C(1, 2);
        it('should negate the imaginary part', () =>
           z.conj().should.eql(C(1, -2)));
    });
    describe('#add(), #sub()', () => {
        it('should add componentwise', () =>
           C(3, 4).add(C(1, 2)).should.eql(C(4, 6)));
        it('should add componentwise with factor', () =>
           C(3, 4).add(C(1, 2), -1).should.eql(C(2, 2)));
        it('should subtract componentwise', () =>
           C(3, 4).sub(C(1, 2)).should.eql(C(2, 2)));
        it('should add real numbers correctly', () =>
           C(3, 4).add(2).should.eql(C(5, 4)));
        it('should add real numbers with factor', () =>
           C(3, 4).add(2, -1).should.eql(C(1, 4)));
        it('should subtract real numbers correctly', () =>
           C(3, 4).sub(2).should.eql(C(1, 4)));
    });
    describe('#mult()', () => {
        it('should multiply to -5 + 10 i', () =>
           C(1, 2).mult(C(3, 4))
                   .should.eql(C(-5, 10)));
        it('should scalar-multiply to 2 + 4 i', () =>
           C(1, 2).mult(2).should.eql(C(2, 4)));
    });
    describe('#pow()', () => {
        it('should raise to the power 1/2', () => {
            let z = C(1, 2);
            z.mult(z);
            let w = z.pow(1/2);
            w.mult(w).should.resemble(z);
        });
        it('should raise to the power 1/3', () => {
            let z = C(1, 2);
            z.mult(z.clone().mult(z));
            let w = z.pow(1/3);
            w.mult(w.clone().mult(w)).should.resemble(z);
        });
    });
    describe('#recip()', () => {
        it('should compute the reciprocal', () =>
           C(3, 4).recip()
                   .should.eql(C(3/25, -4/25)));
        it('should throw when the number is zero', () => {
            let z = C(0);
            z.recip.bind(z).should.throw(/divide by zero/);
        });
    });
    describe('#div()', () => {
        it('should compute the quotient', () =>
           C(1, 2).div(C(3, 4))
                   .should.resemble(C(11/25, 2/25)));
        it('should compute the quotient by a scalar', () =>
           C(1, 2).div(3).should.eql(C(1/3, 2/3)));
        it('should throw when the denominator is zero', () => {
            let z = C(1, 2);
            z.div.bind(z, C(0)).should.throw(/divide by zero/);
        });
        it('should throw when the denominator is the scalar zero', () => {
            let z = C(1, 2);
            z.div.bind(z, 0).should.throw(/divide by zero/);
        });
    });
});
