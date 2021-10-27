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

import { Vector } from "../src/linalg.js";

import './lib/resemble.js';
import should from 'should';

const vec = Vector.create;


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
