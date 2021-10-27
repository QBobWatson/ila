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

import { expect, Assertion } from 'chai';

// import Matrix from '../../src/matrix';
import Complex from '../../src/complex';


declare global {
    export namespace Chai {
        interface Assertion {
            resemble(obj: any, ε?: number): void;
        }
    }
}


function approximate(obj: any, other: any, ε: number): void {
    if (typeof obj === "number")
        expect(obj).to.be.approximately(other, ε);
    /* else if (obj instanceof Matrix) { */
    /*     expect(other).to.be.an.instanceOf(Matrix); */
    /*     expect(obj.m).to.equal(other.m); */
    /*     expect(obj.n).to.equal(other.n); */
    /*     approximate([...obj.rows()], [...other.rows()], ε); */
    /* }  */
    else if (obj instanceof Complex) {
        expect(other).to.be.an.instanceOf(Complex);
        approximate([obj.Re, obj.Im], [other.Re, other.Im], ε);
    } else {
        expect(obj).to.be.an('array');
        expect(other).to.be.an('array');
        expect(other).to.have.length(obj.length);
        for (let i = 0; i < obj.length; ++i)
            approximate(obj[i], other[i], ε);
    }
}

// Custom assertion for comparing arrays of numbers using 'approximately'
Assertion.addMethod(
    'resemble', function(this: typeof Assertion, other: any, ε = 1e-10) {
        approximate(this._obj, other, ε);
    });
