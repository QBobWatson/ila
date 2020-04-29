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

import should from 'should';
import Matrix from '../../src/matrix.js';


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
});
