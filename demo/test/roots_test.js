'use strict'; // -*- js2 -*-

import rpoly from '../lib/roots.js';
import { Complex } from '../lib/linalg.js';

import should from 'should';

function C(a,b=0) {
    if(a instanceof Complex)
        return a;
    return new Complex(a, b);
}
const Exp = Complex.fromPolar;
const π = Math.PI;


function range(a, b=undefined, step=1) {
    let min, max;
    if(b === undefined)
        [min, max] = [0, a];
    else
        [min, max] = [a, b];
    let g = function*() {
        if(step > 0) {
            for(let i = min; i <= max; i += step)
                yield i;
        } else {
            for(let i = max; i >= min; i += step)
                yield i;
        }
    };
    return [...g()];
}

function constant(x, count) {
    return new Array(count).fill(x);
}


class Polynomial extends Array {
    static fromRoots(...roots) {
        let P = [C(1), C(roots.pop()).scale(-1)];
        let n = 2;
        while(roots.length > 0) {
            let c = C(roots.pop());
            P.forEach((a, i) => a.mult(c).scale(-1).add((i < n-1 ? P[i+1] : 0)));
            P.unshift(C(1));
            n++;
        }
        return new Polynomial(...P.map(x => x.Re));
    }

    get deg() {
        return this.length - 1;
    }

    mult(P) {
        return range(0, this.deg + P.deg).map(
            i => range(Math.max(0,i-P.deg), Math.min(i,this.deg)).reduce(
                (a, j) => a + this[j]*P[i-j], 0));
    }

    eval(z) {
        if(z instanceof Complex)
            return this.reduce((a, x) => a.mult(z).add(x), C(0, 0));
        return this.reduce((a, x) => a*z + x, 0);
    }
}

const P = (...coeffs) => new Polynomial(...coeffs);
const PR = (...roots) => Polynomial.fromRoots(...roots);

const jt01 = a => P(1,-1,-(a**2),a**2)
      ({ p: [1,-1,-(a**2),a**2],
                     z: [[a,1],[-a,1],[1,1]] });
const jt07 = a => ({ z: [.001, .01, .1, C(.1, a), C(.1,-a), 1, -10] });
const jt10 = a => ({ z: [a, 1, 1/a] });
const jt11 = m => (
    { z: [...range(1-m, m-1).map(k => Exp( 1, k*π/(2*m))),
          ...range(m, 3*m)  .map(k => Exp(.9, k*π/(2*m)))] });
const uhlig = a => ({ z: [[a, 5], -a, C(0, a), C(0, -a)] });


const TEST_POLYS = {
    // Jenkins--Traub test polynomials
    jt01a: jt01(10**10),
    jt01b: jt01(10**-10),
    jt02: { z: range(1, 17) },
    jt03: { z: range(1, 8).map(x => 1/10**x) },
    jt04: { z: [0.5, 0.6, 0.7, [0.1, 3]] },
    jt05: { z: [0.4, [0.3, 2], [0.2, 3], [0.1, 4]] },
    jt06: { z: [.1,1.001,.998,1.00002,.99999] },
    jt07a: jt07(10**-10),
    jt07b: jt07(10**-9),
    jt07c: jt07(10**-8),
    jt07d: jt07(10**-7),
    jt08: { z: [-1, 5] },
    jt09: {
        p: ([1, ...constant(0, 9), -(10**-20)],
                    [1, ...constant(0, 9), 10**20]),
        z: [...range(0, 9).map(k => Exp(.01, 2*π*k/10)),
            ...range(0, 9).map(k => Exp(100, π*(2*k+1)/10))]
    },
    jt10a: jt10(10**3),
    jt10b: jt10(10**6),
    jt10c: jt10(10**9),
    jt11a: jt11(15),
    jt11b: jt11(20),
    jt11c: jt11(25),

    // Uhlig test polynomials
    uhlig01: uhlig(.01),
    uhlig02: uhlig(.001),
    uhlig05: { z: [...constant(-1, 6), 2, 2] },

    // Goedecker test polynomials
};

