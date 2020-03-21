'use strict';

import { mat, vec, cardano } from './linalg.js';

let test_mats = [
    [
        [10,-7,0],
        [-3, 2,6],
        [ 5,-1,5]
    ], [
        [2, 1, 1, 0],
        [4, 3, 3, 1],
        [8, 7, 9, 5],
        [6, 7, 9, 8]
    ], [
        [-1, 0, 1],
        [ 2, 1, 1],
        [-1, 2, 0]
    ], [
        [ 2, -6,  6],
        [-4,  5, -7],
        [ 3,  5, -1],
        [-6,  4, -8],
        [ 8, -3,  9]
    ]
];

// for(let rows of test_mats) {
//     let A = mat.from_array(rows);
//     let {P, L, U} = A.PLU;
//     console.log(P.toString());
//     console.log(L.toString());
//     console.log(U.toString());
//     console.log(P.mult(A).sub(L.mult(U)).toString());
// }




function test_cardano() {

    let data = [
        [-6, 12, -8], // x^3 - 6x^2 + 12x - 8 = (x-2)^3
        [-4,  5, -2], // x^3 - 4x^2 +  5x - 2 = (x-2) (x-1)^2
        [ 4,  5,  2], // x^3 + 4x^2 +  5x + 2 = (x+2) (x+1)^2
        [ 2, -1, -2], // x^3 + 2x^2 -   x - 2 = (x+2) (x+1) (x-1)
        [ 1, -4, -4], // x^3 +  x^2 -  4x - 4 = (x+2) (x+1) (x-2)
        [-2,  1, -2], // x^2 - 2x^2 +   x - 2 = (x-2) (x^2+1)
    ];

    function polyeval(x, coeffs) {
        let acc = 0;
        let pwr = 1;
        for(let c of coeffs.reverse()) {
            acc += c*pwr;
            pwr *= x;
        }
        return acc + pwr;
    }

    function test(b, c, d) {
        let roots = cardano(b, c, d);
        let str = `f(x) = x^3 + ${b}x^2 + ${c}x + ${d} has `;
        let results = [], strs = [];
        for(let [root, mult] of roots) {
            let type = ['simple', 'double', 'triple'][mult-1];
            strs.push(`a ${type} root at ${root.toFixed(4)}`);
            let res = polyeval(root, [b, c, d]);
            results.push(`f(${root.toFixed(4)}) = ${res.toFixed(4)}`);
        }
        console.log(str + strs.join(' and '));
        console.log(results.join('\n'));
    }

    for(let coeffs of data)
        test(...coeffs);
}

// test_cardano();
