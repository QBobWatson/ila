'use strict'; // -*- js2 -*-

/* Generate random test polynomials for `rpoly`.
 *
 * Usage: `node randpolys.js NUM_POLYS OUT_FILE`
 */

import fs from 'fs';

const { random, floor, min, sin, cos, PI } = Math;

import Polynomial from '../lib/polynomial.js';

const MAX_DEG = 50;
const NUM_POLYS = parseInt(process.argv[2]);
const OUT_FILE = process.argv[3];

console.log(`Generating ${NUM_POLYS} polynomials and saving to ${OUT_FILE}...`);

/* Generate a random multiplicity.
 *
 * The multiplicity has an 80% chance of being 1, a 13% chance of being 2, and a
 * 7% chance of being 3.
 */
function multiplicity() {
    let r = random();
    if(r < .8)  return 1;
    if(r < .93) return 2;
    return 3;
}

const len = (NUM_POLYS-1).toString().length;

let ws = fs.createWriteStream(OUT_FILE);

for(let i = 0; i < NUM_POLYS; ++i) {
    // Random degree between 1 and MAX_DEG
    let deg = floor(random() * MAX_DEG) + 1;
    // Chance of generating real roots
    let pctReal = 0.25 + .75 * random();
    // Roots will have absolute value less than `radius`
    let radius = 0.3 + random() * (30 / deg);

    let P = Polynomial.create(1);
    let roots = [];
    let name = `rand${i.toString().padStart(len, '0')}`;

    while(P.deg < deg) {
        if(random() < pctReal || deg - P.deg == 1) {
            let mult = min(multiplicity(), deg - P.deg);
            // Generate real root
            let root = radius * (2 * random() - 1);
            for(let i = 0; i < mult; ++i) {
                roots.push([root, 0]);
                P = P.mult(Polynomial.create(1, -root));
            }
        } else {
            let mult = floor(min(multiplicity(), (deg - P.deg)/2));
            // Generate complex root
            let arg = random() * PI;
            let mod = random() * radius;
            let re = mod * cos(arg);
            let im = mod * sin(arg);
            for(let i = 0; i < mult; ++i) {
                roots.push([re, im], [re, -im]);
                P = P.mult(Polynomial.create(1, -2*re, mod*mod));
            }
        }
    }

    ws.write(`${name}\n`);
    ws.write(`${deg}\n`);
    for(let coeff of P)
        ws.write(`${coeff}\n`);
    for(let [re, im] of roots)
        ws.write(`${re} ${im}\n`);
}

ws.close();

console.log("Finished.");
