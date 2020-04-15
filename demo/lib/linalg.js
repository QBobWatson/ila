/** @module lib/linalg
 *
 * @file
 * Classes and routines for doing linear algebra.
 */

'use strict';

// TODO:
//  * SVD
//  * pseudo-inverse
//  * chednovsky or whatever
//  * PCA?

import Vector from "./vector.js";
import Complex from "./complex.js";
import Matrix from "./matrix.js";

/** @type {number} */
const π = Math.PI;


/**
 * Roots are counted with multiplicity: `[x, m]` indicates a root of a
 * polynomial at the number `x` with multiplicity `m`
 *
 * @typedef {number[]} Root
 */

/**
 * Find the roots of a quadratic polynomial.
 *
 * This function finds the roots of the polynomial `x^2 + b x + c` using the
 * quadratic formula.
 *
 * @example
 * quadratic(0, -1); // returns [[1, 1], [-1, 1]] (approximately)
 *
 * @param {number} b - The linear coefficient.
 * @param {number} c - The constant coefficient.
 * @param {number} [ε=1e-10] - The roots will be considered equal if the
 *   discriminant `b^2-4c` is smaller than this.
 * @return {Root[]}
 */
export function quadratic(b, c, ε=1e-10) {
    let Δ = b*b - 4*c;
    if(Math.abs(Δ) < ε)
        return [[-b/2, 2]];
    if(Δ < 0) {
        let D = Math.sqrt(-Δ);
        let z = new Complex(-b/2, D/2);
        return [[z, 1], [z.clone().conj(), 1]];
    }
    let D = Math.sqrt(Δ);
    return [[(-b-D)/2, 1], [(-b+D)/2, 1]];
}

/**
 * A cube root of unity.
 *
 * @type {Complex}
 */
const ζ = new Complex(-1/2, Math.sqrt(3)/2);

/**
 * Find the roots of a cubic polynomial.
 *
 * This function finds the roots of the polynomial `x^2 + b x^2 + c x + d` using
 * [Cardano's formula]{@link https://www.encyclopediaofmath.org/index.php/Cardano_formula}.
 *
 * @example
 * cardano(-4, 5, -2); // returns [[1, 2], [2, 1]] (approximately)
 *
 * @param {number} b - The quadratic coefficient.
 * @param {number} c - The linear coefficient.
 * @param {number} d - The constant coefficient.
 * @param {number} [ε=1e-10] - The roots will be considered equal if the
 *   discriminant is smaller than this.
 * @return {Root[]}
 */
export function cardano(b, c, d, ε=1e-10) {
    // Change of variables x --> x-b/3
    let [p, q] = [-1/3*b*b+c, 2/27*b*b*b - 1/3*b*c + d];
    // Discriminant
    let Δ = -27*q*q - 4*p*p*p;
    let ret, cplx = [];
    if(Math.abs(Δ) < ε) {
        if(Math.abs(p) < ε && Math.abs(q) < ε)
            return [[-b/3, 3]]; // Triple root
        // Simple root and double root
        let cr = Math.cbrt(-q/2);
        ret = [[2*cr, 1], [-cr, 2]];
    } else if(Δ > 0) {
        // Three distinct real roots: 2*Re(cube roots of -q/2 + i sqrt(Δ/108))
        let D = Math.sqrt(Δ/108);
        let mod = Math.sqrt(Math.cbrt(q*q/4 + Δ/108)) * 2;
        let arg = Math.atan2(D, -q/2);
        ret = [[mod * Math.cos(arg/3        ), 1],
               [mod * Math.cos(arg/3 + 2*π/3), 1],
               [mod * Math.cos(arg/3 - 2*π/3), 1]
              ];
    } else {
        // Simple real root and conjugate complex roots
        let D = Math.sqrt(-Δ/108);
        let α = Math.cbrt(-q/2 + D), β = Math.cbrt(-q/2 - D);
        ret = [[α + β, 1]];
        let z = ζ.clone().scale(α).add(ζ.clone().conj().scale(β)).sub(b/3);
        if(z.Im > 0) z.conj();
        cplx = [[z, 1], [z.clone().conj(), 1]];
    }
    ret.sort((a, b) => a[0] - b[0]).forEach(a => a[0] -= b/3);
    return ret.concat(cplx);
}


/**
 * Find the roots of a polynomial.
 *
 * This computes the roots of the polynomial
 *    `x^n + a1 x^{n-1} + ... + an`
 * using a numerical algorithm.
 *
 * @param {number[]} coeffs - The coefficients `a1, ..., an`.
 * @param {number} [ε=1e-10] - Roots will be considered equal if their
 *   difference is smaller than this.
 * @return {Root[]}
 */
export function roots(coeffs, ε=1e-10) {
    let n = coeffs.length;
    let degpar = {Degree: n};
    let zeror = new Array(n).fill(0);
    let zeroi = new Array(n).fill(0);
    let coeffs2 = [1, ...coeffs];
    let ret = rpSolve(degpar, coeffs2, zeror, zeroi);
    return [zeror, zeroi];
}


