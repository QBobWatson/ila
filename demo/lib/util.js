'use strict'; // -*- js2 -*-

/** @module lib/matrix
 *
 * @file
 * Utility functions.
 */

/**
 * @summary
 * Iterable over a range of values.
 */
export function* range(a, b=undefined, step=1) {
    let min, max;
    if(b === undefined)
        [min, max] = [0, a];
    else
        [min, max] = [a, b];
    if(step < 0)
        [min, max] = [max, min];
    for(let i = min; i < max; i += step)
        yield i;
}
