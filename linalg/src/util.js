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

'use strict'; // -*- js2 -*-

/** @module util
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
