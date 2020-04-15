'use strict'; // -*- js2 -*-

export function* range(a, b=undefined, step=1, map=i=>i) {
    let min, max;
    if(b === undefined)
        [min, max] = [0, a];
    else
        [min, max] = [a, b];
    if(step > 0) {
        for(let i = min; i < max; i += step)
            yield map(i);
    } else {
        for(let i = min; i > max; i += step)
            yield map(i);
    }
}

export function constant(x, count) {
    return range(0, count, 1, ()=>x);
}

