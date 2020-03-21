'use strict';

// Polynomial root solvers

export function quadratic(b, c, ε=1e-5) {
    let Δ = b*b - 4*c;
    if(Math.abs(Δ) < ε)
        return [[-b/2, 2]];
    if(Δ < 0)
        return [];
    let D = Math.sqrt(Δ);
    return [[(-b+D)/2, 1], [(-b-D)/2, 1]];
}

// Compute roots of x^3 + bx^2 + cx + d
export function cardano(b, c, d, ε=1e-5) {
    // Change of variables x --> x-b/3
    let [p, q] = [-1/3*b*b+c, 2/27*b*b*b - 1/3*b*c + d];
    // Discriminant
    let Δ = -27*q*q - 4*p*p*p;
    let ret;
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
        ret = [[mod * Math.cos(arg/3              ), 1],
               [mod * Math.cos(arg/3 + 2*Math.PI/3), 1],
               [mod * Math.cos(arg/3 - 2*Math.PI/3), 1]
              ];
    } else {
        // Simple real root and conjugate complex roots (ignored)
        let D = Math.sqrt(-Δ/108);
        ret = [[Math.cbrt(-q/2 + D) + Math.cbrt(-q/2 - D), 1]];
    }
    ret.sort((a, b) => a[0] - b[0]).forEach(a => a[0] -= b/3);
    return ret;
}


////////////////////////////////////////////////////////////////////////////////
// * Vector

export class vec {
    static constant(n, c) {
        return new vec(new Array(n).fill(c));
    }
    static zero(n) {
        return vec.constant(n, 0);
    }

    constructor(vals) {
        this.vals = vals;
    }

    clone() {
        return new vec(this.vals.slice());
    }

    toString(precision=4) {
        return "[" + this.vals.map(v => v.toFixed(precision)).join(' ') + "]";
    }

    get size() {
        return this.vals.length;
    }

    get lengthsq() {
        return this.dot(this);
    }
    get length() {
        return Math.sqrt(this.lengthsq);
    }

    normalize() {
        return this.scale(1/this.length);
    }

    // In-place addition and subtraction
    add(other, factor=1, start=0) {
        for(let i = start; i < this.size; ++i)
            this.vals[i] += factor*other.vals[i];
        return this;
    }
    sub(other, start=0) {
        return this.add(other, -1, start);
    }
    scale(c, start=0) {
        for(let i = start; i < this.size; ++i)
            this.vals[i] *= c;
        return this;
    }

    dot(other) {
        return this.vals.reduce((a, v, i) => a + v * other.vals[i], 0);
    }
}


////////////////////////////////////////////////////////////////////////////////
// * Matrix

// This is a simple matrix object that supports matrix factorizations
// Its rows are vec instances.

export class mat {
    static identity(n, λ=1) {
        let rows = new Array(n);
        for(let i = 0; i < n; ++i) {
            rows[i] = vec.zero(n);
            rows[i].vals[i] = λ;
        }
        return new mat(rows);
    }

    static fromArray(arr, m, n) {
        if(m === undefined || n === undefined)
            return new mat(arr.map(row => new vec(row)));
        let rows = [];
        for(let i = 0; i < m; ++i)
            rows.push(new vec(arr.slice(i * n, (i+1) * n)));
        return new mat(rows);
    }

    constructor(rows) {
        this.rows = rows;
        this.m = this.rows.length;
        this.n = this.rows[0].size;

        // Threshold for deciding when a number is zero (settable)
        this.ε = 1e-5;

        // Gaussian elimination
        this._PLU = null;
        this._rref = null;

        // Bases
        this._nullBasis = null;

        // Subspaces
        this._rowSpace = null;
        this._columnSpace = null;
        this._nullSpace = null;
        this._leftNullSpace = null;

        // Eigendata
        this._det = null;
        this._eigendata = null;
    }

    clone() {
        return new mat(this.rows.map(row => row.clone()));
    }

    get PLU() {
        if(!this._PLU) this._calcPLU();
        return this._PLU;
    }

    get pivots() {
        let {pivots} = this.PLU;
        return pivots;
    }

    get rank() {
        return this.pivots.length;
    }

    get nullity() {
        return this.n - this.rank;
    }

    get trace() {
        let acc = 0;
        for(let i = 0; i < this.n; ++i)
            acc += this.entry(i, i);
        return acc;
    }

    get det() {
        if(this._det === null) this._calcDet();
        return this._det;
    }

    get transpose() {
        let rows = new Array(this.n);
        for(let j = 0; j < this.n; ++j)
            rows[j] = new vec(this.col(j));
        return new mat(rows);
    }

    get rref() {
        if(!this._rref) this._calcRREF();
        return this._rref;
    }

    // Four subspaces (TODO)
    get rowSpace() {
    }

    get columnSpace() {
    }

    get nullBasis() {
        if(!this._nullBasis) this._calcNullBasis();
        return this._nullBasis;
    }
    get nullSpace() {
        if(!this._nullSpace)
            this._nullSpace = Subspace.span(this.nullBasis);
        return this._nullSpace;
    }

    get leftNullSpace() {
    }

    // Eigenspaces
    get eigenvalues() {
        if(!this._eigendata) this._calcEigendata();
        return this._eigendata.eigenvalues;
    }

    get charpoly() {
        if(!this._eigendata) this._calcEigendata();
        return this._eigendata.charpoly;
    }

    toString(precision=4) {
        let strings = this.rows.map(
            row => row.vals.map(
                val => val.toFixed(precision)
            )
        );
        let colLens = new Array(this.n);
        for(let j = 0; j < this.n; ++j)
            colLens[j] = Math.max(...strings.map(row => row[j].length));
        let ret = [];
        for(let row of strings) {
            ret.push(row.map((val, j) => val.padStart(colLens[j], ' ')).join(' '));
        }
        return ret.join('\n');
    }

    // In-place addition and subtraction
    add(other, factor=1) {
        this.rows.forEach((row, i) => row.add(other.rows[i], factor));
        return this;
    }
    sub(other) {
        return this.add(other, -1);
    }

    // Multiplication
    mult(other) {
        return new mat(this.rows.map(
            row => {
                let newRow = new Array(other.n);
                for(let i = 0; i < other.n; ++i) {
                    let accum = 0;
                    for(let j = 0; j < this.n; ++j)
                        accum += row.vals[j] * other.entry(j, i);
                    newRow[i] = accum;
                }
                return new vec(newRow);
            }
        ));
    }

    apply(v) {
        return new vec(this.rows.map(row => row.dot(v)));
    }

    row(i) {
        return this.rows[i];
    }

    col(j) {
        return this.rows.map(x => x.vals[j]);
    }

    entry(i, j) {
        return this.rows[i].vals[j];
    }
    setEntry(i, j, c) {
        this.rows[i].vals[j] = c;
        return this;
    }

    // Row operations

    // row[i] *= c
    rowScale(i, c, start=0) {
        this.rows[i].scale(c, start);
        return this;
    }
    // row[i1] += c*row[i2]
    rowReplace(i1, i2, c, start=0) {
        this.rows[i1].add(this.rows[i2], c, start);
        return this;
    }
    // row[i1] <-> row[i2]
    rowSwap(i1, i2) {
        [this.rows[i1], this.rows[i2]] = [this.rows[i2], this.rows[i1]];
        return this;
    }

    // PLU factorization
    _calcPLU() {
        let P = mat.identity(this.m);
        let L = mat.identity(this.m);
        let U = this.clone();
        let {m, n, ε} = this;
        let pivots = [];

        for(let curRow = 0, curCol = 0; curRow < m && curCol < n; ++curCol) {
            // Maximal pivot
            let pivot = U.entry(curRow, curCol), row = curRow;
            for(let i = curRow+1; i < m; ++i) {
                if(Math.abs(U.entry(i, curCol)) > Math.abs(pivot)) {
                    pivot = U.entry(i, curCol);
                    row = i;
                }
            }

            if(Math.abs(pivot) > ε) {
                if(row != curRow) {
                    // Swap
                    P.rowSwap(row, curRow);
                    U.rowSwap(row, curRow);
                    for(let j = 0; j < curCol; ++j)
                        [L.rows[row].vals[j], L.rows[curRow].vals[j]]
                            = [L.rows[curRow].vals[j], L.rows[row].vals[j]];
                }
                // Eliminate
                for(let i = curRow+1; i < m; ++i) {
                    let l = U.entry(i, curCol) / pivot;
                    L.setEntry(i, curRow, l);
                    U.rowReplace(i, curRow, -l, curCol);
                }
                pivots.push([curRow, curCol]);
                curRow++;

            } else {
                // Clear the column so U is really upper-triangular
                for(let i = curRow; i < m; ++i)
                    U.setEntry(i, curCol, 0);
            }
        }

        this._PLU = {P, L, U, pivots};
    }

    // Gaussian elimination
    _calcRREF() {
        let {U, pivots} = this.PLU;
        let rref = U.clone();
        for(let k = pivots.length-1; k >= 0; --k) {
            let [row, col] = pivots[k];
            let pivot = rref.entry(row, col);
            for(let i = 0; i < row; ++i) {
                rref.rowReplace(i, row, -rref.entry(i, col)/pivot, col+1);
                rref.setEntry(i, col, 0);
            }
            rref.rowScale(row, 1/pivot);
        }
        this._rref = rref;
    }

    // Bases
    _calcNullBasis() {
        // Parametric vector form
        let rref = this.rref, pivots = this.pivots, previous = [];
        let basis = [];
        for(let j = 0; j < this.n; ++j) {
            if(pivots.length && pivots[0][1] == j) {
                // Pivot column
                previous.push(pivots.shift());
                continue;
            }
            // Free column
            let v = new Array(this.n).fill(0);
            for(let [row, col] of previous)
                v[col] = -rref.entry(row, j);
            v[j] = 1;
            basis.push(new vec(v));
        }
        this._nullBasis = basis;
    }

    // Eigenvalues and eigenspaces.  Only implemented for 1x1, 2x2, 3x3.
    _calcDet() {
        if(this.m != this.n)
            throw "Determinants only make sense for square matrices";
        switch(this.n) {
        case 0:
            this._det = 1;
            break;
        case 1:
            this._det = this.entry(0,0);
            break;
        case 2: {
            let [a,b] = this.rows[0].vals;
            let [c,d] = this.rows[1].vals;
            this._det = a*d - b*c;
        } break;
        case 3: {
            let [a,b,c] = this.rows[0].vals;
            let [d,e,f] = this.rows[1].vals;
            let [g,h,i] = this.rows[2].vals;
            this._det = a*(e*i-f*h)-b*(d*i-f*g)+c*(d*h-e*g);
        } break;
        default:
            throw "Eigenvalue calculations only implemented for 1 <= n <= 3";
        }
    }

    _calcEigendata() {
        if(this.m != this.n)
            throw "Eigenvalues only make sense for square matrices";
        switch(this.n) {
        case 1: this._calcEigendata1x1(); break;
        case 2: this._calcEigendata2x2(); break;
        case 3: this._calcEigendata3x3(); break;
        default:
            throw "Eigenvalue calculations only implemented for 1 <= n <= 3";
        }
        this._eigendata.eigenvalues.forEach(val => {
            val.eigenspace = this.clone()
                .sub(mat.identity(this.n, val.value))
                .nullSpace;
            // Sanity check
            if(val.eigenspace.dim == 0 || val.eigenspace.dim > val.mult)
                throw "Numerical error computing eigenspaces";
        });
    }

    _calcEigendata1x1() {
        let eigenvalue = this.entry(0, 0);
        this._eigendata = {
            charpoly: [eigenvalue],
            eigenvalues: [{
                value: eigenvalue,
                mult:  1,
            }]
        };
    }

    _calcEigendata2x2() {
        let [b, c] = [-this.trace, this.det];
        let roots = quadratic(b, c, this.ε);
        this._eigendata = {
            charpoly: [b, c],
            eigenvalues: roots.map(([root, mult]) => {
                return {value: root, mult: mult};
            })
        };
    }

    _calcEigendata3x3() {
        let [b, c, d] = [this.trace, (() => {
            let [a,b,c] = this.rows[0].vals;
            let [d,e,f] = this.rows[1].vals;
            let [g,h,i] = this.rows[2].vals;
            return b*d - a*e + c*g + f*h - a*i - e*i;
        })(), this.det];
        let roots = cardano(b, -c, d); // replace x by -x
        this._eigendata = {
            charpoly: [b, c, d],
            eigenvalues: roots.map((root, mult) => {
                return {value: -root, mult: mult};
            })
        };
    }
}


////////////////////////////////////////////////////////////////////////////////
// * Subspace

// Abstract representation of a subspace.
// Internally it is stored as an orthonormal basis.

export class Subspace {
    // Gram--Schmidt
    static span(vecs, n=vecs[0].size, ε=1e-5) {
        let ortho = [];

        for(let v of vecs) {
            let comp = new Subspace(ortho, n).orthogonal(v);
            let length = comp.lengthsq;
            if(length < ε*ε)
                continue;
            length = Math.sqrt(length);
            comp.scale(1/length);
            ortho.push(comp);
        }

        return new Subspace(ortho, n);
    }

    // Construct given an *orthonormal* basis.
    constructor(basis, n=basis[0].size) {
        this.basis = basis;
        this.n = n;
    }

    toString() {
        let str = `Subspace of R^${this.n} of dimension ${this.dim} with basis\n`;
        str += new mat(this.basis).transpose.toString();
        return str;
    }

    get dim() {
        return this.basis.length;
    }

    // Orthogonal projection
    project(v) {
        let coeffs = this.basis.map(b => v.dot(b));
        return this.basis.reduce(
            (ret, val, idx) => ret.add(val, coeffs[idx]),
            vec.zero(this.n)
        );
    }

    // Orthogonal complement
    orthogonal(v) {
        return v.sub(this.project(v));
    }
}
