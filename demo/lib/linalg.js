/** @module lib/linalg
 *
 * @file
 * Classes and routines for doing linear algebra.  Focuses on matrix
 * factorizations, with implementations in dimensions 2 and 3.
 */

'use strict';

/** @type {number} */
const π = Math.PI;

/**
 * Class representing a vector.
 *
 * A vector is a sequence of numbers of a fixed length, called the "length" of
 * the vector.
 *
 * Note that `Vector` uses the `Array` constructor, so `new Vector(3)` will
 * return an "empty" vector of length 3.  Instead use `Vector.constant(1, 3)`.
 *
 * @extends Array
 */
export class Vector extends Array {
    /**
     * Create a Vector with all entries equal to `c`.
     *
     * @param {integer} n - The size of the resulting Vector.
     * @param {number} c - The value of the entries.
     * @return {Vector}
     */
    static constant(n, c) {
        let ret = new Vector(n);
        ret.fill(c);
        return ret;
    }

    /**
     * Create a Vector with all entries equal to 0.
     *
     * @param {integer} n - The size of the resulting Vector.
     * @return {Vector}
     */
    static zero(n) {
        return Vector.constant(n, 0);
    }

    /**
     * Create a Vector with all entries equal to 0, except the ith equals 1.
     *
     * @param {integer} i - The nonzero entry.
     * @param {integer} n - The size of the resulting Vector.
     * @param {number} [λ=1] - The nonzero entry is set to this.
     * @return {Vector}
     */
    static e(i, n, λ=1) {
        let ret = Vector.zero(n);
        ret[i] = λ;
        return ret;
    }

    /**
     * Check if an iterable of vectors is linearly independent.
     *
     * This means that the matrix with columns equal to `vecs` has zero null
     * space.
     *
     * @param {Vector[]} vecs - The vectors to check.
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting.
     * @return {boolean}
     * @throws Will throw an error if the vectors do not have the same length.
     */
    static isLinearlyIndependent(vecs, ε=1e-10) {
        let M = new Matrix(...vecs);
        // Tall matrices never have linearly independent rows.
        if(M.m > M.n) return false;
        M.PLU(ε);
        return M.rank == vecs.length;
    }

    /**
     * Check if an iterable of vectors is linearly dependent.
     *
     * This is an alias for `!Vector.isLinearlyIndependent(vecs)`.
     *
     * @param {Vector[]} vecs - The vectors to check.
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting.
     * @return {boolean}
     * @throws Will throw an error if the vectors do not have the same length.
     */
    static isLinearlyDependent(vecs, ε=1e-10) {
        return !Vector.isLinearlyIndependent(vecs, ε);
    }

    /**
     * Return a linearly independent subset.
     *
     * This deletes some number of vectors from `vecs` until the resulting
     * subset is linearly independent.
     *
     * @param {Vector[]} vecs - The vectors to use.
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting.
     * @return {Vector[]}
     * @throws Will throw an error if the vectors do not have the same length.
     */
    static linearlyIndependentSubset(vecs, ε=1e-10) {
        return new Matrix(...vecs).transpose().colBasis();
    }

    /**
     * Compute a linear combination of vectors.
     *
     * @param {number[]} coeffs - A non-empty array of coefficients.
     * @param {number[]} vecs - A non-empty array of vectors, of the same length
     *   as `coeffs`.
     * @return The sum of the vectors scaled by the coefficients.
     * @throws Will throw an error if the vectors do not have the same length,
     *   or if `vecCoeffs` is empty.
     */
    static linearCombination(coeffs, vecs) {
        let accum = Vector.zero(vecs[0].length);
        let i = 0;
        for(let coeff of coeffs)
            accum.add(vecs[i++].scale(coeff));
        return accum;
    }

    /**
     * Check if this vector is equal to `other`.
     *
     * Two vectors are equal if they have the same number of entries, and all
     * entries are equal.
     *
     * @param {Vector} other - The vector to compare.
     * @param {number} [ε=0] - Entries will test as equal if they are within `ε`
     *   of each other.  This is provided in order to account for rounding
     *   errors.
     * @return {boolean}
     */
    equals(other, ε=0) {
        if(this.length !== other.length)
            return false;
        if(ε === 0)
            return this.every((v, i) => v === other[i]);
        return this.every((v, i) => Math.abs(v - other[i]) < ε);
    }

    /**
     * Create a new Vector with the same entries.
     *
     * @return {Vector}
     */
    clone() {
        return this.slice();
    }

    /**
     * Return a string representation of the vector.
     *
     * @param {integer} [precision=4] - The number of decimal places to include.
     * @return {string}
     */
    toString(precision=4) {
        // this.map() returns a new Vector
        let strings = new Array(this.length), i = 0;
        for(let v of this)
            strings[i++] = v.toFixed(precision);
        return "[" + strings.join(' ') + "]";
    }

    /**
     * Decide if a vector is zero.
     *
     * This means that all entries are at most `ε`.
     *
     * @param {number} [ε=0] - Entries smaller than this will be considered
     *   zero.
     * @return {boolean}
     */
    isZero(ε=0) {
        return this.every(x => Math.abs(x) <= ε);
    }


    /**
     * The square of the geometric length of the vector.
     *
     * @example
     * let v = new Vector([1, 2, 3]);
     * v.sizesq == 1*1 + 2*2 + 3*3
     *
     * @type {number}
     */
    get sizesq() {
        return this.dot(this);
    }

    /**
     * The geometric length of the vector.
     *
     * This is the square root of `this.sizesq`.
     *
     * @example
     * let v = new Vector([1, 2, 3]);
     * v.size == Math.sqrt(1*1 + 2*2 + 3*3)
     *
     * @type {number}
     */
    get size() {
        return Math.sqrt(this.sizesq);
    }

    /**
     * Scale the vector by the reciprocal of its length.
     *
     * This modifies the vector in-place to have length 1.
     *
     * @example
     * let v = new Vector([3, 4]);
     * v.normalize();
     * v.equals(new Vector([3/5, 4/5])) == true;
     *
     * @return {Vector} `this`
     * @throws Will throw an error if `this` is the zero vector.
     */
    normalize() {
        let s = 1/this.size;
        if(!isFinite(s))
            throw new Error("Tried to normalize the zero vector");
        return this.scale(s);
    }

    /**
     * Add a Vector in-place.
     *
     * This modifies the vector in-place by adding the entries of `other`.
     *
     * @example
     * let v = new Vector([1, 2]), w = new Vector([3, 4]);
     * v.add(w);
     * v.equals(new Vector([4, 6])) == true;
     *
     * @param {Vector} other - The vector to add.
     * @param {number} [factor=1] - Add `factor` times `other` instead of just
     *   adding `other`.
     * @param {integer} [start=0] - Only add the entries `start...this.length`.
     *   Provided for optimizations when the entries of `other` before `start`
     *   are known to be zero.
     * @return {Vector} `this`
     * @throws Will throw an error if the vectors have different lengths.
     */
    add(other, factor=1, start=0) {
        if(this.length !== other.length)
            throw new Error(
                'Tried to add vectors of different lengths');
        for(let i = start; i < this.length; ++i)
            this[i] += factor*other[i];
        return this;
    }

    /**
     * Subtract a Vector in-place.
     *
     * This modifies the vector in-place by subtracting the entries of `other`.
     *
     * @example
     * let v = new Vector([1, 2]), w = new Vector([3, 4]);
     * v.sub(w);
     * v.equals(new Vector([-2, -2])) == true;
     *
     * @param {Vector} other - The vector to subtract.
     * @param {integer} [start=0] - Only subtract the entries
     *   `start...this.length`.  Provided for optimizations when the entries of
     *   `other` before `start` are known to be zero.
     * @return {Vector} `this`
     * @throws Will throw an error if the vectors have different lengths.
     */
    sub(other, start=0) {
        return this.add(other, -1, start);
    }

    /**
     * Multiply a Vector by a scalar in-place.
     *
     * This modifies the vector in-place by multiplying all entries by `c`.
     *
     * @example
     * let v = new Vector([1, 2]);
     * v.scale(2);
     * v.equals(new Vector([2, 4])) == true;
     *
     * @param {number} c - The scaling factor.
     * @param {integer} [start=0] - Only scale the entries
     *   `start...this.length`.  Provided for optimizations when the entries
     *   before `start` are known to be zero.
     * @return {Vector} `this`
     */
    scale(c, start=0) {
        for(let i = start; i < this.length; ++i)
            this[i] *= c;
        return this;
    }

    /**
     * Compute the dot product with another vector.
     *
     * The other vector must have the same length as `this`.
     *
     * @example
     * let v = new Vector([1, 2]), w = new Vector([3, 4]);
     * v.dot(w) == 1*3 + 2*4;
     *
     * @param {Vector} other - The vector to dot.
     * @return {number} The dot product.
     * @throws Will throw an error if the vectors have different lengths.
     */
    dot(other) {
        if(this.length !== other.length)
            throw new Error(
                'Tried to take the dot product of vectors of different lengths');
        return this.reduce((a, v, i) => a + v * other[i], 0);
    }
}


/**
 * Class representing a complex number.
 *
 * Complex numbers behave like vectors of length two, in that they can be added
 * and scaled componentwise.  However, complex numbers can also be multiplied
 * together and divided to form new complex numbers.
 *
 * @extends Vector
 */
export class Complex extends Vector {
    /**
     * Create a new Complex number with prescribed polar coordinates.
     *
     * Cartesian and polar coordinates are related by the formula
     * `(x, y) = (r cos(θ), r sin(θ))`.
     *
     * @example
     * Complex.fromPolar(1, Math.PI/2).equals(Complex.i, 1e-10);
     *
     * @param {number} r - The radial coordinate.
     * @param {number} θ - The angular coordinate.
     * @return {Complex}
     */
    static fromPolar(r, θ) {
        return new Complex(r * Math.cos(θ), r * Math.sin(θ));
    }

    /**
     * A copy of `i`, a square root of -1.
     *
     * This returns a new object each time it is accessed.
     *
     * @type {Complex}
     */
    static get i() {
        return new Complex(0, 1);
    }

    /**
     * Create a complex number.
     *
     * Overrides the Vector constructor to require exactly two values.
     *
     * @param {number} a - The real part.
     * @param {number} [b=0] - The imaginary part.
     */
    constructor(a, b=0) {
        super(a, b);
    }

    /**
     * Check if this complex number is equal to `other`.
     *
     * This is the same as {@link Vector#equals}, except that if `other` is a
     * number, it is promoted to a Complex number first.
     *
     * @example
     * new Complex(1, 0).equals(1) == true;
     *
     * @param {(Complex|number)} other - The number to compare.
     * @param {number} [ε=0] - Entries will test as equal if they are within `ε`
     *   of each other.  This is provided in order to account for rounding
     *   errors.
     * @return {boolean}
     */
    equals(other, ε=0) {
        if(typeof other === "number")
            other = new Complex(other, 0);
        return super.equals(other, ε);
    }

    /**
     * Create a new Complex number with the same entries.
     *
     * @return {Complex}
     */
    clone() {
        return new Complex(this.Re, this.Im);
    }

    /**
     * Return a string representation of the complex number.
     *
     * @param {integer} [precision=4] - The number of decimal places to include.
     * @return {string}
     */
    toString(precision=4) {
        return `${this.Re.toFixed(precision)} + ${this.Im.toFixed(precision)} i`;
    }

    /**
     * The real part (first entry) of the complex number.
     *
     * @example
     * new Complex(3, 4).Re == 3;
     *
     * @type {number}
     */
    get Re() { return this[0]; }
    set Re(x) { this[0] = x; }

    /**
     * The imaginary part part (second entry) of the complex number.
     *
     * @example
     * new Complex(3, 4).Re == 4;
     *
     * @type {number}
     */
    get Im() { return this[1]; }
    set Im(x) { this[1] = x; }

    /**
     * The modulus of the complex number.
     *
     * This is an alias for `this.size`.
     *
     * @example
     * new Complex(3, 4).mod == 5;
     *
     * @type {number}
     */
    get mod() { return this.size; }

    /**
     * The argument of the complex number.
     *
     * The value is between -π and π.
     *
     * @type {number}
     */
    get arg() { return Math.atan2(this.Im, this.Re); }

    /**
     * Replace the Complex number with its complex conjugate.
     *
     * This negates the imaginary part.
     *
     * @example
     * new Complex(1, 2).conj().equals(new Complex(1, -2));
     *
     * @return {Complex} `this`
     */
    conj() {
        this.Im *= -1;
        return this;
    }

    /**
     * Add another complex number in-place.
     *
     * This is the same as {@link Vector#add}, except that if `other` is a
     * number, it is promoted to a Complex number first.
     *
     * @example
     * new Complex(1, 2).add(1).equals(new Complex(2, 2)) == true;
     *
     * @param {(Complex|number)} other - The number to add.
     * @param {number} [factor=1] - Add `factor` times `other` instead of just
     *   adding `other`.
     * @return {Complex} `this`
     */
    add(other, factor=1) {
        if(typeof other === "number") {
            this.Re += other * factor;
            return this;
        }
        return super.add(other, factor);
    }

    /**
     * Multiply by a complex number in-place.
     *
     * Multiplication of complex numbers is defined by the formula
     * `(a + b i) (c + d i) = (ac - bd) + (ad + bc) i`.
     *
     * @example
     * let z = new Complex(1, 2), w = new Complex(3, 4);
     * z.mult(w).equals(new Complex(-5, 10)) == true;
     *
     * @param {(Complex|number)} other - The number to multiply.
     * @return {Complex} `this`
     * @throws Will throw an error if `other` is zero.
     */
    mult(other) {
        if(typeof other === "number")
            return this.scale(other);
        [this.Re, this.Im] = [
            this.Re * other.Re - this.Im * other.Im,
            this.Re * other.Im + this.Im * other.Re];
        return this;
    }

    /**
     * Replace the Complex number by its reciprocal.
     *
     * The reciprocal of a nonzero complex number `a + b i` is
     * `(a - b i)/(a^2 + b^2)`.
     *
     * @example
     * let z = new Complex(3, 4);
     * z.recip().equals(new Complex(3/25, -4/25)) == true;
     *
     * @return {Complex} `this`
     * @throws Will throw an error if `this` is zero.
     */
    recip() {
        let s = 1/this.sizesq;
        if(!isFinite(s))
            throw new Error("Tried to divide by zero");
        this.Re *= s;
        this.Im *= -s;
        return this;
    }

    /**
     * Divide by a complex number in-place.
     *
     * This is the same as multiplying by the reciprocal.
     *
     * @example
     * let z = new Complex(1, 2), w = new Complex(3, 4);
     * z.div(w).equals(new Complex(11/25, 2/25)) == true;
     *
     * @param {(Complex|number)} other - The number to divide.
     * @return {Complex} `this`
     * @throws Will throw an error if `other` is zero.
     */
    div(other) {
        if(typeof other === "number") {
            let s = 1/other;
            if(!isFinite(s))
                throw new Error("Tried to divide by zero");
            return this.scale(s);
        }
        return this.mult(other.clone().recip());
    }
}


////////////////////////////////////////////////////////////////////////////////
// * Polynomial root solvers

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
        if(z.Im < 0) z.conj();
        cplx = [[z, 1], [z.clone().conj(), 1]];
    }
    ret.sort((a, b) => a[0] - b[0]).forEach(a => a[0] -= b/3);
    return ret.concat(cplx);
}


/**
 * Class representing a matrix.
 *
 * A matrix is a 2 dimensional array of numbers.  The first dimension indexes
 * the row, written horizontally; the second indexes the column.
 *
 * The rows of a `Matrix` are `Vector` instances.  All rows must have the same
 * length.
 *
 * @extends Array
 */
export class Matrix extends Array {
    /**
     * Create an `n`x`n` identity Matrix.
     *
     * @example
     * Matrix.identity(3).equals(new Matrix([1,0,0],[0,1,0],[0,0,1])) == true;
     *
     * @param {integer} n - The resulting Matrix will have this many rows and columns.
     * @param {number} [λ=1] - The diagonal entries will be equal to `λ`.
     * @return {Matrix}
     */
    static identity(n, λ=1) {
        let rows = new Array(n);
        for(let i = 0; i < n; ++i)
            rows[i] = Vector.e(i, n, λ);
        return new Matrix(...rows);
    }

    /**
     * Create an `m`x`n` zero Matrix.
     *
     * @example
     * Matrix.zero(2, 3).equals(new Matrix([0,0,0],[0,0,0])) == true;
     *
     * @param {integer} m - The resulting Matrix will have this many rows.
     * @param {integer} [n=m] - The resulting Matrix will have this many columns.
     * @return {Matrix}
     */
    static zero(m, n=m) {
        let rows = new Array(m);
        for(let i = 0; i < m; ++i)
            rows[i] = Vector.zero(n);
        return new Matrix(...rows);
    }

    /**
     * Create a permutation Matrix.
     *
     * This is the matrix whose `i`th row is the `j`th unit coordinate
     * vector.
     *
     * @example
     * Matrix.permutation(1, 0).equals(new Matrix([0,1],[1,0])) == true;
     *
     * @param {number[]} vals - If `n` numbers are given, these should be a
     *   permutation of the set {0,1,...,n-1}, although no error is thrown if
     *   each argument is between 0 and n-1.
     * @return {Matrix}
     */
    static permutation(vals) {
        let ret = Matrix.zero(vals.length);
        for(let i = 0; i < vals.length; ++i)
            ret[i][vals[i]] = 1;
        return ret;
    }

    /**
     * Create a Matrix.
     *
     * @param {...(Array<number>|Vector)} rows - The rows of the matrix.  Arrays
     *   will be promoted to Vector instances.  All rows must have the same
     *   length.
     * @throws Will throw an error if the rows do not have the same length.
     */
    constructor(...rows) {
        if(rows.length === 1 && typeof rows[0] === "number") {
            // Need to support Array(n) semantics so this.map() works
            super(rows[0]);
            this.invalidate();
            return;
        }
        if(rows.length === 0) {
            super();
            this.invalidate();
            return;
        }
        super(...rows.map(
            row => row instanceof Vector ? row
                : row.length == 1 ? Vector.constant(1, row[0])
                : new Vector(...row)));
        let n = this[0].length;
        if(this.some(row => row.length !== n))
            throw new Error("Matrix rows must have the same length.");

        this.invalidate();
    }

    /**
     * Invalidate cached computations.
     *
     * Call this after modifying any matrix entries.  It is not called
     * automatically for performance reasons.
     *
     * @return {undefined}
     */
    invalidate() {
        if(this.__cache) delete this.__cache;
    }

    /**
     * @private
     *
     * @type {object}
     * @property {integer} rank - The rank of the matrix.
     * @property {PLUData} PLU - PLU factorization; computed in PLU().
     * @property {Matrix} rref - Reduced row echelon form.
     * @property {Matrix} E - Matrix such that `E*this = rref`.
     * @property {Vector[]} nullBasis - Basis for the null space.
     * @property {Subspace} nullSpace - The null space.
     * @property {Vector[]} colBasis - Basis for the column space.
     * @property {Subspace} nullSpace - The column space.
     * @property {Vector[]} rowBasis - Basis for the row space.
     * @property {Subspace} rowSpace - The row space.
     * @property {Vector[]} leftNullBasis - Basis for the left null space.
     * @property {Subspace} leftNullSpace - The left null space.
     * @property {QLData} QR - QR factorization; computed in QR().
     */
    get _cache() {
        if(!this.__cache) this.__cache = {};
        return this.__cache;
    }

    /**
     * The number of rows of the matrix.
     *
     * @type {integer}
     */
    get m() { return this.length; }

    /**
     * The number of columns of the matrix.
     *
     * @type {integer}
     */
    get n() { return this.length === 0 ? 0 : this[0].length; }

    /**
     * Whether the matrix is square.
     *
     * @type {boolean}
     */
    get isSquare() {
        return this.m == this.n;
    }

    /**
     * Check if this matrix is equal to `other`.
     *
     * Two matrices are equal if they have the same dimensions, and all
     * entries are equal.
     *
     * @param {Matrix} other - The matrix to compare.
     * @param {number} [ε=0] - Entries will test as equal if they are within `ε`
     *   of each other.  This is provided in order to account for rounding
     *   errors.
     * @return {boolean}
     */
    equals(other, ε=0) {
        if(this.m !== other.m || this.n !== other.n)
            return false;
        return this.every((v, i) => v.equals(other[i], ε));
    }

    /**
     * Create a new Matrix with the same entries.
     *
     * @return {Matrix}
     */
    clone() {
        // this.map() calls the Matrix() constructor
        return this.map(row => row.clone());
    }

    /**
     * Return a string representation of the matrix.
     *
     * @param {integer} [precision=4] - The number of decimal places to include.
     * @return {string}
     */
    toString(precision=4) {
        // this.map returns a Matrix
        let strings = new Array(this.m), i = 0, j;
        for(let row of this) {
            let rowStrs = new Array(row.length);
            j = 0;
            for(let v of row)
                rowStrs[j++] = v.toFixed(precision);
            strings[i++] = rowStrs;
        }
        let colLens = new Array(this.n);
        for(j = 0; j < this.n; ++j)
            colLens[j] = Math.max(...strings.map(row => row[j].length));
        let ret = new Array(strings.length);
        i = 0;
        for(let row of strings)
            ret[i++] = row.map((val, j) => val.padStart(colLens[j], ' ')).join(' ');
        return ret.join('\n');
    }

    /**
     * Return the `i`th row.
     *
     * @param {integer} i - The row to return.
     * @return {Vector} The `i`th row of `this`.
     */
    row(i) {
        return this[i];
    }

    /**
     * Return an iterable over the rows.
     *
     * @return {Object}
     */
    rows() {
        return this[Symbol.iterator]();
    }

    /**
     * Return the `j`th column.
     *
     * @param {integer} j - The column to return.
     * @return {Vector} The `j`th column of `this`.
     */
    col(j) {
        let v = Vector.zero(this.m), i = 0;
        for(let row of this)
            v[i++] = row[j];
        return v;
    }

    /**
     * Return an iterable over the columns.
     *
     * @return {Object}
     */
    cols() {
        let self = this;
        return (function*() {
            for(let j = 0; j < self.n; ++j)
                yield self.col(j);
        })();
    }

    /**
     * Add a Matrix in-place.
     *
     * This modifies the matrix in-place by adding the entries of `other`.
     *
     * @param {Matrix} other - The matrix to add.
     * @param {number} [factor=1] - Add `factor` times `other` instead of just
     *   adding `other`.
     * @return {Matrix} `this`
     * @throws Will throw an error if the matrices have different sizes.
     */
    add(other, factor=1) {
        if(this.m !== other.m || this.n !== other.n)
            throw new Error(
                'Tried to add matrices of different sizes');
        this.forEach((row, i) => row.add(other[i], factor));
        return this;
    }

    /**
     * Subtract a Matrix in-place.
     *
     * This modifies the matrix in-place by subtracting the entries of `other`.
     *
     * @param {Matrix} other - The matrix to add.
     * @return {Matrix} `this`
     * @throws Will throw an error if the matrices have different sizes.
     */
    sub(other) {
        return this.add(other, -1);
    }

    /**
     * Scale a Matrix in-place.
     *
     * This modifies the matrix in-place by multiplying all entries by `c`.
     *
     * @param {number} c - The number to multiply.
     * @return {Matrix} `this`
     */
    scale(c) {
        for(let row of this)
            row.scale(c);
        return this;
    }

    /**
     * Create the transpose Matrix.
     *
     * The rows of the transpose are the columns of the matrix.
     *
     * @return {Matrix}
     */
    transpose() {
        let ret = new Matrix(this.n);
        for(let j = 0; j < this.n; ++j)
            ret[j] = this.col(j);
        return ret;
    }

    /**
     * Multiply by another matrix.
     *
     * This creates a new matrix equal to the product of `this` and `other`.
     * This only makes sense when the number of rows of `other` equals the
     * number of columns of `this`.
     *
     * @param {Matrix} other - The matrix to multiply.
     * @return {Matrix} The matrix product.
     * @throws Will throw an error if the matrices have incompatible
     * dimensions.
     */
    mult(other) {
        if(other.m !== this.n)
            throw new Error('Cannot multiply matrices of incompatible dimensions');
        return this.map(
            row => {
                let newRow = Vector.zero(other.n);
                for(let i = 0; i < other.n; ++i) {
                    for(let j = 0; j < this.n; ++j)
                        newRow[i] += row[j] * other[j][i];
                }
                return newRow;
            });
    }

    /**
     * Multiply a matrix times a vector.
     *
     * This creates a new vector equal to the product of `this` and `v`.  This
     * only makes sense when the number of entries of `v` equals the number of
     * columns of `this`.
     *
     * This is an optimized special case of {@link Matrix#mult}.
     *
     * @param {Vector} v - The vector to multiply.
     * @return {Vector} The matrix-vector product.
     * @throws Will throw an error if the matrix and vector have incompatible
     *   dimensions.
     */
    apply(v) {
        if(v.length !== this.n)
            throw new Error('Cannot multiply matrix and vector of incompatible dimensions');
        let ret = Vector.zero(this.m), i = 0;
        for(let row of this)
            ret[i++] = row.dot(v);
        return ret;
    }

    /**
     * Whether the matrix is zero.
     *
     * @param {number} [ε=0] - Entries smaller than this value are assumed
     *   to be zero.
     * @return {boolean}
     */
    isZero(ε=0) {
        return this.every(row => row.isZero(ε));
    }

    /**
     * Whether the matrix is upper-triangular.
     *
     * This means that all entries below the main diagonal are equal to zero.
     *
     * @param {number} [ε=0] - Entries smaller than this value are assumed
     *   to be zero.
     * @return {boolean}
     */
    isUpperTri(ε=0) {
        for(let i = 1; i < this.m; ++i) {
            for(let j = 0; j < i; ++j) {
                if(Math.abs(this[i][j]) > ε)
                    return false;
            }
        }
        return true;
    }

    /**
     * Whether the matrix is upper-triangular with ones on the diagonal.
     *
     * @param {number} [ε=0] - Entries smaller than this value are assumed
     *   to be zero.
     * @return {boolean}
     */
    isUpperUnip(ε=0) {
        for(let i = 1; i < Math.min(this.m, this.n); ++i) {
            if(Math.abs(this[i][i] - 1) > ε)
                return false;
        }
        return this.isUpperTri(ε);
    }

    /**
     * Whether the matrix is lower-triangular.
     *
     * This means that all entries above the main diagonal are equal to zero.
     *
     * @param {number} [ε=0] - Entries smaller than this value are assumed
     *   to be zero.
     * @return {boolean}
     */
    isLowerTri(ε=0) {
        for(let i = 0; i < this.m; ++i) {
            for(let j = i+1; j < this.n; ++j) {
                if(Math.abs(this[i][j]) > ε)
                    return false;
            }
        }
        return true;
    }

    /**
     * Whether the matrix is lower-triangular with ones on the diagonal.
     *
     * @param {number} [ε=0] - Entries smaller than this value are assumed
     *   to be zero.
     * @return {boolean}
     */
    isLowerUnip(ε=0) {
        for(let i = 1; i < Math.min(this.m, this.n); ++i) {
            if(Math.abs(this[i][i] - 1) > ε)
                return false;
        }
        return this.isLowerTri(ε);
    }

    /**
     * Whether the matrix is diagonal.
     *
     * This means the matrix is square, and both upper- and lower-triangular.
     *
     * @param {number} [ε=0] - Entries smaller than this value are assumed
     *   to be zero.
     * @return {boolean}
     */
    isDiagonal(ε=0) {
        if(!this.isSquare) return false;
        return this.isLowerTri(ε) && this.isUpperTri(ε);
    }


    /**
     * A pivot position is recorded as a pair `[row, column]`.
     *
     * @typedef {number[]} Pivot
     */

    /**
     * A list of leading entries of each row.
     *
     * Zero rows do not have a leading entry, and are not included.
     *
     * @param {number} [ε=0] - Entries smaller than this value are assumed
     *   to be zero.
     * @return {Pivot[]}
     */
    leadingEntries(ε=0) {
        let entries = [];
        for(let i = 0; i < this.m; ++i) {
            for(let j = 0; j < this.n; ++j) {
                if(Math.abs(this[i][j]) > ε) {
                    entries.push([i, j]);
                    break;
                }
            }
        }
        return entries;
    }

    /**
     * Whether the matrix is in row-echelon form.
     *
     * This means that the first nonzero entry of each row is to the right of
     * the first nonzero entry of the previous row, which implies that all
     * entries below a pivot are zero and all zero rows are at the bottom.
     *
     * @param {number} [ε=0] - Entries smaller than this value are assumed
     *   to be zero.
     * @return {boolean}
     */
    isEchelon(ε=0) {
        let i0 = -1;
        let j0 = -1;
        for(let [i, j] of this.leadingEntries(ε)) {
            if(i !== i0 + 1 || j <= j0)
                return false;
            i0 = i;
            j0 = j;
        }
        return true;
    }

    /**
     * Whether the matrix is in reduced row-echelon form.
     *
     * This means that it is in row-echelon form, and in addition, all pivots
     * are equal to one, and all entries above a pivot are equal to zero.
     *
     * @param {number} [ε=0] - Entries smaller than this value are assumed
     *   to be zero.
     * @return {boolean}
     */
    isRREF(ε=0) {
        let i0 = -1;
        let j0 = -1;
        for(let [i, j] of this.leadingEntries(ε)) {
            if(i !== i0 + 1 || j <= j0)
                return false;
            if(this[i][j] != 1)
                return false;
            for(let k = 0; k < i; ++k)
                if(this[k][j] !== 0)
                    return false;
            i0 = i;
            j0 = j;
        }
        return true;
    }

    /**
     * The sum of the diagonal entries of the matrix.
     *
     * @type {number}
     */
    get trace() {
        let acc = 0;
        for(let i = 0; i < Math.min(this.m, this.n); ++i)
            acc += this[i][i];
        return acc;
    }

    // Row operations

    /**
     * Scale a row by a constant.
     *
     * This modifies the matrix in-place by multiplying all entries of row `i`
     * by the scalar `c`.
     *
     * @param {integer} i - The row to scale.
     * @param {number} c - The scaling factor.
     * @param {integer} [start=0] - Only scale the entries `start...this.n`.
     *   Provided for optimizations when the entries before `start` are known to
     *   be zero.
     * @return {Matrix} `this`
     */
    rowScale(i, c, start=0) {
        this[i].scale(c, start);
        return this;
    }

    /**
     * Add a constant times one row to another row.
     *
     * This modifies the matrix in-place by adding all entries of row `i2` times
     * the scalar `c` to the entries of row `i1`.
     *
     * @param {integer} i1 - The row to replace.
     * @param {integer} i2 - The row to add.
     * @param {number} c - The scaling factor.
     * @param {integer} [start=0] - Only add the entries `start...this.n`.
     *   Provided for optimizations when the entries before `start` are known to
     *   be zero.
     * @return {Matrix} `this`
     */
    rowReplace(i1, i2, c, start=0) {
        this[i1].add(this[i2], c, start);
        return this;
    }

    /**
     * Swap two rows.
     *
     * This exchanges rows `i1` and `i2`.
     *
     * @param {integer} i1 - The first row.
     * @param {integer} i2 - The second row.
     * @return {Matrix} `this`
     */
    rowSwap(i1, i2) {
        [this[i1], this[i2]] = [this[i2], this[i1]];
        return this;
    }

    // Gaussian elimination

    /**
     * The core data computed by Gaussian elimination.
     *
     * @typedef PLUData
     * @type {object}
     * @property {number[]} P - A permutation of the numbers `1...m-1`.
     * @property {number} detP - The sign of the permutation P.
     * @property {Matrix} L - An `m`x`m` lower-triangular matrix with ones on
     *   the diagonal.
     * @property {Matrix} U - An `m`x`n` matrix in row echelon form.
     * @property {Matrix} E - An `m`x`m` invertible matrix equal to `L^(-1)P`,
     *   so `EA = U`.
     * @property {Pivot[]} pivots - An array of pivot positions.
     */

    /**
     * PA = LU factorization.
     *
     * This computes an `m`x`m` permutation matrix `P`, an `m`x`m`
     * lower-triangular matrix `L` with ones on the diagonal, and a matrix `U`
     * in row echelon form, such that `PA = LU`.
     *
     * This is the core method that implements Gaussian elimination.  It uses
     * maximal partial pivoting.  Many properties of the matrix are derived from
     * the output of this method, which is cached.
     *
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting.
     * @return {PLUData}
     */
    PLU(ε=1e-10) {
        if(this._cache.PLU) return this._cache.PLU;

        let P = new Array(this.m);
        for(let i = 0; i < this.m; ++i) P[i] = i;
        let L = Matrix.identity(this.m);
        let E = Matrix.identity(this.m);
        let U = this.clone();
        let {m, n} = this;
        let pivots = [];
        let detP = 1;

        for(let curRow = 0, curCol = 0; curRow < m && curCol < n; ++curCol) {
            // Find maximal pivot
            let pivot = U[curRow][curCol], row = curRow;
            for(let i = curRow+1; i < m; ++i) {
                if(Math.abs(U[i][curCol]) > Math.abs(pivot)) {
                    pivot = U[i][curCol];
                    row = i;
                }
            }

            if(Math.abs(pivot) > ε) {
                // curCol is a pivot column
                if(row != curRow) {
                    // Row swap
                    [P[row], P[curRow]] = [P[curRow], P[row]];
                    detP *= -1;
                    U.rowSwap(row, curRow);
                    E.rowSwap(row, curRow);
                    for(let j = 0; j < curRow; ++j)
                        [L[row][j], L[curRow][j]]
                            = [L[curRow][j], L[row][j]];
                }
                // Eliminate
                for(let i = curRow+1; i < m; ++i) {
                    let l = U[i][curCol] / pivot;
                    L[i][curRow] = l;
                    U[i][curCol] = 0;
                    U.rowReplace(i, curRow, -l, curCol+1);
                    E.rowReplace(i, curRow, -l);
                }
                pivots.push([curRow, curCol]);
                curRow++;

            } else {
                // Clear the column so U is really upper-triangular
                for(let i = curRow; i < m; ++i)
                    U[i][curCol] = 0;
            }
        }

        this._cache.PLU = {P, L, U, E, pivots, detP};
        return this._cache.PLU;
    }

    /**
     * Find some solution `x` to the equation `Ax=b`.
     *
     * This is the solution obtained by forward- and back-substituting in the
     * `PA=LU` decomposition, taking the free variables equal to zero.  All
     * other solutions can be obtained by adding elements of the null space.
     *
     * @param {Vector} b - The right-hand side of the equation.
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting.
     * @return {(Vector|null)} Returns `null` if no solution exists.
     * @throws Will throw an error if `b.length != this.m`.
     */
    solve(b, ε=1e-10) {
        if(b.length != this.m)
            throw new Error("Incompatible dimensions of matrix and vector");
        let {P, L, U, pivots} = this.PLU(ε);
        let r = this.rank;
        let Pb = Vector.zero(this.m);
        // Solve LUx = PAx = Pb
        for(let i = 0; i < this.m; ++i) Pb[i] = b[P[i]];
        // Solve Ly = Pb by forward-substitution
        let y = Pb.clone();
        for(let i = 1; i < this.m; ++i) {
            for(let ii = 0; ii < i; ++ii)
                y[i] -= L[i][ii] * y[ii];
        }
        // Check if a solution exists
        for(let i = r; i < this.m; ++i) {
            if(Math.abs(y[i]) > ε)
                return null;
        }
        // Solve Ux = y by reverse-substitution
        let x = Vector.zero(this.n);
        for(let p = r-1; p >= 0; --p) {
            let [row, col] = pivots[p];
            x[col] = y[row];
            for(let pp = p+1; pp < r; ++pp) {
                let [, col1] = pivots[pp];
                x[col] -= U[row][col1] * x[col1];
            }
            x[col] /= U[row][col];
        }
        return x;
    }

    /**
     * Compute the pivot positions of the matrix.
     *
     * This requires Gaussian elimination to compute.
     *
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting.
     * @return {Pivot[]}
     */
    pivots(ε=1e-10) {
        let {pivots} = this.PLU(ε);
        return pivots;
    }

    /**
     * The rank of the matrix.
     *
     * This is computed as a side-effect of other computations.
     *
     * @type {integer}
     */
    get rank() {
        if(this._cache.rank === undefined) // Hasn't been computed yet
            this._cache.rank = this.pivots().length;
        return this._cache.rank;
    }

    /**
     * The nullity of the matrix.
     *
     * @type {integer}
     */
    get nullity() {
        return this.n - this.rank;
    }

    /**
     * The reduced row echelon form of the matrix.
     *
     * This completes the Gauss--Jordan elimination on the row echelon form of
     * the matrix, by scaling so all pivots are equal to 1 and performing row
     * replacements so that all entries above a pivot are equal to zero.
     *
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting.
     * @return {Matrix}
     */
    rref(ε=1e-10) {
        if(this._cache.rref) return this._cache.rref;
        let {U, E, pivots} = this.PLU(ε);
        let rowOps = E.clone();
        let rref = U.clone();
        // Start with the last pivot
        for(let k = pivots.length-1; k >= 0; --k) {
            let [row, col] = pivots[k];
            let pivot = rref[row][col];
            for(let i = 0; i < row; ++i) {
                rref.rowReplace(i, row, -rref[i][col]/pivot, col+1);
                rowOps.rowReplace(i, row, -rref[i][col]/pivot);
                rref[i][col] = 0;
            }
            rref.rowScale(row, 1/pivot, col+1);
            rowOps.rowScale(row, 1/pivot);
            rref[row][col] = 1;
        }
        this._cache.rref = rref;
        this._cache.rowOps = rowOps;
        return this._cache.rref;
    }

    /**
     * Return the row operations used in Gauss--Jordan elimination.
     *
     * This returns an invertible `m`x`m` matrix `E` such that `EA = rref`.
     *
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting.
     * @return {Matrix}
     */
    rowOps(ε=1e-10) {
        if(this._cache.rowOps) return this._cache.rowOps;
        this.rref(ε);
        return this._cache.rowOps;
    }

    // Bases

    /**
     * Basis for the null space of the matrix.
     *
     * This is computed using the reduced row echelon form and parametric vector
     * form.
     *
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting.
     * @return {Vector[]}
     */
    nullBasis(ε=1e-10) {
        if(this._cache.nullBasis) return this._cache.nullBasis;
        let rref = this.rref(ε), pivots = this.pivots(ε).slice(), previous = [];
        let basis = [];
        for(let j = 0; j < this.n; ++j) {
            if(pivots.length && pivots[0][1] === j) {
                // Pivot column
                previous.push(pivots.shift());
                continue;
            }
            // Free column
            let v = Vector.zero(this.n);
            for(let [row, col] of previous)
                v[col] = -rref[row][j];
            v[j] = 1;
            basis.push(v);
        }
        this._cache.nullBasis = basis;
        return this._cache.nullBasis;
    }

    /**
     * The null space of the matrix.
     *
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting.
     * @return {Subspace}
     */
    nullSpace(ε=1e-10) {
        if(this._cache.nullSpace) return this._cache.nullSpace;
        this._cache.nullSpace = new Subspace(
            this.nullBasis(ε), {n: this.n, isBasis: true});
        return this._cache.nullSpace;
    }

    /**
     * Basis for the column space of the matrix.
     *
     * This is composed of the pivot columns.
     *
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting.
     * @return {Vector[]}
     */
    colBasis(ε=1e-10) {
        if(this._cache.colBasis) return this._cache.colBasis;
        let {pivots} = this.PLU(ε);
        this._cache.colBasis = pivots.map(([, j]) => this.col(j));
        return this._cache.colBasis;
    }

    /**
     * The column space of the matrix.
     *
     * @return {Subspace}
     */
    colSpace() {
        if(this._cache.colSpace) return this._cache.colSpace;
        this._cache.colSpace = new Subspace(this, {n: this.m});
        return this._cache.colSpace;
    }

    /**
     * Basis for the row space of the matrix.
     *
     * This is composed of the nonzero rows of a REF.
     *
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting.
     * @return {Vector[]}
     */
    rowBasis(ε=1e-10) {
        if(this._cache.rowBasis) return this._cache.rowBasis;
        let {pivots, U} = this.PLU(ε);
        this._cache.rowBasis = pivots.map(([i,]) => U[i].clone());
        return this._cache.rowBasis;
    }

    /**
     * The row space of the matrix.
     *
     * @return {Subspace}
     */
    rowSpace() {
        if(this._cache.rowSpace) return this._cache.rowSpace;
        this._cache.rowSpace = new Subspace(this.transpose(), {n: this.n});
        return this._cache.rowSpace;
    }

    /**
     * Basis for the left null space of the matrix.
     *
     * This is composed of the last `m-r` rows of `L^-1 P`, where `P` and `L`
     * come from the `PA = LU` decomposition and `r` is the rank.
     *
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting.
     * @return {Vector[]}
     */
    leftNullBasis(ε=1e-10) {
        if(this._cache.leftNullBasis) return this._cache.leftNullBasis;
        let {E} = this.PLU(ε);
        let r = this.rank;
        let ret = new Array(this.m - r);
        for(let i = r; i < this.m; ++i)
            ret[i - r] = E[i].clone();
        this._cache.leftNullBasis = ret;
        return this._cache.leftNullBasis;
    }

    /**
     * The left null space of the matrix.
     *
     * @return {Subspace}
     */
    leftNullSpace(ε=1e-10) {
        if(this._cache.leftNullSpace) return this._cache.leftNullSpace;
        this._cache.leftNullSpace = new Subspace(
            this.leftNullBasis(ε), {n: this.n, isBasis: true});
        return this._cache.leftNullSpace;
    }

    /**
     * Whether the matrix has full row rank.
     *
     * This means that there is a pivot in every row.
     *
     * @type {boolean}
     */
    get isFullRowRank() {
        if(this.m > this.n) return false; // Don't need to row reduce
        return this.rank == this.m;
    }

    /**
     * Whether the matrix has full column rank.
     *
     * This means that there is a pivot in every column.
     *
     * @type {boolean}
     */
    get isFullColRank() {
        if(this.m < this.n) return false; // Don't need to row reduce
        return this.rank == this.n;
    }

    /**
     * Whether the matrix is invertible.
     *
     * This means that the matrix is square and has the maximum number of
     * pivots.  Equivalently, there is another matrix (namely `this.inverse()`)
     * such that the product with `this` on both sides is equal to the identity.
     *
     * @type {boolean}
     */
    get isInvertible() {
        return this.isFullRowRank && this.isFullColRank;
    }

    /**
     * Alias for `!this.isInvertible()`.
     *
     * @type {boolean}
     */
    get isSingular() {
        return !this.isInvertible;
    }

    /**
     * Compute the inverse matrix.
     *
     * This only makes sense for square matrices.  Throws an error if the matrix
     * is not square or is not invertible.
     *
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting.
     * @return {Matrix}
     * @throws Will throw an error if the matrix is not square or not
     *   invertible.
     */
    inverse(ε=1e-10) {
        if(!this.isSquare)
            throw new Error("Tried to invert a non-square matrix");
        let E = this.rowOps(ε);
        if(!this.isInvertible)
            throw new Error("Tried to invert a singular matrix");
        return E;
    }

    /**
     * Compute the matrix determinant.
     *
     * This only makes sense for square matrices.  Throws an error if the matrix
     * is not square.
     *
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting.
     * @return {number}
     * @throws Will throw an error if the matrix is not square.
     */
    det(ε=1e-10) {
        if(this._cache.det !== undefined) return this._cache.det;
        if(!this.isSquare)
            throw new Error("Tried to take the determinant of a non-square matrix");
        let {U, detP} = this.PLU(ε);
        let det = 1;
        for(let i = 0; i < this.n; ++i)
            det *= U[i][i];
        this._cache.det = det * detP;
        return this._cache.det;
    }

    // Orthogonality

    /**
     * Check if the matrix is orthogonal.
     *
     * That means that the matrix is square and has orthonormal columns.
     *
     * @param {number} [ε=1e-10] - Numbers smaller than this value are taken
     *   to be zero.
     * @return {QRData}
     */
    isOrthogonal(ε=1e-10) {
        return this.isSquare &&
            this.transpose().mult(this).equals(Matrix.identity(this.n), ε);
    }

    /**
     * The data computed Gram--Schmidt.
     *
     * @typedef QRData
     * @type {object}
     * @property {Matrix} Q - An `m`x`n` matrix with orthogonal columns.
     * @property {Matrix} R - An `n`x`n` upper-triangular matrix.
     * @property {number[]} LD - A list of the zero columns of `Q`.
     */

    /**
     * QR decomposition.
     *
     * This computes a matrix `Q` and an upper-triangular matrix `R` such that
     * `this = QR`.  If `this` has full column rank then `Q` has orthonormal
     * columns and `R` is invertible.  Otherwise `Q` may have zero columns and
     * `R` may have zero rows; the nonzero columns of `Q` are orthonormal.
     *
     * This implements the Gram--Schmidt algorithm, which is not a very
     * efficient or numerically accurate.
     *
     * @param {number} [ε=1e-10] - Vectors smaller than this value are taken
     *   to be zero.
     * @return {QRData}
     */
    QR(ε=1e-10) {
        if(this._cache.QR) return this._cache.QR;
        let {m, n} = this;
        let ui = new Array(n), LD=[];
        let vi = [...this.cols()];
        let R = Matrix.zero(n);
        for(let j = 0; j < n; ++j) {
            let u = vi[j].clone();
            // Compute ui[j] and the jth column of R at the same time
            for(let jj = 0; jj < j; ++jj) {
                let factor = ui[jj].dot(vi[j]);
                u.sub(ui[jj].clone().scale(factor));
                R[jj][j] = factor;
            }
            let l = u.size;
            if(l > ε) {
                R[j][j] = l;
                ui[j] = u.scale(1/l);
            } else {
                R[j][j] = 0;
                ui[j] = Vector.zero(m);
                LD.push(j);
            }
        }
        let Q = new Matrix(...ui).transpose();
        this._cache.QR = {Q, R, LD};
        if(this._cache.rank === undefined)
            this._cache.rank = n - LD.length;
        return this._cache.QR;
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

    // Eigenvalues and eigenspaces.  Only implemented for 1x1, 2x2, 3x3.
    _calcEigendata() {
        if(!this.isSquare)
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
                .sub(Matrix.identity(this.n, val.value))
                .nullSpace;
            // Sanity check
            if(val.eigenspace.dim === 0 || val.eigenspace.dim > val.mult)
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
            let [a,b,c] = this[0];
            let [d,e,f] = this[1];
            let [g,h,i] = this[2];
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


/**
 * Class representing a subspace.
 *
 * A `Subspace` is represented as the column space of a matrix.  Various
 * operations on the subspace will involve computing a basis, or even an
 * orthonormal basis.
 */
export class Subspace {
    /**
     * The subspace `R^n`.
     *
     * @param {integer} n - The dimension of the space.
     * @return {Subspace} The subspace `R^n`.
     */
    static Rn(n) {
        return new Subspace(Matrix.identity(n), {n: n, isBasis: true, isON: true});
    }

    /**
     * The zero subspace of `R^n`.
     *
     * @param {integer} n - The ambient dimension.
     * @return {Subspace} The subspace `{0}`.
     */
    static zero(n) {
        return new Subspace([], {n: n, isBasis: true, isON: true});
    }

    /**
     * Create a Subspace.
     *
     * @param {(Vector[]|Matrix)} generators - A spanning set for the subspace.
     *   If it is a matrix, use the columns of the matrix as the spanning set.
     * @param {Object} hints - Precomputed properties of the subspace.
     * @param {integer} hints.n - The dimension of the ambient R^n.  Only
     *   necessary to provide if `generators` is empty.
     * @param {boolean} hints.isBasis - The generators are linearly independent.
     * @param {boolean} hints.isON - The generators are orthonormal.
     * @throws Will throw an error if the generators do not all have the same
     *   length, or if it can't determine the ambient space R^n.
     */
    constructor(generators, {n, isBasis, isON}={}) {
        /**
         * The dimension of the ambient R^n.
         *
         * @type {integer}
         */
        if(n !== undefined)
            this.n = n;
        else if(generators instanceof Matrix)
            this.n = generators.m;
        else if(generators.length > 0)
            this.n = generators[0].length;
        else
            throw new Error("Cannot determine the ambient R^n from zero generators");

        /**
         * The columns of this matrix span this Subspace.
         *
         * @type {Matrix}
         */
        this.generators = generators instanceof Matrix
            ? generators : new Matrix(...generators).transpose();

        this._basis = isBasis ? this.generators : undefined;

        /** @private */
        this._ONbasis = isON ? this.generators : undefined;
    }

    /**
     * Return a string representation of the Subspace.
     *
     * @param {integer} [precision=4] - The number of decimal places to include.
     * @return {string}
     */
    toString(precision=4) {
        return `Subspace of R^${this.n} spanned by\n`
            + this.generators.toString(precision);
     }


    /**
     * Compute a basis for the subspace.
     *
     * The computed basis will be a subset of `this.generators`.
     *
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting.
     * @return {Matrix} A matrix whose columns form a basis for the Subspace.
     */
    basis(ε=1e-10) {
        if(this._basis) return this._basis;
        this._basis = new Matrix(...this.generators.colBasis(ε)).transpose();
        return this._basis;
    }

    /**
     * The dimension of the subspace.
     *
     * @type {integer}
     */
    get dim() {
        return this.basis().n;
    }

    /**
     * Compute whether the subspace is all of R^n.
     *
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting.
     * @return {boolean}
     */
    isMaximal(ε=1e-10) {
        return this.basis(ε).n == this.n;
    }

    /**
     * Compute whether the subspace is zero.
     *
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting.
     * @return {boolean}
     */
    isZero(ε=1e-10) {
        return this.basis(ε).n == 0;
    }

    /**
     * Check if this Subspace contains a Vector.
     *
     * This means that the vector is a linear combination of the generators.
     *
     * @param {Vector} vec - The vector to test.
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting.
     * @return {boolean}
     * @throws Will throw an error if `vec.length != this.n`.
     */
    contains(vec, ε=1e-10) {
        if(this.n != vec.length)
            throw new Error("Vector has the wrong number of entries");
        if(this.isMaximal(ε))
            return true; // Nothing to check!
        if(this.isZero(ε))
            return vec.isZero(ε); // Nothing to check!
        return !!this.basis(ε).solve(vec, ε);
    }

    /**
     * Check if this Subspace is contained in another Subspace.
     *
     * @param {Subspace} other - The subspace to test.
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting.
     * @return {boolean}
     */
    isSubspaceOf(other, ε=1e-10) {
        if(this.n != other.n)
            return false;
        let B = this.basis(ε);
        return [...B.cols()].every(col => other.contains(col, ε));
    }

    /**
     * Check if this Subspace is equal to `other`.
     *
     * Two subspaces are equal if they have the same ambient dimension, the same
     * dimension, and contain each other.
     *
     * @param {Subspace} other - The subspace to compare.
     * @param {number} [ε=1e-10] - Entries smaller than this value are taken
     *   to be zero for the purposes of pivoting.
     * @return {boolean}
     */
    equals(other, ε=1e-10) {
        this.basis(ε); other.basis(ε);  // Compute dimension
        if(this.dim !== other.dim)
            return false;
        return this.isSubspaceOf(other, ε);
    }

}

