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

/** @module complex
 *
 * @file
 * A class for Complex numbers.
 */


/**
 * @summary
 * Class representing a complex number.
 *
 * @desc
 * Complex numbers behave like vectors of length two, in that they can be added
 * and scaled componentwise.  However, complex numbers can also be multiplied
 * together and divided to form new complex numbers.
 *
 * @example {@lang javascript}
 * new Complex(1, 2).toString(1); // "1.0 + 2.0 i"
 *
 * @extends Vector
 */
class Complex {
    /**
     * @summary
     * The real part (first entry) of the complex number.
     *
     * @example {@lang javascript}
     * let z = new Complex(3, 4);
     * z.Re;           // 3
     * z.Re = 5;
     * z.toString(1);  // "5.0 + 4.0 i"
     */
    Re: number;

    /**
     * @summary
     * The imaginary part part (second entry) of the complex number.
     *
     * @example {@lang javascript}
     * let z = new Complex(3, 4);
     * z.Im;           // 4
     * z.Im = 5;
     * z.toString(1);  // "3.0 + 5.0 i"
     */
    Im: number;

    /**
     * @summary
     * Create a new Complex number with prescribed polar coordinates.
     *
     * @desc
     * Cartesian and polar coordinates are related by the formula
     * `(x, y) = (r cos(θ), r sin(θ))`.  According to Euler's formula, the
     * output of `Complex.fromPolar(r, θ)` is equal to `r e^{i θ}`.
     *
     * @example {@lang javascript}
     * Complex.fromPolar(1, Math.PI/2);  // Complex.i
     *
     * @param r - The radial coordinate.
     * @param θ - The angular coordinate.
     * @return The complex number `r e^{i θ}`.
     */
    static fromPolar(r: number, θ: number): Complex {
        return new Complex(r * Math.cos(θ), r * Math.sin(θ));
    }

    /**
     * @summary
     * A copy of `i = new Complex(0, 1)`, a square root of -1.
     *
     * @desc
     * This returns a new object each time it is accessed.
     */
    static get i(): Complex {
        return new Complex(0, 1);
    }

    /**
     * @summary
     * The square of the modulus of the complex number.
     *
     * @desc
     * This is the square of the real part plus the square of the imaginary
     * part.
     *
     * @example {@lang javascript}
     * new Complex(3, 4).modsq;   // 25
     */
    get modsq(): number {
        return this.Re * this.Re + this.Im * this.Im;
    }

    /**
     * @summary
     * The modulus of the complex number.
     *
     * @desc
     * This is an alias for the inherited property `this.size`.
     *
     * @example {@lang javascript}
     * new Complex(3, 4).mod;  // 5
     */
    get mod(): number {
        return Math.sqrt(this.modsq);
    }

    /**
     * @summary
     * The argument of the complex number.
     *
     * @desc
     * This is the angle component of the polar coordinates for the point
     * `(this.Re, this.Im)`.
     *
     * The value is between -π and π.
     *
     * @example {@lang javascript}
     * new Complex(1, 1).arg;  // Math.PI/4
     */
    get arg(): number { return Math.atan2(this.Im, this.Re); }

    /**
     * @param a - The real part.  If `a` is a Complex number, then this
     *   method clones `a`.
     * @param [b=0] - The imaginary part.
     */
    constructor(a: number | Complex, b: number=0) {
        if(a instanceof Complex) {
            this.Re = a.Re;
            this.Im = a.Im;
        }
        else {
            this.Re = a;
            this.Im = b;
        }
    }

    /**
     * @summary
     * Create a new Complex with the same entries.
     *
     * @return The new complex number.
     */
    clone(): Complex {
        return new Complex(this);
    }

    /**
     * @summary
     * Test if this complex number is equal to `other`.
     *
     * @desc
     * This checks that `this.Re == other.Re` and `this.Im == other.Im`, except
     * that if `other` is a number, it is promoted to a Complex number first.
     *
     * @example {@lang javascript}
     * new Complex(1, 0).equals(1);  // true
     *
     * @param other - The number to compare.
     * @param [ε=0] - Entries will test as equal if they are within `ε`
     *   of each other.  This is provided in order to account for rounding
     *   errors.
     * @return True if `this` equals `other`.
     */
    equals(other: Complex | number, ε: number=0): boolean {
        if(typeof other === "number")
            other = new Complex(other, 0);
        if(ε == 0)
            return this.Re == other.Re && this.Im == other.Im;
        return Math.abs(this.Re - other.Re) < ε &&
            Math.abs(this.Im - other.Im) < ε;
    }

    /**
     * @summary
     * Return a string representation of the complex number.
     *
     * @example {@lang javascript}
     * new Complex(1, 2).toString(2);  // "1.00 + 2.00 i"
     *
     * @param [precision=4] - The number of decimal places to include.
     * @return A string representation of the complex number.
     */
    toString(precision: number=4): string {
        if(this.Im >= 0)
            return `${this.Re.toFixed(precision)} + ${this.Im.toFixed(precision)} i`;
        return `${this.Re.toFixed(precision)} - ${(-this.Im).toFixed(precision)} i`;
    }

    /**
     * @summary
     * Replace the Complex number with its complex conjugate.
     *
     * @desc
     * This modifies `this` in place by negating the imaginary part.
     *
     * @example {@lang javascript}
     * let z = new Complex(1, 2);
     * z.conj();
     * z.toString(1);  // "1.0 - 2.0 i"
     *
     * @return `this`
     */
    conj(): this {
        this.Im *= -1;
        return this;
    }

    /**
     * @summary
     * Add another complex number in-place.
     *
     * @desc
     * Add `other.Re` to `this.Re` and `other.Im` to `this.Im`.  If `other` is a
     * number, it is promoted to a Complex number first.
     *
     * @example {@lang javascript}
     * let z = new Complex(1, 2);
     * z.add(1);
     * z.toString(1);  // "2.0 + 2.0 i"
     *
     * @param other - The number to add.
     * @param [factor=1] - Add `factor` times `other` instead of just
     *   adding `other`.
     * @return `this`
     */
    add(other: Complex | number, factor: number=1): this {
        if(typeof other === "number") {
            this.Re += other * factor;
        } else {
            this.Re += other.Re * factor;
            this.Im += other.Im * factor;
        }
        return this;
    }

    /**
     * @summary
     * Subtract another complex number in-place.
     *
     * @desc
     * Subtract `other.Re` from `this.Re` and `other.Im` from `this.Im`.  If
     * `other` is a number, it is promoted to a Complex number first.
     *
     * @example {@lang javascript}
     * let z = new Complex(1, 2);
     * z.sub(1);
     * z.toString(1);  // "0.0 + 2.0 i"
     *
     * @param other - The number to add.
     * @return `this`
     */
    sub(other: Complex | number): this {
        return this.add(other, -1);
    }

    /**
     * @summary
     * Multiply by a complex number in-place.
     *
     * @desc
     * Multiplication of complex numbers is defined by the formula
     * `(a + b i) (c + d i) = (ac - bd) + (ad + bc) i`.
     *
     * @example {@lang javascript}
     * let z = new Complex(1, 2), w = new Complex(3, 4);
     * z.mult(w);
     * z.toString(1);  // "-5.0 + 10.0 i"
     * z.mult(-2);
     * z.toString(1);  // "10.0 - 20.0 i"
     *
     * @param other - The number to multiply.
     * @return `this`
     */
    mult(other: Complex | number): this {
        if(typeof other === "number") {
            this.Re *= other;
            this.Im *= other;
        } else
            [this.Re, this.Im] = [
                this.Re * other.Re - this.Im * other.Im,
                this.Re * other.Im + this.Im * other.Re];
        return this;
    }

    /**
     * @summary
     * Raise to the power `x`.
     *
     * @desc
     * This returns a new Complex number whose modulus is `this.mod` raised to
     * the power `x` and whose argument is `x*this.arg`.
     *
     * @example {@lang javascript}
     * new Complex( 1, 2).pow(2  ).toString(1);  // "-3.0 + 4.0 i"
     * new Complex(-3, 4).pow(1/2).toString(1);  // "1.0 + 2.0 i"
     *
     * @param x - The exponent.
     * @return A new Complex number equal to the `this` raised to the
     *   power `x`.
     */
    pow(x: number): Complex {
        return Complex.fromPolar(Math.pow(this.mod, x), x * this.arg);
    }

    /**
     * @summary
     * Replace the Complex number by its reciprocal.
     *
     * @desc
     * The reciprocal of a nonzero complex number `a + b i` is
     * `(a - b i)/(a^2 + b^2)`.
     *
     * @example {@lang javascript}
     * let z = new Complex(3, 4);
     * z.recip();
     * z.toString();  // "0.1200 - 0.1600 i"
     *
     * @return `this`
     * @throws Error if `this` is zero.
     */
    recip(): this {
        const s = 1/this.modsq;
        if(!isFinite(s))
            throw new Error("Tried to divide by zero");
        this.Re *= s;
        this.Im *= -s;
        return this;
    }

    /**
     * @summary
     * Divide by a complex number in-place.
     *
     * @desc
     * This is the same as `this.mult(other.recip())`, except `other` is not
     * modified.
     *
     * @example {@lang javascript}
     * let z = new Complex(1, 2), w = new Complex(3, 4);
     * z.div(w);
     * z.toString();  // "0.440 + 0.0800 i"
     *
     * @param other - The number to divide.
     * @return `this`
     * @throws Error if `other` is zero.
     */
    div(other: Complex | number): this {
        if(typeof other === "number") {
            const s = 1/other;
            if(!isFinite(s))
                throw new Error("Tried to divide by zero");
            return this.mult(s);
        }
        const s = 1/other.modsq;
        if(!isFinite(s))
            throw new Error("Tried to divide by zero");
        [this.Re, this.Im] = [
            ( this.Re * other.Re + this.Im * other.Im) * s,
            (-this.Re * other.Im + this.Im * other.Re) * s
        ];
        return this;
    }
}


export default Complex;
