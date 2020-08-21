
# What this library is and is not

This is a self-contained library implementing many basic algorithms and concepts from an introductory linear algebra course.  It was written with two goals in mind:

  1. It should be functional and fast for 3x3 matrices, for use in the interactive demos of the textbook [Interactive Linear Algebra]{@link https://textbooks.math.gatech.edu/ila/}.
  2. It should be suitable as a sandbox for students to play with the basic notions of linear algebra.  As such, it is extensively documented, and can be used from a Javascript console in a browser.

This is *not* a library for serious computation.  Most of the implementations are the naïve schoolbook algorithms: they are meant to illustrate what one learns in the classroom, not what is used in practice.  The implementations are generally not the most efficient or most numerically stable.  If you want to know the state of the art about how matrix computations are performed, see [LAPACK]{@link http://www.netlib.org/lapack/}.

# Use in a browser

While browsing these documentation pages, simply open a [Javascript console]{@link https://webmasters.stackexchange.com/questions/8525} and start typing commands.  The following classes are available in global scope: {@link Matrix}, {@link Vector}, {@link Subspace}, {@link Polynomial}, {@link Complex}. In addition, `mat` is an alias for {@link Matrix.create}, `vec` is an alias for {@link Vector.create}, `poly` is an alias for {@link Polynomial.create}, and `C(a, b)` is an alias for `new Complex(a, b)`.  For example:

![Javascript console example](../static/console.png)

# Getting started

The fundamental class is the {@link Matrix}.  To create a matrix, pass arrays of numbers to {@link Matrix.create}, or its alias `mat`: these are the rows of the matrix.

```javascript
A = mat([1, 2, 3], [4, 5, 6], [7, 8, 9]);
A.toString(0);
  // "[1 2 3]
  //  [4 5 6]
  //  [7 8 9]"
```

There are dozens of operations to perform on matrices, reflecting the basic concepts of linear algebra.  For instance, to compute the reduced row echelon form:

```javascript
A.rref().toString(0);
  // "[1 0 -1]
  //  [0 1  2]
  //  [0 0  0]"
```

To multiply by another matrix:

```javascript
B = mat([1, 1], [2, 2], [3, 4]);
A.mult(B).toString(0);
  // "[14 17]
  //  [32 38]
  //  [50 59]"
```

To work with [Vectors]{@link Vector}, use {@link Vector.create} or its alias `vec`:

```javascript
v = vec(1, -2, 1);
A.apply(v).toString(0);  // [0 0 0]
v.dot(vec(2, 1, 0));     // 0
```

To create a {@link Subspace}, either pass a list of generators, or use Matrix methods:

```javascript
V = A.nullSpace();
V.contains(vec(2, -4, 2));  // true
V.perp().toString(1);
  // "Subspace of R^3 of dimension 2 with basis
  //  [1.0] [0.0]
  //  [0.5] [0.5]
  //  [0.0] [1.0]"
```

# Recommendations

Read the documentation for the individual classes.  The documentation contains links to the source files: click on them!  There is essentially nothing in the source beyond what one would learn in a basic linear algebra class.

# Source code

The latest source can be found on [GitHub]{@link https://github.com/QBobWatson/ila/tree/duke/linalg}.

# License

Copyright (c) 2020 Joseph Rabinoff

This software is released under the [GNU General Public License]{@link https://www.gnu.org/licenses/}.
