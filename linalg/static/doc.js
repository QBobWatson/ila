
let s = document.createElement('script');
s.setAttribute('type', 'module');
s.innerHTML = `
import { Vector, Complex, Matrix, Subspace, Polynomial } from "../src/linalg.js";

window.Complex = Complex;
window.Vector = Vector;
window.Matrix = Matrix;
window.Subspace = Subspace;
window.Polynomial = Polynomial;

window.mat = Matrix.create;
window.vec = Vector.create;
window.poly = Polynomial.create;
window.C = (a, b=0) => new Complex(a, b);
`;
document.body.appendChild(s);
