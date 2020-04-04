
const {abs, max, min, log, exp, round, sqrt, sin, cos, PI} = Math;

const SMALNO = Number.MIN_VALUE;
const INFIN = Number.MAX_VALUE;
const ETA = Number.EPSILON;
const ARE = ETA;
const MRE = ETA;
const LO = SMALNO / ETA;
const LOG2 = log(2);
const COSR = cos(94*PI/180); // ~ -0.069756474
const SINR = sin(94*PI/180); // ~ 0.99756405
const SQRT2INV = sqrt(0.5); // ~ 0.70710678

/**
 * A polynomial is represented as a list of coefficients.
 *
 * The zeroth coefficient is the highest-order term.
 *
 * @typedef {number[]} Polynomial
 */

/**
 * A complex number is represented as an array `[Re, Im]`.
 *
 * @typedef {number[]} Complex
 */

export default function RPOLY(P) {
    // Make a copy of the coefficients
    P = P.slice();
    // Translate so the leading coefficient is nonzero
    while(P[0] === 0) P.shift();

    let NN = P.length;
    let N = NN - 1;

    let zeros = [];

    // Remove the zeros at the origin if any
    while(P[N] === 0) {
        zeros.push([0, 0]);
        P.pop();
        NN -= 1;
        N  -= 1;
    }

    while(N >= 1) {// Main loop
        if(N == 1) {
            zeros.push([-P[1]/P[0], 0]);
            break;
        } else if(N == 2) {
            zeros.push(...QUAD(...P));
            break;
        }

        // Find the largest and smallest moduli of coefficients
        const MAX = max(...P);
        const MIN = min(...P.filter(x => x !== 0));

        // Scale if there are large or very small coefficients.  Computes a
        // scale factor to multiply the coefficients of the polynomial.  The
        // scaling is done to avoid overflow and to avoid undetected
        // underflow interfering with the convergence criterion.
        let SC = LO/MIN;
        if((SC > 1 && INFIN/SC >= MAX) || (SC <= 1 && MAX >= 10)) {
            if(SC === 0) SC = SMALNO;
            let L = round(log(SC)/LOG2);
            if(L != 0) {
                let FACTOR = L >= 0 ? 1 << L : 1/(1 << -L);
                P = P.map(x => x * FACTOR);
            }
        }

        // Compute lower bound on moduli of zeros.
        let PT = P.map(abs);
        PT[N] *= -1;

        // Compute upper estimate of bound
        let X = exp((log(-PT[N])-log(PT[0]))/N);
        if(PT[N-1] !== 0) {
            // If Newton step at the origin is better, use it.
            let XM = -PT[N]/PT[N-1];
            if(XM < X) X = XM;
        }

        // Chop the interval (0,X) until FF <= 0
        while(true) {
            let XM = X * 0.1;
            let FF = PT.reduce((a, x) => a*XM + x, 0);
            if(FF <= 0)
                break;
            X = XM;
        }

        // Do Newton iteration until X converges to two decimal places
        let DX = X;
        while(abs(DX/X) > .005) {
            let FF = PT[0];
            let DF = FF;
            for(let I = 1; I < N; ++I) {
                FF = FF*X + PT[I];
                DF = DF*X + FF;
            }
            FF = FF*X + PT[N];
            DX = FF/DF;
            X -= DX;
        }
        let BND = X;

        // Compute the derivative as the intial K polynomial and do 5 steps
        // with no shift
        let K = P.slice(0, N).map((x, i) => (N-i) * x / N);

        let AA = P[N];
        let BB = P[N-1];
        let ZEROK = K[N-1] === 0;
        for(let JJ = 0; JJ < 5; ++JJ) {
            let CC = K[N-1];
            if(ZEROK) {
                // Use unscaled form of recurrence
                K.unshift(0);
                K.pop();
                ZEROK = K[N-1] === 0;
            } else {
                // Use scaled form of recurrence if value of K at 0 is nonzero
                let T = -AA/CC;
                for(let J = N-1; J >= 1; --J)
                    K[J] = T*K[J-1] + P[J];
                K[0] = P[0];
                ZEROK = abs(K[N-1]) <= abs(BB)*ETA*10;
            }
        }

        // Save K for restarts with new shifts
        const TEMP = K.slice();

        let XX = SQRT2INV;
        let YY = -XX;

        /**
         * The calculation context.
         *
         * @typedef Context
         * @type Object
         * @property {Polynomial} P  - The current polynomial.
         * @property {Polynomial} QP - The quotient by the current factor.
         * @property {Polynomial} K  - The current K-polynomial.
         * @property {Polynomial} QK - The quotient by the current factor.
         * @property {integer}    NN - The length of `P` and `QP`.
         * @property {integer}    N  - The length of `K` and `QK`.
         * @property {number}     U  - The linear coefficient of the denominator.
         * @property {number}     V  - The constant coefficient of the denominator.
         * @property {number}     A  - The linear coefficient of the remainder.
         * @property {number}     B  - The constant coefficient of the remainder.
         * @property {number}     C  - Scalar set in {@link CALCSC}.
         * @property {number}     D  - Scalar set in {@link CALCSC}.
         * @property {number}     F  - Scalar set in {@link CALCSC}.
         * @property {number}     G  - Scalar set in {@link CALCSC}.
         * @property {number}     H  - Scalar set in {@link CALCSC}.
         * @property {number}     A1 - Scalar set in {@link CALCSC}.
         * @property {number}     A3 - Scalar set in {@link CALCSC}.
         * @property {number}     A7 - Scalar set in {@link CALCSC}.
         */
        let context = {P, K, N, NN};
        context.QP = new Array(NN);
        context.QK = new Array(N);

        // Loop to select the quadratic corresponding to each new shift
        let CNT;
        for(CNT = 1; CNT <= 20; ++CNT) {
            // Quadratic corresponds to a double shift to a non-real point and its
            // complex conjugate.  The point has modulus BND and amplitude rotated
            // by 94 degrees from the previous shift
            [XX, YY] = [COSR*XX - SINR*YY, SINR*XX + COSR*YY];
            let SR = BND*XX;
            let U = -2*SR;
            let V = BND;
            Object.assign(context, {U, V});

            let newZeros = FXSHFR(20 * CNT, SR, context);
            if(newZeros.length !== 0) {
                // The second stage jumps directly to one of the third stage
                // iterations and returns here if successful.  Deflate the
                // polynomial, store the zero or zeros and return to the main
                // algorithm.
                zeros.push(...newZeros);
                NN -= newZeros.length;
                N   = NN - 1;
                P   = context.QP.slice(0, NN);
                break;
            } else {
                // If the iteration is unsuccessful another quadratic is chosen
                // after restoring K
                context.K = TEMP.slice();
            }
        }
        if(CNT > 20)
            break; // Return with failure if no convergence with 20 shifts
    }

    return zeros;
}

/**
 * Compute up to L2 fixed shift K-polynomials, testing for convergence in
 * the linear or quadratic case.  Initiates one of the variable shift
 * iterations and returns with the zeros found.
 *
 * @param {integer} L2 - The number of K-polynomials to compute.
 * @param {number} OSS - The starting real value.
 * @param {Context} context - The calculation context.
 * @return {Complex[]} The zeros found.
 */
function FXSHFR(L2, OSS, context) {
    let {N, P} = context;
    let BETAV = 0.25;
    let BETAS = 0.25;
    let OVV = context. V;
    let UI, VI, VV, SS, TV, TS, TYPE, TVV, TSS, OTV, OTS, VPASS, SPASS,
        SVU, SVV, S, VTRY, STRY, IFLAG, SVK;
    // Evaluate polynomial by synthetic division
    [context.A, context.B] = QUADSD(context.U, context.V, P, context.QP);
    TYPE = CALCSC(context);

    for(let J = 0; J < L2; ++J) {
        // Calculate next K polynomial and estimate V
        NEXTK(TYPE, context);
        TYPE = CALCSC(context);
        [UI, VI] = NEWEST(TYPE, context);
        VV = VI;
        // Estimate S
        SS = context.K[N-1] !== 0 ? -P[N]/context.K[N-1] : 0;
        TV = 1;
        TS = 1;

        if(J > 0 && TYPE !== 3) {
            // Compute relative measures of convergence of S and V sequences
            if(VV !== 0) TV = abs((VV-OVV)/VV);
            if(SS !== 0) TS = abs((SS-OSS)/SS);
            // If decreasing, multiply two most recent convergence measures
            TVV = TV < OTV ? TV*OTV : 1;
            TSS = TS < OTS ? TS*OTS : 1;
            // Compare with convergence criteria
            VPASS = TVV < BETAV;
            SPASS = TSS < BETAS;

            if(SPASS || VPASS) {
                // At least one sequence has passed the convergence test.  Store
                // variables before iterating.
                [SVU, SVV, SVK] = [context.U, context.V, context.K.slice()];
                S = SS;
                // Choose iteration according to the fastest converging
                // sequence.
                VTRY = !VPASS;
                STRY = !SPASS;

                while(!VTRY || !STRY) {
                    if(!STRY && (VTRY || TSS < TVV)) {
                        [S, IFLAG] = REALIT(S, context);
                        if(S !== null) return [[S, 0]];
                        // Linear iteration has failed.  Flag that it has been tried
                        // and decrease the convergence criterion
                        STRY = true;
                        BETAS *= 0.25;
                        if(IFLAG) {
                            // If linear iteration signals an almost double real
                            // zero attempt quadratic interation
                            VTRY = false;
                            UI = -(S+S);
                            VI = S*S;
                            continue;
                        }
                    } else if(!VTRY) {
                        let zeros = QUADIT(UI, VI, context);
                        if(zeros.length > 0) return zeros;
                        // Quadratic iteration has failed.  Flag that it has been
                        // tried and decrease the convergence criterion.
                        VTRY = true;
                        BETAV *= 0.25;
                        if(!STRY) {
                            context.K = SVK.slice();
                            continue;
                        }
                    }
                    // Restore variables
                    [context.U, context.V, context.K] = [SVU, SVV, SVK.slice()];
                }

                // Recompute QP and scalar values to continue the second stage
                [context.A, context.B] = QUADSD(context.U, context.V, P, context.QP);
                TYPE = CALCSC(context);
            }
        }

        OVV = VV;
        OSS = SS;
        OTV = TV;
        OTS = TS;
    }

    return [];
}

/**
 * Variable-shift K-polynomial iteration for a quadratic factor converges
 * only if the zeros are equimodular or nearly so.
 *
 * @param {number} UU - The linear coefficient.
 * @param {number} VV - The constant coefficient.
 * @param {Context} context - The calculation context.
 * @return {Complex[]} The zeros found.
 */
function QUADIT(UU, VV, context) {
    let {P, N, QP} = context;
    let TRIED = false;
    let RELSTP, MP, OMP, ZM, EE, T, TYPE, UI, VI;
    let U = UU;
    let V = VV;
    let J = 0;

    while(true) {
        // Main loop
        let roots = QUAD(1, U, V);
        let [[SZR, SZI], [LZR, LZI]] = roots;
        // Return if roots of the quadratic are real and not close to multiple
        // or nearly equal and of opposite sign
        if(abs(abs(SZR)-abs(LZR)) > 0.01*abs(LZR))
            return [];
        // Evaluate polynomial by quadratic synthetic division
        let [A, B] = QUADSD(U, V, P, QP);
        MP = abs(A-SZR*B) + abs(SZI*B);
        // Compute a rigorous bound on the rounding error in evaluting P
        ZM = sqrt(abs(V));
        EE = QP.slice(0, N).reduce((a, x) => a*ZM + abs(x), 0) + abs(QP[0]);
        T = -SZR*B;
        EE = EE*ZM + abs(A+T);
        EE = (5*MRE+4*ARE) * EE
            - (5*MRE+2*ARE) * (abs(A+T)+abs(B)*ZM)
            + 2*ARE*abs(T);
        // Iteration has converged sufficiently if the polynomial value is less
        // than 20 times this bound
        if(MP <= 20*EE)
            return roots;

        // Stop iteration after 20 steps
        if(++J > 20) return [];

        if(J >= 2 && RELSTP <= .01 && MP >= OMP && !TRIED) {
            // A cluster appears to be stalling the convergence.  Five fixed
            // shift steps are taken with a U,V close to the cluster.
            if(RELSTP < ETA) RELSTP = ETA;
            RELSTP = sqrt(RELSTP);
            U -= U*RELSTP;
            V += V*RELSTP;
            [A, B] = QUADSD(U, V, P, QP);
            Object.assign(context, {A, B, U, V, QP});
            for(let I = 0; I < 5; ++I) {
                TYPE = CALCSC(context);
                NEXTK(TYPE, context);
            }
            TRIED = true;
            J = 0;
        }

        // Calculate next K polynomial and new U and V
        OMP = MP;
        Object.assign(context, {A, B, U, V, QP});
        TYPE = CALCSC(context);
        NEXTK(TYPE, context);
        TYPE = CALCSC(context);
        [UI, VI] = NEWEST(TYPE, context);
        // If VI is zero the iteration is not converging
        if(VI === 0) return [];
        RELSTP = abs((VI-V)/VI);
        U = UI;
        V = VI;
    }
}

/**
 * Variable-shift H polynomial iteration for a real zero.
 *
 * @param {number} SSS - The starting point.
 * @param {Context} context - The calculation context.
 * @return {Array} An array `[S, IFLAG]`, where `IFLAG` is a boolean; if `IFLAG`
 *   is true then `S` is the next start position, and if `IFLAG` is false then
 *   `S` is either the root or is `null` if no root was found.
 */
function REALIT(SSS, context) {
    let {N, NN, P, QP, K, QK} = context;
    let NM1 = N - 1;
    let S = SSS;
    let J = 0;
    let T, OMP, KV, PV, MP, MS, EE;
    while(true) {
        // Main loop
        PV = P.reduce((a, x, i) => (QP[i] = a*S + x), 0);
        MP = abs(PV);
        // Compute a rigorous bound on the error in evaluating P
        MS = abs(S);
        EE = (MRE/(ARE+MRE))*abs(QP[0]);
        EE = QP.slice(1).reduce((a, x) => a*MS + abs(x), EE);

        // Iteration has converged sufficiently if the polynomial value is less
        // than 20 times this bound
        if(MP <= 20*((ARE+MRE)*EE-MRE*MP))
            return [S, false];
        // Stop iteration after 10 steps
        if(++J > 10)
            return [null, false];
        if(J >= 2 && abs(T) <= 0.001*abs(S-T) && MP >= OMP)
            // A cluster of zeros near the real axis has been encountered.
            // Return with IFLAG set to initiate a quadratic iteration.
            return [S, true];

        // Return if the polynomial value has increased significantly
        OMP = MP;
        // Compute T, the next polynomial, and the new iterate
        KV = K.reduce((a, x, i) => (QK[i] = a*S + x), 0);
        if(abs(KV) > abs(K[N-1])*10*ETA) {
            // Use the scaled form of the recurrence if the value of K at S is
            // nonzero
            T = -PV/KV;
            context.K = K
                = QP.slice(0, N).map((x, i) => T*(i > 0 ? QK[i-1] : 0) + x);
        } else {
            // Use unscaled form
            context.K = K = [0, ...QK.slice(0, N-1)];
        }
        KV = K.slice(1).reduce((a, x) => a*S + x, KV);
        T = abs(KV) > abs(K[N-1])*10*ETA ? -PV/KV : 0;
        S += T;
    }
}

/**
 * This routine calculates scalar quantities used to compute the next K
 * polynomial and new estimates of the quadratic coefficients.
 *
 * @param {Context} context - The calculation context.
 * @return {integer} Integer indicating how the calculations are normalized
 *   to avoid overflow.
 */
function CALCSC(context) {
    let {N, U, V, A, B, K, QK} = context;
    // Synthetic division of K by the quadratic 1,U,V
    let [C, D] = [context.C, context.D] = QUADSD(U, V, K, QK);

    if(abs(C) <= abs(K[N-1])*100*ETA && abs(D) <= abs(K[N-2])*100*ETA)
        // type 3 indicates the quadratic is almost a factor of K
        return 3;
    let TYPE;
    let E, F, G, H, A1, A3, A7;
    if(abs(D) >= abs(C)) {
        E = A/D;
        F = C/D;
        G = U*B;
        H = V*B;
        A3 = (A+G)*E + H*(B/D);
        A1 = B*F - A;
        A7 = (F+U)*A + H;
        // type 2 indicates that all formulas are divided by D
        TYPE = 2;
    } else {
        E = A/C;
        F = D/C;
        G = U*E;
        H = V*B;
        A3 = A*E + (H/C+G)*B;
        A1 = B - A*(D/C);
        A7 = A + G*D + H*F;
        // type 1 indicates that all formulas are divided by C
        TYPE = 1;
    }
    Object.assign(context, {F, G, H, A1, A3, A7});
    return TYPE;
}

/**
 * Compute the next K polynomials using scalars computed in {@link CALCSC}.
 *
 * @param TYPE - The return value of {@link CALCSC}.
 * @param {Context} context - The calculation context.
 * @return {undefined}
 */
function NEXTK(TYPE, context) {
    let {K, QK, QP, N, A, B, A1, A3, A7} = context;
    if(TYPE == 3) {
        // Use unscaled form of the recurrence if type is 3
        context.K = [0, 0, ...QK];
        return;
    }
    if(abs(A1) > abs(TYPE === 1 ? B : A) * ETA * 10) {
        // Use scaled form of the recurrence
        A7 /= A1;
        context.A7 = A7;
        A3 /= A1;
        context.A3 = A3;
        K[0] = QP[0];
        K[1] = QP[1] - A7*QP[0];
        for(let I = 2; I < N; ++I)
            K[I] = A3*QK[I-2] - A7*QP[I-1] + QP[I];
        return;
    }
    // If A1 is nearly zero then use a special form of the recurrence
    K[0] = 0;
    K[1] = -A7*QP[1];
    for(let I = 2; I < N; ++I)
        K[I] = A3*QK[I-2] - A7*QP[I-1];
}

/**
 * Compute new estimates of the quadratic coefficients using the scalars
 * computed in {@link CALCSC}.
 *
 * @param TYPE - The return value of {@link CALCSC}.
 * @param {Context} context - The calculation context.
 * @return {number[]} An array `[U, V]` containing the new quadratic
 *   coefficients.
 */
function NEWEST(TYPE, context) {
    if(TYPE === 3)
        // If type==3 the quadratic is zeroed
        return [0, 0];
    // Use formulas appropriate to setting of type.
    let {K, P, A1, A3, A7, A, B, C, D, F, G, H, N, NN, U, V} = context;
    let A4, A5;
    if(TYPE === 2) {
        A4 = (A+G)*F + H;
        A5 = (F+U)*C + V*D;
    } else {
        A4 = A + U*B + H*F;
        A5 = C + (U+V*F)*D;
    }
    // Evaluate new quadratic coefficients.
    let B1 = -K[N-1]/P[NN-1];
    let B2 = -(K[N-1-2]+B1*P[N-1])/P[NN-1];
    let C1 = V*B2*A1;
    let C2 = B1*A7;
    let C3 = B1*B1*A3;
    let C4 = C1 - C2 - C3;
    let TEMP = A5 + B1*A4 - C4;
    if(TEMP === 0)
        return [0, 0];
    let UU = U - (U*(C3+C2)+V*(B1*A1+B2*A7))/TEMP;
    let VV = V*(1.+C4/TEMP);
    return [UU, VV];
}

/**
 * Synthetic division of `P` by the quadratic `Z^2 + U*Z + V`.
 *
 * @param {number} U - The linear coefficient.
 * @param {number} V - The constant coefficient.
 * @param {Polynomial} P - The polynomial to divide.
 * @param {Polynomial} Q - The quotient and remainder are placed here.  The
 *   first `deg P - 2 + 1` terms are the quotient, and the last two terms are
 *   the remainder.
 * @return {number[]} An array `[A, B]` defining the remainder.
 */
function QUADSD(U, V, P, Q) {
    return P.reduce(([A, B], x, i) => [(Q[i] = x - U*A - V*B), A], [0, 0]);
}

/**
 * Calculate the zeros of the quadratic A*Z^2+B1*Z+C.  The quadratic formula,
 * modified to avoid overflow, is used to find the larger zero if the zeros are
 * real and both zeros are complex.  The smaller real zero is found directly
 * from the product of the zeros C/A.
 *
 * Returns `[[SR, SI], [LR, LI]]`, where `[SR, SI]` is the smaller root and
 * `[LR, LI]` is the largero.
 *
 * @param {number} A  - The quadratic coefficient.
 * @param {number} B1 - The linear coefficient.
 * @param {number} C  - The constant coefficient.
 * @return {Complex[]} The roots.
 */
function QUAD(A, B1, C) {
    if(A === 0)
        return [[B1 === 0 ? 0 : -C/B1, 0], [0, 0]];
    if(C === 0)
        return [[0, 0], [-B1/A, 0]];
    // Compute discriminant avoiding overflow
    let B = B1/2, D, E;
    if(abs(B) >= abs(C)) {
        E = 1 - (A/B) * (C/B);
        D = sqrt(abs(E)) * abs(B);
    } else {
        E = C < 0 ? -A : A;
        E = B * (B/abs(C)) - E;
        D = sqrt(abs(E)) * sqrt(abs(C));
    }
    if(E >= 0) {
        // Real zeros
        if(B >= 0) D = -D;
        let LR = (-B+D)/A;
        return [[LR !== 0 ? (C/LR)/A : 0, 0], [LR, 0]];
    }
    // Complex conjugate zeros
    let SR = -B/A, SI = abs(D/A);
    return [[SR, SI], [SR, -SI]];
}
