// Import the standard Faust Libraries
import("stdfaust.lib");

// constant values
EPS = ma.EPSILON;
// EPSILON = 2.2204460492503131e-016;
PI = ma.PI;
// PI = 3.141592653589793238462643383279502884197169399375105820974944592307816406;
TWOPI = 2.0 * PI;
SRMAX = 192000;
SR = ma.SR;
NY = SR / 2.0;
T = 1.0 / SR;
PIT = PI * T;
TWOPIT = TWOPI * T;
MAX = ma.MAX;
// MAX = 1.7976931348623158e+308;
MIN = ma.MIN;
// MIN = 2.2250738585072014e-308;

// Global Trigger
globalTrigger = button("TRIG") <: _ > _';

sequenceLenght = 16;
sequenceFrequency = 8;

// Global Counter
globalCounter = (sequenceFrequency : + ~ _ % (SR * sequenceLenght)) / SR <: _, (_ + 1 : vbargraph("globalCounter [style:numerical]", MIN, MAX)) : attach;
//process = globalCounter;

// random float numbers based on a single noise generator
seeds(N, i, trigger) = (abs(noise(1510147 * i)), trigger) : sah <: par(i, N, hash(i + 1));

// Logic Gates
gates(N) = par(i, N, (_ > i) * (_ < (i + 1)));
counterModN(N) = globalCounter <: gates(N);
gateModN(N) = vecOp((si.bus(N), (globalCounter <: gates(N))), *);
// Single Sequence in Matrix
sequence(i, trigger) = vgroup("Matrix", hgroup("Sequence %i", seeds(sequenceLenght, i, trigger) : infos(sequenceLenght) : gateModN(sequenceLenght) :> _));

// Sine Oscillators Banks
oscillators(N, i, fRange, trigger, x) = vgroup("Frequencies", hgroup(" %i", seeds(N, i, trigger) : par(i, N, _ * fRange : int) : infos(N) : par(i, N, sin(x * _ * TWOPI))));

// Network Voice
Voice(oscsN, i, frequencyRange, filterRange, trigger) = tgroup("Discrete States", (sequence(i, trigger) : phasorOffset(1) : oscillators(oscsN, i, frequencyRange, trigger) :> _ * sequence(i + 1, trigger) / oscsN)); // : HPTPT(sequence(i + 2, trigger) * filterRange) : APTPT(sequence(i + 3, trigger) * filterRange));

process = globalTrigger <: Voice(16, 10, 10000, 10000), Voice(16, 20, 10000, 10000);


// signal info
info(i, x) = x <: _, vbargraph(" %i [style:numerical]", MIN, MAX) : attach;
infos(N) = par(i, N, info(i + 1));

// first derivate
derivate(x) = x < x';

// binary selector
selector(sel, x, y) = x * (1 - sel) + y * (sel);

// a classic sample and hold
sah(x, t) = selector(t, _, x) ~ _;

// classic phasor
phasor(f) = f * (0 + 1') : + ~  _ % SR : _ / SR;

// phasor with phase reset
phasorPH(f, reset) = f * (0 + 1') : + * (1 - reset) ~  _ % SR : _ / SR;

// phasor with phase offset : 0 - 1
phasorOffset(f, offset) = f * (0 + 1') : + ~  _ % SR : (_ + offset * SR) % SR : _ / SR;

// pulsetrain
pulsetrain(f) = (1 - 1') + derivate(phasor(f));

// pseudo-random noise with linear congruential generator (LCG)
noise(initSeed) = lcg ~ _ : (_ / m)
with{
    a = 18446744073709551557; c = 12345; m = 2 ^ 31; 
    lcg(seed) = ((a * seed + c) + (initSeed - initSeed') % m);
};

// hash number generator
hash(i, x) = (x <: (_ * (2.718281828 * i) + _ * (3.141592653 * i)) ^ 1.618033988) % 1;

// perform operations on an arbitrary number of vectors
vecOp(vectorsList, op) =
    vectorsList : seq(i, vecDim - 1, vecOp2D , vecBus(vecDim - 2 - i))
    with {
        vecBus(0) = par(i, vecLen, 0 : !);
        vecBus(dim) = par(i, dim, si.bus(vecLen));
        vecOp2D = ro.interleave(vecLen, 2) : par(i, vecLen, op);
        vecDim = outputs(vectorsList) / vecLen;
        vecLen = outputs(ba.take(1, vectorsList));
    };

// Zavalishin Onepole TPT Filter
onePoleTPT(cf, x) = loop ~ _ : ! , si.bus(3)
    with {
        g = tan(cf * PI * ma.T);
        G = g / (1.0 + g);
        loop(s) = u , lp , hp , ap
            with {
            v = (x - s) * G; u = v + lp; lp = v + s; hp = x - lp; ap = lp - hp;
            };
    };

// Lowpass  TPT
LPTPT(cf, x) = onePoleTPT(cf, x) : (_ , ! , !);

// Highpass TPT
HPTPT(cf, x) = onePoleTPT(cf, x) : (! , _ , !);

// Allpass TPT
APTPT(cf, x) = onePoleTPT(cf, x) : (!, !, _);
