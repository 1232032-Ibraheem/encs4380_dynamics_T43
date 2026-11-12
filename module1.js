/* =====================================================================
 * module1.js — ENCS4380 Homework 1, Dynamic Characteristics Toolkit
 * Ibraheem Mhanna (1232032)
 *
 * This file is loaded first and contains, in order:
 *   1. DynMath       — all mathematics shared by the four modules
 *   2. UI            — shared sliders, outputs, plot settings, Predict-then-check
 *   3. DynValidation — registry for the six validation checks
 *   4. Module 1      — first-order time response (+ validation check 1)
 * Parts 1–3 do not touch the page when loaded, so they also run in Node.
 * ===================================================================== */

/* ---------------------- 1. DynMath (shared mathematics) ---------------------- */
(function (root) {
  'use strict';

  const TWO_PI = 2 * Math.PI;

  /* ------------------------------------------------------------------ *
   * Generic helpers
   * ------------------------------------------------------------------ */
  function linspace(a, b, n) {
    const out = new Array(n);
    const step = n > 1 ? (b - a) / (n - 1) : 0;
    for (let i = 0; i < n; i++) out[i] = a + i * step;
    return out;
  }

  function logspace(a, b, n) { // a, b are base-10 exponents
    return linspace(a, b, n).map(e => Math.pow(10, e));
  }

  function clamp(x, lo, hi) { return Math.min(hi, Math.max(lo, x)); }

  function isFinitePositive(x) { return typeof x === 'number' && isFinite(x) && x > 0; }

  const rad2deg = r => r * 180 / Math.PI;
  const deg2rad = d => d * Math.PI / 180;

  /** Least-squares straight line y = a + b x. Returns {slope, intercept, r2}. */
  function linearFit(xs, ys) {
    const n = xs.length;
    if (n < 2) return { slope: NaN, intercept: NaN, r2: NaN };
    let sx = 0, sy = 0, sxx = 0, sxy = 0, syy = 0;
    for (let i = 0; i < n; i++) {
      sx += xs[i]; sy += ys[i];
      sxx += xs[i] * xs[i]; sxy += xs[i] * ys[i]; syy += ys[i] * ys[i];
    }
    const den = n * sxx - sx * sx;
    const slope = (n * sxy - sx * sy) / den;
    const intercept = (sy - slope * sx) / n;
    const ssTot = syy - sy * sy / n;
    let ssRes = 0;
    for (let i = 0; i < n; i++) {
      const r = ys[i] - (intercept + slope * xs[i]);
      ssRes += r * r;
    }
    return { slope, intercept, r2: ssTot > 0 ? 1 - ssRes / ssTot : 1 };
  }

  /** Linear interpolation of the x where the sampled curve crosses `level`
   *  between samples i-1 and i. */
  function interpCross(t, y, i, level) {
    const y0 = y[i - 1], y1 = y[i];
    if (y1 === y0) return t[i];
    return t[i - 1] + (level - y0) * (t[i] - t[i - 1]) / (y1 - y0);
  }

  /* ------------------------------------------------------------------ *
   * Complex numbers (just what the 2x2 eigenvalue problem needs)
   * ------------------------------------------------------------------ */
  function complex(re, im) { return { re, im: im || 0 }; }

  function formatComplex(z, digits) {
    const d = digits === undefined ? 3 : digits;
    const re = Math.abs(z.re) < 5e-13 ? 0 : z.re;
    const im = Math.abs(z.im) < 5e-13 ? 0 : z.im;
    if (im === 0) return re.toFixed(d);
    const sign = im < 0 ? '−' : '+';
    return `${re.toFixed(d)} ${sign} j${Math.abs(im).toFixed(d)}`;
  }

  function complexDistance(a, b) { return Math.hypot(a.re - b.re, a.im - b.im); }

  /** Eigenvalues of a real 2x2 matrix [[a,b],[c,d]] from its characteristic
   *  polynomial  λ² − tr(A) λ + det(A) = 0.
   *  Returned sorted: larger imaginary part first, then larger real part. */
  function eig2x2(A) {
    const a = A[0][0], b = A[0][1], c = A[1][0], d = A[1][1];
    const tr = a + d;
    const det = a * d - b * c;
    const disc = tr * tr / 4 - det;
    if (disc >= 0) {
      const s = Math.sqrt(disc);
      return [complex(tr / 2 + s, 0), complex(tr / 2 - s, 0)];
    }
    const s = Math.sqrt(-disc);
    return [complex(tr / 2, s), complex(tr / 2, -s)];
  }

  /* ------------------------------------------------------------------ *
   * Random numbers — seeded generator so a round can be reproduced
   * ------------------------------------------------------------------ */
  function mulberry32(seed) {
    let a = seed >>> 0;
    return function () {
      a = (a + 0x6D2B79F5) >>> 0;
      let t = a;
      t = Math.imul(t ^ (t >>> 15), t | 1);
      t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
  }

  /** Standard normal sample (Box–Muller) from a uniform generator. */
  function gaussian(rand) {
    let u = 0, v = 0;
    while (u === 0) u = rand();
    while (v === 0) v = rand();
    return Math.sqrt(-2 * Math.log(u)) * Math.cos(TWO_PI * v);
  }

  function uniform(rand, lo, hi) { return lo + (hi - lo) * rand(); }
  function logUniform(rand, lo, hi) { return Math.exp(uniform(rand, Math.log(lo), Math.log(hi))); }

  /* ------------------------------------------------------------------ *
   * FIRST-ORDER SYSTEM   τ dy/dt + y = K x(t),   G(s) = K / (τs + 1)
   * ------------------------------------------------------------------ */
  const firstOrder = {
    /** Step x(t) = A u(t), zero initial condition:  y = KA(1 − e^{−t/τ}) */
    step(t, K, A, tau) {
      if (t <= 0) return 0;
      return K * A * (1 - Math.exp(-t / tau));
    },

    /** Ramp x(t) = R t:  y = KR [ t − τ (1 − e^{−t/τ}) ] */
    ramp(t, K, R, tau) {
      if (t <= 0) return 0;
      return K * R * (t - tau * (1 - Math.exp(-t / tau)));
    },

    /** Ramp error between ideal output K·x(t) and actual y(t):
     *  e(t) = KRτ (1 − e^{−t/τ})  →  steady-state error e_ss = K R τ */
    rampError(t, K, R, tau) {
      if (t <= 0) return 0;
      return K * R * tau * (1 - Math.exp(-t / tau));
    },

    rampSteadyStateError(K, R, tau) { return K * R * tau; },

    /** Time lag of the ramp response once the transient has died: τ */
    rampLag(tau) { return tau; },

    /** Step error fraction  e_m(t) = (y_ss − y(t)) / y_ss = e^{−t/τ}.
     *  ln(e_m) = −t/τ is a straight line of slope −1/τ. */
    errorFraction(t, tau) { return Math.exp(-t / tau); },

    /** Amplitude ratio M(ω) = 1 / sqrt(1 + (ωτ)²)  (normalised by K) */
    magnitude(w, tau) { return 1 / Math.sqrt(1 + (w * tau) * (w * tau)); },

    /** Phase φ(ω) = −tan⁻¹(ωτ)  [rad]  (negative ⇒ output lags) */
    phase(w, tau) { return -Math.atan(w * tau); },

    /** Dynamic error δ(ω) = M(ω) − 1  (negative ⇒ amplitude is reduced) */
    dynamicError(w, tau) { return firstOrder.magnitude(w, tau) - 1; },

    /** Cut-off (corner, −3 dB) angular frequency ω_c = 1/τ */
    cutoff(tau) { return 1 / tau; },

    /** Largest τ that keeps |δ| ≤ δmax for every frequency up to fMax [Hz].
     *  M(ω) ≥ 1 − δmax  ⇔  ωτ ≤ sqrt(1/(1−δmax)² − 1). Since M decreases
     *  with ω the worst case is ω = 2π fMax. */
    tauMaxForDynamicError(deltaMax, fMaxHz) {
      const Mmin = 1 - deltaMax;
      const wtMax = Math.sqrt(1 / (Mmin * Mmin) - 1);
      return wtMax / (TWO_PI * fMaxHz);
    },

    /** Complete response to x = A sin(ωt) starting from rest:
     *  y = KA/(1+(ωτ)²) [ sin ωt − ωτ cos ωt + ωτ e^{−t/τ} ]
     *  The last term is the start-up transient; the rest equals
     *  KA·M·sin(ωt + φ). */
    sine(t, K, A, w, tau, includeTransient) {
      const a = w * tau;
      const g = K * A / (1 + a * a);
      let y = g * (Math.sin(w * t) - a * Math.cos(w * t));
      if (includeTransient && t >= 0) y += g * a * Math.exp(-t / tau);
      return y;
    },
  };

  /** Measure amplitude ratio and phase shift between two sampled sinusoids
   *  of known angular frequency w, by Fourier projection over the last
   *  `periods` whole periods of the record (the transient has died there).
   *  This is "reading off the curve", independent of the closed form. */
  function measureSinePair(t, x, y, w, periods) {
    const P = periods || 3;
    const T = TWO_PI / w;
    const tEnd = t[t.length - 1];
    const tStart = tEnd - P * T;
    let ax = 0, bx = 0, ay = 0, by = 0, span = 0;
    for (let i = 1; i < t.length; i++) {
      if (t[i - 1] < tStart) continue;
      const dt = t[i] - t[i - 1];
      const s0 = Math.sin(w * t[i - 1]), s1 = Math.sin(w * t[i]);
      const c0 = Math.cos(w * t[i - 1]), c1 = Math.cos(w * t[i]);
      ax += 0.5 * dt * (x[i - 1] * s0 + x[i] * s1);
      bx += 0.5 * dt * (x[i - 1] * c0 + x[i] * c1);
      ay += 0.5 * dt * (y[i - 1] * s0 + y[i] * s1);
      by += 0.5 * dt * (y[i - 1] * c0 + y[i] * c1);
      span += dt;
    }
    const ampX = (2 / span) * Math.hypot(ax, bx);
    const ampY = (2 / span) * Math.hypot(ay, by);
    // signal ≈ a sin + b cos = R sin(wt + atan2(b, a))
    const phX = Math.atan2(bx, ax), phY = Math.atan2(by, ay);
    let dphi = phY - phX;
    while (dphi > Math.PI) dphi -= TWO_PI;
    while (dphi < -Math.PI) dphi += TWO_PI;
    return { ampX, ampY, ratio: ampY / ampX, phase: dphi };
  }

  /* ------------------------------------------------------------------ *
   * SECOND-ORDER SYSTEM   G(s) = K ωn² / (s² + 2ζωn s + ωn²)
   * ------------------------------------------------------------------ */
  const CRIT_TOL = 1e-6;

  const secondOrder = {
    classify(z) {
      if (z < 0) return 'unstable';
      if (z === 0) return 'undamped';
      if (Math.abs(z - 1) < CRIT_TOL) return 'critically damped';
      return z < 1 ? 'underdamped' : 'overdamped';
    },

    /** Poles: s = −ζωn ± ωn sqrt(ζ² − 1)  (complex when ζ < 1) */
    poles(z, wn) {
      const disc = z * z - 1;
      if (disc < 0) {
        const wd = wn * Math.sqrt(-disc);
        return [complex(-z * wn, wd), complex(-z * wn, -wd)];
      }
      const s = wn * Math.sqrt(disc);
      return [complex(-z * wn + s, 0), complex(-z * wn - s, 0)];
    },

    dampedFrequency(z, wn) { return z < 1 ? wn * Math.sqrt(1 - z * z) : 0; },

    /** Unit-step response for x = A u(t), zero initial conditions. */
    step(t, K, A, z, wn) {
      if (t <= 0) return 0;
      const KA = K * A;
      if (Math.abs(z - 1) < CRIT_TOL) {                 // critically damped
        return KA * (1 - (1 + wn * t) * Math.exp(-wn * t));
      }
      if (z < 1) {                                       // under/undamped
        const r = Math.sqrt(1 - z * z);
        const wd = wn * r;
        return KA * (1 - Math.exp(-z * wn * t) *
          (Math.cos(wd * t) + (z / r) * Math.sin(wd * t)));
      }
      // overdamped: two real poles p1 (slow) and p2 (fast)
      const s = Math.sqrt(z * z - 1);
      const p1 = -wn * (z - s), p2 = -wn * (z + s);
      return KA * (1 + (p2 * Math.exp(p1 * t) - p1 * Math.exp(p2 * t)) / (p1 - p2));
    },

    /** Closed-form performance figures (for comparison with the values
     *  measured off the curve). */
    formulas(z, wn) {
      const out = { Mp: 0, tp: null, tr0100: null, ts2: null, ts5: null };
      if (z < 1) {
        const r = Math.sqrt(1 - z * z);
        const wd = wn * r;
        out.Mp = Math.exp(-Math.PI * z / r);             // fraction
        out.tp = Math.PI / wd;
        out.tr0100 = (Math.PI - Math.acos(z)) / wd;
      }
      if (z > 0) {
        // classic envelope estimates, valid for underdamped systems
        out.ts2 = z < 1 ? 4 / (z * wn) : null;
        out.ts5 = z < 1 ? 3 / (z * wn) : null;
      }
      return out;
    },

    /** Reasonable plot window that shows the full transient. */
    timeWindow(z, wn) {
      if (z <= 0) return 12 * TWO_PI / wn;                // undamped: 12 cycles
      let ts;
      if (z < 1) ts = 4 / (z * wn);
      else ts = 4 / (wn * (z - Math.sqrt(Math.max(z * z - 1, 0)))); // slow pole
      let tEnd = 1.6 * ts;
      if (z < 1) tEnd = Math.max(tEnd, 2.5 * Math.PI / (wn * Math.sqrt(1 - z * z)));
      return Math.min(tEnd, 40 * TWO_PI / wn);
    },
  };

  /** Read step-response metrics off a sampled curve, as one would with a
   *  ruler on a plot.  yss = final value; band = settling band (0.02 = 2%).
   *  Rise time 10–90 %, rise time 0–100 % (first crossing of yss), peak
   *  time, percent overshoot and settling time. */
  function measureStep(t, y, yss, band) {
    const n = t.length;
    const res = { tr1090: null, tr0100: null, tp: null, yPeak: null,
                  overshoot: 0, ts: null, settled: false, t10: null, t90: null };
    if (!(Math.abs(yss) > 0)) return res;
    const sgn = Math.sign(yss);
    const yn = y.map(v => v * sgn); const Y = Math.abs(yss);
    // rise times
    for (let i = 1; i < n; i++) {
      if (res.t10 === null && yn[i] >= 0.1 * Y) res.t10 = interpCross(t, yn, i, 0.1 * Y);
      if (res.t90 === null && yn[i] >= 0.9 * Y) res.t90 = interpCross(t, yn, i, 0.9 * Y);
      if (res.tr0100 === null && yn[i] >= Y) { res.tr0100 = interpCross(t, yn, i, Y); }
      if (res.t90 !== null && res.tr0100 !== null) break;
    }
    if (res.t10 !== null && res.t90 !== null) res.tr1090 = res.t90 - res.t10;
    // peak (global maximum), refined by parabolic interpolation
    let iMax = 0;
    for (let i = 1; i < n; i++) if (yn[i] > yn[iMax]) iMax = i;
    let tPk = t[iMax], yPk = yn[iMax];
    if (iMax > 0 && iMax < n - 1) {
      const y0 = yn[iMax - 1], y1 = yn[iMax], y2 = yn[iMax + 1];
      const den = y0 - 2 * y1 + y2;
      if (den < 0) {
        const off = 0.5 * (y0 - y2) / den;            // in samples, |off| ≤ 0.5
        const h = t[iMax + 1] - t[iMax];
        tPk = t[iMax] + off * h;
        yPk = y1 - 0.25 * (y0 - y2) * off;
      }
    }
    // an overshoot exists only if the curve rises above the final value and
    // the maximum is an interior peak (not the last sample of a monotone rise)
    const os = (yPk - Y) / Y;
    if (os > 1e-9 && iMax < n - 1) {
      res.overshoot = os;
      res.tp = tPk;
      res.yPeak = yPk * sgn;
    } else {
      res.overshoot = 0;
      res.tp = null;
      res.yPeak = yn[iMax] * sgn;
    }
    // settling: last time the curve is outside ±band·yss
    const lim = band * Y;
    let last = -1;
    for (let i = n - 1; i >= 0; i--) {
      if (Math.abs(yn[i] - Y) > lim) { last = i; break; }
    }
    if (last === -1) { res.ts = 0; res.settled = true; }
    else if (last < n - 1) {
      // interpolate exact exit time of the band between last and last+1
      const lvl = yn[last] > Y ? Y + lim : Y - lim;
      res.ts = interpCross(t, yn, last + 1, lvl);
      res.settled = true;
    }
    return res;
  }

  /* ------------------------------------------------------------------ *
   * STATE SPACE — mass-spring-damper (Example 3 of the SS notes)
   *   x1 = displacement, x2 = velocity
   *   ẋ = [[0,1],[−k/m, −b/m]] x + [0, 1/m]ᵀ F,   y = [1 0] x
   * ------------------------------------------------------------------ */
  const stateSpace = {
    msd(m, b, k) {
      return {
        A: [[0, 1], [-k / m, -b / m]],
        B: [[0], [1 / m]],
        C: [[1, 0]],
        D: [[0]],
      };
    },
    /** Equivalent second-order parameters of the same plant:
     *  X(s)/F(s) = (1/m)/(s² + (b/m)s + k/m) = (1/k)·ωn²/(s²+2ζωn s+ωn²) */
    msdToSecondOrder(m, b, k) {
      const wn = Math.sqrt(k / m);
      return { wn, zeta: b / (2 * Math.sqrt(k * m)), K: 1 / k };
    },
    secondOrderToMsd(m, zeta, wn) {
      const k = m * wn * wn;
      return { m, k, b: 2 * zeta * Math.sqrt(k * m) };
    },
  };

  /* ------------------------------------------------------------------ *
   * ODE solvers (used for the state-space simulation and the bonus page)
   * f(t, x) returns dx/dt as an array.
   * ------------------------------------------------------------------ */
  function eulerStep(f, t, x, h) {
    const k1 = f(t, x);
    return x.map((xi, i) => xi + h * k1[i]);
  }

  function rk4Step(f, t, x, h) {
    const k1 = f(t, x);
    const k2 = f(t + h / 2, x.map((xi, i) => xi + h / 2 * k1[i]));
    const k3 = f(t + h / 2, x.map((xi, i) => xi + h / 2 * k2[i]));
    const k4 = f(t + h, x.map((xi, i) => xi + h * k3[i]));
    return x.map((xi, i) => xi + h / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]));
  }

  /** Integrate from t0 to tEnd with fixed step h. method: 'rk4' | 'euler'. */
  function integrate(f, x0, t0, tEnd, h, method) {
    const stepFn = method === 'euler' ? eulerStep : rk4Step;
    const nSteps = Math.max(1, Math.ceil((tEnd - t0) / h - 1e-9));
    const t = [t0], X = [x0.slice()];
    let x = x0.slice(), tc = t0;
    for (let i = 0; i < nSteps; i++) {
      const hh = Math.min(h, tEnd - tc);
      x = stepFn(f, tc, x, hh);
      tc += hh;
      t.push(tc); X.push(x);
    }
    return { t, X };
  }

  /** Linear state-space model ẋ = A x + B u(t) as an ODE right-hand side. */
  function linearSystem(A, B, u) {
    return function (t, x) {
      const n = x.length, out = new Array(n);
      const ut = u(t);
      for (let i = 0; i < n; i++) {
        let s = B[i][0] * ut;
        for (let j = 0; j < n; j++) s += A[i][j] * x[j];
        out[i] = s;
      }
      return out;
    };
  }

  const api = {
    TWO_PI, linspace, logspace, clamp, isFinitePositive, rad2deg, deg2rad,
    linearFit, interpCross, complex, formatComplex, complexDistance, eig2x2,
    mulberry32, gaussian, uniform, logUniform,
    firstOrder, secondOrder, stateSpace, measureSinePair, measureStep,
    eulerStep, rk4Step, integrate, linearSystem,
  };

  root.DynMath = api;
})(typeof window !== 'undefined' ? window : globalThis);

/* ---------------------- 2. UI (shared interface helpers) ---------------------- */
(function (root) {
  'use strict';

  const $ = id => document.getElementById(id);

  function el(tag, attrs, html) {
    const e = document.createElement(tag);
    if (attrs) Object.entries(attrs).forEach(([k, v]) => {
      if (k === 'class') e.className = v; else if (v !== undefined && v !== null) e.setAttribute(k, v);
    });
    if (html !== undefined) e.innerHTML = html;
    return e;
  }

  /** Format a number to `sig` significant figures without trailing noise. */
  function fmt(v, sig) {
    if (v === null || v === undefined || Number.isNaN(v)) return '—';
    if (!isFinite(v)) return v > 0 ? '∞' : '−∞';
    const s = sig || 4;
    if (v === 0) return '0';
    const a = Math.abs(v);
    let out;
    if (a >= 1e6 || a < 1e-4) out = v.toExponential(s - 1);
    else out = String(Number(v.toPrecision(s)));
    return out.replace('-', '−');
  }
  function fixed(v, d) {
    if (v === null || v === undefined || Number.isNaN(v)) return '—';
    if (!isFinite(v)) return '∞';
    const r = v.toFixed(d);
    return (/^-0\.?0*$/.test(r) ? r.slice(1) : r).replace('-', '−');
  }

  /** Coalesce rapid slider events into one redraw per animation frame. */
  function scheduler(fn) {
    let pending = false;
    return function () {
      if (pending) return;
      pending = true;
      requestAnimationFrame(() => { pending = false; fn(); });
    };
  }

  let uid = 0;

  /** Range slider + number box, kept in sync, with clamping and validation.
   *  opts: {label, unit, min, max, step, value, scale:'lin'|'log', help, digits} */
  function slider(container, opts) {
    const id = opts.id || ('s' + (++uid));
    const log = opts.scale === 'log';
    const wrap = el('div', { class: 'ctl' });
    wrap.innerHTML = `
      <div class="ctl-top"><label for="${id}-num">${opts.label}</label><span class="unit">${opts.unit || ''}</span></div>
      <div class="ctl-row">
        <input type="range" id="${id}" aria-label="${opts.label.replace(/<[^>]+>/g, '')} slider">
        <input type="number" class="num" id="${id}-num" inputmode="decimal" step="any">
      </div>
      ${opts.help ? `<div class="help">${opts.help}</div>` : ''}
      <div class="err" role="status" aria-live="polite"></div>`;
    container.appendChild(wrap);
    const range = wrap.querySelector('input[type=range]');
    const num = wrap.querySelector('.num');
    const err = wrap.querySelector('.err');
    let min = opts.min, max = opts.max;
    const digits = opts.digits || 4;
    let value = opts.value;
    const listeners = [];

    function setRangeBounds() {
      if (log) {
        range.min = Math.log10(min); range.max = Math.log10(max); range.step = 0.001;
      } else {
        range.min = min; range.max = max; range.step = opts.step || 'any';
      }
      num.min = min; num.max = max;
    }
    function round(v) {
      if (!log && opts.step) {
        const k = Math.round((v - min) / opts.step);
        v = min + k * opts.step;
        const dec = (String(opts.step).split('.')[1] || '').length;
        return Number(v.toFixed(dec));
      }
      return Number(v.toPrecision(digits));
    }
    function paint() {
      range.value = log ? Math.log10(value) : value;
      num.value = value;
    }
    function emit() { listeners.forEach(f => f(value)); }

    range.addEventListener('input', () => {
      const r = parseFloat(range.value);
      value = round(log ? Math.pow(10, r) : r);
      value = Math.min(max, Math.max(min, value));
      num.value = value; num.classList.remove('is-invalid'); err.textContent = '';
      emit();
    });
    function fromNum(commit) {
      const raw = num.value.trim();
      const v = parseFloat(raw);
      if (raw === '' || !isFinite(v)) {
        num.classList.add('is-invalid');
        err.textContent = `Enter a number between ${fmt(min)} and ${fmt(max)}.`;
        if (commit) { paint(); num.classList.remove('is-invalid'); }
        return;
      }
      if (v < min || v > max) {
        num.classList.add('is-invalid');
        err.textContent = `Allowed range is ${fmt(min)} to ${fmt(max)}; using the nearest limit.`;
        if (!commit) return;
      } else { num.classList.remove('is-invalid'); err.textContent = ''; }
      value = Math.min(max, Math.max(min, v));
      if (commit) { paint(); num.classList.remove('is-invalid'); }
      else range.value = log ? Math.log10(value) : value;
      emit();
    }
    num.addEventListener('input', () => fromNum(false));
    num.addEventListener('change', () => fromNum(true));

    setRangeBounds(); paint();
    return {
      id, el: wrap,
      get: () => value,
      set(v, silent) {
        if (!isFinite(v)) return;
        value = Math.min(max, Math.max(min, v));
        if (!log && !opts.step) value = Number(value.toPrecision(digits));
        paint(); num.classList.remove('is-invalid'); err.textContent = '';
        if (!silent) emit();
      },
      setBounds(lo, hi) { min = lo; max = hi; setRangeBounds(); value = Math.min(max, Math.max(min, value)); paint(); },
      onChange(f) { listeners.push(f); return this; },
      show(on) { wrap.hidden = !on; },
    };
  }

  /** Segmented toggle (radio group). options: [{value, label}] */
  function segmented(container, opts) {
    const wrap = el('div', { class: 'ctl' });
    const gid = 'g' + (++uid);
    wrap.innerHTML = `<span class="lbl" id="${gid}">${opts.label}</span>
      <div class="seg" role="group" aria-labelledby="${gid}">
        ${opts.options.map(o => `<button type="button" data-v="${o.value}">${o.label}</button>`).join('')}
      </div>`;
    container.appendChild(wrap);
    let value = opts.value;
    const listeners = [];
    const btns = Array.from(wrap.querySelectorAll('button'));
    function paint() { btns.forEach(b => b.setAttribute('aria-pressed', String(b.dataset.v === String(value)))); }
    btns.forEach(b => b.addEventListener('click', () => {
      if (b.dataset.v === String(value)) return;
      value = b.dataset.v; paint(); listeners.forEach(f => f(value));
    }));
    paint();
    return {
      el: wrap, get: () => value,
      set(v, silent) { value = String(v); paint(); if (!silent) listeners.forEach(f => f(value)); },
      onChange(f) { listeners.push(f); return this; },
    };
  }

  /** Readout grid. rows: [{k, v, h, emph}] — v may contain HTML. */
  function readouts(container, rows) {
    container.innerHTML = rows.map(r =>
      `<div class="ro${r.emph ? ' emph' : ''}"><div class="k">${r.k}</div><div class="v">${r.v}</div>${r.h ? `<div class="h">${r.h}</div>` : ''}</div>`
    ).join('');
  }

  /* -------------------------- KaTeX ---------------------------------- */
  function tex(target, latex, display) {
    if (root.katex) {
      try { root.katex.render(latex, target, { displayMode: !!display, throwOnError: false }); return; }
      catch (e) { /* fall through */ }
    }
    target.textContent = latex;
  }
  function renderMathIn(node) {
    if (root.renderMathInElement) {
      root.renderMathInElement(node, {
        delimiters: [{ left: '\\[', right: '\\]', display: true }, { left: '\\(', right: '\\)', display: false }],
        throwOnError: false,
      });
    }
  }

  /* -------------------------- Plotly --------------------------------- */
  // one fixed dark palette, shared by every plot (matches css/style.css)
  const PALETTE = {
    ink: '#e8eaed', ink2: '#c3c9d2', ink3: '#8b94a1', rule: '#2f353f', surface: '#1d2128',
    in: '#9aa4b1', out: '#4fa3e0', third: '#b48cf0', mark: '#f2b347', markSoft: 'rgba(242,179,71,0.14)',
    pole: '#f0706a', accent: '#4fa3e0', pass: '#5cc98a', fail: '#f0706a',
  };
  function colors() { return PALETTE; }
  function axis(c, title, extra) {
    return Object.assign({
      title: { text: title, font: { size: 12, color: c.ink2 }, standoff: 6 },
      gridcolor: c.rule, zerolinecolor: c.ink3, linecolor: c.ink3, tickfont: { size: 11, color: c.ink2 },
      automargin: true, showline: true, mirror: false, zeroline: true, ticks: 'outside', tickcolor: c.ink3,
    }, extra || {});
  }
  function layout(c, opts) {
    const o = opts || {};
    return Object.assign({
      paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(0,0,0,0)',
      font: { family: 'system-ui, -apple-system, Segoe UI, Roboto, Arial, sans-serif', color: c.ink, size: 12 },
      margin: { l: 56, r: 16, t: 10, b: 46 },
      hovermode: 'closest',
      showlegend: true,
      legend: { orientation: 'h', x: 0, y: 1.02, yanchor: 'bottom', font: { size: 11, color: c.ink2 }, bgcolor: 'rgba(0,0,0,0)' },
      xaxis: axis(c, o.x || '', o.xaxis), yaxis: axis(c, o.y || '', o.yaxis),
    }, o.extra || {});
  }
  const plotConfig = {
    responsive: true, displaylogo: false, scrollZoom: false,
    modeBarButtonsToRemove: ['select2d', 'lasso2d', 'autoScale2d', 'toggleSpikelines', 'hoverCompareCartesian', 'hoverClosestCartesian'],
    toImageButtonOptions: { format: 'png', scale: 2 },
  };
  function plot(target, traces, lay) {
    if (!root.Plotly) { target.textContent = 'Plotting library failed to load.'; return; }
    root.Plotly.react(target, traces, lay, plotConfig);
  }
  function resizeAll(scope) {
    if (!root.Plotly) return;
    (scope || document).querySelectorAll('.js-plotly-plot').forEach(p => {
      if (p.offsetParent !== null) root.Plotly.Plots.resize(p);
    });
  }

  /* ---------------- "Predict, then check" component -------------------
   * questions: array of functions (ctx) => {
   *   text, kind:'choice'|'number', options:[...], answer (index | number),
   *   tol (number kind, relative), unit, explain, tryLabel, apply()  }
   * Questions are generated from the current slider values, so the
   * prediction is about what the user is about to see.  */
  function predict(container, questions) {
    let qi = 0, chosen = null, checked = false, q = null, nRight = 0, nDone = 0;
    container.classList.add('predict');
    container.innerHTML = `
      <h3>Predict, then check <span class="counter" aria-live="polite"></span></h3>
      <div class="q"></div><div class="body"></div>
      <div class="btn-row"><button type="button" class="btn small js-check">Check my prediction</button>
      <button type="button" class="btn small ghost js-next">Next question</button></div>
      <div class="result" aria-live="polite"></div>`;
    const qEl = container.querySelector('.q'), body = container.querySelector('.body');
    const res = container.querySelector('.result'), counter = container.querySelector('.counter');
    const checkBtn = container.querySelector('.js-check');

    function load() {
      q = questions[qi % questions.length]();
      chosen = null; checked = false; res.innerHTML = '';
      qEl.innerHTML = q.text;
      if (q.kind === 'number') {
        body.innerHTML = `<div class="numans"><input class="num" type="number" step="any" inputmode="decimal" aria-label="Your prediction"> <span class="muted">${q.unit || ''}</span></div>`;
      } else {
        body.innerHTML = `<div class="opts">${q.options.map((o, i) => `<button type="button" class="opt" data-i="${i}" aria-pressed="false">${o}</button>`).join('')}</div>`;
        body.querySelectorAll('.opt').forEach(b => b.addEventListener('click', () => {
          if (checked) return;
          chosen = +b.dataset.i;
          body.querySelectorAll('.opt').forEach(x => x.setAttribute('aria-pressed', String(x === b)));
        }));
      }
      checkBtn.disabled = false;
      counter.textContent = nDone ? `${nRight} of ${nDone} right so far` : `Question ${qi % questions.length + 1} of ${questions.length}`;
      renderMathIn(container);
    }
    checkBtn.addEventListener('click', () => {
      if (checked) return;
      let right;
      if (q.kind === 'number') {
        const v = parseFloat(body.querySelector('input').value);
        if (!isFinite(v)) { res.innerHTML = '<p>Type a number first, then check.</p>'; return; }
        right = Math.abs(v - q.answer) <= Math.abs(q.answer) * (q.tol || 0.05) + 1e-12;
        res.innerHTML = `<p><span class="badge ${right ? 'badge-pass' : 'badge-fail'}">${right ? 'Good prediction' : 'Not quite'}</span>
          You said ${fmt(v)} ${q.unit || ''}; the toolkit computes <strong>${fmt(q.answer)} ${q.unit || ''}</strong>.</p>`;
      } else {
        if (chosen === null) { res.innerHTML = '<p>Pick an option first, then check.</p>'; return; }
        right = chosen === q.answer;
        body.querySelectorAll('.opt').forEach((b, i) => {
          if (i === q.answer) b.classList.add('correct');
          else if (i === chosen) b.classList.add('wrong');
        });
        res.innerHTML = `<p><span class="badge ${right ? 'badge-pass' : 'badge-fail'}">${right ? 'Good prediction' : 'Not quite'}</span>
          The answer is <strong>${q.options[q.answer]}</strong>.</p>`;
      }
      checked = true; nDone++; if (right) nRight++;
      res.innerHTML += `<p>${q.explain}</p>`;
      if (q.apply) {
        const b = el('button', { type: 'button', class: 'btn small ghost' }, q.tryLabel || 'Show it on the plot');
        b.addEventListener('click', () => { q.apply(); b.disabled = true; b.textContent = 'Applied — look at the plot'; });
        res.appendChild(b);
      }
      counter.textContent = `${nRight} of ${nDone} right so far`;
      checkBtn.disabled = true;
      renderMathIn(res);
    });
    container.querySelector('.js-next').addEventListener('click', () => { qi++; load(); });
    load();
    return { reload: load };
  }

  root.UI = { $, el, fmt, fixed, scheduler, slider, segmented, readouts, tex, renderMathIn,
              colors, layout, axis, plot, resizeAll, predict };
})(typeof window !== 'undefined' ? window : globalThis);

/* ---------------------- 3. Validation registry ---------------------- */
(function (root) {
  'use strict';
  const tests = [];

  function check(quantity, expected, computed, tol, unit, note) {
    return { quantity, expected, computed, tol, unit: unit || '', note: note || '',
             pass: Math.abs(computed - expected) <= tol };
  }
  /** sample fn on [0, tEnd] with n points */
  function sampled(fn, tEnd, n) {
    const t = root.DynMath.linspace(0, tEnd, n);
    return { t, y: t.map(fn) };
  }
  /** value of a sampled curve at time tq (linear interpolation) */
  function readAt(t, y, tq) {
    for (let i = 1; i < t.length; i++) {
      if (t[i] >= tq) return y[i - 1] + (tq - t[i - 1]) / (t[i] - t[i - 1]) * (y[i] - y[i - 1]);
    }
    return y[y.length - 1];
  }
  /** each module registers its own checks: {id, module, title, description, inputs, run() -> [check]} */
  function add(test) { tests.push(test); }
  function runAll() {
    return tests.slice().sort((a, b) => a.id - b.id).map(t => {
      const checks = t.run();
      return Object.assign({}, t, { checks, pass: checks.every(c => c.pass) });
    });
  }
  function fmtV(v) {
    if (typeof v !== 'number' || !isFinite(v)) return String(v);
    if (v === 0) return '0';
    const a = Math.abs(v);
    if (a < 1e-4) return v.toExponential(2);
    if (a >= 100) return v.toFixed(2);
    return Number(v.toPrecision(5)).toString();
  }
  function render(container, summaryEl) {
    const results = runAll();
    const nPass = results.filter(r => r.pass).length;
    container.innerHTML = results.map(r => `
      <article class="vtest ${r.pass ? 'is-pass' : 'is-fail'}">
        <header class="vtest-head">
          <span class="vtest-num">Test ${r.id}</span><h3>${r.title}</h3>
          <span class="badge ${r.pass ? 'badge-pass' : 'badge-fail'}">${r.pass ? 'PASS' : 'FAIL'}</span>
        </header>
        <p class="vtest-desc">${r.description}</p>
        <p class="vtest-inputs"><strong>Inputs:</strong> ${r.inputs} <span class="muted">(${r.module})</span></p>
        <div class="table-wrap"><table>
          <thead><tr><th>Quantity</th><th class="num">Expected</th><th class="num">Calculated</th><th class="num">Tolerance</th><th>Result</th></tr></thead>
          <tbody>${r.checks.map(c => `
            <tr><td>${c.quantity}${c.note ? `<div class="muted small">${c.note}</div>` : ''}</td>
              <td class="num">${fmtV(c.expected)} ${c.unit}</td><td class="num">${fmtV(c.computed)} ${c.unit}</td>
              <td class="num">± ${fmtV(c.tol)} ${c.unit}</td>
              <td><span class="badge ${c.pass ? 'badge-pass' : 'badge-fail'}">${c.pass ? 'PASS' : 'FAIL'}</span></td></tr>`).join('')}
          </tbody></table></div>
      </article>`).join('');
    if (summaryEl) {
      summaryEl.textContent = `${nPass} / ${results.length} tests pass`;
      summaryEl.className = 'badge ' + (nPass === results.length ? 'badge-pass' : 'badge-fail');
    }
    return results;
  }
  root.DynValidation = { add, check, sampled, readAt, runAll, render };
})(typeof window !== 'undefined' ? window : globalThis);

/* ---------------------- 4. Module 1: first-order time response ---------------------- */
(function (root) {
  'use strict';
  const M = root.DynMath, FO = M.firstOrder;
  const { $, fmt, fixed, slider, segmented, readouts, layout, colors, plot } = root.UI;

  const m1 = {
    ready: false,
    init() {
      const box = $('m1-controls');
      box.insertAdjacentHTML('beforeend', '<h3>Parameters</h3>');
      this.input = segmented(box, { label: 'Input', value: 'step', options: [{ value: 'step', label: 'Step' }, { value: 'ramp', label: 'Ramp' }] });
      this.K = slider(box, { label: 'Static sensitivity K', min: 0.1, max: 10, step: 0.1, value: 1 });
      this.tau = slider(box, { label: 'Time constant τ', unit: 's', min: 0.05, max: 50, scale: 'log', value: 2, digits: 3 });
      this.A = slider(box, { label: 'Step height A', min: 0.1, max: 100, step: 0.1, value: 1 });
      this.R = slider(box, { label: 'Ramp rate R', unit: 'per s', min: 0.1, max: 20, step: 0.1, value: 1 });
      const redraw = root.UI.scheduler(() => this.render());
      [this.input, this.K, this.tau, this.A, this.R].forEach(s => s.onChange(redraw));
      this.setupPredict();
      this.ready = true;
    },

    render() {
      const isStep = this.input.get() === 'step';
      this.A.show(isStep);
      this.R.show(!isStep);
      if (isStep) this.renderStep(this.K.get(), this.tau.get(), this.A.get());
      else this.renderRamp(this.K.get(), this.tau.get(), this.R.get());
    },

    renderStep(K, tau, A) {
      const c = colors(), KA = K * A, tEnd = 6 * tau;
      const t = M.linspace(0, tEnd, 1201);
      const y = t.map(x => FO.step(x, K, A, tau));
      const y1 = FO.step(tau, K, A, tau), y3 = FO.step(3 * tau, K, A, tau);

      $('m1-resp-title').textContent = 'Step response';
      $('m1-resp-cap').textContent = 'Dashed: ideal output KA. Solid: sensor output y(t). Markers at t = τ (63.2 %) and t = 3τ (95 %).';
      const shapes = [];
      [[tau, y1], [3 * tau, y3]].forEach(([tx, ty]) => {
        shapes.push({ type: 'line', x0: tx, x1: tx, y0: 0, y1: ty, line: { color: c.mark, width: 1, dash: 'dot' } });
        shapes.push({ type: 'line', x0: 0, x1: tx, y0: ty, y1: ty, line: { color: c.mark, width: 1, dash: 'dot' } });
      });
      plot($('m1-resp'), [
        { x: [0, 0, tEnd], y: [0, KA, KA], mode: 'lines', name: 'Ideal K·x(t)', line: { color: c.in, dash: 'dash' }, hoverinfo: 'skip' },
        { x: t, y, mode: 'lines', name: 'Output y(t)', line: { color: c.out, width: 3 }, hovertemplate: 't = %{x:.4g} s<br>y = %{y:.4g}<extra></extra>' },
        { x: [tau, 3 * tau], y: [y1, y3], mode: 'markers+text', name: 'At τ and 3τ', text: ['63.2 % (t = τ)', '95 % (t = 3τ)'],
          textposition: 'bottom right', textfont: { color: c.mark }, marker: { color: c.mark, size: 9 }, hovertemplate: 't = %{x:.4g} s<br>y = %{y:.4g}<extra></extra>' },
      ], layout(c, { x: 'time t (s)', y: 'output y(t)', yaxis: { range: [0, KA * 1.12] }, xaxis: { range: [0, tEnd] }, extra: { shapes } }));

      // error fraction e_m = (KA − y)/KA, computed from the plotted curve
      const em = y.map(v => (KA - v) / KA);
      const fit = this.renderErr(t, em, tau, 'e_m = (KA − y)/KA');
      $('m1-err-cap').textContent = 'e_m = e^(−t/τ) is a straight line on a log axis. The dashed line is a least-squares fit of ln e_m against t; its slope is −1/τ.';

      readouts($('m1-readouts'), [
        { k: 'Final value KA', v: fmt(KA, 4) },
        { k: 'y(τ)', v: `${fmt(y1, 4)} (${fixed(100 * y1 / KA, 1)} %)`, emph: true },
        { k: 'y(3τ)', v: `${fmt(y3, 4)} (${fixed(100 * y3 / KA, 1)} %)`, emph: true },
        { k: 'Slope of ln e_m (fitted)', v: `${fmt(fit.slope, 4)} s⁻¹`, h: `−1/τ = ${fmt(-1 / tau, 4)} s⁻¹` },
        { k: 'Steady-state error', v: '0', h: 'the output ends exactly at KA' },
      ]);
      $('m1-interpret').innerHTML = `
        <li>The response reaches 63.2 % at t = τ = ${fmt(tau, 3)} s and 95 % at 3τ = ${fmt(3 * tau, 3)} s. Changing K changes the final value, not these times.</li>
        <li>On the log plot the error is a straight line: every τ seconds it shrinks by a factor e ≈ 2.718.</li>`;
    },

    renderRamp(K, tau, R) {
      const c = colors(), tEnd = 8 * tau;
      const t = M.linspace(0, tEnd, 1201);
      const ideal = t.map(x => K * R * x);
      const y = t.map(x => FO.ramp(x, K, R, tau));
      const ess = FO.rampSteadyStateError(K, R, tau);
      // lag read off the curve at t = 6τ: when did the ideal line have the same value?
      const ta = 6 * tau, ya = FO.ramp(ta, K, R, tau), tIdeal = ya / (K * R);
      const xa = 0.97 * tEnd, ylo = FO.ramp(xa, K, R, tau), yhi = K * R * xa;

      $('m1-resp-title').textContent = 'Ramp response';
      $('m1-resp-cap').textContent = 'Dashed: ideal output K·R·t. Solid: sensor output y(t). After the transient the output is the same line shifted τ to the right.';
      const arrow = (x0, y0, x1, y1) => ({ x: x0, y: y0, ax: x1, ay: y1, xref: 'x', yref: 'y', axref: 'x', ayref: 'y',
        showarrow: true, arrowhead: 2, arrowside: 'end+start', arrowcolor: c.mark, arrowwidth: 1.5, text: '' });
      plot($('m1-resp'), [
        { x: t, y: ideal, mode: 'lines', name: 'Ideal K·x(t)', line: { color: c.in, dash: 'dash' }, hovertemplate: 'K·x = %{y:.4g}<extra></extra>' },
        { x: t, y, mode: 'lines', name: 'Output y(t)', line: { color: c.out, width: 3 }, hovertemplate: 't = %{x:.4g} s<br>y = %{y:.4g}<extra></extra>' },
      ], layout(c, { x: 'time t (s)', y: 'output y(t)', xaxis: { range: [0, tEnd] }, yaxis: { range: [0, K * R * tEnd * 1.04] },
        extra: { annotations: [
          arrow(tIdeal, ya, ta, ya),
          { x: (tIdeal + ta) / 2, y: ya, text: `lag = τ = ${fmt(tau, 3)} s`, showarrow: false, yshift: -16, xshift: 10, bgcolor: c.surface, font: { color: c.mark } },
          arrow(xa, ylo, xa, yhi),
          { x: xa, y: (ylo + yhi) / 2, text: `e_ss = KRτ = ${fmt(ess, 3)}`, showarrow: false, xanchor: 'right', xshift: -8, bgcolor: c.surface, font: { color: c.mark } },
        ] } }));

      // the part of the error still to come: (e_ss − e)/e_ss = e^(−t/τ)
      const em = t.map((x, i) => (ess - (ideal[i] - y[i])) / ess);
      const fit = this.renderErr(t, em, tau, '(e_ss − e)/e_ss');
      $('m1-err-cap').textContent = 'The ramp error e = K·x − y grows to e_ss = KRτ. The remaining part (e_ss − e)/e_ss = e^(−t/τ) is again a straight line of slope −1/τ.';

      readouts($('m1-readouts'), [
        { k: 'Time lag', v: `${fmt(tau, 4)} s`, emph: true, h: `read off the curve: ${fmt(ta - tIdeal, 4)} s` },
        { k: 'Steady-state error e_ss = KRτ', v: fmt(ess, 4), emph: true },
        { k: 'Slope of ln error (fitted)', v: `${fmt(fit.slope, 4)} s⁻¹`, h: `−1/τ = ${fmt(-1 / tau, 4)} s⁻¹` },
      ]);
      $('m1-interpret').innerHTML = `
        <li>After the transient the sensor shows the value the input had τ = ${fmt(tau, 3)} s earlier.</li>
        <li>The steady-state error K·R·τ grows with the ramp rate and with τ.</li>`;
    },

    /** draws the error-fraction plot and returns the least-squares fit of ln(e) against t */
    renderErr(t, em, tau, yTitle) {
      const c = colors();
      const idx = t.map((_, i) => i).filter(i => em[i] > 1e-6 && t[i] > 0);
      const fit = M.linearFit(idx.map(i => t[i]), idx.map(i => Math.log(em[i])));
      const t0 = t[1], t1 = t[t.length - 1];
      plot($('m1-err'), [
        { x: t.slice(1), y: em.slice(1), mode: 'lines', name: 'Error fraction', line: { color: c.out, width: 3 }, hovertemplate: 't = %{x:.4g} s<br>e = %{y:.4g}<extra></extra>' },
        { x: [t0, t1], y: [Math.exp(fit.intercept + fit.slope * t0), Math.exp(fit.intercept + fit.slope * t1)], mode: 'lines',
          name: `Fit: slope ${fmt(fit.slope, 3)} s⁻¹ (−1/τ = ${fmt(-1 / tau, 3)})`, line: { color: c.mark, dash: 'dash' }, hoverinfo: 'skip' },
      ], layout(c, { x: 'time t (s)', y: yTitle, xaxis: { range: [0, t1] },
        yaxis: { type: 'log', range: [Math.log10(em[em.length - 1]) - 0.2, 0.08], exponentformat: 'power' } }));
      return fit;
    },

    setupPredict() {
      const self = this;
      root.UI.predict($('m1-predict'), [
        () => {
          const tau = self.tau.get();
          return { kind: 'choice', text: `A step is applied to a sensor with τ = ${fmt(tau, 3)} s. At t = ${fmt(tau, 3)} s, what fraction of the final value will it read?`,
            options: ['50 %', '63.2 %', '86.5 %', '95 %'], answer: 1,
            explain: 'y(τ) = KA(1 − e<sup>−1</sup>) = 0.632 KA for any K, A and τ.',
            tryLabel: 'Show the step response', apply: () => self.input.set('step') };
        },
        () => {
          const tau = self.tau.get(), t2 = Math.min(50, 2 * tau);
          return { kind: 'choice', text: `If τ is doubled from ${fmt(tau, 3)} s to ${fmt(t2, 3)} s, the time to reach 95 % will…`,
            options: ['halve', 'stay the same', 'double', 'quadruple'], answer: 2,
            explain: `The 95 % point is at 3τ, so it moves from ${fmt(3 * tau, 3)} s to ${fmt(3 * t2, 3)} s.`,
            tryLabel: `Set τ = ${fmt(t2, 3)} s`, apply: () => { self.input.set('step'); self.tau.set(t2); } };
        },
        () => {
          const K = self.K.get(), tau = self.tau.get(), R = self.R.get();
          return { kind: 'number', text: `Ramp input with R = ${fmt(R, 3)} per s, K = ${fmt(K, 3)}, τ = ${fmt(tau, 3)} s. What is the steady-state error?`,
            answer: K * R * tau, tol: 0.03, unit: '',
            explain: 'e<sub>ss</sub> = K·R·τ: the output is τ seconds late, and the ideal output rises by K·R·τ in that time.',
            tryLabel: 'Show the ramp response', apply: () => self.input.set('ramp') };
        },
      ]);
    },
  };


  /* validation check 1 (Module 1) */
  root.DynValidation.add({ id: 1, module: 'Module 1', title: 'First-order step response',
    description: 'At t = τ the output must be 63.2 % of its final value KA, and at t = 3τ about 95 %.',
    inputs: '(K = 1, A = 1, τ = 1 s) and (K = 2.5, A = 4, τ = 10 s)',
    run() {
      const V = root.DynValidation, out = [];
      [{ K: 1, A: 1, tau: 1 }, { K: 2.5, A: 4, tau: 10 }].forEach(s => {
        const c = V.sampled(t => FO.step(t, s.K, s.A, s.tau), 6 * s.tau, 2401), KA = s.K * s.A;
        const lbl = `K = ${s.K}, A = ${s.A}, τ = ${s.tau} s`;
        out.push(V.check(`y(τ)/KA (${lbl})`, 0.632, V.readAt(c.t, c.y, s.tau) / KA, 0.0005, '', 'read off the plotted curve'));
        out.push(V.check(`y(3τ)/KA (${lbl})`, 0.95, V.readAt(c.t, c.y, 3 * s.tau) / KA, 0.001, '', 'read off the plotted curve'));
      });
      return out;
    } });

  root.Modules = root.Modules || {};
  root.Modules.m1 = m1;
})(typeof window !== 'undefined' ? window : globalThis);
