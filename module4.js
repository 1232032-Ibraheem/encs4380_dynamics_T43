/* =====================================================================
 * module4.js — Module 4: "Guess the sensor"
 * A hidden sensor with random parameters is generated every round from a
 * seeded PRNG (seeded from crypto randomness), a noisy response is plotted,
 * the user estimates the parameters, and the page reveals and scores.
 * ===================================================================== */
(function (root) {
  'use strict';
  const M = root.DynMath, FO = M.firstOrder, SO = M.secondOrder;
  const { $, fmt, fixed, layout, colors, plot } = root.UI;

  const NOISE = 0.015;                 // noise σ as a fraction of the final value / amplitude
  const FULL = 0.10, ZERO = 0.50;      // full marks within 10 %, zero at 50 % error

  function newSeed() {
    try { const a = new Uint32Array(1); crypto.getRandomValues(a); return a[0]; }
    catch (e) { return (Date.now() ^ (Math.random() * 4294967296)) >>> 0; }
  }
  const sig = (v, n) => Number(v.toPrecision(n));

  /** points for one parameter: 100 within 10 %, linear to 0 at 50 % */
  function paramScore(guess, truth) {
    const err = Math.abs(guess - truth) / Math.abs(truth);
    const pts = err <= FULL ? 100 : Math.max(0, 100 * (1 - (err - FULL) / (ZERO - FULL)));
    return { err, pts };
  }

  /** Build a random hidden sensor and its noisy recorded response. */
  function makeSensor(seed) {
    const r = M.mulberry32(seed);
    const s = { seed };
    const pick = arr => arr[Math.floor(r() * arr.length)];
    if (r() < 0.5) {
      s.order = 1;
      s.K = sig(M.uniform(r, 0.5, 5), 3);
      s.tau = sig(M.logUniform(r, 0.2, 20), 3);
      s.A = pick([1, 2, 5, 10]);
      s.input = r() < 0.35 ? 'sine' : 'step';
      if (s.input === 'step') {
        const tEnd = 6 * s.tau * M.uniform(r, 0.9, 1.3), N = 200;
        s.t = M.linspace(0, tEnd, N);
        s.clean = s.t.map(t => FO.step(t, s.K, s.A, s.tau));
        s.y = s.clean.map(v => v + NOISE * s.K * s.A * M.gaussian(r));
        s.params = [{ key: 'K', label: 'Static sensitivity K', unit: '' }, { key: 'tau', label: 'Time constant τ', unit: 's' }];
      } else {
        s.w = sig(M.uniform(r, 0.5, 2.5) / s.tau, 3);       // ωτ between 0.5 and 2.5
        const T = 2 * Math.PI / s.w;
        const tEnd = Math.min(Math.max(5 * s.tau + 3 * T, 6 * T), 30 * T);
        const N = Math.min(1500, Math.max(300, Math.round(60 * tEnd / T)));
        s.t = M.linspace(0, tEnd, N);
        s.x = s.t.map(t => s.A * Math.sin(s.w * t));
        s.clean = s.t.map(t => FO.sine(t, s.K, s.A, s.w, s.tau, true));
        s.y = s.clean.map(v => v + NOISE * s.K * s.A * M.gaussian(r));
        s.params = [{ key: 'tau', label: 'Time constant τ', unit: 's' }];
      }
    } else {
      s.order = 2;
      s.zeta = sig(M.uniform(r, 0.1, 0.8), 2);
      s.wn = sig(M.logUniform(r, 1, 20), 3);
      s.K = sig(M.uniform(r, 0.5, 5), 3);
      s.A = 1; s.input = 'step';
      const tEnd = SO.timeWindow(s.zeta, s.wn) * M.uniform(r, 0.9, 1.2), N = 300;
      s.t = M.linspace(0, tEnd, N);
      s.clean = s.t.map(t => SO.step(t, s.K, s.A, s.zeta, s.wn));
      s.y = s.clean.map(v => v + NOISE * s.K * M.gaussian(r));
      s.params = [{ key: 'zeta', label: 'Damping ratio ζ', unit: '' }, { key: 'wn', label: 'Natural frequency ω<sub>n</sub>', unit: 'rad/s' }];
    }
    return s;
  }

  function curveFor(s, p) {
    // response of a sensor with parameters p (falls back to the true values)
    const q = Object.assign({ K: s.K, tau: s.tau, zeta: s.zeta, wn: s.wn }, p);
    if (s.order === 1 && s.input === 'step') return s.t.map(t => FO.step(t, q.K, s.A, q.tau));
    if (s.order === 1) return s.t.map(t => FO.sine(t, q.K, s.A, s.w, q.tau, true));
    return s.t.map(t => SO.step(t, q.K, s.A, q.zeta, q.wn));
  }

  const HINTS = {
    step1: [
      'Read the final value first: K = (final value) / A.',
      'The time at which the curve reaches 63.2 % of its final value is τ. Hover over the points to read it.',
      'Cross-check: at 3τ the curve is at 95 %. Or use the log method: ln[(y∞ − y)/y∞] against t has slope −1/τ.',
    ],
    sine1: [
      'Ignore the first few cycles: that is the start-up transient. Use the later cycles, where the output is a clean sine.',
      'Measure the time Δt from an input peak to the next output peak. Phase lag φ = ω·Δt in radians; then τ = tan(φ)/ω.',
      'Or use amplitude: M = (output amplitude)/(K·A), then τ = √(1/M² − 1)/ω.',
    ],
    step2: [
      'Overshoot depends only on ζ. Read y<sub>peak</sub> and the final value and compute M<sub>p</sub> = (y<sub>peak</sub> − y<sub>∞</sub>)/y<sub>∞</sub>.',
      'ζ = −ln M<sub>p</sub> / √(π² + ln² M<sub>p</sub>). For example 16 % overshoot gives ζ ≈ 0.5.',
      'The peak time t<sub>p</sub> = π/ω<sub>d</sub>, and ω<sub>n</sub> = ω<sub>d</sub>/√(1 − ζ²). Noise hides the exact peak, so average a few points around it.',
    ],
  };

  const m4 = {
    ready: false,
    s: null, submitted: false, hintsShown: 0,
    session: { rounds: 0, total: 0, full: 0, last: null },

    init() {
      try { const saved = sessionStorage.getItem('dct-m4-session'); if (saved) this.session = JSON.parse(saved); } catch (e) { /* storage off */ }
      $('m4-form').addEventListener('submit', e => { e.preventDefault(); this.submit(); });
      $('m4-new').addEventListener('click', () => this.newRound());
      $('m4-hint').addEventListener('click', () => this.hint());
      this.newRound();
      this.ready = true;
    },

    hintKey() { return this.s.order === 2 ? 'step2' : (this.s.input === 'sine' ? 'sine1' : 'step1'); },

    newRound() {
      this.s = makeSensor(newSeed());
      this.submitted = false; this.hintsShown = 0;
      const s = this.s;
      $('m4-hints').innerHTML = '';
      $('m4-hint').disabled = false; $('m4-hint').textContent = 'Show a hint';
      $('m4-reveal').hidden = true;
      $('m4-err').textContent = '';
      $('m4-submit').disabled = false;
      let brief;
      if (s.order === 1 && s.input === 'step') {
        brief = `<p><strong>First-order sensor, step input.</strong> A step of height <strong>A = ${s.A}</strong> was applied at t = 0. Estimate the sensor's K and τ.</p>`;
        $('m4-plot-title').textContent = 'Recorded step response';
      } else if (s.order === 1) {
        brief = `<p><strong>First-order sensor, sine input.</strong> Input x(t) = ${s.A}·sin(${fmt(s.w, 3)} t), switched on at t = 0 (ω = ${fmt(s.w, 3)} rad/s, f = ${fmt(s.w / (2 * Math.PI), 3)} Hz). The sensor's K = <strong>${s.K}</strong> is known. Estimate τ.</p>`;
        $('m4-plot-title').textContent = 'Recorded sine response';
      } else {
        brief = `<p><strong>Second-order sensor, step input.</strong> A unit step (A = 1) was applied at t = 0. Estimate the damping ratio ζ and the natural frequency ω<sub>n</sub>.</p>`;
        $('m4-plot-title').textContent = 'Recorded step response';
      }
      $('m4-brief').innerHTML = brief + `<p class="muted small">Measurement noise: ${NOISE * 100} % of full scale.</p>`;
      $('m4-fields').innerHTML = s.params.map(p => `
        <label class="field"><span>${p.label}${p.unit ? ` (${p.unit})` : ''}</span>
          <input type="number" step="any" min="0" inputmode="decimal" name="${p.key}" required></label>`).join('');
      this.renderPlot();
      this.renderScore();
    },

    renderPlot(guessCurve) {
      const c = colors(), s = this.s;
      const traces = [];
      if (s.x) traces.push({ x: s.t, y: s.x, mode: 'lines', name: 'Input x(t)', line: { color: c.in, dash: 'dash', width: 1.5 }, hovertemplate: 't = %{x:.4g} s<br>x = %{y:.4g}<extra></extra>' });
      traces.push({ x: s.t, y: s.y, mode: 'markers', name: 'Measured output (noisy)', marker: { color: c.ink2, size: 4, opacity: 0.75 },
        hovertemplate: 't = %{x:.4g} s<br>y = %{y:.4g}<extra></extra>' });
      if (this.submitted) {
        traces.push({ x: s.t, y: s.clean, mode: 'lines', name: 'True sensor', line: { color: c.out, width: 3 }, hovertemplate: 'true: %{y:.4g}<extra></extra>' });
        traces.push({ x: s.t, y: guessCurve, mode: 'lines', name: 'Your guess', line: { color: c.mark, width: 2.5, dash: 'dash' }, hovertemplate: 'guess: %{y:.4g}<extra></extra>' });
      }
      plot($('m4-plot'), traces, layout(c, { x: 'time t (s)', y: s.x ? 'signal' : 'output y(t)',
        xaxis: { range: [0, s.t[s.t.length - 1]], showspikes: true, spikemode: 'across', spikethickness: 1, spikecolor: c.ink3, spikedash: 'dot' },
        yaxis: { showspikes: true, spikemode: 'across', spikethickness: 1, spikecolor: c.ink3, spikedash: 'dot' } }));
    },

    hint() {
      const list = HINTS[this.hintKey()];
      if (this.hintsShown >= list.length) return;
      $('m4-hints').insertAdjacentHTML('beforeend', `<li>${list[this.hintsShown]}</li>`);
      this.hintsShown++;
      if (this.hintsShown >= list.length) { $('m4-hint').disabled = true; $('m4-hint').textContent = 'No more hints'; }
      else $('m4-hint').textContent = 'Show another hint';
    },

    submit() {
      if (this.submitted) return;
      const s = this.s, err = $('m4-err');
      const guess = {};
      for (const p of s.params) {
        const v = parseFloat($('m4-fields').querySelector(`[name="${p.key}"]`).value);
        if (!isFinite(v) || v <= 0) { err.innerHTML = `Enter a positive number for ${p.label}.`; return; }
        if (p.key === 'zeta' && v > 5) { err.textContent = 'ζ above 5 is outside the range this challenge uses; check your value.'; return; }
        guess[p.key] = v;
      }
      err.textContent = '';
      this.submitted = true;
      $('m4-submit').disabled = true;
      $('m4-fields').querySelectorAll('input').forEach(i => { i.readOnly = true; });

      const rows = s.params.map(p => Object.assign({ p, guess: guess[p.key], truth: s[p.key] }, paramScore(guess[p.key], s[p.key])));
      const score = Math.round(rows.reduce((a, r) => a + r.pts, 0) / rows.length);
      const allFull = rows.every(r => r.err <= FULL);
      const ses = this.session;
      ses.rounds++; ses.total += score; if (allFull) ses.full++; ses.last = score;
      try { sessionStorage.setItem('dct-m4-session', JSON.stringify(ses)); } catch (e) { /* ignore */ }

      // known values that were not estimated are shown as well
      const extra = [];
      if (s.order === 1 && s.input === 'sine') extra.push(['Static sensitivity K (given)', s.K], ['Amplitude ratio M at this ω', FO.magnitude(s.w, s.tau)], ['Phase lag at this ω (°)', -M.rad2deg(FO.phase(s.w, s.tau))]);
      if (s.order === 2) extra.push(['Static sensitivity K', s.K], ['Overshoot (%)', 100 * SO.formulas(s.zeta, s.wn).Mp], ['Peak time (s)', SO.formulas(s.zeta, s.wn).tp]);
      $('m4-table').innerHTML = `<thead><tr><th>Parameter</th><th class="num">Your estimate</th><th class="num">True value</th><th class="num">Error</th><th class="num">Points</th></tr></thead>
        <tbody>${rows.map(r => `<tr><td>${r.p.label}</td><td class="num">${fmt(r.guess, 4)}</td><td class="num">${fmt(r.truth, 4)}</td>
          <td class="num ${r.err <= FULL ? 'good' : 'bad'}">${fixed(100 * r.err, 1)} %</td><td class="num">${Math.round(r.pts)}</td></tr>`).join('')}
        ${extra.map(e => `<tr><td class="muted">${e[0]}</td><td></td><td class="num muted">${fmt(e[1], 4)}</td><td></td><td></td></tr>`).join('')}</tbody>`;
      $('m4-verdict').innerHTML = `<span class="badge ${allFull ? 'badge-pass' : score >= 50 ? 'badge-mark' : 'badge-fail'}">${score} / 100</span> ` +
        (allFull ? 'Every estimate is within 10 % — full marks.' : 'Estimates within 10 % score full marks; the score falls to zero at 50 % error. Compare the dashed guess with the true curve to see which feature you misread.') +
        ' Press <strong>New sensor</strong> to play again.';
      $('m4-reveal').hidden = false;
      const g = Object.assign({}, guess);
      this.renderPlot(curveFor(s, g));
      this.renderScore();
      root.UI.resizeAll($('m4'));
    },

    renderScore() {
      const s = this.session;
      root.UI.readouts($('m4-score'), [
        { k: 'Rounds played', v: String(s.rounds) },
        { k: 'Total score', v: `${s.total} pts` },
        { k: 'Average per round', v: s.rounds ? `${fixed(s.total / s.rounds, 1)}` : '—' },
      ]);
    },
    render() { if (this.s) this.renderPlot(this.submitted ? curveFor(this.s, this.lastGuess()) : undefined); },
    lastGuess() {
      const g = {};
      this.s.params.forEach(p => { g[p.key] = parseFloat($('m4-fields').querySelector(`[name="${p.key}"]`).value); });
      return g;
    },
  };

  // exposed for automated tests
  m4._makeSensor = makeSensor; m4._paramScore = paramScore;
  root.Modules = root.Modules || {};
  root.Modules.m4 = m4;
})(typeof window !== 'undefined' ? window : globalThis);
