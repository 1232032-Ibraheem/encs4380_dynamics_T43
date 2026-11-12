/* module2.js — Module 2: first-order frequency response
 * Input/output sines, amplitude ratio, phase lag, Bode plot with cutoff,
 * and the dynamic-error calculator. */
(function (root) {
  'use strict';
  const M = root.DynMath, FO = M.firstOrder;
  const { $, fmt, fixed, slider, segmented, readouts, layout, colors, plot } = root.UI;
  const TWO_PI = 2 * Math.PI;
  const W_MIN = 0.01, W_MAX = 20000;   // rad/s slider limits

  const m2 = {
    ready: false,
    init() {
      const box = $('m2-controls');
      box.insertAdjacentHTML('beforeend', '<h3>Parameters</h3>');
      this.tau = slider(box, { label: 'Time constant τ', unit: 's', min: 0.0005, max: 20, scale: 'log', value: 0.2, digits: 3 });
      this.unit = segmented(box, { label: 'Frequency unit', value: 'rad', options: [{ value: 'rad', label: 'ω (rad/s)' }, { value: 'hz', label: 'f (Hz)' }] });
      this.freq = slider(box, { label: 'Input frequency', unit: 'rad/s', min: W_MIN, max: W_MAX, scale: 'log', value: 2, digits: 4, help: 'ω = 2πf' });
      box.insertAdjacentHTML('beforeend', '<p class="muted small" style="margin:0">K = 1, A = 1, so the plot shows the amplitude ratio directly.</p>');
      const redraw = root.UI.scheduler(() => this.render());
      [this.tau, this.freq].forEach(s => s.onChange(redraw));
      // switching units converts the slider value, so the frequency itself does not change
      this.unit.onChange(u => {
        const v = this.freq.get();
        this.freq.el.querySelector('.unit').textContent = u === 'hz' ? 'Hz' : 'rad/s';
        if (u === 'hz') { this.freq.setBounds(W_MIN / TWO_PI, W_MAX / TWO_PI); this.freq.set(v / TWO_PI, true); }
        else { this.freq.setBounds(W_MIN, W_MAX); this.freq.set(v * TWO_PI, true); }
        redraw();
      });
      this.setupCalculator();
      this.setupPredict();
      this.ready = true;
    },

    /** current input frequency in rad/s */
    omega() { return this.unit.get() === 'hz' ? this.freq.get() * TWO_PI : this.freq.get(); },
    setOmega(w) { this.freq.set(this.unit.get() === 'hz' ? w / TWO_PI : w); },

    render() {
      const c = colors();
      const tau = this.tau.get(), w = this.omega();
      const Mg = FO.magnitude(w, tau), ph = FO.phase(w, tau), dt = -ph / w, T = TWO_PI / w;

      /* input and output, steady state, three periods */
      const t = M.linspace(0, 3 * T, 600);
      // time lag between an input peak and the matching output peak
      const tIn = 2.25 * T, tOut = tIn + dt;   // sin(ωt) peaks at t = T/4 + nT
      plot($('m2-time'), [
        { x: t, y: t.map(s => Math.sin(w * s)), mode: 'lines', name: 'Input x/A', line: { color: c.in, dash: 'dash' }, hovertemplate: 't = %{x:.4g} s<br>%{y:.3f}<extra></extra>' },
        { x: t, y: t.map(s => Mg * Math.sin(w * s + ph)), mode: 'lines', name: 'Output y/(KA)', line: { color: c.out, width: 3 }, hovertemplate: 't = %{x:.4g} s<br>%{y:.3f}<extra></extra>' },
      ], layout(c, { x: 'time t (s)', y: 'normalised amplitude', yaxis: { range: [-1.35, 1.35] }, xaxis: { range: [0, 3 * T] },
        extra: {
          shapes: [
            { type: 'line', xref: 'paper', x0: 0, x1: 1, y0: Mg, y1: Mg, line: { color: c.out, width: 1, dash: 'dot' } },
            { type: 'line', x0: tIn, x1: tIn, y0: 1, y1: 1.15, line: { color: c.mark, width: 1 } },
            { type: 'line', x0: tOut, x1: tOut, y0: Mg, y1: 1.15, line: { color: c.mark, width: 1 } },
          ],
          annotations: [
            { x: tIn, y: 1.15, ax: tOut, ay: 1.15, xref: 'x', yref: 'y', axref: 'x', ayref: 'y', showarrow: true, arrowhead: 2, arrowside: 'end+start', arrowcolor: c.mark, text: '' },
            { x: (tIn + tOut) / 2, y: 1.15, text: `Δt = ${fmt(dt, 3)} s`, showarrow: false, yshift: 11, font: { color: c.mark } },
            { xref: 'paper', x: 1, y: Mg, text: `M = ${fixed(Mg, 3)}`, showarrow: false, xanchor: 'right', yshift: -10, font: { color: c.out, size: 11 } },
          ],
        } }));

      /* Bode plot over five decades around the cutoff */
      const wc = FO.cutoff(tau);
      const ws = M.logspace(Math.log10(wc) - 2.5, Math.log10(wc) + 2.5, 400);
      const hz = ws.map(v => v / TWO_PI);
      const inRange = w >= ws[0] && w <= ws[ws.length - 1];
      const marker = { name: 'Current ω', mode: 'markers', marker: { color: c.mark, size: 10 }, visible: inRange };
      const cutLine = { type: 'line', x0: wc, x1: wc, yref: 'paper', y0: 0, y1: 1, line: { color: c.pole, dash: 'dash', width: 1 } };
      const xaxis = { type: 'log', exponentformat: 'power', range: [Math.log10(ws[0]), Math.log10(ws[ws.length - 1])] };
      const hover = 'ω = %{x:.4g} rad/s<br>f = %{customdata:.4g} Hz<br>';
      plot($('m2-bode-mag'), [
        { x: ws, y: ws.map(v => 20 * Math.log10(FO.magnitude(v, tau))), customdata: hz, mode: 'lines', name: '|G| (dB)', line: { color: c.out, width: 3 }, hovertemplate: hover + '%{y:.2f} dB<extra></extra>' },
        Object.assign({ x: [w], y: [20 * Math.log10(Mg)] }, marker),
      ], layout(c, { x: 'ω (rad/s)', y: 'magnitude (dB)', xaxis, yaxis: { range: [-52, 5] },
        extra: { shapes: [cutLine], annotations: [{ x: Math.log10(wc), y: -3, text: `ω_c = 1/τ = ${fmt(wc, 3)} rad/s<br>(−3 dB)`, showarrow: true, ax: 50, ay: 30, arrowcolor: c.pole, font: { size: 11, color: c.pole }, xanchor: 'left' }] } }));
      plot($('m2-bode-ph'), [
        { x: ws, y: ws.map(v => M.rad2deg(FO.phase(v, tau))), customdata: hz, mode: 'lines', name: 'φ (°)', line: { color: c.out, width: 3 }, hovertemplate: hover + '%{y:.2f}°<extra></extra>' },
        Object.assign({ x: [w], y: [M.rad2deg(ph)] }, marker),
      ], layout(c, { x: 'ω (rad/s)', y: 'phase (°)', xaxis, yaxis: { range: [-95, 5], dtick: 15 },
        extra: { shapes: [cutLine], annotations: [{ x: Math.log10(wc), y: -45, text: '−45° at ω_c', showarrow: true, ax: 45, ay: -25, arrowcolor: c.pole, font: { size: 11, color: c.pole }, xanchor: 'left' }] } }));

      readouts($('m2-readouts'), [
        { k: 'Input frequency', v: `${fmt(w, 4)} rad/s`, h: `f = ${fmt(w / TWO_PI, 4)} Hz` },
        { k: 'ωτ', v: fmt(w * tau, 4) },
        { k: 'Amplitude ratio M', v: fixed(Mg, 4), emph: true, h: `${fixed(20 * Math.log10(Mg), 2)} dB` },
        { k: 'Phase lag', v: `${fixed(-M.rad2deg(ph), 2)}°`, emph: true, h: `${fixed(-ph, 4)} rad, Δt = ${fmt(dt, 3)} s` },
        { k: 'Dynamic error δ = M − 1', v: `${fixed(100 * (Mg - 1), 2)} %` },
        { k: 'Cutoff ω_c = 1/τ', v: `${fmt(wc, 4)} rad/s`, h: `f_c = ${fmt(wc / TWO_PI, 4)} Hz` },
      ]);
      $('m2-interpret').innerHTML = `
        <li>At ωτ = ${fmt(w * tau, 3)} the output amplitude is ${fixed(100 * Mg, 1)} % of the input and it lags by ${fixed(-M.rad2deg(ph), 1)}°.</li>
        <li>Below the cutoff the sensor is accurate; above it the amplitude falls by 20 dB per decade and the lag approaches 90°.</li>`;
      this.calcA();
    },

    /* dynamic-error calculator: A = error at a frequency, B = largest τ for a design */
    setupCalculator() {
      ['c-tau', 'c-f'].forEach(id => $(id).addEventListener('input', () => this.calcA()));
      ['c-dmax', 'c-fmax', 'c-fph'].forEach(id => $(id).addEventListener('input', () => this.calcB()));
      $('c-load3').addEventListener('click', () => {
        $('c-dmax').value = 5; $('c-fmax').value = 100; $('c-fph').value = 50;
        this.calcB();
      });
      this.calcB();
    },
    calcA() {
      const tau = parseFloat($('c-tau').value), f = parseFloat($('c-f').value);
      if (!(tau > 0) || !(f > 0)) { $('c-out-a').innerHTML = '<p class="small muted">Enter a positive τ and f.</p>'; return; }
      const w = TWO_PI * f, Mg = FO.magnitude(w, tau);
      readouts($('c-out-a'), [
        { k: 'M', v: fixed(Mg, 4), h: `ω = ${fmt(w, 4)} rad/s` },
        { k: 'δ = M − 1', v: `${fixed(100 * (Mg - 1), 3)} %`, emph: true },
        { k: 'Phase lag', v: `${fixed(-M.rad2deg(FO.phase(w, tau)), 2)}°` },
      ]);
    },
    calcB() {
      const d = parseFloat($('c-dmax').value) / 100, fmax = parseFloat($('c-fmax').value), fph = parseFloat($('c-fph').value);
      if (!(d > 0 && d < 1) || !(fmax > 0) || !(fph > 0)) { $('c-out-b').innerHTML = '<p class="small muted">Use 0 &lt; |δ| &lt; 100 % and positive frequencies.</p>'; return; }
      const tm = FO.tauMaxForDynamicError(d, fmax);
      readouts($('c-out-b'), [
        { k: 'τ_max', v: tm < 0.1 ? `${fmt(tm * 1000, 4)} ms` : `${fmt(tm, 4)} s`, emph: true },
        { k: `Phase lag at ${fmt(fph, 4)} Hz`, v: `${fixed(-M.rad2deg(FO.phase(TWO_PI * fph, tm)), 2)}°`, emph: true },
      ]);
    },

    setupPredict() {
      const self = this;
      root.UI.predict($('m2-predict'), [
        () => {
          const tau = self.tau.get();
          return { kind: 'choice', text: `τ = ${fmt(tau, 3)} s. At the cutoff ω = 1/τ = ${fmt(1 / tau, 3)} rad/s, what is the amplitude ratio M?`,
            options: ['0.5', '0.707', '0.9', '1.0'], answer: 1,
            explain: 'At ωτ = 1, M = 1/√2 = 0.707 (−3 dB) and the phase lag is 45°.',
            tryLabel: 'Set ω = 1/τ', apply: () => self.setOmega(1 / tau) };
        },
        () => {
          const tau = self.tau.get(), w = self.omega();
          return { kind: 'number', text: `With τ = ${fmt(tau, 3)} s and ω = ${fmt(w, 3)} rad/s, predict the amplitude ratio M.`,
            answer: FO.magnitude(w, tau), tol: 0.03, unit: '',
            explain: `M = 1/√(1 + (ωτ)²) with ωτ = ${fmt(w * tau, 3)}.` };
        },
        () => {
          const tau = self.tau.get(), th = Math.max(0.0005, tau / 2);
          return { kind: 'choice', text: `If τ is halved to ${fmt(th, 3)} s, the cutoff frequency will…`,
            options: ['halve', 'stay the same', 'double'], answer: 2,
            explain: 'ω_c = 1/τ: a faster sensor follows signals of higher frequency.',
            tryLabel: `Set τ = ${fmt(th, 3)} s`, apply: () => self.tau.set(th) };
        },
      ]);
    },
  };


  /* validation checks 2 and 3 (Module 2) */
  root.DynValidation.add({ id: 2, module: 'Module 2', title: 'First-order sine response',
    description: 'τ = 0.2 s with input sin 2t: output amplitude 0.93 of the input, lagging by 21.8°.',
    inputs: 'τ = 0.2 s, ω = 2 rad/s (f = 0.318 Hz), K = 1',
    run() {
      const V = root.DynValidation, tau = 0.2, w = 2;
      // simulate 5 periods from rest and measure amplitude and phase from the waveforms
      const t = M.linspace(0, 10 * Math.PI, 8001);
      const meas = M.measureSinePair(t, t.map(s => Math.sin(w * s)), t.map(s => FO.sine(s, 1, 1, w, tau, true)), w, 3);
      return [
        V.check('Amplitude ratio M (formula)', 0.93, FO.magnitude(w, tau), 0.005),
        V.check('Amplitude ratio M (measured from the waveforms)', 0.93, meas.ratio, 0.005),
        V.check('Phase lag (formula)', 21.8, -M.rad2deg(FO.phase(w, tau)), 0.05, '°'),
        V.check('Phase lag (measured from the waveforms)', 21.8, -M.rad2deg(meas.phase), 0.05, '°'),
      ];
    } });
  root.DynValidation.add({ id: 3, module: 'Module 2', title: 'Dynamic-error design',
    description: 'Signals up to 100 Hz with at most 5 % dynamic error need τ_max ≈ 0.52 ms; the phase lag at 50 Hz is then 9.3°.',
    inputs: '|δ| ≤ 0.05, f_max = 100 Hz, phase at 50 Hz with τ = τ_max',
    run() {
      const V = root.DynValidation, tauMax = FO.tauMaxForDynamicError(0.05, 100);
      return [
        V.check('τ_max', 0.52, tauMax * 1000, 0.01, 'ms'),
        V.check('Dynamic error δ at 100 Hz', -5.0, FO.dynamicError(TWO_PI * 100, tauMax) * 100, 0.01, '%'),
        V.check('Phase lag at 50 Hz', 9.3, -M.rad2deg(FO.phase(TWO_PI * 50, tauMax)), 0.1, '°'),
      ];
    } });

  root.Modules = root.Modules || {};
  root.Modules.m2 = m2;
})(typeof window !== 'undefined' ? window : globalThis);
