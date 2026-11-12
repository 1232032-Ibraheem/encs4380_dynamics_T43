/* module3.js — Module 3: second-order response and state-space panel.
 * (ζ, ωn) and (m, b, k) describe the same mass-spring-damper and stay in sync. */
(function (root) {
  'use strict';
  const M = root.DynMath, SO = M.secondOrder, SS = M.stateSpace;
  const { $, fmt, fixed, slider, readouts, layout, colors, plot, tex } = root.UI;
  const BAND = 0.02;   // ±2 % settling band
  const F0 = 1;        // step force for the state-space simulation [N]

  /** exact time at which the closed-form response first reaches `level` (bisection) */
  function exactCrossing(fn, level, tEnd) {
    let prev = 0;
    for (let i = 1; i <= 4000; i++) {
      const t = tEnd * i / 4000;
      if (fn(t) >= level) {
        let lo = prev, hi = t;
        for (let k = 0; k < 60; k++) { const mid = (lo + hi) / 2; if (fn(mid) >= level) hi = mid; else lo = mid; }
        return (lo + hi) / 2;
      }
      prev = t;
    }
    return null;
  }

  const m3 = {
    ready: false,
    p: { m: 1, b: 2, k: 4 },   // the plant; ζ and ωn are computed from it

    init() {
      const box = $('m3-controls');
      box.insertAdjacentHTML('beforeend', '<h3>ζ and ω<sub>n</sub></h3>');
      this.zeta = slider(box, { label: 'Damping ratio ζ', min: 0, max: 3, step: 0.01, value: 0.5 });
      this.wn = slider(box, { label: 'Natural frequency ω<sub>n</sub>', unit: 'rad/s', min: 0.2, max: 20, scale: 'log', value: 2, digits: 4 });
      const presets = root.UI.el('div', { class: 'btn-row' });
      presets.innerHTML = [0.2, 1, 2].map(z => `<button type="button" class="btn small ghost" data-z="${z}">ζ = ${z}</button>`).join('');
      box.appendChild(presets);
      presets.querySelectorAll('button').forEach(b => b.addEventListener('click', () => this.zeta.set(+b.dataset.z)));

      box.insertAdjacentHTML('beforeend', '<h3 style="margin-top:6px">Mass, damper, spring</h3>');
      this.m = slider(box, { label: 'Mass m', unit: 'kg', min: 0.1, max: 5, step: 0.1, value: 1 });
      this.b = slider(box, { label: 'Damping b', unit: 'N·s/m', min: 0.001, max: 600, scale: 'log', value: 2, digits: 4, help: 'For b = 0 set ζ = 0.' });
      this.k = slider(box, { label: 'Stiffness k', unit: 'N/m', min: 0.004, max: 2000, scale: 'log', value: 4, digits: 4 });

      const redraw = root.UI.scheduler(() => this.render());
      // ζ or ωn moved: keep m, compute b and k
      const fromZetaWn = () => {
        const q = SS.secondOrderToMsd(this.p.m, this.zeta.get(), this.wn.get());
        this.p = { m: q.m, b: q.b, k: q.k };
        this.b.set(q.b, true); this.k.set(q.k, true);
        redraw();
      };
      // m, b or k moved: compute ζ and ωn
      const fromMBK = () => {
        this.p = { m: this.m.get(), b: this.b.get(), k: this.k.get() };
        const q = SS.msdToSecondOrder(this.p.m, this.p.b, this.p.k);
        this.zeta.set(q.zeta, true); this.wn.set(q.wn, true);
        if (Math.abs(this.zeta.get() - q.zeta) > 1e-9 || Math.abs(this.wn.get() - q.wn) > 1e-9 * q.wn) {
          // outside the ζ/ωn slider range: use the limited values for both views
          this.p = SS.secondOrderToMsd(this.p.m, this.zeta.get(), this.wn.get());
          this.b.set(this.p.b, true); this.k.set(this.p.k, true);
        }
        redraw();
      };
      this.zeta.onChange(fromZetaWn); this.wn.onChange(fromZetaWn);
      this.m.onChange(fromMBK); this.b.onChange(fromMBK); this.k.onChange(fromMBK);

      tex($('ss-eq'), String.raw`\dot{\mathbf x}=\begin{bmatrix}0&1\\-\tfrac{k}{m}&-\tfrac{b}{m}\end{bmatrix}\mathbf x+\begin{bmatrix}0\\ \tfrac1m\end{bmatrix}F,\quad y=\begin{bmatrix}1&0\end{bmatrix}\mathbf x`, true);
      this.setupPredict();
      this.ready = true;
    },

    setMBK(m, b, k) { this.m.set(m, true); this.b.set(b, true); this.k.set(k); },

    render() {
      const c = colors();
      const { m, b, k } = this.p;
      const { zeta: z, wn } = SS.msdToSecondOrder(m, b, k);
      const cls = SO.classify(z);
      $('m3-class').textContent = cls;

      /* step response (K = 1, unit step) and metrics read off the sampled curve */
      const tEnd = SO.timeWindow(z, wn);
      const t = M.linspace(0, tEnd, 3001);
      const fn = s => SO.step(s, 1, 1, z, wn);
      const y = t.map(fn);
      const meas = M.measureStep(t, y, 1, BAND);
      const form = SO.formulas(z, wn);

      const pts = { x: [], y: [], text: [], pos: [] };
      const add = (x, yv, label, pos) => { pts.x.push(x); pts.y.push(yv); pts.text.push(label); pts.pos.push(pos); };
      if (meas.t10 !== null) { add(meas.t10, 0.1, '10 %', 'middle right'); add(meas.t90, 0.9, '90 %', 'middle left'); }
      if (meas.tp !== null) add(meas.tp, meas.yPeak, `peak ${fixed(100 * meas.overshoot, 1)} %`, 'top right');
      if (meas.ts > 0) add(meas.ts, fn(meas.ts), 't_s', 'bottom right');
      plot($('m3-step'), [
        { x: t, y, mode: 'lines', name: `y(t), ζ = ${fixed(z, 2)}`, line: { color: c.out, width: 3 }, hovertemplate: 't = %{x:.4g} s<br>y = %{y:.4g}<extra></extra>' },
        { x: pts.x, y: pts.y, text: pts.text, textposition: pts.pos, mode: 'markers+text', name: 'Read off the curve',
          marker: { color: c.mark, size: 9 }, textfont: { color: c.mark, size: 11 }, hovertemplate: 't = %{x:.4g} s<br>y = %{y:.4g}<extra></extra>' },
      ], layout(c, { x: 'time t (s)', y: 'output y(t)', xaxis: { range: [0, tEnd] },
        yaxis: { range: [-0.1, z === 0 ? 2.1 : Math.max(1.15, Math.max.apply(null, y) * 1.12)] },
        extra: { shapes: [
          { type: 'rect', xref: 'paper', x0: 0, x1: 1, y0: 1 - BAND, y1: 1 + BAND, fillcolor: c.markSoft, line: { width: 0 }, layer: 'below' },
          { type: 'line', xref: 'paper', x0: 0, x1: 1, y0: 1, y1: 1, line: { color: c.in, width: 1, dash: 'dash' } },
        ] } }));

      // metrics: measured on the curve vs formula
      const t10 = exactCrossing(fn, 0.1, tEnd), t90 = exactCrossing(fn, 0.9, tEnd);
      const sec = v => (v === null || v === undefined) ? '—' : `${fmt(v, 4)} s`;
      const under = z > 0 && z < 1;
      const rows = [
        ['Rise time (10 % → 90 %)', sec(meas.tr1090), sec(t10 !== null && t90 !== null ? t90 - t10 : null)],
        ['Rise time (0 → 100 %)', under ? sec(meas.tr0100) : '—', under ? sec(form.tr0100) : '—'],
        ['Peak time t_p', meas.tp !== null ? sec(meas.tp) : 'no peak', z < 1 ? sec(form.tp) : '—'],
        ['Overshoot M_p', `${fixed(100 * meas.overshoot, 2)} %`, `${fixed(100 * form.Mp, 2)} %`],
        ['Settling time t_s (±2 %)', z === 0 ? 'never' : (meas.settled ? sec(meas.ts) : 'beyond plot'), under ? `≈ ${sec(form.ts2)}` : '—'],
      ];
      $('m3-metrics').innerHTML = `<thead><tr><th>Quantity</th><th class="num">Read off the curve</th><th class="num">Formula</th></tr></thead>
        <tbody>${rows.map(r => `<tr><td>${r[0]}</td><td class="num"><strong>${r[1]}</strong></td><td class="num">${r[2]}</td></tr>`).join('')}</tbody>`;

      this.renderStateSpace(c, z, wn);
      $('m3-interpret').innerHTML = `
        <li>ζ = ${fixed(z, 2)}: the system is <strong>${cls}</strong>. ${z < 1 ? 'Complex poles give an oscillating response.' : 'Real poles give a response without oscillation.'}</li>
        <li>Overshoot depends only on ζ (formula e<sup>−πζ/√(1−ζ²)</sup>); changing ω<sub>n</sub> only scales the time axis.</li>`;
    },

    renderStateSpace(c, z, wn) {
      const { m, b, k } = this.p;
      const { A, B } = SS.msd(m, b, k);
      const num = v => (Math.abs(v) < 1e-12 ? '0' : String(Number(v.toPrecision(4))));
      tex($('ss-mat'), String.raw`A=\begin{bmatrix}0&1\\${num(A[1][0])}&${num(A[1][1])}\end{bmatrix},\quad B=\begin{bmatrix}0\\${num(B[1][0])}\end{bmatrix}`, true);
      const eig = M.eig2x2(A), poles = SO.poles(z, wn);
      const diff = Math.max(M.complexDistance(eig[0], poles[0]), M.complexDistance(eig[1], poles[1]));
      readouts($('ss-readouts'), [
        { k: 'Eigenvalues of A', v: `${M.formatComplex(eig[0], 3)}<br>${M.formatComplex(eig[1], 3)}`, emph: true },
        { k: 'Poles of G(s)', v: `${M.formatComplex(poles[0], 3)}<br>${M.formatComplex(poles[1], 3)}`, emph: true },
        { k: 'max |λ − p|', v: fmt(diff, 2), h: diff < 1e-9 ? 'equal' : 'different' },
        { k: 'det(sI − A)', v: `s² + ${fmt(b / m, 4)}s + ${fmt(k / m, 4)}` },
      ]);

      /* s-plane: poles (×) and eigenvalues (○) */
      const R = Math.max(Math.hypot(poles[0].re, poles[0].im), Math.hypot(poles[1].re, poles[1].im), wn) * 1.25;
      plot($('m3-splane'), [
        { x: poles.map(p => p.re), y: poles.map(p => p.im), mode: 'markers', name: 'Poles', marker: { symbol: 'x-thin', size: 15, color: c.pole, line: { width: 3, color: c.pole } }, hovertemplate: '%{x:.4f} %{y:+.4f}j<extra>pole</extra>' },
        { x: eig.map(p => p.re), y: eig.map(p => p.im), mode: 'markers', name: 'Eigenvalues of A', marker: { symbol: 'circle-open', size: 21, color: c.out, line: { width: 2.5 } }, hovertemplate: '%{x:.4f} %{y:+.4f}j<extra>eigenvalue</extra>' },
      ], layout(c, { x: 'Re(s)', y: 'Im(s)',
        xaxis: { range: [Math.min(-R, poles[1].re * 1.1), R * 0.35] },
        yaxis: { range: [-R, R], scaleanchor: 'x', scaleratio: 1 } }));

      /* phase plane: RK4 simulation of ẋ = Ax + B·F0 from rest */
      const tEnd = SO.timeWindow(z, wn);
      const h = Math.min(tEnd / 2000, 0.02 * 2 * Math.PI / wn);
      const sim = M.integrate(M.linearSystem(A, B, () => F0), [0, 0], 0, tEnd, h, 'rk4');
      const step = Math.max(1, Math.floor(sim.t.length / 2500));
      const X = sim.X.filter((_, i) => i % step === 0);
      plot($('m3-phase'), [
        { x: X.map(v => v[0]), y: X.map(v => v[1]), mode: 'lines', name: 'Trajectory', line: { color: c.out, width: 2.5 }, hovertemplate: 'x₁ = %{x:.4g} m<br>x₂ = %{y:.4g} m/s<extra></extra>' },
        { x: [0], y: [0], mode: 'markers', name: 'Start', marker: { color: c.in, size: 9 }, hoverinfo: 'skip' },
        { x: [F0 / k], y: [0], mode: 'markers', name: 'Equilibrium F₀/k', marker: { color: c.pole, size: 12, symbol: 'x-thin', line: { width: 3, color: c.pole } }, hoverinfo: 'skip' },
      ], layout(c, { x: 'x₁ displacement (m)', y: 'x₂ velocity (m/s)' }));
    },

    setupPredict() {
      const self = this;
      root.UI.predict($('m3-predict'), [
        () => ({ kind: 'choice', text: 'If ζ is increased from 0.2 to 0.7 (same ω<sub>n</sub>), the overshoot will…',
          options: ['increase', 'decrease', 'stay the same'], answer: 1,
          explain: 'M<sub>p</sub> = e<sup>−πζ/√(1−ζ²)</sup>: 52.7 % at ζ = 0.2, 4.6 % at ζ = 0.7.',
          tryLabel: 'Set ζ = 0.7', apply: () => self.zeta.set(0.7) }),
        () => ({ kind: 'choice', text: 'For m = 1, b = 2, k = 4, what are the eigenvalues of A?',
          options: ['−1 ± j1.732', '−2 ± j2', '−1 and −4', '± j2'], answer: 0,
          explain: 'det(sI − A) = s² + 2s + 4 = 0 gives s = −1 ± j√3. These are also the poles (ζ = 0.5, ω<sub>n</sub> = 2).',
          tryLabel: 'Load m = 1, b = 2, k = 4', apply: () => self.setMBK(1, 2, 4) }),
        () => ({ kind: 'choice', text: 'With ζ = 2, will the phase-plane trajectory spiral around the equilibrium?',
          options: ['yes', 'no'], answer: 1,
          explain: 'With ζ > 1 the eigenvalues are real, so there is no oscillation and no spiral.',
          tryLabel: 'Set ζ = 2', apply: () => self.zeta.set(2) }),
      ]);
    },
  };


  /* validation checks 4, 5 and 6 (Module 3) */
  root.DynValidation.add({ id: 4, module: 'Module 3', title: 'Second-order overshoot, ζ = 0.5',
    description: 'With ζ = 0.5 the step response overshoots by 16.3 %.',
    inputs: 'ζ = 0.5, ωn = 2 rad/s, K = 1, unit step',
    run() {
      const V = root.DynValidation, c = V.sampled(t => SO.step(t, 1, 1, 0.5, 2), SO.timeWindow(0.5, 2), 4001);
      return [
        V.check('Overshoot (read off the curve)', 16.3, M.measureStep(c.t, c.y, 1, 0.02).overshoot * 100, 0.05, '%'),
        V.check('Overshoot (formula e^(−πζ/√(1−ζ²)))', 16.3, SO.formulas(0.5, 2).Mp * 100, 0.05, '%'),
      ];
    } });
  root.DynValidation.add({ id: 5, module: 'Module 3', title: 'Critically damped, ζ = 1',
    description: 'With ζ = 1 the step response never exceeds its final value.',
    inputs: 'ζ = 1, ωn = 2 rad/s, K = 1, unit step',
    run() {
      const V = root.DynValidation, c = V.sampled(t => SO.step(t, 1, 1, 1, 2), 3 * SO.timeWindow(1, 2), 4001);
      return [
        V.check('Overshoot (read off the curve)', 0, M.measureStep(c.t, c.y, 1, 0.02).overshoot * 100, 1e-6, '%'),
        V.check('max y(t) − final value (must not be positive)', 0, Math.max(0, Math.max.apply(null, c.y) - 1), 1e-9),
      ];
    } });
  root.DynValidation.add({ id: 6, module: 'Module 3', title: 'State-space eigenvalues = poles',
    description: 'm = 1, b = 2, k = 4 gives ζ = 0.5, ωn = 2 rad/s and eigenvalues −1 ± j1.732, equal to the poles.',
    inputs: 'm = 1 kg, b = 2 N·s/m, k = 4 N/m, A = [[0, 1], [−4, −2]]',
    run() {
      const V = root.DynValidation, eq = SS.msdToSecondOrder(1, 2, 4);
      const eig = M.eig2x2(SS.msd(1, 2, 4).A), poles = SO.poles(eq.zeta, eq.wn);
      const diff = Math.max(M.complexDistance(eig[0], poles[0]), M.complexDistance(eig[1], poles[1]));
      return [
        V.check('ζ = b / (2√(km))', 0.5, eq.zeta, 1e-9),
        V.check('ωn = √(k/m)', 2, eq.wn, 1e-9, 'rad/s'),
        V.check('Re λ (eigenvalues of A)', -1, eig[0].re, 1e-9),
        V.check('|Im λ| (eigenvalues of A)', 1.732, Math.abs(eig[0].im), 0.0005),
        V.check('max |λ − p| (eigenvalues vs poles)', 0, diff, 1e-9, '',
          `λ = ${M.formatComplex(eig[0], 4)}, ${M.formatComplex(eig[1], 4)}; p = ${M.formatComplex(poles[0], 4)}, ${M.formatComplex(poles[1], 4)}`),
      ];
    } });

  root.Modules = root.Modules || {};
  root.Modules.m3 = m3;
})(typeof window !== 'undefined' ? window : globalThis);

/* ---------------- Bonus: Euler and RK4 against the exact solution ---------------- */
(function (root) {
  'use strict';
  const M = root.DynMath, FO = M.firstOrder, SO = M.secondOrder;
  const { $, fmt, slider, segmented, readouts, layout, colors, plot } = root.UI;

  const bonus = {
    ready: false,
    init() {
      const box = $('b-controls');
      box.insertAdjacentHTML('beforeend', '<h3>Parameters</h3>');
      this.sys = segmented(box, { label: 'System', value: '1', options: [{ value: '1', label: 'First order' }, { value: '2', label: 'Second order' }] });
      this.tau = slider(box, { label: 'Time constant τ', unit: 's', min: 0.1, max: 10, scale: 'log', value: 1, digits: 3 });
      this.zeta = slider(box, { label: 'Damping ratio ζ', unit: '', min: 0.05, max: 2, step: 0.01, value: 0.3 });
      this.wn = slider(box, { label: 'Natural frequency ω<sub>n</sub>', unit: 'rad/s', min: 0.5, max: 10, scale: 'log', value: 2, digits: 3 });
      this.h = slider(box, { label: 'Solver step h', unit: 's', min: 0.005, max: 3, scale: 'log', value: 0.25, digits: 3,
        help: 'Euler is stable only for h < 2τ (first order) or h < 2ζ/ω<sub>n</sub> (second order).' });
      const r = root.UI.scheduler(() => this.render());
      [this.sys, this.tau, this.zeta, this.wn, this.h].forEach(s => s.onChange(r));
      this.ready = true;
    },
    render() {
      const c = colors();
      const first = this.sys.get() === '1';
      this.tau.show(first); this.zeta.show(!first); this.wn.show(!first);
      const h = this.h.get();
      let f, exact, tEnd, x0, limit;
      if (first) {
        const tau = this.tau.get();
        f = (t, x) => [(1 - x[0]) / tau];
        exact = t => FO.step(t, 1, 1, tau);
        tEnd = 6 * tau; x0 = [0]; limit = 2 * tau;
      } else {
        const z = this.zeta.get(), wn = this.wn.get();
        f = (t, x) => [x[1], wn * wn * (1 - x[0]) - 2 * z * wn * x[1]];
        exact = t => SO.step(t, 1, 1, z, wn);
        tEnd = SO.timeWindow(z, wn); x0 = [0, 0];
        // Euler is stable when |1 + hλ| < 1 for both poles λ
        limit = z < 1 ? 2 * z / wn : 2 / (wn * (z + Math.sqrt(z * z - 1)));
      }
      const hh = Math.min(h, tEnd / 4);
      const eu = M.integrate(f, x0, 0, tEnd, hh, 'euler');
      const rk = M.integrate(f, x0, 0, tEnd, hh, 'rk4');
      const te = M.linspace(0, tEnd, 800);
      const ye = te.map(exact);
      const eE = eu.t.map((t, i) => Math.abs(eu.X[i][0] - exact(t)));
      const eR = rk.t.map((t, i) => Math.abs(rk.X[i][0] - exact(t)));
      const maxE = Math.max.apply(null, eE), maxR = Math.max.apply(null, eR);
      plot($('b-resp'), [
        { x: te, y: ye, mode: 'lines', name: 'Closed form', line: { color: c.in, width: 2 }, hovertemplate: 'exact %{y:.5g}<extra></extra>' },
        { x: eu.t, y: eu.X.map(v => v[0]), mode: 'lines+markers', name: 'Euler', line: { color: c.third, width: 1.3 }, marker: { size: 5 }, hovertemplate: 'Euler %{y:.5g}<extra></extra>' },
        { x: rk.t, y: rk.X.map(v => v[0]), mode: 'markers', name: 'RK4', marker: { color: c.out, size: 8, symbol: 'diamond' }, hovertemplate: 'RK4 %{y:.5g}<extra></extra>' },
      ], layout(c, { x: 'time t (s)', y: 'normalised output y/(KA)', xaxis: { range: [0, tEnd] }, yaxis: { range: [-0.5, 2.2] } }));
      const pos = (t, e) => { const o = { x: [], y: [] }; t.forEach((v, i) => { if (e[i] > 0) { o.x.push(v); o.y.push(e[i]); } }); return o; };
      const pE = pos(eu.t, eE), pR = pos(rk.t, eR);
      plot($('b-err'), [
        { x: pE.x, y: pE.y, mode: 'lines+markers', name: 'Euler |error|', line: { color: c.third }, marker: { size: 4 } },
        { x: pR.x, y: pR.y, mode: 'lines+markers', name: 'RK4 |error|', line: { color: c.out }, marker: { size: 4 } },
      ], layout(c, { x: 'time t (s)', y: '|numerical − exact|', yaxis: { type: 'log', exponentformat: 'power' }, xaxis: { range: [0, tEnd] } }));
      readouts($('b-readouts'), [
        { k: 'Step h used', v: `${fmt(hh, 3)} s`, h: `${eu.t.length - 1} steps` },
        { k: 'Max |error|, Euler', v: fmt(maxE, 3) },
        { k: 'Max |error|, RK4', v: fmt(maxR, 3), emph: true },
        { k: 'RK4 is more accurate by', v: maxR > 0 ? `${fmt(maxE / maxR, 3)}×` : '—' },
        { k: 'Euler stability limit', v: `h < ${fmt(limit, 3)} s`, h: limit && hh > limit ? 'exceeded: Euler diverges' : '' },
      ]);
    },
  };
  root.Modules = root.Modules || {};
  root.Modules.bonus = bonus;
})(typeof window !== 'undefined' ? window : globalThis);
