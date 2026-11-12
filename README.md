# Dynamic Characteristics Interactive Toolkit

ENCS4380 Interfacing Techniques — Homework 1, Birzeit University
Instructor: Dr. Wasel Ghanem

**Submitted by:** Ibraheem Mhanna (1232032) — individual submission, all four modules.

**Live site:** https://1232032-ibraheem.github.io/encs4380_dynamics_T43/

An interactive website that shows the dynamic characteristics of first- and second-order sensors. All graphs are computed in the browser and update as the sliders move.

## Modules

1. **First-order time response** (`module1.js`) — step and ramp inputs with K and τ sliders; markers at 63.2 % (t = τ) and 95 % (t = 3τ); ramp lag and steady-state error KRτ; error fraction on a log axis with a fitted slope of −1/τ.
2. **First-order frequency response** (`module2.js`) — input and output sines on one plot; amplitude ratio and phase lag; Bode magnitude and phase plots with the cutoff ω = 1/τ; frequency in rad/s or Hz; dynamic-error calculator (δ = M − 1, and the largest τ for a given error and bandwidth).
3. **Second-order response** (`module3.js`) — ζ and ωn sliders, or m, b and k (kept in sync); under-, critically and overdamped step responses; rise time, peak time, overshoot and settling time read off the curve and compared with the formulas; poles in the s-plane. The state-space panel shows A and B updating live, the eigenvalues of A on top of the poles, and the (x₁, x₂) phase-plane trajectory.
4. **Guess the sensor** (`module4.js`) — a random first- or second-order sensor with hidden parameters and measurement noise. Estimate τ (and K) or ζ and ωn from the graph; the page reveals the true values, overlays the true and guessed curves and scores the guess (full marks within 10 %), with hints, a **New sensor** button and a session score.

Each module has a **Predict, then check** question. `module1.js` also contains the calculations, interface helpers and validation code shared by all modules. The **Bonus** page (in `module3.js`) compares Euler and RK4 numerical solutions with the exact response.

## How to use

Open the live site, choose a module from the navigation bar, and move the sliders or type values in the boxes. Hover over a graph to read values; drag to zoom and double-click to reset. It works on a phone.

To run it locally, serve the folder with any static web server (for example `python -m http.server`) and open `index.html`. Plotly.js 4.3.0 and KaTeX 0.19.1 are loaded from the jsDelivr CDN, so an internet connection is needed.

## Validation

The **Validation** page runs the six checks from the homework sheet with the same functions that draw the plots, and shows the expected value, the calculated value and PASS/FAIL:

| Check | Expected | Calculated |
|---|---|---|
| 1. First-order step | y(τ) = 0.632 KA, y(3τ) = 0.95 KA | 0.6321 KA, 0.9502 KA |
| 2. τ = 0.2 s, input sin 2t | M = 0.93, lag 21.8° | 0.9285, 21.80° |
| 3. 100 Hz, ≤ 5 % dynamic error | τmax ≈ 0.52 ms, 9.3° at 50 Hz | 0.5231 ms, 9.33° |
| 4. ζ = 0.5 | overshoot 16.3 % | 16.30 % |
| 5. ζ = 1 | no overshoot | 0 % |
| 6. m = 1, b = 2, k = 4 | eigenvalues −1 ± j1.732 = poles | −1 ± j1.7321, difference 0 |

All six pass. Before submission the calculations were also checked with an automated test script (58 checks), and the website was tested in Chrome at desktop and phone sizes, using every module, slider and button (115 checks, no errors).

## Files

| File | Contents |
|---|---|
| `index.html` | page layout and navigation |
| `style.css` | styling |
| `module1.js` | shared maths, interface helpers and validation registry; Module 1 |
| `module2.js` | Module 2 |
| `module3.js` | Module 3, state-space panel and bonus |
| `module4.js` | Module 4 |

## AI usage disclosure

This is an individual submission by Ibraheem Mhanna (1232032). I reviewed the homework requirements, set the scope and direction of the toolkit, and guided its development through detailed instructions and feedback on the results. The implementation was done with Claude (Anthropic) as an AI coding assistant: it wrote the HTML, CSS and JavaScript and this README, and ran the automated tests described above. The equations are standard first- and second-order results, and the validation page checks them against the reference values in the homework sheet. I take responsibility for the submitted work and will explain its code and equations in the demo.
