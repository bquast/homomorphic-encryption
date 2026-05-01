/**
 * Math Core - Complex Numbers and Linear Algebra
 */
class Complex {
    constructor(re, im = 0) { this.re = re; this.im = im; }
    add(c) { return new Complex(this.re + c.re, this.im + c.im); }
    sub(c) { return new Complex(this.re - c.re, this.im - c.im); }
    mul(c) {
        if (typeof c === 'number') return new Complex(this.re * c, this.im * c);
        return new Complex(this.re * c.re - this.im * c.im, this.re * c.im + this.im * c.re);
    }
    div(c) {
        if (typeof c === 'number') return new Complex(this.re / c, this.im / c);
        const denom = c.re * c.re + c.im * c.im;
        return new Complex((this.re * c.re + this.im * c.im) / denom, (this.im * c.re - this.re * c.im) / denom);
    }
    conjugate() { return new Complex(this.re, -this.im); }
    magnitude() { return Math.sqrt(this.re * this.re + this.im * this.im); }
    pow(n) {
        let res = new Complex(1, 0), base = this;
        for (let i = 0; i < n; i++) res = res.mul(base);
        return res;
    }
    static exp(theta) { return new Complex(Math.cos(theta), Math.sin(theta)); }
}

function solve(A, b) {
    const n = A.length;
    let mat = A.map((row, i) => [...row, b[i]]);
    for (let i = 0; i < n; i++) {
        let max = i;
        for (let j = i + 1; j < n; j++) if (mat[j][i].magnitude() > mat[max][i].magnitude()) max = j;
        let temp = mat[i]; mat[i] = mat[max]; mat[max] = temp;
        let pivot = mat[i][i];
        for (let j = i; j <= n; j++) mat[i][j] = mat[i][j].div(pivot);
        for (let j = 0; j < n; j++) {
            if (i !== j) {
                let factor = mat[j][i];
                for (let k = i; k <= n; k++) mat[j][k] = mat[j][k].sub(factor.mul(mat[i][k]));
            }
        }
    }
    return mat.map(row => row[n]);
}

function transpose(mat) { return mat[0].map((_, col) => mat.map(row => row[col])); }
function vdot(a, b) {
    let sum = new Complex(0, 0);
    for (let i = 0; i < a.length; i++) sum = sum.add(a[i].conjugate().mul(b[i]));
    return sum;
}
function matrixVectorMul(mat, vec) {
    let res = [];
    for (let i = 0; i < mat.length; i++) {
        let sum = new Complex(0, 0);
        for (let j = 0; j < vec.length; j++) {
            let val = vec[j];
            if (typeof val === 'number') val = new Complex(val, 0);
            sum = sum.add(mat[i][j].mul(val));
        }
        res.push(sum);
    }
    return res;
}

/**
 * Pipeline State & RLWE Constants
 */
let encodedPolyCoeffs = null;
let sk = null;
let pk = null;
let ciphertext = null;
let decryptedPolyCoeffs = null;

let N_ring = 4;
const Q = 1073741824n; // 2^30 Modulus for polynomial ring

/**
 * RLWE Math helpers
 */
function modQ(val) { let v = val % Q; if (v < 0n) v += Q; return v; }
function centerModQ(val) { let v = modQ(val); if (v > Q / 2n) v -= Q; return v; }

function polyAdd(p1, p2) { return p1.map((c, i) => modQ(c + p2[i])); }
function polyMul(p1, p2) {
    const res = new Array(N_ring * 2).fill(0n);
    for (let i = 0; i < N_ring; i++) {
        for (let j = 0; j < N_ring; j++) res[i + j] += p1[i] * p2[j];
    }
    for (let i = N_ring; i < 2 * N_ring; i++) {
        res[i - N_ring] -= res[i]; // Reduction X^N = -1
        res[i] = 0n;
    }
    return res.slice(0, N_ring).map(modQ);
}

function sampleUniformPoly() { return Array.from({length: N_ring}, () => BigInt(Math.floor(Math.random() * Number(Q)))); }
function sampleSmallPoly() {
    return Array.from({length: N_ring}, () => {
        let r = Math.random(); return r < 0.33 ? 0n : r < 0.66 ? 1n : Q - 1n; // Ternary [-1, 0, 1]
    });
}

/**
 * CKKS Encoder Core Logic
 */
function vandermonde(xi, M) {
    const N = Math.floor(M / 2);
    let matrix = [];
    for (let i = 0; i < N; i++) {
        let root = xi.pow(2 * i + 1), row = [];
        for (let j = 0; j < N; j++) row.push(root.pow(j));
        matrix.push(row);
    }
    return matrix;
}
function pi_inverse(z) { return z.concat([...z].reverse().map(x => x.conjugate())); }
function pi(M, z) { return z.slice(0, Math.floor(M / 4)); }
function sigma(xi, M, p_coeffs) {
    const N = Math.floor(M / 2);
    return Array.from({length: N}, (_, i) => {
        let root = xi.pow(2 * i + 1), sum = new Complex(0, 0);
        for (let j = 0; j < p_coeffs.length; j++) sum = sum.add(p_coeffs[j].mul(root.pow(j)));
        return sum;
    });
}

/**
 * UI Log Utilities
 */
function log(title, data, highlight = false) {
    const logBox = document.getElementById('log-output');
    const entry = document.createElement('div');
    entry.className = 'log-entry' + (highlight ? ' highlight-entry' : '');
    
    let displayData = data;
    if (Array.isArray(data) && data.length > 0 && data[0] instanceof Complex) displayData = formatVector(data);
    entry.innerHTML = `<div class="log-title">${title}</div><div class="log-data">${displayData}</div>`;
    logBox.appendChild(entry);
    logBox.scrollTop = logBox.scrollHeight;
}

function formatVector(vec) {
    if (vec.length === 0) return '[]';
    return `[\n  ${vec.map(c => `${c.re.toFixed(4).padStart(8)} ${c.im >= 0 ? '+' : '-'} ${Math.abs(c.im).toFixed(4).padStart(8)}i`).join(',\n  ')}\n]`;
}

function formatPoly(coeffs, isBigInt = false) {
    let terms = [];
    for (let i = 0; i < coeffs.length; i++) {
        let val = isBigInt ? centerModQ(coeffs[i]) : coeffs[i];
        let c = isBigInt ? Number(val) : val;
        if (c === 0 && terms.length > 0) continue;
        let term = i === 0 ? `${c}` : i === 1 ? `${c}X` : `${c}X^${i}`;
        terms.push(term);
    }
    return terms.length === 0 ? '0' : terms.join(' + ').replace(/\+ -/g, '- ');
}

/**
 * Pipeline Execution Steps
 */
function doEncode() {
    const logBox = document.getElementById('log-output');
    if (logBox.querySelector('p')) logBox.innerHTML = ''; 

    const M = parseInt(document.getElementById('param-m').value);
    const scale = parseFloat(document.getElementById('param-scale').value);
    N_ring = M / 2;

    const N_vec = Math.max(1, Math.floor(M / 4));
    let z = [];
    for(let i = 0; i < N_vec; i++) {
        z.push(new Complex(parseFloat(document.getElementById(`z-re-${i}`).value), parseFloat(document.getElementById(`z-im-${i}`).value)));
    }

    log("--- PHASE 1: CKKS ENCODING ---", "Transforming original complex vector to polynomial via discretized coordinate mapping.");
    log("1. Original Vector z", formatVector(z));
    
    let pi_z = pi_inverse(z);
    log("2. π⁻¹(z) (Symmetric Extension)", formatVector(pi_z));
    
    let scaled_pi_z = pi_z.map(v => v.mul(scale));
    log(`3. Scaled by Δ (${scale})`, formatVector(scaled_pi_z));
    
    const theta = 2 * Math.PI / M;
    const xi = Complex.exp(theta);
    
    const sigma_R_basis = transpose(vandermonde(xi, M));
    const coordinates = sigma_R_basis.map(b => vdot(scaled_pi_z, b).div(vdot(b, b).re).re);
    log("4. Canonical Basis Coordinates", `[${coordinates.map(c => c.toFixed(4)).join(', ')}]`);
    
    const rounded_coordinates = coordinates.map(c => {
        let r = c - Math.floor(c);
        return Math.random() < (1 - r) ? Math.floor(c) : Math.ceil(c);
    });
    log("5. Randomized Int Rounding", `[${rounded_coordinates.join(', ')}]`);
    
    const rounded_scale_pi_zi = matrixVectorMul(transpose(sigma_R_basis), rounded_coordinates);
    const p_coeffs_complex = solve(vandermonde(xi, M), rounded_scale_pi_zi);
    encodedPolyCoeffs = p_coeffs_complex.map(c => Math.round(c.re));
    
    log("6. Final Encoded Polynomial m(X)", formatPoly(encodedPolyCoeffs), true);

    document.getElementById('step-keygen').classList.remove('disabled-overlay');
    document.getElementById('btn-encode').disabled = true;
}

function doKeyGen() {
    log("--- PHASE 2: RLWE KEYGEN ---", `Generating ring polynomials in ℤq[X]/(X^${N_ring}+1). Modulus Q=${Q}`);
    
    sk = sampleSmallPoly();
    const a = sampleUniformPoly();
    const e = sampleSmallPoly();
    
    const as_e = polyAdd(polyMul(a, sk), e);
    const pk0 = as_e.map(c => modQ(-c));
    pk = [pk0, a];
    
    log("Secret Key (s) — ternary noise", formatPoly(sk, true));
    log("Error (e) — ternary noise", formatPoly(e, true));
    log("Public Key (pk0)", "[-a * s + e mod Q] -> (Uniform random poly)");
    log("Public Key (pk1)", "[a] -> (Uniform random poly)", true);
    
    document.getElementById('step-encrypt').classList.remove('disabled-overlay');
    document.getElementById('btn-keygen').disabled = true;
}

function doEncrypt() {
    log("--- PHASE 3: RLWE ENCRYPTION ---", "Masking m(X) with public key and ephemeral noise.");
    
    const m = encodedPolyCoeffs.map(c => modQ(BigInt(c)));
    const u = sampleSmallPoly(), e1 = sampleSmallPoly(), e2 = sampleSmallPoly();
    
    let c0 = polyAdd(polyAdd(polyMul(pk[0], u), e1), m);
    let c1 = polyAdd(polyMul(pk[1], u), e2);
    ciphertext = [c0, c1];
    
    log("Ephemeral key & Noise (u, e1, e2)", "Sampled from small distribution.");
    log("Ciphertext c0", `pk0 * u + e1 + m(X) mod Q`);
    log("Ciphertext c1", `pk1 * u + e2 mod Q`, true);
    
    document.getElementById('step-decrypt').classList.remove('disabled-overlay');
    document.getElementById('btn-encrypt').disabled = true;
}

function doDecrypt() {
    log("--- PHASE 4: RLWE DECRYPTION ---", "Recovering approximate polynomial using secret key.");
    
    let m_noisy = polyAdd(ciphertext[0], polyMul(ciphertext[1], sk));
    let m_recovered = m_noisy.map(centerModQ); // Fix negative wrapping

    decryptedPolyCoeffs = m_recovered.map(c => Number(c));
    
    let noise = decryptedPolyCoeffs.map((c, i) => c - encodedPolyCoeffs[i]);
    
    log("Decrypted Noisy Polynomial m'(X)", formatPoly(decryptedPolyCoeffs));
    log("RLWE Accumulated Noise (m' - m)", formatPoly(noise), true);
    
    document.getElementById('step-decode').classList.remove('disabled-overlay');
    document.getElementById('btn-decrypt').disabled = true;
}

function doDecode() {
    log("--- PHASE 5: CKKS DECODING ---", "Evaluating m'(X) on canonical roots and reversing encoding scaling.");
    
    const M = parseInt(document.getElementById('param-m').value);
    const scale = parseFloat(document.getElementById('param-scale').value);
    const xi = Complex.exp(2 * Math.PI / M);
    
    let rescaled_p = decryptedPolyCoeffs.map(c => new Complex(c / scale, 0));
    log(`1. Scale Down by 1/Δ (1/${scale})`, `m(X) coefficients divided by ${scale}`);
    
    let z_ext = sigma(xi, M, rescaled_p);
    log("2. Matrix Sigma Eval (σ)", formatVector(z_ext));
    
    let z_final = pi(M, z_ext);
    log("3. Extracted Final Vector (π(z'))", formatVector(z_final), true);
    
    document.getElementById('btn-decode').disabled = true;
}

/**
 * Event Listeners and Parameter Handlers
 */
function updateVectorInputs() {
    const M = parseInt(document.getElementById('param-m').value) || 8;
    const N_vec = Math.max(1, Math.floor(M / 4));
    const container = document.getElementById('vector-inputs');
    container.innerHTML = '';
    
    const defaultZ = [new Complex(3, 4), new Complex(2, -1), new Complex(1.5, 0.5), new Complex(-2, -3)];
    
    for(let i = 0; i < N_vec; i++) {
        const val = defaultZ[i] || new Complex(1, 1);
        const div = document.createElement('div');
        div.className = 'input-group-row';
        div.innerHTML = `
            <div class="input-group">
                <label>z[${i}] Real</label>
                <input type="number" id="z-re-${i}" class="input z-input" value="${val.re}">
            </div>
            <div class="input-group">
                <label>z[${i}] Imag</label>
                <input type="number" id="z-im-${i}" class="input z-input" value="${val.im}">
            </div>
        `;
        container.appendChild(div);
    }
    resetPipeline();
}

function resetPipeline() {
    document.getElementById('step-keygen').classList.add('disabled-overlay');
    document.getElementById('step-encrypt').classList.add('disabled-overlay');
    document.getElementById('step-decrypt').classList.add('disabled-overlay');
    document.getElementById('step-decode').classList.add('disabled-overlay');
    
    document.getElementById('btn-encode').disabled = false;
    document.getElementById('btn-keygen').disabled = false;
    document.getElementById('btn-encrypt').disabled = false;
    document.getElementById('btn-decrypt').disabled = false;
    document.getElementById('btn-decode').disabled = false;
    
    document.getElementById('log-output').innerHTML = '<p class="text-secondary">Adjust parameters and click "Encode Vector" to begin...</p>';
}

document.addEventListener('DOMContentLoaded', () => {
    document.getElementById('param-m').addEventListener('change', updateVectorInputs);
    document.getElementById('param-scale').addEventListener('change', resetPipeline);
    
    document.getElementById('btn-encode').addEventListener('click', doEncode);
    document.getElementById('btn-keygen').addEventListener('click', doKeyGen);
    document.getElementById('btn-encrypt').addEventListener('click', doEncrypt);
    document.getElementById('btn-decrypt').addEventListener('click', doDecrypt);
    document.getElementById('btn-decode').addEventListener('click', doDecode);
    
    updateVectorInputs(); // initialize dynamically based on default params
});