/**
 * Complex Number Math Class
 */
class Complex {
    constructor(re, im = 0) {
        this.re = re;
        this.im = im;
    }

    add(c) {
        return new Complex(this.re + c.re, this.im + c.im);
    }

    sub(c) {
        return new Complex(this.re - c.re, this.im - c.im);
    }

    mul(c) {
        if (typeof c === 'number') {
            return new Complex(this.re * c, this.im * c);
        }
        return new Complex(
            this.re * c.re - this.im * c.im,
            this.re * c.im + this.im * c.re
        );
    }

    div(c) {
        if (typeof c === 'number') {
            return new Complex(this.re / c, this.im / c);
        }
        const denom = c.re * c.re + c.im * c.im;
        return new Complex(
            (this.re * c.re + this.im * c.im) / denom,
            (this.im * c.re - this.re * c.im) / denom
        );
    }

    conjugate() {
        return new Complex(this.re, -this.im);
    }

    magnitude() {
        return Math.sqrt(this.re * this.re + this.im * this.im);
    }

    pow(n) {
        let res = new Complex(1, 0);
        let base = this;
        for (let i = 0; i < n; i++) {
            res = res.mul(base);
        }
        return res;
    }

    static exp(theta) {
        return new Complex(Math.cos(theta), Math.sin(theta));
    }
}

/**
 * Linear Algebra Tools
 */

// Gaussian elimination for complex matrices: Solves A * x = b
function solve(A, b) {
    const n = A.length;
    let mat = A.map((row, i) => [...row, b[i]]); // Create Augmented matrix

    for (let i = 0; i < n; i++) {
        // Find pivot
        let max = i;
        for (let j = i + 1; j < n; j++) {
            if (mat[j][i].magnitude() > mat[max][i].magnitude()) {
                max = j;
            }
        }

        // Swap rows
        let temp = mat[i];
        mat[i] = mat[max];
        mat[max] = temp;

        // Eliminate
        let pivot = mat[i][i];
        for (let j = i; j <= n; j++) {
            mat[i][j] = mat[i][j].div(pivot);
        }
        for (let j = 0; j < n; j++) {
            if (i !== j) {
                let factor = mat[j][i];
                for (let k = i; k <= n; k++) {
                    mat[j][k] = mat[j][k].sub(factor.mul(mat[i][k]));
                }
            }
        }
    }
    return mat.map(row => row[n]);
}

function transpose(mat) {
    return mat[0].map((_, colIndex) => mat.map(row => row[colIndex]));
}

// Complex dot product (np.vdot behavior: sum(conj(a_i) * b_i))
function vdot(a, b) {
    let sum = new Complex(0, 0);
    for (let i = 0; i < a.length; i++) {
        sum = sum.add(a[i].conjugate().mul(b[i]));
    }
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
 * CKKS Encoder Core logic
 */

function vandermonde(xi, M) {
    const N = Math.floor(M / 2);
    let matrix = [];
    for (let i = 0; i < N; i++) {
        let root = xi.pow(2 * i + 1);
        let row = [];
        for (let j = 0; j < N; j++) {
            row.push(root.pow(j));
        }
        matrix.push(row);
    }
    return matrix;
}

function sigma_inverse(xi, M, b) {
    const A = vandermonde(xi, M);
    return solve(A, b); // Returns array of Complex polynomials coefficients
}

function sigma(xi, M, p_coeffs) {
    const N = Math.floor(M / 2);
    let outputs = [];
    for (let i = 0; i < N; i++) {
        let root = xi.pow(2 * i + 1);
        let sum = new Complex(0, 0);
        for (let j = 0; j < p_coeffs.length; j++) {
            sum = sum.add(p_coeffs[j].mul(root.pow(j)));
        }
        outputs.push(sum);
    }
    return outputs;
}

function create_sigma_R_basis(xi, M) {
    return transpose(vandermonde(xi, M));
}

function compute_basis_coordinates(sigma_R_basis, z) {
    return sigma_R_basis.map(b => {
        let num = vdot(z, b);
        let den = vdot(b, b); // Always real
        return num.div(den).re; // Extract real part
    });
}

function coordinate_wise_random_rounding(coordinates) {
    return coordinates.map(c => {
        let r = c - Math.floor(c);
        let rand = Math.random();
        return rand < (1 - r) ? Math.floor(c) : Math.ceil(c);
    });
}

function sigma_R_discretization(xi, M, z) {
    const sigma_R_basis = create_sigma_R_basis(xi, M);
    const coordinates = compute_basis_coordinates(sigma_R_basis, z);
    const rounded_coordinates = coordinate_wise_random_rounding(coordinates);
    // sigma_R_basis is Vandermonde^T. Transposing it back gives Vandermonde.
    const V = vandermonde(xi, M);
    return matrixVectorMul(V, rounded_coordinates);
}

function pi(M, z) {
    const N = Math.floor(M / 4);
    return z.slice(0, N);
}

function pi_inverse(z) {
    let z_conj_rev = [...z].reverse().map(x => x.conjugate());
    return z.concat(z_conj_rev);
}

function encode(xi, M, scale, z) {
    let pi_z = pi_inverse(z);
    let scaled_pi_z = pi_z.map(val => val.mul(scale));
    let rounded_scale_pi_zi = sigma_R_discretization(xi, M, scaled_pi_z);
    let p_coeffs_complex = sigma_inverse(xi, M, rounded_scale_pi_zi);
    return p_coeffs_complex.map(c => Math.round(c.re)); // Integer polynomial coefficients
}

function decode(xi, M, scale, p_coeffs) {
    let rescaled_p = p_coeffs.map(c => new Complex(c / scale, 0));
    let z = sigma(xi, M, rescaled_p);
    return pi(M, z);
}


/**
 * DOM Interaction & Formatting
 */

function formatPoly(coeffs) {
    let terms = [];
    for (let i = 0; i < coeffs.length; i++) {
        if (coeffs[i] === 0) continue;
        let term = '';
        if (i === 0) term = `${coeffs[i]}`;
        else if (i === 1) term = `${coeffs[i]}X`;
        else term = `${coeffs[i]}X^${i}`;
        terms.push(term);
    }
    if (terms.length === 0) return '0';
    return terms.join(' + ').replace(/\+ -/g, '- ');
}

function formatVector(vec) {
    return `[\n  ${vec.map(c => `${c.re.toFixed(4).padStart(8)} ${c.im >= 0 ? '+' : '-'} ${Math.abs(c.im).toFixed(4).padStart(8)}i`).join(',\n  ')}\n]`;
}

// Initial defaults to match Python example
const defaultZ = [
    new Complex(3, 4),
    new Complex(2, -1)
];

function updateVectorInputs() {
    const M = parseInt(document.getElementById('input-m').value) || 8;
    const N = Math.max(1, Math.floor(M / 4));
    
    document.getElementById('vector-size-badge').innerText = `N = ${N}`;
    const container = document.getElementById('vector-inputs-container');
    container.innerHTML = '';

    for(let i = 0; i < N; i++) {
        const val = defaultZ[i] || new Complex(0, 0);
        
        const row = document.createElement('div');
        row.className = 'vector-row';
        row.innerHTML = `
            <div class="vector-index">z[${i}]</div>
            <div class="vector-part">
                <input type="number" step="any" class="input z-re" value="${val.re}" placeholder="Real">
            </div>
            <div style="color: var(--neutral); font-weight: 500;">+</div>
            <div class="vector-part">
                <input type="number" step="any" class="input z-im" value="${val.im}" placeholder="Imag">
                <span style="color: var(--neutral); font-family: 'JetBrains Mono', monospace;">i</span>
            </div>
        `;
        container.appendChild(row);
    }
}

function executeCKKS() {
    const M = parseInt(document.getElementById('input-m').value) || 8;
    const scale = parseInt(document.getElementById('input-scale').value) || 64;
    
    // Parse vector inputs
    const reInputs = document.querySelectorAll('.z-re');
    const imInputs = document.querySelectorAll('.z-im');
    let z = [];
    for(let i=0; i<reInputs.length; i++) {
        z.push(new Complex(
            parseFloat(reInputs[i].value) || 0,
            parseFloat(imInputs[i].value) || 0
        ));
    }

    // Mathematical Constants setup
    const theta = 2 * Math.PI / M;
    const xi = Complex.exp(theta);

    // 1. Encode
    const startEncode = performance.now();
    const p_coeffs = encode(xi, M, scale, z);
    const endEncode = performance.now();

    // 2. Decode
    const startDecode = performance.now();
    const decoded_z = decode(xi, M, scale, p_coeffs);
    const endDecode = performance.now();

    // Update UI
    document.getElementById('output-poly').innerText = 
        `// Encoding took ${(endEncode - startEncode).toFixed(2)}ms\np(X) = ${formatPoly(p_coeffs)}`;
        
    document.getElementById('output-decoded').innerText = 
        `// Decoding took ${(endDecode - startDecode).toFixed(2)}ms\nDecoded z = ${formatVector(decoded_z)}`;

    // Calculate Error
    let errorVec = [];
    for(let i=0; i<z.length; i++) {
        errorVec.push(new Complex(
            Math.abs(z[i].re - decoded_z[i].re),
            Math.abs(z[i].im - decoded_z[i].im)
        ));
    }
    document.getElementById('output-error').innerText = 
        `Δ = ${formatVector(errorVec)}`;
}

// Event Listeners
document.getElementById('input-m').addEventListener('change', updateVectorInputs);
document.getElementById('run-btn').addEventListener('click', executeCKKS);

// Initialize App
document.addEventListener('DOMContentLoaded', () => {
    updateVectorInputs();
    executeCKKS(); // run immediately with defaults
});