// =========================================
// BFV parameters (Tiny insecure demo params)
// =========================================
// Polynomial modulus: x^N + 1
const N = 8n; 
// Ciphertext modulus (a prime number)
const Q = 12289n; 
// Plaintext modulus (allows integers 0 to 255)
const T = 256n; 

// Scaling factor Delta = floor(Q/T)
const DELTA = Q / T;

// Global state for demo purposes
let secretKey = null;
let publicKey = null;
let cipherA = null;
let cipherB = null;
let cipherResult = null;

// =========================================
// Math & Polynomial Helpers (Ring Arith)
// =========================================

// Correct modulo for negative numbers: (-5) mod 12289 -> 12284
const modQ = (val) => ((val % Q) + Q) % Q;

// Generate a random polynomial with coeffs in [0, Q-1]
const sampleUniformPoly = () => {
    const poly = [];
    for (let i = 0n; i < N; i++) {
        // Insecure random for demo only
        poly.push(BigInt(Math.floor(Math.random() * Number(Q))));
    }
    return poly;
};

// Generate "small" noise polynomial (ternary distribution {-1, 0, 1})
const sampleSmallPoly = () => {
    const poly = [];
    for (let i = 0n; i < N; i++) {
        const rand = Math.random();
        if (rand < 0.33) poly.push(0n);
        else if (rand < 0.66) poly.push(1n);
        else poly.push(Q - 1n); // Represents -1 mod Q
    }
    return poly;
};

// Polynomial Addition: (p1 + p2) mod Q
const polyAdd = (p1, p2) => {
    return p1.map((coeff, i) => modQ(coeff + p2[i]));
};

// Polynomial Subtraction: (p1 - p2) mod Q
const polySub = (p1, p2) => {
    return p1.map((coeff, i) => modQ(coeff - p2[i]));
};

// Polynomial Multiplication in ring Zq[x] / (x^N + 1)
// Standard schoolbook multiplication O(N^2) followed by reduction.
const polyMul = (p1, p2) => {
    // double degree temp buffer
    const res = new Array(Number(N) * 2).fill(0n);
    
    // Convolution
    for (let i = 0; i < Number(N); i++) {
        for (let j = 0; j < Number(N); j++) {
            res[i + j] = res[i + j] + (p1[i] * p2[j]);
        }
    }

    // Reduction modulo (x^N + 1).
    // Rule: x^N = -1, x^(N+1) = -x, etc.
    // Coeffs at indices >= N wrap around and flip sign.
    for (let i = Number(N); i < 2 * Number(N); i++) {
        res[i - Number(N)] = res[i - Number(N)] - res[i];
        res[i] = 0n;
    }

    // Final modulo Q on first N coefficients
    return res.slice(0, Number(N)).map(c => modQ(c));
};


// =========================================
// BFV Core Functions
// =========================================

// 1. Key Generation
// Generates Secret Key (sk) and Public Key (pk = (pk0, pk1))
const keyGen = () => {
    const s = sampleSmallPoly(); // Secret key
    const a = sampleUniformPoly(); // Public parameter part 1
    const e = sampleSmallPoly(); // Error term

    // pk0 = -(a*s + e) mod Q
    const as_plus_e = polyAdd(polyMul(a, s), e);
    const pk0 = as_plus_e.map(c => modQ(-c));
    const pk1 = a;

    return { sk: s, pk: [pk0, pk1] };
};

// 2. Encryption
// Encrypts integer message 'm' into ciphertext ct = (c0, c1)
const encrypt = (m, pk) => {
    const mBig = BigInt(m);
    const [pk0, pk1] = pk;

    // Encode message: map integer m to constant polynomial scaled by Delta
    const scaledM = new Array(Number(N)).fill(0n);
    scaledM[0] = modQ(DELTA * mBig);

    const u = sampleSmallPoly(); // Random ephemeral key
    const e1 = sampleSmallPoly(); // Error 1
    const e2 = sampleSmallPoly(); // Error 2

    // c0 = pk0 * u + e1 + Delta*m
    const c0_temp = polyAdd(polyMul(pk0, u), e1);
    const c0 = polyAdd(c0_temp, scaledM);

    // c1 = pk1 * u + e2
    const c1 = polyAdd(polyMul(pk1, u), e2);

    // Ciphertext is a pair of polynomials
    return [c0, c1];
};

// 3. Decryption
// Decrypts ciphertext ct = (c0, c1) using sk
const decrypt = (ct, sk) => {
    const [c0, c1] = ct;

    // Calculate noisy plaintext polynomial: c0 + c1 * s mod Q
    const noisyPoly = polyAdd(c0, polyMul(c1, sk));

    // Decode: Extract constant coefficient, scale down by T/Q, and round.
    // Rounding formula for integer arithmetic: floor( (noisy_val * T + Q/2) / Q ) mod T
    const noisyConst = noisyPoly[0];
    const numerator = (noisyConst * T) + (Q / 2n);
    const decryptedBig = (numerator / Q) % T;

    return Number(decryptedBig);
};

// 4. Homomorphic Addition
// Adds two ciphertexts: (c0+d0, c1+d1)
const addCiphertexts = (ctA, ctB) => {
    const newC0 = polyAdd(ctA[0], ctB[0]);
    const newC1 = polyAdd(ctA[1], ctB[1]);
    return [newC0, newC1];
};


// =========================================
// UI Glue Code
// =========================================

// Helper for logging to HTML
const log = (title, data) => {
    const logBox = document.getElementById('logOutput');
    const entry = document.createElement('div');
    entry.className = 'log-entry';
    
    // Simple formatter for BigInt arrays for display
    let displayData = data;
    if (Array.isArray(data) && typeof data[0] === 'bigint') {
        displayData = `[${data.slice(0,4).map(b => b.toString()).join(', ')}...] (Poly deg ${N-1n})`;
    } else if (Array.isArray(data) && Array.isArray(data[0])) {
        displayData = `(Poly Pair) C0: [${data[0].slice(0,3).join(',')}...], C1: [${data[1].slice(0,3).join(',')}...]`;
    }

    entry.innerHTML = `<div class="log-title">${title}</div><div class="log-data">${displayData}</div>`;
    
    // Clear default msg if exists
    if (logBox.querySelector('.text-secondary')) {
        logBox.innerHTML = '';
    }
    logBox.appendChild(entry);
    logBox.scrollTop = logBox.scrollHeight;
};

// DOM Elements
const btnKeygen = document.getElementById('btnKeygen');
const btnEncrypt = document.getElementById('btnEncrypt');
const btnAddCipher = document.getElementById('btnAddCipher');
const btnDecrypt = document.getElementById('btnDecrypt');

const stepEnc = document.getElementById('encryptionStep');
const stepComp = document.getElementById('computationStep');
const stepDec = document.getElementById('decryptionStep');
const keyStatus = document.getElementById('keyStatus');
const num1Input = document.getElementById('num1');
const num2Input = document.getElementById('num2');

// Event Listeners
btnKeygen.addEventListener('click', () => {
    log('Action', 'Generating Keys...');
    const keys = keyGen();
    secretKey = keys.sk;
    publicKey = keys.pk;
    
    log('Keys Generated', `Public Key (pk0, pk1) and Secret Key (s) created. N=${N}, Q=${Q}`);
    
    // Update UI state
    keyStatus.classList.remove('hidden');
    btnKeygen.disabled = true;
    btnKeygen.innerText = "Keys Generated";
    stepEnc.classList.add('active');
    btnEncrypt.disabled = false;
});

btnEncrypt.addEventListener('click', () => {
    const valA = parseInt(num1Input.value) || 0;
    const valB = parseInt(num2Input.value) || 0;

    log('Action', `Encrypting values: ${valA} and ${valB}`);

    cipherA = encrypt(valA, publicKey);
    cipherB = encrypt(valB, publicKey);

    log('Encryption Complete', 'Generated Ciphertext A and Ciphertext B.');
    log('Ciphertext A Sample', cipherA);

    // Update UI State
    stepComp.classList.add('active');
    btnAddCipher.disabled = false;
    btnEncrypt.disabled = true;
});

btnAddCipher.addEventListener('click', () => {
    log('Action', 'Computing Homomorphic Addition (CtxA + CtxB) in encrypted domain...');
    
    // The magic happens here. Adding polynomials, not numbers.
    cipherResult = addCiphertexts(cipherA, cipherB);

    log('Computation Complete', 'New resultant ciphertext created.');
    log('Result Ciphertext Sample', cipherResult);

    // Update UI State
    stepDec.classList.add('active');
    btnDecrypt.disabled = false;
    btnAddCipher.disabled = true;
});

btnDecrypt.addEventListener('click', () => {
    log('Action', 'Decrypting resultant ciphertext using Secret Key...');

    const result = decrypt(cipherResult, secretKey);

    const expected = (parseInt(num1Input.value) + parseInt(num2Input.value)) % Number(T);
    
    log('Decryption Successful', `<span class="log-highlight">Final Result: ${result}</span>`);
    log('Verification', `Expected value ( (${num1Input.value} + ${num2Input.value}) mod ${T} ) is: ${expected}`);
    
    btnDecrypt.disabled = true;

    if (result === expected) {
         log('SUCCESS', 'Homomorphic addition matched expected plaintext sum.');
    } else {
         log('FAILURE', 'Decryption did not match. Noise budget may have been exceeded.');
    }
});

// Initial state
btnEncrypt.disabled = true;
btnAddCipher.disabled = true;
btnDecrypt.disabled = true;