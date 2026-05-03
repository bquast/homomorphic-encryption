// =========================================
// BGV parameters (From Python reference)
// =========================================
// Polynomial modulus: x^N + 1
const N = 16n; 
// Ciphertext modulus
const Q = 868n; 
// Plaintext modulus (tiny prime)
const P = 7n; 


// Global state for demo purposes
let secretKey = null;
let publicKey = null;
let cipherA = null;
let cipherB = null;
let cipherResult = null;

// =========================================
// Math & Polynomial Helpers (Ring Arith)
// =========================================

// Correct modulo Q (ciphertext space): (-5) mod 868 -> 863
const modQ = (val) => ((val % Q) + Q) % Q;

// Correct modulo P (plaintext space): (-2) mod 7 -> 5
const modP = (val) => ((val % P) + P) % P;

// Generate a random polynomial with coeffs in [0, Q-1]
const sampleUniformPoly = () => {
    const poly = [];
    for (let i = 0n; i < N; i++) {
        poly.push(BigInt(Math.floor(Math.random() * Number(Q))));
    }
    return poly;
};

// Generate ternary polynomial {-1, 0, 1} for keys
const sampleTernaryPoly = () => {
    const poly = [];
    for (let i = 0n; i < N; i++) {
        const rand = Math.random();
        if (rand < 0.33) poly.push(0n);
        else if (rand < 0.66) poly.push(1n);
        else poly.push(Q - 1n); // Represents -1 mod Q
    }
    return poly;
};

// Generate error polynomial.
// Python used Gaussian normal(0, N/3). For N=16, sigma is ~5.
// We use a uniform distribution over approx [-3, 3] for simplicity in JS.
const sampleErrorPoly = () => {
    const poly = [];
    for (let i = 0n; i < N; i++) {
        // random int between 0 and 6, shift down by 3 to get [-3, 3]
        let val = BigInt(Math.floor(Math.random() * 7)) - 3n;
        poly.push(modQ(val));
    }
    return poly;
};

// Polynomial Addition: (p1 + p2) mod Q
const polyAdd = (p1, p2) => {
    return p1.map((coeff, i) => modQ(coeff + p2[i]));
};

// Polynomial Multiplication in ring Zq[x] / (x^N + 1)
// Standard O(N^2) multiplication followed by reduction.
const polyMul = (p1, p2) => {
    const res = new Array(Number(N) * 2).fill(0n);
    
    // Convolution
    for (let i = 0; i < Number(N); i++) {
        for (let j = 0; j < Number(N); j++) {
            res[i + j] = res[i + j] + (p1[i] * p2[j]);
        }
    }

    // Reduction modulo (x^N + 1).
    // Rule: x^N = -1. Coeffs at indices >= N flip sign and wrap.
    for (let i = Number(N); i < 2 * Number(N); i++) {
        res[i - Number(N)] = res[i - Number(N)] - res[i];
        res[i] = 0n;
    }

    // Final modulo Q on first N coefficients
    return res.slice(0, Number(N)).map(c => modQ(c));
};


// =========================================
// BGV Core Functions
// =========================================

// 1. Key Generation
// Generates Secret Key (s) and Public Key pk = (pk1, pk2)
const keyGen = () => {
    const s = sampleTernaryPoly(); // Secret key
    const a = sampleUniformPoly(); // Public parameter 'a' (pk2)
    const e = sampleErrorPoly();   // Error term

    // BGV Public Key structure:
    // pk1 = -(a*s + P*e) mod Q
    // pk2 = a

    const as = polyMul(a, s);
    // Scale error by plaintext modulus P
    const Pe = e.map(c => modQ(c * P));

    const as_plus_Pe = polyAdd(as, Pe);
    
    const pk1 = as_plus_Pe.map(c => modQ(-c));
    const pk2 = a;

    // In this demo implementation mapping:
    // pk[0] is the encrypted secret part (pk1 in python)
    // pk[1] is the uniform part 'a' (pk2 in python)
    return { sk: s, pk: [pk1, pk2] };
};

// 2. Encryption
// Encrypts integer message 'm' (0 to P-1) into ciphertext ct = (ct1, ct2)
const encrypt = (m, pk) => {
    const mBig = BigInt(m);
    // pk[0] is pk1, pk[1] is pk2 (a)
    const [pk1, pk2] = pk;

    const u = sampleTernaryPoly(); // Ephemeral key
    const e1 = sampleErrorPoly(); 
    const e2 = sampleErrorPoly();

    // Embed message into constant term of a polynomial
    const mPoly = new Array(Number(N)).fill(0n);
    mPoly[0] = modQ(mBig);

    // Scale errors by P
    const Pe1 = e1.map(c => modQ(c * P));
    const Pe2 = e2.map(c => modQ(c * P));

    // BGV Encryption:
    // ct1 = pk1*u + P*e1 + m mod Q  (contains message)
    // ct2 = pk2*u + P*e2 mod Q      (contains 'a' part)

    // Calculate ct1
    const pk1u = polyMul(pk1, u);
    const ct1_temp = polyAdd(pk1u, Pe1);
    const ct1 = polyAdd(ct1_temp, mPoly);

    // Calculate ct2
    const pk2u = polyMul(pk2, u);
    const ct2 = polyAdd(pk2u, Pe2);

    return [ct1, ct2];
};

// 3. Decryption
// Decrypts ciphertext ct = (ct1, ct2) using sk
const decrypt = (ct, sk) => {
    const [ct1, ct2] = ct;

    // BGV Decryption attempt:
    // noisy_poly = ct1 + ct2*s mod Q
    const ct2s = polyMul(ct2, sk);
    const noisyPoly = polyAdd(ct1, ct2s);

    // The result is (m + P*error) mod Q.
    // Since Q (868) is much larger than P (7) times small errors, 
    // we can recover m by taking the constant coefficient modulo P.
    
    // Extract constant term. It's currently in [0, Q-1].
    const noisyConst = noisyPoly[0];

    // Final reduction modulo P to recover message.
    const result = modP(noisyConst);

    return Number(result);
};

// 4. Homomorphic Addition
// Adds two ciphertexts component-wise: (ctA1+ctB1, ctA2+ctB2)
const addCiphertexts = (ctA, ctB) => {
    const newCt1 = polyAdd(ctA[0], ctB[0]);
    const newCt2 = polyAdd(ctA[1], ctB[1]);
    return [newCt1, newCt2];
};


// =========================================
// UI Glue Code (Identical structure to BFV)
// =========================================

// Helper for logging to HTML
const log = (title, data) => {
    const logBox = document.getElementById('logOutput');
    const entry = document.createElement('div');
    entry.className = 'log-entry';
    
    let displayData = data;
    if (Array.isArray(data) && typeof data[0] === 'bigint') {
        displayData = `[${data.slice(0,4).map(b => b.toString()).join(', ')}...] (Poly deg ${N-1n})`;
    } else if (Array.isArray(data) && Array.isArray(data[0])) {
        // Displaying a ciphertext pair (ct1, ct2)
        displayData = `(Poly Pair)<br>CT1 (msg part): [${data[0].slice(0,4).join(',')}...]<br>CT2 ('a' part): [${data[1].slice(0,4).join(',')}...]`;
    }

    entry.innerHTML = `<div class="log-title">${title}</div><div class="log-data">${displayData}</div>`;
    
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
    
    log('Keys Generated', `Public Key (pk1, pk2) and Secret Key (s) created. N=${N}, Q=${Q}, P=${P}`);
    
    keyStatus.classList.remove('hidden');
    btnKeygen.disabled = true;
    btnKeygen.innerText = "Keys Generated";
    stepEnc.classList.add('active');
    btnEncrypt.disabled = false;
});

btnEncrypt.addEventListener('click', () => {
    // Ensure inputs are within 0-(P-1)
    let valA = parseInt(num1Input.value) || 0;
    let valB = parseInt(num2Input.value) || 0;
    valA = Math.max(0, Math.min(Number(P)-1, valA));
    valB = Math.max(0, Math.min(Number(P)-1, valB));
    num1Input.value = valA;
    num2Input.value = valB;

    log('Action', `Encrypting values: ${valA} and ${valB}`);

    cipherA = encrypt(valA, publicKey);
    cipherB = encrypt(valB, publicKey);

    log('Encryption Complete', 'Generated Ciphertext A and Ciphertext B.');
    log('Ciphertext A Sample', cipherA);

    stepComp.classList.add('active');
    btnAddCipher.disabled = false;
    btnEncrypt.disabled = true;
});

btnAddCipher.addEventListener('click', () => {
    log('Action', 'Computing Homomorphic Addition (CtxA + CtxB) in encrypted domain...');
    
    // Performing addition on polynomials.
    cipherResult = addCiphertexts(cipherA, cipherB);

    log('Computation Complete', 'New resultant ciphertext created.');
    log('Result Ciphertext Sample', cipherResult);

    stepDec.classList.add('active');
    btnDecrypt.disabled = false;
    btnAddCipher.disabled = true;
});

btnDecrypt.addEventListener('click', () => {
    log('Action', 'Decrypting resultant ciphertext using Secret Key...');

    const result = decrypt(cipherResult, secretKey);

    const inputA = parseInt(num1Input.value);
    const inputB = parseInt(num2Input.value);
    // Expected result is (A + B) mod P
    const expected = (inputA + inputB) % Number(P);
    
    log('Decryption Successful', `<span class="log-highlight">Final Result: ${result}</span>`);
    log('Verification', `Expected value ( (${inputA} + ${inputB}) mod ${P} ) is: ${expected}`);
    
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