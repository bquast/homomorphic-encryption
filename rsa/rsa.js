function modpow(base, exp, mod) {
  if (mod === 1n) return 0n;
  let result = 1n;
  base = base % mod;
  while (exp > 0n) {
    if (exp % 2n === 1n) result = result * base % mod;
    exp = exp >> 1n;
    base = base * base % mod;
  }
  return result;
}

function gcd(a, b) { return b === 0n ? a : gcd(b, a % b); }

function modInv(a, m) {
  let [old_r, r] = [a, m];
  let [old_s, s] = [1n, 0n];
  while (r !== 0n) {
    const q = old_r / r;
    [old_r, r] = [r, old_r - q * r];
    [old_s, s] = [s, old_s - q * s];
  }
  return ((old_s % m) + m) % m;
}

function findE(phi) {
  for (let e = 3n; e < phi; e += 2n) {
    if (gcd(e, phi) === 1n) return e;
  }
  return 3n;
}

export function compute(P, Q, m1, m2) {
  const Pb = BigInt(P);
  const Qb = BigInt(Q);
  const m1b = BigInt(m1);
  const m2b = BigInt(m2);

  const n = Pb * Qb;
  const phi = (Pb - 1n) * (Qb - 1n);
  const e = findE(phi);
  let d;
  try { d = modInv(e, phi); } catch { d = 1n; }

  const n2 = n * n;
  const c1 = modpow(m1b, e, n2);
  const c2 = modpow(m2b, e, n2);
  const cp = (c1 * c2) % n2;
  const dc1 = modpow(c1, d, n);
  const dc2 = modpow(c2, d, n);
  const dcp = modpow(cp, d, n);
  const plainprod = (m1b * m2b) % n;

  return {
    n, phi, e, d, n2,
    c1, c2, cp,
    dc1, dc2, dcp,
    plainprod,
    holds: dcp === plainprod,
    warn: m1b >= n || m2b >= n
  };
}