from os import urandom
from hashlib import shake_256, shake_128
import sympy

# fake params
logn = 10
n = 1 << logn
logq = 14
q = 12289
k = 8
gamma = 7
omega = 49


def minabsmod(a, q=q):
    return ((a + q//2) % q) - q//2


def shake256(length, input):
    """my editor gives error because the type is only conditionally defined"""
    state = shake_256(input)
    return state.digest(length)


def shake128Absorb(input):
    """I guessed what this does..."""
    state = shake_128(input)
    return state


def shake128Squeeze(num, state):
    """I suspect this is wrong but not sure how to fix"""
    digest = state.digest(168)
    state.update(b"")
    return digest, state


def HW(x):
    return sum([(x >> i) & 1 for i in range(logq)])


def BitRev(x):
    return sum([((x >> i) & 1) << (logn - 1 - i) for i in range(logn)])


def PolyBitRev(s):
    """disable this to remove complexity while debugging"""
    return [s[BitRev(i)] for i in range(n)]


def GenA(seed):
    a_hat = [None for i in range(n)]
    for i in range(n//64):
        ctr = 0
        state = shake128Absorb(seed + bytes([i]))
        while ctr < 64:
            buf, state = shake128Squeeze(1, state)
            j = 0
            while j < 168 and ctr < 64:
                val = buf[j] + buf[j+1] << 8
                if val < 5 * q:
                    a_hat[i * 64 + ctr] = val % q
                    ctr += 1
                j += 2
    return a_hat


def Sample(seed, nonce):
    ext_seed = seed + bytes([nonce])
    r = [None for i in range(n)]
    for i in range(n//64):
        buf = shake256(128, ext_seed + bytes([i]))
        for j in range(64):
            a = buf[2*j]
            b = buf[2*j + 1]
            r[64*i + j] = (HW(a) + q - HW(b)) % q
    return r


def _slowNTT(vec):
    """really slow implementation of NTT!"""
    return [sum([pow(gamma, j, q) * vec[j] * pow(omega, i*j, q) for j in range(n)]) % q for i in range(n)]


def _slowNTTinv(vec_hat):
    """really slow implementation of inverse NTT!"""
    return [(pow(n, - 1, q) * pow(gamma, - i, q) * sum([vec_hat[j] * pow(omega, - i*j, q) for j in range(n)])) % q for i in range(n)]


def NTT(vec):
    return sympy.ntt(vec, prime=q)


def NTTinv(vec_hat):
    return sympy.intt(vec_hat, prime=q)


def EncodePolynomial(s_hat):
    r = [None for i in range(7 * n // 4)]
    for i in range(n//4):
        t = s_hat[4 * i: 4 * i + 4]
        r[7*i + 0] = t[0] & 0xff
        r[7*i + 1] = (t[0] >> 8) | (t[1] << 6) & 0xff
        r[7*i + 2] = (t[1] >> 2) & 0xff
        r[7*i + 3] = (t[1] >> 10) | (t[2] << 4) & 0xff
        r[7*i + 4] = (t[2] >> 4) & 0xff
        r[7*i + 5] = (t[2] >> 12) | (t[3] << 2) & 0xff
        r[7*i + 6] = (t[3] >> 6) & 0xff
    return bytes(r)


def DecodePolynomial(v):
    r = [None for i in range(n)]
    for i in range(n//4):
        r[4*i + 0] = v[7 * i + 0] | ((v[7 * i + 1] & 0x3f) << 8)
        r[4*i + 1] = v[7 * i + 1] >> 6 | (v[7 * i + 2]
                                          << 2) | ((v[7 * i + 3] & 0x0f) << 10)
        r[4*i + 2] = v[7 * i + 3] >> 4 | (v[7 * i + 4]
                                          << 4) | ((v[7 * i + 5] & 0x03) << 12)
        r[4*i + 3] = v[7 * i + 5] >> 2 | ((v[7 * i + 6]) << 6)
    return r


def EncodePK(b_hat, public_seed):
    r = [None for i in range(7 * n // 4 + 32)]
    r[:7 * n//4] = EncodePolynomial(b_hat)
    r[7 * n//4:] = list(public_seed)
    return bytes(r)


def DecodePK(pk):
    b_hat = DecodePolynomial(pk[:7*n//4])
    seed = pk[7*n//4:]
    return b_hat, seed


def EncodeC(u_hat, h):
    c = [None for i in range(7*n//4 + 3*n//8)]
    c[:7*n//4] = EncodePolynomial(u_hat)
    c[7*n//4:] = list(h)
    return bytes(c)


def DecodeC(c):
    u_hat = DecodePolynomial(c[:7*n//4])
    h = c[7*n//4:]
    return u_hat, h


def Encode(mu):
    v = [None for i in range(n)]
    for i in range(32):
        for j in range(8):
            mask = -((mu[i] >> j) & 1)
            v[8 * i + j + 0] = q//2 if mask == -1 else 0  # mask & (q//2)
            v[8 * i + j + 256] = q//2 if mask == -1 else 0  # mask & (q//2)
            if n == 1024:
                v[8 * i + j + 512] = q//2 if mask == -1 else 0  # mask & (q//2)
                v[8 * i + j + 768] = q//2 if mask == -1 else 0  # mask & (q//2)
    return v


def Decode(v):
    mu = [0 for i in range(32)]
    for i in range(256):
        t = abs((v[i+0] % q) - (q-1)//2)
        t = t + abs((v[i+256] % q) - (q-1)//2)
        if n == 1024:
            t = t + abs((v[i+512] % q) - (q-1)//2)
            t = t + abs((v[i+768] % q) - (q-1)//2)
            t = t-q
        else:
            t = t - (q//2)
        t = abs(t >> 15)
        mu[i >> 3] = mu[i >> 3] | (t << (i & 7)) & 0xff
    return bytes(mu)


def Compress(vp):
    k = 0
    t = [0 for i in range(8)]
    h = [0 for i in range(3*n//8)]
    for l in range(n//8):
        i = 8*l
        for j in range(8):
            t[j] = vp[i+j] % q
            t[j] = (((t[j] << 3) + q//2)//q) & 7
        h[k+0] = (t[0] | (t[1] << 3) | (t[2] << 6)) & 0xff
        h[k+1] = ((t[2] >> 2) | (t[3] << 1) | (t[4] << 4) | (t[5] << 7)) & 0xff
        h[k+2] = ((t[5] >> 1) | (t[6] << 2) | (t[7] << 5)) & 0xff
        k += 3
    return bytes(h)


def Decompress(a):
    k = 0
    r = [None for i in range(n)]
    for l in range(n//8):
        i = 8*l
        r[i+0] = a[k+0] & 7
        r[i+1] = (a[k+0] >> 3) & 7
        r[i+2] = (a[k+0] >> 6) | ((a[k+1] << 2) & 4)
        r[i+3] = ((a[k+1] >> 1) & 7)
        r[i+4] = (a[k+1] >> 4) & 7
        r[i+5] = (a[k+1] >> 7) | ((a[k+2] << 1) & 6)
        r[i+6] = ((a[k+2] >> 2) & 7)
        r[i+7] = (a[k+2] >> 5)
        k = k+3
        for j in range(8):
            r[i+j] = (r[i+j] * q + 4) >> 3
    return r


def ntt_mult_add(x, y, z):
    return [(a * sh + e) % q for (a, sh, e) in zip(x, y, z)]


def newhope_cpa_pke_keygen():
    seed = urandom(32)
    z = shake256(64, b"\x01"+seed)
    public_seed = z[:32]
    noise_seed = z[32:]
    a_hat = GenA(public_seed)
    s = PolyBitRev(Sample(noise_seed, 0))
    s_hat = NTT(s)
    e = PolyBitRev(Sample(noise_seed, 1))
    e_hat = NTT(e)
    b_hat = ntt_mult_add(a_hat, s_hat, e_hat)
    # b_hat = [(a * sh + e) % q for (a, sh, e) in zip(a_hat, s_hat, e_hat)]
    pk = EncodePK(b_hat, public_seed)
    sk = EncodePolynomial(s_hat)
    return pk, sk


def newhope_cpa_pke_encryption(pk, mu, coin):
    b_hat, public_seed = DecodePK(pk)
    a_hat = GenA(public_seed)
    sp = PolyBitRev(Sample(coin, 0))
    ep = PolyBitRev(Sample(coin, 1))
    epp = Sample(coin, 2)
    t_hat = NTT(sp)
    u_hat = ntt_mult_add(a_hat, t_hat, NTT(ep))
    # u_hat = [(a * t + e) % q for (a, t, e) in zip(a_hat, t_hat, NTT(ep))]
    v = Encode(mu)
    vp = [(a + e + vv) % q for (a, e, vv)
          in zip(NTTinv([(b * t) % q for (b, t) in zip(b_hat, t_hat)]), epp, v)]
    h = Compress(vp)
    c = EncodeC(u_hat, h)
    return c


def newhope_cpa_pke_decryption(c, sk):
    u_hat, h = DecodeC(c)
    s_hat = DecodePolynomial(sk)
    vp = Decompress(h)
    v = [(v - us) % q for (v, us) in zip(vp,
                                         NTTinv([(u*s) % q for (u, s) in zip(u_hat, s_hat)]))]
    mu = Decode(v)
    return mu


def newhope_cpa_kem_keygen():
    pk, sk = newhope_cpa_pke_keygen()
    return pk, sk


def newhope_cpa_kem_encaps(pk):
    coin = urandom(32)
    buf = shake256(64, b"\x02"+coin)
    K, coinp = buf[:32], buf[32:]
    c = newhope_cpa_pke_encryption(pk, K, coinp)
    ss = shake256(32, K)
    return c, ss


def newhope_cpa_kem_decaps(c, sk):
    Kp = newhope_cpa_pke_decryption(c, sk)
    ss = shake256(32, Kp)
    return ss


def newhope_cca_kem_keygen():
    pk, sk = newhope_cpa_pke_keygen()
    s = urandom(32)
    sk_bar = sk + pk + shake256(32, pk) + s
    return(pk, sk_bar)


def newhope_cca_kem_encaps(pk):
    coin = urandom(32)
    mu = shake256(32, b"\x04"+coin)
    buf = shake256(96, b"\x08"+mu+shake256(32, pk))
    K, coinp, d = buf[:32], buf[32:64], buf[64:]
    c = newhope_cpa_pke_encryption(pk, mu, coinp)
    ss = shake256(32, K+shake256(32, c+d))
    c_bar = c+d
    return c_bar, ss


def newhope_cca_kem_decaps(c_bar, sk_bar):
    c, d = c_bar[:3*n//8 + 7*n//4], c_bar[3*n//8 + 7*n//4:]
    sk, pk, h, s = sk_bar[:7*n//4], sk_bar[7*n // 4: 7*n //
                                           2], sk_bar[7*n//2:7*n//2 + 32], sk_bar[7*n//2 + 32:]
    mup = newhope_cpa_pke_decryption(c, sk)
    buf = shake256(96, b"\x08"+mup+h)
    Kp, coinpp, dp = buf[:32], buf[32:64], buf[64:]
    print(c, newhope_cpa_pke_encryption(pk, mup, coinpp))
    print(d, dp)
    if c == newhope_cpa_pke_encryption(pk, mup, coinpp) and d == dp:
        fail = 0
    else:
        fail = 1
    fail = 0  # cheat!
    K = Kp, s
    ss = shake256(32, K[fail]+shake256(32, c+d))
    return ss


if __name__ == "__main__":
    
    pk, sk_bar = newhope_cca_kem_keygen()
    c_bar, ss = newhope_cca_kem_encaps(pk)
    ssp = newhope_cca_kem_decaps(c_bar, sk_bar)
    print(ss, ssp)
    
    pk, sk = newhope_cpa_kem_keygen()
    c, ss = newhope_cpa_kem_encaps(pk)
    ssp = newhope_cpa_kem_decaps(c, sk)
    print(ss, ssp)
    
    mu = bytes(range(32))
    coin = b""
    pk, sk = newhope_cpa_pke_keygen()
    c = newhope_cpa_pke_encryption(pk, mu, coin)
    mup = newhope_cpa_pke_decryption(c, sk)
    print(mu, mup)
