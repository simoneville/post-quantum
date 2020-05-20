import newHope
import unittest
import random


class TestShake(unittest.TestCase):
    def test_shake_256(self):
        x = newHope.shake256(5, b"123")
        self.assertEqual(len(x), 5)
        self.assertIsInstance(x, bytes)

    def test_shake_absorb(self):
        x = newHope.shake128Absorb(b"asdf")
        self.assertIsInstance(x, newHope.shake_128)

    def test_shake_squeeze(self):
        x = newHope.shake128Absorb(b"asdf")
        out, y = newHope.shake128Squeeze(1, x)
        self.assertIsInstance(y, newHope.shake_128)
        self.assertIsInstance(out, bytes)
        self.assertGreaterEqual(len(out), 168)


class TestRev(unittest.TestCase):
    def test_bitrev(self):
        x = 0b0001101100
        x_rev = 0b0011011000
        self.assertEqual(newHope.BitRev(x), x_rev)

    def test_bitrev_self_inverse(self):
        x = 0x167
        self.assertEqual(x, newHope.BitRev(newHope.BitRev(x)))

    def test_polybitrev(self):
        x = [1, 2, 3, 4] + [0 for _ in range(1024 - 4)]
        x_rev = [0 for _ in range(1024)]
        x_rev[0] = 1
        x_rev[512] = 2
        x_rev[256] = 3
        x_rev[768] = 4
        self.assertEqual(newHope.PolyBitRev(x), x_rev)

    def test_polybitrev_self_inverse(self):
        x = list(range(1024))
        self.assertEqual(x, newHope.PolyBitRev(newHope.PolyBitRev(x)))


class TestHW(unittest.TestCase):
    def test_zero_hw(self):
        self.assertEqual(newHope.HW(0), 0)

    def test_nonzero_hw(self):
        self.assertEqual(newHope.HW(0x654), 5)


class TestNTT(unittest.TestCase):
    def test_NTTinv_is_NTT_inv(self):
        x = list(range(1024))
        self.assertEqual(x, newHope.NTTinv(newHope.NTT(x)))

    def test_NTT_is_NTTinv_inv(self):
        x = list(range(1024))
        self.assertEqual(x, newHope.NTT(newHope.NTTinv(x)))

    def test_multiply_by_x(self):
        x = list(range(1024))
        y = [1 if _ == 1 else 0 for _ in range(1024)]
        x_hat = newHope.NTT(x)
        y_hat = newHope.NTT(y)
        z_hat = [a*b for (a, b) in zip(x_hat, y_hat)]
        z = newHope.NTTinv(z_hat)
        self.assertEqual(z, [newHope.q - 1023] + list(range(1023)))


class TestEncoding(unittest.TestCase):
    def test_poly_codec(self):
        x = list(range(1024))
        self.assertEqual(x, newHope.DecodePolynomial(
            newHope.EncodePolynomial(x)))

    def test_poly_decod(self):
        x = bytes([_ % 256 for _ in range(7 * 1024 // 4)])
        self.assertEqual(x, newHope.EncodePolynomial(
            newHope.DecodePolynomial(x)))

    def test_pk_codec(self):
        b_hat = list(range(1024))
        seed = bytes(range(32))
        self.assertEqual((b_hat, seed), newHope.DecodePK(
            newHope.EncodePK(b_hat, seed)))

    def test_msg_codec(self):
        x = bytes(range(32))
        self.assertEqual(x, newHope.Decode(newHope.Encode(x)))

    def test_c_codec(self):
        # TODO
        u_hat = list(range(1024))
        h = bytes(i % 256 for i in range(1024))
        self.assertEqual((u_hat, h), newHope.DecodeC(
            newHope.EncodeC(u_hat, h)))

    def test_comp_decomp(self):
        vp = [((newHope.q * _ // 16) + i) % newHope.q for (_, i)
              in zip(range(0, 1 << 16, 1 << 6), range(1023, -1, -1))]
        h = newHope.Compress(vp)
        self.assertLessEqual(
            max([abs(newHope.minabsmod(a - b)) for (a, b) in zip(vp, newHope.Decompress(h))]), newHope.q//4)
        v = [random.randint(0, newHope.q - 1) for _ in range(1024)]
        vp = [v[_] + random.randint(- newHope.q//16, newHope.q//16)
              for _ in range(1024)]
        h = newHope.Compress(vp)
        self.assertLessEqual(
            max([abs(newHope.minabsmod(a - b)) for (a, b) in zip(vp, newHope.Decompress(h))]), newHope.q//4)
