/**
 * 
 */
package linearalgebra;

import static org.junit.Assert.*;

import org.junit.Test;

/**
 * Decomposition test.
 * 
 * @author Jacob Malter
 *
 */
public class DecompositionTest {

	/** acceptable margin of error */
	private static final double DELTA = 0.0001d;

	// Related to LU Decomposition

	/**
	 * Test method for {@link Decomposition#lu(double[][])}.
	 */
	@Test
	public void testLu() {
		// variables required for testing
		double[][] a;
		double[][] actualL;
		double[][] actualU;
		double[][] actualA;
		double[][] expectedLu;
		double[][] actualLu;
		double[][] testLu = null;

		// first test on 2 by 2 array
		a = new double[][] { new double[] { 2, -1 }, new double[] { -4, -1 } };
		expectedLu = new double[][] { new double[] { 1, 0 },
				new double[] { -2, 1 }, new double[] { 2, -1 },
				new double[] { 0, -3 } };
		actualLu = Decomposition.lu(a);

		assertEquals(expectedLu.length, actualLu.length);
		for (int row = 0; row < expectedLu.length; row++) {
			assertEquals(expectedLu[row].length, actualLu[row].length);
			for (int col = 0; col < expectedLu[row].length; col++)
				assertEquals(expectedLu[row][col], actualLu[row][col], DELTA);
		}

		actualL = new double[a.length][a.length];
		actualU = new double[a.length][a.length];
		for (int row = 0; row < actualLu.length; row++)
			for (int col = 0; col < actualLu[row].length; col++)
				if (row < a.length)
					actualL[row][col] = actualLu[row][col];

				else
					actualU[row - a.length][col] = actualLu[row][col];
		actualA = Matrices.multiply(actualL, actualU);

		assertEquals(a.length, actualA.length);
		for (int row = 0; row < a.length; row++) {
			assertEquals(a[row].length, actualA[row].length);
			for (int col = 0; col < a[row].length; col++)
				assertEquals(a[row][col], actualA[row][col], DELTA);
		}

		// second test on 3 by 3 array
		a = new double[][] { new double[] { 2, 4, -1 },
				new double[] { -6, -15, 5 }, new double[] { 8, 22, -10 } };
		expectedLu = new double[][] { new double[] { 1, 0, 0 },
				new double[] { -3, 1, 0, }, new double[] { 4, -2, 1 },
				new double[] { 2, 4, -1 }, new double[] { 0, -3, 2 },
				new double[] { 0, 0, -2 } };
		actualLu = Decomposition.lu(a);

		assertEquals(expectedLu.length, actualLu.length);
		for (int row = 0; row < expectedLu.length; row++) {
			assertEquals(expectedLu[row].length, actualLu[row].length);
			for (int col = 0; col < expectedLu[row].length; col++)
				assertEquals(expectedLu[row][col], actualLu[row][col], DELTA);
		}

		actualL = new double[a.length][a.length];
		actualU = new double[a.length][a.length];
		for (int row = 0; row < actualLu.length; row++)
			for (int col = 0; col < actualLu[row].length; col++)
				if (row < a.length)
					actualL[row][col] = actualLu[row][col];

				else
					actualU[row - a.length][col] = actualLu[row][col];
		actualA = Matrices.multiply(actualL, actualU);

		assertEquals(a.length, actualA.length);
		for (int row = 0; row < a.length; row++) {
			assertEquals(a[row].length, actualA[row].length);
			for (int col = 0; col < a[row].length; col++)
				assertEquals(a[row][col], actualA[row][col], DELTA);
		}

		// third test on 4 by 4 array
		a = new double[][] { new double[] { 2, -1, 2, -4 },
				new double[] { 4, -5, 5, -12 }, new double[] { 2, 8, -5, 5 },
				new double[] { -8, 16, 4, 40 } };
		expectedLu = new double[][] { new double[] { 1, 0, 0, 0 },
				new double[] { 2, 1, 0, 0 }, new double[] { 1, -3, 1, 0 },
				new double[] { -4, -4, -4, 1 }, new double[] { 2, -1, 2, -4 },
				new double[] { 0, -3, 1, -4 }, new double[] { 0, 0, -4, -3 },
				new double[] { 0, 0, 0, -4 } };
		actualLu = Decomposition.lu(a);

		assertEquals(expectedLu.length, actualLu.length);
		for (int row = 0; row < expectedLu.length; row++) {
			assertEquals(expectedLu[row].length, actualLu[row].length);
			for (int col = 0; col < expectedLu[row].length; col++)
				assertEquals(expectedLu[row][col], actualLu[row][col], DELTA);
		}

		actualL = new double[a.length][a.length];
		actualU = new double[a.length][a.length];
		for (int row = 0; row < actualLu.length; row++)
			for (int col = 0; col < actualLu[row].length; col++)
				if (row < a.length)
					actualL[row][col] = actualLu[row][col];

				else
					actualU[row - a.length][col] = actualLu[row][col];
		actualA = Matrices.multiply(actualL, actualU);

		assertEquals(a.length, actualA.length);
		for (int row = 0; row < a.length; row++) {
			assertEquals(a[row].length, actualA[row].length);
			for (int col = 0; col < a[row].length; col++)
				assertEquals(a[row][col], actualA[row][col], DELTA);
		}

		try {
			// matrix a is null
			testLu = Decomposition.lu(null);
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testLu);
		}

		try {
			// matrix a does contain null column
			testLu = Decomposition.lu(new double[][] { null });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testLu);
		}

		try {
			// matrix a must does not have as many columns as rows
			testLu = Decomposition.lu(new double[][] { new double[] { 0, 0 },
					new double[] { 0, 0 }, new double[] { 0, 0 } });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testLu);
		}
	}

	// Related to Cholesky Decomposition

	/**
	 * Test method for {@link Decomposition#cholesky(double[][])}.
	 */
	@Test
	public void testCholesky() {
		// variables required for testing
		double[][] a;
		double[][] actualA;
		double[][] expectedL;
		double[][] actualL;
		double[][] testL = null;

		// first test on 2 by 2 array
		a = new double[][] { new double[] { 16, 8 }, new double[] { 8, 20 } };
		expectedL = new double[][] { new double[] { 4, 0 },
				new double[] { 2, 4 } };
		actualL = Decomposition.cholesky(a);

		assertEquals(expectedL.length, actualL.length);
		for (int row = 0; row < expectedL.length; row++) {
			assertEquals(expectedL[row].length, actualL[row].length);
			for (int col = 0; col < expectedL[row].length; col++)
				assertEquals(expectedL[row][col], actualL[row][col], DELTA);
		}

		actualA = Matrices.multiply(actualL, Matrices.transpose(actualL));

		assertEquals(a.length, actualA.length);
		for (int row = 0; row < a.length; row++) {
			assertEquals(a[row].length, actualA[row].length);
			for (int col = 0; col < a[row].length; col++)
				assertEquals(a[row][col], actualA[row][col], DELTA);
		}

		// second test on 3 by 3 array
		a = new double[][] { new double[] { 4, 4, -6 },
				new double[] { 4, 20, 10 }, new double[] { -6, 10, 29 } };
		expectedL = new double[][] { new double[] { 2, 0, 0 },
				new double[] { 2, 4, 0 }, new double[] { -3, 4, 2 } };
		actualL = Decomposition.cholesky(a);

		assertEquals(expectedL.length, actualL.length);
		for (int row = 0; row < expectedL.length; row++) {
			assertEquals(expectedL[row].length, actualL[row].length);
			for (int col = 0; col < expectedL[row].length; col++)
				assertEquals(expectedL[row][col], actualL[row][col], DELTA);
		}

		actualA = Matrices.multiply(actualL, Matrices.transpose(actualL));

		assertEquals(a.length, actualA.length);
		for (int row = 0; row < a.length; row++) {
			assertEquals(a[row].length, actualA[row].length);
			for (int col = 0; col < a[row].length; col++)
				assertEquals(a[row][col], actualA[row][col], DELTA);
		}

		// third test on 4 by 4 array
		a = new double[][] { new double[] { 4, 8, -6, 6 },
				new double[] { 8, 17, -14, 14 },
				new double[] { -6, -14, 17, -5 },
				new double[] { 6, 14, -5, 38 } };
		expectedL = new double[][] { new double[] { 2, 0, 0, 0 },
				new double[] { 4, 1, 0, 0 }, new double[] { -3, -2, 2, 0 },
				new double[] { 3, 2, 4, 3 } };
		actualL = Decomposition.cholesky(a);

		assertEquals(expectedL.length, actualL.length);
		for (int row = 0; row < expectedL.length; row++) {
			assertEquals(expectedL[row].length, actualL[row].length);
			for (int col = 0; col < expectedL[row].length; col++)
				assertEquals(expectedL[row][col], actualL[row][col], DELTA);
		}

		actualA = Matrices.multiply(actualL, Matrices.transpose(actualL));

		assertEquals(a.length, actualA.length);
		for (int row = 0; row < a.length; row++) {
			assertEquals(a[row].length, actualA[row].length);
			for (int col = 0; col < a[row].length; col++)
				assertEquals(a[row][col], actualA[row][col], DELTA);
		}

		try {
			// matrix a is null
			testL = Decomposition.cholesky(null);
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testL);
		}

		try {
			// matrix a does contain null column
			testL = Decomposition.cholesky(new double[][] { null });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testL);
		}

		try {
			// matrix a must does not have as many columns as rows
			testL = Decomposition.cholesky(new double[][] {
					new double[] { 0, 0 }, new double[] { 0, 0 },
					new double[] { 0, 0 } });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testL);
		}

		try {
			// matrix a must does not have symmetry across diagonal
			testL = Decomposition.cholesky(new double[][] {
					new double[] { 1, 2 }, new double[] { 3, 4 } });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testL);
		}
	}

	// Related to QR Decomposition

	/**
	 * Test method for {@link Decomposition#qr(double[][])}.
	 */
	@Test
	public void testQr() {
		// variables required for testing
		double[][] a;
		double[][] actualQ;
		double[][] actualR;
		double[][] actualA;
		double[][] actualI;
		double[][] expectedQr;
		double[][] actualQr;
		double[][] testQr = null;

		// first test on 2 by 2 array
		a = new double[][] { new double[] { 14, -6 }, new double[] { 14, -9 },
				new double[] { 7, -6 } };
		expectedQr = new double[][] { new double[] { 2d / 3, 2d / 3 },
				new double[] { 2d / 3, -1d / 3 },
				new double[] { 1d / 3, -2d / 3 }, new double[] { 21, -12 },
				new double[] { 0, 3 } };
		actualQr = Decomposition.qr(a);

		assertEquals(expectedQr.length, actualQr.length);
		for (int row = 0; row < expectedQr.length; row++) {
			assertEquals(expectedQr[row].length, actualQr[row].length);
			for (int col = 0; col < expectedQr[row].length; col++)
				assertEquals(expectedQr[row][col], actualQr[row][col], DELTA);
		}

		actualQ = new double[a.length][a[0].length];
		actualR = new double[a[0].length][a[0].length];
		for (int row = 0; row < actualQr.length; row++)
			for (int col = 0; col < actualQr[row].length; col++)
				if (row < a.length)
					actualQ[row][col] = actualQr[row][col];

				else
					actualR[row - a.length][col] = actualQr[row][col];
		actualA = Matrices.multiply(actualQ, actualR);

		assertEquals(a.length, actualA.length);
		for (int row = 0; row < a.length; row++) {
			assertEquals(a[row].length, actualA[row].length);
			for (int col = 0; col < a[row].length; col++)
				assertEquals(a[row][col], actualA[row][col], DELTA);
		}

		actualI = Matrices.multiply(Matrices.transpose(actualQ), actualQ);
		for (int row = 0; row < actualI.length; row++)
			for (int col = 0; col < actualI[row].length; col++)
				if (row == col)
					assertEquals(1, actualI[row][col], DELTA);

				else
					assertEquals(0, actualI[row][col], DELTA);

		// second test on 4 by 3 array
		a = new double[][] { new double[] { -3, -8, -9 },
				new double[] { -3, -10, 15 }, new double[] { -3, -8, 7 },
				new double[] { -3, -10, -1 } };
		expectedQr = new double[][] { new double[] { -0.5, 0.5, -0.5 },
				new double[] { -0.5, -0.5, 0.5 },
				new double[] { -0.5, 0.5, 0.5 },
				new double[] { -0.5, -0.5, -0.5 }, new double[] { 6, 18, -6 },
				new double[] { 0, 2, -8 }, new double[] { 0, 0, 16 } };
		actualQr = Decomposition.qr(a);

		assertEquals(expectedQr.length, actualQr.length);
		for (int row = 0; row < expectedQr.length; row++) {
			assertEquals(expectedQr[row].length, actualQr[row].length);
			for (int col = 0; col < expectedQr[row].length; col++)
				assertEquals(expectedQr[row][col], actualQr[row][col], DELTA);
		}

		actualQ = new double[a.length][a[0].length];
		actualR = new double[a[0].length][a[0].length];
		for (int row = 0; row < actualQr.length; row++)
			for (int col = 0; col < actualQr[row].length; col++)
				if (row < a.length)
					actualQ[row][col] = actualQr[row][col];

				else
					actualR[row - a.length][col] = actualQr[row][col];
		actualA = Matrices.multiply(actualQ, actualR);

		assertEquals(a.length, actualA.length);
		for (int row = 0; row < a.length; row++) {
			assertEquals(a[row].length, actualA[row].length);
			for (int col = 0; col < a[row].length; col++)
				assertEquals(a[row][col], actualA[row][col], DELTA);
		}

		actualI = Matrices.multiply(Matrices.transpose(actualQ), actualQ);
		for (int row = 0; row < actualI.length; row++)
			for (int col = 0; col < actualI[row].length; col++)
				if (row == col)
					assertEquals(1, actualI[row][col], DELTA);

				else
					assertEquals(0, actualI[row][col], DELTA);

		try {
			// matrix a is null
			testQr = Decomposition.qr(null);
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testQr);
		}

		try {
			// matrix a does contain null column
			testQr = Decomposition.qr(new double[][] { null });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testQr);
		}

		try {
			// matrix a contains as fewer rows than columns
			testQr = Decomposition.qr(new double[][] {
					new double[] { 0, 0, 0 }, new double[] { 0, 0, 0 } });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testQr);
		}

		try {
			// matrix a row width does not remain constant
			testQr = Decomposition.qr(new double[][] { new double[] { 0, 0 },
					new double[] { 0 } });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testQr);
		}
	}

	// Related to Eigen Decomposition

	/**
	 * Test method for {@link Decomposition#eigen(double[][])}.
	 */
	@Test
	public void testEigen() {
		// variables required for testing
		double[][] a;
		double[][] actualL;
		double[][] actualV;
		double[][] actualAV;
		double[][] actualLV;
		double[][] expectedEig;
		double[][] actualEig;
		double[][] testEig = null;

		// first test on 1 by 1 array
		a = new double[][] { new double[] { 7 } };
		expectedEig = new double[][] { new double[] { 7 }, new double[] { 1 } };
		actualEig = Decomposition.eigen(a);

		assertEquals(expectedEig.length, actualEig.length);
		for (int row = 0; row < expectedEig.length; row++) {
			assertEquals(expectedEig[row].length, actualEig[row].length);
			for (int col = 0; col < expectedEig[row].length; col++)
				assertEquals(expectedEig[row][col], actualEig[row][col], DELTA);
		}

		actualL = new double[a.length][a[0].length];
		actualV = new double[a[0].length][a[0].length];
		for (int row = 0; row < actualEig.length; row++)
			for (int col = 0; col < actualEig[row].length; col++)
				if (row < a.length)
					actualL[row][col] = actualEig[row][col];

				else
					actualV[row - a.length][col] = actualEig[row][col];
		actualAV = Matrices.multiply(a, actualV);
		actualLV = Matrices.multiply(actualV, actualL);

		assertEquals(actualAV.length, actualLV.length);
		for (int row = 0; row < actualAV.length; row++) {
			assertEquals(actualAV[row].length, actualLV[row].length);
			for (int col = 0; col < actualAV[row].length; col++)
				assertEquals(actualAV[row][col], actualLV[row][col], DELTA);
		}

		// second test on 2 by 2 array
		a = new double[][] { new double[] { -1, 6 }, new double[] { -2, 6 } };
		expectedEig = new double[][] { new double[] { 3, 0 },
				new double[] { 0, 2 },
				new double[] { 3 / Math.sqrt(13d), 2 / Math.sqrt(5d) },
				new double[] { 2 / Math.sqrt(13d), 1 / Math.sqrt(5d) } };
		actualEig = Decomposition.eigen(a);

		assertEquals(expectedEig.length, actualEig.length);
		for (int row = 0; row < expectedEig.length; row++) {
			assertEquals(expectedEig[row].length, actualEig[row].length);
			for (int col = 0; col < expectedEig[row].length; col++) {
				assertEquals(expectedEig[row][col], actualEig[row][col], DELTA);
			}
		}

		actualL = new double[a.length][a[0].length];
		actualV = new double[a[0].length][a[0].length];
		for (int row = 0; row < actualEig.length; row++)
			for (int col = 0; col < actualEig[row].length; col++)
				if (row < a.length)
					actualL[row][col] = actualEig[row][col];

				else
					actualV[row - a.length][col] = actualEig[row][col];
		actualAV = Matrices.multiply(a, actualV);
		actualLV = Matrices.multiply(actualV, actualL);

		assertEquals(actualAV.length, actualLV.length);
		for (int row = 0; row < actualAV.length; row++) {
			assertEquals(actualAV[row].length, actualLV[row].length);
			for (int col = 0; col < actualAV[row].length; col++)
				assertEquals(actualAV[row][col], actualLV[row][col], DELTA);
		}

		// third test on 2 by 2 array
		a = new double[][] { new double[] { -2, -2 }, new double[] { 2, 3 } };
		expectedEig = new double[][] { new double[] { 2, 0 },
				new double[] { 0, -1 },
				new double[] { -1 / Math.sqrt(5), -2 / Math.sqrt(5) },
				new double[] { 2 / Math.sqrt(5), 1 / Math.sqrt(5) } };
		actualEig = Decomposition.eigen(a);

		assertEquals(expectedEig.length, actualEig.length);
		for (int row = 0; row < expectedEig.length; row++) {
			assertEquals(expectedEig[row].length, actualEig[row].length);
			for (int col = 0; col < expectedEig[row].length; col++)
				assertEquals(expectedEig[row][col], actualEig[row][col], DELTA);
		}

		actualL = new double[a.length][a[0].length];
		actualV = new double[a[0].length][a[0].length];
		for (int row = 0; row < actualEig.length; row++)
			for (int col = 0; col < actualEig[row].length; col++)
				if (row < a.length)
					actualL[row][col] = actualEig[row][col];

				else
					actualV[row - a.length][col] = actualEig[row][col];
		actualAV = Matrices.multiply(a, actualV);
		actualLV = Matrices.multiply(actualV, actualL);

		assertEquals(actualAV.length, actualLV.length);
		for (int row = 0; row < actualAV.length; row++) {
			assertEquals(actualAV[row].length, actualLV[row].length);
			for (int col = 0; col < actualAV[row].length; col++)
				assertEquals(actualAV[row][col], actualLV[row][col], DELTA);
		}

		// fourth test on 3 by 3 array
		a = new double[][] { new double[] { 2, 3, -3 },
				new double[] { -2, -1, 2 }, new double[] { 2, 4, -3 } };
		expectedEig = new double[][] {
				new double[] { 1, 0, 0 },
				new double[] { 0, -1, 0 },
				new double[] { 0, 0, -2 },
				new double[] { 0, 1 / Math.sqrt(2), 9 / Math.sqrt(185) },
				new double[] { 1 / Math.sqrt(2), 0, -2 / Math.sqrt(185) },
				new double[] { 1 / Math.sqrt(2), 1 / Math.sqrt(2),
						10 / Math.sqrt(185) } };
		actualEig = Decomposition.eigen(a);

		assertEquals(expectedEig.length, actualEig.length);
		for (int row = 0; row < expectedEig.length; row++) {
			assertEquals(expectedEig[row].length, actualEig[row].length);
			for (int col = 0; col < expectedEig[row].length; col++)
				assertEquals(expectedEig[row][col], actualEig[row][col], DELTA);
		}

		actualL = new double[a.length][a[0].length];
		actualV = new double[a[0].length][a[0].length];
		for (int row = 0; row < actualEig.length; row++)
			for (int col = 0; col < actualEig[row].length; col++)
				if (row < a.length)
					actualL[row][col] = actualEig[row][col];

				else
					actualV[row - a.length][col] = actualEig[row][col];
		actualAV = Matrices.multiply(a, actualV);
		actualLV = Matrices.multiply(actualV, actualL);

		assertEquals(actualAV.length, actualLV.length);
		for (int row = 0; row < actualAV.length; row++) {
			assertEquals(actualAV[row].length, actualLV[row].length);
			for (int col = 0; col < actualAV[row].length; col++)
				assertEquals(actualAV[row][col], actualLV[row][col], DELTA);
		}

		// fifth test on 3 by 3 array
		a = new double[][] { new double[] { 3, 2, 2 },
				new double[] { -1, -4, 2 }, new double[] { -2, -2, -1 } };
		expectedEig = new double[][] {
				new double[] { 1, 0, 0 },
				new double[] { 0, -1, 0 },
				new double[] { 0, 0, -2 },
				new double[] { -7 / Math.sqrt(74), -1 / Math.sqrt(3),
						-2 / Math.sqrt(17) },
				new double[] { 3 / Math.sqrt(74), 1 / Math.sqrt(3),
						3 / Math.sqrt(17) },
				new double[] { 4 / Math.sqrt(74), 1 / Math.sqrt(3),
						2 / Math.sqrt(17) } };
		actualEig = Decomposition.eigen(a);

		assertEquals(expectedEig.length, actualEig.length);
		for (int row = 0; row < expectedEig.length; row++) {
			assertEquals(expectedEig[row].length, actualEig[row].length);
			for (int col = 0; col < expectedEig[row].length; col++)
				assertEquals(expectedEig[row][col], actualEig[row][col], DELTA);
		}

		actualL = new double[a.length][a[0].length];
		actualV = new double[a[0].length][a[0].length];
		for (int row = 0; row < actualEig.length; row++)
			for (int col = 0; col < actualEig[row].length; col++)
				if (row < a.length)
					actualL[row][col] = actualEig[row][col];

				else
					actualV[row - a.length][col] = actualEig[row][col];
		actualAV = Matrices.multiply(a, actualV);
		actualLV = Matrices.multiply(actualV, actualL);

		assertEquals(actualAV.length, actualLV.length);
		for (int row = 0; row < actualAV.length; row++) {
			assertEquals(actualAV[row].length, actualLV[row].length);
			for (int col = 0; col < actualAV[row].length; col++)
				assertEquals(actualAV[row][col], actualLV[row][col], DELTA);
		}

		try {
			// matrix a is null
			testEig = Decomposition.eigen(null);
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testEig);
		}

		try {
			// matrix a does contain null column
			testEig = Decomposition.eigen(new double[][] { null });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testEig);
		}

		try {
			// matrix a must does not have as many columns as rows
			testEig = Decomposition.eigen(new double[][] {
					new double[] { 0, 0 }, new double[] { 0, 0 },
					new double[] { 0, 0 } });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testEig);
		}

		try {
			// matrix a is too large
			testEig = Decomposition.eigen(new double[][] {
					new double[] { 0, 0, 0, 0 }, new double[] { 0, 0, 0, 0 },
					new double[] { 0, 0, 0, 0 }, new double[] { 0, 0, 0, 0 } });
			fail();
		} catch (UnsupportedOperationException ex) {
			assertNull(testEig);
		}
	}

	// Related to Singular Value Decomposition

	/**
	 * Test method for {@link Decomposition#singularValue(double[][])}.
	 */
	@Test
	public void testSingularValue() {
		// variables required for testing
		double[][] a;
		double[][] actualA;
		double[][] actualU;
		double[][] actualS;
		double[][] actualV;
		double[][] actualI;
		double[][] expectedSing;
		double[][] actualSing;
		double[][] testSing = null;

		// first test on 2 by 2 array
		a = new double[][] { new double[] { 8, 6 }, new double[] { 3, -4 } };
		expectedSing = new double[][] { new double[] { 1, 0 },
				new double[] { 0, -1 }, new double[] { 10, 0 },
				new double[] { 0, 5 }, new double[] { 0.8, -0.6 },
				new double[] { 0.6, 0.8 } };
		actualSing = Decomposition.singularValue(a);

		assertEquals(expectedSing.length, actualSing.length);
		for (int row = 0; row < expectedSing.length; row++) {
			assertEquals(expectedSing[row].length, actualSing[row].length);
			for (int col = 0; col < expectedSing[row].length; col++)
				assertEquals(expectedSing[row][col], actualSing[row][col],
						DELTA);
		}

		actualU = new double[a.length][a[0].length];
		actualS = new double[a[0].length][a[0].length];
		actualV = new double[a[0].length][a[0].length];
		for (int row = 0; row < actualSing.length; row++)
			for (int col = 0; col < actualSing[row].length; col++)
				if (row < a.length)
					actualU[row][col] = actualSing[row][col];

				else if (row < a.length + a[0].length)
					actualS[row - a.length][col] = actualSing[row][col];

				else
					actualV[row - a.length - a[0].length][col] = actualSing[row][col];

		actualA = Matrices.multiply(Matrices.multiply(actualU, actualS),
				Matrices.transpose(actualV));

		assertEquals(a.length, actualA.length);
		for (int row = 0; row < actualA.length; row++) {
			assertEquals(a[row].length, actualA[row].length);
			for (int col = 0; col < actualA[row].length; col++)
				assertEquals(a[row][col], actualA[row][col], DELTA);
		}

		actualI = Matrices.multiply(Matrices.transpose(actualU), actualU);

		for (int row = 0; row < actualI.length; row++)
			for (int col = 0; col < actualI[row].length; col++)
				if (row == col)
					assertEquals(1, actualI[row][col], DELTA);

				else
					assertEquals(0, actualI[row][col], DELTA);

		actualI = Matrices.multiply(Matrices.transpose(actualV), actualV);

		for (int row = 0; row < actualI.length; row++)
			for (int col = 0; col < actualI[row].length; col++)
				if (row == col)
					assertEquals(1, actualI[row][col], DELTA);

				else
					assertEquals(0, actualI[row][col], DELTA);

		for (int row = 0; row < actualS.length; row++)
			for (int col = 0; col < actualS[row].length; col++)
				if (row == col)
					assertNotEquals(0, actualS[row][col], DELTA);

				else
					assertEquals(0, actualS[row][col], DELTA);

		// second test on 2 by 3 array
		a = new double[][] { new double[] { 18, 26 }, new double[] { 18, 1 },
				new double[] { -27, -14 } };
		expectedSing = new double[][] { new double[] { 2d / 3, 2d / 3 },
				new double[] { 1d / 3, -2d / 3 },
				new double[] { -2d / 3, 1d / 3 }, new double[] { 45, 0 },
				new double[] { 0, 15 }, new double[] { 0.8, -0.6 },
				new double[] { 0.6, 0.8 } };
		actualSing = Decomposition.singularValue(a);

		assertEquals(expectedSing.length, actualSing.length);
		for (int row = 0; row < expectedSing.length; row++) {
			assertEquals(expectedSing[row].length, actualSing[row].length);
			for (int col = 0; col < expectedSing[row].length; col++)
				assertEquals(expectedSing[row][col], actualSing[row][col],
						DELTA);
		}

		actualU = new double[a.length][a[0].length];
		actualS = new double[a[0].length][a[0].length];
		actualV = new double[a[0].length][a[0].length];
		for (int row = 0; row < actualSing.length; row++)
			for (int col = 0; col < actualSing[row].length; col++)
				if (row < a.length)
					actualU[row][col] = actualSing[row][col];

				else if (row < a.length + a[0].length)
					actualS[row - a.length][col] = actualSing[row][col];

				else
					actualV[row - a.length - a[0].length][col] = actualSing[row][col];

		actualA = Matrices.multiply(Matrices.multiply(actualU, actualS),
				Matrices.transpose(actualV));

		assertEquals(a.length, actualA.length);
		for (int row = 0; row < actualA.length; row++) {
			assertEquals(a[row].length, actualA[row].length);
			for (int col = 0; col < actualA[row].length; col++)
				assertEquals(a[row][col], actualA[row][col], DELTA);
		}

		actualI = Matrices.multiply(Matrices.transpose(actualU), actualU);

		for (int row = 0; row < actualI.length; row++)
			for (int col = 0; col < actualI[row].length; col++)
				if (row == col)
					assertEquals(1, actualI[row][col], DELTA);

				else
					assertEquals(0, actualI[row][col], DELTA);

		actualI = Matrices.multiply(Matrices.transpose(actualV), actualV);

		for (int row = 0; row < actualI.length; row++)
			for (int col = 0; col < actualI[row].length; col++)
				if (row == col)
					assertEquals(1, actualI[row][col], DELTA);

				else
					assertEquals(0, actualI[row][col], DELTA);

		for (int row = 0; row < actualS.length; row++)
			for (int col = 0; col < actualS[row].length; col++)
				if (row == col)
					assertNotEquals(0, actualS[row][col], DELTA);

				else
					assertEquals(0, actualS[row][col], DELTA);

		// third test on 3 by 3 array
		a = new double[][] { new double[] { 6, 6, 3 },
				new double[] { -1, 2, -2 }, new double[] { 4, -2, -4 } };
		expectedSing = new double[][] { new double[] { 1, 0, 0 },
				new double[] { 0, 0, -1 }, new double[] { 0, -1, 0 },
				new double[] { 9, 0, 0 }, new double[] { 0, 6, 0 },
				new double[] { 0, 0, 3 },
				new double[] { 2d / 3, -2d / 3, 1d / 3 },
				new double[] { 2d / 3, 1d / 3, -2d / 3 },
				new double[] { 1d / 3, 2d / 3, 2d / 3 } };
		actualSing = Decomposition.singularValue(a);

		assertEquals(expectedSing.length, actualSing.length);
		for (int row = 0; row < expectedSing.length; row++) {
			assertEquals(expectedSing[row].length, actualSing[row].length);
			for (int col = 0; col < expectedSing[row].length; col++)
				assertEquals(expectedSing[row][col], actualSing[row][col],
						DELTA);
		}

		actualU = new double[a.length][a[0].length];
		actualS = new double[a[0].length][a[0].length];
		actualV = new double[a[0].length][a[0].length];
		for (int row = 0; row < actualSing.length; row++)
			for (int col = 0; col < actualSing[row].length; col++)
				if (row < a.length)
					actualU[row][col] = actualSing[row][col];

				else if (row < a.length + a[0].length)
					actualS[row - a.length][col] = actualSing[row][col];

				else
					actualV[row - a.length - a[0].length][col] = actualSing[row][col];

		actualA = Matrices.multiply(Matrices.multiply(actualU, actualS),
				Matrices.transpose(actualV));

		assertEquals(a.length, actualA.length);
		for (int row = 0; row < actualA.length; row++) {
			assertEquals(a[row].length, actualA[row].length);
			for (int col = 0; col < actualA[row].length; col++)
				assertEquals(a[row][col], actualA[row][col], DELTA);
		}

		actualI = Matrices.multiply(Matrices.transpose(actualU), actualU);

		for (int row = 0; row < actualI.length; row++)
			for (int col = 0; col < actualI[row].length; col++)
				if (row == col)
					assertEquals(1, actualI[row][col], DELTA);

				else
					assertEquals(0, actualI[row][col], DELTA);

		actualI = Matrices.multiply(Matrices.transpose(actualV), actualV);

		for (int row = 0; row < actualI.length; row++)
			for (int col = 0; col < actualI[row].length; col++)
				if (row == col)
					assertEquals(1, actualI[row][col], DELTA);

				else
					assertEquals(0, actualI[row][col], DELTA);

		for (int row = 0; row < actualS.length; row++)
			for (int col = 0; col < actualS[row].length; col++)
				if (row == col)
					assertNotEquals(0, actualS[row][col], DELTA);

				else
					assertEquals(0, actualS[row][col], DELTA);

		// fourth test on 4 by 2 array
		a = new double[][] { new double[] { 18, 1 }, new double[] { 18, 1 },
				new double[] { 6, 17 }, new double[] { -6, -17 } };
		expectedSing = new double[][] { new double[] { 0.5, -0.5 },
				new double[] { 0.5, -0.5 }, new double[] { 0.5, 0.5 },
				new double[] { -0.5, -0.5 }, new double[] { 30, 0 },
				new double[] { 0, 20 }, new double[] { 0.8, -0.6 },
				new double[] { 0.6, 0.8 } };
		actualSing = Decomposition.singularValue(a);

		assertEquals(expectedSing.length, actualSing.length);
		for (int row = 0; row < expectedSing.length; row++) {
			assertEquals(expectedSing[row].length, actualSing[row].length);
			for (int col = 0; col < expectedSing[row].length; col++)
				assertEquals(expectedSing[row][col], actualSing[row][col],
						DELTA);
		}

		actualU = new double[a.length][a[0].length];
		actualS = new double[a[0].length][a[0].length];
		actualV = new double[a[0].length][a[0].length];
		for (int row = 0; row < actualSing.length; row++)
			for (int col = 0; col < actualSing[row].length; col++)
				if (row < a.length)
					actualU[row][col] = actualSing[row][col];

				else if (row < a.length + a[0].length)
					actualS[row - a.length][col] = actualSing[row][col];

				else
					actualV[row - a.length - a[0].length][col] = actualSing[row][col];

		actualA = Matrices.multiply(Matrices.multiply(actualU, actualS),
				Matrices.transpose(actualV));

		assertEquals(a.length, actualA.length);
		for (int row = 0; row < actualA.length; row++) {
			assertEquals(a[row].length, actualA[row].length);
			for (int col = 0; col < actualA[row].length; col++)
				assertEquals(a[row][col], actualA[row][col], DELTA);
		}

		actualI = Matrices.multiply(Matrices.transpose(actualU), actualU);

		for (int row = 0; row < actualI.length; row++)
			for (int col = 0; col < actualI[row].length; col++)
				if (row == col)
					assertEquals(1, actualI[row][col], DELTA);

				else
					assertEquals(0, actualI[row][col], DELTA);

		actualI = Matrices.multiply(Matrices.transpose(actualV), actualV);

		for (int row = 0; row < actualI.length; row++)
			for (int col = 0; col < actualI[row].length; col++)
				if (row == col)
					assertEquals(1, actualI[row][col], DELTA);

				else
					assertEquals(0, actualI[row][col], DELTA);

		for (int row = 0; row < actualS.length; row++)
			for (int col = 0; col < actualS[row].length; col++)
				if (row == col)
					assertNotEquals(0, actualS[row][col], DELTA);

				else
					assertEquals(0, actualS[row][col], DELTA);

		// fifth test on 4 by 3 array
		a = new double[][] { new double[] { 11, 8, 2 },
				new double[] { 13, 4, -2 }, new double[] { 3, 12, -6 },
				new double[] { -5, -8, 10 } };
		expectedSing = new double[][] { new double[] { -0.5, 0.5, 0.5 },
				new double[] { -0.5, 0.5, -0.5 },
				new double[] { -0.5, -0.5, 0.5 },
				new double[] { 0.5, 0.5, 0.5 }, new double[] { 24, 0, 0 },
				new double[] { 0, 12, 0 }, new double[] { 0, 0, 6 },
				new double[] { -2d / 3, 2d / 3, -1d / 3 },
				new double[] { -2d / 3, -1d / 3, 2d / 3 },
				new double[] { 1d / 3, 2d / 3, 2d / 3 } };
		actualSing = Decomposition.singularValue(a);

		assertEquals(expectedSing.length, actualSing.length);
		for (int row = 0; row < expectedSing.length; row++) {
			assertEquals(expectedSing[row].length, actualSing[row].length);
			for (int col = 0; col < expectedSing[row].length; col++)
				assertEquals(expectedSing[row][col], actualSing[row][col],
						DELTA);
		}

		actualU = new double[a.length][a[0].length];
		actualS = new double[a[0].length][a[0].length];
		actualV = new double[a[0].length][a[0].length];
		for (int row = 0; row < actualSing.length; row++)
			for (int col = 0; col < actualSing[row].length; col++)
				if (row < a.length)
					actualU[row][col] = actualSing[row][col];

				else if (row < a.length + a[0].length)
					actualS[row - a.length][col] = actualSing[row][col];

				else
					actualV[row - a.length - a[0].length][col] = actualSing[row][col];

		actualA = Matrices.multiply(Matrices.multiply(actualU, actualS),
				Matrices.transpose(actualV));

		assertEquals(a.length, actualA.length);
		for (int row = 0; row < actualA.length; row++) {
			assertEquals(a[row].length, actualA[row].length);
			for (int col = 0; col < actualA[row].length; col++)
				assertEquals(a[row][col], actualA[row][col], DELTA);
		}

		actualI = Matrices.multiply(Matrices.transpose(actualU), actualU);

		for (int row = 0; row < actualI.length; row++)
			for (int col = 0; col < actualI[row].length; col++)
				if (row == col)
					assertEquals(1, actualI[row][col], DELTA);

				else
					assertEquals(0, actualI[row][col], DELTA);

		actualI = Matrices.multiply(Matrices.transpose(actualV), actualV);

		for (int row = 0; row < actualI.length; row++)
			for (int col = 0; col < actualI[row].length; col++)
				if (row == col)
					assertEquals(1, actualI[row][col], DELTA);

				else
					assertEquals(0, actualI[row][col], DELTA);

		for (int row = 0; row < actualS.length; row++)
			for (int col = 0; col < actualS[row].length; col++)
				if (row == col)
					assertNotEquals(0, actualS[row][col], DELTA);

				else
					assertEquals(0, actualS[row][col], DELTA);

		try {
			// matrix a is null
			testSing = Decomposition.singularValue(null);
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testSing);
		}

		try {
			// matrix a does contain null column
			testSing = Decomposition.singularValue(new double[][] { null });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testSing);
		}

		try {
			// matrix a contains as fewer rows than columns
			testSing = Decomposition.singularValue(new double[][] {
					new double[] { 0, 0, 0 }, new double[] { 0, 0, 0 } });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testSing);
		}

		try {
			// matrix a row width does not remain constant
			testSing = Decomposition.singularValue(new double[][] {
					new double[] { 0, 0 }, new double[] { 0 } });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testSing);
		}
	}

}