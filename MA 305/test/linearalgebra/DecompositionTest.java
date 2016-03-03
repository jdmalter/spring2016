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

	/**
	 * Test method for {@link Decomposition#lu(double[][])}.
	 */
	@Test
	public void testLu() {
		// variables required for testing
		double[][] a;
		double[][] actualA;
		double[][] expectedL;
		double[][] actualL;
		double[][] expectedU;
		double[][] actualU;
		double[][][] actualLU;
		double[][][] testLU = null;

		// first test on 2 by 2 array
		a = new double[][] { new double[] { 2, -1 }, new double[] { -4, -1 } };
		expectedL = new double[][] { new double[] { 1, 0 },
				new double[] { -2, 1 } };
		expectedU = new double[][] { new double[] { 2, -1 },
				new double[] { 0, -3 } };
		actualLU = Decomposition.lu(a);
		actualL = actualLU[0];
		actualU = actualLU[1];

		assertEquals(expectedL.length, actualL.length);
		for (int row = 0; row < expectedL.length; row++) {
			assertEquals(expectedL[row].length, actualL[row].length);
			for (int col = 0; col < expectedL[row].length; col++)
				assertEquals(expectedL[row][col], actualL[row][col], DELTA);
		}

		assertEquals(expectedU.length, actualU.length);
		for (int row = 0; row < expectedU.length; row++) {
			assertEquals(expectedU[row].length, actualU[row].length);
			for (int col = 0; col < expectedU[row].length; col++)
				assertEquals(expectedU[row][col], actualU[row][col], DELTA);
		}

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
		expectedL = new double[][] { new double[] { 1, 0, 0 },
				new double[] { -3, 1, 0, }, new double[] { 4, -2, 1 } };
		expectedU = new double[][] { new double[] { 2, 4, -1 },
				new double[] { 0, -3, 2 }, new double[] { 0, 0, -2 } };
		actualLU = Decomposition.lu(a);
		actualL = actualLU[0];
		actualU = actualLU[1];

		assertEquals(expectedL.length, actualL.length);
		for (int row = 0; row < expectedL.length; row++) {
			assertEquals(expectedL[row].length, actualL[row].length);
			for (int col = 0; col < expectedL[row].length; col++)
				assertEquals(expectedL[row][col], actualL[row][col], DELTA);
		}

		assertEquals(expectedU.length, actualU.length);
		for (int row = 0; row < expectedU.length; row++) {
			assertEquals(expectedU[row].length, actualU[row].length);
			for (int col = 0; col < expectedU[row].length; col++)
				assertEquals(expectedU[row][col], actualU[row][col], DELTA);
		}

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
		expectedL = new double[][] { new double[] { 1, 0, 0, 0 },
				new double[] { 2, 1, 0, 0 }, new double[] { 1, -3, 1, 0 },
				new double[] { -4, -4, -4, 1 } };
		expectedU = new double[][] { new double[] { 2, -1, 2, -4 },
				new double[] { 0, -3, 1, -4 }, new double[] { 0, 0, -4, -3 },
				new double[] { 0, 0, 0, -4 } };
		actualLU = Decomposition.lu(a);
		actualL = actualLU[0];
		actualU = actualLU[1];

		assertEquals(expectedL.length, actualL.length);
		for (int row = 0; row < expectedL.length; row++) {
			assertEquals(expectedL[row].length, actualL[row].length);
			for (int col = 0; col < expectedL[row].length; col++)
				assertEquals(expectedL[row][col], actualL[row][col], DELTA);
		}

		assertEquals(expectedU.length, actualU.length);
		for (int row = 0; row < expectedU.length; row++) {
			assertEquals(expectedU[row].length, actualU[row].length);
			for (int col = 0; col < expectedU[row].length; col++)
				assertEquals(expectedU[row][col], actualU[row][col], DELTA);
		}

		actualA = Matrices.multiply(actualL, actualU);
		assertEquals(a.length, actualA.length);
		for (int row = 0; row < a.length; row++) {
			assertEquals(a[row].length, actualA[row].length);
			for (int col = 0; col < a[row].length; col++)
				assertEquals(a[row][col], actualA[row][col], DELTA);
		}

		try {
			// matrix a is null
			testLU = Decomposition.lu(null);
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testLU);
		}

		try {
			// matrix a does contain null column
			testLU = Decomposition.lu(new double[][] { null });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testLU);
		}

		try {
			// matrix a must does not have as many columns as rows
			testLU = Decomposition.lu(new double[][] { new double[] { 0, 0 },
					new double[] { 0, 0 }, new double[] { 0, 0 } });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testLU);
		}
	}

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

	/**
	 * Test method for {@link Decomposition#qr(double[][])}.
	 */
	@Test
	public void testQr() {
		// variables required for testing
		double[][] a;
		double[][] actualA;
		double[][] expectedQ;
		double[][] actualQ;
		double[][] expectedR;
		double[][] actualR;
		double[][] actualI;
		double[][][] actualQR;
		double[][][] testQR = null;

		// first test on 2 by 2 array
		a = new double[][] { new double[] { 14, -6 }, new double[] { 14, -9 },
				new double[] { 7, -6 } };
		expectedQ = new double[][] { new double[] { 2d / 3, 2d / 3 },
				new double[] { 2d / 3, -1d / 3 },
				new double[] { 1d / 3, -2d / 3 } };
		expectedR = new double[][] { new double[] { 21, -12 },
				new double[] { 0, 3 } };
		actualQR = Decomposition.qr(a);
		actualQ = actualQR[0];
		actualR = actualQR[1];

		assertEquals(expectedQ.length, actualQ.length);
		for (int row = 0; row < expectedQ.length; row++) {
			assertEquals(expectedQ[row].length, actualQ[row].length);
			for (int col = 0; col < expectedQ[row].length; col++)
				assertEquals(expectedQ[row][col], actualQ[row][col], DELTA);
		}

		assertEquals(expectedR.length, actualR.length);
		for (int row = 0; row < expectedR.length; row++) {
			assertEquals(expectedR[row].length, actualR[row].length);
			for (int col = 0; col < expectedR[row].length; col++)
				assertEquals(expectedR[row][col], actualR[row][col], DELTA);
		}

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
		expectedQ = new double[][] { new double[] { -0.5, 0.5, -0.5 },
				new double[] { -0.5, -0.5, 0.5 },
				new double[] { -0.5, 0.5, 0.5 },
				new double[] { -0.5, -0.5, -0.5 } };
		expectedR = new double[][] { new double[] { 6, 18, -6 },
				new double[] { 0, 2, -8 }, new double[] { 0, 0, 16 } };
		actualQR = Decomposition.qr(a);
		actualQ = actualQR[0];
		actualR = actualQR[1];

		assertEquals(expectedQ.length, actualQ.length);
		for (int row = 0; row < expectedQ.length; row++) {
			assertEquals(expectedQ[row].length, actualQ[row].length);
			for (int col = 0; col < expectedQ[row].length; col++)
				assertEquals(expectedQ[row][col], actualQ[row][col], DELTA);
		}

		assertEquals(expectedR.length, actualR.length);
		for (int row = 0; row < expectedR.length; row++) {
			assertEquals(expectedR[row].length, actualR[row].length);
			for (int col = 0; col < expectedR[row].length; col++)
				assertEquals(expectedR[row][col], actualR[row][col], DELTA);
		}

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
			testQR = Decomposition.qr(null);
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testQR);
		}

		try {
			// matrix a does contain null column
			testQR = Decomposition.qr(new double[][] { null });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testQR);
		}

		try {
			// matrix a contains as fewer rows than columns
			testQR = Decomposition.qr(new double[][] {
					new double[] { 0, 0, 0 }, new double[] { 0, 0, 0 } });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testQR);
		}

		try {
			// matrix a row width does not remain constant
			testQR = Decomposition.qr(new double[][] { new double[] { 0, 0 },
					new double[] { 0 } });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testQR);
		}
	}

	/**
	 * Test method for {@link Decomposition#eigen(double[][])}.
	 */
	@Test
	public void testEigen() {
		// variables required for testing
		double[][] a;
		double[][] expectedL;
		double[][] actualL;
		double[][] expectedV;
		double[][] actualV;
		double[][] actualAV;
		double[][] actualLV;
		double[][][] actualEigen;
		double[][][] testEigen = null;

		// first test on 1 by 1 array
		a = new double[][] { new double[] { 7 } };
		expectedL = new double[][] { new double[] { 7 } };
		expectedV = new double[][] { new double[] { 1 } };
		actualEigen = Decomposition.eigen(a);
		actualL = actualEigen[0];
		actualV = actualEigen[1];

		assertEquals(expectedL.length, actualL.length);
		for (int row = 0; row < expectedL.length; row++) {
			assertEquals(expectedL[row].length, actualL[row].length);
			for (int col = 0; col < expectedL[row].length; col++)
				assertEquals(expectedL[row][col], actualL[row][col], DELTA);
		}

		assertEquals(expectedV.length, actualV.length);
		for (int row = 0; row < expectedV.length; row++) {
			assertEquals(expectedV[row].length, actualV[row].length);
			for (int col = 0; col < expectedV[row].length; col++)
				assertEquals(expectedV[row][col], actualV[row][col], DELTA);
		}

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
		expectedL = new double[][] { new double[] { 3, 0 },
				new double[] { 0, 2 } };
		expectedV = new double[][] {
				new double[] { 3 / Math.sqrt(13d), 2 / Math.sqrt(5d) },
				new double[] { 2 / Math.sqrt(13d), 1 / Math.sqrt(5d) } };
		actualEigen = Decomposition.eigen(a);
		actualL = actualEigen[0];
		actualV = actualEigen[1];

		assertEquals(expectedL.length, actualL.length);
		for (int row = 0; row < expectedL.length; row++) {
			assertEquals(expectedL[row].length, actualL[row].length);
			for (int col = 0; col < expectedL[row].length; col++) {
				assertEquals(expectedL[row][col], actualL[row][col], DELTA);
			}
		}

		assertEquals(expectedV.length, actualV.length);
		for (int row = 0; row < expectedV.length; row++) {
			assertEquals(expectedV[row].length, actualV[row].length);
			for (int col = 0; col < expectedV[row].length; col++)
				assertEquals(expectedV[row][col], actualV[row][col], DELTA);
		}

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
		expectedL = new double[][] { new double[] { 2, 0 },
				new double[] { 0, -1 } };
		expectedV = new double[][] {
				new double[] { -1 / Math.sqrt(5), -2 / Math.sqrt(5) },
				new double[] { 2 / Math.sqrt(5), 1 / Math.sqrt(5) } };
		actualEigen = Decomposition.eigen(a);
		actualL = actualEigen[0];
		actualV = actualEigen[1];

		assertEquals(expectedL.length, actualL.length);
		for (int row = 0; row < expectedL.length; row++) {
			assertEquals(expectedL[row].length, actualL[row].length);
			for (int col = 0; col < expectedL[row].length; col++)
				assertEquals(expectedL[row][col], actualL[row][col], DELTA);
		}

		assertEquals(expectedV.length, actualV.length);
		for (int row = 0; row < expectedV.length; row++) {
			assertEquals(expectedV[row].length, actualV[row].length);
			for (int col = 0; col < expectedV[row].length; col++)
				assertEquals(expectedV[row][col], actualV[row][col], DELTA);
		}

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
		expectedL = new double[][] { new double[] { 1, 0, 0 },
				new double[] { 0, -1, 0 }, new double[] { 0, 0, -2 } };
		expectedV = new double[][] {
				new double[] { 0, 1 / Math.sqrt(2), 9 / Math.sqrt(185) },
				new double[] { 1 / Math.sqrt(2), 0, -2 / Math.sqrt(185) },
				new double[] { 1 / Math.sqrt(2), 1 / Math.sqrt(2),
						10 / Math.sqrt(185) } };
		actualEigen = Decomposition.eigen(a);
		actualL = actualEigen[0];
		actualV = actualEigen[1];

		assertEquals(expectedL.length, expectedL.length);
		for (int row = 0; row < expectedL.length; row++) {
			assertEquals(expectedL[row].length, expectedL[row].length);
			for (int col = 0; col < expectedL[row].length; col++)
				assertEquals(expectedL[row][col], expectedL[row][col], DELTA);
		}

		assertEquals(expectedV.length, actualV.length);
		for (int row = 0; row < expectedV.length; row++) {
			assertEquals(expectedV[row].length, actualV[row].length);
			for (int col = 0; col < expectedV[row].length; col++)
				assertEquals(expectedV[row][col], actualV[row][col], DELTA);
		}

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
		expectedL = new double[][] { new double[] { 1, 0, 0 },
				new double[] { 0, -1, 0 }, new double[] { 0, 0, -2 } };
		expectedV = new double[][] {
				new double[] { -7 / Math.sqrt(74), -1 / Math.sqrt(3),
						-2 / Math.sqrt(17) },
				new double[] { 3 / Math.sqrt(74), 1 / Math.sqrt(3),
						3 / Math.sqrt(17) },
				new double[] { 4 / Math.sqrt(74), 1 / Math.sqrt(3),
						2 / Math.sqrt(17) } };
		actualEigen = Decomposition.eigen(a);
		actualL = actualEigen[0];
		actualV = actualEigen[1];

		assertEquals(expectedL.length, expectedL.length);
		for (int row = 0; row < expectedL.length; row++) {
			assertEquals(expectedL[row].length, expectedL[row].length);
			for (int col = 0; col < expectedL[row].length; col++)
				assertEquals(expectedL[row][col], expectedL[row][col], DELTA);
		}

		assertEquals(expectedV.length, actualV.length);
		for (int row = 0; row < expectedV.length; row++) {
			assertEquals(expectedV[row].length, actualV[row].length);
			for (int col = 0; col < expectedV[row].length; col++)
				assertEquals(expectedV[row][col], actualV[row][col], DELTA);
		}

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
			testEigen = Decomposition.eigen(null);
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testEigen);
		}

		try {
			// matrix a does contain null column
			testEigen = Decomposition.eigen(new double[][] { null });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testEigen);
		}

		try {
			// matrix a must does not have as many columns as rows
			testEigen = Decomposition.eigen(new double[][] {
					new double[] { 0, 0 }, new double[] { 0, 0 },
					new double[] { 0, 0 } });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testEigen);
		}

		try {
			// matrix a is too large
			testEigen = Decomposition.eigen(new double[][] {
					new double[] { 0, 0, 0, 0 }, new double[] { 0, 0, 0, 0 },
					new double[] { 0, 0, 0, 0 }, new double[] { 0, 0, 0, 0 } });
			fail();
		} catch (UnsupportedOperationException ex) {
			assertNull(testEigen);
		}
	}

	/**
	 * Test method for {@link Decomposition#singularValue(double[][])}.
	 */
	@Test
	public void testSingularValue() {
		// variables required for testing
		double[][] a;
		double[][] actualA;
		double[][] expectedU;
		double[][] actualU;
		double[][] expectedS;
		double[][] actualS;
		double[][] expectedV;
		double[][] actualV;
		double[][] actualI;
		double[][][] actualSingular;
		double[][][] testSingular = null;

		// first test on 2 by 2 array
		a = new double[][] { new double[] { 8, 6 }, new double[] { 3, -4 } };
		expectedU = new double[][] { new double[] { 1, 0 },
				new double[] { 0, -1 } };
		expectedS = new double[][] { new double[] { 10, 0 },
				new double[] { 0, 5 } };
		expectedV = new double[][] { new double[] { 0.8, -0.6 },
				new double[] { 0.6, 0.8 } };
		actualSingular = Decomposition.singularValue(a);
		actualU = actualSingular[0];
		actualS = actualSingular[1];
		actualV = actualSingular[2];

		assertEquals(expectedU.length, actualU.length);
		for (int row = 0; row < expectedU.length; row++) {
			assertEquals(expectedU[row].length, actualU[row].length);
			for (int col = 0; col < expectedU[row].length; col++)
				assertEquals(expectedU[row][col], actualU[row][col], DELTA);
		}

		assertEquals(expectedS.length, actualS.length);
		for (int row = 0; row < expectedS.length; row++) {
			assertEquals(expectedS[row].length, actualS[row].length);
			for (int col = 0; col < expectedS[row].length; col++)
				assertEquals(expectedS[row][col], actualS[row][col], DELTA);
		}

		assertEquals(expectedV.length, actualV.length);
		for (int row = 0; row < expectedV.length; row++) {
			assertEquals(expectedV[row].length, actualV[row].length);
			for (int col = 0; col < expectedV[row].length; col++)
				assertEquals(expectedV[row][col], actualV[row][col], DELTA);
		}

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
		expectedU = new double[][] { new double[] { 2d / 3, 2d / 3 },
				new double[] { 1d / 3, -2d / 3 },
				new double[] { -2d / 3, 1d / 3 } };
		expectedS = new double[][] { new double[] { 45, 0 },
				new double[] { 0, 15 } };
		expectedV = new double[][] { new double[] { 0.8, -0.6 },
				new double[] { 0.6, 0.8 } };
		actualSingular = Decomposition.singularValue(a);
		actualU = actualSingular[0];
		actualS = actualSingular[1];
		actualV = actualSingular[2];

		assertEquals(expectedU.length, actualU.length);
		for (int row = 0; row < expectedU.length; row++) {
			assertEquals(expectedU[row].length, actualU[row].length);
			for (int col = 0; col < expectedU[row].length; col++)
				assertEquals(expectedU[row][col], actualU[row][col], DELTA);
		}

		assertEquals(expectedS.length, actualS.length);
		for (int row = 0; row < expectedS.length; row++) {
			assertEquals(expectedS[row].length, actualS[row].length);
			for (int col = 0; col < expectedS[row].length; col++)
				assertEquals(expectedS[row][col], actualS[row][col], DELTA);
		}

		assertEquals(expectedV.length, actualV.length);
		for (int row = 0; row < expectedV.length; row++) {
			assertEquals(expectedV[row].length, actualV[row].length);
			for (int col = 0; col < expectedV[row].length; col++)
				assertEquals(expectedV[row][col], actualV[row][col], DELTA);
		}

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
		expectedU = new double[][] { new double[] { 1, 0, 0 },
				new double[] { 0, 0, -1 }, new double[] { 0, -1, 0 } };
		expectedS = new double[][] { new double[] { 9, 0, 0 },
				new double[] { 0, 6, 0 }, new double[] { 0, 0, 3 } };
		expectedV = new double[][] { new double[] { 2d / 3, -2d / 3, 1d / 3 },
				new double[] { 2d / 3, 1d / 3, -2d / 3 },
				new double[] { 1d / 3, 2d / 3, 2d / 3 } };
		actualSingular = Decomposition.singularValue(a);
		actualU = actualSingular[0];
		actualS = actualSingular[1];
		actualV = actualSingular[2];

		assertEquals(expectedU.length, actualU.length);
		for (int row = 0; row < expectedU.length; row++) {
			assertEquals(expectedU[row].length, actualU[row].length);
			for (int col = 0; col < expectedU[row].length; col++)
				assertEquals(expectedU[row][col], actualU[row][col], DELTA);
		}

		assertEquals(expectedS.length, actualS.length);
		for (int row = 0; row < expectedS.length; row++) {
			assertEquals(expectedS[row].length, actualS[row].length);
			for (int col = 0; col < expectedS[row].length; col++)
				assertEquals(expectedS[row][col], actualS[row][col], DELTA);
		}

		assertEquals(expectedV.length, actualV.length);
		for (int row = 0; row < expectedV.length; row++) {
			assertEquals(expectedV[row].length, actualV[row].length);
			for (int col = 0; col < expectedV[row].length; col++)
				assertEquals(expectedV[row][col], actualV[row][col], DELTA);
		}

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
		expectedU = new double[][] { new double[] { 0.5, -0.5 },
				new double[] { 0.5, -0.5 }, new double[] { 0.5, 0.5 },
				new double[] { -0.5, -0.5 } };
		expectedS = new double[][] { new double[] { 30, 0 },
				new double[] { 0, 20 } };
		expectedV = new double[][] { new double[] { 0.8, -0.6 },
				new double[] { 0.6, 0.8 } };
		actualSingular = Decomposition.singularValue(a);
		actualU = actualSingular[0];
		actualS = actualSingular[1];
		actualV = actualSingular[2];

		assertEquals(expectedU.length, actualU.length);
		for (int row = 0; row < expectedU.length; row++) {
			assertEquals(expectedU[row].length, actualU[row].length);
			for (int col = 0; col < expectedU[row].length; col++)
				assertEquals(expectedU[row][col], actualU[row][col], DELTA);
		}

		assertEquals(expectedS.length, actualS.length);
		for (int row = 0; row < expectedS.length; row++) {
			assertEquals(expectedS[row].length, actualS[row].length);
			for (int col = 0; col < expectedS[row].length; col++)
				assertEquals(expectedS[row][col], actualS[row][col], DELTA);
		}

		assertEquals(expectedV.length, actualV.length);
		for (int row = 0; row < expectedV.length; row++) {
			assertEquals(expectedV[row].length, actualV[row].length);
			for (int col = 0; col < expectedV[row].length; col++)
				assertEquals(expectedV[row][col], actualV[row][col], DELTA);
		}

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
		expectedU = new double[][] { new double[] { -0.5, 0.5, 0.5 },
				new double[] { -0.5, 0.5, -0.5 },
				new double[] { -0.5, -0.5, 0.5 },
				new double[] { 0.5, 0.5, 0.5 } };
		expectedS = new double[][] { new double[] { 24, 0, 0 },
				new double[] { 0, 12, 0 }, new double[] { 0, 0, 6 } };
		expectedV = new double[][] { new double[] { -2d / 3, 2d / 3, -1d / 3 },
				new double[] { -2d / 3, -1d / 3, 2d / 3 },
				new double[] { 1d / 3, 2d / 3, 2d / 3 } };
		actualSingular = Decomposition.singularValue(a);
		actualU = actualSingular[0];
		actualS = actualSingular[1];
		actualV = actualSingular[2];

		assertEquals(expectedU.length, actualU.length);
		for (int row = 0; row < expectedU.length; row++) {
			assertEquals(expectedU[row].length, actualU[row].length);
			for (int col = 0; col < expectedU[row].length; col++)
				assertEquals(expectedU[row][col], actualU[row][col], DELTA);
		}

		assertEquals(expectedS.length, actualS.length);
		for (int row = 0; row < expectedS.length; row++) {
			assertEquals(expectedS[row].length, actualS[row].length);
			for (int col = 0; col < expectedS[row].length; col++)
				assertEquals(expectedS[row][col], actualS[row][col], DELTA);
		}

		assertEquals(expectedV.length, actualV.length);
		for (int row = 0; row < expectedV.length; row++) {
			assertEquals(expectedV[row].length, actualV[row].length);
			for (int col = 0; col < expectedV[row].length; col++)
				assertEquals(expectedV[row][col], actualV[row][col], DELTA);
		}

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
			testSingular = Decomposition.singularValue(null);
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testSingular);
		}

		try {
			// matrix a does contain null column
			testSingular = Decomposition.singularValue(new double[][] { null });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testSingular);
		}

		try {
			// matrix a contains as fewer rows than columns
			testSingular = Decomposition.singularValue(new double[][] {
					new double[] { 0, 0, 0 }, new double[] { 0, 0, 0 } });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testSingular);
		}

		try {
			// matrix a row width does not remain constant
			testSingular = Decomposition.singularValue(new double[][] {
					new double[] { 0, 0 }, new double[] { 0 } });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testSingular);
		}
	}

}