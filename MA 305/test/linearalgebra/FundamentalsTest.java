package linearalgebra;

import static org.junit.Assert.*;
import linearalgebra.Fundamentals;
import linearalgebra.Matrices;

import org.junit.Test;

/**
 * Fundamentals test.
 * 
 * @author Jacob Malter
 *
 */
public class FundamentalsTest {

	/** acceptable margin of error */
	private static final double DELTA = 0.0001d;

	// Related to Gauss Jordan method

	/**
	 * Test method for {@link Fundamentals#augment(double[][], double[][])}.
	 */
	@Test
	public void testAugment() {
		// variables required for testing
		double[][] a;
		double[][] b;
		double[][] expectedAB;
		double[][] actualAB;
		double[][] testAB = null;

		// first test on 3 by 3 and 3 by 1 arrays
		a = new double[][] { new double[] { -3, 1, 2 },
				new double[] { 9, -4, -8 }, new double[] { -9, 5, 8 } };
		b = new double[][] { new double[] { -3 }, new double[] { 3 },
				new double[] { 1 } };
		expectedAB = new double[][] { new double[] { -3, 1, 2, -3 },
				new double[] { 9, -4, -8, 3 }, new double[] { -9, 5, 8, 1 } };
		actualAB = Fundamentals.augment(a, b);

		assertEquals(expectedAB.length, actualAB.length);
		for (int row = 0; row < expectedAB.length; row++) {
			assertEquals(expectedAB[row].length, actualAB[row].length);
			for (int col = 0; col < expectedAB[row].length; col++)
				assertEquals(expectedAB[row][col], actualAB[row][col], DELTA);
		}

		// second test on 2 by 2 and 2 by 1 arrays
		a = new double[][] { new double[] { -3, 2 }, new double[] { 3, 1 } };
		b = new double[][] { new double[] { -5 }, new double[] { 4 } };
		expectedAB = new double[][] { new double[] { -3, 2, -5 },
				new double[] { 3, 1, 4 } };
		actualAB = Fundamentals.augment(a, b);

		assertEquals(expectedAB.length, actualAB.length);
		for (int row = 0; row < expectedAB.length; row++) {
			assertEquals(expectedAB[row].length, actualAB[row].length);
			for (int col = 0; col < expectedAB[row].length; col++)
				assertEquals(expectedAB[row][col], actualAB[row][col], DELTA);
		}

		try {
			// matrix a is null
			testAB = Fundamentals.augment(null, new double[][] {});
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testAB);
		}

		try {
			// matrix b is null
			testAB = Fundamentals.augment(new double[][] {}, null);
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testAB);
		}

		try {
			// matrix a and b do not have equal number of rows
			testAB = Fundamentals.augment(new double[][] { new double[] {} },
					new double[][] {});
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testAB);
		}

		try {
			// matrix a does contain null column
			testAB = Fundamentals.augment(new double[][] { null },
					new double[][] { null });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testAB);
		}

		try {
			// matrix b does contain null column
			testAB = Fundamentals.augment(new double[][] { new double[] {} },
					new double[][] { null });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testAB);
		}

		try {
			// matrix a row width does not remain constant
			testAB = Fundamentals.augment(new double[][] {
					new double[] { 0, 0 }, new double[] { 0 } },
					new double[][] { new double[] {}, new double[] {} });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testAB);
		}

		try {
			// matrix b row width does not remain constant
			testAB = Fundamentals.augment(new double[][] { new double[] { 0 },
					new double[] { 0 } }, new double[][] {
					new double[] { 0, 0 }, new double[] { 0 } });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testAB);
		}
	}

	/**
	 * Test method for {@link Fundamentals#forwardEliminate(double[][])}.
	 */
	@Test
	public void testForwardEliminate() {
		// variables required for testing
		double[][] a;
		double[][] expectedFor;
		double[][] actualFor;
		double[][] testFor = null;

		// first test on 3 by 4 matrix
		a = new double[][] { new double[] { 3, -6, 6, -6 },
				new double[] { -3, 9, -8, 0 }, new double[] { 3, -9, 6, 6 } };
		expectedFor = new double[][] { new double[] { 3, -6, 6, -6 },
				new double[] { 0, 3, -2, -6 }, new double[] { 0, 0, -2, 6 } };
		actualFor = Fundamentals.forwardEliminate(a);

		assertEquals(expectedFor.length, actualFor.length);
		for (int row = 0; row < expectedFor.length; row++) {
			assertEquals(expectedFor[row].length, actualFor[row].length);
			for (int col = 0; col < expectedFor[row].length; col++)
				assertEquals(expectedFor[row][col], actualFor[row][col], DELTA);
		}

		// second test on 3 by 4 matrix
		a = new double[][] { new double[] { 4, -6, 8, -2 },
				new double[] { -4, 3, -2, -1 }, new double[] { 4, 0, -2, 0 } };
		expectedFor = new double[][] { new double[] { 4, -6, 8, -2 },
				new double[] { 0, -3, 6, -3 }, new double[] { 0, 0, 2, -4 } };
		actualFor = Fundamentals.forwardEliminate(a);

		assertEquals(expectedFor.length, actualFor.length);
		for (int row = 0; row < expectedFor.length; row++) {
			assertEquals(expectedFor[row].length, actualFor[row].length);
			for (int col = 0; col < expectedFor[row].length; col++)
				assertEquals(expectedFor[row][col], actualFor[row][col], DELTA);
		}

		try {
			// matrix a is null
			testFor = Fundamentals.forwardEliminate(null);
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testFor);
		}

		try {
			// matrix a does contain null column
			testFor = Fundamentals.forwardEliminate(new double[][] { null });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testFor);
		}

		try {
			// matrix a must does not have as many or more columns than rows
			testFor = Fundamentals.forwardEliminate(new double[][] {
					new double[] { 0, 0 }, new double[] { 0, 0 },
					new double[] { 0, 0 } });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testFor);
		}

		try {
			// matrix a row width does not remain constant
			testFor = Fundamentals.forwardEliminate(new double[][] {
					new double[] { 0, 0 }, new double[] { 0 } });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testFor);
		}

		try {
			// matrix a row width does not remain constant
			testFor = Fundamentals.forwardEliminate(new double[][] {
					new double[] { 0, 0 }, new double[] { 0 } });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testFor);
		}
	}

	/**
	 * Test method for {@link Fundamentals#backwardEliminate(double[][])}.
	 */
	@Test
	public void testBackwardEliminate() {
		// variables required for testing
		double[][] a;
		double[][] expectedBack;
		double[][] actualBack;
		double[][] testBack = null;

		// first test on 3 by 4 matrix
		a = new double[][] { new double[] { 3, -6, 6, -6 },
				new double[] { 0, 3, -2, -6 }, new double[] { 0, 0, -2, 6 } };
		expectedBack = new double[][] { new double[] { 3, 0, 0, -12 },
				new double[] { 0, 3, 0, -12 }, new double[] { 0, 0, -2, 6 } };
		actualBack = Fundamentals.backwardEliminate(a);

		assertEquals(expectedBack.length, actualBack.length);
		for (int row = 0; row < expectedBack.length; row++) {
			assertEquals(expectedBack[row].length, actualBack[row].length);
			for (int col = 0; col < expectedBack[row].length; col++)
				assertEquals(expectedBack[row][col], actualBack[row][col],
						DELTA);
		}

		// second test on 3 by 4 matrix
		a = new double[][] { new double[] { 4, -6, 8, -2 },
				new double[] { 0, -3, 6, -3 }, new double[] { 0, 0, 2, -4 } };
		expectedBack = new double[][] { new double[] { 4, 0, 0, -4 },
				new double[] { 0, -3, 0, 9 }, new double[] { 0, 0, 2, -4 } };
		actualBack = Fundamentals.backwardEliminate(a);

		assertEquals(expectedBack.length, actualBack.length);
		for (int row = 0; row < expectedBack.length; row++) {
			assertEquals(expectedBack[row].length, actualBack[row].length);
			for (int col = 0; col < expectedBack[row].length; col++)
				assertEquals(expectedBack[row][col], actualBack[row][col],
						DELTA);
		}

		try {
			// matrix a is null
			testBack = Fundamentals.backwardEliminate(null);
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testBack);
		}

		try {
			// matrix a does contain null column
			testBack = Fundamentals.backwardEliminate(new double[][] { null });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testBack);
		}

		try {
			// matrix a must does not have as many or more columns than rows
			testBack = Fundamentals.backwardEliminate(new double[][] {
					new double[] { 0, 0 }, new double[] { 0, 0 },
					new double[] { 0, 0 } });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testBack);
		}

		try {
			// matrix a row width does not remain constant
			testBack = Fundamentals.backwardEliminate(new double[][] {
					new double[] { 0, 0 }, new double[] { 0 } });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testBack);
		}

		try {
			// matrix a row width does not remain constant
			testBack = Fundamentals.backwardEliminate(new double[][] {
					new double[] { 0, 0 }, new double[] { 0 } });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testBack);
		}
	}

	/**
	 * Test method for {@link Fundamentals#gaussJordan(double[][], double[][])}
	 * .
	 */
	@Test
	public void testGaussJordan() {
		// variables required for testing
		double[][] a;
		double[][] b;
		double[][] actualB;
		double[][] expectedX;
		double[][] actualX;
		double[][] testX = null;

		// first test on 3 by 3 and 3 by 1 arrays
		a = new double[][] { new double[] { -3, 1, 2 },
				new double[] { 9, -4, -8 }, new double[] { -9, 5, 8 } };
		b = new double[][] { new double[] { -3 }, new double[] { 3 },
				new double[] { 1 } };
		expectedX = new double[][] { new double[] { 3 }, new double[] { 4 },
				new double[] { 1 } };
		actualX = Fundamentals.gaussJordan(a, b);

		assertEquals(expectedX.length, actualX.length);
		for (int row = 0; row < expectedX.length; row++) {
			assertEquals(expectedX[row].length, actualX[row].length);
			for (int col = 0; col < expectedX[row].length; col++)
				assertEquals(expectedX[row][col], actualX[row][col], DELTA);
		}

		actualB = Matrices.multiply(a, actualX);
		assertEquals(b.length, actualB.length);
		for (int row = 0; row < b.length; row++) {
			assertEquals(b[row].length, actualB[row].length);
			for (int col = 0; col < b[row].length; col++)
				assertEquals(b[row][col], actualB[row][col], DELTA);
		}

		// second test on 2 by 2 and 2 by 1 arrays
		a = new double[][] { new double[] { 4, -6, 8 },
				new double[] { -4, 3, -2 }, new double[] { 4, 0, -2 } };
		b = new double[][] { new double[] { -2 }, new double[] { -1 },
				new double[] { 0 } };
		expectedX = new double[][] { new double[] { -1 }, new double[] { -3 },
				new double[] { -2 } };
		actualX = Fundamentals.gaussJordan(a, b);

		assertEquals(expectedX.length, actualX.length);
		for (int row = 0; row < expectedX.length; row++) {
			assertEquals(expectedX[row].length, actualX[row].length);
			for (int col = 0; col < expectedX[row].length; col++)
				assertEquals(expectedX[row][col], actualX[row][col], DELTA);
		}

		actualB = Matrices.multiply(a, actualX);
		assertEquals(b.length, actualB.length);
		for (int row = 0; row < b.length; row++) {
			assertEquals(b[row].length, actualB[row].length);
			for (int col = 0; col < b[row].length; col++)
				assertEquals(b[row][col], actualB[row][col], DELTA);
		}

		try {
			// matrix a is null
			testX = Fundamentals.gaussJordan(null, new double[][] {});
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testX);
		}

		try {
			// matrix b is null
			testX = Fundamentals.gaussJordan(new double[][] {}, null);
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testX);
		}

		try {
			// matrix a and b do not have equal number of rows
			testX = Fundamentals.gaussJordan(
					new double[][] { new double[] {} }, new double[][] {});
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testX);
		}

		try {
			// matrix a does contain null column
			testX = Fundamentals.gaussJordan(new double[][] { null },
					new double[][] { null });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testX);
		}

		try {
			// matrix b does contain null column
			testX = Fundamentals
					.gaussJordan(new double[][] { new double[] {} },
							new double[][] { null });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testX);
		}

		try {
			// matrix a row width does not remain constant
			testX = Fundamentals.gaussJordan(new double[][] {
					new double[] { 0, 0 }, new double[] { 0 } },
					new double[][] { new double[] {}, new double[] {} });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testX);
		}

		try {
			// matrix b row width does not remain constant
			testX = Fundamentals.gaussJordan(new double[][] {
					new double[] { 0 }, new double[] { 0 } }, new double[][] {
					new double[] { 0, 0 }, new double[] { 0 } });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testX);
		}
	}

	// Related to Cramer's method

	/**
	 * Test method for {@link Fundamentals#cramer(double[][], double[][])}.
	 */
	@Test
	public void testCramer() {
		// variables required for testing
		double[][] a;
		double[][] b;
		double[][] actualB;
		double[][] expectedX;
		double[][] actualX;
		double[][] testX = null;

		// first test on 3 by 3 and 3 by 1 arrays
		a = new double[][] { new double[] { -3, 1, 2 },
				new double[] { 9, -4, -8 }, new double[] { -9, 5, 8 } };
		b = new double[][] { new double[] { -3 }, new double[] { 3 },
				new double[] { 1 } };
		expectedX = new double[][] { new double[] { 3 }, new double[] { 4 },
				new double[] { 1 } };
		actualX = Fundamentals.cramer(a, b);

		assertEquals(expectedX.length, actualX.length);
		for (int row = 0; row < expectedX.length; row++) {
			assertEquals(expectedX[row].length, actualX[row].length);
			for (int col = 0; col < expectedX[row].length; col++)
				assertEquals(expectedX[row][col], actualX[row][col], DELTA);
		}

		actualB = Matrices.multiply(a, actualX);
		assertEquals(b.length, actualB.length);
		for (int row = 0; row < b.length; row++) {
			assertEquals(b[row].length, actualB[row].length);
			for (int col = 0; col < b[row].length; col++)
				assertEquals(b[row][col], actualB[row][col], DELTA);
		}

		// second test on 2 by 2 and 2 by 1 arrays
		a = new double[][] { new double[] { 4, -6, 8 },
				new double[] { -4, 3, -2 }, new double[] { 4, 0, -2 } };
		b = new double[][] { new double[] { -2 }, new double[] { -1 },
				new double[] { 0 } };
		expectedX = new double[][] { new double[] { -1 }, new double[] { -3 },
				new double[] { -2 } };
		actualX = Fundamentals.cramer(a, b);

		assertEquals(expectedX.length, actualX.length);
		for (int row = 0; row < expectedX.length; row++) {
			assertEquals(expectedX[row].length, actualX[row].length);
			for (int col = 0; col < expectedX[row].length; col++)
				assertEquals(expectedX[row][col], actualX[row][col], DELTA);
		}

		actualB = Matrices.multiply(a, actualX);
		assertEquals(b.length, actualB.length);
		for (int row = 0; row < b.length; row++) {
			assertEquals(b[row].length, actualB[row].length);
			for (int col = 0; col < b[row].length; col++)
				assertEquals(b[row][col], actualB[row][col], DELTA);
		}

		try {
			// matrix a is null
			testX = Fundamentals.cramer(null, new double[][] {});
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testX);
		}

		try {
			// matrix b is null
			testX = Fundamentals.cramer(new double[][] {}, null);
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testX);
		}

		try {
			// matrix a and b do not have equal number of rows
			testX = Fundamentals.cramer(new double[][] { new double[] {} },
					new double[][] {});
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testX);
		}

		try {
			// matrix a does contain null column
			testX = Fundamentals.cramer(new double[][] { null },
					new double[][] { null });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testX);
		}

		try {
			// matrix b does contain null column
			testX = Fundamentals.cramer(new double[][] { new double[] {} },
					new double[][] { null });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testX);
		}

		try {
			// matrix a row width does not remain constant
			testX = Fundamentals.cramer(new double[][] { new double[] { 0, 0 },
					new double[] { 0 } }, new double[][] { new double[] {},
					new double[] {} });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testX);
		}

		try {
			// matrix b row width does not remain constant
			testX = Fundamentals.cramer(new double[][] { new double[] { 0 },
					new double[] { 0 } }, new double[][] {
					new double[] { 0, 0 }, new double[] { 0 } });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testX);
		}
	}

	// Related to matrix inversion

	/**
	 * Test method for {@link Fundamentals#invert(double[][])}.
	 */
	@Test
	public void testInvert() {
		// variables required for testing
		double[][] a;
		double[][] actualI;
		double[][] expectedInv;
		double[][] actualInv;
		double[][] testInv = null;

		// first test on 2 by 2 array
		a = new double[][] { new double[] { -1, 4 }, new double[] { 1, -5 } };
		expectedInv = new double[][] { new double[] { -5, -4 },
				new double[] { -1, -1 } };
		actualInv = Fundamentals.invert(a);

		assertEquals(expectedInv.length, actualInv.length);
		for (int row = 0; row < expectedInv.length; row++) {
			assertEquals(expectedInv[row].length, actualInv[row].length);
			for (int col = 0; col < expectedInv[row].length; col++)
				assertEquals(expectedInv[row][col], actualInv[row][col], DELTA);
		}

		actualI = Matrices.multiply(a, actualInv);
		for (int row = 0; row < actualI.length; row++)
			for (int col = 0; col < actualI[row].length; col++)
				if (row == col)
					assertEquals(1, actualI[row][col], DELTA);

				else
					assertEquals(0, actualI[row][col], DELTA);

		// second test on 3 by 3 array
		a = new double[][] { new double[] { -1, 4, -4 },
				new double[] { -1, 0, -2 }, new double[] { -3, 0, -4 } };
		expectedInv = new double[][] { new double[] { 0, 2, -1 },
				new double[] { 0.25, -1, 0.25 }, new double[] { 0, -1.5, 0.5 } };
		actualInv = Fundamentals.invert(a);

		assertEquals(expectedInv.length, actualInv.length);
		for (int row = 0; row < expectedInv.length; row++) {
			assertEquals(expectedInv[row].length, actualInv[row].length);
			for (int col = 0; col < expectedInv[row].length; col++)
				assertEquals(expectedInv[row][col], actualInv[row][col], DELTA);
		}

		actualI = Matrices.multiply(a, actualInv);
		for (int row = 0; row < actualI.length; row++)
			for (int col = 0; col < actualI[row].length; col++)
				if (row == col)
					assertEquals(1, actualI[row][col], DELTA);

				else
					assertEquals(0, actualI[row][col], DELTA);

		try {
			// matrix a is null
			testInv = Fundamentals.invert(null);
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testInv);
		}

		try {
			// matrix a does contain null column
			testInv = Fundamentals.invert(new double[][] { null });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testInv);
		}

		try {
			// matrix a row width does not remain constant
			testInv = Fundamentals.invert(new double[][] {
					new double[] { 0, 0 }, new double[] { 0 } });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testInv);
		}

		try {
			// divide by 0
			testInv = Fundamentals.invert(new double[][] {
					new double[] { 0, 0 }, new double[] { 0, 0 } });
			fail();
		} catch (ArithmeticException ex) {
			assertNull(testInv);
		}
	}

}