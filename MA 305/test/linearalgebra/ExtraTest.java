/**
 * 
 */
package linearalgebra;

import static org.junit.Assert.*;

import org.junit.Test;

/**
 * Extra test.
 * 
 * @author Jacob Malter
 *
 */
public class ExtraTest {

	/** acceptable margin of error */
	private static final double DELTA = 0.001d;

	/**
	 * Test method for {@link Extra#overconstrainedPI(double[][], double[][])}.
	 */
	@Test
	public void testOverconstrainedPI() {
		double[][] a;
		double[][] b;
		double[][][] actualXE;
		double[][] expectedX;
		double[][] actualX;
		double[][] expectedE;
		double[][] actualE;
		double[][] actualAX;
		double[][][] testXE = null;

		// first test on 3 by 2 and 3 by 1 arrays
		a = new double[][] { new double[] { 4, -3 }, new double[] { -1, 1 },
				new double[] { -4, -3 } };
		b = new double[][] { new double[] { -1 }, new double[] { -3 },
				new double[] { 4 } };
		expectedX = new double[][] { new double[] { -0.535 },
				new double[] { -0.660 } };
		expectedE = new double[][] { new double[] { 0.839 },
				new double[] { 2.875 }, new double[] { 0.120 } };
		actualXE = Extra.overconstrainedPI(a, b);
		actualX = actualXE[0];
		actualE = actualXE[1];

		assertEquals(expectedX.length, actualX.length);
		for (int row = 0; row < expectedX.length; row++) {
			assertEquals(expectedX[row].length, actualX[row].length);
			for (int col = 0; col < expectedX[row].length; col++)
				assertEquals(expectedX[row][col], actualX[row][col], DELTA);
		}

		assertEquals(expectedE.length, actualE.length);
		for (int row = 0; row < expectedE.length; row++) {
			assertEquals(expectedE[row].length, actualE[row].length);
			for (int col = 0; col < expectedE[row].length; col++)
				assertEquals(expectedE[row][col], actualE[row][col], DELTA);
		}

		actualAX = Matrices.multiply(a, actualX);
		assertEquals(expectedE.length, actualAX.length);
		assertEquals(expectedE.length, b.length);
		for (int row = 0; row < expectedE.length; row++) {
			assertEquals(expectedE[row].length, actualAX[row].length);
			assertEquals(expectedE[row].length, b[row].length);
			for (int col = 0; col < expectedE[row].length; col++)
				assertEquals(expectedE[row][col], actualAX[row][col]
						- b[row][col], DELTA);
		}

		// second test on 4 by 3 and 4 by 1 arrays
		a = new double[][] { new double[] { -3, 4, 3 },
				new double[] { -3, -2, -4 }, new double[] { 4, 2, -1 },
				new double[] { 3, 1, -3 } };
		b = new double[][] { new double[] { 2 }, new double[] { -3 },
				new double[] { 1 }, new double[] { 1 } };
		expectedX = new double[][] { new double[] { 0.239 },
				new double[] { 0.473 }, new double[] { 0.266 } };
		expectedE = new double[][] { new double[] { -0.030 },
				new double[] { 0.273 }, new double[] { 0.637 },
				new double[] { -0.607 } };
		actualXE = Extra.overconstrainedPI(a, b);
		actualX = actualXE[0];
		actualE = actualXE[1];

		assertEquals(expectedX.length, actualX.length);
		for (int row = 0; row < expectedX.length; row++) {
			assertEquals(expectedX[row].length, actualX[row].length);
			for (int col = 0; col < expectedX[row].length; col++)
				assertEquals(expectedX[row][col], actualX[row][col], DELTA);
		}

		assertEquals(expectedE.length, actualE.length);
		for (int row = 0; row < expectedE.length; row++) {
			assertEquals(expectedE[row].length, actualE[row].length);
			for (int col = 0; col < expectedE[row].length; col++)
				assertEquals(expectedE[row][col], actualE[row][col], DELTA);
		}

		actualAX = Matrices.multiply(a, actualX);
		assertEquals(expectedE.length, actualAX.length);
		assertEquals(expectedE.length, b.length);
		for (int row = 0; row < expectedE.length; row++) {
			assertEquals(expectedE[row].length, actualAX[row].length);
			assertEquals(expectedE[row].length, b[row].length);
			for (int col = 0; col < expectedE[row].length; col++)
				assertEquals(expectedE[row][col], actualAX[row][col]
						- b[row][col], DELTA);
		}

		// third test on 5 by 3 and 5 by 1 arrays
		a = new double[][] { new double[] { 3, -1, -1 },
				new double[] { -3, 2, -3 }, new double[] { 2, -4, -2 },
				new double[] { -4, -3, -1 }, new double[] { 3, 1, -1 } };
		b = new double[][] { new double[] { 3 }, new double[] { 3 },
				new double[] { 1 }, new double[] { 1 }, new double[] { 4 } };
		expectedX = new double[][] { new double[] { 0.307 },
				new double[] { 0.232 }, new double[] { -1.318 } };
		expectedE = new double[][] { new double[] { -0.994 },
				new double[] { 0.497 }, new double[] { 1.319 },
				new double[] { -1.606 }, new double[] { -1.530 } };
		actualXE = Extra.overconstrainedPI(a, b);
		actualX = actualXE[0];
		actualE = actualXE[1];

		assertEquals(expectedX.length, actualX.length);
		for (int row = 0; row < expectedX.length; row++) {
			assertEquals(expectedX[row].length, actualX[row].length);
			for (int col = 0; col < expectedX[row].length; col++)
				assertEquals(expectedX[row][col], actualX[row][col], DELTA);
		}

		assertEquals(expectedE.length, actualE.length);
		for (int row = 0; row < expectedE.length; row++) {
			assertEquals(expectedE[row].length, actualE[row].length);
			for (int col = 0; col < expectedE[row].length; col++)
				assertEquals(expectedE[row][col], actualE[row][col], DELTA);
		}

		actualAX = Matrices.multiply(a, actualX);
		assertEquals(expectedE.length, actualAX.length);
		assertEquals(expectedE.length, b.length);
		for (int row = 0; row < expectedE.length; row++) {
			assertEquals(expectedE[row].length, actualAX[row].length);
			assertEquals(expectedE[row].length, b[row].length);
			for (int col = 0; col < expectedE[row].length; col++)
				assertEquals(expectedE[row][col], actualAX[row][col]
						- b[row][col], DELTA);
		}

		try {
			// matrix a is null
			testXE = Extra.overconstrainedPI(null, b);
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testXE);
		}

		try {
			// matrix b is null
			testXE = Extra.overconstrainedPI(a, null);
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testXE);
		}

		try {
			// matrix a and b do not have equal number of rows
			testXE = Extra.overconstrainedPI(
					new double[][] { new double[] {} }, new double[][] {});
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testXE);
		}

		try {
			// matrix a does contain null column
			testXE = Extra.overconstrainedPI(new double[][] { null },
					new double[][] { null });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testXE);
		}

		try {
			// matrix b does contain null column
			testXE = Extra
					.overconstrainedPI(new double[][] { new double[] {} },
							new double[][] { null });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testXE);
		}

		try {
			// matrix a must does not have more rows than columns
			testXE = Extra.overconstrainedPI(new double[][] {
					new double[] { 0, 0 }, new double[] { 0, 0 } },
					new double[][] { new double[] { 0 }, new double[] { 0 } });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testXE);
		}

		try {
			// matrix a row width does not remain constant
			testXE = Extra.overconstrainedPI(new double[][] {
					new double[] { 0, 0 }, new double[] { 0, 0 },
					new double[] { 0, 0, 0 } },
					new double[][] { new double[] { 0 }, new double[] { 0 },
							new double[] { 0 } });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testXE);
		}

		try {
			// matrix b row width does not remain constant
			testXE = Extra.overconstrainedPI(new double[][] {
					new double[] { 0, 0 }, new double[] { 0, 0 },
					new double[] { 0, 0 } }, new double[][] {
					new double[] { 0 }, new double[] { 0 },
					new double[] { 0, 0 } });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testXE);
		}
	}

	/**
	 * Test method for {@link Extra#underconstrainedPI(double[][], double[][])}.
	 */
	@Test
	public void testUnderconstrainedPI() {
		double[][] a;
		double[][] b;
		double[][] actualAX;
		double[][] expectedX;
		double[][] actualX;
		double[][] testX = null;

		// first test on 2 by 3 and 2 by 1 arrays
		a = new double[][] { new double[] { 4, -3, -1 },
				new double[] { 1, -4, -3 } };
		b = new double[][] { new double[] { -1 }, new double[] { -3 } };
		expectedX = new double[][] { new double[] { 0.206 },
				new double[] { 0.454 }, new double[] { 0.463 } };
		actualX = Extra.underconstrainedPI(a, b);

		assertEquals(expectedX.length, actualX.length);
		for (int row = 0; row < expectedX.length; row++) {
			assertEquals(expectedX[row].length, actualX[row].length);
			for (int col = 0; col < expectedX[row].length; col++)
				assertEquals(expectedX[row][col], actualX[row][col], DELTA);
		}

		actualAX = Matrices.multiply(a, actualX);
		assertEquals(b.length, actualAX.length);
		assertEquals(b.length, b.length);
		for (int row = 0; row < b.length; row++) {
			assertEquals(b[row].length, actualAX[row].length);
			assertEquals(b[row].length, b[row].length);
			for (int col = 0; col < b[row].length; col++)
				assertEquals(b[row][col], actualAX[row][col], DELTA);
		}

		// second test on 3 by 4 and 3 by 1 arrays
		a = new double[][] { new double[] { 4, -3, -1, 1 },
				new double[] { -4, -3, -1, -3 }, new double[] { 4, 2, 3, 3 } };
		b = new double[][] { new double[] { 3 }, new double[] { 1 },
				new double[] { -1 } };
		expectedX = new double[][] { new double[] { 0.279 },
				new double[] { -0.555 }, new double[] { -0.277 },
				new double[] { -0.058 } };
		actualX = Extra.underconstrainedPI(a, b);

		assertEquals(expectedX.length, actualX.length);
		for (int row = 0; row < expectedX.length; row++) {
			assertEquals(expectedX[row].length, actualX[row].length);
			for (int col = 0; col < expectedX[row].length; col++)
				assertEquals(expectedX[row][col], actualX[row][col], DELTA);
		}

		actualAX = Matrices.multiply(a, actualX);
		assertEquals(b.length, actualAX.length);
		assertEquals(b.length, b.length);
		for (int row = 0; row < b.length; row++) {
			assertEquals(b[row].length, actualAX[row].length);
			assertEquals(b[row].length, b[row].length);
			for (int col = 0; col < b[row].length; col++)
				assertEquals(b[row][col], actualAX[row][col], DELTA);
		}

		try {
			// matrix a is null
			testX = Extra.underconstrainedPI(null, b);
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testX);
		}

		try {
			// matrix b is null
			testX = Extra.underconstrainedPI(a, null);
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testX);
		}

		try {
			// matrix a and b do not have equal number of rows
			testX = Extra.underconstrainedPI(
					new double[][] { new double[] {} }, new double[][] {});
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testX);
		}

		try {
			// matrix a does contain null column
			testX = Extra.underconstrainedPI(new double[][] { null },
					new double[][] { null });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testX);
		}

		try {
			// matrix b does contain null column
			testX = Extra
					.underconstrainedPI(new double[][] { new double[] {} },
							new double[][] { null });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testX);
		}

		try {
			// matrix a must does not have more rows than columns
			testX = Extra.underconstrainedPI(new double[][] {
					new double[] { 0, 0 }, new double[] { 0, 0 } },
					new double[][] { new double[] { 0 }, new double[] { 0 } });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testX);
		}

		try {
			// matrix a row width does not remain constant
			testX = Extra.underconstrainedPI(new double[][] {
					new double[] { 0, 0, 0 }, new double[] { 0, 0, 0, 0 } },
					new double[][] { new double[] { 0 }, new double[] { 0 } });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testX);
		}

		try {
			// matrix b row width does not remain constant
			testX = Extra
					.underconstrainedPI(
							new double[][] { new double[] { 0, 0, 0 },
									new double[] { 0, 0, 0 } }, new double[][] {
									new double[] { 0 }, new double[] { 0, 0 } });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testX);
		}
	}

	/**
	 * Test method for {@link Extra#overconstrainedCD(double[][], double[][])}.
	 */
	@Test
	public void testOverconstrainedCD() {
		double[][] a;
		double[][] b;
		double[][][] actualXE;
		double[][] expectedX;
		double[][] actualX;
		double[][] expectedE;
		double[][] actualE;
		double[][] actualAX;
		double[][][] testXE = null;

		// first test on 3 by 2 and 3 by 1 arrays
		a = new double[][] { new double[] { 4, -3 }, new double[] { -1, 1 },
				new double[] { -4, -3 } };
		b = new double[][] { new double[] { -1 }, new double[] { -3 },
				new double[] { 4 } };
		expectedX = new double[][] { new double[] { -0.535 },
				new double[] { -0.660 } };
		expectedE = new double[][] { new double[] { 0.839 },
				new double[] { 2.875 }, new double[] { 0.120 } };
		actualXE = Extra.overconstrainedCD(a, b);
		actualX = actualXE[0];
		actualE = actualXE[1];

		assertEquals(expectedX.length, actualX.length);
		for (int row = 0; row < expectedX.length; row++) {
			assertEquals(expectedX[row].length, actualX[row].length);
			for (int col = 0; col < expectedX[row].length; col++)
				assertEquals(expectedX[row][col], actualX[row][col], DELTA);
		}

		assertEquals(expectedE.length, actualE.length);
		for (int row = 0; row < expectedE.length; row++) {
			assertEquals(expectedE[row].length, actualE[row].length);
			for (int col = 0; col < expectedE[row].length; col++)
				assertEquals(expectedE[row][col], actualE[row][col], DELTA);
		}

		actualAX = Matrices.multiply(a, actualX);
		assertEquals(expectedE.length, actualAX.length);
		assertEquals(expectedE.length, b.length);
		for (int row = 0; row < expectedE.length; row++) {
			assertEquals(expectedE[row].length, actualAX[row].length);
			assertEquals(expectedE[row].length, b[row].length);
			for (int col = 0; col < expectedE[row].length; col++)
				assertEquals(expectedE[row][col], actualAX[row][col]
						- b[row][col], DELTA);
		}

		// second test on 4 by 3 and 4 by 1 arrays
		a = new double[][] { new double[] { -3, 4, 3 },
				new double[] { -3, -2, -4 }, new double[] { 4, 2, -1 },
				new double[] { 3, 1, -3 } };
		b = new double[][] { new double[] { 2 }, new double[] { -3 },
				new double[] { 1 }, new double[] { 1 } };
		expectedX = new double[][] { new double[] { 0.239 },
				new double[] { 0.473 }, new double[] { 0.266 } };
		expectedE = new double[][] { new double[] { -0.030 },
				new double[] { 0.273 }, new double[] { 0.637 },
				new double[] { -0.607 } };
		actualXE = Extra.overconstrainedCD(a, b);
		actualX = actualXE[0];
		actualE = actualXE[1];

		assertEquals(expectedX.length, actualX.length);
		for (int row = 0; row < expectedX.length; row++) {
			assertEquals(expectedX[row].length, actualX[row].length);
			for (int col = 0; col < expectedX[row].length; col++)
				assertEquals(expectedX[row][col], actualX[row][col], DELTA);
		}

		assertEquals(expectedE.length, actualE.length);
		for (int row = 0; row < expectedE.length; row++) {
			assertEquals(expectedE[row].length, actualE[row].length);
			for (int col = 0; col < expectedE[row].length; col++)
				assertEquals(expectedE[row][col], actualE[row][col], DELTA);
		}

		actualAX = Matrices.multiply(a, actualX);
		assertEquals(expectedE.length, actualAX.length);
		assertEquals(expectedE.length, b.length);
		for (int row = 0; row < expectedE.length; row++) {
			assertEquals(expectedE[row].length, actualAX[row].length);
			assertEquals(expectedE[row].length, b[row].length);
			for (int col = 0; col < expectedE[row].length; col++)
				assertEquals(expectedE[row][col], actualAX[row][col]
						- b[row][col], DELTA);
		}

		// third test on 5 by 3 and 5 by 1 arrays
		a = new double[][] { new double[] { 3, -1, -1 },
				new double[] { -3, 2, -3 }, new double[] { 2, -4, -2 },
				new double[] { -4, -3, -1 }, new double[] { 3, 1, -1 } };
		b = new double[][] { new double[] { 3 }, new double[] { 3 },
				new double[] { 1 }, new double[] { 1 }, new double[] { 4 } };
		expectedX = new double[][] { new double[] { 0.307 },
				new double[] { 0.232 }, new double[] { -1.318 } };
		expectedE = new double[][] { new double[] { -0.994 },
				new double[] { 0.497 }, new double[] { 1.319 },
				new double[] { -1.606 }, new double[] { -1.530 } };
		actualXE = Extra.overconstrainedCD(a, b);
		actualX = actualXE[0];
		actualE = actualXE[1];

		assertEquals(expectedX.length, actualX.length);
		for (int row = 0; row < expectedX.length; row++) {
			assertEquals(expectedX[row].length, actualX[row].length);
			for (int col = 0; col < expectedX[row].length; col++)
				assertEquals(expectedX[row][col], actualX[row][col], DELTA);
		}

		assertEquals(expectedE.length, actualE.length);
		for (int row = 0; row < expectedE.length; row++) {
			assertEquals(expectedE[row].length, actualE[row].length);
			for (int col = 0; col < expectedE[row].length; col++)
				assertEquals(expectedE[row][col], actualE[row][col], DELTA);
		}

		actualAX = Matrices.multiply(a, actualX);
		assertEquals(expectedE.length, actualAX.length);
		assertEquals(expectedE.length, b.length);
		for (int row = 0; row < expectedE.length; row++) {
			assertEquals(expectedE[row].length, actualAX[row].length);
			assertEquals(expectedE[row].length, b[row].length);
			for (int col = 0; col < expectedE[row].length; col++)
				assertEquals(expectedE[row][col], actualAX[row][col]
						- b[row][col], DELTA);
		}

		try {
			// matrix a is null
			testXE = Extra.overconstrainedCD(null, b);
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testXE);
		}

		try {
			// matrix b is null
			testXE = Extra.overconstrainedCD(a, null);
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testXE);
		}

		try {
			// matrix a and b do not have equal number of rows
			testXE = Extra.overconstrainedCD(
					new double[][] { new double[] {} }, new double[][] {});
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testXE);
		}

		try {
			// matrix a does contain null column
			testXE = Extra.overconstrainedCD(new double[][] { null },
					new double[][] { null });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testXE);
		}

		try {
			// matrix b does contain null column
			testXE = Extra
					.overconstrainedCD(new double[][] { new double[] {} },
							new double[][] { null });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testXE);
		}

		try {
			// matrix a must does not have more rows than columns
			testXE = Extra.overconstrainedCD(new double[][] {
					new double[] { 0, 0 }, new double[] { 0, 0 } },
					new double[][] { new double[] { 0 }, new double[] { 0 } });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testXE);
		}

		try {
			// matrix a row width does not remain constant
			testXE = Extra.overconstrainedCD(new double[][] {
					new double[] { 0, 0 }, new double[] { 0, 0 },
					new double[] { 0, 0, 0 } },
					new double[][] { new double[] { 0 }, new double[] { 0 },
							new double[] { 0 } });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testXE);
		}

		try {
			// matrix b row width does not remain constant
			testXE = Extra.overconstrainedCD(new double[][] {
					new double[] { 0, 0 }, new double[] { 0, 0 },
					new double[] { 0, 0 } }, new double[][] {
					new double[] { 0 }, new double[] { 0 },
					new double[] { 0, 0 } });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testXE);
		}
	}

	/**
	 * Test method for {@link Extra#overconstrainedQR(double[][], double[][])}.
	 */
	@Test
	public void testOverconstrainedQR() {
		double[][] a;
		double[][] b;
		double[][][] actualXE;
		double[][] expectedX;
		double[][] actualX;
		double[][] expectedE;
		double[][] actualE;
		double[][] actualAX;
		double[][][] testXE = null;

		// first test on 3 by 2 and 3 by 1 arrays
		a = new double[][] { new double[] { 4, -3 }, new double[] { -1, 1 },
				new double[] { -4, -3 } };
		b = new double[][] { new double[] { -1 }, new double[] { -3 },
				new double[] { 4 } };
		expectedX = new double[][] { new double[] { -0.535 },
				new double[] { -0.660 } };
		expectedE = new double[][] { new double[] { 0.839 },
				new double[] { 2.875 }, new double[] { 0.120 } };
		actualXE = Extra.overconstrainedQR(a, b);
		actualX = actualXE[0];
		actualE = actualXE[1];

		assertEquals(expectedX.length, actualX.length);
		for (int row = 0; row < expectedX.length; row++) {
			assertEquals(expectedX[row].length, actualX[row].length);
			for (int col = 0; col < expectedX[row].length; col++)
				assertEquals(expectedX[row][col], actualX[row][col], DELTA);
		}

		assertEquals(expectedE.length, actualE.length);
		for (int row = 0; row < expectedE.length; row++) {
			assertEquals(expectedE[row].length, actualE[row].length);
			for (int col = 0; col < expectedE[row].length; col++)
				assertEquals(expectedE[row][col], actualE[row][col], DELTA);
		}

		actualAX = Matrices.multiply(a, actualX);
		assertEquals(expectedE.length, actualAX.length);
		assertEquals(expectedE.length, b.length);
		for (int row = 0; row < expectedE.length; row++) {
			assertEquals(expectedE[row].length, actualAX[row].length);
			assertEquals(expectedE[row].length, b[row].length);
			for (int col = 0; col < expectedE[row].length; col++)
				assertEquals(expectedE[row][col], actualAX[row][col]
						- b[row][col], DELTA);
		}

		// second test on 4 by 3 and 4 by 1 arrays
		a = new double[][] { new double[] { -3, 4, 3 },
				new double[] { -3, -2, -4 }, new double[] { 4, 2, -1 },
				new double[] { 3, 1, -3 } };
		b = new double[][] { new double[] { 2 }, new double[] { -3 },
				new double[] { 1 }, new double[] { 1 } };
		expectedX = new double[][] { new double[] { 0.239 },
				new double[] { 0.473 }, new double[] { 0.266 } };
		expectedE = new double[][] { new double[] { -0.030 },
				new double[] { 0.273 }, new double[] { 0.637 },
				new double[] { -0.607 } };
		actualXE = Extra.overconstrainedQR(a, b);
		actualX = actualXE[0];
		actualE = actualXE[1];

		assertEquals(expectedX.length, actualX.length);
		for (int row = 0; row < expectedX.length; row++) {
			assertEquals(expectedX[row].length, actualX[row].length);
			for (int col = 0; col < expectedX[row].length; col++)
				assertEquals(expectedX[row][col], actualX[row][col], DELTA);
		}

		assertEquals(expectedE.length, actualE.length);
		for (int row = 0; row < expectedE.length; row++) {
			assertEquals(expectedE[row].length, actualE[row].length);
			for (int col = 0; col < expectedE[row].length; col++)
				assertEquals(expectedE[row][col], actualE[row][col], DELTA);
		}

		actualAX = Matrices.multiply(a, actualX);
		assertEquals(expectedE.length, actualAX.length);
		assertEquals(expectedE.length, b.length);
		for (int row = 0; row < expectedE.length; row++) {
			assertEquals(expectedE[row].length, actualAX[row].length);
			assertEquals(expectedE[row].length, b[row].length);
			for (int col = 0; col < expectedE[row].length; col++)
				assertEquals(expectedE[row][col], actualAX[row][col]
						- b[row][col], DELTA);
		}

		// third test on 5 by 3 and 5 by 1 arrays
		a = new double[][] { new double[] { 3, -1, -1 },
				new double[] { -3, 2, -3 }, new double[] { 2, -4, -2 },
				new double[] { -4, -3, -1 }, new double[] { 3, 1, -1 } };
		b = new double[][] { new double[] { 3 }, new double[] { 3 },
				new double[] { 1 }, new double[] { 1 }, new double[] { 4 } };
		expectedX = new double[][] { new double[] { 0.307 },
				new double[] { 0.232 }, new double[] { -1.318 } };
		expectedE = new double[][] { new double[] { -0.994 },
				new double[] { 0.497 }, new double[] { 1.319 },
				new double[] { -1.606 }, new double[] { -1.530 } };
		actualXE = Extra.overconstrainedQR(a, b);
		actualX = actualXE[0];
		actualE = actualXE[1];

		assertEquals(expectedX.length, actualX.length);
		for (int row = 0; row < expectedX.length; row++) {
			assertEquals(expectedX[row].length, actualX[row].length);
			for (int col = 0; col < expectedX[row].length; col++)
				assertEquals(expectedX[row][col], actualX[row][col], DELTA);
		}

		assertEquals(expectedE.length, actualE.length);
		for (int row = 0; row < expectedE.length; row++) {
			assertEquals(expectedE[row].length, actualE[row].length);
			for (int col = 0; col < expectedE[row].length; col++)
				assertEquals(expectedE[row][col], actualE[row][col], DELTA);
		}

		actualAX = Matrices.multiply(a, actualX);
		assertEquals(expectedE.length, actualAX.length);
		assertEquals(expectedE.length, b.length);
		for (int row = 0; row < expectedE.length; row++) {
			assertEquals(expectedE[row].length, actualAX[row].length);
			assertEquals(expectedE[row].length, b[row].length);
			for (int col = 0; col < expectedE[row].length; col++)
				assertEquals(expectedE[row][col], actualAX[row][col]
						- b[row][col], DELTA);
		}

		try {
			// matrix a is null
			testXE = Extra.overconstrainedQR(null, b);
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testXE);
		}

		try {
			// matrix b is null
			testXE = Extra.overconstrainedQR(a, null);
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testXE);
		}

		try {
			// matrix a and b do not have equal number of rows
			testXE = Extra.overconstrainedQR(
					new double[][] { new double[] {} }, new double[][] {});
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testXE);
		}

		try {
			// matrix a does contain null column
			testXE = Extra.overconstrainedQR(new double[][] { null },
					new double[][] { null });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testXE);
		}

		try {
			// matrix b does contain null column
			testXE = Extra
					.overconstrainedQR(new double[][] { new double[] {} },
							new double[][] { null });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testXE);
		}

		try {
			// matrix a must does not have more rows than columns
			testXE = Extra.overconstrainedQR(new double[][] {
					new double[] { 0, 0 }, new double[] { 0, 0 } },
					new double[][] { new double[] { 0 }, new double[] { 0 } });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testXE);
		}

		try {
			// matrix a row width does not remain constant
			testXE = Extra.overconstrainedQR(new double[][] {
					new double[] { 0, 0 }, new double[] { 0, 0 },
					new double[] { 0, 0, 0 } },
					new double[][] { new double[] { 0 }, new double[] { 0 },
							new double[] { 0 } });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testXE);
		}

		try {
			// matrix b row width does not remain constant
			testXE = Extra.overconstrainedQR(new double[][] {
					new double[] { 0, 0 }, new double[] { 0, 0 },
					new double[] { 0, 0 } }, new double[][] {
					new double[] { 0 }, new double[] { 0 },
					new double[] { 0, 0 } });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testXE);
		}
	}

	/**
	 * Test method for {@link Extra#overconstrainedSV(double[][], double[][])}.
	 */
	@Test
	public void testOverconstrainedSV() {
		double[][] a;
		double[][] b;
		double[][][] actualXE;
		double[][] expectedX;
		double[][] actualX;
		double[][] expectedE;
		double[][] actualE;
		double[][] actualAX;
		double[][][] testXE = null;

		// first test on 3 by 2 and 3 by 1 arrays
		a = new double[][] { new double[] { 4, -3 }, new double[] { -1, 1 },
				new double[] { -4, -3 } };
		b = new double[][] { new double[] { -1 }, new double[] { -3 },
				new double[] { 4 } };
		expectedX = new double[][] { new double[] { -0.535 },
				new double[] { -0.660 } };
		expectedE = new double[][] { new double[] { 0.839 },
				new double[] { 2.875 }, new double[] { 0.120 } };
		actualXE = Extra.overconstrainedSV(a, b);
		actualX = actualXE[0];
		actualE = actualXE[1];

		assertEquals(expectedX.length, actualX.length);
		for (int row = 0; row < expectedX.length; row++) {
			assertEquals(expectedX[row].length, actualX[row].length);
			for (int col = 0; col < expectedX[row].length; col++)
				assertEquals(expectedX[row][col], actualX[row][col], DELTA);
		}

		assertEquals(expectedE.length, actualE.length);
		for (int row = 0; row < expectedE.length; row++) {
			assertEquals(expectedE[row].length, actualE[row].length);
			for (int col = 0; col < expectedE[row].length; col++)
				assertEquals(expectedE[row][col], actualE[row][col], DELTA);
		}

		actualAX = Matrices.multiply(a, actualX);
		assertEquals(expectedE.length, actualAX.length);
		assertEquals(expectedE.length, b.length);
		for (int row = 0; row < expectedE.length; row++) {
			assertEquals(expectedE[row].length, actualAX[row].length);
			assertEquals(expectedE[row].length, b[row].length);
			for (int col = 0; col < expectedE[row].length; col++)
				assertEquals(expectedE[row][col], actualAX[row][col]
						- b[row][col], DELTA);
		}

		// second test on 4 by 3 and 4 by 1 arrays
		a = new double[][] { new double[] { -3, 4, 3 },
				new double[] { -3, -2, -4 }, new double[] { 4, 2, -1 },
				new double[] { 3, 1, -3 } };
		b = new double[][] { new double[] { 2 }, new double[] { -3 },
				new double[] { 1 }, new double[] { 1 } };
		expectedX = new double[][] { new double[] { 0.239 },
				new double[] { 0.473 }, new double[] { 0.266 } };
		expectedE = new double[][] { new double[] { -0.030 },
				new double[] { 0.273 }, new double[] { 0.637 },
				new double[] { -0.607 } };
		actualXE = Extra.overconstrainedSV(a, b);
		actualX = actualXE[0];
		actualE = actualXE[1];

		assertEquals(expectedX.length, actualX.length);
		for (int row = 0; row < expectedX.length; row++) {
			assertEquals(expectedX[row].length, actualX[row].length);
			for (int col = 0; col < expectedX[row].length; col++)
				assertEquals(expectedX[row][col], actualX[row][col], DELTA);
		}

		assertEquals(expectedE.length, actualE.length);
		for (int row = 0; row < expectedE.length; row++) {
			assertEquals(expectedE[row].length, actualE[row].length);
			for (int col = 0; col < expectedE[row].length; col++)
				assertEquals(expectedE[row][col], actualE[row][col], DELTA);
		}

		actualAX = Matrices.multiply(a, actualX);
		assertEquals(expectedE.length, actualAX.length);
		assertEquals(expectedE.length, b.length);
		for (int row = 0; row < expectedE.length; row++) {
			assertEquals(expectedE[row].length, actualAX[row].length);
			assertEquals(expectedE[row].length, b[row].length);
			for (int col = 0; col < expectedE[row].length; col++)
				assertEquals(expectedE[row][col], actualAX[row][col]
						- b[row][col], DELTA);
		}

		// third test on 5 by 3 and 5 by 1 arrays
		a = new double[][] { new double[] { 3, -1, -1 },
				new double[] { -3, 2, -3 }, new double[] { 2, -4, -2 },
				new double[] { -4, -3, -1 }, new double[] { 3, 1, -1 } };
		b = new double[][] { new double[] { 3 }, new double[] { 3 },
				new double[] { 1 }, new double[] { 1 }, new double[] { 4 } };
		expectedX = new double[][] { new double[] { 0.307 },
				new double[] { 0.232 }, new double[] { -1.318 } };
		expectedE = new double[][] { new double[] { -0.994 },
				new double[] { 0.497 }, new double[] { 1.319 },
				new double[] { -1.606 }, new double[] { -1.530 } };
		actualXE = Extra.overconstrainedSV(a, b);
		actualX = actualXE[0];
		actualE = actualXE[1];

		assertEquals(expectedX.length, actualX.length);
		for (int row = 0; row < expectedX.length; row++) {
			assertEquals(expectedX[row].length, actualX[row].length);
			for (int col = 0; col < expectedX[row].length; col++)
				assertEquals(expectedX[row][col], actualX[row][col], DELTA);
		}

		assertEquals(expectedE.length, actualE.length);
		for (int row = 0; row < expectedE.length; row++) {
			assertEquals(expectedE[row].length, actualE[row].length);
			for (int col = 0; col < expectedE[row].length; col++)
				assertEquals(expectedE[row][col], actualE[row][col], DELTA);
		}

		actualAX = Matrices.multiply(a, actualX);
		assertEquals(expectedE.length, actualAX.length);
		assertEquals(expectedE.length, b.length);
		for (int row = 0; row < expectedE.length; row++) {
			assertEquals(expectedE[row].length, actualAX[row].length);
			assertEquals(expectedE[row].length, b[row].length);
			for (int col = 0; col < expectedE[row].length; col++)
				assertEquals(expectedE[row][col], actualAX[row][col]
						- b[row][col], DELTA);
		}

		try {
			// matrix a is null
			testXE = Extra.overconstrainedSV(null, b);
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testXE);
		}

		try {
			// matrix b is null
			testXE = Extra.overconstrainedSV(a, null);
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testXE);
		}

		try {
			// matrix a and b do not have equal number of rows
			testXE = Extra.overconstrainedSV(
					new double[][] { new double[] {} }, new double[][] {});
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testXE);
		}

		try {
			// matrix a does contain null column
			testXE = Extra.overconstrainedSV(new double[][] { null },
					new double[][] { null });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testXE);
		}

		try {
			// matrix b does contain null column
			testXE = Extra
					.overconstrainedSV(new double[][] { new double[] {} },
							new double[][] { null });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testXE);
		}

		try {
			// matrix a must does not have more rows than columns
			testXE = Extra.overconstrainedSV(new double[][] {
					new double[] { 0, 0 }, new double[] { 0, 0 } },
					new double[][] { new double[] { 0 }, new double[] { 0 } });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testXE);
		}

		try {
			// matrix a row width does not remain constant
			testXE = Extra.overconstrainedSV(new double[][] {
					new double[] { 0, 0 }, new double[] { 0, 0 },
					new double[] { 0, 0, 0 } },
					new double[][] { new double[] { 0 }, new double[] { 0 },
							new double[] { 0 } });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testXE);
		}

		try {
			// matrix b row width does not remain constant
			testXE = Extra.overconstrainedSV(new double[][] {
					new double[] { 0, 0 }, new double[] { 0, 0 },
					new double[] { 0, 0 } }, new double[][] {
					new double[] { 0 }, new double[] { 0 },
					new double[] { 0, 0 } });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testXE);
		}
	}

	/**
	 * Test method for {@link Extra#underconstrainedSV(double[][], double[][])}.
	 */
	@Test
	public void testUnderconstrainedSV() {
		double[][] a;
		double[][] b;
		double[][] actualAX;
		double[][] expectedX;
		double[][] actualX;
		double[][] testX = null;

		// first test on 2 by 3 and 2 by 1 arrays
		a = new double[][] { new double[] { 4, -3, -1 },
				new double[] { 1, -4, -3 } };
		b = new double[][] { new double[] { -1 }, new double[] { -3 } };
		expectedX = new double[][] { new double[] { 0.206 },
				new double[] { 0.454 }, new double[] { 0.463 } };
		actualX = Extra.underconstrainedSV(a, b);

		assertEquals(expectedX.length, actualX.length);
		for (int row = 0; row < expectedX.length; row++) {
			assertEquals(expectedX[row].length, actualX[row].length);
			for (int col = 0; col < expectedX[row].length; col++)
				assertEquals(expectedX[row][col], actualX[row][col], DELTA);
		}

		actualAX = Matrices.multiply(a, actualX);
		assertEquals(b.length, actualAX.length);
		assertEquals(b.length, b.length);
		for (int row = 0; row < b.length; row++) {
			assertEquals(b[row].length, actualAX[row].length);
			assertEquals(b[row].length, b[row].length);
			for (int col = 0; col < b[row].length; col++)
				assertEquals(b[row][col], actualAX[row][col], DELTA);
		}

		// second test on 3 by 4 and 3 by 1 arrays
		a = new double[][] { new double[] { 4, -3, -1, 1 },
				new double[] { -4, -3, -1, -3 }, new double[] { 4, 2, 3, 3 } };
		b = new double[][] { new double[] { 3 }, new double[] { 1 },
				new double[] { -1 } };
		expectedX = new double[][] { new double[] { 0.279 },
				new double[] { -0.555 }, new double[] { -0.277 },
				new double[] { -0.058 } };
		actualX = Extra.underconstrainedSV(a, b);

		assertEquals(expectedX.length, actualX.length);
		for (int row = 0; row < expectedX.length; row++) {
			assertEquals(expectedX[row].length, actualX[row].length);
			for (int col = 0; col < expectedX[row].length; col++)
				assertEquals(expectedX[row][col], actualX[row][col], DELTA);
		}

		actualAX = Matrices.multiply(a, actualX);
		assertEquals(b.length, actualAX.length);
		assertEquals(b.length, b.length);
		for (int row = 0; row < b.length; row++) {
			assertEquals(b[row].length, actualAX[row].length);
			assertEquals(b[row].length, b[row].length);
			for (int col = 0; col < b[row].length; col++)
				assertEquals(b[row][col], actualAX[row][col], DELTA);
		}

		try {
			// matrix a is null
			testX = Extra.underconstrainedSV(null, b);
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testX);
		}

		try {
			// matrix b is null
			testX = Extra.underconstrainedSV(a, null);
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testX);
		}

		try {
			// matrix a and b do not have equal number of rows
			testX = Extra.underconstrainedSV(
					new double[][] { new double[] {} }, new double[][] {});
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testX);
		}

		try {
			// matrix a does contain null column
			testX = Extra.underconstrainedSV(new double[][] { null },
					new double[][] { null });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testX);
		}

		try {
			// matrix b does contain null column
			testX = Extra
					.underconstrainedSV(new double[][] { new double[] {} },
							new double[][] { null });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testX);
		}

		try {
			// matrix a must does not have more rows than columns
			testX = Extra.underconstrainedSV(new double[][] {
					new double[] { 0, 0 }, new double[] { 0, 0 } },
					new double[][] { new double[] { 0 }, new double[] { 0 } });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testX);
		}

		try {
			// matrix a row width does not remain constant
			testX = Extra.underconstrainedSV(new double[][] {
					new double[] { 0, 0, 0 }, new double[] { 0, 0, 0, 0 } },
					new double[][] { new double[] { 0 }, new double[] { 0 } });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testX);
		}

		try {
			// matrix b row width does not remain constant
			testX = Extra
					.underconstrainedSV(
							new double[][] { new double[] { 0, 0, 0 },
									new double[] { 0, 0, 0 } }, new double[][] {
									new double[] { 0 }, new double[] { 0, 0 } });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testX);
		}
	}

}