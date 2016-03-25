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
	 * Test method for {@link Extra#overconstrainted(double[][], double[][])}.
	 */
	@Test
	public void testOverconstrainted() {
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
		actualXE = Extra.overconstrainted(a, b);
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
		actualXE = Extra.overconstrainted(a, b);
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
		actualXE = Extra.overconstrainted(a, b);
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
			testXE = Extra.overconstrainted(null, b);
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testXE);
		}

		try {
			// matrix b is null
			testXE = Extra.overconstrainted(a, null);
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testXE);
		}

		try {
			// matrix a and b do not have equal number of rows
			testXE = Extra.overconstrainted(new double[][] { new double[] {} },
					new double[][] {});
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testXE);
		}

		try {
			// matrix a does contain null column
			testXE = Extra.overconstrainted(new double[][] { null },
					new double[][] { null });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testXE);
		}

		try {
			// matrix b does contain null column
			testXE = Extra.overconstrainted(new double[][] { new double[] {} },
					new double[][] { null });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testXE);
		}

		try {
			// matrix a must does not have more rows than columns
			testXE = Extra.overconstrainted(new double[][] {
					new double[] { 0, 0 }, new double[] { 0, 0 } },
					new double[][] { new double[] { 0 }, new double[] { 0 } });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testXE);
		}

		try {
			// matrix a row width does not remain constant
			testXE = Extra.overconstrainted(new double[][] {
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
			testXE = Extra.overconstrainted(new double[][] {
					new double[] { 0, 0 }, new double[] { 0, 0 },
					new double[] { 0, 0 } }, new double[][] {
					new double[] { 0 }, new double[] { 0 },
					new double[] { 0, 0 } });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testXE);
		}
	}

}