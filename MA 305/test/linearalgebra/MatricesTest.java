package linearalgebra;

import static org.junit.Assert.*;
import linearalgebra.Matrices;

import org.junit.Test;

/**
 * Matrices test.
 * 
 * @author Jacob Malter
 *
 */
public class MatricesTest {

	/** acceptable margin of error */
	private static final double DELTA = 0.0001d;

	/**
	 * Test method for {@link Matrices#determinant(double[][])}.
	 */
	@Test
	public void testDeterminant() {
		// variables required for testing
		double[][] a;
		double expectedDet;
		double actualDet;
		double testDet = 0;

		// first test on 0 component matrix
		a = new double[][] {};
		expectedDet = 0;
		actualDet = Matrices.determinant(a);

		assertEquals(expectedDet, actualDet, DELTA);

		// first test on 1 component matrix
		a = new double[][] { new double[] { 3 } };
		expectedDet = 3;
		actualDet = Matrices.determinant(a);

		assertEquals(expectedDet, actualDet, DELTA);

		// second test on 2 component matrix
		a = new double[][] { new double[] { 3, 8 }, new double[] { 4, 6 } };
		expectedDet = -14;
		actualDet = Matrices.determinant(a);

		assertEquals(expectedDet, actualDet, DELTA);

		// third test on 3 component matrix
		a = new double[][] { new double[] { 6, 1, 1 }, new double[] { 4, -2, 5 }, new double[] { 2, 8, 7 } };
		expectedDet = -306;
		actualDet = Matrices.determinant(a);

		assertEquals(expectedDet, actualDet, DELTA);

		// fourth test on 4 component matrix
		a = new double[][] { new double[] { 0, 2, 0, 0 }, new double[] { 6, 0, 1, 1 }, new double[] { 4, 0, -2, 5 },
				new double[] { 2, 0, 8, 7 } };
		expectedDet = 612;
		actualDet = Matrices.determinant(a);

		assertEquals(expectedDet, actualDet, DELTA);

		// fifth test on 4 component matrix
		a = new double[][] { new double[] { 2, 5, 3, 5 }, new double[] { 4, 6, 6, 3 }, new double[] { 11, 3, 2, -2 },
				new double[] { 4, -7, 9, 3 } };
		expectedDet = 2960;
		actualDet = Matrices.determinant(a);

		assertEquals(expectedDet, actualDet, DELTA);

		// sixth test on 5 component matrix
		a = new double[][] { new double[] { 276, 1, 179, 23, 9387 }, new double[] { 0, 0, 78, 0, 0 },
				new double[] { 0, 0, -1, 0, 1 }, new double[] { 0, 0, 1994, -1, 1089 },
				new double[] { 1, 0, 212, 726, -378 } };
		expectedDet = 78;
		actualDet = Matrices.determinant(a);

		assertEquals(expectedDet, actualDet, DELTA);

		try {
			// array a is null
			testDet = Matrices.determinant(null);
			fail();
		} catch (IllegalArgumentException ex) {
			assertEquals(0, testDet, DELTA);
		}

		try {
			// matrix a does contain null column
			testDet = Matrices.determinant(new double[][] { null });
			fail();
		} catch (IllegalArgumentException ex) {
			assertEquals(0, testDet, DELTA);
		}

		try {
			// matrix a must does not have as many columns as rows
			testDet = Matrices.determinant(
					new double[][] { new double[] { 0, 0 }, new double[] { 0, 0 }, new double[] { 0, 0 } });
			fail();
		} catch (IllegalArgumentException ex) {
			assertEquals(0, testDet, DELTA);
		}
	}

	/**
	 * Test method for {@link Matrices#dotProduct(double[], double[])}.
	 */
	@Test
	public void testDotProduct() {
		// variables required for testing
		double[] a;
		double[] b;
		double expectedDot;
		double actualDot;
		double testDot = 0;

		// first test on 3 component arrays
		a = new double[] { 1, 2, 3 };
		b = new double[] { 7, 9, 11 };
		expectedDot = 58;
		actualDot = Matrices.dotProduct(a, b);

		assertEquals(expectedDot, actualDot, DELTA);

		// second test on 2 component arrays
		a = new double[] { -4, -9 };
		b = new double[] { -1, 2 };
		expectedDot = -14;
		actualDot = Matrices.dotProduct(a, b);

		assertEquals(expectedDot, actualDot, DELTA);

		try {
			// array a is null
			testDot = Matrices.dotProduct(null, new double[] {});
			fail();
		} catch (IllegalArgumentException ex) {
			assertEquals(0, testDot, DELTA);
		}

		try {
			// array b is null
			testDot = Matrices.dotProduct(new double[] {}, null);
			fail();
		} catch (IllegalArgumentException ex) {
			assertEquals(0, testDot, DELTA);
		}

		try {
			// array a and b do not have equal number of components
			testDot = Matrices.dotProduct(new double[] { 0 }, new double[] {});
			fail();
		} catch (IllegalArgumentException ex) {
			assertEquals(0, testDot, DELTA);
		}
	}

	/**
	 * Test method for {@link Matrices#minor(double[][], int, int)}
	 */
	@Test
	public void testMinor() {
		// variables required for testing
		double[][] a;
		double[][] expectedMin;
		double[][] actualMin;
		double[][] testMin = null;

		// first test on 1 by 1 array
		a = new double[][] { new double[] { 0 } };
		expectedMin = new double[][] {};
		actualMin = Matrices.minor(a, 0, 0);

		assertEquals(expectedMin.length, actualMin.length);
		for (int row = 0; row < expectedMin.length; row++) {
			assertEquals(expectedMin[row].length, actualMin[row].length);
			for (int col = 0; col < expectedMin[row].length; col++)
				assertEquals(expectedMin[row][col], actualMin[row][col], DELTA);
		}

		// second test on 2 by 2 array
		a = new double[][] { new double[] { 0, 1 }, new double[] { 2, 3 } };
		expectedMin = new double[][] { new double[] { 1 } };
		actualMin = Matrices.minor(a, 1, 0);

		assertEquals(expectedMin.length, actualMin.length);
		for (int row = 0; row < expectedMin.length; row++) {
			assertEquals(expectedMin[row].length, actualMin[row].length);
			for (int col = 0; col < expectedMin[row].length; col++)
				assertEquals(expectedMin[row][col], actualMin[row][col], DELTA);
		}

		// third test on 3 by 3 array
		a = new double[][] { new double[] { 3, 0, 2 }, new double[] { 2, 0, -2 }, new double[] { 0, 1, 1 } };
		expectedMin = new double[][] { new double[] { 3, 2 }, new double[] { 0, 1 } };
		actualMin = Matrices.minor(a, 1, 1);

		assertEquals(expectedMin.length, actualMin.length);
		for (int row = 0; row < expectedMin.length; row++) {
			assertEquals(expectedMin[row].length, actualMin[row].length);
			for (int col = 0; col < expectedMin[row].length; col++)
				assertEquals(expectedMin[row][col], actualMin[row][col], DELTA);
		}

		try {
			// matrix a is null
			testMin = Matrices.minor(null, 0, 0);
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testMin);
		}

		try {
			// skipRow is less than 0
			testMin = Matrices.minor(new double[][] {}, -1, 0);
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testMin);
		}

		try {
			// skipCol is less than 0
			testMin = Matrices.minor(new double[][] {}, 0, -1);
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testMin);
		}

		try {
			// skipRow is equal to number of rows
			testMin = Matrices.minor(new double[][] {}, 0, 0);
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testMin);
		}

		try {
			// matrix a does contain null column
			testMin = Matrices.minor(new double[][] { null }, 0, 0);
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testMin);
		}

		try {
			// matrix a row width does not remain constant
			testMin = Matrices.minor(new double[][] { new double[] { 0, 0 }, new double[] { 0 } }, 0, 0);
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testMin);
		}

		try {
			// skipCol is equal to number of columns
			testMin = Matrices.minor(new double[][] { new double[] { 0, 0 }, new double[] { 0, 0 } }, 0, 2);
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testMin);
		}
	}

	/**
	 * Test method for {@link Matrices#multiply(double[][], double[][])}.
	 */
	@Test
	public void testMultiply() {
		// variables required for testing
		double[][] a;
		double[][] b;
		double[][] expectedMul;
		double[][] actualMul;
		double[][] testMul = null;

		// first test on 2 by 3 and 3 by 2 arrays
		a = new double[][] { new double[] { 1, 2, 3 }, new double[] { 4, 5, 6 } };
		b = new double[][] { new double[] { 7, 8 }, new double[] { 9, 10 }, new double[] { 11, 12 } };
		expectedMul = new double[][] { new double[] { 58, 64 }, new double[] { 139, 154 } };
		actualMul = Matrices.multiply(a, b);

		assertEquals(expectedMul.length, actualMul.length);
		for (int row = 0; row < expectedMul.length; row++) {
			assertEquals(expectedMul[row].length, actualMul[row].length);
			for (int col = 0; col < expectedMul[row].length; col++)
				assertEquals(expectedMul[row][col], actualMul[row][col], DELTA);
		}

		// first test on 1 by 3 and 3 by 4 arrays
		a = new double[][] { new double[] { 3, 4, 2 } };
		b = new double[][] { new double[] { 13, 9, 7, 15 }, new double[] { 8, 7, 4, 6 }, new double[] { 6, 4, 0, 3 } };
		expectedMul = new double[][] { new double[] { 83, 63, 37, 75 } };
		actualMul = Matrices.multiply(a, b);

		assertEquals(expectedMul.length, actualMul.length);
		for (int row = 0; row < expectedMul.length; row++) {
			assertEquals(expectedMul[row].length, actualMul[row].length);
			for (int col = 0; col < expectedMul[row].length; col++)
				assertEquals(expectedMul[row][col], actualMul[row][col], DELTA);
		}

		try {
			// matrix a is null
			testMul = Matrices.multiply(null, new double[][] {});
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testMul);
		}

		try {
			// matrix b is null
			testMul = Matrices.multiply(new double[][] {}, null);
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testMul);
		}

		try {
			// matrix a does contain null column
			testMul = Matrices.multiply(new double[][] { null }, new double[][] { null });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testMul);
		}

		try {
			// matrix b does contain null column
			testMul = Matrices.multiply(new double[][] { new double[] {} }, new double[][] { null });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testMul);
		}

		try {
			// matrix a row width does not remain constant
			testMul = Matrices.multiply(new double[][] { new double[] { 0, 0 }, new double[] { 0 } },
					new double[][] { new double[] {}, new double[] {} });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testMul);
		}

		try {
			// matrix b number of rows does not equal matrix a number of columns
			testMul = Matrices.multiply(new double[][] { new double[] { 0, 0, 0 }, new double[] { 0, 0, 0 } },
					new double[][] { new double[] { 0 }, new double[] { 0 } });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testMul);
		}
	}

	/**
	 * Test method for {@link Matrices#transpose(double[][])}.
	 */
	@Test
	public void testTranspose() {
		// variables required for testing
		double[][] a;
		double[][] expectedAT;
		double[][] actualAT;
		double[][] testAT = null;

		// first test on 1 by 1 array
		a = new double[][] { new double[] { 0 } };
		expectedAT = new double[][] { new double[] { 0 } };
		actualAT = Matrices.transpose(a);

		assertEquals(expectedAT.length, actualAT.length);
		for (int row = 0; row < expectedAT.length; row++) {
			assertEquals(expectedAT[row].length, actualAT[row].length);
			for (int col = 0; col < expectedAT[row].length; col++)
				assertEquals(expectedAT[row][col], actualAT[row][col], DELTA);
		}

		// second test on 2 by 2 array
		a = new double[][] { new double[] { 0, 1 }, new double[] { 2, 3 } };
		expectedAT = new double[][] { new double[] { 0, 2 }, new double[] { 1, 3 } };
		actualAT = Matrices.transpose(a);

		assertEquals(expectedAT.length, actualAT.length);
		for (int row = 0; row < expectedAT.length; row++) {
			assertEquals(expectedAT[row].length, actualAT[row].length);
			for (int col = 0; col < expectedAT[row].length; col++)
				assertEquals(expectedAT[row][col], actualAT[row][col], DELTA);
		}

		// third test on 3 by 3 array
		a = new double[][] { new double[] { 3, 0, 2 }, new double[] { 2, 0, -2 }, new double[] { 0, 1, 1 } };
		expectedAT = new double[][] { new double[] { 3, 2, 0 }, new double[] { 0, 0, 1 }, new double[] { 2, -2, 1 } };
		actualAT = Matrices.transpose(a);

		assertEquals(expectedAT.length, actualAT.length);
		for (int row = 0; row < expectedAT.length; row++) {
			assertEquals(expectedAT[row].length, actualAT[row].length);
			for (int col = 0; col < expectedAT[row].length; col++)
				assertEquals(expectedAT[row][col], actualAT[row][col], DELTA);
		}

		try {
			// matrix a is null
			testAT = Matrices.transpose(null);
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testAT);
		}

		try {
			// matrix a does contain null column
			testAT = Matrices.transpose(new double[][] { null });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testAT);
		}

		try {
			// matrix a row width does not remain constant
			testAT = Matrices.transpose(new double[][] { new double[] { 0, 0 }, new double[] { 0 } });
			fail();
		} catch (IllegalArgumentException ex) {
			assertNull(testAT);
		}
	}

	/**
	 * Test method for {@link Matrices#sofu(double[][])}.
	 */
	@Test
	public void testSofu() {
		// variables required for testing
		double[][] a;
		double expectedSofu;
		double actualSofu;
		double testSofu = 0;

		// first test on 2 by 1 array
		a = new double[][] { new double[] { 3 }, new double[] { -4 } };
		expectedSofu = 5;
		actualSofu = Matrices.sofu(a);

		assertEquals(expectedSofu, actualSofu, DELTA);

		// second test on 3 by 2 array
		a = new double[][] { new double[] { -2, -4 }, new double[] { 3, 1 }, new double[] { 2, -1 } };
		expectedSofu = 15;
		actualSofu = Matrices.sofu(a);

		assertEquals(expectedSofu, actualSofu, DELTA);

		// third test on 4 by 5 array
		a = new double[][] { new double[] { -2, 2, -4, 1 }, new double[] { -1, -4, -1, -1 },
				new double[] { 5, 4, -4, 4 }, new double[] { -2, 1, -2, 1 }, new double[] { -3, 3, 5, -5 } };
		expectedSofu = 526.49786324352732362438329844498;
		actualSofu = Matrices.sofu(a);

		assertEquals(expectedSofu, actualSofu, DELTA);

		try {
			// matrix a is null
			testSofu = Matrices.sofu(null);
			fail();
		} catch (IllegalArgumentException ex) {
			assertEquals(0, testSofu, DELTA);
		}

		try {
			// matrix a does contain null column
			testSofu = Matrices.sofu(new double[][] { null });
			fail();
		} catch (IllegalArgumentException ex) {
			assertEquals(0, testSofu, DELTA);
		}

		try {
			// matrix a row width does not remain constant
			testSofu = Matrices.sofu(new double[][] { new double[] { 0, 0 }, new double[] { 0 } });
			fail();
		} catch (IllegalArgumentException ex) {
			assertEquals(0, testSofu, DELTA);
		}
	}

}