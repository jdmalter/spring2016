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
   public void testLU() {
      // variables required for testing
      double[][] a;
      double[][] expectedL;
      double[][] expectedU;
      double[][][] test = null;

      // first test on 2 by 2 array
      a = new double[][] { { 2, -1 }, { -4, -1 } };
      expectedL = new double[][] { { 1, 0 }, { -2, 1 } };
      expectedU = new double[][] { { 2, -1 }, { 0, -3 } };
      testcaseLU(a, expectedL, expectedU);

      // second test on 3 by 3 array
      a = new double[][] { { 2, 4, -1 }, { -6, -15, 5 }, { 8, 22, -10 } };
      expectedL = new double[][] { { 1, 0, 0 }, { -3, 1, 0, }, { 4, -2, 1 } };
      expectedU = new double[][] { { 2, 4, -1 }, { 0, -3, 2 }, { 0, 0, -2 } };
      testcaseLU(a, expectedL, expectedU);

      // third test on 4 by 4 array
      a = new double[][] { { 2, -1, 2, -4 }, { 4, -5, 5, -12 }, { 2, 8, -5, 5 },
         { -8, 16, 4, 40 } };
      expectedL = new double[][] { { 1, 0, 0, 0 }, { 2, 1, 0, 0 },
         { 1, -3, 1, 0 }, { -4, -4, -4, 1 } };
      expectedU = new double[][] { { 2, -1, 2, -4 }, { 0, -3, 1, -4 },
         { 0, 0, -4, -3 }, { 0, 0, 0, -4 } };
      testcaseLU(a, expectedL, expectedU);

      try {
         // matrix a is null
         test = Decomposition.lu(null);
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix a does contain null column
         test = Decomposition.lu(new double[][] { null });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix a must does not have as many columns as rows
         test =
            Decomposition.lu(new double[][] { { 0, 0 }, { 0, 0 }, { 0, 0 } });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }
   }

   private void testcaseLU(double[][] a, double[][] expectedL,
      double[][] expectedU) {
      double[][][] actualLU = Decomposition.lu(a);
      double[][] actual = actualLU[0];

      assertEquals(expectedL.length, actual.length);
      for (int row = 0; row < expectedL.length; row++) {
         assertEquals(expectedL[row].length, actual[row].length);
         for (int col = 0; col < expectedL[row].length; col++)
            assertEquals(expectedL[row][col], actual[row][col], DELTA);
      }

      actual = actualLU[1];
      assertEquals(expectedU.length, actual.length);
      for (int row = 0; row < expectedU.length; row++) {
         assertEquals(expectedU[row].length, actual[row].length);
         for (int col = 0; col < expectedU[row].length; col++)
            assertEquals(expectedU[row][col], actual[row][col], DELTA);
      }

      actual = Matrices.multiply(actualLU[0], actualLU[1]);
      assertEquals(a.length, actual.length);
      for (int row = 0; row < a.length; row++) {
         assertEquals(a[row].length, actual[row].length);
         for (int col = 0; col < a[row].length; col++)
            assertEquals(a[row][col], actual[row][col], DELTA);
      }
   }

   /**
    * Test method for {@link Decomposition#cholesky(double[][])}.
    */
   @Test
   public void testCholesky() {
      // variables required for testing
      double[][] a;
      double[][] expected;
      double[][] test = null;

      // first test on 2 by 2 array
      a = new double[][] { { 16, 8 }, { 8, 20 } };
      expected = new double[][] { { 4, 0 }, { 2, 4 } };
      testcaseCholesky(a, expected);

      // second test on 3 by 3 array
      a = new double[][] { { 4, 4, -6 }, { 4, 20, 10 }, { -6, 10, 29 } };
      expected = new double[][] { { 2, 0, 0 }, { 2, 4, 0 }, { -3, 4, 2 } };
      testcaseCholesky(a, expected);

      // third test on 4 by 4 array
      a = new double[][] { { 4, 8, -6, 6 }, { 8, 17, -14, 14 },
         { -6, -14, 17, -5 }, { 6, 14, -5, 38 } };
      expected = new double[][] { { 2, 0, 0, 0 }, { 4, 1, 0, 0 },
         { -3, -2, 2, 0 }, { 3, 2, 4, 3 } };
      testcaseCholesky(a, expected);

      try {
         // matrix a is null
         test = Decomposition.cholesky(null);
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix a does contain null column
         test = Decomposition.cholesky(new double[][] { null });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix a must does not have as many columns as rows
         test = Decomposition
            .cholesky(new double[][] { { 0, 0 }, { 0, 0 }, { 0, 0 } });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix a must does not have symmetry across diagonal
         test = Decomposition.cholesky(new double[][] { { 1, 2 }, { 3, 4 } });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }
   }

   private void testcaseCholesky(double[][] a, double[][] expected) {
      double[][] actual = Decomposition.cholesky(a);

      assertEquals(expected.length, actual.length);
      for (int row = 0; row < expected.length; row++) {
         assertEquals(expected[row].length, actual[row].length);
         for (int col = 0; col < expected[row].length; col++)
            assertEquals(expected[row][col], actual[row][col], DELTA);
      }

      actual = Matrices.multiply(actual, Matrices.transpose(actual));
      assertEquals(a.length, actual.length);
      for (int row = 0; row < a.length; row++) {
         assertEquals(a[row].length, actual[row].length);
         for (int col = 0; col < a[row].length; col++)
            assertEquals(a[row][col], actual[row][col], DELTA);
      }
   }

   /**
    * Test method for {@link Decomposition#qr(double[][])}.
    */
   @Test
   public void testQR() {
      // variables required for testing
      double[][] a;
      double[][] expectedQ;
      double[][] expectedR;
      double[][][] test = null;

      // first test on 2 by 2 array
      a = new double[][] { { 14, -6 }, { 14, -9 }, { 7, -6 } };
      expectedQ = new double[][] { { 2d / 3, 2d / 3 }, { 2d / 3, -1d / 3 },
         { 1d / 3, -2d / 3 } };
      expectedR = new double[][] { { 21, -12 }, { 0, 3 } };
      testcaseQR(a, expectedQ, expectedR);

      // second test on 4 by 3 array
      a = new double[][] { new double[] { -3, -8, -9 },
         new double[] { -3, -10, 15 }, new double[] { -3, -8, 7 },
         new double[] { -3, -10, -1 } };
      expectedQ = new double[][] { new double[] { -0.5, 0.5, -0.5 },
         new double[] { -0.5, -0.5, 0.5 }, new double[] { -0.5, 0.5, 0.5 },
         new double[] { -0.5, -0.5, -0.5 } };
      expectedR = new double[][] { new double[] { 6, 18, -6 },
         new double[] { 0, 2, -8 }, new double[] { 0, 0, 16 } };
      testcaseQR(a, expectedQ, expectedR);

      try {
         // matrix a is null
         test = Decomposition.qr(null);
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix a does contain null column
         test = Decomposition.qr(new double[][] { null });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix a contains as fewer rows than columns
         test = Decomposition.qr(new double[][] { { 0, 0, 0 }, { 0, 0, 0 } });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix a row width does not remain constant
         test = Decomposition.qr(new double[][] { { 0, 0 }, { 0 } });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }
   }

   private void testcaseQR(double[][] a, double[][] expectedQ,
      double[][] expectedR) {
      double[][][] actualQR = Decomposition.qr(a);
      double[][] actual = actualQR[0];

      assertEquals(expectedQ.length, actual.length);
      for (int row = 0; row < expectedQ.length; row++) {
         assertEquals(expectedQ[row].length, actual[row].length);
         for (int col = 0; col < expectedQ[row].length; col++)
            assertEquals(expectedQ[row][col], actual[row][col], DELTA);
      }

      actual = actualQR[1];
      assertEquals(expectedR.length, actual.length);
      for (int row = 0; row < expectedR.length; row++) {
         assertEquals(expectedR[row].length, actual[row].length);
         for (int col = 0; col < expectedR[row].length; col++)
            assertEquals(expectedR[row][col], actual[row][col], DELTA);
      }

      actual = Matrices.multiply(actualQR[0], actualQR[1]);
      assertEquals(a.length, actual.length);
      for (int row = 0; row < a.length; row++) {
         assertEquals(a[row].length, actual[row].length);
         for (int col = 0; col < a[row].length; col++)
            assertEquals(a[row][col], actual[row][col], DELTA);
      }

      actual = Matrices.multiply(Matrices.transpose(actualQR[0]), actualQR[0]);
      for (int row = 0; row < actual.length; row++)
         for (int col = 0; col < actual[row].length; col++)
            if (row == col) assertEquals(1, actual[row][col], DELTA);
            else assertEquals(0, actual[row][col], DELTA);
   }

   /**
    * Test method for {@link Decomposition#eigen(double[][])}.
    */
   @Test
   public void testEigen() {
      // variables required for testing
      double[][] a;
      double[][] expectedL;
      double[][] expectedV;
      double[][][] test = null;

      // first test on 1 by 1 array
      a = new double[][] { { 7 } };
      expectedL = new double[][] { { 7 } };
      expectedV = new double[][] { { 1 } };
      testcaseEigen(a, expectedL, expectedV);

      // second test on 2 by 2 array
      a = new double[][] { { -1, 6 }, { -2, 6 } };
      expectedL = new double[][] { { 3, 0 }, { 0, 2 } };
      expectedV = new double[][] { { 3 / Math.sqrt(13d), 2 / Math.sqrt(5d) },
         { 2 / Math.sqrt(13d), 1 / Math.sqrt(5d) } };
      testcaseEigen(a, expectedL, expectedV);

      // third test on 2 by 2 array
      a = new double[][] { { -2, -2 }, { 2, 3 } };
      expectedL = new double[][] { { 2, 0 }, { 0, -1 } };
      expectedV = new double[][] { { -1 / Math.sqrt(5), -2 / Math.sqrt(5) },
         { 2 / Math.sqrt(5), 1 / Math.sqrt(5) } };
      testcaseEigen(a, expectedL, expectedV);

      // fourth test on 3 by 3 array
      a = new double[][] { { 2, 3, -3 }, { -2, -1, 2 }, { 2, 4, -3 } };
      expectedL = new double[][] { { 1, 0, 0 }, { 0, -1, 0 }, { 0, 0, -2 } };
      expectedV = new double[][] { { 0, 1 / Math.sqrt(2), 9 / Math.sqrt(185) },
         { 1 / Math.sqrt(2), 0, -2 / Math.sqrt(185) },
         { 1 / Math.sqrt(2), 1 / Math.sqrt(2), 10 / Math.sqrt(185) } };
      testcaseEigen(a, expectedL, expectedV);

      // fifth test on 3 by 3 array
      a = new double[][] { { 3, 2, 2 }, { -1, -4, 2 }, { -2, -2, -1 } };
      expectedL = new double[][] { { 1, 0, 0 }, { 0, -1, 0 }, { 0, 0, -2 } };
      expectedV = new double[][] {
         { -7 / Math.sqrt(74), -1 / Math.sqrt(3), -2 / Math.sqrt(17) },
         { 3 / Math.sqrt(74), 1 / Math.sqrt(3), 3 / Math.sqrt(17) },
         { 4 / Math.sqrt(74), 1 / Math.sqrt(3), 2 / Math.sqrt(17) } };
      testcaseEigen(a, expectedL, expectedV);

      try {
         // matrix a is null
         test = Decomposition.eigen(null);
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix a does contain null column
         test = Decomposition.eigen(new double[][] { null });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix a must does not have as many columns as rows
         test = Decomposition
            .eigen(new double[][] { { 0, 0 }, { 0, 0 }, { 0, 0 } });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix a is too large
         test = Decomposition.eigen(new double[][] { { 0, 0, 0, 0 },
            { 0, 0, 0, 0 }, { 0, 0, 0, 0 }, { 0, 0, 0, 0 } });
         fail();
      } catch (UnsupportedOperationException ex) {
         assertNull(test);
      }
   }

   private void testcaseEigen(double[][] a, double[][] expectedL,
      double[][] expectedV) {
      double[][][] actualEigen = Decomposition.eigen(a);

      assertEquals(actualEigen[0].length, expectedL.length);
      for (int row = 0; row < expectedL.length; row++) {
         assertEquals(actualEigen[0][row].length, expectedL[row].length);
         for (int col = 0; col < expectedL[row].length; col++)
            assertEquals(actualEigen[0][row][col], expectedL[row][col], DELTA);
      }

      assertEquals(actualEigen[1].length, expectedV.length);
      for (int row = 0; row < expectedV.length; row++) {
         assertEquals(actualEigen[1][row].length, expectedV[row].length);
         for (int col = 0; col < expectedV[row].length; col++)
            assertEquals(actualEigen[1][row][col], expectedV[row][col], DELTA);
      }

      double[][] actualAV = Matrices.multiply(a, actualEigen[1]);
      double[][] actualVL = Matrices.multiply(actualEigen[1], actualEigen[0]);
      assertEquals(actualAV.length, actualVL.length);
      for (int row = 0; row < actualAV.length; row++) {
         assertEquals(actualAV[row].length, actualVL[row].length);
         for (int col = 0; col < actualAV[row].length; col++)
            assertEquals(actualAV[row][col], actualVL[row][col], DELTA);
      }
   }

   /**
    * Test method for {@link Decomposition#singularValue(double[][])}.
    */
   @Test
   public void testSingularValue() {
      // variables required for testing
      double[][] a;
      double[][] expectedU;
      double[][] expectedS;
      double[][] expectedV;
      double[][][] testSingular = null;

      // first test on 2 by 2 array
      a = new double[][] { { 8, 6 }, { 3, -4 } };
      expectedU = new double[][] { { 1, 0 }, { 0, -1 } };
      expectedS = new double[][] { { 10, 0 }, { 0, 5 } };
      expectedV = new double[][] { { 0.8, -0.6 }, { 0.6, 0.8 } };
      testcaseSingularValue(a, expectedU, expectedS, expectedV);

      // second test on 2 by 3 array
      a = new double[][] { { 18, 26 }, { 18, 1 }, { -27, -14 } };
      expectedU = new double[][] { { 2d / 3, 2d / 3 }, { 1d / 3, -2d / 3 },
         { -2d / 3, 1d / 3 } };
      expectedS = new double[][] { { 45, 0 }, { 0, 15 } };
      expectedV = new double[][] { { 0.8, -0.6 }, { 0.6, 0.8 } };
      testcaseSingularValue(a, expectedU, expectedS, expectedV);

      // third test on 3 by 3 array
      a = new double[][] { { 6, 6, 3 }, { -1, 2, -2 }, { 4, -2, -4 } };
      expectedU = new double[][] { { 1, 0, 0 }, { 0, 0, -1 }, { 0, -1, 0 } };
      expectedS = new double[][] { { 9, 0, 0 }, { 0, 6, 0 }, { 0, 0, 3 } };
      expectedV = new double[][] { { 2d / 3, -2d / 3, 1d / 3 },
         { 2d / 3, 1d / 3, -2d / 3 }, { 1d / 3, 2d / 3, 2d / 3 } };
      testcaseSingularValue(a, expectedU, expectedS, expectedV);

      // fourth test on 4 by 2 array
      a = new double[][] { { 18, 1 }, { 18, 1 }, { 6, 17 }, { -6, -17 } };
      expectedU = new double[][] { { 0.5, -0.5 }, { 0.5, -0.5 }, { 0.5, 0.5 },
         { -0.5, -0.5 } };
      expectedS = new double[][] { { 30, 0 }, { 0, 20 } };
      expectedV = new double[][] { { 0.8, -0.6 }, { 0.6, 0.8 } };
      testcaseSingularValue(a, expectedU, expectedS, expectedV);

      // fifth test on 4 by 3 array
      a = new double[][] { { 11, 8, 2 }, { 13, 4, -2 }, { 3, 12, -6 },
         { -5, -8, 10 } };
      expectedU = new double[][] { { -0.5, 0.5, 0.5 }, { -0.5, 0.5, -0.5 },
         { -0.5, -0.5, 0.5 }, { 0.5, 0.5, 0.5 } };
      expectedS = new double[][] { { 24, 0, 0 }, { 0, 12, 0 }, { 0, 0, 6 } };
      expectedV = new double[][] { { -2d / 3, 2d / 3, -1d / 3 },
         { -2d / 3, -1d / 3, 2d / 3 }, { 1d / 3, 2d / 3, 2d / 3 } };
      testcaseSingularValue(a, expectedU, expectedS, expectedV);

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
         testSingular = Decomposition
            .singularValue(new double[][] { { 0, 0, 0 }, { 0, 0, 0 } });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(testSingular);
      }

      try {
         // matrix a row width does not remain constant
         testSingular =
            Decomposition.singularValue(new double[][] { { 0, 0 }, { 0 } });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(testSingular);
      }
   }

   private void testcaseSingularValue(double[][] a, double[][] expectedU,
      double[][] expectedS, double[][] expectedV) {
      double[][][] actualSingular = Decomposition.singularValue(a);

      assertEquals(expectedU.length, actualSingular[0].length);
      for (int row = 0; row < expectedU.length; row++) {
         assertEquals(expectedU[row].length, actualSingular[0][row].length);
         for (int col = 0; col < expectedU[row].length; col++)
            assertEquals(expectedU[row][col], actualSingular[0][row][col],
               DELTA);
      }

      assertEquals(expectedS.length, actualSingular[1].length);
      for (int row = 0; row < expectedS.length; row++) {
         assertEquals(expectedS[row].length, actualSingular[1][row].length);
         for (int col = 0; col < expectedS[row].length; col++)
            assertEquals(expectedS[row][col], actualSingular[1][row][col],
               DELTA);
      }

      assertEquals(expectedV.length, actualSingular[2].length);
      for (int row = 0; row < expectedV.length; row++) {
         assertEquals(expectedV[row].length, actualSingular[2][row].length);
         for (int col = 0; col < expectedV[row].length; col++)
            assertEquals(expectedV[row][col], actualSingular[2][row][col],
               DELTA);
      }

      double[][] actual = Matrices.multiply(
         Matrices.multiply(actualSingular[0], actualSingular[1]),
         Matrices.transpose(actualSingular[2]));
      assertEquals(a.length, actual.length);
      for (int row = 0; row < actual.length; row++) {
         assertEquals(a[row].length, actual[row].length);
         for (int col = 0; col < actual[row].length; col++)
            assertEquals(a[row][col], actual[row][col], DELTA);
      }

      actual = Matrices.multiply(Matrices.transpose(actualSingular[0]),
         actualSingular[0]);
      for (int row = 0; row < actual.length; row++)
         for (int col = 0; col < actual[row].length; col++)
            if (row == col) assertEquals(1, actual[row][col], DELTA);

            else assertEquals(0, actual[row][col], DELTA);

      actual = Matrices.multiply(Matrices.transpose(actualSingular[2]),
         actualSingular[2]);
      for (int row = 0; row < actual.length; row++)
         for (int col = 0; col < actual[row].length; col++)
            if (row == col) assertEquals(1, actual[row][col], DELTA);

            else assertEquals(0, actual[row][col], DELTA);

      for (int row = 0; row < actualSingular[1].length; row++)
         for (int col = 0; col < actualSingular[1][row].length; col++)
            if (row == col)
               assertNotEquals(0, actualSingular[1][row][col], DELTA);
            else assertEquals(0, actualSingular[1][row][col], DELTA);
   }

}