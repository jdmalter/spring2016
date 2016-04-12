/**
 * 
 */
package linearalgebra;

import static org.junit.Assert.*;

import org.junit.Test;

/**
 * Constrained test.
 * 
 * @author Jacob Malter
 *
 */
public class ContrainedTest {

   /** acceptable margin of error */
   private static final double DELTA = 0.001d;

   /**
    * Test method for
    * {@link Constrained#overconstrainedPI(double[][], double[][])}.
    */
   @Test
   public void testOverconstrainedPI() {
      double[][] a;
      double[][] b;
      double[][] expectedX;
      double[][] expectedE;
      double[][][] test = null;

      // first test on 3 by 2 and 3 by 1 arrays
      a = new double[][] { { 4, -3 }, { -1, 1 }, { -4, -3 } };
      b = new double[][] { { -1 }, { -3 }, { 4 } };
      expectedX = new double[][] { { -0.535 }, { -0.660 } };
      expectedE = new double[][] { { 0.839 }, { 2.875 }, { 0.120 } };
      testcaseOverconstrainedPI(a, b, expectedX, expectedE);

      // second test on 4 by 3 and 4 by 1 arrays
      a = new double[][] { { -3, 4, 3 }, { -3, -2, -4 }, { 4, 2, -1 },
         { 3, 1, -3 } };
      b = new double[][] { { 2 }, { -3 }, { 1 }, { 1 } };
      expectedX = new double[][] { { 0.239 }, { 0.473 }, { 0.266 } };
      expectedE =
         new double[][] { { -0.030 }, { 0.273 }, { 0.637 }, { -0.607 } };
      testcaseOverconstrainedPI(a, b, expectedX, expectedE);

      // third test on 5 by 3 and 5 by 1 arrays
      a = new double[][] { { 3, -1, -1 }, { -3, 2, -3 }, { 2, -4, -2 },
         { -4, -3, -1 }, { 3, 1, -1 } };
      b = new double[][] { { 3 }, { 3 }, { 1 }, { 1 }, { 4 } };
      expectedX = new double[][] { { 0.307 }, { 0.232 }, { -1.318 } };
      expectedE = new double[][] { { -0.994 }, { 0.497 }, { 1.319 }, { -1.606 },
         { -1.530 } };
      testcaseOverconstrainedPI(a, b, expectedX, expectedE);

      try {
         // matrix a is null
         test = Constrained.overconstrainedPI(null, b);
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix b is null
         test = Constrained.overconstrainedPI(a, null);
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix a and b do not have equal number of rows
         test = Constrained.overconstrainedPI(new double[][] { {} },
            new double[][] {});
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix a does contain null column
         test = Constrained.overconstrainedPI(new double[][] { null },
            new double[][] { null });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix b does contain null column
         test = Constrained.overconstrainedPI(new double[][] { {} },
            new double[][] { null });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix a must does not have more rows than columns
         test =
            Constrained.overconstrainedPI(new double[][] { { 0, 0 }, { 0, 0 } },
               new double[][] { { 0 }, { 0 } });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix a row width does not remain constant
         test = Constrained.overconstrainedPI(
            new double[][] { { 0, 0 }, { 0, 0 }, { 0, 0, 0 } },
            new double[][] { { 0 }, { 0 }, { 0 } });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix b row width does not remain constant
         test = Constrained.overconstrainedPI(
            new double[][] { { 0, 0 }, { 0, 0 }, { 0, 0 } },
            new double[][] { { 0 }, { 0 }, { 0, 0 } });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }
   }

   private void testcaseOverconstrainedPI(double[][] a, double[][] b,
      double[][] expectedX, double[][] expectedE) {
      double[][][] actual = Constrained.overconstrainedPI(a, b);

      assertEquals(expectedX.length, actual[0].length);
      for (int row = 0; row < expectedX.length; row++) {
         assertEquals(expectedX[row].length, actual[0][row].length);
         for (int col = 0; col < expectedX[row].length; col++)
            assertEquals(expectedX[row][col], actual[0][row][col], DELTA);
      }

      assertEquals(expectedE.length, actual[1].length);
      for (int row = 0; row < expectedE.length; row++) {
         assertEquals(expectedE[row].length, actual[1][row].length);
         for (int col = 0; col < expectedE[row].length; col++)
            assertEquals(expectedE[row][col], actual[1][row][col], DELTA);
      }

      double[][] actualAX = Matrices.multiply(a, actual[0]);
      assertEquals(expectedE.length, actualAX.length);
      assertEquals(expectedE.length, b.length);
      for (int row = 0; row < expectedE.length; row++) {
         assertEquals(expectedE[row].length, actualAX[row].length);
         assertEquals(expectedE[row].length, b[row].length);
         for (int col = 0; col < expectedE[row].length; col++)
            assertEquals(expectedE[row][col], actualAX[row][col] - b[row][col],
               DELTA);
      }
   }

   /**
    * Test method for
    * {@link Constrained#underconstrainedPI(double[][], double[][])}.
    */
   @Test
   public void testUnderconstrainedPI() {
      double[][] a;
      double[][] b;
      double[][] expected;
      double[][] test = null;

      // first test on 2 by 3 and 2 by 1 arrays
      a = new double[][] { { 4, -3, -1 }, { 1, -4, -3 } };
      b = new double[][] { { -1 }, { -3 } };
      expected = new double[][] { { 0.206 }, { 0.454 }, { 0.463 } };
      testcaseUnderconstrainedPI(a, b, expected);

      // second test on 3 by 4 and 3 by 1 arrays
      a = new double[][] { { 4, -3, -1, 1 }, { -4, -3, -1, -3 },
         { 4, 2, 3, 3 } };
      b = new double[][] { { 3 }, { 1 }, { -1 } };
      expected = new double[][] { { 0.279 }, { -0.555 },
         new double[] { -0.277 }, { -0.058 } };
      testcaseUnderconstrainedPI(a, b, expected);

      try {
         // matrix a is null
         test = Constrained.underconstrainedPI(null, b);
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix b is null
         test = Constrained.underconstrainedPI(a, null);
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix a and b do not have equal number of rows
         test = Constrained.underconstrainedPI(new double[][] { {} },
            new double[][] {});
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix a does contain null column
         test = Constrained.underconstrainedPI(new double[][] { null },
            new double[][] { null });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix b does contain null column
         test = Constrained.underconstrainedPI(new double[][] { {} },
            new double[][] { null });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix a must does not have more rows than columns
         test = Constrained.underconstrainedPI(
            new double[][] { { 0, 0 }, { 0, 0 } },
            new double[][] { { 0 }, { 0 } });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix a row width does not remain constant
         test = Constrained.underconstrainedPI(
            new double[][] { { 0, 0, 0 }, { 0, 0, 0, 0 } },
            new double[][] { { 0 }, { 0 } });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix b row width does not remain constant
         test = Constrained.underconstrainedPI(
            new double[][] { { 0, 0, 0 }, { 0, 0, 0 } },
            new double[][] { { 0 }, { 0, 0 } });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }
   }

   private void testcaseUnderconstrainedPI(double[][] a, double[][] b,
      double[][] expectedX) {
      double[][] actual = Constrained.underconstrainedPI(a, b);

      assertEquals(expectedX.length, actual.length);
      for (int row = 0; row < expectedX.length; row++) {
         assertEquals(expectedX[row].length, actual[row].length);
         for (int col = 0; col < expectedX[row].length; col++)
            assertEquals(expectedX[row][col], actual[row][col], DELTA);
      }

      double[][] actualAX = Matrices.multiply(a, actual);
      assertEquals(b.length, actualAX.length);
      assertEquals(b.length, b.length);
      for (int row = 0; row < b.length; row++) {
         assertEquals(b[row].length, actualAX[row].length);
         assertEquals(b[row].length, b[row].length);
         for (int col = 0; col < b[row].length; col++)
            assertEquals(b[row][col], actualAX[row][col], DELTA);
      }
   }

   /**
    * Test method for
    * {@link Constrained#overconstrainedCD(double[][], double[][])}.
    */
   @Test
   public void testOverconstrainedCD() {
      double[][] a;
      double[][] b;
      double[][] expectedX;
      double[][] expectedE;
      double[][][] test = null;

      // first test on 3 by 2 and 3 by 1 arrays
      a = new double[][] { { 4, -3 }, { -1, 1 }, { -4, -3 } };
      b = new double[][] { { -1 }, { -3 }, { 4 } };
      expectedX = new double[][] { { -0.535 }, { -0.660 } };
      expectedE = new double[][] { { 0.839 }, { 2.875 }, { 0.120 } };
      testcaseOverconstrainedCD(a, b, expectedX, expectedE);

      // second test on 4 by 3 and 4 by 1 arrays
      a = new double[][] { { -3, 4, 3 }, { -3, -2, -4 }, { 4, 2, -1 },
         { 3, 1, -3 } };
      b = new double[][] { { 2 }, { -3 }, { 1 }, { 1 } };
      expectedX = new double[][] { { 0.239 }, { 0.473 }, { 0.266 } };
      expectedE =
         new double[][] { { -0.030 }, { 0.273 }, { 0.637 }, { -0.607 } };
      testcaseOverconstrainedCD(a, b, expectedX, expectedE);

      // third test on 5 by 3 and 5 by 1 arrays
      a = new double[][] { { 3, -1, -1 }, { -3, 2, -3 }, { 2, -4, -2 },
         { -4, -3, -1 }, { 3, 1, -1 } };
      b = new double[][] { { 3 }, { 3 }, { 1 }, { 1 }, { 4 } };
      expectedX = new double[][] { { 0.307 }, { 0.232 }, { -1.318 } };
      expectedE = new double[][] { { -0.994 }, { 0.497 }, { 1.319 }, { -1.606 },
         { -1.530 } };
      testcaseOverconstrainedCD(a, b, expectedX, expectedE);

      try {
         // matrix a is null
         test = Constrained.overconstrainedCD(null, b);
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix b is null
         test = Constrained.overconstrainedCD(a, null);
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix a and b do not have equal number of rows
         test = Constrained.overconstrainedCD(new double[][] { {} },
            new double[][] {});
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix a does contain null column
         test = Constrained.overconstrainedCD(new double[][] { null },
            new double[][] { null });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix b does contain null column
         test = Constrained.overconstrainedCD(new double[][] { {} },
            new double[][] { null });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix a must does not have more rows than columns
         test =
            Constrained.overconstrainedCD(new double[][] { { 0, 0 }, { 0, 0 } },
               new double[][] { { 0 }, { 0 } });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix a row width does not remain constant
         test = Constrained.overconstrainedCD(
            new double[][] { { 0, 0 }, { 0, 0 }, { 0, 0, 0 } },
            new double[][] { { 0 }, { 0 }, { 0 } });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix b row width does not remain constant
         test = Constrained.overconstrainedCD(
            new double[][] { { 0, 0 }, { 0, 0 }, { 0, 0 } },
            new double[][] { { 0 }, { 0 }, { 0, 0 } });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }
   }

   private void testcaseOverconstrainedCD(double[][] a, double[][] b,
      double[][] expectedX, double[][] expectedE) {
      double[][][] actual = Constrained.overconstrainedCD(a, b);

      assertEquals(expectedX.length, actual[0].length);
      for (int row = 0; row < expectedX.length; row++) {
         assertEquals(expectedX[row].length, actual[0][row].length);
         for (int col = 0; col < expectedX[row].length; col++)
            assertEquals(expectedX[row][col], actual[0][row][col], DELTA);
      }

      assertEquals(expectedE.length, actual[1].length);
      for (int row = 0; row < expectedE.length; row++) {
         assertEquals(expectedE[row].length, actual[1][row].length);
         for (int col = 0; col < expectedE[row].length; col++)
            assertEquals(expectedE[row][col], actual[1][row][col], DELTA);
      }

      double[][] actualAX = Matrices.multiply(a, actual[0]);
      assertEquals(expectedE.length, actualAX.length);
      assertEquals(expectedE.length, b.length);
      for (int row = 0; row < expectedE.length; row++) {
         assertEquals(expectedE[row].length, actualAX[row].length);
         assertEquals(expectedE[row].length, b[row].length);
         for (int col = 0; col < expectedE[row].length; col++)
            assertEquals(expectedE[row][col], actualAX[row][col] - b[row][col],
               DELTA);
      }
   }

   /**
    * Test method for
    * {@link Constrained#overconstrainedQR(double[][], double[][])}.
    */
   @Test
   public void testOverconstrainedQR() {
      double[][] a;
      double[][] b;
      double[][] expectedX;
      double[][] expectedE;
      double[][][] test = null;

      // first test on 3 by 2 and 3 by 1 arrays
      a = new double[][] { { 4, -3 }, { -1, 1 }, { -4, -3 } };
      b = new double[][] { { -1 }, { -3 }, { 4 } };
      expectedX = new double[][] { { -0.535 }, { -0.660 } };
      expectedE = new double[][] { { 0.839 }, { 2.875 }, { 0.120 } };
      testcaseOverconstrainedQR(a, b, expectedX, expectedE);

      // second test on 4 by 3 and 4 by 1 arrays
      a = new double[][] { { -3, 4, 3 }, { -3, -2, -4 }, { 4, 2, -1 },
         { 3, 1, -3 } };
      b = new double[][] { { 2 }, { -3 }, { 1 }, { 1 } };
      expectedX = new double[][] { { 0.239 }, { 0.473 }, { 0.266 } };
      expectedE =
         new double[][] { { -0.030 }, { 0.273 }, { 0.637 }, { -0.607 } };
      testcaseOverconstrainedQR(a, b, expectedX, expectedE);

      // third test on 5 by 3 and 5 by 1 arrays
      a = new double[][] { { 3, -1, -1 }, { -3, 2, -3 }, { 2, -4, -2 },
         { -4, -3, -1 }, { 3, 1, -1 } };
      b = new double[][] { { 3 }, { 3 }, { 1 }, { 1 }, { 4 } };
      expectedX = new double[][] { { 0.307 }, { 0.232 }, { -1.318 } };
      expectedE = new double[][] { { -0.994 }, { 0.497 }, { 1.319 }, { -1.606 },
         { -1.530 } };
      testcaseOverconstrainedQR(a, b, expectedX, expectedE);

      try {
         // matrix a is null
         test = Constrained.overconstrainedQR(null, b);
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix b is null
         test = Constrained.overconstrainedQR(a, null);
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix a and b do not have equal number of rows
         test = Constrained.overconstrainedQR(new double[][] { {} },
            new double[][] {});
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix a does contain null column
         test = Constrained.overconstrainedQR(new double[][] { null },
            new double[][] { null });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix b does contain null column
         test = Constrained.overconstrainedQR(new double[][] { {} },
            new double[][] { null });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix a must does not have more rows than columns
         test =
            Constrained.overconstrainedQR(new double[][] { { 0, 0 }, { 0, 0 } },
               new double[][] { { 0 }, { 0 } });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix a row width does not remain constant
         test = Constrained.overconstrainedQR(
            new double[][] { { 0, 0 }, { 0, 0 }, { 0, 0, 0 } },
            new double[][] { { 0 }, { 0 }, { 0 } });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix b row width does not remain constant
         test = Constrained.overconstrainedQR(
            new double[][] { { 0, 0 }, { 0, 0 }, { 0, 0 } },
            new double[][] { { 0 }, { 0 }, { 0, 0 } });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }
   }

   private void testcaseOverconstrainedQR(double[][] a, double[][] b,
      double[][] expectedX, double[][] expectedE) {
      double[][][] actual = Constrained.overconstrainedQR(a, b);

      assertEquals(expectedX.length, actual[0].length);
      for (int row = 0; row < expectedX.length; row++) {
         assertEquals(expectedX[row].length, actual[0][row].length);
         for (int col = 0; col < expectedX[row].length; col++)
            assertEquals(expectedX[row][col], actual[0][row][col], DELTA);
      }

      assertEquals(expectedE.length, actual[1].length);
      for (int row = 0; row < expectedE.length; row++) {
         assertEquals(expectedE[row].length, actual[1][row].length);
         for (int col = 0; col < expectedE[row].length; col++)
            assertEquals(expectedE[row][col], actual[1][row][col], DELTA);
      }

      double[][] actualAX = Matrices.multiply(a, actual[0]);
      assertEquals(expectedE.length, actualAX.length);
      assertEquals(expectedE.length, b.length);
      for (int row = 0; row < expectedE.length; row++) {
         assertEquals(expectedE[row].length, actualAX[row].length);
         assertEquals(expectedE[row].length, b[row].length);
         for (int col = 0; col < expectedE[row].length; col++)
            assertEquals(expectedE[row][col], actualAX[row][col] - b[row][col],
               DELTA);
      }
   }

   /**
    * Test method for
    * {@link Constrained#overconstrainedSV(double[][], double[][])}.
    */
   @Test
   public void testOverconstrainedSV() {
      double[][] a;
      double[][] b;
      double[][] expectedX;
      double[][] expectedE;
      double[][][] test = null;

      // first test on 3 by 2 and 3 by 1 arrays
      a = new double[][] { { 4, -3 }, { -1, 1 }, { -4, -3 } };
      b = new double[][] { { -1 }, { -3 }, { 4 } };
      expectedX = new double[][] { { -0.535 }, { -0.660 } };
      expectedE = new double[][] { { 0.839 }, { 2.875 }, { 0.120 } };
      testcaseOverconstrainedSV(a, b, expectedX, expectedE);

      // second test on 4 by 3 and 4 by 1 arrays
      a = new double[][] { { -3, 4, 3 }, { -3, -2, -4 }, { 4, 2, -1 },
         new double[] { 3, 1, -3 } };
      b = new double[][] { { 2 }, { -3 }, { 1 }, { 1 } };
      expectedX = new double[][] { { 0.239 }, { 0.473 }, { 0.266 } };
      expectedE =
         new double[][] { { -0.030 }, { 0.273 }, { 0.637 }, { -0.607 } };
      testcaseOverconstrainedSV(a, b, expectedX, expectedE);

      // third test on 5 by 3 and 5 by 1 arrays
      a = new double[][] { { 3, -1, -1 }, { -3, 2, -3 }, { 2, -4, -2 },
         { -4, -3, -1 }, { 3, 1, -1 } };
      b = new double[][] { { 3 }, { 3 }, { 1 }, { 1 }, { 4 } };
      expectedX = new double[][] { { 0.307 }, { 0.232 }, { -1.318 } };
      expectedE = new double[][] { { -0.994 }, { 0.497 }, { 1.319 }, { -1.606 },
         { -1.530 } };
      testcaseOverconstrainedSV(a, b, expectedX, expectedE);

      try {
         // matrix a is null
         test = Constrained.overconstrainedSV(null, b);
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix b is null
         test = Constrained.overconstrainedSV(a, null);
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix a and b do not have equal number of rows
         test = Constrained.overconstrainedSV(new double[][] { {} },
            new double[][] {});
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix a does contain null column
         test = Constrained.overconstrainedSV(new double[][] { null },
            new double[][] { null });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix b does contain null column
         test = Constrained.overconstrainedSV(new double[][] { {} },
            new double[][] { null });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix a must does not have more rows than columns
         test =
            Constrained.overconstrainedSV(new double[][] { { 0, 0 }, { 0, 0 } },
               new double[][] { { 0 }, { 0 } });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix a row width does not remain constant
         test = Constrained.overconstrainedSV(
            new double[][] { { 0, 0 }, { 0, 0 }, { 0, 0, 0 } },
            new double[][] { { 0 }, { 0 }, { 0 } });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix b row width does not remain constant
         test = Constrained.overconstrainedSV(
            new double[][] { { 0, 0 }, { 0, 0 }, { 0, 0 } },
            new double[][] { { 0 }, { 0 }, { 0, 0 } });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }
   }

   private void testcaseOverconstrainedSV(double[][] a, double[][] b,
      double[][] expectedX, double[][] expectedE) {
      double[][][] actual = Constrained.overconstrainedSV(a, b);

      assertEquals(expectedX.length, actual[0].length);
      for (int row = 0; row < expectedX.length; row++) {
         assertEquals(expectedX[row].length, actual[0][row].length);
         for (int col = 0; col < expectedX[row].length; col++)
            assertEquals(expectedX[row][col], actual[0][row][col], DELTA);
      }

      assertEquals(expectedE.length, actual[1].length);
      for (int row = 0; row < expectedE.length; row++) {
         assertEquals(expectedE[row].length, actual[1][row].length);
         for (int col = 0; col < expectedE[row].length; col++)
            assertEquals(expectedE[row][col], actual[1][row][col], DELTA);
      }

      double[][] actualAX = Matrices.multiply(a, actual[0]);
      assertEquals(expectedE.length, actualAX.length);
      assertEquals(expectedE.length, b.length);
      for (int row = 0; row < expectedE.length; row++) {
         assertEquals(expectedE[row].length, actualAX[row].length);
         assertEquals(expectedE[row].length, b[row].length);
         for (int col = 0; col < expectedE[row].length; col++)
            assertEquals(expectedE[row][col], actualAX[row][col] - b[row][col],
               DELTA);
      }
   }

   /**
    * Test method for
    * {@link Constrained#underconstrainedSV(double[][], double[][])}.
    */
   @Test
   public void testUnderconstrainedSV() {
      double[][] a;
      double[][] b;
      double[][] expected;
      double[][] test = null;

      // first test on 2 by 3 and 2 by 1 arrays
      a = new double[][] { { 4, -3, -1 }, { 1, -4, -3 } };
      b = new double[][] { { -1 }, { -3 } };
      expected = new double[][] { { 0.206 }, { 0.454 }, { 0.463 } };
      testcaseUnderconstrainedSV(a, b, expected);

      // second test on 3 by 4 and 3 by 1 arrays
      a = new double[][] { { 4, -3, -1, 1 }, { -4, -3, -1, -3 },
         { 4, 2, 3, 3 } };
      b = new double[][] { { 3 }, { 1 }, { -1 } };
      expected =
         new double[][] { { 0.279 }, { -0.555 }, { -0.277 }, { -0.058 } };
      testcaseUnderconstrainedSV(a, b, expected);

      try {
         // matrix a is null
         test = Constrained.underconstrainedSV(null, b);
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix b is null
         test = Constrained.underconstrainedSV(a, null);
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix a and b do not have equal number of rows
         test = Constrained.underconstrainedSV(new double[][] { {} },
            new double[][] {});
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix a does contain null column
         test = Constrained.underconstrainedSV(new double[][] { null },
            new double[][] { null });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix b does contain null column
         test = Constrained.underconstrainedSV(new double[][] { {} },
            new double[][] { null });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix a must does not have more rows than columns
         test = Constrained.underconstrainedSV(
            new double[][] { { 0, 0 }, { 0, 0 } },
            new double[][] { { 0 }, { 0 } });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix a row width does not remain constant
         test = Constrained.underconstrainedSV(
            new double[][] { { 0, 0, 0 }, { 0, 0, 0, 0 } },
            new double[][] { { 0 }, { 0 } });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix b row width does not remain constant
         test = Constrained.underconstrainedSV(
            new double[][] { { 0, 0, 0 }, { 0, 0, 0 } },
            new double[][] { { 0 }, { 0, 0 } });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }
   }

   private void testcaseUnderconstrainedSV(double[][] a, double[][] b,
      double[][] expectedX) {
      double[][] actual = Constrained.underconstrainedSV(a, b);

      assertEquals(expectedX.length, actual.length);
      for (int row = 0; row < expectedX.length; row++) {
         assertEquals(expectedX[row].length, actual[row].length);
         for (int col = 0; col < expectedX[row].length; col++)
            assertEquals(expectedX[row][col], actual[row][col], DELTA);
      }

      double[][] actualAX = Matrices.multiply(a, actual);
      assertEquals(b.length, actualAX.length);
      assertEquals(b.length, b.length);
      for (int row = 0; row < b.length; row++) {
         assertEquals(b[row].length, actualAX[row].length);
         assertEquals(b[row].length, b[row].length);
         for (int col = 0; col < b[row].length; col++)
            assertEquals(b[row][col], actualAX[row][col], DELTA);
      }
   }

}