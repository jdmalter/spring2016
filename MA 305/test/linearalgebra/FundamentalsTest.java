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

   /**
    * Test method for {@link Fundamentals#forwardEliminate(double[][])}.
    */
   @Test
   public void testForwardEliminate() {
      // variables required for testing
      double[][] a;
      double[][] expected;
      double[][] test = null;

      // first test on 3 by 4 matrix
      a = new double[][] { { 3, -6, 6, -6 }, { -3, 9, -8, 0 },
         { 3, -9, 6, 6 } };
      expected =
         new double[][] { { 3, -6, 6, -6 }, { 0, 3, -2, -6 }, { 0, 0, -2, 6 } };
      testcaseForwardEliminate(a, expected);

      // second test on 3 by 4 matrix
      a = new double[][] { { 4, -6, 8, -2 }, { -4, 3, -2, -1 },
         { 4, 0, -2, 0 } };
      expected =
         new double[][] { { 4, -6, 8, -2 }, { 0, -3, 6, -3 }, { 0, 0, 2, -4 } };
      testcaseForwardEliminate(a, expected);

      try {
         // matrix a is null
         test = Fundamentals.forwardEliminate(null);
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix a does contain null column
         test = Fundamentals.forwardEliminate(new double[][] { null });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix a must does not have as many or more columns than rows
         test = Fundamentals
            .forwardEliminate(new double[][] { { 0, 0 }, { 0, 0 }, { 0, 0 } });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix a row width does not remain constant
         test =
            Fundamentals.forwardEliminate(new double[][] { { 0, 0 }, { 0 } });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix a row width does not remain constant
         test =
            Fundamentals.forwardEliminate(new double[][] { { 0, 0 }, { 0 } });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }
   }

   private void testcaseForwardEliminate(double[][] a, double[][] expected) {
      double[][] actual = Fundamentals.forwardEliminate(a);

      assertEquals(expected.length, actual.length);
      for (int row = 0; row < expected.length; row++) {
         assertEquals(expected[row].length, actual[row].length);
         for (int col = 0; col < expected[row].length; col++)
            assertEquals(expected[row][col], actual[row][col], DELTA);
      }
   }

   /**
    * Test method for {@link Fundamentals#backwardEliminate(double[][])}.
    */
   @Test
   public void testBackwardEliminate() {
      // variables required for testing
      double[][] a;
      double[][] expected;
      double[][] test = null;

      // first test on 3 by 4 matrix
      a = new double[][] { { 3, -6, 6, -6 }, { 0, 3, -2, -6 },
         { 0, 0, -2, 6 } };
      expected =
         new double[][] { { 3, 0, 0, -12 }, { 0, 3, 0, -12 }, { 0, 0, -2, 6 } };
      testcaseBackwardEliminate(a, expected);

      // second test on 3 by 4 matrix
      a = new double[][] { { 4, -6, 8, -2 }, { 0, -3, 6, -3 },
         { 0, 0, 2, -4 } };
      expected =
         new double[][] { { 4, 0, 0, -4 }, { 0, -3, 0, 9 }, { 0, 0, 2, -4 } };
      testcaseBackwardEliminate(a, expected);

      try {
         // matrix a is null
         test = Fundamentals.backwardEliminate(null);
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix a does contain null column
         test = Fundamentals.backwardEliminate(new double[][] { null });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix a must does not have as many or more columns than rows
         test = Fundamentals
            .backwardEliminate(new double[][] { { 0, 0 }, { 0, 0 }, { 0, 0 } });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix a row width does not remain constant
         test =
            Fundamentals.backwardEliminate(new double[][] { { 0, 0 }, { 0 } });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix a row width does not remain constant
         test =
            Fundamentals.backwardEliminate(new double[][] { { 0, 0 }, { 0 } });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }
   }

   private void testcaseBackwardEliminate(double[][] a, double[][] expected) {
      double[][] actual = Fundamentals.backwardEliminate(a);

      assertEquals(expected.length, actual.length);
      for (int row = 0; row < expected.length; row++) {
         assertEquals(expected[row].length, actual[row].length);
         for (int col = 0; col < expected[row].length; col++)
            assertEquals(expected[row][col], actual[row][col], DELTA);
      }
   }

   /**
    * Test method for {@link Fundamentals#rref(double[][])}.
    */
   @Test
   public void testRREF() {
      double[][] a;
      double[][] expected;
      double[][] test = null;

      // first test on 1 by 1 array
      a = new double[][] { { 7 } };
      expected = new double[][] { { 1 } };
      testcaseRREF(a, expected);

      // second test on 2 by 2 array
      a = new double[][] { { 10, -60 }, { 15, 60 } };
      expected = new double[][] { { 1, -6 }, { 0, 0 } };
      testcaseRREF(a, expected);

      // third test on 3 by 3 array
      a = new double[][] { { -1, 4, -4 }, { -1, 0, -2 }, { -3, 0, -4 } };
      expected = new double[][] { { 1, 0, 2 }, { 0, 1, -0.5 }, { 0, 0, 0 } };
      testcaseRREF(a, expected);

      // fourth test on 3 by 4 array
      a = new double[][] { { 3, -6, 6, -6 }, { -3, 9, -8, 0 },
         { 3, -9, 6, 6 } };
      expected =
         new double[][] { { 1, 0, 0, -4 }, { 0, 1, 0, -4 }, { 0, 0, 1, -3 } };
      testcaseRREF(a, expected);

      try {
         // matrix a is null
         test = Fundamentals.rref(null);
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix a does contain null column
         test = Fundamentals.rref(new double[][] { null });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix a must does not have as many or more columns than rows
         test =
            Fundamentals.rref(new double[][] { { 0, 0 }, { 0, 0 }, { 0, 0 } });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix a row width does not remain constant
         test = Fundamentals.rref(new double[][] { { 0, 0 }, { 0 } });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix a row width does not remain constant
         test = Fundamentals.rref(new double[][] { { 0, 0 }, { 0 } });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }
   }

   private void testcaseRREF(double[][] a, double[][] expected) {
      double[][] actual = Fundamentals.rref(a);

      assertEquals(expected.length, actual.length);
      for (int row = 0; row < expected.length; row++) {
         assertEquals(expected[row].length, actual[row].length);
         for (int col = 0; col < expected[row].length; col++)
            assertEquals(expected[row][col], actual[row][col], DELTA);
      }
   }

   /**
    * Test method for {@link Fundamentals#gaussJordan(double[][], double[][])}.
    */
   @Test
   public void testGaussJordan() {
      // variables required for testing
      double[][] a;
      double[][] b;
      double[][] expected;
      double[][] test = null;

      // first test on 3 by 3 and 3 by 1 arrays
      a = new double[][] { { -3, 1, 2 }, { 9, -4, -8 }, { -9, 5, 8 } };
      b = new double[][] { { -3 }, { 3 }, { 1 } };
      expected = new double[][] { { 3 }, { 4 }, { 1 } };
      testcaseGaussJordan(a, b, expected);

      // second test on 2 by 2 and 2 by 1 arrays
      a = new double[][] { { 4, -6, 8 }, { -4, 3, -2 }, { 4, 0, -2 } };
      b = new double[][] { { -2 }, { -1 }, { 0 } };
      expected = new double[][] { { -1 }, { -3 }, { -2 } };
      testcaseGaussJordan(a, b, expected);

      try {
         // matrix a is null
         test = Fundamentals.gaussJordan(null, new double[][] {});
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix b is null
         test = Fundamentals.gaussJordan(new double[][] {}, null);
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix a and b do not have equal number of rows
         test =
            Fundamentals.gaussJordan(new double[][] { {} }, new double[][] {});
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix a does contain null column
         test = Fundamentals.gaussJordan(new double[][] { null },
            new double[][] { null });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix b does contain null column
         test = Fundamentals.gaussJordan(new double[][] { {} },
            new double[][] { null });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix a row width does not remain constant
         test = Fundamentals.gaussJordan(new double[][] { { 0, 0 }, { 0 } },
            new double[][] { {}, {} });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix b row width does not remain constant
         test = Fundamentals.gaussJordan(new double[][] { { 0 }, { 0 } },
            new double[][] { { 0, 0 }, { 0 } });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // no solution exists
         test = Fundamentals.gaussJordan(
            new double[][] { { -2, 6, -4 }, { 2, -4, 2 }, { 2, -4, 2 } },
            new double[][] { { -8 }, { -6 }, { -2 } });
         fail();
      } catch (ArithmeticException ex) {
         assertNull(test);
      }

      try {
         // many solutions exist
         test = Fundamentals.gaussJordan(
            new double[][] { { 4, -2, 2 }, { 4, 0, 4 }, { 16, -16, 0 } },
            new double[][] { { 4 }, { 8 }, { 0 } });
         fail();
      } catch (ArithmeticException ex) {
         assertNull(test);
      }
   }

   private void testcaseGaussJordan(double[][] a, double[][] b,
      double[][] expected) {
      double[][] actual = Fundamentals.gaussJordan(a, b);

      assertEquals(expected.length, actual.length);
      for (int row = 0; row < expected.length; row++) {
         assertEquals(expected[row].length, actual[row].length);
         for (int col = 0; col < expected[row].length; col++)
            assertEquals(expected[row][col], actual[row][col], DELTA);
      }

      actual = Matrices.multiply(a, actual);
      assertEquals(b.length, actual.length);
      for (int row = 0; row < b.length; row++) {
         assertEquals(b[row].length, actual[row].length);
         for (int col = 0; col < b[row].length; col++)
            assertEquals(b[row][col], actual[row][col], DELTA);
      }
   }

   /**
    * Test method for
    * {@link Fundamentals#gaussJordanMany(double[][], double[][])}.
    */
   @Test
   public void testGaussJordanMany() {
      // variables required for testing
      double[][] a;
      double[][] b;
      double[][] expectedX;
      double[][] expectedAB;
      double[][][] test = null;

      // first test on 3 by 3 and 3 by 1 arrays
      a = new double[][] { { -3, 1, 2 }, { 9, -4, -8 }, { -9, 5, 8 } };
      b = new double[][] { { -3 }, { 3 }, { 1 } };
      expectedAB = new double[][] { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
      expectedX = new double[][] { { 3 }, { 4 }, { 1 } };
      testcaseGaussJordanMany(a, b, expectedAB, expectedX);

      // second test on 2 by 2 and 2 by 1 arrays
      a = new double[][] { { 4, -6, 8 }, { -4, 3, -2 }, { 4, 0, -2 } };
      b = new double[][] { { -2 }, { -1 }, { 0 } };
      expectedAB = new double[][] { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
      expectedX = new double[][] { { -1 }, { -3 }, { -2 } };
      testcaseGaussJordanMany(a, b, expectedAB, expectedX);

      // third test on 3 by 3 and 3 by 1 arrays
      a = new double[][] { { 4, -2, 2 }, { 4, 0, 4 }, { 16, -16, 0 } };
      b = new double[][] { { 4 }, { 8 }, { 0 } };
      expectedAB = new double[][] { { 1, 0, -1 }, { 0, 1, -1 }, { 0, 0, 1 } };
      expectedX = new double[][] { { 2 }, { 2 }, { 0 } };
      testcaseGaussJordanMany(a, b, expectedAB, expectedX);

      // fourth test on 3 by 3 and 3 by 1 arrays
      a = new double[][] { { -2, -6, -4 }, { 6, 18, 12 }, { 2, 6, 4 } };
      b = new double[][] { { -4 }, { 12 }, { 4 } };
      expectedAB = new double[][] { { 1, -3, -2 }, { 0, 1, 0 }, { 0, 0, 1 } };
      expectedX = new double[][] { { 2 }, { 0 }, { 0 } };
      testcaseGaussJordanMany(a, b, expectedAB, expectedX);

      // fifth test on 4 by 4 and 4 by 1 arrays
      a = new double[][] { { 2, 1, 2, 10 }, { -2, -2, 2, -12 },
         { -2, -3, 6, -14 }, { 4, 0, 12, 16 } };
      b = new double[][] { { -7 }, { 10 }, { 13 }, { -8 } };
      expectedAB = new double[][] { { 1, 0, -3, -4 }, { 0, 1, 4, -2 },
         { 0, 0, 1, 0 }, { 0, 0, 0, 1 } };
      expectedX = new double[][] { { -2 }, { -3 }, { 0 }, { 0 } };
      testcaseGaussJordanMany(a, b, expectedAB, expectedX);

      // sixth test on 4 by 4 and 4 by 1 arrays
      a = new double[][] { { 2, 6, 2, 2 }, { 4, 12, 4, 4 }, { -6, -18, -6, -6 },
         { 2, 6, 2, 2 } };
      b = new double[][] { { -4 }, { -8 }, { 12 }, { -4 } };
      expectedAB = new double[][] { { 1, -3, -1, -1 }, { 0, 1, 0, 0 },
         { 0, 0, 1, 0 }, { 0, 0, 0, 1 } };
      expectedX = new double[][] { { -2 }, { 0 }, { 0 }, { 0 } };
      testcaseGaussJordanMany(a, b, expectedAB, expectedX);

      try {
         // matrix a is null
         test = Fundamentals.gaussJordanMany(null, new double[][] {});
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix b is null
         test = Fundamentals.gaussJordanMany(new double[][] {}, null);
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix a and b do not have equal number of rows
         test = Fundamentals.gaussJordanMany(new double[][] { {} },
            new double[][] {});
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix a does contain null column
         test = Fundamentals.gaussJordanMany(new double[][] { null },
            new double[][] { null });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix b does contain null column
         test = Fundamentals.gaussJordanMany(new double[][] { {} },
            new double[][] { null });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix a row width does not remain constant
         test = Fundamentals.gaussJordanMany(new double[][] { { 0, 0 }, { 0 } },
            new double[][] { {}, {} });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix b row width does not remain constant
         test = Fundamentals.gaussJordanMany(new double[][] { { 0 }, { 0 } },
            new double[][] { { 0, 0 }, { 0 } });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // no solution exists
         test = Fundamentals.gaussJordanMany(
            new double[][] { { -2, 6, -4 }, { 2, -4, 2 }, { 2, -4, 2 } },
            new double[][] { { -8 }, { -6 }, { -2 } });
         fail();
      } catch (ArithmeticException ex) {
         assertNull(test);
      }

      try {
         // no solution exists
         test = Fundamentals.gaussJordanMany(
            new double[][] { { 3, 2, -2, 8 }, { 12, 7, -12, 19 },
               { -12, -6, 18, 0 }, { -9, -4, 18, 14 } },
            new double[][] { { 3 }, { 13 }, { -6 }, { 2 } });
         fail();
      } catch (ArithmeticException ex) {
         assertNull(test);
      }

      try {
         // no solution exists
         test = Fundamentals.gaussJordanMany(
            new double[][] { { -2, 2, 2, -4 }, { 2, -2, -2, 4 },
               { 8, -8, -8, 16 }, { -4, 4, 4, -8 } },
            new double[][] { { 6 }, { -10 }, { -8 }, { -4 } });
         fail();
      } catch (ArithmeticException ex) {
         assertNull(test);
      }
   }

   private void testcaseGaussJordanMany(double[][] a, double[][] b,
      double[][] expectedAB, double[][] expectedX) {
      double[][][] actualABX = Fundamentals.gaussJordanMany(a, b);

      assertEquals(expectedAB.length, actualABX[0].length);
      for (int row = 0; row < expectedAB.length; row++) {
         assertEquals(expectedAB[row].length, actualABX[0][row].length);
         for (int col = 0; col < expectedAB[row].length; col++)
            assertEquals(expectedAB[row][col], actualABX[0][row][col], DELTA);
      }

      assertEquals(expectedX.length, actualABX[1].length);
      for (int row = 0; row < expectedX.length; row++) {
         assertEquals(expectedX[row].length, actualABX[1][row].length);
         for (int col = 0; col < expectedX[row].length; col++)
            assertEquals(expectedX[row][col], actualABX[1][row][col], DELTA);
      }
   }

   /**
    * Test method for {@link Fundamentals#cramer(double[][], double[][])}.
    */
   @Test
   public void testCramer() {
      // variables required for testing
      double[][] a;
      double[][] b;
      double[][] expected;
      double[][] test = null;

      // first test on 3 by 3 and 3 by 1 arrays
      a = new double[][] { { -3, 1, 2 }, { 9, -4, -8 }, { -9, 5, 8 } };
      b = new double[][] { { -3 }, { 3 }, { 1 } };
      expected = new double[][] { { 3 }, { 4 }, { 1 } };
      testcaseCramer(a, b, expected);

      // second test on 2 by 2 and 2 by 1 arrays
      a = new double[][] { { 4, -6, 8 }, { -4, 3, -2 }, { 4, 0, -2 } };
      b = new double[][] { { -2 }, { -1 }, { 0 } };
      expected = new double[][] { { -1 }, { -3 }, { -2 } };
      testcaseCramer(a, b, expected);

      try {
         // matrix a is null
         test = Fundamentals.cramer(null, new double[][] {});
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix b is null
         test = Fundamentals.cramer(new double[][] {}, null);
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix a and b do not have equal number of rows
         test = Fundamentals.cramer(new double[][] { {} }, new double[][] {});
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix a does contain null column
         test = Fundamentals.cramer(new double[][] { null },
            new double[][] { null });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix b does contain null column
         test =
            Fundamentals.cramer(new double[][] { {} }, new double[][] { null });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix a row width does not remain constant
         test = Fundamentals.cramer(new double[][] { { 0, 0 }, { 0 } },
            new double[][] { {}, {} });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix b row width does not remain constant
         test = Fundamentals.cramer(new double[][] { { 0 }, { 0 } },
            new double[][] { { 0, 0 }, { 0 } });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }
   }

   private void testcaseCramer(double[][] a, double[][] b,
      double[][] expected) {
      double[][] actual = Fundamentals.cramer(a, b);

      assertEquals(expected.length, actual.length);
      for (int row = 0; row < expected.length; row++) {
         assertEquals(expected[row].length, actual[row].length);
         for (int col = 0; col < expected[row].length; col++)
            assertEquals(expected[row][col], actual[row][col], DELTA);
      }

      actual = Matrices.multiply(a, actual);
      assertEquals(b.length, actual.length);
      for (int row = 0; row < b.length; row++) {
         assertEquals(b[row].length, actual[row].length);
         for (int col = 0; col < b[row].length; col++)
            assertEquals(b[row][col], actual[row][col], DELTA);
      }
   }

   /**
    * Test method for {@link Fundamentals#invert(double[][])}.
    */
   @Test
   public void testInvert() {
      // variables required for testing
      double[][] a;
      double[][] expected;
      double[][] test = null;

      // first test on 2 by 2 array
      a = new double[][] { { -1, 4 }, { 1, -5 } };
      expected = new double[][] { { -5, -4 }, { -1, -1 } };
      testcaseInvert(a, expected);

      // second test on 3 by 3 array
      a = new double[][] { { -1, 4, -4 }, { -1, 0, -2 }, { -3, 0, -4 } };
      expected =
         new double[][] { { 0, 2, -1 }, { 0.25, -1, 0.25 }, { 0, -1.5, 0.5 } };
      testcaseInvert(a, expected);

      try {
         // matrix a is null
         test = Fundamentals.invert(null);
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix a does contain null column
         test = Fundamentals.invert(new double[][] { null });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix a row width does not remain constant
         test = Fundamentals.invert(new double[][] { { 0, 0 }, { 0 } });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // divide by 0
         test = Fundamentals.invert(new double[][] { { 0, 0 }, { 0, 0 } });
         fail();
      } catch (ArithmeticException ex) {
         assertNull(test);
      }
   }

   private void testcaseInvert(double[][] a, double[][] expected) {
      double[][] actual = Fundamentals.invert(a);

      assertEquals(expected.length, actual.length);
      for (int row = 0; row < expected.length; row++) {
         assertEquals(expected[row].length, actual[row].length);
         for (int col = 0; col < expected[row].length; col++)
            assertEquals(expected[row][col], actual[row][col], DELTA);
      }

      actual = Matrices.multiply(a, actual);
      for (int row = 0; row < actual.length; row++)
         for (int col = 0; col < actual[row].length; col++)
            if (row == col) assertEquals(1, actual[row][col], DELTA);
            else assertEquals(0, actual[row][col], DELTA);
   }

   /**
    * Test method for {@link Fundamentals#invertFaddeev(double[][])}.
    */
   @Test
   public void testInvertFaddeev() {
      // variables required for testing
      double[][] a;
      double[][] expected;
      double[][] test = null;

      // first test on 2 by 2 array
      a = new double[][] { { -1, 4 }, { 1, -5 } };
      expected = new double[][] { { -5, -4 }, { -1, -1 } };
      testcaseInvertFaddeev(a, expected);

      // second test on 3 by 3 array
      a = new double[][] { { -1, 4, -4 }, { -1, 0, -2 }, { -3, 0, -4 } };
      expected =
         new double[][] { { 0, 2, -1 }, { 0.25, -1, 0.25 }, { 0, -1.5, 0.5 } };
      testcaseInvertFaddeev(a, expected);

      try {
         // matrix a is null
         test = Fundamentals.invertFaddeev(null);
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix a does contain null column
         test = Fundamentals.invertFaddeev(new double[][] { null });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix a row width does not remain constant
         test = Fundamentals.invertFaddeev(new double[][] { { 0, 0 }, { 0 } });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // divide by 0
         test =
            Fundamentals.invertFaddeev(new double[][] { { 0, 0 }, { 0, 0 } });
         fail();
      } catch (ArithmeticException ex) {
         assertNull(test);
      }
   }

   private void testcaseInvertFaddeev(double[][] a, double[][] expected) {
      double[][] actual = Fundamentals.invertFaddeev(a);

      assertEquals(expected.length, actual.length);
      for (int row = 0; row < expected.length; row++) {
         assertEquals(expected[row].length, actual[row].length);
         for (int col = 0; col < expected[row].length; col++)
            assertEquals(expected[row][col], actual[row][col], DELTA);
      }

      actual = Matrices.multiply(a, actual);
      for (int row = 0; row < actual.length; row++)
         for (int col = 0; col < actual[row].length; col++)
            if (row == col) assertEquals(1, actual[row][col], DELTA);
            else assertEquals(0, actual[row][col], DELTA);
   }

}