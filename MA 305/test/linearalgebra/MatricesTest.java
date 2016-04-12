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
      double test = 0;

      // property test on identity matrx
      a = new double[][] { { 1, 0, 0, 0, 0, 0 }, { 0, 1, 0, 0, 0, 0 },
         { 0, 0, 1, 0, 0, 0 }, { 0, 0, 0, 1, 0, 0 }, { 0, 0, 0, 0, 1, 0 },
         { 0, 0, 0, 0, 0, 1 } };
      testcaseDeterminant(a, 1);

      // property test on tranpose
      a = new double[][] { { 1, 2, 3 }, { 4, 5, 6 }, { 7, 8, 9 } };
      testcaseDeterminant(a, Matrices.determinant(Matrices.transpose(a)));

      // first test on 0 component matrix
      a = new double[][] {};
      testcaseDeterminant(a, 0);

      // second test on 1 component matrix
      a = new double[][] { { 3 } };
      testcaseDeterminant(a, 3);

      // third test on 2 component matrix
      a = new double[][] { { 3, 8 }, { 4, 6 } };
      testcaseDeterminant(a, -14);

      // fourth test on 3 component matrix
      a = new double[][] { { 6, 1, 1 }, { 4, -2, 5 }, { 2, 8, 7 } };
      testcaseDeterminant(a, -306);

      // fifth test on 4 component matrix
      a = new double[][] { { 0, 2, 0, 0 }, { 6, 0, 1, 1 }, { 4, 0, -2, 5 },
         { 2, 0, 8, 7 } };
      testcaseDeterminant(a, 612);

      // sixth test on 4 component matrix
      a = new double[][] { { 2, 5, 3, 5 }, { 4, 6, 6, 3 }, { 11, 3, 2, -2 },
         { 4, -7, 9, 3 } };
      testcaseDeterminant(a, 2960);

      // seventh test on 5 component matrix
      a = new double[][] { { 276, 1, 179, 23, 9387 }, { 0, 0, 78, 0, 0 },
         { 0, 0, -1, 0, 1 }, { 0, 0, 1994, -1, 1089 },
         { 1, 0, 212, 726, -378 } };
      testcaseDeterminant(a, 78);

      try {
         // array a is null
         test = Matrices.determinant(null);
         fail();
      } catch (IllegalArgumentException ex) {
         assertEquals(0, test, DELTA);
      }

      try {
         // matrix a does contain null column
         test = Matrices.determinant(new double[][] { null });
         fail();
      } catch (IllegalArgumentException ex) {
         assertEquals(0, test, DELTA);
      }

      try {
         // matrix a must does not have as many columns as rows
         test = Matrices
            .determinant(new double[][] { { 0, 0 }, { 0, 0 }, { 0, 0 } });
         fail();
      } catch (IllegalArgumentException ex) {
         assertEquals(0, test, DELTA);
      }
   }

   private void testcaseDeterminant(double[][] a, double expected) {
      double actual = Matrices.determinant(a);

      assertEquals(expected, actual, DELTA);
   }

   /**
    * Test method for {@link Matrices#minor(double[][], int, int)}
    */
   @Test
   public void testMinor() {
      // variables required for testing
      double[][] a;
      double[][] expected;
      double[][] test = null;

      // first test on 1 by 1 array
      a = new double[][] { { 0 } };
      expected = new double[][] {};
      testcaseMinor(a, expected, 0, 0);

      // second test on 2 by 2 array
      a = new double[][] { { 0, 1 }, { 2, 3 } };
      expected = new double[][] { { 1 } };
      testcaseMinor(a, expected, 1, 0);

      // third test on 3 by 3 array
      a = new double[][] { { 3, 0, 2 }, { 2, 0, -2 }, { 0, 1, 1 } };
      expected = new double[][] { { 3, 2 }, { 0, 1 } };
      testcaseMinor(a, expected, 1, 1);

      try {
         // matrix a is null
         test = Matrices.minor(null, 0, 0);
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // skipRow is less than 0
         test = Matrices.minor(new double[][] {}, -1, 0);
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // skipCol is less than 0
         test = Matrices.minor(new double[][] {}, 0, -1);
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // skipRow is equal to number of rows
         test = Matrices.minor(new double[][] {}, 0, 0);
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix a does contain null column
         test = Matrices.minor(new double[][] { null }, 0, 0);
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix a row width does not remain constant
         test = Matrices.minor(new double[][] { { 0, 0 }, { 0 } }, 0, 0);
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // skipCol is equal to number of columns
         test = Matrices.minor(new double[][] { { 0, 0 }, { 0, 0 } }, 0, 2);
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }
   }

   private void testcaseMinor(double[][] a, double[][] expected, int skipRow,
      int skipCol) {
      double[][] actual = Matrices.minor(a, skipRow, skipCol);

      assertEquals(expected.length, actual.length);
      for (int row = 0; row < expected.length; row++) {
         assertEquals(expected[row].length, actual[row].length);
         for (int col = 0; col < expected[row].length; col++)
            assertEquals(expected[row][col], actual[row][col], DELTA);
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
      double[][] expected;
      double[][] test = null;

      // first test on 2 by 3 and 3 by 2 arrays
      a = new double[][] { { 1, 2, 3 }, { 4, 5, 6 } };
      b = new double[][] { { 7, 8 }, { 9, 10 }, { 11, 12 } };
      expected = new double[][] { { 58, 64 }, { 139, 154 } };
      testcaseMultiply(a, b, expected);

      // second test on 1 by 3 and 3 by 4 arrays
      a = new double[][] { { 3, 4, 2 } };
      b = new double[][] { { 13, 9, 7, 15 }, { 8, 7, 4, 6 }, { 6, 4, 0, 3 } };
      expected = new double[][] { { 83, 63, 37, 75 } };
      testcaseMultiply(a, b, expected);

      // third test on 4 by 4 and 4 by 4 arrays
      a = new double[][] { { 0, 1, 0, 3 }, { 3, 1, 8, 9 }, { 7, 3, 6, 8 },
         { 4, 8, 3, 1 } };
      b = new double[][] { { 0, 3, 7, 4 }, { 1, 1, 3, 8 }, { 0, 8, 6, 3 },
         { 3, 9, 8, 1 } };
      expected = new double[][] { { 10, 28, 27, 11 }, { 28, 155, 144, 53 },
         { 27, 144, 158, 78 }, { 11, 53, 78, 90 } };
      testcaseMultiply(a, b, expected);

      try {
         // matrix a is null
         test = Matrices.multiply(null, new double[][] {});
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix b is null
         test = Matrices.multiply(new double[][] {}, null);
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix a does contain null column
         test =
            Matrices.multiply(new double[][] { null }, new double[][] { null });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix b does contain null column
         test =
            Matrices.multiply(new double[][] { {} }, new double[][] { null });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix a row width does not remain constant
         test = Matrices.multiply(new double[][] { { 0, 0 }, { 0 } },
            new double[][] { {}, {} });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix b number of rows does not equal matrix a number of columns
         test = Matrices.multiply(new double[][] { { 0, 0, 0 }, { 0, 0, 0 } },
            new double[][] { { 0 }, { 0 } });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }
   }

   private void testcaseMultiply(double[][] a, double[][] b,
      double[][] expected) {
      double[][] actual = Matrices.multiply(a, b);

      assertEquals(expected.length, actual.length);
      for (int row = 0; row < expected.length; row++) {
         assertEquals(expected[row].length, actual[row].length);
         for (int col = 0; col < expected[row].length; col++)
            assertEquals(expected[row][col], actual[row][col], DELTA);
      }
   }

   /**
    * Test method for {@link Matrices#multiplyATA(double[][])}.
    */
   @Test
   public void testMultiplyATA() {
      // variables required for testing
      double[][] a;
      double[][] expected;
      double[][] test = null;

      // first test on 2 by 3 array
      a = new double[][] { { 1, 2, 3 }, { 4, 5, 6 } };
      expected =
         new double[][] { { 17, 22, 27 }, { 22, 29, 36 }, { 27, 36, 45 } };
      testcaseMultiplyATA(a, expected);

      // second test on 4 by 3 array
      a = new double[][] { { 11, 8, 2 }, { 13, 4, -2 }, { 3, 12, -6 },
         { -5, -8, 10 } };
      expected = new double[][] { { 324, 216, -72 }, { 216, 288, -144 },
         { -72, -144, 144 } };
      testcaseMultiplyATA(a, expected);

      // third test on 4 by 4 array
      a = new double[][] { { 0, 3, 7, 4 }, { 1, 1, 3, 8 }, { 0, 8, 6, 3 },
         { 3, 9, 8, 1 } };
      expected = new double[][] { { 10, 28, 27, 11 }, { 28, 155, 144, 53 },
         { 27, 144, 158, 78 }, { 11, 53, 78, 90 } };
      testcaseMultiplyATA(a, expected);

      try {
         // matrix a is null
         test = Matrices.multiplyATA(null);
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix a does contain null column
         test = Matrices.multiplyATA(new double[][] { null });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix a row width does not remain constant
         test = Matrices.multiplyATA(new double[][] { { 0, 0 }, { 0 } });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

   }

   private void testcaseMultiplyATA(double[][] a, double[][] expected) {
      double[][] actual = Matrices.multiplyATA(a);

      assertEquals(expected.length, actual.length);
      for (int row = 0; row < expected.length; row++) {
         assertEquals(expected[row].length, actual[row].length);
         for (int col = 0; col < expected[row].length; col++)
            assertEquals(expected[row][col], actual[row][col], DELTA);
      }
   }

   /**
    * Test method for {@link Matrices#transpose(double[][])}.
    */
   @Test
   public void testTranspose() {
      // variables required for testing
      double[][] a;
      double[][] expected;
      double[][] test = null;

      // first test on 1 by 1 array
      a = new double[][] { { 0 } };
      expected = new double[][] { { 0 } };
      testcaseTranspose(a, expected);

      // second test on 2 by 2 array
      a = new double[][] { { 0, 1 }, { 2, 3 } };
      expected = new double[][] { { 0, 2 }, { 1, 3 } };
      testcaseTranspose(a, expected);

      // third test on 3 by 3 array
      a = new double[][] { { 3, 0, 2 }, { 2, 0, -2 }, { 0, 1, 1 } };
      expected = new double[][] { { 3, 2, 0 }, { 0, 0, 1 }, { 2, -2, 1 } };
      testcaseTranspose(a, expected);

      // fourth test on 4 by 6 array
      a = new double[][] { { 1, 2, 3, 4, 5, 6 }, { 1, 2, 3, 4, 5, 6 },
         { 1, 2, 3, 4, 5, 6 }, { 1, 2, 3, 4, 5, 6 } };
      expected = new double[][] { { 1, 1, 1, 1 }, { 2, 2, 2, 2 },
         { 3, 3, 3, 3 }, { 4, 4, 4, 4 }, { 5, 5, 5, 5 }, { 6, 6, 6, 6 } };
      testcaseTranspose(a, expected);

      try {
         // matrix a is null
         test = Matrices.transpose(null);
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix a does contain null column
         test = Matrices.transpose(new double[][] { null });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix a row width does not remain constant
         test = Matrices.transpose(new double[][] { { 0, 0 }, { 0 } });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }
   }

   private void testcaseTranspose(double[][] a, double[][] expected) {
      double[][] actual = Matrices.transpose(a);

      assertEquals(expected.length, actual.length);
      for (int row = 0; row < expected.length; row++) {
         assertEquals(expected[row].length, actual[row].length);
         for (int col = 0; col < expected[row].length; col++)
            assertEquals(expected[row][col], actual[row][col], DELTA);
      }
   }

   /**
    * Test method for {@link Matrices#sofu(double[][])}.
    */
   @Test
   public void testSofu() {
      // variables required for testing
      double[][] a;
      double test = 0;

      // first test on 2 by 1 array
      a = new double[][] { { 3 }, { -4 } };
      testcaseSofu(a, 5);

      // second test on 3 by 2 array
      a = new double[][] { { -2, -4 }, { 3, 1 }, { 2, -1 } };
      testcaseSofu(a, 15);

      // third test on 4 by 5 array
      a = new double[][] { { -2, 2, -4, 1 }, { -1, -4, -1, -1 },
         { 5, 4, -4, 4 }, { -2, 1, -2, 1 }, { -3, 3, 5, -5 } };
      testcaseSofu(a, 526.4979);

      try {
         // matrix a is null
         test = Matrices.sofu(null);
         fail();
      } catch (IllegalArgumentException ex) {
         assertEquals(0, test, DELTA);
      }

      try {
         // matrix a does contain null column
         test = Matrices.sofu(new double[][] { null });
         fail();
      } catch (IllegalArgumentException ex) {
         assertEquals(0, test, DELTA);
      }

      try {
         // matrix a row width does not remain constant
         test = Matrices.sofu(new double[][] { { 0, 0 }, { 0 } });
         fail();
      } catch (IllegalArgumentException ex) {
         assertEquals(0, test, DELTA);
      }
   }

   private void testcaseSofu(double[][] a, double expected) {
      double actual = Matrices.sofu(a);

      assertEquals(expected, actual, DELTA);
   }

   /**
    * Test method for {@link Matrices#faddeev(double[][])}.
    */
   @Test
   public void testFaddeev() {
      // variables required for testing
      double[][] a;
      double[] test = null;

      // first test on 0 by 0 array
      a = new double[][] {};
      testcaseFaddeev(a, new double[] { 1 });

      // second test on 1 by 1 array
      a = new double[][] { { 5 } };
      testcaseFaddeev(a, new double[] { 1, -5 });

      // third test on 2 by 2 array
      a = new double[][] { { -1, 6 }, { -2, 6 } };
      testcaseFaddeev(a, new double[] { 1, -5, 6 });

      // fourth test on 3 by 3 array
      a = new double[][] { { 2, -1, 1 }, { -1, 2, 1 }, { 1, -1, 2 } };
      testcaseFaddeev(a, new double[] { 1, -6, 11, -6 });

      // fifth test on 4 by 4 array
      a = new double[][] { { 8, -1, 3, -1 }, { -1, 6, 2, 0 }, { 3, 2, 9, 1 },
         { -1, 0, 1, 7 } };
      testcaseFaddeev(a, new double[] { 1, -30, 319, -1410, 2138 });

      // sixth test on 3 by 3 array
      a = new double[][] { { 2, 3, -3 }, { -2, -1, 2 }, { 2, 4, -3 } };
      testcaseFaddeev(a, new double[] { 1, 2, -1, -2 });

      try {
         // array a is null
         test = Matrices.faddeev(null);
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix a does contain null column
         test = Matrices.faddeev(new double[][] { null });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }

      try {
         // matrix a must does not have as many columns as rows
         test =
            Matrices.faddeev(new double[][] { { 0, 0 }, { 0, 0 }, { 0, 0 } });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(test);
      }
   }

   private void testcaseFaddeev(double[][] a, double[] expected) {
      double[] actual = Matrices.faddeev(a);

      assertEquals(expected.length, actual.length);
      for (int row = 0; row < expected.length; row++)
         assertEquals(expected[row], actual[row], DELTA);
   }

}