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
      a = new double[][] { new double[] { -3, 1, 2 }, new double[] { 9, -4, -8 },
         new double[] { -9, 5, 8 } };
      b = new double[][] { new double[] { -3 }, new double[] { 3 }, new double[] { 1 } };
      expectedAB = new double[][] { new double[] { -3, 1, 2, -3 }, new double[] { 9, -4, -8, 3 },
         new double[] { -9, 5, 8, 1 } };
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
      expectedAB = new double[][] { new double[] { -3, 2, -5 }, new double[] { 3, 1, 4 } };
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
         testAB = Fundamentals.augment(new double[][] { new double[] {} }, new double[][] {});
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(testAB);
      }

      try {
         // matrix a does contain null column
         testAB = Fundamentals.augment(new double[][] { null }, new double[][] { null });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(testAB);
      }

      try {
         // matrix b does contain null column
         testAB = Fundamentals.augment(new double[][] { new double[] {} }, new double[][] { null });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(testAB);
      }

      try {
         // matrix a row width does not remain constant
         testAB = Fundamentals.augment(new double[][] { new double[] { 0, 0 }, new double[] { 0 } },
            new double[][] { new double[] {}, new double[] {} });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(testAB);
      }

      try {
         // matrix b row width does not remain constant
         testAB = Fundamentals.augment(new double[][] { new double[] { 0 }, new double[] { 0 } },
            new double[][] { new double[] { 0, 0 }, new double[] { 0 } });
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
      a = new double[][] { new double[] { 3, -6, 6, -6 }, new double[] { -3, 9, -8, 0 },
         new double[] { 3, -9, 6, 6 } };
      expectedFor = new double[][] { new double[] { 3, -6, 6, -6 }, new double[] { 0, 3, -2, -6 },
         new double[] { 0, 0, -2, 6 } };
      actualFor = Fundamentals.forwardEliminate(a);

      assertEquals(expectedFor.length, actualFor.length);
      for (int row = 0; row < expectedFor.length; row++) {
         assertEquals(expectedFor[row].length, actualFor[row].length);
         for (int col = 0; col < expectedFor[row].length; col++)
            assertEquals(expectedFor[row][col], actualFor[row][col], DELTA);
      }

      // second test on 3 by 4 matrix
      a = new double[][] { new double[] { 4, -6, 8, -2 }, new double[] { -4, 3, -2, -1 },
         new double[] { 4, 0, -2, 0 } };
      expectedFor = new double[][] { new double[] { 4, -6, 8, -2 }, new double[] { 0, -3, 6, -3 },
         new double[] { 0, 0, 2, -4 } };
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
         testFor = Fundamentals.forwardEliminate(
            new double[][] { new double[] { 0, 0 }, new double[] { 0, 0 }, new double[] { 0, 0 } });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(testFor);
      }

      try {
         // matrix a row width does not remain constant
         testFor = Fundamentals
            .forwardEliminate(new double[][] { new double[] { 0, 0 }, new double[] { 0 } });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(testFor);
      }

      try {
         // matrix a row width does not remain constant
         testFor = Fundamentals
            .forwardEliminate(new double[][] { new double[] { 0, 0 }, new double[] { 0 } });
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
      a = new double[][] { new double[] { 3, -6, 6, -6 }, new double[] { 0, 3, -2, -6 },
         new double[] { 0, 0, -2, 6 } };
      expectedBack = new double[][] { new double[] { 3, 0, 0, -12 }, new double[] { 0, 3, 0, -12 },
         new double[] { 0, 0, -2, 6 } };
      actualBack = Fundamentals.backwardEliminate(a);

      assertEquals(expectedBack.length, actualBack.length);
      for (int row = 0; row < expectedBack.length; row++) {
         assertEquals(expectedBack[row].length, actualBack[row].length);
         for (int col = 0; col < expectedBack[row].length; col++)
            assertEquals(expectedBack[row][col], actualBack[row][col], DELTA);
      }

      // second test on 3 by 4 matrix
      a = new double[][] { new double[] { 4, -6, 8, -2 }, new double[] { 0, -3, 6, -3 },
         new double[] { 0, 0, 2, -4 } };
      expectedBack = new double[][] { new double[] { 4, 0, 0, -4 }, new double[] { 0, -3, 0, 9 },
         new double[] { 0, 0, 2, -4 } };
      actualBack = Fundamentals.backwardEliminate(a);

      assertEquals(expectedBack.length, actualBack.length);
      for (int row = 0; row < expectedBack.length; row++) {
         assertEquals(expectedBack[row].length, actualBack[row].length);
         for (int col = 0; col < expectedBack[row].length; col++)
            assertEquals(expectedBack[row][col], actualBack[row][col], DELTA);
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
         testBack = Fundamentals.backwardEliminate(
            new double[][] { new double[] { 0, 0 }, new double[] { 0, 0 }, new double[] { 0, 0 } });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(testBack);
      }

      try {
         // matrix a row width does not remain constant
         testBack = Fundamentals
            .backwardEliminate(new double[][] { new double[] { 0, 0 }, new double[] { 0 } });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(testBack);
      }

      try {
         // matrix a row width does not remain constant
         testBack = Fundamentals
            .backwardEliminate(new double[][] { new double[] { 0, 0 }, new double[] { 0 } });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(testBack);
      }
   }

   /**
    * Test method for {@link Fundamentals#rref(double[][])}.
    */
   @Test
   public void testRref() {
      double[][] a;
      double[][] expectedRref;
      double[][] actualRref;
      double[][] testRref = null;

      // first test on 1 by 1 array
      a = new double[][] { new double[] { 7 } };
      expectedRref = new double[][] { new double[] { 1 } };
      actualRref = Fundamentals.rref(a);

      assertEquals(expectedRref.length, actualRref.length);
      for (int row = 0; row < expectedRref.length; row++) {
         assertEquals(expectedRref[row].length, actualRref[row].length);
         for (int col = 0; col < expectedRref[row].length; col++)
            assertEquals(expectedRref[row][col], actualRref[row][col], DELTA);
      }

      // second test on 2 by 2 array
      a = new double[][] { new double[] { 10, -60 }, new double[] { 15, 60 } };
      expectedRref = new double[][] { new double[] { 1, -6 }, new double[] { 0, 0 } };
      actualRref = Fundamentals.rref(a);

      assertEquals(expectedRref.length, actualRref.length);
      for (int row = 0; row < expectedRref.length; row++) {
         assertEquals(expectedRref[row].length, actualRref[row].length);
         for (int col = 0; col < expectedRref[row].length; col++)
            assertEquals(expectedRref[row][col], actualRref[row][col], DELTA);
      }

      // third test on 3 by 3 array
      a = new double[][] { new double[] { -1, 4, -4 }, new double[] { -1, 0, -2 },
         new double[] { -3, 0, -4 } };
      expectedRref = new double[][] { new double[] { 1, 0, 2 }, new double[] { 0, 1, -0.5 },
         new double[] { 0, 0, 0 } };
      actualRref = Fundamentals.rref(a);

      assertEquals(expectedRref.length, actualRref.length);
      for (int row = 0; row < expectedRref.length; row++) {
         assertEquals(expectedRref[row].length, actualRref[row].length);
         for (int col = 0; col < expectedRref[row].length; col++)
            assertEquals(expectedRref[row][col], actualRref[row][col], DELTA);
      }

      // fourth test on 3 by 4 array
      a = new double[][] { new double[] { 3, -6, 6, -6 }, new double[] { -3, 9, -8, 0 },
         new double[] { 3, -9, 6, 6 } };
      expectedRref = new double[][] { new double[] { 1, 0, 0, -4 }, new double[] { 0, 1, 0, -4 },
         new double[] { 0, 0, 1, -3 } };
      actualRref = Fundamentals.rref(a);

      assertEquals(expectedRref.length, actualRref.length);
      for (int row = 0; row < expectedRref.length; row++) {
         assertEquals(expectedRref[row].length, actualRref[row].length);
         for (int col = 0; col < expectedRref[row].length; col++)
            assertEquals(expectedRref[row][col], actualRref[row][col], DELTA);
      }

      try {
         // matrix a is null
         testRref = Fundamentals.rref(null);
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(testRref);
      }

      try {
         // matrix a does contain null column
         testRref = Fundamentals.rref(new double[][] { null });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(testRref);
      }

      try {
         // matrix a must does not have as many or more columns than rows
         testRref = Fundamentals.rref(
            new double[][] { new double[] { 0, 0 }, new double[] { 0, 0 }, new double[] { 0, 0 } });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(testRref);
      }

      try {
         // matrix a row width does not remain constant
         testRref = Fundamentals.rref(new double[][] { new double[] { 0, 0 }, new double[] { 0 } });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(testRref);
      }

      try {
         // matrix a row width does not remain constant
         testRref = Fundamentals.rref(new double[][] { new double[] { 0, 0 }, new double[] { 0 } });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(testRref);
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
      double[][] actualB;
      double[][] expectedX;
      double[][] actualX;
      double[][] testX = null;

      // first test on 3 by 3 and 3 by 1 arrays
      a = new double[][] { new double[] { -3, 1, 2 }, new double[] { 9, -4, -8 },
         new double[] { -9, 5, 8 } };
      b = new double[][] { new double[] { -3 }, new double[] { 3 }, new double[] { 1 } };
      expectedX = new double[][] { new double[] { 3 }, new double[] { 4 }, new double[] { 1 } };
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
      a = new double[][] { new double[] { 4, -6, 8 }, new double[] { -4, 3, -2 },
         new double[] { 4, 0, -2 } };
      b = new double[][] { new double[] { -2 }, new double[] { -1 }, new double[] { 0 } };
      expectedX = new double[][] { new double[] { -1 }, new double[] { -3 }, new double[] { -2 } };
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
         testX = Fundamentals.gaussJordan(new double[][] { new double[] {} }, new double[][] {});
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(testX);
      }

      try {
         // matrix a does contain null column
         testX = Fundamentals.gaussJordan(new double[][] { null }, new double[][] { null });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(testX);
      }

      try {
         // matrix b does contain null column
         testX =
            Fundamentals.gaussJordan(new double[][] { new double[] {} }, new double[][] { null });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(testX);
      }

      try {
         // matrix a row width does not remain constant
         testX =
            Fundamentals.gaussJordan(new double[][] { new double[] { 0, 0 }, new double[] { 0 } },
               new double[][] { new double[] {}, new double[] {} });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(testX);
      }

      try {
         // matrix b row width does not remain constant
         testX = Fundamentals.gaussJordan(new double[][] { new double[] { 0 }, new double[] { 0 } },
            new double[][] { new double[] { 0, 0 }, new double[] { 0 } });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(testX);
      }

      try {
         // no solution exists
         testX = Fundamentals.gaussJordan(
            new double[][] { new double[] { -2, 6, -4 }, new double[] { 2, -4, 2 },
               new double[] { 2, -4, 2 } },
            new double[][] { new double[] { -8 }, new double[] { -6 }, new double[] { -2 } });
         fail();
      } catch (ArithmeticException ex) {
         assertNull(testX);
      }

      try {
         // many solutions exist
         testX = Fundamentals.gaussJordan(
            new double[][] { new double[] { 4, -2, 2 }, new double[] { 4, 0, 4 },
               new double[] { 16, -16, 0 } },
            new double[][] { new double[] { 4 }, new double[] { 8 }, new double[] { 0 } });
         fail();
      } catch (ArithmeticException ex) {
         assertNull(testX);
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
      double[][][] actualABX;
      double[][] expectedX;
      double[][] actualX;
      double[][] expectedAB;
      double[][] actualAB;
      double[][][] testABX = null;

      // first test on 3 by 3 and 3 by 1 arrays
      a = new double[][] { new double[] { -3, 1, 2 }, new double[] { 9, -4, -8 },
         new double[] { -9, 5, 8 } };
      b = new double[][] { new double[] { -3 }, new double[] { 3 }, new double[] { 1 } };
      expectedAB = new double[][] { new double[] { 1, 0, 0 }, new double[] { 0, 1, 0 },
         new double[] { 0, 0, 1 } };
      expectedX = new double[][] { new double[] { 3 }, new double[] { 4 }, new double[] { 1 } };
      actualABX = Fundamentals.gaussJordanMany(a, b);
      actualAB = actualABX[0];
      actualX = actualABX[1];

      assertEquals(expectedAB.length, actualAB.length);
      for (int row = 0; row < expectedAB.length; row++) {
         assertEquals(expectedAB[row].length, actualAB[row].length);
         for (int col = 0; col < expectedAB[row].length; col++)
            assertEquals(expectedAB[row][col], actualAB[row][col], DELTA);
      }

      assertEquals(expectedX.length, actualX.length);
      for (int row = 0; row < expectedX.length; row++) {
         assertEquals(expectedX[row].length, actualX[row].length);
         for (int col = 0; col < expectedX[row].length; col++)
            assertEquals(expectedX[row][col], actualX[row][col], DELTA);
      }

      // second test on 2 by 2 and 2 by 1 arrays
      a = new double[][] { new double[] { 4, -6, 8 }, new double[] { -4, 3, -2 },
         new double[] { 4, 0, -2 } };
      b = new double[][] { new double[] { -2 }, new double[] { -1 }, new double[] { 0 } };
      expectedAB = new double[][] { new double[] { 1, 0, 0 }, new double[] { 0, 1, 0 },
         new double[] { 0, 0, 1 } };
      expectedX = new double[][] { new double[] { -1 }, new double[] { -3 }, new double[] { -2 } };
      actualABX = Fundamentals.gaussJordanMany(a, b);
      actualAB = actualABX[0];
      actualX = actualABX[1];

      assertEquals(expectedAB.length, actualAB.length);
      for (int row = 0; row < expectedAB.length; row++) {
         assertEquals(expectedAB[row].length, actualAB[row].length);
         for (int col = 0; col < expectedAB[row].length; col++)
            assertEquals(expectedAB[row][col], actualAB[row][col], DELTA);
      }

      assertEquals(expectedX.length, actualX.length);
      for (int row = 0; row < expectedX.length; row++) {
         assertEquals(expectedX[row].length, actualX[row].length);
         for (int col = 0; col < expectedX[row].length; col++)
            assertEquals(expectedX[row][col], actualX[row][col], DELTA);
      }

      // third test on 3 by 3 and 3 by 1 arrays
      a = new double[][] { new double[] { 4, -2, 2 }, new double[] { 4, 0, 4 },
         new double[] { 16, -16, 0 } };
      b = new double[][] { new double[] { 4 }, new double[] { 8 }, new double[] { 0 } };
      expectedAB = new double[][] { new double[] { 1, 0, -1 }, new double[] { 0, 1, -1 },
         new double[] { 0, 0, 1 } };
      expectedX = new double[][] { new double[] { 2 }, new double[] { 2 }, new double[] { 0 } };
      actualABX = Fundamentals.gaussJordanMany(a, b);
      actualAB = actualABX[0];
      actualX = actualABX[1];

      assertEquals(expectedAB.length, actualAB.length);
      for (int row = 0; row < expectedAB.length; row++) {
         assertEquals(expectedAB[row].length, actualAB[row].length);
         for (int col = 0; col < expectedAB[row].length; col++)
            assertEquals(expectedAB[row][col], actualAB[row][col], DELTA);
      }

      assertEquals(expectedX.length, actualX.length);
      for (int row = 0; row < expectedX.length; row++) {
         assertEquals(expectedX[row].length, actualX[row].length);
         for (int col = 0; col < expectedX[row].length; col++)
            assertEquals(expectedX[row][col], actualX[row][col], DELTA);
      }

      // fourth test on 3 by 3 and 3 by 1 arrays
      a = new double[][] { new double[] { -2, -6, -4 }, new double[] { 6, 18, 12 },
         new double[] { 2, 6, 4 } };
      b = new double[][] { new double[] { -4 }, new double[] { 12 }, new double[] { 4 } };
      expectedAB = new double[][] { new double[] { 1, -3, -2 }, new double[] { 0, 1, 0 },
         new double[] { 0, 0, 1 } };
      expectedX = new double[][] { new double[] { 2 }, new double[] { 0 }, new double[] { 0 } };
      actualABX = Fundamentals.gaussJordanMany(a, b);
      actualAB = actualABX[0];
      actualX = actualABX[1];

      assertEquals(expectedAB.length, actualAB.length);
      for (int row = 0; row < expectedAB.length; row++) {
         assertEquals(expectedAB[row].length, actualAB[row].length);
         for (int col = 0; col < expectedAB[row].length; col++)
            assertEquals(expectedAB[row][col], actualAB[row][col], DELTA);
      }

      assertEquals(expectedX.length, actualX.length);
      for (int row = 0; row < expectedX.length; row++) {
         assertEquals(expectedX[row].length, actualX[row].length);
         for (int col = 0; col < expectedX[row].length; col++)
            assertEquals(expectedX[row][col], actualX[row][col], DELTA);
      }

      // fifth test on 4 by 4 and 4 by 1 arrays
      a = new double[][] { new double[] { 2, 1, 2, 10 }, new double[] { -2, -2, 2, -12 },
         new double[] { -2, -3, 6, -14 }, new double[] { 4, 0, 12, 16 } };
      b = new double[][] { new double[] { -7 }, new double[] { 10 }, new double[] { 13 },
         new double[] { -8 } };
      expectedAB = new double[][] { new double[] { 1, 0, -3, -4 }, new double[] { 0, 1, 4, -2 },
         new double[] { 0, 0, 1, 0 }, new double[] { 0, 0, 0, 1 } };
      expectedX = new double[][] { new double[] { -2 }, new double[] { -3 }, new double[] { 0 },
         new double[] { 0 } };
      actualABX = Fundamentals.gaussJordanMany(a, b);
      actualAB = actualABX[0];
      actualX = actualABX[1];

      assertEquals(expectedAB.length, actualAB.length);
      for (int row = 0; row < expectedAB.length; row++) {
         assertEquals(expectedAB[row].length, actualAB[row].length);
         for (int col = 0; col < expectedAB[row].length; col++)
            assertEquals(expectedAB[row][col], actualAB[row][col], DELTA);
      }

      assertEquals(expectedX.length, actualX.length);
      for (int row = 0; row < expectedX.length; row++) {
         assertEquals(expectedX[row].length, actualX[row].length);
         for (int col = 0; col < expectedX[row].length; col++)
            assertEquals(expectedX[row][col], actualX[row][col], DELTA);
      }

      // sixth test on 4 by 4 and 4 by 1 arrays
      a = new double[][] { new double[] { 2, 6, 2, 2 }, new double[] { 4, 12, 4, 4 },
         new double[] { -6, -18, -6, -6 }, new double[] { 2, 6, 2, 2 } };
      b = new double[][] { new double[] { -4 }, new double[] { -8 }, new double[] { 12 },
         new double[] { -4 } };
      expectedAB = new double[][] { new double[] { 1, -3, -1, -1 }, new double[] { 0, 1, 0, 0 },
         new double[] { 0, 0, 1, 0 }, new double[] { 0, 0, 0, 1 } };
      expectedX = new double[][] { new double[] { -2 }, new double[] { 0 }, new double[] { 0 },
         new double[] { 0 } };
      actualABX = Fundamentals.gaussJordanMany(a, b);
      actualAB = actualABX[0];
      actualX = actualABX[1];

      assertEquals(expectedAB.length, actualAB.length);
      for (int row = 0; row < expectedAB.length; row++) {
         assertEquals(expectedAB[row].length, actualAB[row].length);
         for (int col = 0; col < expectedAB[row].length; col++)
            assertEquals(expectedAB[row][col], actualAB[row][col], DELTA);
      }

      assertEquals(expectedX.length, actualX.length);
      for (int row = 0; row < expectedX.length; row++) {
         assertEquals(expectedX[row].length, actualX[row].length);
         for (int col = 0; col < expectedX[row].length; col++)
            assertEquals(expectedX[row][col], actualX[row][col], DELTA);
      }

      try {
         // matrix a is null
         testABX = Fundamentals.gaussJordanMany(null, new double[][] {});
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(testABX);
      }

      try {
         // matrix b is null
         testABX = Fundamentals.gaussJordanMany(new double[][] {}, null);
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(testABX);
      }

      try {
         // matrix a and b do not have equal number of rows
         testABX =
            Fundamentals.gaussJordanMany(new double[][] { new double[] {} }, new double[][] {});
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(testABX);
      }

      try {
         // matrix a does contain null column
         testABX = Fundamentals.gaussJordanMany(new double[][] { null }, new double[][] { null });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(testABX);
      }

      try {
         // matrix b does contain null column
         testABX = Fundamentals.gaussJordanMany(new double[][] { new double[] {} },
            new double[][] { null });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(testABX);
      }

      try {
         // matrix a row width does not remain constant
         testABX = Fundamentals.gaussJordanMany(
            new double[][] { new double[] { 0, 0 }, new double[] { 0 } },
            new double[][] { new double[] {}, new double[] {} });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(testABX);
      }

      try {
         // matrix b row width does not remain constant
         testABX =
            Fundamentals.gaussJordanMany(new double[][] { new double[] { 0 }, new double[] { 0 } },
               new double[][] { new double[] { 0, 0 }, new double[] { 0 } });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(testABX);
      }

      try {
         // no solution exists
         testABX = Fundamentals.gaussJordanMany(
            new double[][] { new double[] { -2, 6, -4 }, new double[] { 2, -4, 2 },
               new double[] { 2, -4, 2 } },
            new double[][] { new double[] { -8 }, new double[] { -6 }, new double[] { -2 } });
         fail();
      } catch (ArithmeticException ex) {
         assertNull(testABX);
      }

      try {
         // no solution exists
         testABX = Fundamentals.gaussJordanMany(
            new double[][] { new double[] { 3, 2, -2, 8 }, new double[] { 12, 7, -12, 19 },
               new double[] { -12, -6, 18, 0 }, new double[] { -9, -4, 18, 14 } },
            new double[][] { new double[] { 3 }, new double[] { 13 }, new double[] { -6 },
               new double[] { 2 } });
         fail();
      } catch (ArithmeticException ex) {
         assertNull(testABX);
      }

      try {
         // no solution exists
         testABX = Fundamentals.gaussJordanMany(
            new double[][] { new double[] { -2, 2, 2, -4 }, new double[] { 2, -2, -2, 4 },
               new double[] { 8, -8, -8, 16 }, new double[] { -4, 4, 4, -8 } },
            new double[][] { new double[] { 6 }, new double[] { -10 }, new double[] { -8 },
               new double[] { -4 } });
         fail();
      } catch (ArithmeticException ex) {
         assertNull(testABX);
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
      double[][] actualB;
      double[][] expectedX;
      double[][] actualX;
      double[][] testX = null;

      // first test on 3 by 3 and 3 by 1 arrays
      a = new double[][] { new double[] { -3, 1, 2 }, new double[] { 9, -4, -8 },
         new double[] { -9, 5, 8 } };
      b = new double[][] { new double[] { -3 }, new double[] { 3 }, new double[] { 1 } };
      expectedX = new double[][] { new double[] { 3 }, new double[] { 4 }, new double[] { 1 } };
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
      a = new double[][] { new double[] { 4, -6, 8 }, new double[] { -4, 3, -2 },
         new double[] { 4, 0, -2 } };
      b = new double[][] { new double[] { -2 }, new double[] { -1 }, new double[] { 0 } };
      expectedX = new double[][] { new double[] { -1 }, new double[] { -3 }, new double[] { -2 } };
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
         testX = Fundamentals.cramer(new double[][] { new double[] {} }, new double[][] {});
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(testX);
      }

      try {
         // matrix a does contain null column
         testX = Fundamentals.cramer(new double[][] { null }, new double[][] { null });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(testX);
      }

      try {
         // matrix b does contain null column
         testX = Fundamentals.cramer(new double[][] { new double[] {} }, new double[][] { null });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(testX);
      }

      try {
         // matrix a row width does not remain constant
         testX = Fundamentals.cramer(new double[][] { new double[] { 0, 0 }, new double[] { 0 } },
            new double[][] { new double[] {}, new double[] {} });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(testX);
      }

      try {
         // matrix b row width does not remain constant
         testX = Fundamentals.cramer(new double[][] { new double[] { 0 }, new double[] { 0 } },
            new double[][] { new double[] { 0, 0 }, new double[] { 0 } });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(testX);
      }
   }

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
      expectedInv = new double[][] { new double[] { -5, -4 }, new double[] { -1, -1 } };
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
            if (row == col) assertEquals(1, actualI[row][col], DELTA);

            else assertEquals(0, actualI[row][col], DELTA);

      // second test on 3 by 3 array
      a = new double[][] { new double[] { -1, 4, -4 }, new double[] { -1, 0, -2 },
         new double[] { -3, 0, -4 } };
      expectedInv = new double[][] { new double[] { 0, 2, -1 }, new double[] { 0.25, -1, 0.25 },
         new double[] { 0, -1.5, 0.5 } };
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
            if (row == col) assertEquals(1, actualI[row][col], DELTA);

            else assertEquals(0, actualI[row][col], DELTA);

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
         testInv =
            Fundamentals.invert(new double[][] { new double[] { 0, 0 }, new double[] { 0 } });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(testInv);
      }

      try {
         // divide by 0
         testInv =
            Fundamentals.invert(new double[][] { new double[] { 0, 0 }, new double[] { 0, 0 } });
         fail();
      } catch (ArithmeticException ex) {
         assertNull(testInv);
      }
   }

   /**
    * Test method for {@link Fundamentals#invertFaddeev(double[][])}.
    */
   @Test
   public void testInvertFaddeev() {
      // variables required for testing
      double[][] a;
      double[][] actualI;
      double[][] expectedInv;
      double[][] actualInv;
      double[][] testInv = null;

      // first test on 2 by 2 array
      a = new double[][] { new double[] { -1, 4 }, new double[] { 1, -5 } };
      expectedInv = new double[][] { new double[] { -5, -4 }, new double[] { -1, -1 } };
      actualInv = Fundamentals.invertFaddeev(a);

      assertEquals(expectedInv.length, actualInv.length);
      for (int row = 0; row < expectedInv.length; row++) {
         assertEquals(expectedInv[row].length, actualInv[row].length);
         for (int col = 0; col < expectedInv[row].length; col++)
            assertEquals(expectedInv[row][col], actualInv[row][col], DELTA);
      }

      actualI = Matrices.multiply(a, actualInv);
      for (int row = 0; row < actualI.length; row++)
         for (int col = 0; col < actualI[row].length; col++)
            if (row == col) assertEquals(1, actualI[row][col], DELTA);

            else assertEquals(0, actualI[row][col], DELTA);

      // second test on 3 by 3 array
      a = new double[][] { new double[] { -1, 4, -4 }, new double[] { -1, 0, -2 },
         new double[] { -3, 0, -4 } };
      expectedInv = new double[][] { new double[] { 0, 2, -1 }, new double[] { 0.25, -1, 0.25 },
         new double[] { 0, -1.5, 0.5 } };
      actualInv = Fundamentals.invertFaddeev(a);

      assertEquals(expectedInv.length, actualInv.length);
      for (int row = 0; row < expectedInv.length; row++) {
         assertEquals(expectedInv[row].length, actualInv[row].length);
         for (int col = 0; col < expectedInv[row].length; col++)
            assertEquals(expectedInv[row][col], actualInv[row][col], DELTA);
      }

      actualI = Matrices.multiply(a, actualInv);
      for (int row = 0; row < actualI.length; row++)
         for (int col = 0; col < actualI[row].length; col++)
            if (row == col) assertEquals(1, actualI[row][col], DELTA);

            else assertEquals(0, actualI[row][col], DELTA);

      try {
         // matrix a is null
         testInv = Fundamentals.invertFaddeev(null);
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(testInv);
      }

      try {
         // matrix a does contain null column
         testInv = Fundamentals.invertFaddeev(new double[][] { null });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(testInv);
      }

      try {
         // matrix a row width does not remain constant
         testInv = Fundamentals
            .invertFaddeev(new double[][] { new double[] { 0, 0 }, new double[] { 0 } });
         fail();
      } catch (IllegalArgumentException ex) {
         assertNull(testInv);
      }

      try {
         // divide by 0
         testInv = Fundamentals
            .invertFaddeev(new double[][] { new double[] { 0, 0 }, new double[] { 0, 0 } });
         fail();
      } catch (ArithmeticException ex) {
         assertNull(testInv);
      }
   }

}