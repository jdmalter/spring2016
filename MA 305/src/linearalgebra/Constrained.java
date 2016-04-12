package linearalgebra;

/**
 * Contains functions which compute overconstrained and underconstrained systems
 * through multiple implementations.
 * 
 * @author Jacob Malter
 *
 */
public class Constrained {

   /**
    * Does nothing. Never called.
    */
   private Constrained() {
   }

   /**
    * Assumes row by column array representation of matrix. Given that x is the
    * first element of the resulting array and e is the second element of the
    * resulting array, a*x - b = e where the magnitude of e is minimized.
    * Returns array of array of n by 1 and m by 1 array of double,
    * 
    * Given that at is tranpose of given matrix a, this implementation finds the
    * pseudo-inverse of matrix a by solving ((at*a)^-1)*at.
    * 
    * Throws IllegalArgumentException if any of the following is true:
    * {@code (null == a)}, {@code (null == b)}, {@code (a.length != b.length)}
    * where length is number of rows, {@code (null == a[row])} where a[row] is
    * any row, {@code (null == b[row])} where b[row] is any row,
    * {@code (a.length <= acol)} where length is number rows and acol is number
    * of columns on first row, {@code (acol != a[row].length)} where acol is
    * number of columns on first row and a[row].length is number of columns, or
    * {@code (1 != b[row].length)} where b[row].length is number of columns on
    * any row.
    * 
    * There is a special condition where divide by 0 occurs when the determinant
    * of the resulting matrix from at*a equals 0.
    * 
    * @param a
    *        m by n array of double where m > n
    * @param b
    *        m by 1 array of double
    * @return array of n by 1 and m by 1 array of double
    */
   public static double[][][] overconstrainedPI(double[][] a, double[][] b) {
      if (null == a)
         throw new IllegalArgumentException("matrix a must not be null");

      else if (null == b)
         throw new IllegalArgumentException("matrix b must not be null");

      else if (a.length != b.length) throw new IllegalArgumentException(
         "matrix a and b number of rows must be equal");

      int acol = -1;

      for (int row = 0; row < a.length; row++)
         if (null == a[row]) throw new IllegalArgumentException(
            "matrix a must not contain null column");

         else if (null == b[row]) throw new IllegalArgumentException(
            "matrix b must not contain null column");

         else if (row == 0) {
            acol = a[0].length;
            if (a.length <= acol) throw new IllegalArgumentException(
               "matrix a must contain more rows than columns");

         } else if (acol != a[row].length) throw new IllegalArgumentException(
            "matrix a row width must remain constant");

         else if (1 != b[row].length) throw new IllegalArgumentException(
            "matrix b must only contain one column");

      double[][] x = Matrices.multiply(
         Matrices.multiply(Fundamentals.invert(Matrices.multiplyATA(a)),
            Matrices.transpose(a)),
         b);

      // inline matrix multiply on vector

      double[][] e = new double[a.length][1];

      // find dot products of every row and column combination

      for (int row = 0; row < e.length; row++) {

         // inline dot product

         double sum = 0;

         for (int brow = 0; brow < x.length; brow++)
            sum += a[row][brow] * x[brow][0];

         // change to matrix multiply to satisfy specifications

         e[row][0] = sum - b[row][0];
      }

      return new double[][][] { x, e };
   }

   /**
    * Assumes row by column array representation of matrix. Given that x is
    * resulting matrix, a*x - b = 0 where the magnitude of x is minimized.
    * Returns n by 1 array of double.
    * 
    * Given that at is tranpose of given matrix a, this implementation finds the
    * pseudo-inverse of matrix a by solving at*((a*at)^-1).
    * 
    * Throws IllegalArgumentException if any of the following is true:
    * {@code (null == a)}, {@code (null == b)}, {@code (a.length != b.length)}
    * where length is number of rows, {@code (null == a[row])} where a[row] is
    * any row, {@code (null == b[row])} where b[row] is any row,
    * {@code (a.length >= acol)} where length is number rows and acol is number
    * of columns on first row, {@code (acol != a[row].length)} where acol is
    * number of columns on first row and a[row].length is number of columns, or
    * {@code (1 != b[row].length)} where b[row].length is number of columns on
    * any row.
    * 
    * There is a special condition where divide by 0 occurs when the determinant
    * of the resulting matrix from a*at equal 0.
    * 
    * @param a
    *        m by n array of double where m < n
    * @param b
    *        m by 1 array of double
    * @return n by 1 array of double
    */
   public static double[][] underconstrainedPI(double[][] a, double[][] b) {
      if (null == a)
         throw new IllegalArgumentException("matrix a must not be null");

      else if (null == b)
         throw new IllegalArgumentException("matrix b must not be null");

      else if (a.length != b.length) throw new IllegalArgumentException(
         "matrix a and b number of rows must be equal");

      int acol = -1;

      for (int row = 0; row < a.length; row++)
         if (null == a[row]) throw new IllegalArgumentException(
            "matrix a must not contain null column");

         else if (null == b[row]) throw new IllegalArgumentException(
            "matrix b must not contain null column");

         else if (row == 0) {
            acol = a[0].length;
            if (a.length >= acol) throw new IllegalArgumentException(
               "matrix a must contain less rows than columns");

         } else if (acol != a[row].length) throw new IllegalArgumentException(
            "matrix a row width must remain constant");

         else if (1 != b[row].length) throw new IllegalArgumentException(
            "matrix b must only contain one column");

      double[][] at = Matrices.transpose(a);
      return Matrices.multiply(
         Matrices.multiply(at, Fundamentals.invert(Matrices.multiplyATA(at))),
         b);
   }

   /**
    * Assumes row by column array representation of matrix. Given that x is the
    * first element of the resulting array and e is the second element of the
    * resulting array, a*x - b = e where the magnitude of e is minimized.
    * Returns array of array of n by 1 and m by 1 array of double,
    * 
    * This implementation finds x in two steps. It is given that at is tranpose
    * of given matrix a, l is the cholesky decomposition of a, and lt is
    * tranpose of l. First, y is found by solving l*y = at*b. Second, x is found
    * by solving lt*x = y.
    * 
    * Throws IllegalArgumentException if any of the following is true:
    * {@code (null == a)}, {@code (null == b)}, {@code (a.length != b.length)}
    * where length is number of rows, {@code (null == a[row])} where a[row] is
    * any row, {@code (null == b[row])} where b[row] is any row,
    * {@code (a.length <= acol)} where length is number rows and acol is number
    * of columns on first row, {@code (acol != a[row].length)} where acol is
    * number of columns on first row and a[row].length is number of columns, or
    * {@code (1 != b[row].length)} where b[row].length is number of columns on
    * any row.
    * 
    * There is a special condition where divide by 0 occurs when the determinant
    * of the resulting matrix from at*a equals 0.
    * 
    * @param a
    *        m by n array of double where m > n
    * @param b
    *        m by 1 array of double
    * @return array of n by 1 and m by 1 array of double
    */
   public static double[][][] overconstrainedCD(double[][] a, double[][] b) {
      if (null == a)
         throw new IllegalArgumentException("matrix a must not be null");

      else if (null == b)
         throw new IllegalArgumentException("matrix b must not be null");

      else if (a.length != b.length) throw new IllegalArgumentException(
         "matrix a and b number of rows must be equal");

      int acol = -1;

      for (int row = 0; row < a.length; row++)
         if (null == a[row]) throw new IllegalArgumentException(
            "matrix a must not contain null column");

         else if (null == b[row]) throw new IllegalArgumentException(
            "matrix b must not contain null column");

         else if (row == 0) {
            acol = a[0].length;
            if (a.length <= acol) throw new IllegalArgumentException(
               "matrix a must contain more rows than columns");

         } else if (acol != a[row].length) throw new IllegalArgumentException(
            "matrix a row width must remain constant");

         else if (1 != b[row].length) throw new IllegalArgumentException(
            "matrix b must only contain one column");

      double[][] l = Decomposition.cholesky(Matrices.multiplyATA(a));

      double[][] x =
         Fundamentals.gaussJordan(Matrices.transpose(l), Fundamentals
            .gaussJordan(l, Matrices.multiply(Matrices.transpose(a), b)));

      // inline matrix multiply on vector

      double[][] e = new double[a.length][1];

      // find dot products of every row and column combination

      for (int row = 0; row < e.length; row++) {

         // inline dot product

         double sum = 0;

         for (int brow = 0; brow < x.length; brow++)
            sum += a[row][brow] * x[brow][0];

         // change to matrix multiply to satisfy specifications

         e[row][0] = sum - b[row][0];
      }

      return new double[][][] { x, e };
   }

   /**
    * Assumes row by column array representation of matrix. Given that x is the
    * first element of the resulting array and e is the second element of the
    * resulting array, a*x - b = e where the magnitude of e is minimized.
    * Returns array of array of n by 1 and m by 1 array of double,
    * 
    * This implementation finds x in two steps. First, array qr is the resulting
    * array from qr decomposition on given matrix a where q is the first element
    * of qr and r is the second element of qr. Second, given that qt is tranpose
    * of q, x is found by solving r*x = qt*b.
    * 
    * Throws IllegalArgumentException if any of the following is true:
    * {@code (null == a)}, {@code (null == b)}, {@code (a.length != b.length)}
    * where length is number of rows, {@code (null == a[row])} where a[row] is
    * any row, {@code (null == b[row])} where b[row] is any row,
    * {@code (a.length <= acol)} where length is number rows and acol is number
    * of columns on first row, {@code (acol != a[row].length)} where acol is
    * number of columns on first row and a[row].length is number of columns, or
    * {@code (1 != b[row].length)} where b[row].length is number of columns on
    * any row.
    * 
    * There is a special condition where divide by 0 occurs when the determinant
    * of the resulting matrix from at*a equals 0.
    * 
    * @param a
    *        m by n array of double where m > n
    * @param b
    *        m by 1 array of double
    * @return array of n by 1 and m by 1 array of double
    */
   public static double[][][] overconstrainedQR(double[][] a, double[][] b) {
      if (null == a)
         throw new IllegalArgumentException("matrix a must not be null");

      else if (null == b)
         throw new IllegalArgumentException("matrix b must not be null");

      else if (a.length != b.length) throw new IllegalArgumentException(
         "matrix a and b number of rows must be equal");

      int acol = -1;

      for (int row = 0; row < a.length; row++)
         if (null == a[row]) throw new IllegalArgumentException(
            "matrix a must not contain null column");

         else if (null == b[row]) throw new IllegalArgumentException(
            "matrix b must not contain null column");

         else if (row == 0) {
            acol = a[0].length;
            if (a.length <= acol) throw new IllegalArgumentException(
               "matrix a must contain more rows than columns");

         } else if (acol != a[row].length) throw new IllegalArgumentException(
            "matrix a row width must remain constant");

         else if (1 != b[row].length) throw new IllegalArgumentException(
            "matrix b must only contain one column");

      double[][][] qr = Decomposition.qr(a);
      double[][] q = qr[0];
      double[][] r = qr[1];

      double[][] x = Fundamentals.gaussJordan(r,
         Matrices.multiply(Matrices.transpose(q), b));

      // inline matrix multiply on vector

      double[][] e = new double[a.length][1];

      // find dot products of every row and column combination

      for (int row = 0; row < e.length; row++) {

         // inline dot product

         double sum = 0;

         for (int brow = 0; brow < x.length; brow++)
            sum += a[row][brow] * x[brow][0];

         // change to matrix multiply to satisfy specifications

         e[row][0] = sum - b[row][0];
      }

      return new double[][][] { x, e };
   }

   /**
    * Assumes row by column array representation of matrix. Given that x is the
    * first element of the resulting array and e is the second element of the
    * resulting array, a*x - b = e where the magnitude of e is minimized.
    * Returns array of array of n by 1 and m by 1 array of double,
    * 
    * This implementation finds x in one step. Given that array usv is the
    * resulting array of sv decomposition on given matrix a where ut is tranpose
    * of the first element of usv, (s^-1) is the inverse of the second element
    * of usv, v is the third element of usv, x is equal to v*(s^-1)*ut*b.
    * 
    * Throws IllegalArgumentException if any of the following is true:
    * {@code (null == a)}, {@code (null == b)}, {@code (a.length != b.length)}
    * where length is number of rows, {@code (null == a[row])} where a[row] is
    * any row, {@code (null == b[row])} where b[row] is any row,
    * {@code (a.length <= acol)} where length is number rows and acol is number
    * of columns on first row, {@code (acol != a[row].length)} where acol is
    * number of columns on first row and a[row].length is number of columns, or
    * {@code (1 != b[row].length)} where b[row].length is number of columns on
    * any row.
    * 
    * There is a special condition where divide by 0 occurs when the determinant
    * of the resulting matrix from at*a equals 0.
    * 
    * @param a
    *        m by n array of double where m > n
    * @param b
    *        m by 1 array of double
    * @return array of n by 1 and m by 1 array of double
    */
   public static double[][][] overconstrainedSV(double[][] a, double[][] b) {
      if (null == a)
         throw new IllegalArgumentException("matrix a must not be null");

      else if (null == b)
         throw new IllegalArgumentException("matrix b must not be null");

      else if (a.length != b.length) throw new IllegalArgumentException(
         "matrix a and b number of rows must be equal");

      int acol = -1;

      for (int row = 0; row < a.length; row++)
         if (null == a[row]) throw new IllegalArgumentException(
            "matrix a must not contain null column");

         else if (null == b[row]) throw new IllegalArgumentException(
            "matrix b must not contain null column");

         else if (row == 0) {
            acol = a[0].length;
            if (a.length <= acol) throw new IllegalArgumentException(
               "matrix a must contain more rows than columns");

         } else if (acol != a[row].length) throw new IllegalArgumentException(
            "matrix a row width must remain constant");

         else if (1 != b[row].length) throw new IllegalArgumentException(
            "matrix b must only contain one column");

      double[][][] usv = Decomposition.singularValue(a);
      double[][] u = usv[0];
      double[][] s = usv[1];
      double[][] v = usv[2];

      // simple inverse because sigma is diagonal array
      double[][] inverse = new double[acol][acol];
      for (int row = 0; row < acol; row++)
         inverse[row][row] = 1 / s[row][row];

      double[][] x = Matrices.multiply(Matrices
         .multiply(Matrices.multiply(v, inverse), Matrices.transpose(u)), b);

      // inline matrix multiply on vector

      double[][] e = new double[a.length][1];

      // find dot products of every row and column combination

      for (int row = 0; row < e.length; row++) {

         // inline dot product

         double sum = 0;

         for (int brow = 0; brow < x.length; brow++)
            sum += a[row][brow] * x[brow][0];

         // change to matrix multiply to satisfy specifications

         e[row][0] = sum - b[row][0];
      }

      return new double[][][] { x, e };
   }

   /**
    * Assumes row by column array representation of matrix. Given that x is
    * resulting matrix, a*x - b = 0 where the magnitude of x is minimized.
    * Returns n by 1 array of double.
    * 
    * This implementation finds x in one step. Given that array usv is the
    * resulting array of sv decomposition on tranpose of the given matrix a
    * where u is the first element of usv, (s^-1) is the inverse of the second
    * element of usv, vt is tranpose of the third element of usv, x is found by
    * solving u*(s^-1)*vt*b.
    * 
    * Throws IllegalArgumentException if any of the following is true:
    * {@code (null == a)}, {@code (null == b)}, {@code (a.length != b.length)}
    * where length is number of rows, {@code (null == a[row])} where a[row] is
    * any row, {@code (null == b[row])} where b[row] is any row,
    * {@code (a.length >= acol)} where length is number rows and acol is number
    * of columns on first row, {@code (acol != a[row].length)} where acol is
    * number of columns on first row and a[row].length is number of columns, or
    * {@code (1 != b[row].length)} where b[row].length is number of columns on
    * any row.
    * 
    * There is a special condition where divide by 0 occurs when the determinant
    * of the resulting matrix from a*at equals 0.
    * 
    * @param a
    *        m by n array of double where m < n
    * @param b
    *        m by 1 array of double
    * @return n by 1 array of double
    */
   public static double[][] underconstrainedSV(double[][] a, double[][] b) {
      if (null == a)
         throw new IllegalArgumentException("matrix a must not be null");

      else if (null == b)
         throw new IllegalArgumentException("matrix b must not be null");

      else if (a.length != b.length) throw new IllegalArgumentException(
         "matrix a and b number of rows must be equal");

      int acol = -1;

      for (int row = 0; row < a.length; row++)
         if (null == a[row]) throw new IllegalArgumentException(
            "matrix a must not contain null column");

         else if (null == b[row]) throw new IllegalArgumentException(
            "matrix b must not contain null column");

         else if (row == 0) {
            acol = a[0].length;
            if (a.length >= acol) throw new IllegalArgumentException(
               "matrix a must contain less rows than columns");

         } else if (acol != a[row].length) throw new IllegalArgumentException(
            "matrix a row width must remain constant");

         else if (1 != b[row].length) throw new IllegalArgumentException(
            "matrix b must only contain one column");

      double[][][] usv = Decomposition.singularValue(Matrices.transpose(a));
      double[][] u = usv[0];
      double[][] s = usv[1];
      double[][] v = usv[2];

      // simple inverse because sigma is diagonal array
      double[][] inverse = new double[a.length][a.length];
      for (int row = 0; row < a.length; row++)
         inverse[row][row] = 1 / s[row][row];

      return Matrices.multiply(Matrices.multiply(Matrices.multiply(u, inverse),
         Matrices.transpose(v)), b);
   }

}