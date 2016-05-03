package linearalgebra;

/**
 * Combines complex functions with detailed specifications. Functions take one
 * matrix and return multiple matrices.
 * 
 * @author Jacob Malter
 *
 */
public class Decomposition {

   /**
    * Does nothing. Never called.
    */
   private Decomposition() {
   }

   /**
    * Assumes row by column array representation of matrix. Given that l is the
    * first element of the resulting array and u is the second element of the
    * resulting array, l*u - a = 0. Returns array of n by n and n by n array of
    * double.
    * 
    * Throws IllegalArgumentException if any of the following is true:
    * {@code (null == a)}, {@code (null == a[row])} where a[row] is any row, or
    * {@code (acol != a[row].length)} where acol is number of columns on first
    * row and a[row].length is number of columns.
    * 
    * There is a special condition where ArithmeticException if divide by 0
    * occurs.
    * 
    * @param a
    *        n by n array of double
    * @return array of n by n and n by n array of double
    */
   public static double[][][] lu(double[][] a) {
      if (null == a)
         throw new IllegalArgumentException("matrix a must not be null");

      for (int row = 0; row < a.length; row++)
         if (null == a[row]) throw new IllegalArgumentException(
            "matrix a must not contain null column");

         else if (a.length != a[row].length) throw new IllegalArgumentException(
            "matrix a number of rows and columns must be equal");

      double[][] l = new double[a.length][a.length];
      double[][] u = new double[a.length][a.length];

      // set l to identity matrix for now

      for (int row = 0; row < l.length; row++)
         for (int col = row; col < l.length; col++)
            if (row == col) l[row][col] = 1;

      // copy array a into array u

      for (int row = 0; row < u.length; row++)
         for (int col = 0; col < u[row].length; col++)
            u[row][col] = a[row][col];

      // inline forward elimination

      for (int diag = 0; diag < u.length - 1; diag++) {
         if (u[diag][diag] == 0) throw new ArithmeticException("divide by 0");

         for (int row = diag + 1; row < u.length; row++) {
            // inline is necessary because array l holds quotients
            l[row][diag] = u[row][diag] / u[diag][diag];

            for (int col = 0; col < u[row].length; col++)
               u[row][col] -= (l[row][diag] * u[diag][col]);
         }
      }

      return new double[][][] { l, u };
   }

   /**
    * Assumes row by column array representation of matrix. Given that l is
    * resulting matrix and lt is transpose of resulting matrix, l*lt - a = 0.
    * Returns n by n array of double.
    * 
    * Throws IllegalArgumentException if any of the following is true:
    * {@code (null == a)}, {@code (null == a[row])} where a[row] is any row,
    * {@code (acol != a[row].length)} where acol is number of columns on first
    * row and a[row].length is number of columns, or
    * {@code (a[row][col] != a[col][row])}.
    * 
    * There is a special condition where ArithmeticException if square root of
    * negative value occurs.
    * 
    * There is a special condition where ArithmeticException if divide by 0
    * occurs.
    * 
    * Credit to http://www.sci.utah.edu/~wallstedt/LU.htm for being a good
    * resource.
    * 
    * @param a
    *        n by n array of double
    * @return n by n array of double
    */
   public static double[][] cholesky(double[][] a) {
      if (null == a)
         throw new IllegalArgumentException("matrix a must not be null");

      for (int row = 0; row < a.length; row++)
         if (null == a[row]) throw new IllegalArgumentException(
            "matrix a must not contain null column");

         else if (a.length != a[row].length) throw new IllegalArgumentException(
            "matrix a number of rows and columns must be equal");

      for (int row = 0; row < a.length; row++)
         for (int col = 0; col < a.length; col++)
            if (a[row][col] != a[col][row]) throw new IllegalArgumentException(
               "matrix a must have symmetry across diagonal");

      double[][] l = new double[a.length][a.length];

      for (int row = 0; row < a.length; row++) {
         double sum = 0;

         for (int col = 0; col < row; col++)
            sum += l[row][col] * l[row][col];
         // l[row][row] is not exactly determinant on vector
         l[row][row] = Math.sqrt(a[row][row] - sum);

         for (int lrow = row + 1; lrow < a.length; lrow++) {
            double lsum = 0;

            for (int lcol = 0; lcol < row; ++lcol)
               lsum += l[lrow][lcol] * l[row][lcol];

            if (l[row][row] == 0) throw new ArithmeticException("divide by 0");

            l[lrow][row] = (a[lrow][row] - lsum) / l[row][row];
         }
      }

      return l;
   }

   /**
    * Assumes row by column array representation of matrix. Given that q is the
    * first element of the resulting array and r is the second element of the
    * resulting array, q*r - a = 0. Given that qt is transpose of q, qt*q - I =
    * 0 where I is identity matrix. r is be upper triangular. Returns array of m
    * by n and n by n array of double.
    * 
    * Throws IllegalArgumentException if any of the following is true:
    * {@code (null == a)}, {@code (null == a[row])} where a[row] is any row,
    * {@code (a.length < acol)} where length is number rows and acol is number
    * of columns on first row, or {@code (acol != a[row].length)} where acol is
    * number of columns on first row and a[row].length is number of columns.
    * 
    * @param a
    *        n by m array of double where m >= n
    * @return array of m by n and n by n array of double
    */
   public static double[][][] qr(double[][] a) {
      if (null == a)
         throw new IllegalArgumentException("matrix a must not be null");

      int acol = -1;

      for (int row = 0; row < a.length; row++)
         if (null == a[row]) throw new IllegalArgumentException(
            "matrix a must not contain null column");

         else if (row == 0) {
            acol = a[0].length;
            if (a.length < acol) throw new IllegalArgumentException(
               "matrix a must contain as many or more rows than columns");
         }

         else if (acol != a[row].length) throw new IllegalArgumentException(
            "matrix a row width must remain constant");

      double[][] q = new double[a.length][acol];
      double[][] r = new double[acol][acol];

      for (int col = 0; col < acol; col++) {
         for (int rrow = 0; rrow < col; rrow++)
            for (int row = 0; row < a.length; row++)
               r[rrow][col] += q[row][rrow] * a[row][col];

         double[] w = new double[a.length];
         for (int row = 0; row < a.length; row++) {
            w[row] = a[row][col];

            for (int wrow = 0; wrow < col; wrow++)
               w[row] -= (r[wrow][col] * q[row][wrow]);
         }

         for (int row = 0; row < a.length; row++)
            r[col][col] += w[row] * w[row];
         r[col][col] = Math.sqrt(r[col][col]);

         for (int row = 0; row < a.length; row++)
            q[row][col] = w[row] / r[col][col];
      }

      return new double[][][] { q, r };
   }

   /**
    * Assumes row by column array representation of matrix. Given that l is the
    * first element of the resulting array and v is the second element of the
    * resulting array, a*v - v*l = 0, or given that (v^-1) is the inverse of v,
    * a - v*l*(v^-1) = 0. Returns array of n by n and n by n array of double.
    * 
    * Throws IllegalArgumentException if any of the following is true:
    * {@code (null == a)}, {@code (null == a[row])} where a[row] is any row, or
    * {@code (acol != a[row].length)} where acol is number of columns on first
    * row and a[row].length is number of columns.
    * 
    * Credit to http://www.1728.org/cubic2.htm for solving eigen values on 3 by
    * 3 matrix.
    * 
    * @param a
    *        n by n array of double
    * @return array of n by n and n by n array of double
    */
   public static double[][][] eigen(double[][] a) {
      if (null == a)
         throw new IllegalArgumentException("matrix a must not be null");

      for (int row = 0; row < a.length; row++)
         if (null == a[row]) throw new IllegalArgumentException(
            "matrix a must not contain null column");

         else if (a.length != a[row].length) throw new IllegalArgumentException(
            "matrix a number of rows and columns must be equal");

      double[][] l = new double[a.length][a.length];
      double[][] v = new double[a.length][a.length];

      if (a.length == 0) return new double[][][] { l, v };

      else if (a.length == 1) {

         // x - det(a) = 0
         l[0][0] = a[0][0];
         v[0][0] = 1;
         return new double[][][] { l, v };

      } else if (a.length == 2) {

         // x^2 + bx + det(a) = 0
         double b = a[0][0] + a[1][1];
         double sqrt =
            Math.sqrt(b * b - 4 * (a[0][0] * a[1][1] - a[0][1] * a[1][0]));

         // quadratic formula
         l[0][0] = (b + sqrt) / 2;
         l[1][1] = (b - sqrt) / 2;

         // simplified reduced row echelon form
         for (int row = 0; row < a.length; row++) {
            v[0][row] = -a[0][1] / (a[0][0] - l[row][row]);
            v[1][row] = 1;
         }

         // divide elements in each row by magnitude of vector
         for (int row = 0; row < a.length; row++) {
            double mag = Math.sqrt((v[0][row] * v[0][row]) + 1);
            v[0][row] /= mag;
            v[1][row] /= mag;
         }

         return new double[][][] { l, v };

      } else if (a.length == 3) {

         // find characteristic equation of the form -x^3 + trace(a)*x^2 -
         // (determinant of traces)*x^1 + determinant(a)*x^0
         double[] polynomial = Matrices.faddeev(a);
         double b = -polynomial[1];
         double c = -polynomial[2];
         double d = -polynomial[3];

         double f = (-(3 * c) - (b * b)) / 3;
         double g = (-(2 * b * b * b) - (9 * b * c) - (27 * d)) / 27;
         double h = (g * g / 4) + (f * f * f / 27);
         double i = Math.sqrt((g * g / 4) - h);
         double j = Math.pow(i, 1d / 3);
         double k = Math.acos(-g / (2 * i));

         double m = Math.cos(k / 3);
         double n = Math.sqrt(3) * Math.sin(k / 3);
         double p = (b / 3);
         double q = 2 * j * Math.cos(k / 3) + (b / 3);
         double r = -j * (m + n) + p;
         double s = -j * (m - n) + p;

         if (q > r && q > s) {

            // q is greatest root

            l[0][0] = q;

            if (r > s) {

               // s is least root

               l[1][1] = r;
               l[2][2] = s;

            } else {

               // r is least root

               l[1][1] = s;
               l[2][2] = r;
            }

         } else if (r > s) {

            // r is greatest root

            l[0][0] = r;

            if (q > s) {

               // s is least root

               l[1][1] = q;
               l[2][2] = s;

            } else {

               // q is least root

               l[1][1] = s;
               l[2][2] = q;
            }

         } else {

            // s is greatest root

            l[0][0] = s;

            if (q > r) {

               // r is least root

               l[1][1] = q;
               l[2][2] = r;

            } else {

               // q is least root

               l[1][1] = r;
               l[2][2] = q;
            }
         }

      } else {
         throw new UnsupportedOperationException("Not yet");
      }

      // for each eigen value

      for (int row = 0; row < a.length; row++) {

         // copy array a into array x

         double[][] x = new double[a.length - 1][a.length];
         for (int xrow = 0; xrow < x.length; xrow++)
            for (int xcol = 0; xcol < a.length; xcol++) {
               x[xrow][xcol] = a[xrow][xcol];
               if (xrow == xcol)
                  // if on the diagonal, then subtract lambda
                  x[xrow][xcol] -= l[row][row];
            }

         x = Fundamentals.rref(x);

         // solve for "z" where z is on last row

         v[0][row] = -x[0][2] / x[0][0];
         v[1][row] = -x[1][2] / x[1][1];
         v[x.length][row] = 1;
      }

      // normalize by row

      for (int row = 0; row < a.length; row++) {
         double determinant =
            Math.sqrt((v[0][row] * v[0][row]) + (v[1][row] * v[1][row]) + 1);
         v[0][row] /= determinant;
         v[1][row] /= determinant;
         v[2][row] /= determinant;
      }

      return new double[][][] { l, v };
   }

   /**
    * Assumes row by column array representation of matrix. Given the resulting
    * array, the first element is matrix u, the second element is matrix s, and
    * the third element is matrix v. Given that ut is transpose of u, ut*u - I =
    * 0 where I is identity matrix. Given that vt is transpose of v, vt*v - I =
    * 0 where I is identity matrix. s is diagonal such that elements on diagonal
    * are non-zero and elements not on diagonal equal zero. From the
    * specifications, it is given that matrix a equals u*s*vt. Returns array of
    * m by n, n by n, and n by n array of double.
    * 
    * Throws IllegalArgumentException if any of the following is true:
    * {@code (null == a)}, {@code (null == a[row])} where a[row] is any row,
    * {@code (a.length < acol)} where length is number rows and acol is number
    * of columns on first row, or {@code (acol != a[row].length)} where acol is
    * number of columns on first row and a[row].length is number of columns.
    * 
    * There is a special condition where divide by 0 occurs when l matrix from
    * eigen decomposition of (at*a given at is transpose of a matrix multiplied
    * by a) contains 0 along diagonal.
    * 
    * @param a
    *        m by n array of double where m >= n
    * @return array of m by n, n by n, and n by n array of double
    */
   public static double[][][] singularValue(double[][] a) {
      if (null == a)
         throw new IllegalArgumentException("matrix a must not be null");

      int acol = -1;

      for (int row = 0; row < a.length; row++)
         if (null == a[row]) throw new IllegalArgumentException(
            "matrix a must not contain null column");

         else if (row == 0) {
            acol = a[0].length;
            if (a.length < acol) throw new IllegalArgumentException(
               "matrix a must contain as many or more rows than columns");
         }

         else if (acol != a[row].length) throw new IllegalArgumentException(
            "matrix a row width must remain constant");

      double[][][] lv = eigen(Matrices.multiplyATA(a));
      double[][] l = lv[0];
      double[][] v = lv[1];

      double[][] sigma = new double[acol][acol];

      for (int row = 0; row < acol; row++)
         sigma[row][row] = Math.sqrt(l[row][row]);

      // simple inverse because sigma is diagonal array
      double[][] inverse = new double[acol][acol];
      for (int row = 0; row < acol; row++)
         inverse[row][row] = 1 / sigma[row][row];

      return new double[][][] {
         Matrices.multiply(Matrices.multiply(a, v), inverse), sigma, v };
   }

}