package linearalgebra;

/**
 * Efficiently computes problems.
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
	 * Assumes row by columnn array representation of matrix. Given that l is
	 * the first element of the resulting array and u is the second element of
	 * the resulting array, l*u (through matrix multiplication) equals a.
	 * Returns array of n by n and n by n array of double.
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
	 *            n by n array of double
	 * @return array of n by n and n by n array of double
	 */
	public static double[][][] lu(double[][] a) {
		if (null == a)
			throw new IllegalArgumentException("matrix a must not be null");

		for (int row = 0; row < a.length; row++)
			if (null == a[row])
				throw new IllegalArgumentException(
						"matrix a must not contain null column");

			else if (a.length != a[row].length)
				throw new IllegalArgumentException(
						"matrix a number of rows and columns must be equal");

		double[][] l = new double[a.length][a.length];
		double[][] u = new double[a.length][a.length];

		// set l to identity matrix for now

		for (int row = 0; row < l.length; row++)
			for (int col = row; col < l.length; col++)
				if (row == col)
					l[row][col] = 1;

		// copy array a into array u

		for (int row = 0; row < u.length; row++)
			for (int col = 0; col < u[row].length; col++)
				u[row][col] = a[row][col];

		// inline forward elimination

		for (int diag = 0; diag < u.length - 1; diag++) {
			if (u[diag][diag] == 0)
				throw new ArithmeticException("divide by 0");

			for (int row = diag + 1; row < u.length; row++) {
				l[row][diag] = u[row][diag] / u[diag][diag];

				for (int col = 0; col < u[row].length; col++)
					u[row][col] -= (l[row][diag] * u[diag][col]);
			}
		}

		return new double[][][] { l, u };
	}

	/**
	 * Assumes row by columnn array representation of matrix. Given that l is
	 * resulting matrix and lt is transpose of resulting matrix, l*lt (through
	 * matrix multiplication) equals a. Returns n by n array of double.
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
	 *            n by n array of double
	 * @return n by n array of double
	 */
	public static double[][] cholesky(double[][] a) {
		if (null == a)
			throw new IllegalArgumentException("matrix a must not be null");

		for (int row = 0; row < a.length; row++)
			if (null == a[row])
				throw new IllegalArgumentException(
						"matrix a must not contain null column");

			else if (a.length != a[row].length)
				throw new IllegalArgumentException(
						"matrix a number of rows and columns must be equal");

		for (int row = 0; row < a.length; row++)
			for (int col = 0; col < a.length; col++)
				if (a[row][col] != a[col][row])
					throw new IllegalArgumentException(
							"matrix a must have symmetry across diagonal");

		double[][] l = new double[a.length][a.length];

		for (int row = 0; row < a.length; row++) {
			double sum = 0;

			// inline determinant on vector

			for (int col = 0; col < row; col++)
				sum += l[row][col] * l[row][col];
			l[row][row] = Math.sqrt(a[row][row] - sum);

			for (int lrow = row + 1; lrow < a.length; lrow++) {
				double lsum = 0;

				for (int lcol = 0; lcol < row; ++lcol)
					lsum += l[lrow][lcol] * l[row][lcol];

				if (l[row][row] == 0)
					throw new ArithmeticException("divide by 0");

				l[lrow][row] = (a[lrow][row] - lsum) / l[row][row];
			}
		}

		return l;
	}

	/**
	 * Assumes row by columnn array representation of matrix. Given that q is
	 * the first element of the resulting array and r is the second element of
	 * the resulting array, q*r (through matrix multiplication) equals a. Given
	 * that qt is tranpose of q, qt*q equals I where I is identity matrix. r
	 * must be upper triangular. Returns array of m by n and n by n array of
	 * double.
	 * 
	 * Throws IllegalArgumentException if any of the following is true:
	 * {@code (null == a)}, {@code (null == a[row])} where a[row] is any row,
	 * {@code (a.length < acol)} where length is number rows and acol is number
	 * of columns on first row, or {@code (acol != a[row].length)} where acol is
	 * number of columns on first row and a[row].length is number of columns.
	 * 
	 * @param a
	 *            n by m array of double where m >= n
	 * @return array of m by n and n by n array of double
	 */
	public static double[][][] qr(double[][] a) {
		if (null == a)
			throw new IllegalArgumentException("matrix a must not be null");

		int acol = -1;

		for (int row = 0; row < a.length; row++)
			if (null == a[row])
				throw new IllegalArgumentException(
						"matrix a must not contain null column");

			else if (row == 0) {
				acol = a[0].length;
				if (a.length < acol)
					throw new IllegalArgumentException(
							"matrix a must contain as many or more rows than columns");
			}

			else if (acol != a[row].length)
				throw new IllegalArgumentException(
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

			// inline determinant on vector

			for (int row = 0; row < a.length; row++)
				r[col][col] += w[row] * w[row];
			r[col][col] = Math.sqrt(r[col][col]);

			for (int row = 0; row < a.length; row++)
				q[row][col] = w[row] / r[col][col];
		}

		return new double[][][] { q, r };
	}

	/**
	 * Assumes row by columnn array representation of matrix. Given that l is
	 * the first element of the resulting array and v is the second element of
	 * the resulting array, a*v (through matrix multiplication) equals v*l
	 * (through matrix multiplication). Returns array of n by n and n by n array
	 * of double.
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
	 *            n by n array of double
	 * @return array of n by n and n by n array of double
	 */
	public static double[][][] eigen(double[][] a) {
		if (null == a)
			throw new IllegalArgumentException("matrix a must not be null");

		for (int row = 0; row < a.length; row++)
			if (null == a[row])
				throw new IllegalArgumentException(
						"matrix a must not contain null column");

			else if (a.length != a[row].length)
				throw new IllegalArgumentException(
						"matrix a number of rows and columns must be equal");

		double[][] l = new double[a.length][a.length];
		double[][] v = new double[a.length][a.length];

		if (a.length == 1) {

			// assumed double a = -1;
			double b = a[0][0];

			l[0][0] = b;

			v[0][0] = 1;

		} else if (a.length == 2) {

			// assumed double a = 1;
			double b = -a[0][0] - a[1][1];
			double c = a[0][0] * a[1][1] - a[0][1] * a[1][0];
			double sqrt = Math.sqrt((b * b) - (4 * c));

			l[0][0] = (-b + sqrt) / 2;
			l[1][1] = (-b - sqrt) / 2;

			for (int row = 0; row < a.length; row++) {
				v[0][row] = -a[0][1] / (a[0][0] - l[row][row]);
				v[1][row] = 1;
			}

			for (int row = 0; row < a.length; row++) {
				double mag = Math.sqrt((v[0][row] * v[0][row]) + 1);
				v[0][row] /= mag;
				v[1][row] /= mag;
			}

		} else if (a.length == 3) {

			// assumed double a = -1;
			double b = a[0][0] + a[1][1] + a[2][2];
			double c = -(a[0][0] * a[1][1]) - (a[0][0] * a[2][2])
					- (a[1][1] * a[2][2]) + (a[0][1] * a[1][0])
					+ (a[0][2] * a[2][0]) + (a[1][2] * a[2][1]);
			double d = (a[2][0] * a[0][1] * a[1][2])
					+ (a[1][0] * a[2][1] * a[0][2])
					- (a[0][0] * a[2][1] * a[1][2])
					- (a[1][1] * a[2][0] * a[0][2])
					- (a[2][2] * a[1][0] * a[0][1])
					+ (a[0][0] * a[1][1] * a[2][2]);

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

			// 19 variables...

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

				// inline forward elimination

				for (int diag = 0; diag < x.length - 1; diag++) {
					if (x[diag][diag] == 0)
						throw new ArithmeticException("divide by 0");

					for (int xrow = diag + 1; xrow < x.length; xrow++) {
						double quotient = x[xrow][diag] / x[diag][diag];

						for (int xcol = 0; xcol < x[xrow].length; xcol++)
							x[xrow][xcol] -= (quotient * x[diag][xcol]);
					}
				}

				// inline backward elimination

				for (int diag = x.length - 1; diag > 0; diag--) {
					if (x[diag][diag] == 0)
						throw new ArithmeticException("divide by 0");

					for (int xrow = diag - 1; xrow >= 0; xrow--) {
						double quotient = x[xrow][diag] / x[diag][diag];

						for (int xcol = 0; xcol < x[xrow].length; xcol++)
							x[xrow][xcol] -= (quotient * x[diag][xcol]);
					}
				}

				// solve for "z" where z is on row 5

				v[0][row] = -x[0][2] / x[0][0];
				v[1][row] = -x[1][2] / x[1][1];
				v[2][row] = 1;
			}

			// normalize by row

			for (int row = 0; row < a.length; row++) {

				// determinant on vector

				double det = Math.sqrt((v[0][row] * v[0][row])
						+ (v[1][row] * v[1][row]) + 1);

				v[0][row] /= det;
				v[1][row] /= det;
				v[2][row] /= det;
			}

			return new double[][][] { l, v };

		} else
			// because finding zeros of larger polynomials is difficult
			throw new UnsupportedOperationException(
					"Cannot handle matrices larger than 3 by 3");

		return new double[][][] { l, v };
	}

	/**
	 * Assumes row by columnn array representation of matrix. Given the
	 * resulting array,the first element is matrix u, the second element is
	 * matrix s, and the third element is matrix v. Input array a equals
	 * u*s*(v^-1) (through matrix multiplication and inverse of matrix). Given
	 * that ut is tranpose of u, ut*u (through matrix multiplication) equals I
	 * where I is identity matrix, and given that vt is tranpose of v, vt*v
	 * (through matrix multiplication) equals I where I is identity matrix.
	 * Resulting matrix s is diagonal such that elements on diagonal are
	 * non-zero and element not on diagonal equal zero. Returns array of m by n,
	 * n by n, and n by n array of double.
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
	 *            m by n array of double where m >= n
	 * @return array of m by n, n by n, and n by n array of double
	 */
	public static double[][][] singularValue(double[][] a) {
		if (null == a)
			throw new IllegalArgumentException("matrix a must not be null");

		int acol = -1;

		for (int row = 0; row < a.length; row++)
			if (null == a[row])
				throw new IllegalArgumentException(
						"matrix a must not contain null column");

			else if (row == 0) {
				acol = a[0].length;
				if (a.length < acol)
					throw new IllegalArgumentException(
							"matrix a must contain as many or more rows than columns");
			}

			else if (acol != a[row].length)
				throw new IllegalArgumentException(
						"matrix a row width must remain constant");

		// inline tranpose

		double[][] at = new double[acol][a.length];

		for (int row = 0; row < a.length; row++)
			for (int col = 0; col < acol; col++)
				// row becomes column; column becomes row
				at[col][row] = a[row][col];

		double[][] ata = new double[acol][acol];

		// inline matrix multiplication

		for (int row = 0; row < ata.length; row++) {
			for (int col = 0; col < ata[row].length; col++) {

				// inline dot product

				double sum = 0;

				for (int brow = 0; brow < a.length; brow++)
					sum += at[row][brow] * a[brow][col];

				ata[row][col] = sum;
			}
		}

		// inline eigen

		double[][] l = new double[ata.length][ata.length];
		double[][] v = new double[ata.length][ata.length];

		if (ata.length == 1) {

			// assumed double a = -1;
			double b = ata[0][0];

			l[0][0] = b;

			v[1][0] = 1;

		} else if (ata.length == 2) {

			// assumed double a = 1;
			double b = -ata[0][0] - ata[1][1];
			double c = ata[0][0] * ata[1][1] - ata[0][1] * ata[1][0];
			double sqrt = Math.sqrt((b * b) - (4 * c));

			l[0][0] = (-b + sqrt) / 2;
			l[1][1] = (-b - sqrt) / 2;

			for (int row = 0; row < ata.length; row++) {
				v[0][row] = -ata[0][1] / (ata[0][0] - l[row][row]);
				v[1][row] = 1;
			}

			for (int row = 0; row < ata.length; row++) {
				double mag = Math.sqrt((v[0][row] * v[0][row]) + 1);
				v[0][row] /= mag;
				v[1][row] /= mag;
			}

		} else if (ata.length == 3) {

			// assumed double a = -1;
			double b = ata[0][0] + ata[1][1] + ata[2][2];
			double c = -(ata[0][0] * ata[1][1]) - (ata[0][0] * ata[2][2])
					- (ata[1][1] * ata[2][2]) + (ata[0][1] * ata[1][0])
					+ (ata[0][2] * ata[2][0]) + (ata[1][2] * ata[2][1]);
			double d = (ata[2][0] * ata[0][1] * ata[1][2])
					+ (ata[1][0] * ata[2][1] * ata[0][2])
					- (ata[0][0] * ata[2][1] * ata[1][2])
					- (ata[1][1] * ata[2][0] * ata[0][2])
					- (ata[2][2] * ata[1][0] * ata[0][1])
					+ (ata[0][0] * ata[1][1] * ata[2][2]);

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

			// more than 19 variables

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

			for (int row = 0; row < ata.length; row++) {

				// copy array sig into array x

				double[][] x = new double[ata.length - 1][ata.length];
				for (int xrow = 0; xrow < x.length; xrow++)
					for (int xcol = 0; xcol < ata.length; xcol++) {
						x[xrow][xcol] = ata[xrow][xcol];
						if (xrow == xcol)
							// if on the diagonal, then subtract lambda
							x[xrow][xcol] -= l[row][row];
					}

				// inline forward elimination

				for (int diag = 0; diag < x.length - 1; diag++) {
					if (x[diag][diag] == 0)
						throw new ArithmeticException("divide by 0");

					for (int xrow = diag + 1; xrow < x.length; xrow++) {
						double quotient = x[xrow][diag] / x[diag][diag];

						for (int xcol = 0; xcol < x[xrow].length; xcol++)
							x[xrow][xcol] -= (quotient * x[diag][xcol]);
					}
				}

				// inline backward elimination

				for (int diag = x.length - 1; diag > 0; diag--) {
					if (x[diag][diag] == 0)
						throw new ArithmeticException("divide by 0");

					for (int xrow = diag - 1; xrow >= 0; xrow--) {
						double quotient = x[xrow][diag] / x[diag][diag];

						for (int xcol = 0; xcol < x[xrow].length; xcol++)
							x[xrow][xcol] -= (quotient * x[diag][xcol]);
					}
				}

				// solve for "z" where z is on row 5

				v[0][row] = -x[0][2] / x[0][0];
				v[1][row] = -x[1][2] / x[1][1];
				v[2][row] = 1;
			}

			// normalize by row

			for (int row = 0; row < ata.length; row++) {

				// determinant on vector

				double det = Math.sqrt((v[0][row] * v[0][row])
						+ (v[1][row] * v[1][row]) + 1);
				v[0][row] /= det;
				v[1][row] /= det;
				v[2][row] /= det;
			}

		} else
			// because finding zeros of larger polynomials is difficult
			throw new UnsupportedOperationException(
					"Cannot handle matrices larger than 3 by 3");

		// copy array v into array usv

		double[][] sigma = new double[ata.length][ata.length];

		for (int row = 0; row < ata.length; row++)
			for (int col = 0; col < ata.length; col++)
				sigma[row][col] = Math.sqrt(l[row][col]);

		double[][] sigmaInv = new double[ata.length][ata.length];

		for (int row = 0; row < ata.length; row++)
			for (int col = 0; col < ata.length; col++) {
				if (sigma[row][col] == 0)
					// prevent divides by 0
					sigmaInv[row][col] = 0;

				else
					// inverse is simple inverse since sigma is diagonal
					sigmaInv[row][col] = 1 / sigma[row][col];
			}

		double[][] av = new double[a.length][v[0].length];

		// inline matrix multiply

		for (int row = 0; row < av.length; row++) {
			for (int col = 0; col < av[row].length; col++) {

				// inline dot product

				double sum = 0;

				for (int brow = 0; brow < v.length; brow++)
					sum += a[row][brow] * v[brow][col];

				av[row][col] = sum;
			}
		}

		double[][] u = new double[av.length][sigmaInv[0].length];

		// inline matrix multiply

		for (int row = 0; row < u.length; row++) {
			for (int col = 0; col < u[row].length; col++) {

				// inline dot product

				double sum = 0;

				for (int brow = 0; brow < sigmaInv.length; brow++)
					sum += av[row][brow] * sigmaInv[brow][col];

				u[row][col] = sum;
			}
		}

		// 27 variables NOT counting any variable declared inside for loop

		return new double[][][] { u, sigma, v };
	}

}