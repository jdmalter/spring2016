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
	 * rows 0 through n of resulting matrix and u is rows n through n+n of
	 * resulting matrix, l*u (through matrix multiplication) equals a. Returns
	 * n+n by n array of double.
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
	 * @return n+n by n array of double
	 */
	public static double[][] lu(double[][] a) {
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

		double[][] lu = new double[a.length + a.length][a.length];

		// copy array l and array u into array lu

		for (int row = 0; row < lu.length; row++)
			for (int col = 0; col < a.length; col++)
				if (row < a.length)
					lu[row][col] = l[row][col];

				else
					lu[row][col] = u[row - a.length][col];

		return lu;
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
	 * rows 0 through m of resulting matrix and r is rows m through m+n of
	 * resulting matrix, q*r (through matrix multiplication) equals a. Given
	 * that qt is tranpose of q, qt*q equals I where I is identity matrix. r
	 * must be upper triangular. Returns m+n by n array of double.
	 * 
	 * Throws IllegalArgumentException if any of the following is true:
	 * {@code (null == a)}, {@code (null == a[row])} where a[row] is any row,
	 * {@code (a.length < acol)} where length is number rows and acol is number
	 * of columns on first row, or {@code (acol != a[row].length)} where acol is
	 * number of columns on first row and a[row].length is number of columns.
	 * 
	 * @param a
	 *            n by m array of double where m >= n
	 * @return m+n by n array of double
	 */
	public static double[][] qr(double[][] a) {
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

		double[][] qr = new double[a.length + acol][acol];

		for (int row = 0; row < qr.length; row++)
			for (int col = 0; col < acol; col++)
				if (row < a.length)
					qr[row][col] = q[row][col];

				else
					qr[row][col] = r[row - a.length][col];

		return qr;
	}

	/**
	 * Assumes row by columnn array representation of matrix. Given that l is
	 * rows 0 through n of the resulting matrix, and v is rows n through n+n of
	 * the resulting matrix, a*v (through matrix multiplication) equals v*l
	 * (through matrix multiplication). Returns n+n by n array of double.
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
	 * @return n+n by n array of double
	 */
	public static double[][] eigen(double[][] a) {
		if (null == a)
			throw new IllegalArgumentException("matrix a must not be null");

		for (int row = 0; row < a.length; row++)
			if (null == a[row])
				throw new IllegalArgumentException(
						"matrix a must not contain null column");

			else if (a.length != a[row].length)
				throw new IllegalArgumentException(
						"matrix a number of rows and columns must be equal");

		double[][] lv = new double[a.length + a.length][a.length];

		if (a.length == 1) {

			// assumed double a = -1;
			double b = a[0][0];

			lv[0][0] = b;

			lv[1][0] = 1;

		} else if (a.length == 2) {

			// assumed double a = 1;
			double b = -a[0][0] - a[1][1];
			double c = a[0][0] * a[1][1] - a[0][1] * a[1][0];
			double sqrt = Math.sqrt((b * b) - (4 * c));

			lv[0][0] = (-b + sqrt) / 2;
			lv[1][1] = (-b - sqrt) / 2;

			for (int row = 0; row < a.length; row++) {
				lv[2][row] = -a[0][1] / (a[0][0] - lv[row][row]);
				lv[3][row] = 1;
			}

			for (int row = 0; row < a.length; row++) {
				double mag = Math.sqrt((lv[2][row] * lv[2][row]) + 1);
				lv[2][row] /= mag;
				lv[3][row] /= mag;
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
			double l = -j;
			double m = Math.cos(k / 3);
			double n = Math.sqrt(3) * Math.sin(k / 3);
			double p = (b / 3);
			double q = 2 * j * Math.cos(k / 3) + (b / 3);
			double r = l * (m + n) + p;
			double s = l * (m - n) + p;

			// 19 variables...

			if (q > r && q > s) {

				// q is greatest root

				lv[0][0] = q;

				if (r > s) {

					// s is least root

					lv[1][1] = r;
					lv[2][2] = s;

				} else {

					// r is least root

					lv[1][1] = s;
					lv[2][2] = r;
				}

			} else if (r > s) {

				// r is greatest root

				lv[0][0] = r;

				if (q > s) {

					// s is least root

					lv[1][1] = q;
					lv[2][2] = s;

				} else {

					// q is least root

					lv[1][1] = s;
					lv[2][2] = q;
				}

			} else {

				// s is greatest root

				lv[0][0] = s;

				if (q > r) {

					// r is least root

					lv[1][1] = q;
					lv[2][2] = r;

				} else {

					// q is least root

					lv[1][1] = r;
					lv[2][2] = q;
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
							x[xrow][xcol] -= lv[row][row];
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

				lv[3][row] = -x[0][2] / x[0][0];
				lv[4][row] = -x[1][2] / x[1][1];
				lv[5][row] = 1;
			}

			// normalize by row

			for (int row = 0; row < a.length; row++) {

				// determinant on vector

				double det = Math.sqrt((lv[3][row] * lv[3][row])
						+ (lv[4][row] * lv[4][row]) + 1);

				lv[3][row] /= det;
				lv[4][row] /= det;
				lv[5][row] /= det;
			}

			return lv;

		} else
			// because finding zeros of larger polynomials is difficult
			throw new UnsupportedOperationException(
					"Cannot handle matrices larger than 3 by 3");

		return lv;
	}

	/**
	 * Assumes row by columnn array representation of matrix. Given the
	 * resulting array, rows 0 to m belong to matrix u, rows m to m+n belong to
	 * matrix s, and rows m+n to m+n+n belong to matrix v. Input array a equals
	 * u*s*(v^-1) (through matrix multiplication and inverse of matrix). Given
	 * that ut is tranpose of u, ut*u (through matrix multiplication) equals I
	 * where I is identity matrix, and given that vt is tranpose of v, vt*v
	 * (through matrix multiplication) equals I where I is identity matrix.
	 * Resulting matrix s is diagonal such that elements on diagonal are
	 * non-zero and element not on diagonal equal zero. Returns m+n+n by n array
	 * of double.
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
	 * @return m+n+n by n array of double
	 */
	public static double[][] singularValue(double[][] a) {
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

		double[][] sig = new double[acol][acol];

		// inline matrix multiplication

		for (int row = 0; row < sig.length; row++) {
			for (int col = 0; col < sig[row].length; col++) {

				// inline dot product

				double sum = 0;

				for (int brow = 0; brow < a.length; brow++)
					sum += at[row][brow] * a[brow][col];

				sig[row][col] = sum;
			}
		}

		// inline eigen

		double[][] lv = new double[sig.length + sig.length][sig.length];

		if (sig.length == 1) {

			// assumed double a = -1;
			double b = sig[0][0];

			lv[0][0] = b;

			lv[1][0] = 1;

		} else if (sig.length == 2) {

			// assumed double a = 1;
			double b = -sig[0][0] - sig[1][1];
			double c = sig[0][0] * sig[1][1] - sig[0][1] * sig[1][0];
			double sqrt = Math.sqrt((b * b) - (4 * c));

			lv[0][0] = (-b + sqrt) / 2;
			lv[1][1] = (-b - sqrt) / 2;

			for (int row = 0; row < sig.length; row++) {
				lv[2][row] = -sig[0][1] / (sig[0][0] - lv[row][row]);
				lv[3][row] = 1;
			}

			for (int row = 0; row < sig.length; row++) {
				double mag = Math.sqrt((lv[2][row] * lv[2][row]) + 1);
				lv[2][row] /= mag;
				lv[3][row] /= mag;
			}

		} else if (sig.length == 3) {

			// assumed double a = -1;
			double b = sig[0][0] + sig[1][1] + sig[2][2];
			double c = -(sig[0][0] * sig[1][1]) - (sig[0][0] * sig[2][2])
					- (sig[1][1] * sig[2][2]) + (sig[0][1] * sig[1][0])
					+ (sig[0][2] * sig[2][0]) + (sig[1][2] * sig[2][1]);
			double d = (sig[2][0] * sig[0][1] * sig[1][2])
					+ (sig[1][0] * sig[2][1] * sig[0][2])
					- (sig[0][0] * sig[2][1] * sig[1][2])
					- (sig[1][1] * sig[2][0] * sig[0][2])
					- (sig[2][2] * sig[1][0] * sig[0][1])
					+ (sig[0][0] * sig[1][1] * sig[2][2]);

			double f = (-(3 * c) - (b * b)) / 3;
			double g = (-(2 * b * b * b) - (9 * b * c) - (27 * d)) / 27;
			double h = (g * g / 4) + (f * f * f / 27);
			double i = Math.sqrt((g * g / 4) - h);
			double j = Math.pow(i, 1d / 3);
			double k = Math.acos(-g / (2 * i));
			double l = -j;
			double m = Math.cos(k / 3);
			double n = Math.sqrt(3) * Math.sin(k / 3);
			double p = (b / 3);
			double q = 2 * j * Math.cos(k / 3) + (b / 3);
			double r = l * (m + n) + p;
			double s = l * (m - n) + p;

			// more than 19 variables

			if (q > r && q > s) {

				// q is greatest root

				lv[0][0] = q;

				if (r > s) {

					// s is least root

					lv[1][1] = r;
					lv[2][2] = s;

				} else {

					// r is least root

					lv[1][1] = s;
					lv[2][2] = r;
				}

			} else if (r > s) {

				// r is greatest root

				lv[0][0] = r;

				if (q > s) {

					// s is least root

					lv[1][1] = q;
					lv[2][2] = s;

				} else {

					// q is least root

					lv[1][1] = s;
					lv[2][2] = q;
				}

			} else {

				// s is greatest root

				lv[0][0] = s;

				if (q > r) {

					// r is least root

					lv[1][1] = q;
					lv[2][2] = r;

				} else {

					// q is least root

					lv[1][1] = r;
					lv[2][2] = q;
				}
			}

			for (int row = 0; row < sig.length; row++) {

				// copy array sig into array x

				double[][] x = new double[sig.length - 1][sig.length];
				for (int xrow = 0; xrow < x.length; xrow++)
					for (int xcol = 0; xcol < sig.length; xcol++) {
						x[xrow][xcol] = sig[xrow][xcol];
						if (xrow == xcol)
							// if on the diagonal, then subtract lambda
							x[xrow][xcol] -= lv[row][row];
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

				lv[3][row] = -x[0][2] / x[0][0];
				lv[4][row] = -x[1][2] / x[1][1];
				lv[5][row] = 1;
			}

			// normalize by row

			for (int row = 0; row < sig.length; row++) {

				// determinant on vector

				double det = Math.sqrt((lv[3][row] * lv[3][row])
						+ (lv[4][row] * lv[4][row]) + 1);
				lv[3][row] /= det;
				lv[4][row] /= det;
				lv[5][row] /= det;
			}

		} else
			// because finding zeros of larger polynomials is difficult
			throw new UnsupportedOperationException(
					"Cannot handle matrices larger than 3 by 3");

		double[][] usv = new double[a.length + acol + acol][acol];

		// copy array v into array usv

		for (int row = a.length + acol; row < a.length + acol + acol; row++)
			for (int col = 0; col < acol; col++)
				usv[row][col] = lv[row - a.length][col];

		// copy array l into array usv

		for (int row = a.length; row < a.length + acol; row++)
			for (int col = 0; col < acol; col++)
				usv[row][col] = Math.sqrt(lv[row - a.length][col]);

		double[][] v = new double[acol][acol];

		// find array v inside array lv

		for (int row = 0; row < acol; row++)
			for (int col = 0; col < acol; col++)
				v[row][col] = lv[row + acol][col];

		double[][] sigInv = new double[acol][acol];

		for (int row = a.length; row < a.length + acol; row++)
			for (int col = 0; col < acol; col++) {
				if (usv[row][col] == 0)
					// prevent divides by 0
					sigInv[row - a.length][col] = 0;

				else
					// inverse is simple inverse since sigma is diagonal
					sigInv[row - a.length][col] = 1 / usv[row][col];
			}

		double[][] mul = new double[a.length][v[0].length];

		// inline matrix multiply

		for (int row = 0; row < mul.length; row++) {
			for (int col = 0; col < mul[row].length; col++) {

				// inline dot product

				double sum = 0;

				for (int brow = 0; brow < v.length; brow++)
					sum += a[row][brow] * v[brow][col];

				mul[row][col] = sum;
			}
		}

		double[][] u = new double[mul.length][sigInv[0].length];

		// inline matrix multiply

		for (int row = 0; row < u.length; row++) {
			for (int col = 0; col < u[row].length; col++) {

				// inline dot product

				double sum = 0;

				for (int brow = 0; brow < sigInv.length; brow++)
					sum += mul[row][brow] * sigInv[brow][col];

				u[row][col] = sum;
			}
		}

		// copy array u into array usv

		for (int row = 0; row < a.length; row++)
			for (int col = 0; col < acol; col++)
				usv[row][col] = u[row][col];

		// 27 variables NOT counting any variable declared inside for loop

		return usv;
	}

}