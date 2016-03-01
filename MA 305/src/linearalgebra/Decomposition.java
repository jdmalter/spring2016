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

	// Related to LU Decomposition

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

		for (int row = 0; row < l.length; row++)
			for (int col = row; col < l.length; col++)
				if (row == col)
					l[row][col] = 1;

		for (int row = 0; row < u.length; row++)
			for (int col = 0; col < u[row].length; col++)
				u[row][col] = a[row][col];

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

		for (int row = 0; row < lu.length; row++)
			for (int col = 0; col < a.length; col++)
				if (row < a.length)
					lu[row][col] = l[row][col];

				else
					lu[row][col] = u[row - a.length][col];

		return lu;
	}

	// Related to Cholesky Decomposition

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

	// Related to QR Decomposition

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

	// Related to Eigen Decomposition

	/**
	 * Assumes row by columnn array representation of matrix. Given that l is
	 * rows 0 through n of the resulting matrix, and v is rows n through n+n of
	 * the resulting matrix, v being some basis vector, a*v (through matrix
	 * multiplication) equals l*v (through matrix multiplication), or the
	 * determinant of (l*I - A) equals 0 where I is identity matrix. Returns n+n
	 * by n array of double.
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

			double x1 = 2 * j * Math.cos(k / 3) + (b / 3);
			double x2 = l * (m + n) + p;
			double x3 = l * (m - n) + p;

			if (x1 > x2 && x1 > x3) {
				lv[0][0] = x1;

				if (x2 > x3) {
					lv[1][1] = x2;
					lv[2][2] = x3;

				} else {
					lv[1][1] = x3;
					lv[2][2] = x2;
				}

			} else if (x2 > x3) {
				lv[0][0] = x2;

				if (x1 > x3) {
					lv[1][1] = x1;
					lv[2][2] = x3;

				} else {
					lv[1][1] = x3;
					lv[2][2] = x1;
				}

			} else {
				lv[0][0] = x3;

				if (x1 > x2) {
					lv[1][1] = x1;
					lv[2][2] = x2;

				} else {
					lv[1][1] = x2;
					lv[2][2] = x1;
				}
			}

			for (int row = 0; row < a.length; row++) {
				double[][] x = new double[a.length - 1][a.length];
				for (int xrow = 0; xrow < x.length; xrow++)
					for (int xcol = 0; xcol < a.length; xcol++) {
						x[xrow][xcol] = a[xrow][xcol];
						if (xrow == xcol)
							x[xrow][xcol] -= lv[row][row];
					}

				for (int diag = 0; diag < x.length - 1; diag++) {
					if (x[diag][diag] == 0)
						throw new ArithmeticException("divide by 0");

					for (int xrow = diag + 1; xrow < x.length; xrow++) {
						double quotient = x[xrow][diag] / x[diag][diag];

						for (int xcol = 0; xcol < x[xrow].length; xcol++)
							x[xrow][xcol] -= (quotient * x[diag][xcol]);
					}
				}

				for (int diag = x.length - 1; diag > 0; diag--) {
					if (x[diag][diag] == 0)
						throw new ArithmeticException("divide by 0");

					for (int xrow = diag - 1; xrow >= 0; xrow--) {
						double quotient = x[xrow][diag] / x[diag][diag];

						for (int xcol = 0; xcol < x[xrow].length; xcol++)
							x[xrow][xcol] -= (quotient * x[diag][xcol]);
					}
				}

				lv[3][row] = -x[0][2] / x[0][0];
				lv[4][row] = -x[1][2] / x[1][1];
				lv[5][row] = 1;
			}

			for (int row = 0; row < a.length; row++) {
				double mag = Math.sqrt((lv[3][row] * lv[3][row])
						+ (lv[4][row] * lv[4][row]) + 1);
				lv[3][row] /= mag;
				lv[4][row] /= mag;
				lv[5][row] /= mag;
			}

			return lv;

		} else
			throw new UnsupportedOperationException(
					"Cannot handle matrices larger than 3 by 3");

		return lv;
	}

	// Related to Singular Value Decomposition

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

		double[][] usv = new double[a.length + acol + acol][acol];

		double[][] lv = Decomposition.eigen(Matrices.multiply(
				Matrices.transpose(a), a));

		for (int row = a.length + acol; row < a.length + acol + acol; row++)
			for (int col = 0; col < acol; col++)
				usv[row][col] = lv[row - a.length][col];

		for (int row = a.length; row < a.length + acol; row++)
			for (int col = 0; col < acol; col++)
				usv[row][col] = Math.sqrt(lv[row - a.length][col]);

		double[][] v = new double[acol][acol];

		for (int row = 0; row < acol; row++)
			for (int col = 0; col < acol; col++)
				v[row][col] = lv[row + acol][col];

		double[][] sInv = new double[acol][acol];

		for (int row = a.length; row < a.length + acol; row++)
			for (int col = 0; col < acol; col++) {
				if (usv[row][col] == 0)
					sInv[row - a.length][col] = 0;
				else
					sInv[row - a.length][col] = 1 / usv[row][col];
			}

		double[][] u = Matrices.multiply(Matrices.multiply(a, v), sInv);

		for (int row = 0; row < a.length; row++)
			for (int col = 0; col < acol; col++)
				usv[row][col] = u[row][col];

		return usv;
	}

}