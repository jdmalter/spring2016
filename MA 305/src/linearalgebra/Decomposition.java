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

				for (int p = 0; p < row; ++p)
					lsum += l[lrow][p] * l[row][p];

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
				for (int arow = 0; arow < a.length; arow++)
					r[rrow][col] += q[arow][rrow] * a[arow][col];

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
	 * Assumes row by columnn array representation of matrix. Given that v is
	 * resulting matrix, v is some basis vector, and l is lambda, a*v (through
	 * matrix multiplication) equals l*v (through matrix multiplication), or the
	 * determinant of (l*I - A) equals 0 where I is identity matrix. Returns n
	 * by n array of double.
	 * 
	 * Throws IllegalArgumentException if any of the following is true:
	 * {@code (null == a)}, {@code (null == a[row])} where a[row] is any row, or
	 * {@code (acol != a[row].length)} where acol is number of columns on first
	 * row and a[row].length is number of columns.
	 * 
	 * @param a
	 *            n by n array of double
	 * @return n by n array of double
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

		double[][] v = new double[a.length][a.length];

		// TODO learn more about eigen decomposition

		return v;
	}

}