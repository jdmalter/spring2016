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
	 * Assumes row by columnn array representation of matrix. Given that x is
	 * the first element of the resulting array and e is the second element of
	 * the resulting array, a*x - b = e where the magnitude of e is minimized.
	 * Returns array of array of n by 1 and m by 1 array of double,
	 * 
	 * Given that at is tranpose of given matrix a, this implementation finds
	 * the pseudo-inverse of matrix a by solving ((at*a)^-1)*at.
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
	 * There is a special condition where divide by 0 occurs when the
	 * determinant of the resulting matrix from at*a equals 0.
	 * 
	 * @param a
	 *            m by n array of double where m > n
	 * @param b
	 *            m by 1 array of double
	 * @return array of n by 1 and m by 1 array of double
	 */
	public static double[][][] overconstrainedPI(double[][] a, double[][] b) {
		if (null == a)
			throw new IllegalArgumentException("matrix a must not be null");

		else if (null == b)
			throw new IllegalArgumentException("matrix b must not be null");

		else if (a.length != b.length)
			throw new IllegalArgumentException(
					"matrix a and b number of rows must be equal");

		int acol = -1;

		for (int row = 0; row < a.length; row++)
			if (null == a[row])
				throw new IllegalArgumentException(
						"matrix a must not contain null column");

			else if (null == b[row])
				throw new IllegalArgumentException(
						"matrix b must not contain null column");

			else if (row == 0) {
				acol = a[0].length;
				if (a.length <= acol)
					throw new IllegalArgumentException(
							"matrix a must contain more rows than columns");

			} else if (acol != a[row].length)
				throw new IllegalArgumentException(
						"matrix a row width must remain constant");

			else if (1 != b[row].length)
				throw new IllegalArgumentException(
						"matrix b must only contain one column");

		// inline matrix multiply and tranpose

		double[][] ata = new double[acol][acol];

		// find dot products of every row and column combination

		for (int row = 0; row < acol; row++) {
			for (int col = 0; col < acol; col++) {

				// inline dot product

				double sum = 0;

				for (int brow = 0; brow < a.length; brow++)
					sum += a[brow][row] * a[brow][col];

				ata[row][col] = sum;
			}
		}

		// inline invert

		double[][] adjugate = new double[ata.length][ata.length];
		double determinant = 0;

		if (ata.length == 0) {

			// inline matrix determinant on point

			// matrix is singular
			throw new ArithmeticException("divide by 0");

		} else if (ata.length == 1) {

			// inline matrix determinant on line

			determinant = ata[0][0];
			if (determinant == 0)
				// matrix is singular
				throw new ArithmeticException("divide by 0");

			adjugate[0][0] = 1;

		} else if (ata.length == 2) {

			// inline matrix determinant on vectors

			determinant = ata[0][0] * ata[1][1] - ata[0][1] * ata[1][0];
			if (determinant == 0)
				// matrix is singular
				throw new ArithmeticException("divide by 0");

			adjugate[0][0] = ata[1][1];
			adjugate[0][1] = -ata[0][1];
			adjugate[1][0] = -ata[1][0];
			adjugate[1][1] = ata[0][0];

		} else if (ata.length == 3) {

			// inline matrix determinant on triple product of vectors

			determinant = (ata[0][0] * (ata[1][1] * ata[2][2] - ata[1][2]
					* ata[2][1]))
					- (ata[0][1] * (ata[1][0] * ata[2][2] - ata[1][2]
							* ata[2][0]))
					+ (ata[0][2] * (ata[1][0] * ata[2][1] - ata[1][1]
							* ata[2][0]));
			if (determinant == 0)
				// matrix is singular
				throw new ArithmeticException("divide by 0");

			// multiple inline matrix determinants on vectors

			adjugate[0][0] = (ata[1][1] * ata[2][2] - ata[1][2] * ata[2][1]);
			adjugate[0][1] = -(ata[0][1] * ata[2][2] - ata[0][2] * ata[2][1]);
			adjugate[0][2] = (ata[0][1] * ata[1][2] - ata[0][2] * ata[1][1]);
			adjugate[1][0] = -(ata[1][0] * ata[2][2] - ata[1][2] * ata[2][0]);
			adjugate[1][1] = (ata[0][0] * ata[2][2] - ata[0][2] * ata[2][0]);
			adjugate[1][2] = -(ata[0][0] * ata[1][2] - ata[0][2] * ata[1][0]);
			adjugate[2][0] = (ata[1][0] * ata[2][1] - ata[1][1] * ata[2][0]);
			adjugate[2][1] = -(ata[0][0] * ata[2][1] - ata[0][1] * ata[2][0]);
			adjugate[2][2] = (ata[0][0] * ata[1][1] - ata[0][1] * ata[1][0]);

		} else {

			determinant = Matrices.determinant(ata);
			if (determinant == 0)
				// matrix is singular
				throw new ArithmeticException("divide by 0");

			for (int row = 0; row < ata.length; row++)
				for (int col = 0; col < ata.length; col++) {

					// inline matrix minor

					double[][] minor = new double[ata.length - 1][ata.length - 1];

					for (int mrow = 0, srow = 0; mrow < minor.length; mrow++, srow++) {
						if (row == mrow)
							// ignore row being skipped and access new row
							srow++;

						for (int mcol = 0, scol = 0; mcol < minor.length; mcol++, scol++) {
							if (col == mcol)
								// ignore column being skipped and access new
								// col
								scol++;

							minor[mrow][mcol] = ata[srow][scol];
						}
					}

					double det = Matrices.determinant(minor);

					if ((row + col) % 2 == 0)
						adjugate[col][row] = det;

					else
						adjugate[col][row] = -det;
				}
		}

		// divide adjugate elements by determinant

		for (int row = 0; row < ata.length; row++)
			for (int col = 0; col < ata.length; col++)
				adjugate[row][col] /= determinant;

		// inline matrix multiply

		double[][] adjugateAt = new double[acol][a.length];

		// find dot products of every row and column combination

		for (int row = 0; row < adjugateAt.length; row++) {
			for (int col = 0; col < adjugateAt[row].length; col++) {

				// inline dot product

				double sum = 0;

				for (int brow = 0; brow < acol; brow++)
					sum += adjugate[row][brow] * a[col][brow];

				adjugateAt[row][col] = sum;
			}
		}

		// inline matrix multiply on vector

		double[][] x = new double[acol][1];

		// find dot products of every row and column combination

		for (int row = 0; row < x.length; row++) {

			// inline dot product

			double sum = 0;

			for (int brow = 0; brow < b.length; brow++)
				sum += adjugateAt[row][brow] * b[brow][0];

			x[row][0] = sum;
		}

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
	 * Assumes row by columnn array representation of matrix. Given that x is
	 * resulting matrix, a*x - b = 0 where the magnitude of x is minimized.
	 * Returns n by 1 array of double.
	 * 
	 * Given that at is tranpose of given matrix a, this implementation finds
	 * the pseudo-inverse of matrix a by solving at*((a*at)^-1).
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
	 * There is a special condition where divide by 0 occurs when the
	 * determinant of the resulting matrix from a*at equal 0.
	 * 
	 * @param a
	 *            m by n array of double where m < n
	 * @param b
	 *            m by 1 array of double
	 * @return n by 1 array of double
	 */
	public static double[][] underconstrainedPI(double[][] a, double[][] b) {
		if (null == a)
			throw new IllegalArgumentException("matrix a must not be null");

		else if (null == b)
			throw new IllegalArgumentException("matrix b must not be null");

		else if (a.length != b.length)
			throw new IllegalArgumentException(
					"matrix a and b number of rows must be equal");

		int acol = -1;

		for (int row = 0; row < a.length; row++)
			if (null == a[row])
				throw new IllegalArgumentException(
						"matrix a must not contain null column");

			else if (null == b[row])
				throw new IllegalArgumentException(
						"matrix b must not contain null column");

			else if (row == 0) {
				acol = a[0].length;
				if (a.length >= acol)
					throw new IllegalArgumentException(
							"matrix a must contain less rows than columns");

			} else if (acol != a[row].length)
				throw new IllegalArgumentException(
						"matrix a row width must remain constant");

			else if (1 != b[row].length)
				throw new IllegalArgumentException(
						"matrix b must only contain one column");

		// inline matrix multiply and tranpose

		double[][] aat = new double[a.length][a.length];

		// find dot products of every row and column combination

		for (int row = 0; row < aat.length; row++) {
			for (int col = 0; col < aat[row].length; col++) {

				// inline dot product

				double sum = 0;

				for (int brow = 0; brow < a[0].length; brow++)
					sum += a[row][brow] * a[col][brow];

				aat[row][col] = sum;
			}
		}

		// inline invert

		double[][] adjugate = new double[aat.length][aat.length];
		double determinant = 0;

		if (aat.length == 0) {

			// inline matrix determinant on point

			// matrix is singular
			throw new ArithmeticException("divide by 0");

		} else if (aat.length == 1) {

			// inline matrix determinant on line

			determinant = aat[0][0];
			if (determinant == 0)
				// matrix is singular
				throw new ArithmeticException("divide by 0");

			adjugate[0][0] = 1;

		} else if (aat.length == 2) {

			// inline matrix determinant on vectors

			determinant = aat[0][0] * aat[1][1] - aat[0][1] * aat[1][0];
			if (determinant == 0)
				// matrix is singular
				throw new ArithmeticException("divide by 0");

			adjugate[0][0] = aat[1][1];
			adjugate[0][1] = -aat[0][1];
			adjugate[1][0] = -aat[1][0];
			adjugate[1][1] = aat[0][0];

		} else if (aat.length == 3) {

			// inline matrix determinant on triple product of vectors

			determinant = (aat[0][0] * (aat[1][1] * aat[2][2] - aat[1][2]
					* aat[2][1]))
					- (aat[0][1] * (aat[1][0] * aat[2][2] - aat[1][2]
							* aat[2][0]))
					+ (aat[0][2] * (aat[1][0] * aat[2][1] - aat[1][1]
							* aat[2][0]));
			if (determinant == 0)
				// matrix is singular
				throw new ArithmeticException("divide by 0");

			// multiple inline matrix determinants on vectors

			adjugate[0][0] = (aat[1][1] * aat[2][2] - aat[1][2] * aat[2][1]);
			adjugate[0][1] = -(aat[0][1] * aat[2][2] - aat[0][2] * aat[2][1]);
			adjugate[0][2] = (aat[0][1] * aat[1][2] - aat[0][2] * aat[1][1]);
			adjugate[1][0] = -(aat[1][0] * aat[2][2] - aat[1][2] * aat[2][0]);
			adjugate[1][1] = (aat[0][0] * aat[2][2] - aat[0][2] * aat[2][0]);
			adjugate[1][2] = -(aat[0][0] * aat[1][2] - aat[0][2] * aat[1][0]);
			adjugate[2][0] = (aat[1][0] * aat[2][1] - aat[1][1] * aat[2][0]);
			adjugate[2][1] = -(aat[0][0] * aat[2][1] - aat[0][1] * aat[2][0]);
			adjugate[2][2] = (aat[0][0] * aat[1][1] - aat[0][1] * aat[1][0]);

		} else {

			determinant = Matrices.determinant(aat);
			if (determinant == 0)
				// matrix is singular
				throw new ArithmeticException("divide by 0");

			for (int row = 0; row < aat.length; row++)
				for (int col = 0; col < aat.length; col++) {

					// inline matrix minor

					double[][] minor = new double[aat.length - 1][aat.length - 1];

					for (int mrow = 0, srow = 0; mrow < minor.length; mrow++, srow++) {
						if (row == mrow)
							// ignore row being skipped and access new row
							srow++;

						for (int mcol = 0, scol = 0; mcol < minor.length; mcol++, scol++) {
							if (col == mcol)
								// ignore column being skipped and access new
								// col
								scol++;

							minor[mrow][mcol] = aat[srow][scol];
						}
					}

					double det = Matrices.determinant(minor);

					if ((row + col) % 2 == 0)
						adjugate[col][row] = det;

					else
						adjugate[col][row] = -det;
				}
		}

		// divide adjugate elements by determinant

		for (int row = 0; row < aat.length; row++)
			for (int col = 0; col < aat.length; col++)
				adjugate[row][col] /= determinant;

		// inline matrix multiply

		double[][] atAdjugate = new double[acol][aat.length];

		// find dot products of every row and column combination

		for (int row = 0; row < atAdjugate.length; row++) {
			for (int col = 0; col < atAdjugate[row].length; col++) {

				// inline dot product

				double sum = 0;

				for (int brow = 0; brow < adjugate.length; brow++)
					sum += a[brow][row] * adjugate[brow][col];

				atAdjugate[row][col] = sum;
			}
		}

		// inline matrix multiply on vector

		double[][] x = new double[atAdjugate.length][1];

		// find dot products of every row and column combination

		for (int row = 0; row < x.length; row++) {

			// inline dot product

			double sum = 0;

			for (int brow = 0; brow < b.length; brow++)
				sum += atAdjugate[row][brow] * b[brow][0];

			x[row][0] = sum;
		}

		return x;
	}

	/**
	 * Assumes row by columnn array representation of matrix. Given that x is
	 * the first element of the resulting array and e is the second element of
	 * the resulting array, a*x - b = e where the magnitude of e is minimized.
	 * Returns array of array of n by 1 and m by 1 array of double,
	 * 
	 * This implementation finds x in two steps. It is given that at is tranpose
	 * of given matrix a, l is the cholesky decomposition of a, and lt is
	 * tranpose of l. First, y is found by solving l*y = at*b. Second, x is
	 * found by solving lt*x = y.
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
	 * There is a special condition where divide by 0 occurs when the
	 * determinant of the resulting matrix from at*a equals 0.
	 * 
	 * @param a
	 *            m by n array of double where m > n
	 * @param b
	 *            m by 1 array of double
	 * @return array of n by 1 and m by 1 array of double
	 */
	public static double[][][] overconstrainedCD(double[][] a, double[][] b) {
		if (null == a)
			throw new IllegalArgumentException("matrix a must not be null");

		else if (null == b)
			throw new IllegalArgumentException("matrix b must not be null");

		else if (a.length != b.length)
			throw new IllegalArgumentException(
					"matrix a and b number of rows must be equal");

		int acol = -1;

		for (int row = 0; row < a.length; row++)
			if (null == a[row])
				throw new IllegalArgumentException(
						"matrix a must not contain null column");

			else if (null == b[row])
				throw new IllegalArgumentException(
						"matrix b must not contain null column");

			else if (row == 0) {
				acol = a[0].length;
				if (a.length <= acol)
					throw new IllegalArgumentException(
							"matrix a must contain more rows than columns");

			} else if (acol != a[row].length)
				throw new IllegalArgumentException(
						"matrix a row width must remain constant");

			else if (1 != b[row].length)
				throw new IllegalArgumentException(
						"matrix b must only contain one column");

		// inline matrix multiply and tranpose

		double[][] ata = new double[acol][acol];

		// find dot products of every row and column combination

		for (int row = 0; row < acol; row++) {
			for (int col = 0; col < acol; col++) {

				// inline dot product

				double sum = 0;

				for (int brow = 0; brow < a.length; brow++)
					sum += a[brow][row] * a[brow][col];

				ata[row][col] = sum;
			}
		}

		// inline cholesy

		double[][] l = new double[acol][acol];

		for (int row = 0; row < acol; row++) {
			double sum = 0;

			// inline determinant on vector

			for (int col = 0; col < row; col++)
				sum += l[row][col] * l[row][col];
			l[row][row] = Math.sqrt(ata[row][row] - sum);

			for (int lrow = row + 1; lrow < acol; lrow++) {
				double lsum = 0;

				for (int lcol = 0; lcol < row; ++lcol)
					lsum += l[lrow][lcol] * l[row][lcol];

				if (l[row][row] == 0)
					throw new ArithmeticException("divide by 0");

				l[lrow][row] = (ata[lrow][row] - lsum) / l[row][row];
			}
		}

		// inline matrix multiply and tranpose

		double[][] atb = new double[acol][1];

		// find dot products of every row and column combination

		for (int row = 0; row < acol; row++) {

			// inline dot product

			double sum = 0;

			for (int brow = 0; brow < b.length; brow++)
				sum += a[brow][row] * b[brow][0];

			atb[row][0] = sum;
		}

		// inline gauss jordan

		// inline augment

		double[][] latb = new double[acol][acol + 1];

		for (int row = 0; row < acol; row++)
			for (int col = 0; col < acol; col++)
				latb[row][col] = l[row][col];

		for (int row = 0; row < acol; row++)
			latb[row][acol] = atb[row][0];

		// slight modification to the order of declarations
		double[][] y = new double[acol][1];

		// inline rref

		if (latb.length == 0) {
			// if l has 0 rows, then x has 0 rows
			y = new double[][] {};

		} else if (latb.length == 1) {
			if (latb[0][0] == 0)
				y = new double[][] { new double[] { 0 } };
			else
				y = new double[][] { new double[] { 1 } };

		} else {
			// copy array ab into array c

			double[][] c = new double[latb[0].length - 1][latb[0].length];
			for (int row = 0; row < c.length; row++)
				for (int col = 0; col < latb[0].length; col++)
					c[row][col] = latb[row][col];

			// inline forward elimination

			for (int diag = 0; diag < c.length - 1; diag++) {
				if (c[diag][diag] == 0)
					throw new ArithmeticException("divide by 0");

				for (int row = diag + 1; row < c.length; row++) {
					double quotient = c[row][diag] / c[diag][diag];

					for (int col = 0; col < c[row].length; col++)
						c[row][col] -= (quotient * c[diag][col]);
				}
			}

			// inline backward elimination

			for (int diag = c.length - 1; diag > 0; diag--) {
				if (c[diag][diag] == 0)
					throw new ArithmeticException("divide by 0");

				for (int row = diag - 1; row >= 0; row--) {
					double quotient = c[row][diag] / c[diag][diag];

					for (int col = 0; col < c[row].length; col++)
						c[row][col] -= (quotient * c[diag][col]);
				}
			}

			// normalize c

			double[][] rref = new double[acol][latb[0].length];

			for (int row = 0; row < c.length; row++) {
				rref[row][latb[0].length - 1] = c[row][latb[0].length - 1]
						/ c[row][row];
				rref[row][row] = 1;
			}

			// copy last column of results into array y

			for (int row = 0; row < y.length; row++)
				y[row][0] = rref[row][rref[row].length - 1] / rref[row][row];
		}

		// inline gauss jordan

		// inline augment

		double[][] lty = new double[acol][acol + 1];

		for (int row = 0; row < acol; row++)
			for (int col = 0; col < acol; col++)
				lty[row][col] = l[col][row];

		for (int row = 0; row < y.length; row++)
			lty[row][acol] = y[row][0];

		// slight modification to the order of declarations
		double[][] x = new double[acol][1];

		// inline rref

		if (lty.length == 0) {
			x = new double[][] {};

		} else if (lty.length == 1) {
			if (lty[0][0] == 0)
				x = new double[][] { new double[] { 0 } };
			else
				x = new double[][] { new double[] { 1 } };

		} else {
			// copy array ab into array c

			double[][] c = new double[lty[0].length - 1][lty[0].length];
			for (int row = 0; row < c.length; row++)
				for (int col = 0; col < lty[0].length; col++)
					c[row][col] = lty[row][col];

			// inline forward elimination

			for (int diag = 0; diag < c.length - 1; diag++) {
				if (c[diag][diag] == 0)
					throw new ArithmeticException("divide by 0");

				for (int row = diag + 1; row < c.length; row++) {
					double quotient = c[row][diag] / c[diag][diag];

					for (int col = 0; col < c[row].length; col++)
						c[row][col] -= (quotient * c[diag][col]);
				}
			}

			// inline backward elimination

			for (int diag = c.length - 1; diag > 0; diag--) {
				if (c[diag][diag] == 0)
					throw new ArithmeticException("divide by 0");

				for (int row = diag - 1; row >= 0; row--) {
					double quotient = c[row][diag] / c[diag][diag];

					for (int col = 0; col < c[row].length; col++)
						c[row][col] -= (quotient * c[diag][col]);
				}
			}

			// normalize c

			double[][] rref = new double[a.length][lty[0].length];

			for (int row = 0; row < c.length; row++) {
				rref[row][lty[0].length - 1] = c[row][lty[0].length - 1]
						/ c[row][row];
				rref[row][row] = 1;
			}

			// copy last column of results into array x

			for (int row = 0; row < x.length; row++)
				x[row][0] = rref[row][rref[row].length - 1] / rref[row][row];
		}

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
	 * Assumes row by columnn array representation of matrix. Given that x is
	 * the first element of the resulting array and e is the second element of
	 * the resulting array, a*x - b = e where the magnitude of e is minimized.
	 * Returns array of array of n by 1 and m by 1 array of double,
	 * 
	 * This implementation finds x in two steps. First, array qr is the
	 * resulting array from qr decomposition on given matrix a where q is the
	 * first element of qr and r is the second element of qr. Second, given that
	 * qt is tranpose of q, x is found by solving r*x = qt*b.
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
	 * There is a special condition where divide by 0 occurs when the
	 * determinant of the resulting matrix from at*a equals 0.
	 * 
	 * @param a
	 *            m by n array of double where m > n
	 * @param b
	 *            m by 1 array of double
	 * @return array of n by 1 and m by 1 array of double
	 */
	public static double[][][] overconstrainedQR(double[][] a, double[][] b) {
		if (null == a)
			throw new IllegalArgumentException("matrix a must not be null");

		else if (null == b)
			throw new IllegalArgumentException("matrix b must not be null");

		else if (a.length != b.length)
			throw new IllegalArgumentException(
					"matrix a and b number of rows must be equal");

		int acol = -1;

		for (int row = 0; row < a.length; row++)
			if (null == a[row])
				throw new IllegalArgumentException(
						"matrix a must not contain null column");

			else if (null == b[row])
				throw new IllegalArgumentException(
						"matrix b must not contain null column");

			else if (row == 0) {
				acol = a[0].length;
				if (a.length <= acol)
					throw new IllegalArgumentException(
							"matrix a must contain more rows than columns");

			} else if (acol != a[row].length)
				throw new IllegalArgumentException(
						"matrix a row width must remain constant");

			else if (1 != b[row].length)
				throw new IllegalArgumentException(
						"matrix b must only contain one column");

		// inline qr

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

		// inline matrix multiply and tranpose

		double[][] qtb = new double[acol][1];

		// find dot products of every row and column combination

		for (int row = 0; row < qtb.length; row++) {

			// inline dot product

			double sum = 0;

			for (int brow = 0; brow < b.length; brow++)
				sum += q[brow][row] * b[brow][0];

			qtb[row][0] = sum;
		}

		// inline gauss jordan

		// inline augment

		double[][] rqtb = new double[acol][acol + 1];
		for (int row = 0; row < acol; row++)
			for (int col = 0; col < acol; col++)
				rqtb[row][col] = r[row][col];

		for (int row = 0; row < acol; row++)
			rqtb[row][acol] = qtb[row][0];

		// slight modification to the order of declarations
		double[][] x = new double[acol][1];

		// inline rref

		if (rqtb.length == 0) {
			x = new double[][] {};

		} else if (rqtb.length == 1) {
			if (rqtb[0][0] == 0)
				x = new double[][] { new double[] { 0 } };
			else
				x = new double[][] { new double[] { 1 } };

		} else {
			// copy array ab into array c

			double[][] c = new double[rqtb[0].length - 1][rqtb[0].length];
			for (int row = 0; row < c.length; row++)
				for (int col = 0; col < rqtb[0].length; col++)
					c[row][col] = rqtb[row][col];

			// inline forward elimination

			for (int diag = 0; diag < c.length - 1; diag++) {
				if (c[diag][diag] == 0)
					throw new ArithmeticException("divide by 0");

				for (int row = diag + 1; row < c.length; row++) {
					double quotient = c[row][diag] / c[diag][diag];

					for (int col = 0; col < c[row].length; col++)
						c[row][col] -= (quotient * c[diag][col]);
				}
			}

			// inline backward elimination

			for (int diag = c.length - 1; diag > 0; diag--) {
				if (c[diag][diag] == 0)
					throw new ArithmeticException("divide by 0");

				for (int row = diag - 1; row >= 0; row--) {
					double quotient = c[row][diag] / c[diag][diag];

					for (int col = 0; col < c[row].length; col++)
						c[row][col] -= (quotient * c[diag][col]);
				}
			}

			// normalize c

			double[][] rref = new double[acol][rqtb[0].length];

			for (int row = 0; row < c.length; row++) {
				rref[row][rqtb[0].length - 1] = c[row][rqtb[0].length - 1]
						/ c[row][row];
				rref[row][row] = 1;
			}

			// copy last column of results into array x

			for (int row = 0; row < x.length; row++)
				x[row][0] = rref[row][rref[row].length - 1] / rref[row][row];
		}

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
	 * Assumes row by columnn array representation of matrix. Given that x is
	 * the first element of the resulting array and e is the second element of
	 * the resulting array, a*x - b = e where the magnitude of e is minimized.
	 * Returns array of array of n by 1 and m by 1 array of double,
	 * 
	 * This implementation finds x in one step. Given that array usv is the
	 * resulting array of sv decomposition on given matrix a where ut is
	 * tranpose of the first element of usv, (s^-1) is the inverse of the second
	 * element of usv, v is the third element of usv, x is equal to
	 * v*(s^-1)*ut*b.
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
	 * There is a special condition where divide by 0 occurs when the
	 * determinant of the resulting matrix from at*a equals 0.
	 * 
	 * @param a
	 *            m by n array of double where m > n
	 * @param b
	 *            m by 1 array of double
	 * @return array of n by 1 and m by 1 array of double
	 */
	public static double[][][] overconstrainedSV(double[][] a, double[][] b) {
		if (null == a)
			throw new IllegalArgumentException("matrix a must not be null");

		else if (null == b)
			throw new IllegalArgumentException("matrix b must not be null");

		else if (a.length != b.length)
			throw new IllegalArgumentException(
					"matrix a and b number of rows must be equal");

		int acol = -1;

		for (int row = 0; row < a.length; row++)
			if (null == a[row])
				throw new IllegalArgumentException(
						"matrix a must not contain null column");

			else if (null == b[row])
				throw new IllegalArgumentException(
						"matrix b must not contain null column");

			else if (row == 0) {
				acol = a[0].length;
				if (a.length <= acol)
					throw new IllegalArgumentException(
							"matrix a must contain more rows than columns");

			} else if (acol != a[row].length)
				throw new IllegalArgumentException(
						"matrix a row width must remain constant");

			else if (1 != b[row].length)
				throw new IllegalArgumentException(
						"matrix b must only contain one column");

		// inline singular value decomposition

		// inline matrix multiplication and tranpose

		double[][] ata = new double[acol][acol];

		for (int row = 0; row < ata.length; row++) {
			for (int col = 0; col < ata[row].length; col++) {

				// inline dot product

				double sum = 0;

				for (int brow = 0; brow < a.length; brow++)
					sum += a[brow][row] * a[brow][col];

				ata[row][col] = sum;
			}
		}

		// inline eigen

		double[][] l = new double[ata.length][ata.length];
		double[][] v = new double[ata.length][ata.length];

		if (ata.length == 1) {

			// assumed double a = -1;
			double bb = ata[0][0];

			l[0][0] = bb;

			v[1][0] = 1;

		} else if (ata.length == 2) {

			// assumed double a = 1;
			double bb = -ata[0][0] - ata[1][1];
			double c = ata[0][0] * ata[1][1] - ata[0][1] * ata[1][0];
			double sqrt = Math.sqrt((bb * bb) - (4 * c));

			l[0][0] = (-bb + sqrt) / 2;
			l[1][1] = (-bb - sqrt) / 2;

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
			double bb = ata[0][0] + ata[1][1] + ata[2][2];
			double c = -(ata[0][0] * ata[1][1]) - (ata[0][0] * ata[2][2])
					- (ata[1][1] * ata[2][2]) + (ata[0][1] * ata[1][0])
					+ (ata[0][2] * ata[2][0]) + (ata[1][2] * ata[2][1]);
			double d = (ata[2][0] * ata[0][1] * ata[1][2])
					+ (ata[1][0] * ata[2][1] * ata[0][2])
					- (ata[0][0] * ata[2][1] * ata[1][2])
					- (ata[1][1] * ata[2][0] * ata[0][2])
					- (ata[2][2] * ata[1][0] * ata[0][1])
					+ (ata[0][0] * ata[1][1] * ata[2][2]);

			double f = (-(3 * c) - (bb * bb)) / 3;
			double g = (-(2 * bb * bb * bb) - (9 * bb * c) - (27 * d)) / 27;
			double h = (g * g / 4) + (f * f * f / 27);
			double i = Math.sqrt((g * g / 4) - h);
			double j = Math.pow(i, 1d / 3);
			double k = Math.acos(-g / (2 * i));

			double m = Math.cos(k / 3);
			double n = Math.sqrt(3) * Math.sin(k / 3);
			double p = (bb / 3);
			double q = 2 * j * Math.cos(k / 3) + (bb / 3);
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

		double[][] inverse = new double[ata.length][ata.length];

		for (int row = 0; row < ata.length; row++)
			for (int col = 0; col < ata.length; col++) {
				if (l[row][col] == 0)
					// prevent divides by 0
					inverse[row][col] = 0;

				else
					// inverse is simple inverse since sigma is diagonal
					inverse[row][col] = 1 / Math.sqrt(l[row][col]);
			}

		// inline matrix multiply

		double[][] av = new double[a.length][v[0].length];

		for (int row = 0; row < av.length; row++) {
			for (int col = 0; col < av[row].length; col++) {

				// inline dot product

				double sum = 0;

				for (int brow = 0; brow < v.length; brow++)
					sum += a[row][brow] * v[brow][col];

				av[row][col] = sum;
			}
		}

		// inline matrix multiply and tranpose

		double[][] ut = new double[inverse[0].length][av.length];

		for (int col = 0; col < ut.length; col++) {
			for (int row = 0; row < ut[col].length; row++) {

				// inline dot product

				double sum = 0;

				for (int brow = 0; brow < inverse.length; brow++)
					sum += av[row][brow] * inverse[brow][col];

				ut[col][row] = sum;
			}
		}

		// 27 variables NOT counting any variable declared inside for loop

		// inline matrix multiply

		double[][] vsInv = new double[acol][acol];

		// find dot products of every row and column combination

		for (int row = 0; row < vsInv.length; row++) {
			for (int col = 0; col < vsInv[row].length; col++) {

				// inline dot product

				double sum = 0;

				for (int brow = 0; brow < acol; brow++)
					sum += v[row][brow] * inverse[brow][col];

				vsInv[row][col] = sum;
			}
		}

		// inline matrix multiply

		double[][] vsInvUt = new double[acol][a.length];

		// find dot products of every row and column combination

		for (int row = 0; row < vsInvUt.length; row++) {
			for (int col = 0; col < vsInvUt[row].length; col++) {

				// inline dot product

				double sum = 0;

				for (int brow = 0; brow < ut.length; brow++)
					sum += vsInv[row][brow] * ut[brow][col];

				vsInvUt[row][col] = sum;
			}
		}

		// inline matrix multiply

		double[][] x = new double[acol][1];

		// find dot products of every row and column combination

		for (int row = 0; row < x.length; row++) {
			for (int col = 0; col < x[row].length; col++) {

				// inline dot product

				double sum = 0;

				for (int brow = 0; brow < b.length; brow++)
					sum += vsInvUt[row][brow] * b[brow][col];

				x[row][col] = sum;
			}
		}

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
	 * Assumes row by columnn array representation of matrix. Given that x is
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
	 * There is a special condition where divide by 0 occurs when the
	 * determinant of the resulting matrix from a*at equals 0.
	 * 
	 * @param a
	 *            m by n array of double where m < n
	 * @param b
	 *            m by 1 array of double
	 * @return n by 1 array of double
	 */
	public static double[][] underconstrainedSV(double[][] a, double[][] b) {
		if (null == a)
			throw new IllegalArgumentException("matrix a must not be null");

		else if (null == b)
			throw new IllegalArgumentException("matrix b must not be null");

		else if (a.length != b.length)
			throw new IllegalArgumentException(
					"matrix a and b number of rows must be equal");

		int acol = -1;

		for (int row = 0; row < a.length; row++)
			if (null == a[row])
				throw new IllegalArgumentException(
						"matrix a must not contain null column");

			else if (null == b[row])
				throw new IllegalArgumentException(
						"matrix b must not contain null column");

			else if (row == 0) {
				acol = a[0].length;
				if (a.length >= acol)
					throw new IllegalArgumentException(
							"matrix a must contain less rows than columns");

			} else if (acol != a[row].length)
				throw new IllegalArgumentException(
						"matrix a row width must remain constant");

			else if (1 != b[row].length)
				throw new IllegalArgumentException(
						"matrix b must only contain one column");

		// inline matrix multiplication and tranpose (different tranpose)

		double[][] aat = new double[a.length][a.length];

		for (int row = 0; row < a.length; row++) {
			for (int col = 0; col < a.length; col++) {

				// inline dot product

				double sum = 0;

				for (int brow = 0; brow < acol; brow++)
					sum += a[row][brow] * a[col][brow];

				aat[row][col] = sum;
			}
		}

		// inline eigen

		double[][] l = new double[aat.length][aat.length];
		double[][] v = new double[aat.length][aat.length];

		if (aat.length == 1) {

			// assumed double a = -1;
			double bb = aat[0][0];

			l[0][0] = bb;

			v[1][0] = 1;

		} else if (aat.length == 2) {

			// assumed double a = 1;
			double bb = -aat[0][0] - aat[1][1];
			double c = aat[0][0] * aat[1][1] - aat[0][1] * aat[1][0];
			double sqrt = Math.sqrt((bb * bb) - (4 * c));

			l[0][0] = (-bb + sqrt) / 2;
			l[1][1] = (-bb - sqrt) / 2;

			for (int row = 0; row < aat.length; row++) {
				v[0][row] = -aat[0][1] / (aat[0][0] - l[row][row]);
				v[1][row] = 1;
			}

			for (int row = 0; row < aat.length; row++) {
				double mag = Math.sqrt((v[0][row] * v[0][row]) + 1);
				v[0][row] /= mag;
				v[1][row] /= mag;
			}

		} else if (aat.length == 3) {

			// assumed double a = -1;
			double bb = aat[0][0] + aat[1][1] + aat[2][2];
			double c = -(aat[0][0] * aat[1][1]) - (aat[0][0] * aat[2][2])
					- (aat[1][1] * aat[2][2]) + (aat[0][1] * aat[1][0])
					+ (aat[0][2] * aat[2][0]) + (aat[1][2] * aat[2][1]);
			double d = (aat[2][0] * aat[0][1] * aat[1][2])
					+ (aat[1][0] * aat[2][1] * aat[0][2])
					- (aat[0][0] * aat[2][1] * aat[1][2])
					- (aat[1][1] * aat[2][0] * aat[0][2])
					- (aat[2][2] * aat[1][0] * aat[0][1])
					+ (aat[0][0] * aat[1][1] * aat[2][2]);

			double f = (-(3 * c) - (bb * bb)) / 3;
			double g = (-(2 * bb * bb * bb) - (9 * bb * c) - (27 * d)) / 27;
			double h = (g * g / 4) + (f * f * f / 27);
			double i = Math.sqrt((g * g / 4) - h);
			double j = Math.pow(i, 1d / 3);
			double k = Math.acos(-g / (2 * i));

			double m = Math.cos(k / 3);
			double n = Math.sqrt(3) * Math.sin(k / 3);
			double p = (bb / 3);
			double q = 2 * j * Math.cos(k / 3) + (bb / 3);
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

			for (int row = 0; row < aat.length; row++) {

				// copy array sig into array x

				double[][] x = new double[aat.length - 1][aat.length];
				for (int xrow = 0; xrow < x.length; xrow++)
					for (int xcol = 0; xcol < aat.length; xcol++) {
						x[xrow][xcol] = aat[xrow][xcol];
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

			for (int row = 0; row < aat.length; row++) {

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

		double[][] inverse = new double[aat.length][aat.length];

		for (int row = 0; row < aat.length; row++)
			for (int col = 0; col < aat.length; col++) {
				if (l[row][col] == 0)
					// prevent divides by 0
					inverse[row][col] = 0;

				else
					// inverse is simple inverse since sigma is diagonal
					inverse[row][col] = 1 / Math.sqrt(l[row][col]);
			}

		// inline matrix multiply

		double[][] av = new double[acol][a.length];

		for (int row = 0; row < av.length; row++) {
			for (int col = 0; col < av[row].length; col++) {

				// inline dot product

				double sum = 0;

				for (int brow = 0; brow < v.length; brow++)
					sum += a[brow][row] * v[brow][col];

				av[row][col] = sum;
			}
		}

		// inline matrix multiply

		double[][] u = new double[acol][a.length];

		for (int row = 0; row < u.length; row++) {
			for (int col = 0; col < u[row].length; col++) {

				// inline dot product

				double sum = 0;

				for (int brow = 0; brow < inverse.length; brow++)
					sum += av[row][brow] * inverse[brow][col];

				u[row][col] = sum;
			}
		}

		// 27 variables NOT counting any variable declared inside for loop

		// inline matrix multiply

		double[][] usInv = new double[acol][a.length];

		// find dot products of every row and column combination

		for (int row = 0; row < usInv.length; row++) {
			for (int col = 0; col < usInv[row].length; col++) {

				// inline dot product

				double sum = 0;

				for (int brow = 0; brow < inverse.length; brow++)
					sum += u[row][brow] * inverse[brow][col];

				usInv[row][col] = sum;
			}
		}

		// inline matrix multiply and tranpose

		double[][] usInvVt = new double[acol][a.length];

		// find dot products of every row and column combination

		for (int row = 0; row < usInvVt.length; row++) {
			for (int col = 0; col < usInvVt[row].length; col++) {

				// inline dot product

				double sum = 0;

				for (int brow = 0; brow < a.length; brow++)
					sum += usInv[row][brow] * v[col][brow];

				usInvVt[row][col] = sum;
			}
		}

		// inline matrix multiply

		double[][] x = new double[acol][1];

		// find dot products of every row and column combination

		for (int row = 0; row < x.length; row++) {

			// inline dot product

			double sum = 0;

			for (int brow = 0; brow < b.length; brow++)
				sum += usInvVt[row][brow] * b[brow][0];

			x[row][0] = sum;
		}

		return x;
	}

}