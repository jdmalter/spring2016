package linearalgebra;

/**
 * Builds more matrix operations with detailed specifications.
 * 
 * @author Jacob Malter
 *
 */
public class Fundamentals {

	/**
	 * Does nothing. Never called.
	 */
	private Fundamentals() {
	}

	/**
	 * Assumes row by columnn array representation of matrix. Columns 0 through
	 * x of the resulting matrix contain array a, and columns x through y of the
	 * resulting matrix contain array b. Returns n by x+y array of double.
	 * 
	 * Throws IllegalArgumentException if any of the following is true:
	 * {@code (null == a)}, {@code (null == b)}, {@code (a.length != b.length)}
	 * where length is number of rows, {@code (null == a[row])} where a[row] is
	 * any row, {@code (null == b[row])} where b[row] is any row,
	 * {@code (acol != a[row].length)} where acol is number of columns on first
	 * row and a[row].length is number of columns, or
	 * {@code (bcol != b[row].length)} where bcol is number of columns on first
	 * row and b[row].length is number of columns on any row.
	 * 
	 * @param a
	 *            n by x array of double
	 * @param b
	 *            n by y array of double
	 * @return n by x+y array of double
	 */
	public static double[][] augment(double[][] a, double[][] b) {
		if (null == a)
			throw new IllegalArgumentException("matrix a must not be null");

		else if (null == b)
			throw new IllegalArgumentException("matrix b must not be null");

		else if (a.length != b.length)
			throw new IllegalArgumentException("matrix a and b number of rows must be equal");

		int acol = -1;
		int bcol = -1;

		for (int row = 0; row < a.length; row++)
			if (null == a[row])
				throw new IllegalArgumentException("matrix a must not contain null column");

			else if (null == b[row])
				throw new IllegalArgumentException("matrix b must not contain null column");

			else if (row == 0) {
				acol = a[0].length;
				bcol = b[0].length;
			}

			else if (acol != a[row].length)
				throw new IllegalArgumentException("matrix a row width must remain constant");

			else if (bcol != b[row].length)
				throw new IllegalArgumentException("matrix b row width must remain constant");

		double[][] aug = new double[a.length][acol + bcol];

		// fill in left rows of augment

		for (int row = 0; row < a.length; row++)
			for (int col = 0; col < acol; col++)
				aug[row][col] = a[row][col];

		// fill in right rows of augment

		for (int row = 0; row < b.length; row++)
			for (int col = 0; col < bcol; col++)
				aug[row][col + acol] = b[row][col];

		return aug;
	}

	/**
	 * Assumes row by columnn array representation of matrix. Creates new matrix
	 * equal to given matrix a except the lower triangular is set to 0. Returns
	 * n by m array of double.
	 * 
	 * Throws IllegalArgumentException if any of the following is true:
	 * {@code (null == a)}, {@code (null == a[row])} where a[row] is any row,
	 * {@code (a.length > acol)} where length is number rows and acol is number
	 * of columns on first row, or {@code (acol != a[row].length)} where acol is
	 * number of columns on first row and a[row].length is number of columns.
	 * 
	 * Throws ArithmeticException if there is 0 on diagonal.
	 * 
	 * @param a
	 *            n by m array of double
	 * @return n by m array of double
	 */
	public static double[][] forwardEliminate(double[][] a) {
		if (null == a)
			throw new IllegalArgumentException("matrix a must not be null");

		int acol = -1;

		for (int row = 0; row < a.length; row++)
			if (null == a[row])
				throw new IllegalArgumentException("matrix a must not contain null column");

			else if (row == 0) {
				acol = a[0].length;
				if (a.length > acol)
					throw new IllegalArgumentException("matrix a must contain as many or more columns than rows");
			}

			else if (acol != a[row].length)
				throw new IllegalArgumentException("matrix a row width must remain constant");

		// copy array a into forward

		double[][] forward = new double[a.length][acol];
		for (int row = 0; row < forward.length; row++)
			for (int col = 0; col < forward[row].length; col++)
				forward[row][col] = a[row][col];

		for (int diag = 0; diag < a.length - 1; diag++) {
			if (forward[diag][diag] == 0)
				// prevent divide by 0 a few lines later
				throw new ArithmeticException("divide by 0");

			for (int row = diag + 1; row < a.length; row++) {
				double quotient = forward[row][diag] / forward[diag][diag];

				for (int col = 0; col < forward[row].length; col++)
					forward[row][col] -= (quotient * forward[diag][col]);
			}
		}

		return forward;
	}

	/**
	 * Assumes row by columnn array representation of matrix. Creates new matrix
	 * equal to given matrix a except the upper triangular is set to 0. Returns
	 * n by m array of double.
	 * 
	 * Throws IllegalArgumentException if any of the following is true:
	 * {@code (null == a)}, {@code (null == a[row])} where a[row] is any row,
	 * {@code (a.length > acol)} where length is number rows and acol is number
	 * of columns on first row, or {@code (acol != a[row].length)} where acol is
	 * number of columns on first row and a[row].length is number of columns.
	 * 
	 * Throws ArithmeticException if there is 0 on diagonal.
	 * 
	 * @param a
	 *            n by m array of double
	 * @return n by m array of double
	 */
	public static double[][] backwardEliminate(double[][] a) {
		if (null == a)
			throw new IllegalArgumentException("matrix a must not be null");

		int acol = -1;

		for (int row = 0; row < a.length; row++)
			if (null == a[row])
				throw new IllegalArgumentException("matrix a must not contain null column");

			else if (row == 0) {
				acol = a[0].length;
				if (a.length > acol)
					throw new IllegalArgumentException("matrix a must contain as many or more columns than rows");
			}

			else if (acol != a[row].length)
				throw new IllegalArgumentException("matrix a row width must remain constant");

		// copy a into backward

		double[][] backward = new double[a.length][acol];
		for (int row = 0; row < backward.length; row++)
			for (int col = 0; col < backward[row].length; col++)
				backward[row][col] = a[row][col];

		for (int diag = a.length - 1; diag > 0; diag--) {
			if (backward[diag][diag] == 0)
				// prevent divide by 0 a few lines later
				throw new ArithmeticException("divide by 0");

			for (int row = diag - 1; row >= 0; row--) {
				double quotient = backward[row][diag] / backward[diag][diag];

				for (int col = 0; col < backward[row].length; col++)
					backward[row][col] -= (quotient * backward[diag][col]);
			}
		}

		return backward;
	}

	/**
	 * Assumes row by columnn array representation of matrix. Given that row
	 * echlon form is where all zeros belong to the highest rows and lowest
	 * columns and the leftmost nonzero number is in a higher column that the
	 * leftmost nonzero number in the column above, the resulting matrix is in
	 * row echelon form and every row's leftmost number is 1 which is the only
	 * element in its column. Returns n by m array of double.
	 * 
	 * Throws IllegalArgumentException if any of the following is true:
	 * {@code (null == a)}, {@code (null == a[row])} where a[row] is any row,
	 * {@code (a.length > acol)} where length is number rows and acol is number
	 * of columns on first row, or {@code (acol != a[row].length)} where acol is
	 * number of columns on first row and a[row].length is number of columns.
	 * 
	 * @param a
	 *            n by m array of double
	 * @return n by m array of double
	 */
	public static double[][] rref(double[][] a) {
		if (null == a)
			throw new IllegalArgumentException("matrix a must not be null");

		int acol = -1;

		for (int row = 0; row < a.length; row++)
			if (null == a[row])
				throw new IllegalArgumentException("matrix a must not contain null column");

			else if (row == 0) {
				acol = a[0].length;
				if (a.length > acol)
					throw new IllegalArgumentException("matrix a must contain as many or more columns than rows");

			} else if (acol != a[row].length)
				throw new IllegalArgumentException("matrix a row width must remain constant");

		if (a.length == 0)
			// there are no numbers to work with
			return new double[][] {};

		else if (a.length == 1)
			if (a[0][0] == 0)
				// prevent 0 from dividing by itself
				return new double[][] { new double[] { 0 } };
			else
				// otherwise, result is always 1
				return new double[][] { new double[] { 1 } };

		// copy array a into array c

		double[][] c = new double[acol - 1][acol];
		for (int row = 0; row < c.length; row++)
			for (int col = 0; col < acol; col++)
				c[row][col] = a[row][col];

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

		double[][] rref = new double[a.length][acol];

		for (int row = 0; row < c.length; row++) {
			rref[row][acol - 1] = c[row][acol - 1] / c[row][row];
			rref[row][row] = 1;
		}

		return rref;
	}

	/**
	 * Assumes row by columnn array representation of matrix. Given that x is
	 * the resulting matrix, a*x - b = 0. Returns n by 1 array of double.
	 * 
	 * Throws IllegalArgumentException if any of the following is true:
	 * {@code (null == a)}, {@code (null == b)}, {@code (a.length != b.length)}
	 * where length is number of rows, {@code (null == a[row])} where a[row] is
	 * any row, {@code (null == b[row])} where b[row] is any row,
	 * {@code (acol != a[row].length)} where acol is number of columns on first
	 * row and a[row].length is number of columns,
	 * {@code (bcol != b[row].length)} where bcol is number of columns on first
	 * row and b[row].length is number of columns on any row,
	 * {@code (a.length != a[row].length)} where length is number of rows and
	 * a[row].length is number of columns on any row, or
	 * {@code (1 != b[row].length)} where b[row].length is number of columns on
	 * any row.
	 * 
	 * Throws ArithmeticException if there are no solutions or many solutions
	 * for x.
	 * 
	 * @param a
	 *            n by n array of double
	 * @param b
	 *            n by 1 array of double
	 * @return n by 1 array of double
	 */
	public static double[][] gaussJordan(double[][] a, double[][] b) {
		if (null == a)
			throw new IllegalArgumentException("matrix a must not be null");

		else if (null == b)
			throw new IllegalArgumentException("matrix b must not be null");

		else if (a.length != b.length)
			throw new IllegalArgumentException("matrix a and b number of rows must be equal");

		for (int row = 0; row < a.length; row++)
			if (null == a[row])
				throw new IllegalArgumentException("matrix a must not contain null column");

			else if (null == b[row])
				throw new IllegalArgumentException("matrix b must not contain null column");

			else if (a.length != a[row].length)
				throw new IllegalArgumentException("matrix a number of rows and columns must be equal");

			else if (1 != b[row].length)
				throw new IllegalArgumentException("matrix b must contain only one column");

		// inline augment

		double[][] ab = new double[a.length][a[0].length + b[0].length];

		for (int row = 0; row < a.length; row++)
			for (int col = 0; col < a[0].length; col++)
				ab[row][col] = a[row][col];

		for (int row = 0; row < b.length; row++)
			for (int col = 0; col < b[0].length; col++)
				ab[row][col + a[0].length] = b[row][col];

		// inline rref

		if (ab.length == 0)
			return new double[][] {};

		else if (ab.length == 1)
			if (ab[0][0] == 0)
				return new double[][] { new double[] { 0 } };
			else
				return new double[][] { new double[] { 1 } };

		// copy array ab into array c

		double[][] c = new double[ab[0].length - 1][ab[0].length];
		for (int row = 0; row < c.length; row++)
			for (int col = 0; col < ab[0].length; col++)
				c[row][col] = ab[row][col];

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

		double[][] rref = new double[a.length][ab[0].length];

		for (int row = 0; row < c.length; row++) {
			rref[row][ab[0].length - 1] = c[row][ab[0].length - 1] / c[row][row];
			rref[row][row] = 1;
		}

		// copy last column of results into array x

		double[][] x = new double[a.length][1];
		for (int row = 0; row < x.length; row++)
			x[row][0] = rref[row][rref[row].length - 1] / rref[row][row];

		return x;
	}

	/**
	 * Assumes row by columnn array representation of matrix. In the first
	 * element in the resulting array, each column, or variable, besides the
	 * furthest right column, has a lead variable of 1. Above the 1 in each
	 * column, any variable may have non-zero numbers which can be multipled by
	 * any value of their respective variable while maintaining row echlon form.
	 * Supports multiple solutions for linear equations. Given that x is the
	 * combination of the two matrices in the resulting array, a*x - b = 0.
	 * Returns array of n by n and n by 1 array of double.
	 * 
	 * Throws IllegalArgumentException if any of the following is true:
	 * {@code (null == a)}, {@code (null == b)}, {@code (a.length != b.length)}
	 * where length is number of rows, {@code (null == a[row])} where a[row] is
	 * any row, {@code (null == b[row])} where b[row] is any row,
	 * {@code (acol != a[row].length)} where acol is number of columns on first
	 * row and a[row].length is number of columns,
	 * {@code (bcol != b[row].length)} where bcol is number of columns on first
	 * row and b[row].length is number of columns on any row,
	 * {@code (a.length != a[row].length)} where length is number of rows and
	 * a[row].length is number of columns on any row, or
	 * {@code (1 != b[row].length)} where b[row].length is number of columns on
	 * any row.
	 * 
	 * Throws ArithmeticException if there are no solutions for x.
	 * 
	 * @param a
	 *            n by n array of double
	 * @param b
	 *            n by 1 array of double
	 * @return array of n by n and n by 1 array of double
	 */
	public static double[][][] gaussJordanMany(double[][] a, double[][] b) {
		if (null == a)
			throw new IllegalArgumentException("matrix a must not be null");

		else if (null == b)
			throw new IllegalArgumentException("matrix b must not be null");

		else if (a.length != b.length)
			throw new IllegalArgumentException("matrix a and b number of rows must be equal");

		for (int row = 0; row < a.length; row++)
			if (null == a[row])
				throw new IllegalArgumentException("matrix a must not contain null column");

			else if (null == b[row])
				throw new IllegalArgumentException("matrix b must not contain null column");

			else if (a.length != a[row].length)
				throw new IllegalArgumentException("matrix a number of rows and columns must be equal");

			else if (1 != b[row].length)
				throw new IllegalArgumentException("matrix b must contain only one column");

		// inline augment

		double[][] ab = new double[a.length][a[0].length + b[0].length];

		for (int row = 0; row < a.length; row++)
			for (int col = 0; col < a[0].length; col++)
				ab[row][col] = a[row][col];

		for (int row = 0; row < b.length; row++)
			for (int col = 0; col < b[0].length; col++)
				ab[row][col + a[0].length] = b[row][col];

		// keep track of columns ignore in computation

		int zeroCol = 0;

		// inline forward elimination

		for (int diag = 0; diag < a.length - 1; diag++) {
			if (ab[diag][diag] == 0) {
				if (ab[diag][a[0].length + b[0].length - 1] != 0)
					// implies 0 != 0 therefore throw a fit
					throw new ArithmeticException("No solution");

				zeroCol++;
				continue;
			}

			for (int row = diag + 1; row < a.length; row++) {
				double quotient = ab[row][diag] / ab[diag][diag];

				for (int col = 0; col < ab[row].length; col++)
					ab[row][col] -= (quotient * ab[diag][col]);
			}
		}

		// inline backward elimination

		for (int diag = a.length - 1 - zeroCol; diag > 0; diag--) {
			if (ab[diag][diag] == 0) {
				if (ab[diag][a[0].length + b[0].length - 1] != 0)
					// implies 0 != 0 therefore throw a fit
					throw new ArithmeticException("No solution");

				zeroCol++;
				continue;
			}

			for (int row = diag - 1; row >= 0; row--) {
				double quotient = ab[row][diag] / ab[diag][diag];

				for (int col = 0; col < ab[row].length; col++)
					ab[row][col] -= (quotient * ab[diag][col]);
			}
		}

		// normalize and deal with ignored columns

		// copy last column of results into array x

		double[][] x = new double[a.length][1];

		for (int row = 0; row < ab.length; row++) {
			if (row < ab.length - zeroCol) {

				// columns not being ignored

				for (int col = 0; col < a.length; col++)
					if (row != col)
						// avoid changing diagonal or last element
						ab[row][col] /= -ab[row][row];

				x[row][0] = ab[row][a.length] / ab[row][row];
				ab[row][row] = 1;

			} else {

				// columns being ignored

				for (int col = 0; col < ab[row].length; col++)
					if (row != col)
						ab[row][col] = 0;
					else
						ab[row][row] = 1;
			}
		}

		// Remove last column of ab

		double[][] copyAB = new double[a.length][a.length];

		for (int row = 0; row < a.length; row++)
			for (int col = 0; col < a.length; col++)
				copyAB[row][col] = ab[row][col];

		return new double[][][] { copyAB, x };
	}

	/**
	 * Assumes row by columnn array representation of matrix. Given that x is
	 * the resulting matrix, a*x - b = 0. Returns n by 1 array of double.
	 * 
	 * Throws IllegalArgumentException if any of the following is true:
	 * {@code (null == a)}, {@code (null == b)}, {@code (a.length != b.length)}
	 * where length is number of rows, {@code (null == a[row])} where a[row] is
	 * any row, {@code (null == b[row])} where b[row] is any row,
	 * {@code (acol != a[row].length)} where acol is number of columns on first
	 * row and a[row].length is number of columns, or
	 * {@code (bcol != b[row].length)} where bcol is number of columns on first
	 * row and b[row].length is number of columns on any row.
	 * 
	 * @param a
	 *            n by n array of double
	 * @param b
	 *            n by 1 array of double
	 * @return n by 1 array of double
	 */
	public static double[][] cramer(double[][] a, double[][] b) {
		if (null == a)
			throw new IllegalArgumentException("matrix a must not be null");

		else if (null == b)
			throw new IllegalArgumentException("matrix b must not be null");

		else if (a.length != b.length)
			throw new IllegalArgumentException("matrix a and b number of rows must be equal");

		for (int row = 0; row < a.length; row++)
			if (null == a[row])
				throw new IllegalArgumentException("matrix a must not contain null column");

			else if (null == b[row])
				throw new IllegalArgumentException("matrix b must not contain null column");

			else if (a.length != a[row].length)
				throw new IllegalArgumentException("matrix a number of rows and columns must be equal");

			else if (1 != b[row].length)
				throw new IllegalArgumentException("matrix b must contain only one column");

		double denominator = Matrices.determinant(a);

		double[][] cramer = new double[a.length][1];

		for (int col = 0; col < a.length; col++) {
			double[][] copy = new double[a.length][a.length];

			// copy array a into array copy
			for (int drow = 0; drow < copy.length; drow++)
				for (int dcol = 0; dcol < copy.length; dcol++)
					copy[drow][dcol] = a[drow][dcol];

			// copy array b into array copy at appropriate row
			for (int row = 0; row < copy.length; row++)
				copy[row][col] = b[row][0];

			cramer[col][0] = Matrices.determinant(copy) / denominator;
		}

		return cramer;
	}

	/**
	 * Assumes row by columnn array representation of matrix. Given that x is
	 * the resulting matrix, a*x - I = 0 where I is identity matrix. Returns n
	 * by n array of double.
	 * 
	 * Throws IllegalArgumentException if any of the following is true:
	 * {@code (null == a)}, {@code (null == a[row])} where a[row] is any row, or
	 * {@code (a.length != a[row].length)} where length is number of rows and
	 * a[row].length is number of columns on any row.
	 * 
	 * Throws ArithmeticException if any of the following is true:
	 * {@code (determinant == 0)} where determinant is computed from matrix a.
	 * 
	 * @param a
	 *            n by n array of double
	 * @return n by n array of double
	 */
	public static double[][] invert(double[][] a) {
		if (null == a)
			throw new IllegalArgumentException("matrix a must not be null");

		for (int row = 0; row < a.length; row++)
			if (null == a[row])
				throw new IllegalArgumentException("matrix a must not contain null column");

			else if (a.length != a[row].length)
				throw new IllegalArgumentException("matrix a number of rows and columns must be equal");

		double[][] adjugate = new double[a.length][a.length];
		double determinant = 0;

		if (a.length == 0) {

			// inline matrix determinant on point

			// matrix is singular
			throw new ArithmeticException("divide by 0");

		} else if (a.length == 1) {

			// inline matrix determinant on line

			determinant = a[0][0];
			if (determinant == 0)
				// matrix is singular
				throw new ArithmeticException("divide by 0");

			adjugate[0][0] = 1;

		} else if (a.length == 2) {

			// inline matrix determinant on vectors

			determinant = a[0][0] * a[1][1] - a[0][1] * a[1][0];
			if (determinant == 0)
				// matrix is singular
				throw new ArithmeticException("divide by 0");

			adjugate[0][0] = a[1][1];
			adjugate[0][1] = -a[0][1];
			adjugate[1][0] = -a[1][0];
			adjugate[1][1] = a[0][0];

		} else if (a.length == 3) {

			// inline matrix determinant on triple product of vectors

			determinant = (a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1]))
					- (a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0]))
					+ (a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]));
			if (determinant == 0)
				// matrix is singular
				throw new ArithmeticException("divide by 0");

			// multiple inline matrix determinants on vectors

			adjugate[0][0] = (a[1][1] * a[2][2] - a[1][2] * a[2][1]);
			adjugate[0][1] = -(a[0][1] * a[2][2] - a[0][2] * a[2][1]);
			adjugate[0][2] = (a[0][1] * a[1][2] - a[0][2] * a[1][1]);
			adjugate[1][0] = -(a[1][0] * a[2][2] - a[1][2] * a[2][0]);
			adjugate[1][1] = (a[0][0] * a[2][2] - a[0][2] * a[2][0]);
			adjugate[1][2] = -(a[0][0] * a[1][2] - a[0][2] * a[1][0]);
			adjugate[2][0] = (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
			adjugate[2][1] = -(a[0][0] * a[2][1] - a[0][1] * a[2][0]);
			adjugate[2][2] = (a[0][0] * a[1][1] - a[0][1] * a[1][0]);

		} else {

			determinant = Matrices.determinant(a);
			if (determinant == 0)
				// matrix is singular
				throw new ArithmeticException("divide by 0");

			for (int row = 0; row < a.length; row++)
				for (int col = 0; col < a.length; col++) {

					// inline matrix minor

					double[][] minor = new double[a.length - 1][a.length - 1];

					for (int mrow = 0, srow = 0; mrow < minor.length; mrow++, srow++) {
						if (row == mrow)
							// ignore row being skipped and access new row
							srow++;

						for (int mcol = 0, scol = 0; mcol < minor.length; mcol++, scol++) {
							if (col == mcol)
								// ignore column being skipped and access new
								// col
								scol++;

							minor[mrow][mcol] = a[srow][scol];
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

		for (int row = 0; row < a.length; row++)
			for (int col = 0; col < a.length; col++)
				adjugate[row][col] /= determinant;

		return adjugate;
	}

}