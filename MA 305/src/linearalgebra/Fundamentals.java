package linearalgebra;

/**
 * Efficiently computes problems.
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

	// Related to Gauss Jordan method

	/**
	 * Assumes row by columnn array representation of matrix. Rows of elements
	 * remain unchanged. Columns 0 (inclusive) through x (exclusive) of
	 * resulting matrix contain elements of array a, and columns x (inclusive)
	 * through y (exclusive) of resulting matrix contain elements of array b.
	 * Returns n by x+y array of double.
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
			throw new IllegalArgumentException(
					"matrix a and b number of rows must be equal");

		int acol = -1;
		int bcol = -1;

		for (int row = 0; row < a.length; row++)
			if (null == a[row])
				throw new IllegalArgumentException(
						"matrix a must not contain null column");

			else if (null == b[row])
				throw new IllegalArgumentException(
						"matrix b must not contain null column");

			else if (row == 0) {
				acol = a[0].length;
				bcol = b[0].length;
			}

			else if (acol != a[row].length)
				throw new IllegalArgumentException(
						"matrix a row width must remain constant");

			else if (bcol != b[row].length)
				throw new IllegalArgumentException(
						"matrix b row width must remain constant");

		double[][] aug = new double[a.length][acol + bcol];

		for (int row = 0; row < a.length; row++)
			for (int col = 0; col < acol; col++)
				aug[row][col] = a[row][col];

		for (int row = 0; row < b.length; row++)
			for (int col = 0; col < bcol; col++)
				aug[row][col + acol] = b[row][col];

		return aug;
	}

	/**
	 * Assumes row by columnn array representation of matrix. Sets lower
	 * triangular of a to 0. Returns n by m array of double.
	 * 
	 * Throws IllegalArgumentException if any of the following is true:
	 * {@code (null == a)}, {@code (null == a[row])} where a[row] is any row,
	 * {@code (a.length > acol)} where length is number rows and acol is number
	 * of columns on first row, or {@code (acol != a[row].length)} where acol is
	 * number of columns on first row and a[row].length is number of columns.
	 * 
	 * There is a special condition where ArithmeticException if divide by 0
	 * occurs.
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
				throw new IllegalArgumentException(
						"matrix a must not contain null column");

			else if (row == 0) {
				acol = a[0].length;
				if (a.length > acol)
					throw new IllegalArgumentException(
							"matrix a must contain as many or more columns than rows");
			}

			else if (acol != a[row].length)
				throw new IllegalArgumentException(
						"matrix a row width must remain constant");

		double[][] forward = new double[a.length][acol];
		for (int row = 0; row < forward.length; row++)
			for (int col = 0; col < forward[row].length; col++)
				forward[row][col] = a[row][col];

		int diagonal = a.length;
		for (int diag = 0; diag < diagonal - 1; diag++) {
			double divisor = forward[diag][diag];
			if (divisor == 0)
				throw new ArithmeticException("divide by 0");

			for (int row = diag + 1; row < diagonal; row++) {
				double dividend = forward[row][diag];
				double quotient = dividend / divisor;

				for (int col = 0; col < forward[row].length; col++)
					forward[row][col] -= (quotient * forward[diag][col]);
			}
		}

		return forward;
	}

	/**
	 * Assumes row by columnn array representation of matrix. Sets upper
	 * triangular of a to 0. Returns n by m array of double.
	 * 
	 * Throws IllegalArgumentException if any of the following is true:
	 * {@code (null == a)}, {@code (null == a[row])} where a[row] is any row,
	 * {@code (a.length > acol)} where length is number rows and acol is number
	 * of columns on first row, or {@code (acol != a[row].length)} where acol is
	 * number of columns on first row and a[row].length is number of columns.
	 * 
	 * There is a special condition where ArithmeticException if divide by 0
	 * occurs.
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
				throw new IllegalArgumentException(
						"matrix a must not contain null column");

			else if (row == 0) {
				acol = a[0].length;
				if (a.length > acol)
					throw new IllegalArgumentException(
							"matrix a must contain as many or more columns than rows");
			}

			else if (acol != a[row].length)
				throw new IllegalArgumentException(
						"matrix a row width must remain constant");

		double[][] backward = new double[a.length][acol];
		for (int row = 0; row < backward.length; row++)
			for (int col = 0; col < backward[row].length; col++)
				backward[row][col] = a[row][col];

		int diagonal = a.length;
		for (int diag = diagonal - 1; diag > 0; diag--) {
			double divisor = backward[diag][diag];
			if (divisor == 0)
				throw new ArithmeticException("divide by 0");

			for (int row = diag - 1; row >= 0; row--) {
				double dividend = backward[row][diag];
				double quotient = dividend / divisor;

				for (int col = 0; col < backward[row].length; col++)
					backward[row][col] -= (quotient * backward[diag][col]);
			}
		}

		return backward;
	}

	/**
	 * Assumes row by columnn array representation of matrix. Given that x is
	 * resulting matrix, a*x (through matrix multiplication) equals b. Returns n
	 * by 1 array of double.
	 * 
	 * Throws IllegalArgumentException if any of the following is true:
	 * {@link GaussJordan#augment(double[][], double[][])},
	 * {@code (a.length != a[row].length)} where length is number of rows and
	 * a[row].length is number of columns on any row, or
	 * {@code (1 != b[row].length)} where b[row].length is number of columns on
	 * any row.
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
			throw new IllegalArgumentException(
					"matrix a and b number of rows must be equal");

		for (int row = 0; row < a.length; row++)
			if (null == a[row])
				throw new IllegalArgumentException(
						"matrix a must not contain null column");

			else if (null == b[row])
				throw new IllegalArgumentException(
						"matrix b must not contain null column");

			else if (a.length != a[row].length)
				throw new IllegalArgumentException(
						"matrix a number of rows and columns must be equal");

			else if (1 != b[row].length)
				throw new IllegalArgumentException(
						"matrix b must contain only one column");

		double[][] ab = backwardEliminate(forwardEliminate(augment(a, b)));

		double[][] x = new double[a.length][1];
		for (int row = 0; row < x.length; row++)
			x[row][0] = ab[row][ab[row].length - 1] / ab[row][row];

		return x;
	}

	// Related to Cramer's method

	/**
	 * Assumes row by columnn array representation of matrix. Given that x is
	 * resulting matrix, a*x (through matrix multiplication) equals b. Returns n
	 * by 1 array of double.
	 * 
	 * Throws IllegalArgumentException if any of the following is true:
	 * {@link GaussJordan#augment(double[][], double[][])}.
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
			throw new IllegalArgumentException(
					"matrix a and b number of rows must be equal");

		for (int row = 0; row < a.length; row++)
			if (null == a[row])
				throw new IllegalArgumentException(
						"matrix a must not contain null column");

			else if (null == b[row])
				throw new IllegalArgumentException(
						"matrix b must not contain null column");

			else if (a.length != a[row].length)
				throw new IllegalArgumentException(
						"matrix a number of rows and columns must be equal");

			else if (1 != b[row].length)
				throw new IllegalArgumentException(
						"matrix b must contain only one column");

		double denominator = Matrices.determinant(a);

		double[][] cramer = new double[a.length][1];
		for (int col = 0; col < a.length; col++) {

			double[][] duplicate = new double[a.length][a.length];
			for (int drow = 0; drow < duplicate.length; drow++)
				for (int dcol = 0; dcol < duplicate.length; dcol++)
					duplicate[drow][dcol] = a[drow][dcol];

			for (int row = 0; row < duplicate.length; row++)
				duplicate[row][col] = b[row][0];

			double numerator = Matrices.determinant(duplicate);
			cramer[col][0] = numerator / denominator;
		}

		return cramer;
	}

	// Related to matrix inversion

	/**
	 * Assumes row by columnn array representation of matrix. Given that x is
	 * resulting matrix, a*x (through matrix multiplication) equals I where I is
	 * identity matrix. Returns n by n array of double.
	 * 
	 * Throws IllegalArgumentException if any of the following is true:
	 * {@link Matrices#determinant(double[][])}.
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
				throw new IllegalArgumentException(
						"matrix a must not contain null column");

			else if (a.length != a[row].length)
				throw new IllegalArgumentException(
						"matrix a number of rows and columns must be equal");

		double determinant = Matrices.determinant(a);
		if (determinant == 0)
			throw new ArithmeticException("divide by 0");

		double[][] adjugate = new double[a.length][a.length];
		for (int row = 0; row < a.length; row++)
			for (int col = 0; col < a.length; col++) {
				double[][] minor = Matrices.minor(a, row, col);

				if ((row + col) % 2 == 0)
					adjugate[col][row] = Matrices.determinant(minor);

				else
					adjugate[col][row] = -Matrices.determinant(minor);
			}

		for (int row = 0; row < a.length; row++)
			for (int col = 0; col < a.length; col++)
				adjugate[row][col] *= (1 / determinant);

		return adjugate;
	}

}