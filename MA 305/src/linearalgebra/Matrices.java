package linearalgebra;

/**
 * Efficiently computes problems.
 * 
 * @author Jacob Malter
 *
 */
public class Matrices {

	/**
	 * Does nothing. Never called.
	 */
	private Matrices() {
	}

	/**
	 * Assumes row by columnn array representation of matrix. Returns amount of
	 * "sofu" enclosed in linear object defined by matrix a.
	 * 
	 * This implementation does not use recursion when number of rows is less
	 * than four.
	 * 
	 * Throws IllegalArgumentException if any of the following is true:
	 * {@code (null == a)}, {@code (null == a[row])} where a[row] is any row, or
	 * {@code (a.length != a[row].length)} where length is number of rows and
	 * a[row].length is number of columns on any row.
	 * 
	 * @param a
	 *            n by n array of double
	 * @return amount of "sofu" enclosed in linear object defined by matrix a
	 */
	public static double determinant(double[][] a) {
		if (null == a)
			throw new IllegalArgumentException("matrix a must not be null");

		for (int row = 0; row < a.length; row++)
			if (null == a[row])
				throw new IllegalArgumentException(
						"matrix a must not contain null column");

			else if (a.length != a[row].length)
				throw new IllegalArgumentException(
						"matrix a number of rows and columns must be equal");

		if (a.length == 0)
			// point with no dimensions
			return 0;

		else if (a.length == 1)
			// length of line
			return a[0][0];

		else if (a.length == 2)
			// area between vectors
			return a[0][0] * a[1][1] - a[0][1] * a[1][0];

		else if (a.length == 3)
			// volume inside triple product of vectors
			return (a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1]))
					- (a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0]))
					+ (a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]));

		double determinant = 0;

		for (int col = 0; col < a.length; col++) {
			if (col % 2 == 0)
				determinant += a[0][col] * determinant(minor(a, 0, col));

			else
				determinant -= a[0][col] * determinant(minor(a, 0, col));
		}

		// other dimensional "sofu"
		return determinant;
	}

	/**
	 * Returns sum of products of respective components.
	 * 
	 * Throws IllegalArgumentException if any of the following is true:
	 * {@code (null == a)}, {@code (null == b)}, or
	 * {@code (a.length != b.length)} where length is number of components.
	 * 
	 * @param a
	 *            n component array
	 * @param b
	 *            n component array
	 * @return sum of products of respective components
	 */
	public static double dotProduct(double[] a, double[] b) {
		if (null == a)
			throw new IllegalArgumentException("matrix a must not be null");

		else if (null == b)
			throw new IllegalArgumentException("matrix b must not be null");

		else if (a.length != b.length)
			throw new IllegalArgumentException(
					"matrix a and b number of components must be equal");

		double sum = 0;

		// multiply respective components and accumulate

		for (int cmpnt = 0; cmpnt < a.length; cmpnt++)
			sum += (a[cmpnt] * b[cmpnt]);

		return sum;
	}

	/**
	 * Assumes row by columnn array representation of matrix. Resulting matrix
	 * does not contain row at skipRow or column at skipCol of matrix a. Returns
	 * n-1 by m-1 array of double.
	 * 
	 * Throws IllegalArgumentException if any of the following is true:
	 * {@code (null == a)}, {@code (skipRow < 0)}, {@code (skipCol < 0)},
	 * {@code (skipRow >= a.length)} where length is number of rows,
	 * {@code (null == a[row])} where a[row] is any row,
	 * {@code (acol != a[row].length)} where acol is number of columns on first
	 * row and a[row].length is number of columns, or {@code (skipCol >= acol)}
	 * where acol is number of columns on first row.
	 * 
	 * @param a
	 *            n by m array of double
	 * @param skipRow
	 *            removed row from a
	 * @param skipCol
	 *            removed column from a
	 * @return n-1 by m-1 array of double
	 */
	public static double[][] minor(double[][] a, int skipRow, int skipCol) {
		if (null == a)
			throw new IllegalArgumentException("matrix a must not be null");

		else if (skipRow < 0)
			throw new IllegalArgumentException(
					"skipRow must not be less than zero");

		else if (skipCol < 0)
			throw new IllegalArgumentException(
					"skipCol must not be less than zero");

		else if (skipRow >= a.length)
			throw new IllegalArgumentException(
					"skipRow must not be greater than or equal to number of rows");

		int acol = -1;

		for (int row = 0; row < a.length; row++)
			if (null == a[row])
				throw new IllegalArgumentException(
						"matrix a must not contain null column");

			else if (row == 0)
				acol = a[0].length;

			else if (acol != a[row].length)
				throw new IllegalArgumentException(
						"matrix a row width must remain constant");

			else if (skipCol >= acol)
				throw new IllegalArgumentException(
						"skipCol must not be greater than or equal to number of columns");

		// minor has one less row and column than original

		double[][] minor = new double[a.length - 1][a.length - 1];

		for (int row = 0, srow = 0; row < minor.length; row++, srow++) {
			if (skipRow == row)
				// ignore row being skipped and access a from new row
				srow++;

			for (int col = 0, scol = 0; col < minor.length; col++, scol++) {
				if (skipCol == col)
					// ignore column being skipped and access a from new col
					scol++;
				minor[row][col] = a[srow][scol];
			}
		}

		return minor;
	}

	/**
	 * Assumes row by columnn array representation of matrix. Returns n by p
	 * array of double.
	 * 
	 * Throws IllegalArgumentException if any of the following is true:
	 * {@code (null == a)}, {@code (null == b)}, {@code (null == a[row])} where
	 * a[row] is any row, {@code (null == b[row])} where b[row] is any row,
	 * {@code (acol != a[row].length)} where acol is number of columns on first
	 * row and a[row].length is number of columns, or
	 * {@code (b.length != a[row].length)} where length is number of rows and
	 * a[row].length is number of columns.
	 * 
	 * @param a
	 *            n by m array of double
	 * @param b
	 *            m by p array of double
	 * @return n by p array of double
	 */
	public static double[][] multiply(double[][] a, double[][] b) {
		if (null == a)
			throw new IllegalArgumentException("matrix a must not be null");

		else if (null == b)
			throw new IllegalArgumentException("matrix b must not be null");

		int acol = -1;
		int bcol = -1;

		for (int row = 0; row < a.length; row++)
			if (null == a[row])
				throw new IllegalArgumentException(
						"matrix a must not contain null column");

			else if (row < b.length && null == b[row])
				throw new IllegalArgumentException(
						"matrix b must not contain null column");

			else if (row == 0) {
				acol = a[0].length;
				bcol = b[0].length;
			}

			else if (acol != a[row].length)
				throw new IllegalArgumentException(
						"matrix a row width must remain constant");

			else if (b.length != a[row].length)
				throw new IllegalArgumentException(
						"matrix b number of rows must equal matrix a number of columns");

		double[][] mul = new double[a.length][bcol];

		// find dot products of every row and column combination

		for (int row = 0; row < mul.length; row++) {
			for (int col = 0; col < mul[row].length; col++) {

				// inline dot product

				double sum = 0;

				for (int brow = 0; brow < b.length; brow++)
					sum += a[row][brow] * b[brow][col];

				mul[row][col] = sum;
			}
		}

		return mul;
	}

	/**
	 * Assumes row by columnn array representation of matrix. Rows become
	 * respective columns, and columns become respective rows. Alternatively,
	 * flip element of matrix a over its diagonal. Returns m by n array of
	 * double.
	 * 
	 * Throws IllegalArgumentException if any of the following is true:
	 * {@code (null == a)}, {@code (null == a[row])} where a[row] is any row, or
	 * {@code (acol != a[row].length)} where acol is number of columns on first
	 * row and a[row].length is number of columns.
	 * 
	 * @param a
	 *            n by m array of double
	 * @return m by n array of double
	 */
	public static double[][] transpose(double[][] a) {
		if (null == a)
			throw new IllegalArgumentException("matrix a must not be null");

		int acol = -1;

		for (int row = 0; row < a.length; row++)
			if (null == a[row])
				throw new IllegalArgumentException(
						"matrix a must not contain null column");

			else if (row == 0)
				acol = a[0].length;

			else if (acol != a[row].length)
				throw new IllegalArgumentException(
						"matrix a row width must remain constant");

		double[][] at = new double[acol][a.length];

		for (int row = 0; row < a.length; row++)
			for (int col = 0; col < acol; col++)
				// row becomes column; column becomes row
				at[col][row] = a[row][col];

		return at;
	}

	/**
	 * Assumes row by columnn array representation of matrix. Returns amount of
	 * "sofu" defined by n vectors enclosing a linear object in m-th dimensional
	 * space with m components.
	 * 
	 * Throws IllegalArgumentException if any of the following is true:
	 * {@link Matrices#transpose(double[][])}.
	 * 
	 * @param a
	 *            n by m array of double
	 * @return amount of "sofu" defined by n vectors enclosing a linear object
	 *         in m-th dimensional space with m components
	 */
	public static double sofu(double[][] a) {
		if (null == a)
			throw new IllegalArgumentException("matrix a must not be null");

		int acol = -1;

		for (int row = 0; row < a.length; row++)
			if (null == a[row])
				throw new IllegalArgumentException(
						"matrix a must not contain null column");

			else if (row == 0)
				acol = a[0].length;

			else if (acol != a[row].length)
				throw new IllegalArgumentException(
						"matrix a row width must remain constant");

		// inline matrix multiplication

		double[][] mul = new double[a[0].length][a[0].length];

		for (int row = 0; row < mul.length; row++) {
			for (int col = 0; col < mul[row].length; col++) {
				double sum = 0;

				for (int brow = 0; brow < a.length; brow++)
					sum += a[brow][row] * a[brow][col];

				mul[row][col] = sum;
			}
		}

		// inline determinant whose result is placed under square root

		if (mul.length == 0)
			return 0;

		else if (mul.length == 1)
			return Math.sqrt(mul[0][0]);

		else if (mul.length == 2)
			return Math.sqrt(mul[0][0] * mul[1][1] - mul[0][1] * mul[1][0]);

		else if (a.length == 3)
			return Math.sqrt((mul[0][0] * (mul[1][1] * mul[2][2] - mul[1][2]
					* mul[2][1]))
					- (mul[0][1] * (mul[1][0] * mul[2][2] - mul[1][2]
							* mul[2][0]))
					+ (mul[0][2] * (mul[1][0] * mul[2][1] - mul[1][1]
							* mul[2][0])));

		double determinant = 0;

		for (int col = 0; col < mul.length; col++) {
			if (col % 2 == 0)
				determinant += mul[0][col] * determinant(minor(mul, 0, col));

			else
				determinant -= mul[0][col] * determinant(minor(mul, 0, col));
		}

		return Math.sqrt(determinant);
	}

}