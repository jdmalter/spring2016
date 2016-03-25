package linearalgebra;

/**
 * Efficiently computes problems.
 * 
 * @author Jacob Malter
 *
 */
public class Extra {

	/**
	 * Does nothing. Never called.
	 */
	private Extra() {
	}

	/**
	 * Assumes row by columnn array representation of matrix. Given that x is
	 * the first element of the resulting array and e is the second element of
	 * the resulting array, a*x (through matrix multiplication) - b equals e
	 * where the magnitude of e is minimized. Returns array of array of n by 1
	 * and m by 1 array of double,
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
	 * determinant of the resulting matrix from at*a (through matrix
	 * multiplication) eqauls 0.
	 * 
	 * @param a
	 *            m by n array of double where m > n
	 * @param b
	 *            m by 1 array of double
	 * @return array of n by 1 and m by 1 array of double
	 */
	public static double[][][] overconstrainted(double[][] a, double[][] b) {
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
	 * resulting matrix, a*x (through matrix multiplication) equals b where the
	 * magnitude of x is minimized. Returns n by 1 array of double.
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
	 * determinant of the resulting matrix from a*at (through matrix
	 * multiplication) eqauls 0.
	 * 
	 * @param a
	 *            m by n array of double where m < n
	 * @param b
	 *            m by 1 array of double
	 * @return n by 1 array of double
	 */
	public static double[][] underconstrained(double[][] a, double[][] b) {
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

		double[][] atAdjugate = new double[a[0].length][aat.length];

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

}