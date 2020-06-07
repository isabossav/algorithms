/*
**
** SeamCarver class implements a content-aware resizing technique where the image is reduced in size
** by one pixel of height (or width) at a time. SeamCarver is based on the seam-carving content-aware resizing technique,
** which unlike standard content-agnostic resizing techniques (such as cropping and scaling), preserves the most 
** interesting features of the image (such as aspect ratio, set of objects presents, etc.)
** 
*/

import edu.princeton.cs.algs4.Picture;
import java.awt.Color;
import java.util.Arrays;

public class SeamCarver {
	private int[][] color;
	private double[][] energy;
	private int width;
	private int height;

	// create a seam carver object based on the given picture
	public SeamCarver(Picture picture) {
		if (picture == null) throw new IllegalArgumentException("Invalid argument");
		Picture pic = new Picture(picture);
		width = pic.width();
		height = pic.height();
		
		// Temporarily store color information
		color = new int[height][width];
		for (int i = 0; i < height; i++) 
			for (int j = 0; j < width; j++) 
				color[i][j] = pic.getRGB(j, i);
			
		// Calculate energy of each pixel
		energy = new double[height][width];
		for (int i = 0; i < width; i++) {
			for (int j = 0; j < height; j++) {
				double delta = delta(i, j, width, height);
				energy[j][i] = delta;
			}
		}
	}

	// Return (Dx)^2
	private double delta(int col, int row, int w, int h) {

		// default delta - use it for border pixels
		double delta = 1000;

		// Check if it's an inner or a border pixel
		boolean xborder = col == 0 || col == w - 1;
		boolean yborder = row == 0 || row == h - 1;

		if (xborder || yborder) 
			return delta;
		
		// Calculate dx sq
		int rgb1 = color[row][col-1];
		int rgb2 = color[row][col+1];

		int rx = ((rgb1 >> 16) & 0xFF) - ((rgb2 >> 16) & 0xFF);
		int gx = ((rgb1 >> 8) & 0xFF) - ((rgb2 >> 8) & 0xFF);
		int bx = ((rgb1 >> 0) & 0xFF) - ((rgb2 >> 0) & 0xFF);

		rx *= rx;
		gx *= gx;
		bx *= bx;
		int dx = rx + gx + bx;

		// Calculate dy sq
		rgb1 = color[row-1][col];
		rgb2 = color[row+1][col];
		int ry = ((rgb1 >> 16) & 0xFF) - ((rgb2 >> 16) & 0xFF);
		int gy = ((rgb1 >> 8) & 0xFF) - ((rgb2 >> 8) & 0xFF);
		int by = ((rgb1 >> 0) & 0xFF) - ((rgb2 >> 0) & 0xFF);

		ry *= ry;
		gy *= gy;
		by *= by;
		int dy = ry + gy + by;

		// Calculate delta
		delta = Math.sqrt(dx + dy);
		return delta;
	}


	// current picture
	public Picture picture() {

		Picture pic = new Picture(width, height);
		for (int i = 0; i < height; i++) 
			for (int j = 0; j < width; j++) 
				pic.setRGB(j, i, color[i][j]);

		return pic;
	}

	// width of current picture
	public int width() {
		return width;
	}

	// height of current picture
	public int height() {
		return height;
	}

	// energy of pixel at column x and row y
	public double energy(int x, int y) {
		if (x < 0 || x >= width) throw new IllegalArgumentException("Invalid column\n");
		if (y < 0 || y >= height) throw new IllegalArgumentException("Invalid row\n");
		return energy[y][x];
	}


	// sequence of indices for horizontal seam
	public int[] findHorizontalSeam() {

		// Special case, there is only one column
		if (width == 1) return new int[]{0}; 

		double[][] distTo = new double[height][width];
		int[][] comingFrom = new int[height][width];

		// Min energy of horizontal seam and the row where it's at
		double minEnd = Double.POSITIVE_INFINITY;
		int minRow = -1;

		// Find minimum path using topological order
		// Left to right, top to bottom
		int borderval = 1000;
		for (int j = 0; j < width; j++) {
			for (int i = 0; i < height; i++) {

				// Special cases: first, second and last column
				// If this is the first col, distTo is always 1000
				if (j == 0) {
					distTo[i][j] = borderval;
					comingFrom[i][j] = i;
				}

				// If this is the last col, simply add 1000 to the value of the row
				// before it (since last col will always have energy 1000)
				else if (j == width-1) {
					distTo[i][j] = distTo[i][j-1] + borderval;
					comingFrom[i][j] = i;

					if (distTo[i][j] < minEnd) {
						minEnd = distTo[i][j];
						minRow = i;
					}
				}

				// If this is the second col, simply add own energy to the value
				// of the col before it
				else if (j == 1) {
					distTo[i][j] = energy[i][j] + borderval;
					comingFrom[i][j] = i;
				}

				// Otherwise, find the path that minimizes distance to pixel
				else 
					findMinSeamH(i, j, distTo, comingFrom);
				
			}
		}

		int[] sequence = new int[width];
		int lastCol = width-1;
		sequence[lastCol] = minRow;
		for (int i = lastCol; i > 0; i--) 
			sequence[i-1] = comingFrom[sequence[i]][i];

		return sequence;
	}

	private void findMinSeamH(int row, int col, double[][] distTo, int[][] comingFrom) {

		// To reach (row, col), we can come from three pixels:
		// (row-1, col-1) (row, col-1) and (row+1, col-1)

		// Assume the min if (row, col-1) - the pixel right before it (to it's left)
		int minPrevious = row;
		double minDistance = distTo[row][col-1];

		// If there is a pixel in (row-1, col-1), find if distance would be 
		// smaller coming from there
		if (row > 0 && distTo[row-1][col-1] < minDistance) {
			minPrevious = row-1;
			minDistance = distTo[row-1][col-1];
		}

		// If there is a pixel in (row+1, col-1), find if distance would be
		// smaller coming from there
		if (row < height-1 && distTo[row+1][col-1] < minDistance) {
			minPrevious = row+1;
			minDistance = distTo[row+1][col-1];
		}

		double totalDistance = energy[row][col] + minDistance;
		distTo[row][col] = totalDistance;
		comingFrom[row][col] = minPrevious;
	}


	// sequence of indices for vertical seam {
	public int[] findVerticalSeam() {

		// Special case, there is only one row
		if (height == 1) return new int[]{0}; 

		double[][] distTo = new double[height][width];
		int[][] comingFrom = new int[height][width];

		// Min energy of vertical seam and the col where it's at
		double minEnd = Double.POSITIVE_INFINITY;
		int minCol = -1;

		// Find minimum path using topological order
		// Top to bottom, left to right
		int borderval = 1000;
		for (int i = 0; i < height; i++) {
			for (int j = 0; j < width; j++) {

				// Special cases: first, second and last row
				// If this is the first row, distTo is always 1000
				if (i == 0) {
					distTo[i][j] = borderval;
					comingFrom[i][j] = j; 
				}

				// If this is the last row, simply add 1000 to the value of the col 
				// above it (since last row will always have energy 1000)
				else if (i == height-1) {
					distTo[i][j] = distTo[i-1][j] + borderval;
					comingFrom[i][j] = j;

					if (distTo[i][j] < minEnd) {
						minEnd = distTo[i][j];
						minCol = j;
					}
				}

				// If this is the second row, simply add own energy to the value 
				// of the col above it (which is always 1000)
				else if (i == 1) {
					distTo[i][j] = energy[i][j] + borderval;
					comingFrom[i][j] = j;
				}

				// Otherwise, find the path that minimizes distance to pixel
				else 
					findMinSeam(i, j, distTo, comingFrom);
			}
		}

		int[] sequence = new int[height];
		int lastRow = height-1;
		sequence[lastRow] = minCol;
		for (int i = lastRow; i > 0; i--) 
			sequence[i-1] = comingFrom[i][sequence[i]];

		return sequence;
	}

	private void findMinSeam(int row, int col, double[][] distTo, int[][] comingFrom) {

		// To reach (i, j), we can come from three pixels:
		// (i-1, j-1), (i-1, j) and (i-1, j+1)

		// Assume the min is (i-1, j) - the pixel right above it
		int minPrevious = col;
		double minDistance = distTo[row-1][col];

		// If there is a pixel in (i-1, j-1), find if distance would be
		// smaller coming from there
		if (col > 0 && distTo[row-1][col-1] < minDistance) {
			minPrevious = col-1;
			minDistance = distTo[row-1][col-1];
		}

		// If there is a pixel in (i-1, j+1), find if distance would be
		// smaller coming from there
		if (col < width-1 && distTo[row-1][col+1] < minDistance) {
			minPrevious = col+1;
			minDistance = distTo[row-1][col+1];
		}

		double totalDistance = energy[row][col] + minDistance;
		distTo[row][col] = totalDistance;
		comingFrom[row][col] = minPrevious;
	}

	// remove horizontal seam from current picture
	public void removeHorizontalSeam(int[] seam) {
		if (height <= 1) throw new IllegalArgumentException("Cannot remove horizontal seams from picture");
		if (seam == null) throw new IllegalArgumentException("Null argument\n");
		if (seam.length != width) throw new IllegalArgumentException("Invalid seam length\n");

		int[] s = new int[seam.length];
		int previous = seam[0];
		for (int i = 0; i < seam.length; i++) {
			int row = seam[i];
			if (row < 0 || row >= height || Math.abs(row-previous) > 1) 
				throw new IllegalArgumentException("Invalid column");
			s[i] = row;
			previous = row;
		}

		// Update color
		for (int j = 0; j < width; j++) 
			for (int i = s[j]; i < height-1; i++)
				color[i][j] = color[i+1][j];

		// Update energy
		for (int j = 0; j < width; j++) {
			int row = s[j];

			// Update previous row
			if (row >= 1)
				energy[row-1][j] = delta(j, row-1, width, height-1);

			// Update next row (which should now be in current row position)
			if (row+1 < height)
				energy[row][j] = delta(j, row, width, height-1);

			// Shift the remaining rows
			for (int i = row+1; i < height-1; i++) 
				energy[i][j] = energy[i+1][j];
		}

		// Update height (width stays the same)
		height--;

	}


	// remove vertical seam from current picture
	public void removeVerticalSeam(int[] seam) {
		if (width <= 1) throw new IllegalArgumentException("Cannot remove vertical seams from picture\n");
		if (seam == null) throw new IllegalArgumentException("Null argument\n");
		if (seam.length != height) throw new IllegalArgumentException("Invalid seam length\n");
		
		int[] s = new int[seam.length];
		int previous = seam[0];
		for (int i = 0; i < seam.length; i++) {
			int col = seam[i];
			if (col < 0 || col >= width || Math.abs(col-previous) > 1) 
				throw new IllegalArgumentException("Invalid column\n");
			s[i] = col;
			previous = col;
		}

		// Update color
		for (int i = 0; i < height; i++) 
			for (int j = s[i]; j < width-1; j++) 
				color[i][j] = color[i][j+1];
		
		// Update energy
		for (int i = 0; i < height; i++) {
			int col = s[i];

			// Update previous column
			if (col >= 1) 
				energy[i][col-1] = delta(col-1, i, width-1, height);

			// Update next column (which should now be in current col position)
			if (col+1 < width)
				energy[i][col] = delta(col, i, width-1, height);

			// Shift the remaining columns
			for (int j = col+1; j < width-1; j++) 
				energy[i][j] = energy[i][j+1];
		}

		// Update width (height stays the same)
		width--;
	}
}