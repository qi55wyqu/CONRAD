package edu.stanford.rsl.tutorial.qi55wyqu;

import edu.stanford.rsl.conrad.data.numeric.Grid2D;

public class Phantom extends Grid2D {

	public Phantom(int[] size, double[] spacing, double[] origin) {
		super(size[0], size[1]);
		this.setSpacing(spacing);
		this.setOrigin(origin);
		this.drawEllipse(new int[] {(int)(0.4*size[0]), (int)(0.5*size[1])}, new int[] {(int)(0.2*size[0]), (int)(0.2*size[1])}, 0.3f);
		this.drawRectangle(new int[] {(int)(0.25*size[0]), (int)(0.25*size[1])}, new int[] {(int)(0.5*size[0]), (int)(0.4*size[1])}, 0.7f);
		this.drawEllipse(new int[] {(int)(0.6*size[0]), (int)(0.55*size[1])}, new int[] {(int)(0.15*size[0]), (int)(0.3*size[1])}, 0.5f);
	}
		
	public Phantom(int[] size) {
		this(size, new double[] {1.0, 1.0}, new double[] {0.0, 0.0});
	}
	
	public Phantom() {
		this(new int[] {512, 512});
	}
	
	private void drawRectangle(int[] pos, int[] size, float greyVal) {
		assert pos[0] + size[0] <= this.getWidth() && pos[1] + size[1] <= this.getHeight() 
			: "The given Rectangle is too large to fit into the image!";
		for (int i = pos[0]; i <= size[0]; i++) {
			for (int j = pos[1]; j <= size[1]; j++) {
				this.setAtIndex(i, j, greyVal);
			}
		}
	}
	
	private void drawEllipse(int[] center, int[] radii, float greyVal) {
		assert center[0] - radii[0] >= 0 
			&& center[0] + radii[0] <= this.getWidth() 
			&& center[1] - radii[1] >= 0
			&& center[1] + radii[1] <= this.getHeight()
			: "The given ellipse is too large to fit into the image!";
		/* Approximate ellipse using rectangles
		   Calculate edge points for one quarter of the ellipse and 
		   mirror it to the other quadrants
		*/ 
		for (int alpha = 0; alpha <= 90; alpha++) {
			// x = r_x * cos(a)
			// y = r_y * sin(a)
			int x = (int) (radii[0] * Math.cos(Math.toRadians(alpha)));
			int y = (int) (radii[1] * Math.sin(Math.toRadians(alpha)));
			for (int i = -x; i <= x; i++) {
				for (int j = -y; j <= y; j++) {
					this.setAtIndex(center[0]+i, center[1]+j, greyVal);
				}
			}
		}
	}
	
	public static void main(String[] args) {
		Phantom testPhantom = new Phantom(
			new int[] {512, 512}, new double[] {1.0, 1.0}, new double[] {0.0, 0.0}
		);
		testPhantom.show();
	}

}
