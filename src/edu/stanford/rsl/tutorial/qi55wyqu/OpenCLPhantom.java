package edu.stanford.rsl.tutorial.qi55wyqu;

import edu.stanford.rsl.conrad.data.numeric.Grid2D;
import edu.stanford.rsl.conrad.data.numeric.opencl.OpenCLGrid2D;

public class OpenCLPhantom extends OpenCLGrid2D {
	
	public OpenCLPhantom(int[] size, double[] spacing) {
		super((Grid2D) new Phantom(size, spacing));
//		Grid2D phantom2D = new Phantom(size, spacing);
//		Grid1D phantom1D = new Grid1D(size[0]*size[1]);
//		for (int y = 0; y < phantom2D.getHeight(); y++) {
//			for (int x = 0; x < phantom2D.getWidth(); x++) {
//				phantom1D.setAtIndex(y * phantom2D.getWidth() + x, phantom2D.getAtIndex(x, y));
//			}
//		}
//		super(phantom1D);
	}
	
}
