package edu.stanford.rsl.tutorial.qi55wyqu;

import edu.stanford.rsl.conrad.data.numeric.Grid2D;
import edu.stanford.rsl.conrad.data.numeric.opencl.OpenCLGrid2D;

public class OpenCLPhantom extends OpenCLGrid2D {
	
	public OpenCLPhantom(int[] size, double[] spacing) {
		super((Grid2D) (new Phantom(size, spacing)));
	}

}
