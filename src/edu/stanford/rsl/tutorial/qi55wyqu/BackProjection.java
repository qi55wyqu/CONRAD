package edu.stanford.rsl.tutorial.qi55wyqu;

import edu.stanford.rsl.conrad.data.numeric.Grid1D;
import edu.stanford.rsl.conrad.data.numeric.Grid2D;
import edu.stanford.rsl.conrad.utils.VisualizationUtil;
import ij.ImageJ;
import ij.ImagePlus;

public class BackProjection {
	
	public static Grid2D backProjection(Grid2D sinogram, int[] size, double[] spacing) {
		Grid2D backProjection = new Grid2D(size[0], size[1]);
		backProjection.setSpacing(spacing);
		final double angularSpacing = sinogram.getSpacing()[1];
		for (int thetaIndex = 0; thetaIndex < sinogram.getHeight(); thetaIndex++) {
			double theta = thetaIndex * angularSpacing;
			double sinTheta = Math.sin(Math.toRadians(theta));
			double cosTheta = Math.cos(Math.toRadians(theta));
			Grid1D projectionLine = sinogram.getSubGrid(thetaIndex);
			for (int s = 0; s < backProjection.getWidth(); s++) {
				int x = (int) (s * cosTheta);
				int y = (int) (s * sinTheta);
				backProjection.setAtIndex(x, y, projectionLine.getAtIndex(s));
			}
		}
		return backProjection;
	}

	public static void main(String[] args) {

		int[] size = new int[] { 256, 256 };
		double[] spacing = new double[] {1.0, 1.0};
		int numProjections = 180;
		int numDetectorPixels = size[0];
		double detectorSpacing = 1.0;
		
		new ImageJ();
		
		Phantom phantom = new Phantom(size, spacing);
		ImagePlus phan = VisualizationUtil.showGrid2D(phantom, "Test Phantom");
		phan.show();
		
		Grid2D sinogram = RadonTransform.radonTransform(phantom, numProjections, numDetectorPixels, detectorSpacing);
		ImagePlus sino = VisualizationUtil.showGrid2D(sinogram, "Meins");
		sino.show();
		
//		ParallelBackprojector2D backproj = new ParallelBackprojector2D(256, 256, 1, 1);
//		backproj.backprojectPixelDriven(sinogram).show("The Reconstruction");
		
		Grid2D backProjection = backProjection(sinogram, size, spacing);
		backProjection.show();
		

	}

}
