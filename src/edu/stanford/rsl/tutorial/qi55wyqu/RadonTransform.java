package edu.stanford.rsl.tutorial.qi55wyqu;

import java.util.ArrayList;

import edu.stanford.rsl.apps.gui.opengl.PointCloudViewer;
import edu.stanford.rsl.conrad.data.numeric.Grid2D;
import edu.stanford.rsl.conrad.geometry.shapes.simple.Box;
import edu.stanford.rsl.conrad.geometry.shapes.simple.Point2D;
import edu.stanford.rsl.conrad.geometry.shapes.simple.PointND;

public class RadonTransform {

	public static float[][] radonTransform(
		Grid2D phantom, 
		int numProjections, 
		int numDetectorPixels, 
		double[] detectorSpacing, 
		double angleStart) {
		
			Box box = new Box(phantom.getWidth(), phantom.getHeight(), 0);
			double[] lowerPhantomCornerPhysical = phantom.indexToPhysical(0, phantom.getHeight());
			double[] upperPhantomCornerPhysical = phantom.indexToPhysical(phantom.getWidth(), 0);
			Point2D lowerCorner = new Point2D(lowerPhantomCornerPhysical);
			Point2D upperCorner = new Point2D(upperPhantomCornerPhysical);
			box.setLowerCorner(lowerCorner);
			box.setUpperCorner(upperCorner);
		
			float[][] radon = new float[numProjections][numDetectorPixels];
			final double rotationDegrees = 180;
			for (double theta = angleStart; theta <= angleStart + rotationDegrees; theta += rotationDegrees / numProjections) {
				double bottomX = 0.5 * phantom.getSpacing()[1] * Math.sin(theta);
				
				for (int s = 0; s < numDetectorPixels; s++) {
					
				}
			}
			
			return radon;
	}
	

	public static void main(String[] args) {
		
		
		
	}

}
