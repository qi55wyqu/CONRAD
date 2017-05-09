package edu.stanford.rsl.tutorial.qi55wyqu;

import java.util.ArrayList;

import edu.stanford.rsl.apps.gui.opengl.PointCloudViewer;
import edu.stanford.rsl.conrad.data.numeric.Grid2D;
import edu.stanford.rsl.conrad.data.numeric.InterpolationOperators;
import edu.stanford.rsl.conrad.geometry.shapes.simple.Box;
import edu.stanford.rsl.conrad.geometry.shapes.simple.Point2D;
import edu.stanford.rsl.conrad.geometry.shapes.simple.PointND;
import edu.stanford.rsl.conrad.geometry.shapes.simple.StraightLine;
import edu.stanford.rsl.conrad.numerics.SimpleVector;

public class RadonTransform {

	public static Grid2D radonTransform(
		Grid2D phantom, 
		int numProjections, 
		int numDetectorPixels, 
		double detectorSpacing) {
		
			double rotationDegrees = 180;
			double angularIncrement = rotationDegrees / numProjections;
		
			Grid2D sinogram = new Grid2D(numDetectorPixels, numProjections);
			sinogram.setOrigin(new double[] {-0.5*numDetectorPixels*detectorSpacing, 0.0});
			sinogram.setSpacing(new double[] {detectorSpacing, angularIncrement});
		
			Box box = new Box(phantom.getWidth(), phantom.getHeight(), 1);
			double[] lowerPhantomCornerPhysical = phantom.indexToPhysical(-0.5*phantom.getWidth(), -0.5*phantom.getHeight()); //spacing
			double[] upperPhantomCornerPhysical = phantom.indexToPhysical(0.5*phantom.getWidth(), 0.5*phantom.getHeight());
			PointND lowerCorner = new PointND(lowerPhantomCornerPhysical[0], lowerPhantomCornerPhysical[1], 0);
			PointND upperCorner = new PointND(upperPhantomCornerPhysical[0], upperPhantomCornerPhysical[1], 0);
			box.setLowerCorner(lowerCorner);
			box.setUpperCorner(upperCorner);
		
			for (double theta = 0; theta < rotationDegrees; theta += angularIncrement) {
				for (double s = -0.5*numDetectorPixels*detectorSpacing; s < 0.5*numDetectorPixels; s += detectorSpacing) {
					double x = s * Math.cos(Math.toRadians(theta));
					double y = s * Math.sin(Math.toRadians(theta));
					PointND currDetPix = new PointND(new double[] {x, y, 0});
					SimpleVector vec = new SimpleVector(new double[] {-s * Math.sin(Math.toRadians(theta)), s * Math.cos(Math.toRadians(theta)), 0});
					StraightLine line = new StraightLine(currDetPix, vec);
					ArrayList<PointND> intersections = box.intersect(line);
					double[] firstPoint = intersections.get(0).getCoordinates();
					double[] secondPoint = intersections.get(1).getCoordinates();
					float sum = 0f;
					for (double i = 0; i < phantom.getHeight(); i += 0.5) {
						double xWorld = firstPoint[0] + i*Math.sin(Math.toRadians(theta));
						double yWorld = firstPoint[1] + i*Math.cos(Math.toRadians(theta));
						double[] idx = phantom.physicalToIndex(xWorld, yWorld);
						sum += InterpolationOperators.interpolateLinear(phantom, idx[0], idx[1]);
					}
				}
			}
			
			return sinogram;
	}
	

	public static void main(String[] args) {
		
		Phantom phantom = new Phantom(new int[] {512, 512}, new double[] {1.0, 1.0});
		Grid2D sinogram = radonTransform(phantom, 360, 512, 1.0);
		
	}

}
