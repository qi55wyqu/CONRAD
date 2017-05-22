package edu.stanford.rsl.tutorial.qi55wyqu;

import ij.ImageJ;
import ij.ImagePlus;

import java.util.ArrayList;

import edu.stanford.rsl.conrad.data.numeric.Grid2D;
import edu.stanford.rsl.conrad.data.numeric.InterpolationOperators;
import edu.stanford.rsl.conrad.geometry.shapes.simple.Box;
import edu.stanford.rsl.conrad.geometry.shapes.simple.PointND;
import edu.stanford.rsl.conrad.geometry.shapes.simple.StraightLine;
import edu.stanford.rsl.conrad.numerics.SimpleVector;
import edu.stanford.rsl.conrad.utils.VisualizationUtil;
import edu.stanford.rsl.tutorial.parallel.ParallelProjector2D;

public class RadonTransform {

	public static Grid2D radonTransform(Grid2D phantom, int numProjections, int numDetectorPixels, double detectorSpacing) {

		final double rotationDegrees = 180;
		final double sampleSpacing = 0.5;
		final double angularIncrement = rotationDegrees / numProjections;
		final double detectorLength = numDetectorPixels * detectorSpacing;

		Grid2D sinogram = new Grid2D(numDetectorPixels, numProjections);
		sinogram.setOrigin(new double[] { -0.5 * detectorLength, 0.0 });
		sinogram.setSpacing(new double[] { detectorSpacing, angularIncrement });

		Box box = new Box();
		box.setLowerCorner(new PointND(-(phantom.getWidth()-1) * phantom.getSpacing()[0] / 2, -(phantom.getHeight()-1) * phantom.getSpacing()[1] / 2, 0));
		box.setUpperCorner(new PointND((phantom.getWidth()-1) * phantom.getSpacing()[0] / 2, (phantom.getHeight()-1) * phantom.getSpacing()[1] / 2, 0));

		for (int thetaIndex = 0; thetaIndex < rotationDegrees; thetaIndex++) {
			
			double theta = thetaIndex * angularIncrement;
			double sinTheta = Math.sin(Math.toRadians(theta));
			double cosTheta = Math.cos(Math.toRadians(theta));
			System.out.println("\u03B8 = " + theta + "°");

			for (int sIndex = 0; sIndex <= numDetectorPixels; sIndex++) {

				double s = -detectorLength / 2 + sIndex * detectorSpacing;

				PointND currDetPix = new PointND(new double[] { s * cosTheta, s * sinTheta, 0 });
				SimpleVector dir = new SimpleVector(-sinTheta, cosTheta, 0);
				StraightLine line = new StraightLine(currDetPix, dir);

				ArrayList<PointND> intersections = box.intersect(line);
				if (intersections.size() < 2) continue;

				double distance = intersections.get(0).euclideanDistance(intersections.get(1));
				SimpleVector sampleIncrement = new SimpleVector(intersections.get(1).getAbstractVector());
				sampleIncrement.subtract(intersections.get(0).getAbstractVector());
				sampleIncrement.divideBy(distance);
				sampleIncrement.multiplyBy(sampleSpacing);
				float sum = 0f;
				PointND currentPoint = new PointND(intersections.get(0));
				for (double i = 0; i <= distance; i += sampleSpacing) {
			
					double[] currentIdx = phantom.physicalToIndex(currentPoint.get(0), currentPoint.get(1));
					sum += InterpolationOperators.interpolateLinear(phantom, currentIdx[0], currentIdx[1]);
					currentPoint.getAbstractVector().add(sampleIncrement);
					
				}
				sum *= sampleSpacing;
				sinogram.setAtIndex(sIndex, thetaIndex, sum);
			}
		}
		return sinogram;
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
		
		Grid2D sinogram = radonTransform(phantom, numProjections, numDetectorPixels, detectorSpacing);
		ImagePlus sino1 = VisualizationUtil.showGrid2D(sinogram, "Meins");
		sino1.show();

		ParallelProjector2D projector = new ParallelProjector2D(Math.toRadians(180), Math.toRadians(180/numProjections), numDetectorPixels*		detectorSpacing, detectorSpacing);
		Grid2D sinogram2 = projector.projectRayDriven(phantom);
		ImagePlus sino2 = VisualizationUtil.showGrid2D(sinogram2, "Offiziell");
		sino2.show();
	}

}