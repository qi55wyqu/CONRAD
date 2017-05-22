package edu.stanford.rsl.tutorial.qi55wyqu;

import java.util.ArrayList;

import edu.stanford.rsl.conrad.data.numeric.Grid2D;
import edu.stanford.rsl.conrad.data.numeric.InterpolationOperators;
import edu.stanford.rsl.conrad.geometry.shapes.simple.Box;
import edu.stanford.rsl.conrad.geometry.shapes.simple.PointND;
import edu.stanford.rsl.conrad.geometry.shapes.simple.StraightLine;
import edu.stanford.rsl.conrad.numerics.SimpleVector;
import edu.stanford.rsl.conrad.utils.VisualizationUtil;
import edu.stanford.rsl.tutorial.parallel.ParallelBackprojector2D;
import ij.ImageJ;
import ij.ImagePlus;

public class BackProjection {
	
	public static Grid2D backProjection(Grid2D sinogram, int[] size, double[] spacing) {
		Grid2D backProjection = new Grid2D(size[0], size[1]);
		backProjection.setSpacing(spacing);
		backProjection.setOrigin(new double[] {-(size[0] - 1) * spacing[0] / 2, -(size[1] - 1) * spacing[1] / 2});
		final double angularSpacing = sinogram.getSpacing()[1];
		final int numDetectorPixels = size[0];
		final double detectorSpacing = spacing[0];
		final double detectorLength = numDetectorPixels * detectorSpacing;
		final double sampleSpacing = 1.0;

		Box box = new Box();
		box.setLowerCorner(new PointND(-(size[0]-1)*spacing[0] / 2, -(size[1]-1)*spacing[1] / 2, 0));
		box.setUpperCorner(new PointND((size[0]-1)*spacing[0] / 2, (size[1]-1)*spacing[1] / 2, 0));

		for (int thetaIndex = 0; thetaIndex < sinogram.getHeight(); thetaIndex++) {
			
			double theta = thetaIndex * angularSpacing;
			double sinTheta = Math.sin(Math.toRadians(theta));
			double cosTheta = Math.cos(Math.toRadians(theta));
			System.out.println("\u03B8 = " + theta + "°");

			for (int sIndex = 0; sIndex < sinogram.getWidth(); sIndex++) {

				double s = -detectorLength / 2 + sIndex * detectorSpacing;
				float currSinoVal = sinogram.getAtIndex(sIndex, thetaIndex);

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
				PointND currentPoint = new PointND(intersections.get(0));
				for (double i = 1; i < distance; i += sampleSpacing) {

					currentPoint.getAbstractVector().add(sampleIncrement);
					double[] currentIdx = backProjection.physicalToIndex(currentPoint.get(0), currentPoint.get(1));
//					if (currentIdx[0] < 1 || currentIdx[0] > backProjection.getWidth()-1 || currentIdx[1] < 1 || currentIdx[1] < backProjection.getHeight()-1) {
//						currentPoint.getAbstractVector().add(sampleIncrement);
//						continue;
//					}
					InterpolationOperators.addInterpolateLinear(backProjection, currentIdx[0], currentIdx[1], currSinoVal);
					
				}
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
		ImagePlus sino = VisualizationUtil.showGrid2D(sinogram, "Radon Transform");
		sino.show();
		
		Grid2D backProjection = backProjection(sinogram, size, spacing);
		ImagePlus backProj = VisualizationUtil.showGrid2D(backProjection, "Unfiltered Backprojection");
		backProj.show();

//		ParallelBackprojector2D backProjection2 = new ParallelBackprojector2D(256, 256, 1, 1);
//		backProjection2.backprojectPixelDriven(sinogram).show("The Reconstruction");
		

	}

}
