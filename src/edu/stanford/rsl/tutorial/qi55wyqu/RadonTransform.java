package edu.stanford.rsl.tutorial.qi55wyqu;

import ij.ImageJ;
import ij.ImagePlus;

import java.util.ArrayList;

import edu.stanford.rsl.conrad.data.numeric.Grid2D;
import edu.stanford.rsl.conrad.data.numeric.InterpolationOperators;
import edu.stanford.rsl.conrad.geometry.shapes.simple.Box;
import edu.stanford.rsl.conrad.geometry.shapes.simple.PointND;
import edu.stanford.rsl.conrad.geometry.shapes.simple.StraightLine;
import edu.stanford.rsl.conrad.geometry.transforms.Transform;
import edu.stanford.rsl.conrad.geometry.transforms.Translation;
import edu.stanford.rsl.conrad.numerics.SimpleVector;
import edu.stanford.rsl.conrad.utils.VisualizationUtil;

public class RadonTransform {

	public static Grid2D radonTransform(
		Grid2D phantom, 
		int numProjections, 
		int numDetectorPixels, 
		double detectorSpacing) {
		
			final double rotationDegrees = 180;
			final double sampleSpacing = 0.5;
			final double angularIncrement = rotationDegrees / numProjections;
			final double detectorLength = numDetectorPixels * detectorSpacing;
			final double[] phantomSpacing = phantom.getSpacing();
			final double[] phantomPhysicalSize = new double[] {phantom.getWidth() * phantomSpacing[0], phantom.getHeight() * phantomSpacing[1]};
					
			Grid2D sinogram = new Grid2D(numDetectorPixels, numProjections);
			sinogram.setOrigin(new double[] {-0.5*detectorLength, 0.0});
			sinogram.setSpacing(new double[] {detectorSpacing, angularIncrement});
			
			int method = 1;
		
			if (method == 1) {
				// Translation translation = new Translation(-0.5*phantomPhysicalSize[0], -0.5*phantomPhysicalSize[1], 0);
				Box box = new Box(phantomPhysicalSize[0], phantomPhysicalSize[1], 1);
				box.setLowerCorner(new PointND(-phantomPhysicalSize[0]/2, -phantomPhysicalSize[1]/2, 0));
				box.setUpperCorner(new PointND(phantomPhysicalSize[0]/2, phantomPhysicalSize[1]/2, 0));
				//box.applyTransform(translation);

				for (int thetaIndex = 0; thetaIndex < rotationDegrees; thetaIndex++) { 
					double theta = thetaIndex * angularIncrement;
					double sinTheta = Math.sin(Math.toRadians(theta));
					double cosTheta = Math.cos(Math.toRadians(theta));
					System.out.println("Theta = " + theta);

					for (int sIndex = 0; sIndex <= numDetectorPixels; sIndex++) {
						
						double s = -0.5 * detectorLength + sIndex * detectorSpacing;
						
						PointND currDetPix = new PointND(new double[] {s * cosTheta, s * sinTheta, 0});
						//PointND perpPoint = new PointND(new double[] {s * cosTheta - sinTheta, s * sinTheta + cosTheta, 0});
						//StraightLine line = new StraightLine(currDetPix, perpPoint);
						SimpleVector dir = new SimpleVector(-sinTheta, cosTheta, 0);
						StraightLine line = new StraightLine(currDetPix, dir);						
						
						
						/*
						 * tests
						 * */
						PointND currDetPix2 = new PointND(new double[] {0, 0, 0});
						StraightLine line2 = new StraightLine(currDetPix2, dir);
						ArrayList<PointND> intersections2 = box.intersect(line2);
					
						
						ArrayList<PointND> intersections = box.intersect(line);
						if (intersections.size() < 2){
							continue;
						}
						//System.out.print(intersections.get(0));
						
						//System.out.println(intersections.get(1));
						SimpleVector sampleIncrement = new SimpleVector(intersections.get(1).getAbstractVector());
						
						sampleIncrement.subtract(intersections.get(0).getAbstractVector());
						double distance = intersections.get(0).euclideanDistance(intersections.get(1));
						sampleIncrement.divideBy(distance / sampleSpacing);
						float sum = 0f;
						PointND currentPoint = new PointND(intersections.get(0));
						for (double i = 0; i < distance; i += sampleSpacing) {
							double[] currentIdx = phantom.physicalToIndex(currentPoint.get(0), currentPoint.get(1));
							sum += InterpolationOperators.interpolateLinear(phantom, currentIdx[0], currentIdx[1]);
							currentPoint.getAbstractVector().add(sampleIncrement);
						}
						sum *= sampleSpacing;
						sinogram.setAtIndex(sIndex, thetaIndex, sum);
						
						/*
						double[] firstPoint = intersections.get(0).getCoordinates();
						double[] secondPoint = intersections.get(1).getCoordinates();
						double xWorld = firstPoint[0];
						double yWorld = firstPoint[1];
						double step = 0d;
						float sum = 0f;
						while (xWorld <= secondPoint[0] && yWorld <= secondPoint[1]) {
							System.out.println("first: {"+firstPoint[0]+", "+firstPoint[1]+"}\tcurr: {"+xWorld+", "+yWorld+"}\n");
							xWorld = firstPoint[0] + step * cosTheta;
							yWorld = firstPoint[1] + step * sinTheta;
							
							//double[] idx = phantom.physicalToIndex(xWorld, yWorld);
							//sum += InterpolationOperators.interpolateLinear(phantom, idx[0], idx[1]);
							sum += InterpolationOperators.interpolateLinear(phantom, xWorld, yWorld);
							step += sampleSpacing;
						}
						double sIdx = s / detectorSpacing + 0.5 * detectorLength;
						sinogram.setAtIndex((int) sIdx, (int) thetaIdx, sum);
						*/
					}
				}
			} else {
				for (double theta = 0; theta < rotationDegrees; theta += angularIncrement) {
					System.out.println("Theta = " + theta);
					for (double s = -0.5*detectorLength; s < 0.5*detectorLength; s += detectorSpacing) {
						float sum = 0;
						for (int x = 0; x < phantom.getWidth(); x++) {
							for (int y = 0; y < phantom.getHeight(); y++) {
								double[] idx = phantom.indexToPhysical(x, y);
								int condition = (int) (idx[0]*Math.cos(Math.toRadians(theta)) + idx[1]*Math.sin(Math.toRadians(theta)) - s);
								if (condition == 0) {
									sum += phantom.getAtIndex(x, y);
								}
							}
						}
						double[] pos = sinogram.physicalToIndex(s, theta);
						sinogram.addAtIndex((int) pos[0], (int) pos[1], sum);
					}
				}
			}
			
			return sinogram;
	}
	

	public static void main(String[] args) {
		
		Phantom phantom = new Phantom(new int[] {256, 256}, new double[] {1.0, 1.5});
		//phantom.show();
		Grid2D sinogram = radonTransform(phantom, 180, 256, 1.0);
		//ParallelProjector2D projector = new ParallelProjector2D(Math.PI, Math.PI/180.0, 128, 0.5);
		//Grid2D sinogram = projector.projectRayDriven(phantom);
		new ImageJ();
		ImagePlus phan = VisualizationUtil.showGrid2D(phantom, "Test Phantom");
		phan.show();
		ImagePlus sino = VisualizationUtil.showGrid2D(sinogram, "Test Phantom");
		sino.show();
	}

}
