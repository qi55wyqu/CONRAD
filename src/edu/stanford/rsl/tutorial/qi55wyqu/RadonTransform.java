package edu.stanford.rsl.tutorial.qi55wyqu;

import ij.ImageJ;
import ij.ImagePlus;

import java.util.ArrayList;

import edu.stanford.rsl.apps.gui.opengl.PointCloudViewer;
import edu.stanford.rsl.conrad.data.numeric.Grid2D;
import edu.stanford.rsl.conrad.data.numeric.InterpolationOperators;
import edu.stanford.rsl.conrad.geometry.shapes.simple.Box;
import edu.stanford.rsl.conrad.geometry.shapes.simple.Point2D;
import edu.stanford.rsl.conrad.geometry.shapes.simple.PointND;
import edu.stanford.rsl.conrad.geometry.shapes.simple.StraightLine;
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
			final double[] phantomPhysicalSize = new double[2];
			phantomPhysicalSize[0] = phantom.getWidth() * phantomSpacing[0];
			phantomPhysicalSize[1] = phantom.getHeight() * phantomSpacing[1];
					
			Grid2D sinogram = new Grid2D(numDetectorPixels, numProjections);
			sinogram.setOrigin(new double[] {-0.5*detectorLength, 0.0});
			sinogram.setSpacing(new double[] {detectorSpacing, angularIncrement});
			
			int method = 1;
		
			if (method == 1) {
				Box box = new Box(phantomPhysicalSize[0], phantomPhysicalSize[1], 1);
				double[] lowerPhantomCornerPhysical = phantom.indexToPhysical(
					-0.5*phantomPhysicalSize[0], -0.5*phantomPhysicalSize[1]);
				double[] upperPhantomCornerPhysical = phantom.indexToPhysical(
					0.5*phantomPhysicalSize[0], 0.5*phantomPhysicalSize[1]);
				box.setLowerCorner(new PointND(lowerPhantomCornerPhysical[0], lowerPhantomCornerPhysical[1], 0));
				box.setUpperCorner(new PointND(upperPhantomCornerPhysical[0], upperPhantomCornerPhysical[1], 0));
				
				for (double theta = 0; theta < rotationDegrees; theta += angularIncrement) {
					System.out.println("Theta = " + theta);
					for (double s = -0.5*detectorLength; s <= 0.5*detectorLength; s += detectorSpacing) {
						double x = s * Math.cos(Math.toRadians(theta));
						double y = s * Math.sin(Math.toRadians(theta));
						PointND currDetPix = new PointND(new double[] {x, y, 0});
						SimpleVector vec = new SimpleVector(new double[] {x-(s/Math.cos(Math.toRadians(theta))), -y, 0});
						StraightLine line = new StraightLine(currDetPix, vec);
						ArrayList<PointND> intersections = box.intersect(line);
						if (s == -0.5*detectorLength) {
							System.out.println(line.toString());
							if (intersections.size() < 2) {
								System.out.println("error");
								continue;
							}
						}
						if (intersections.size() < 2) {
							continue;
						}
						double[] firstPoint = intersections.get(0).getCoordinates();
						double[] secondPoint = intersections.get(1).getCoordinates();
						
						float sum = 0f;
						
						double xWorld = firstPoint[0];
						double yWorld = firstPoint[1];
						double step = 0;
						
						while (xWorld <= secondPoint[0] && yWorld <= secondPoint[1]) {
							xWorld = firstPoint[0] + step * Math.sin(Math.toRadians(theta));
							yWorld = firstPoint[1] + step * Math.cos(Math.toRadians(theta));
							double[] idx = phantom.physicalToIndex(xWorld, yWorld);
							sum += InterpolationOperators.interpolateLinear(phantom, idx[0], idx[1]);
							step += sampleSpacing;
						}
						double[] pointSinogram = sinogram.physicalToIndex(s, theta);
						sinogram.putPixelValue((int) pointSinogram[0], (int) pointSinogram[1], sum);
					}
				}
			} else {
				for (double theta = 0; theta < rotationDegrees; theta += angularIncrement) {
					System.out.println("Theta = " + theta);
					for (double s = -0.5*detectorLength; s < 0.5*detectorLength; s += detectorSpacing) {
						float pixVal = 0;
						for (int x = 0; x < phantom.getWidth(); x++) {
							for (int y = 0; y < phantom.getHeight(); y++) {
								double[] idx = phantom.indexToPhysical(x, y);
								int condition = (int) (idx[0]*Math.cos(Math.toRadians(theta)) + idx[1]*Math.sin(Math.toRadians(theta)) - s);
								if (condition == 0) {
									pixVal += phantom.getAtIndex(x, y);
								}
							}
						}
						double[] pos = sinogram.physicalToIndex(s, theta);
						sinogram.addAtIndex((int) pos[0], (int) pos[1], pixVal);
					}
				}
			}
			
			return sinogram;
	}
	

	public static void main(String[] args) {
		
		Phantom phantom = new Phantom(new int[] {128, 128}, new double[] {1.0, 1.0});
		Grid2D sinogram = radonTransform(phantom, 180, 128, 1.0);
		new ImageJ();
		ImagePlus img = VisualizationUtil.showGrid2D(sinogram, "Test Phantom");
		img.show();
	}

}
