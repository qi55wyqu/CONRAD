package edu.stanford.rsl.tutorial.qi55wyqu;

import java.util.ArrayList;

import edu.stanford.rsl.conrad.data.numeric.Grid2D;
import edu.stanford.rsl.conrad.data.numeric.InterpolationOperators;
import edu.stanford.rsl.conrad.geometry.shapes.simple.Box;
import edu.stanford.rsl.conrad.geometry.shapes.simple.PointND;
import edu.stanford.rsl.conrad.geometry.shapes.simple.StraightLine;
import edu.stanford.rsl.conrad.numerics.SimpleVector;
import edu.stanford.rsl.conrad.utils.VisualizationUtil;
import ij.ImageJ;
import ij.ImagePlus;

public class FanBeamProjection {
	
	protected int numProjections, numDetectorPixels;
	protected double detectorSpacing, angularIncrement, distSourceIso, distSourceDet;
	private double sampleSpacing, detectorLength, rotationDegress, distDetIso, delta;
	private double[] beta;
	private PointND[] sourcePos;
	private SimpleVector[] dirSourceDet;
	private SimpleVector[] dirDet;

	public FanBeamProjection(
		double detectorSpacing, 
		int numDetectorPixels, 
		double angularIncrement, 
		int numProjections, 
		double distSourceIso, 
		double distSourceDet) {
			
		this.detectorSpacing   = detectorSpacing;
		this.numDetectorPixels = numDetectorPixels;
		this.angularIncrement  = angularIncrement;
		this.numProjections    = numProjections;
		this.distSourceIso     = distSourceIso;
		this.distSourceDet     = distSourceDet;
		
		this.initialise();
	}
	
	private void initialise() {
		this.sampleSpacing = 0.5;
		this.detectorLength = this.numDetectorPixels * this.detectorSpacing;
		this.rotationDegress = this.numProjections * this.angularIncrement;
		this.distDetIso = this.distSourceDet - this.distSourceIso;
		
		this.delta = Math.atan(0.5 * this.detectorLength / this.distSourceDet);
		
		this.beta = new double[this.numProjections];
		this.sourcePos = new PointND[this.numProjections];
		this.dirSourceDet = new SimpleVector[this.numProjections];
		this.dirDet = new SimpleVector[this.numProjections];
		for (int betaIndex = 0; betaIndex < this.numProjections; betaIndex++) {	
			double beta = betaIndex * this.angularIncrement;
			this.beta[betaIndex] = beta;
			double cosBeta = Math.cos(Math.toRadians(beta));
			double sinBeta = Math.sin(Math.toRadians(beta));
			this.sourcePos[betaIndex] = new PointND(new double[] { cosBeta * this.distSourceIso, sinBeta * this.distSourceIso, 0 });
			this.dirSourceDet[betaIndex] = new SimpleVector(cosBeta, sinBeta, 0);
			this.dirDet[betaIndex] = new SimpleVector(-sinBeta, cosBeta, 0);
		}
		
	}
	
	private boolean isValidGeometry(Grid2D inputImage) {
		double diag = Math.sqrt(
				Math.pow(0.5 * inputImage.getWidth()  * inputImage.getSpacing()[0], 2) 
			  + Math.pow(0.5 * inputImage.getHeight() * inputImage.getSpacing()[1], 2)
		);
		if (this.distSourceIso <= diag) {
			System.err.println("Tube will hit object!");
			return false;
		}
		else if (this.distDetIso <= diag) {
			System.err.println("Detector will hit object!");
			return false;
		}
		return true;
	}
	
	private Box getBoundingBox(Grid2D inputImage) {
		Box box = new Box();
		box.setLowerCorner(new PointND(-(inputImage.getWidth()-1) * inputImage.getSpacing()[0] / 2, -(inputImage.getHeight()-1) * inputImage.getSpacing()[1] / 2, 0));
		box.setUpperCorner(new PointND((inputImage.getWidth()-1) * inputImage.getSpacing()[0] / 2, (inputImage.getHeight()-1) * inputImage.getSpacing()[1] / 2, 0));
		return box;
	}
	
	public Grid2D project(Grid2D inputImage) {

		Grid2D fanogram = new Grid2D(this.numDetectorPixels, this.numProjections);
		if (!isValidGeometry(inputImage)) return fanogram;
		
		fanogram.setSpacing(new double[] { this.detectorSpacing, this.angularIncrement });
		fanogram.setOrigin(new double[] { -(this.numDetectorPixels - 1) * this.detectorSpacing / 2, 0.0 });
		
		Box box = getBoundingBox(inputImage);
		
		for (int betaIndex = 0; betaIndex < this.numProjections; betaIndex++) {			
			System.out.println("\u03B2" +  " = " + this.beta[betaIndex] + "°");
			
			for (int tIndex = 0; tIndex < this.numDetectorPixels; tIndex++) {

				SimpleVector currDetPixVec = new SimpleVector(this.sourcePos[betaIndex].getAbstractVector());
				currDetPixVec.subtract(this.dirSourceDet[betaIndex].multipliedBy(this.distSourceDet));
				currDetPixVec.subtract(this.dirDet[betaIndex].multipliedBy(fanogram.indexToPhysical(tIndex, 0)[0]));
				
				StraightLine line = new StraightLine(this.sourcePos[betaIndex], new PointND(currDetPixVec));
				ArrayList<PointND> intersections = box.intersect(line);
				if (intersections.size() < 2) continue;

				double distance = intersections.get(0).euclideanDistance(intersections.get(1));
				SimpleVector sampleIncrement = new SimpleVector(intersections.get(1).getAbstractVector());
				sampleIncrement.subtract(intersections.get(0).getAbstractVector());
				sampleIncrement.divideBy(distance);
				sampleIncrement.multiplyBy(this.sampleSpacing);

				float sum = 0f;
				PointND currentPoint = new PointND(intersections.get(0));
				for (double i = 0; i <= distance; i += this.sampleSpacing) {
					double[] currentIdx = inputImage.physicalToIndex(currentPoint.get(0), currentPoint.get(1));
					sum += InterpolationOperators.interpolateLinear(inputImage, currentIdx[0], currentIdx[1]);
					currentPoint.getAbstractVector().add(sampleIncrement);
				}
				fanogram.setAtIndex(tIndex, betaIndex, (float) (sum * this.sampleSpacing));
				
			}
		}
				
		return fanogram;
		
	}
	
	public Grid2D rebin(Grid2D fanogram) {
		Grid2D sinogram = new Grid2D(fanogram.getWidth(), fanogram.getHeight());
		sinogram.setSpacing(fanogram.getSpacing());
		sinogram.setOrigin(fanogram.getOrigin());
		for (int betaIndex = 0; betaIndex < fanogram.getHeight(); betaIndex++) {
			for (int tIndex = 0; tIndex < fanogram.getWidth(); tIndex++) {
				double fanoWorldT = fanogram.indexToPhysical(tIndex, betaIndex)[0];
				double gamma = Math.atan(fanoWorldT / this.distSourceDet);
				double beta = betaIndex * this.angularIncrement;
//				if (beta >= Math.PI + 2 * gamma && beta <= Math.PI + 2 * this.delta) {
//					beta = beta + 2 * gamma + Math.PI;
//					gamma = -gamma;
//				}
				double s = this.distSourceIso * Math.sin(gamma);
				double theta = beta + gamma;
				double[] sinoIdx = sinogram.physicalToIndex(s, theta);				
				float val = InterpolationOperators.interpolateLinear(fanogram, sinoIdx[0], sinoIdx[1]);
				sinogram.setAtIndex((int) (sinoIdx[0]), (int) (sinoIdx[1]), val);
			}
		}
		
		return sinogram;
	}
	
	public static void main(String[] args) {

		int[] size = new int[] { 256, 256 };
		double[] spacing = new double[] {1.0, 1.0};
		int numProjections = 180;
		int numDetectorPixels = 2 * size[0];
		double detectorSpacing = 1.0;
		double angularIncrement = 1.0;
		double distSourceIso = 10 * size[0];
		double distSourceDet = 20 * size[0];
		
		new ImageJ();
		
		Phantom phantom = new Phantom(size, spacing);
		ImagePlus phan = VisualizationUtil.showGrid2D(phantom, "Test Phantom");
		phan.show();
		
		FanBeamProjection fanBeamProjection = new FanBeamProjection(detectorSpacing, numDetectorPixels, angularIncrement, numProjections, distSourceIso, distSourceDet);
		Grid2D fanogram = fanBeamProjection.project(phantom);
		ImagePlus fano1 = VisualizationUtil.showGrid2D(fanogram, "Fanogram");
		fano1.show();
		
		Grid2D sinogram = fanBeamProjection.rebin(fanogram);
		ImagePlus sino1 = VisualizationUtil.showGrid2D(sinogram, "Rebinned Sinogram");
		sino1.show();
		
//		BackProjection backProjection = new BackProjection(size, spacing);
//		Grid2D backProjected = backProjection.backProject(sinogram);
//		ImagePlus backProj = VisualizationUtil.showGrid2D(backProjected, "Unfiltered Backprojection");
//		backProj.show();
//		
//		Grid2D sinoRamLakFiltered = BackProjection.ramLakFilter(sinogram);
//		Grid2D ramLakFilteredBackProjection = backProjection.backProject(sinoRamLakFiltered);
//		ImagePlus ramLakFilteredBackProj = VisualizationUtil.showGrid2D(ramLakFilteredBackProjection, "RamLak-filtered Backprojection");
//		ramLakFilteredBackProj.show();

	}

}