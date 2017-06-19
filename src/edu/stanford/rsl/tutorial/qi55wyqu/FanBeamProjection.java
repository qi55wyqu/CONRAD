package edu.stanford.rsl.tutorial.qi55wyqu;

import edu.stanford.rsl.conrad.data.numeric.Grid2D;
import edu.stanford.rsl.conrad.geometry.shapes.simple.Point2D;
import edu.stanford.rsl.conrad.utils.VisualizationUtil;
import ij.ImageJ;
import ij.ImagePlus;

public class FanBeamProjection {
	
	protected int numProjections, numDetectorPixels;
	protected double detectorSpacing, angularIncrement, distSourceIso, distSourceDet;
	protected double sampleSpacing, detectorLength, gamma, rotationDegress, distDetIso;

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
		this.gamma = Math.atan(0.5 * this.detectorLength / this.distSourceDet);
		this.rotationDegress = this.numProjections * this.angularIncrement;
		this.distDetIso = this.distSourceDet - this.distSourceIso;
	}
	
	public Grid2D project(Grid2D inputImage) {

		double diag = Math.sqrt(
			  0.5 * inputImage.getWidth()  * inputImage.getSpacing()[0] 
			+ 0.5 * inputImage.getHeight() * inputImage.getSpacing()[1]
		);
		if (this.distSourceIso <= diag) {
			System.err.println("Tube will hit object!");
		}
		else if (this.distDetIso <= diag) {
			System.err.println("Detector will hit object!");
		}
		
		Grid2D fanogram = new Grid2D(this.numDetectorPixels, this.numProjections);
		fanogram.setSpacing(new double[] { this.detectorSpacing, this.angularIncrement });
		fanogram.setOrigin(new double[] { -(this.numDetectorPixels - 1) * this.detectorSpacing / 2, 0.0 });
		
		for (int betaIndex = 0; betaIndex < this.rotationDegress; betaIndex++) {
			
			double beta = betaIndex * this.angularIncrement;
			double sourcePosX = Math.cos(Math.toRadians(beta)) * this.distSourceIso;
			double sourcePosY = Math.sin(Math.toRadians(beta)) * this.distSourceIso;
			Point2D sourcePos = new Point2D(new double[] { sourcePosX, sourcePosY });
			
			for (int tIndex = 0; tIndex < this.numDetectorPixels; tIndex++) {
				
				
				
			}
		}
				
		return fanogram;
		
	}
	
	public static void main(String[] args) {

		int[] size = new int[] { 256, 256 };
		double[] spacing = new double[] {1.0, 1.0};
		int numProjections = 180;
		int numDetectorPixels = size[0];
		double detectorSpacing = 1.0;
		double angularIncrement = 1.0;
		double distSourceIso = 1.5 * size[0] * spacing[0];
		double distSourceDet = distSourceIso;
		
		new ImageJ();
		
		Phantom phantom = new Phantom(size, spacing);
		ImagePlus phan = VisualizationUtil.showGrid2D(phantom, "Test Phantom");
		phan.show();
		
		FanBeamProjection fanBeamProjection = new FanBeamProjection(detectorSpacing, numDetectorPixels, angularIncrement, numProjections, distSourceIso, distSourceDet);
		Grid2D fanogram = fanBeamProjection.project(phantom);
		ImagePlus fano1 = VisualizationUtil.showGrid2D(fanogram, "Fanogram");
		fano1.show();

	}

}