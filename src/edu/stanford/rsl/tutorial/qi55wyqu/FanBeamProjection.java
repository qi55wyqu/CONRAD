package edu.stanford.rsl.tutorial.qi55wyqu;

import edu.stanford.rsl.conrad.data.numeric.Grid2D;
import edu.stanford.rsl.conrad.utils.VisualizationUtil;
import ij.ImageJ;
import ij.ImagePlus;

public class FanBeamProjection {
	
	protected int numProjections, numDetectorPixels;
	protected double detectorSpacing, numDetectorElements, sampleSpacing, angularIncrement, distSourceIso, distSourceDet;

	public FanBeamProjection(
		double detectorSpacing, 
		int numDetectorPixels, 
		double angularIncrement, 
		int numProjections, 
		double distSourceIso, 
		double distSourceDet) {
			
		this.detectorSpacing = detectorSpacing;
		this.numDetectorPixels = numDetectorPixels;
		this.angularIncrement = angularIncrement;
		this.numProjections = numProjections;
		this.distSourceIso = distSourceIso;
		this.distSourceDet = distSourceDet;
		this.sampleSpacing = 0.5;
	}
	
	public Grid2D project(Grid2D inputImage) {
		
		Grid2D fanogram = new Grid2D(this.numDetectorPixels, this.numProjections);
		fanogram.setSpacing(new double[] { this.detectorSpacing, this.angularIncrement });
		fanogram.setOrigin(new double[] { -(this.numDetectorPixels - 1) * this.detectorSpacing / 2, 0.0 });
		
		
		
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