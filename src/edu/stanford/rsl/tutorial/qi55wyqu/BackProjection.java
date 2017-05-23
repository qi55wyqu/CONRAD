package edu.stanford.rsl.tutorial.qi55wyqu;

import java.util.ArrayList;

import edu.stanford.rsl.conrad.data.numeric.Grid1DComplex;
import edu.stanford.rsl.conrad.data.numeric.Grid2D;
import edu.stanford.rsl.conrad.data.numeric.Grid2DComplex;
import edu.stanford.rsl.conrad.data.numeric.InterpolationOperators;
import edu.stanford.rsl.conrad.utils.VisualizationUtil;
import edu.stanford.rsl.tutorial.parallel.ParallelBackprojector2D;
import ij.ImageJ;
import ij.ImagePlus;

public class BackProjection {
	
	public static Grid2D backProjection(Grid2D sinogram, int[] size, double[] spacing) {
		
		Grid2D backProjection = new Grid2D(size[0], size[1]);
		backProjection.setSpacing(spacing);
		backProjection.setOrigin(new double[] {-(size[0] - 1) * spacing[0] / 2, -(size[1] - 1) * spacing[1] / 2});
		
		for (int y = 0; y < backProjection.getHeight(); y++) {
			if (y % 50 == 0) {
				System.out.println("Backprojection running... line " + y);
			}
			for (int x = 0; x < backProjection.getWidth(); x++) {
				double[] backProjIdx = backProjection.indexToPhysical(x, y);
				
				for (int thetaIndex = 0; thetaIndex < sinogram.getHeight(); thetaIndex++) {
					
					double theta = thetaIndex * sinogram.getSpacing()[1];
					double s = backProjIdx[0] * Math.cos(Math.toRadians(theta)) + backProjIdx[1] * Math.sin(Math.toRadians(theta));
					
					double[] sinoIdx = sinogram.physicalToIndex(s, theta);
					float sinoVal = InterpolationOperators.interpolateLinear(sinogram, sinoIdx[0], sinoIdx[1]);
					backProjection.addAtIndex(x, y, sinoVal);
					
				}
			}
		}
		
		return backProjection;
	}
	
	public static Grid2D rampFilter(Grid2D sinogram) {
		Grid2DComplex sinoFourier = new Grid2DComplex(sinogram);
//		sinoFourier.show();
		sinoFourier.transformForward();
		Grid1DComplex rampFilter = new Grid1DComplex(sinoFourier.getWidth());
		rampFilter.setSpacing(1 / (sinogram.getSpacing()[0] * (sinoFourier.getWidth()-sinogram.getWidth())));
		for (int x = 0; x < rampFilter.getSize()[0] / 2; x++) {
			rampFilter.setAtIndex(x, x);
			rampFilter.setAtIndex(rampFilter.getSize()[0]-1-x, x);
		}
		rampFilter.show();
		for (int y = 0; y < sinoFourier.getHeight(); y++) {
			for (int x = 0; x < sinoFourier.getWidth(); x++) {
				sinoFourier.multiplyAtIndex(x, y, rampFilter.getRealAtIndex(x), rampFilter.getImagAtIndex(x));
			}
		}
//		sinoFourier.show();
		sinoFourier.transformInverse();
//		sinoFourier.show();
		Grid2D sinoFiltered = sinoFourier.getRealSubGrid(0, 0, sinogram.getWidth(), sinogram.getHeight());
		sinoFiltered.setSpacing(sinogram.getSpacing());
		sinoFiltered.setOrigin(sinogram.getOrigin());
//		sinoFiltered.show();
		return sinoFiltered;
	}

	public static void main(String[] args) {

		int[] size = new int[] { 512, 512 };
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
		
		Grid2D sinoFiltered = rampFilter(sinogram);
		ImagePlus sinoFil = VisualizationUtil.showGrid2D(sinoFiltered, "Ramp Filtered");
		sinoFil.show();
		
		Grid2D backProjection = backProjection(sinogram, size, spacing);
		ImagePlus backProj = VisualizationUtil.showGrid2D(backProjection, "Unfiltered Backprojection");
		backProj.show();
		
		Grid2D filteredBackProjection = backProjection(sinoFiltered, size, spacing);
		ImagePlus filteredBackProj = VisualizationUtil.showGrid2D(filteredBackProjection, "Filtered Backprojection");
		filteredBackProj.show();


//		ParallelBackprojector2D backProjection2 = new ParallelBackprojector2D(256, 256, 1, 1);
//		backProjection2.backprojectPixelDriven(sinogram).show("The Reconstruction");
		

	}

}
