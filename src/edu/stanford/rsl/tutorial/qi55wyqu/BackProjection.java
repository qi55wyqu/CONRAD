package edu.stanford.rsl.tutorial.qi55wyqu;

import edu.stanford.rsl.conrad.data.numeric.Grid1D;
import edu.stanford.rsl.conrad.data.numeric.Grid1DComplex;
import edu.stanford.rsl.conrad.data.numeric.Grid2D;
import edu.stanford.rsl.conrad.data.numeric.Grid2DComplex;
import edu.stanford.rsl.conrad.data.numeric.InterpolationOperators;
import edu.stanford.rsl.conrad.utils.VisualizationUtil;
import edu.stanford.rsl.tutorial.parallel.ParallelBackprojector2D;
import ij.ImageJ;
import ij.ImagePlus;

public class BackProjection {
	
	protected int[] size;
	protected double[] spacing;
	
	public BackProjection(int[] size, double[] spacing) {
		this.size = size;
		this.spacing = spacing;
	}
	
	public Grid2D backProject(Grid2D sinogram) {
		
		Grid2D backProjection = new Grid2D(this.size[0], this.size[1]);
		backProjection.setSpacing(this.spacing);
		backProjection.setOrigin(new double[] {-(this.size[0] - 1) * this.spacing[0] / 2, -(this.size[1] - 1) * this.spacing[1] / 2});
		
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
		sinoFourier.transformForward();
		final int width = sinoFourier.getWidth();
		Grid1DComplex rampFilter = new Grid1DComplex(width);
		rampFilter.setSpacing(1 / (sinogram.getSpacing()[0] * (sinoFourier.getWidth() - sinogram.getWidth())));
//		for (int x = 0; x <= rampFilter.getSize()[0] / 2; x++) {
//			rampFilter.setAtIndex((int) (x/2), (float) (x/2));
//			rampFilter.setAtIndex((int) ((rampFilter.getSize()[0]-x)/2), (float) (x/2));
//		}
		for (int x = 0; x < sinogram.getWidth() / 2; x++) {
			rampFilter.setAtIndex(x, x);
		}
		for (int x = sinogram.getWidth() / 2; x < sinogram.getWidth(); x++) {
			rampFilter.setAtIndex(x, sinogram.getWidth() - x);
		}
		rampFilter.show();
		for (int y = 0; y < sinoFourier.getHeight(); y++) {
			for (int x = 0; x < sinoFourier.getWidth(); x++) {
				sinoFourier.multiplyAtIndex(x, y, rampFilter.getRealAtIndex(x), rampFilter.getImagAtIndex(x));
			}
		}
		sinoFourier.transformInverse();
		Grid2D sinoFiltered = sinoFourier.getRealSubGrid(0, 0, sinogram.getWidth(), sinogram.getHeight());
		sinoFiltered.setSpacing(sinogram.getSpacing());
		sinoFiltered.setOrigin(sinogram.getOrigin());
		return sinoFiltered;
	}
	
	public static Grid2D ramLakFilter(Grid2D sinogram) {
		Grid2DComplex sinoFourier = new Grid2DComplex(sinogram);
		sinoFourier.transformForward();
		final int width = sinoFourier.getWidth();
		Grid1DComplex ramLak = new Grid1DComplex(width);
		ramLak.setSpacing(1 / sinogram.getSpacing()[0] * (sinoFourier.getWidth() - sinogram.getWidth()));
		ramLak.setAtIndex(0, 0.25f);
		for (int x = 0; x < width / 2; x++) {
			if (x % 2 != 0) {
				ramLak.setRealAtIndex(x, (float) (-1/Math.pow(Math.PI * x, 2)));
			}
		}
		for (int x = width / 2; x < width; x++) {
			if (x % 2 != 0) {
				ramLak.setRealAtIndex(x, (float) (-1/Math.pow(Math.PI * (width - x), 2)));
			}
		}
		ramLak.show();
		ramLak.transformForward();
		ramLak.show();
		for (int y = 0; y < sinoFourier.getHeight(); y++) {
			for (int x = 0; x < sinoFourier.getWidth(); x++) {
				sinoFourier.multiplyAtIndex(x, y, ramLak.getRealAtIndex(x), ramLak.getImagAtIndex(x));
			}
		}
		sinoFourier.transformInverse();
		Grid2D sinoFiltered = sinoFourier.getRealSubGrid(0, 0, sinogram.getWidth(), sinogram.getHeight());
		sinoFiltered.setSpacing(sinogram.getSpacing());
		sinoFiltered.setOrigin(sinogram.getOrigin());
		return sinoFiltered;
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
		
		RadonTransform radonTransform = new RadonTransform(numProjections, numDetectorPixels, detectorSpacing);
		Grid2D sinogram = radonTransform.transform(phantom);
		ImagePlus sino = VisualizationUtil.showGrid2D(sinogram, "Radon Transform");
		sino.show();
				
		BackProjection backProjection = new BackProjection(size, spacing);
		Grid2D backProjected = backProjection.backProject(sinogram);
		ImagePlus backProj = VisualizationUtil.showGrid2D(backProjected, "Unfiltered Backprojection");
		backProj.show();
		
		Grid2D sinoRampFiltered = rampFilter(sinogram);
		ImagePlus sinoRampFil = VisualizationUtil.showGrid2D(sinoRampFiltered, "Ramp Filtered");
		sinoRampFil.show();

		Grid2D rampFilteredBackProjection = backProjection.backProject(sinoRampFiltered);
		ImagePlus rampFilteredBackProj = VisualizationUtil.showGrid2D(rampFilteredBackProjection, "Ramp-filtered Backprojection");
		rampFilteredBackProj.show();
		
		Grid2D sinoRamLakFiltered = ramLakFilter(sinogram);
		ImagePlus sinoRamLakFil = VisualizationUtil.showGrid2D(sinoRamLakFiltered, "RamLak Filtered");
		sinoRamLakFil.show();
		
		Grid2D ramLakFilteredBackProjection = backProjection.backProject(sinoRamLakFiltered);
		ImagePlus ramLakFilteredBackProj = VisualizationUtil.showGrid2D(ramLakFilteredBackProjection, "RamLak-filtered Backprojection");
		ramLakFilteredBackProj.show();	


//		ParallelBackprojector2D backProjection2 = new ParallelBackprojector2D(256, 256, 1, 1);
//		backProjection2.backprojectPixelDriven(sinogram).show("The Reconstruction");
		

	}

}
