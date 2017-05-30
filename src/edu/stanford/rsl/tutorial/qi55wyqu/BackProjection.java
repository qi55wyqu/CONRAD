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
	
	final static boolean debug = true;
	
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
		
		double[][] thetaArray = new double[sinogram.getHeight()][3];
		for (int thetaIndex = 0; thetaIndex < thetaArray.length; thetaIndex++) {
			double theta = thetaIndex * sinogram.getSpacing()[1];
			thetaArray[thetaIndex][0] = Math.cos(Math.toRadians(theta));
			thetaArray[thetaIndex][1] = Math.sin(Math.toRadians(theta));
			thetaArray[thetaIndex][2] = theta;
		}
		
		for (int y = 0; y < backProjection.getHeight(); y++) {
			if (y % 50 == 0) System.out.println("Backprojection running... line " + y);
			for (int x = 0; x < backProjection.getWidth(); x++) {
				double[] backProjPhysical = backProjection.indexToPhysical(x, y);
				float sum = 0.f;
				for (int thetaIndex = 0; thetaIndex < sinogram.getHeight(); thetaIndex++) {
					double s = backProjPhysical[0] * thetaArray[thetaIndex][0] + backProjPhysical[1] * thetaArray[thetaIndex][1];
					double[] sinoIdx = sinogram.physicalToIndex(s, thetaArray[thetaIndex][2]);
					sum += InterpolationOperators.interpolateLinear(sinogram, sinoIdx[0], sinoIdx[1]);
				}
				backProjection.setAtIndex(x, y, sum);
			}
		}
		return backProjection;
	}
	
	public static Grid2D rampFilter(Grid2D sinogram) {
//		Grid2DComplex sinoFourier = new Grid2DComplex(sinogram, false);
//		sinoFourier.transformForward();
		Grid1DComplex rampFilter = new Grid1DComplex(sinogram.getWidth(), false);
		
//		rampFilter.setSpacing(1 / (sinogram.getSpacing()[0] * (sinoFourier.getWidth() - sinogram.getWidth())));
//		int width = sinogram.getWidth();
		int width = rampFilter.getSize()[0];
		for (int x = 0; x < width / 2; x++) {
			rampFilter.setAtIndex(x, x);
		}
		for (int x = width / 2; x < width; x++) {
			rampFilter.setAtIndex(x, width - x);
		}
		if (debug) rampFilter.show();
		for (int y = 0; y < sinogram.getHeight(); y++) {
			Grid1DComplex projection = new Grid1DComplex(sinogram.getSubGrid(y), false);
			projection.transformForward();
			for (int x = 0; x < sinogram.getWidth(); x++) {
				projection.multiplyAtIndex(x, rampFilter.getRealAtIndex(x), rampFilter.getImagAtIndex(x));
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
		Grid1DComplex ramLak = new Grid1DComplex(sinogram.getWidth());
		ramLak.setSpacing(1 / sinogram.getSpacing()[0] * (sinoFourier.getWidth() - sinogram.getWidth()));
		ramLak.setAtIndex(0, 0.25f);
		final int width = ramLak.getSize()[0];
		for (int x = 1; x < width / 2; x++) {
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
		
		int[] size = new int[] { 128, 128 };
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
		if (debug) {
			ImagePlus sinoRampFil = VisualizationUtil.showGrid2D(sinoRampFiltered, "Ramp Filtered");
			sinoRampFil.show();
		}

		Grid2D rampFilteredBackProjection = backProjection.backProject(sinoRampFiltered);
		ImagePlus rampFilteredBackProj = VisualizationUtil.showGrid2D(rampFilteredBackProjection, "Ramp-filtered Backprojection");
		rampFilteredBackProj.show();
		
		Grid2D sinoRamLakFiltered = ramLakFilter(sinogram);
		if (debug) {
			ImagePlus sinoRamLakFil = VisualizationUtil.showGrid2D(sinoRamLakFiltered, "RamLak Filtered");
			sinoRamLakFil.show();
		}
		
		Grid2D ramLakFilteredBackProjection = backProjection.backProject(sinoRamLakFiltered);
		ImagePlus ramLakFilteredBackProj = VisualizationUtil.showGrid2D(ramLakFilteredBackProjection, "RamLak-filtered Backprojection");
		ramLakFilteredBackProj.show();	


//		ParallelBackprojector2D backProjection2 = new ParallelBackprojector2D(256, 256, 1, 1);
//		backProjection2.backprojectPixelDriven(sinogram).show("The Reconstruction");
		

	}

}
