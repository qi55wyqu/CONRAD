package edu.stanford.rsl.tutorial.qi55wyqu;

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
		Grid1DComplex rampFilter = new Grid1DComplex(sinoFourier.getWidth());
		rampFilter.setSpacing(1 / (sinogram.getSpacing()[0] * (sinoFourier.getWidth()-sinogram.getWidth())));
		for (int x = 0; x <= rampFilter.getSize()[0] / 2; x++) {
			rampFilter.setAtIndex((int) (x/2), x);
			rampFilter.setAtIndex((int)((rampFilter.getSize()[0]-1-x)/2), x);
		}
		rampFilter.show();
		for (int y = 0; y < sinoFourier.getHeight(); y++) {
			for (int x = 0; x <= sinoFourier.getWidth() / 2; x++) {
				sinoFourier.multiplyAtIndex(x, y, rampFilter.getRealAtIndex(x), rampFilter.getImagAtIndex(x));
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
		
		Grid2D sinoFiltered = rampFilter(sinogram);
		ImagePlus sinoFil = VisualizationUtil.showGrid2D(sinoFiltered, "Ramp Filtered");
		sinoFil.show();
		
		BackProjection backProjection = new BackProjection(size, spacing);
		Grid2D backProjected = backProjection.backProject(sinogram);
		ImagePlus backProj = VisualizationUtil.showGrid2D(backProjected, "Unfiltered Backprojection");
		backProj.show();
		
		Grid2D filteredBackProjection = backProjection.backProject(sinoFiltered);
		ImagePlus filteredBackProj = VisualizationUtil.showGrid2D(filteredBackProjection, "Filtered Backprojection");
		filteredBackProj.show();


//		ParallelBackprojector2D backProjection2 = new ParallelBackprojector2D(256, 256, 1, 1);
//		backProjection2.backprojectPixelDriven(sinogram).show("The Reconstruction");
		

	}

}
