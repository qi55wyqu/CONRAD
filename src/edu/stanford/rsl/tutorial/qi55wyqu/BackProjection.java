package edu.stanford.rsl.tutorial.qi55wyqu;

import java.io.IOException;
import java.nio.FloatBuffer;

import com.jogamp.opencl.CLBuffer;
import com.jogamp.opencl.CLCommandQueue;
import com.jogamp.opencl.CLContext;
import com.jogamp.opencl.CLDevice;
import com.jogamp.opencl.CLImage2d;
import com.jogamp.opencl.CLImageFormat;
import com.jogamp.opencl.CLKernel;
import com.jogamp.opencl.CLProgram;
import com.jogamp.opencl.CLImageFormat.ChannelOrder;
import com.jogamp.opencl.CLImageFormat.ChannelType;
import com.jogamp.opencl.CLMemory.Mem;
import edu.stanford.rsl.conrad.data.numeric.Grid1D;
import edu.stanford.rsl.conrad.data.numeric.Grid1DComplex;
import edu.stanford.rsl.conrad.data.numeric.Grid2D;
import edu.stanford.rsl.conrad.data.numeric.Grid2DComplex;
import edu.stanford.rsl.conrad.data.numeric.InterpolationOperators;
import edu.stanford.rsl.conrad.data.numeric.opencl.OpenCLGrid2D;
import edu.stanford.rsl.conrad.opencl.OpenCLUtil;
import edu.stanford.rsl.conrad.utils.VisualizationUtil;
import edu.stanford.rsl.tutorial.filters.RamLakKernel;
import edu.stanford.rsl.tutorial.parallel.ParallelBackprojector2D;
import ij.ImageJ;
import ij.ImagePlus;

public class BackProjection {
	
	final static boolean debug = false;
	
	protected int[] size;
	protected double[] spacing;
	protected double[] origin;
	
	public BackProjection(int[] size, double[] spacing) {
		this.size = size;
		this.spacing = spacing;
		this.origin = new double[] {-(this.size[0] - 1) * this.spacing[0] / 2, -(this.size[1] - 1) * this.spacing[1] / 2};
	}
		
	public Grid2D backProject(Grid2D sinogram) {
		Grid2D backProjection = new Grid2D(this.size[0], this.size[1]);
		backProjection.setSpacing(this.spacing);
		backProjection.setOrigin(this.origin);
		for (int thetaIndex = 0; thetaIndex < sinogram.getHeight(); thetaIndex++) {
			double theta = thetaIndex * sinogram.getSpacing()[1];
//			System.out.println("Backprojection: \u03B8 = " + theta);
			double sinTheta = Math.sin(Math.toRadians(theta));
			double cosTheta = Math.cos(Math.toRadians(theta));
			for (int y = 0; y < backProjection.getHeight(); y++) {
				for (int x = 0; x < backProjection.getWidth(); x++) {
					double[] backProjPhysical = backProjection.indexToPhysical(x, y);
					double s = backProjPhysical[0] * cosTheta + backProjPhysical[1] * sinTheta;
					double[] sinoIdx = sinogram.physicalToIndex(s, theta);
					backProjection.addAtIndex(x, y, InterpolationOperators.interpolateLinear(sinogram, sinoIdx[0], sinoIdx[1]));
				}
			}
		}
		return backProjection;
	}
	
	public Grid2D backProjectNoHelp(Grid2D sinogram) {
		Grid2D backProjection = new Grid2D(this.size[0], this.size[1]);
		int[] backProjSize = {backProjection.getWidth(), backProjection.getHeight()};
		for (int thetaIndex = 0; thetaIndex < sinogram.getHeight(); thetaIndex++) {
			backProjection = backProjStep(sinogram, thetaIndex, sinogram.getSpacing(), sinogram.getOrigin(), this.spacing, this.origin, backProjection, backProjSize);
		}
		backProjection.setSpacing(this.spacing);
		backProjection.setOrigin(this.origin);
		return backProjection;
	}
	
	private static Grid2D backProjStep(Grid2D sinogram, double thetaIndex, double[] sinogramSpacing, double[] sinogramOrigin, double[] backProjSpacing, double[] backProjOrigin, Grid2D backProjection, int[] backProjSize) {
		double theta = thetaIndex * sinogramSpacing[1];
		double sinTheta = Math.sin(Math.toRadians(theta));
		double cosTheta = Math.cos(Math.toRadians(theta));
		for (int y = 0; y < backProjSize[1]; y++) {
			for (int x = 0; x < backProjSize[0]; x++) {
				double s = (x * backProjSpacing[0] + backProjOrigin[0]) * cosTheta + 
						   (y * backProjSpacing[1] + backProjOrigin[1]) * sinTheta;
				backProjection.addAtIndex(x, y, 
					InterpolationOperators.interpolateLinear(
						sinogram, 
						(double) (1 / sinogramSpacing[0] * (s - sinogramOrigin[0])), 
						(double) (1 / sinogramSpacing[1] * (theta - sinogramOrigin[1]))));
			}
		}
		return backProjection;
	}
	
	public Grid2D backProjectTest(Grid2D sinogram) {
		Grid2D backProjection = new Grid2D(this.size[0], this.size[1]);
		backProjection.setSpacing(this.spacing);
		backProjection.setOrigin(this.origin);
		double[][][] backProjPhysical = new double[backProjection.getWidth()][backProjection.getHeight()][2];
		for (int y = 0; y < backProjection.getHeight(); y++) {
			for (int x = 0; x < backProjection.getWidth(); x++) {
				double[] physical = backProjection.indexToPhysical(x, y);
				for (int i = 0; i < 2; i++) {
					backProjPhysical[y][x][i] = physical[i];
				}
			}
		}
		for (int thetaIndex = 0; thetaIndex < sinogram.getHeight(); thetaIndex++) {
			double theta = thetaIndex * sinogram.getSpacing()[1];
//			System.out.println("Backprojection: \u03B8 = " + theta);
			double sinTheta = Math.sin(Math.toRadians(theta));
			double cosTheta = Math.cos(Math.toRadians(theta));
			for (int y = 0; y < backProjection.getHeight(); y++) {
				for (int x = 0; x < backProjection.getWidth(); x++) {
					double s = backProjPhysical[y][x][0] * cosTheta + backProjPhysical[y][x][1] * sinTheta;
					double[] sinoIdx = sinogram.physicalToIndex(s, theta);
					backProjection.addAtIndex(x, y, InterpolationOperators.interpolateLinear(sinogram, sinoIdx[0], sinoIdx[1]));
				}
			}
		}
		return backProjection;
	}
	
	public Grid2D backProjectOnGPU(Grid2D sinogram) {
		
		OpenCLGrid2D gSinogram = new OpenCLGrid2D(sinogram);
		OpenCLGrid2D gBackProjection = new OpenCLGrid2D(new Grid2D(this.size[0], this.size[1]));
		
		float[] thetaIndexMax = {sinogram.getHeight(), sinogram.getHeight()};
		float[] sinogramSpacing = {(float) sinogram.getSpacing()[0], (float) sinogram.getSpacing()[1]}; 
		float[] sinogramOrigin = {(float) sinogram.getOrigin()[0], (float) sinogram.getOrigin()[1]};
		float[] backProjOrigin = {(float) (-(this.size[0] - 1) * this.spacing[0] / 2), (float) (-(this.size[1] - 1) * this.spacing[1] / 2)};
		float[] backProjSpacing = {(float) this.spacing[0], (float) this.spacing[1]};		
		float[] backProjSize = {(float) this.size[0], (float) this.size[1]};

		CLContext context = OpenCLUtil.getStaticContext();
		CLDevice device = context.getMaxFlopsDevice();
		CLCommandQueue commandQueue = device.createCommandQueue();
		java.io.InputStream is = BackProjection.class.getResourceAsStream("backProjectionKernel.cl");
		CLProgram program = null;
		try {
			program = context.createProgram(is).build();
		} catch (IOException e) {
			e.printStackTrace();
		}
		if (program == null) {
			System.err.println("Build failed");
			return new Grid2D(this.size[0], this.size[1]);
		}
		CLKernel kernelFunction = program.createCLKernel("backProjectionKernel");
		
		CLBuffer<FloatBuffer> gThetaIndexMax = context.createFloatBuffer(thetaIndexMax.length, Mem.READ_ONLY);
		gThetaIndexMax.getBuffer().put(thetaIndexMax);
		gThetaIndexMax.getBuffer().rewind();
		
		CLBuffer<FloatBuffer> gSinogramSpacing = context.createFloatBuffer(sinogramSpacing.length, Mem.READ_ONLY);
		gSinogramSpacing.getBuffer().put(sinogramSpacing);
		gSinogramSpacing.getBuffer().rewind();
		
		CLBuffer<FloatBuffer> gSinogramOrigin = context.createFloatBuffer(sinogramOrigin.length, Mem.READ_ONLY);
		gSinogramOrigin.getBuffer().put(sinogramOrigin);
		gSinogramOrigin.getBuffer().rewind();
		
		CLBuffer<FloatBuffer> gBackProjSpacing = context.createFloatBuffer(backProjSpacing.length, Mem.READ_ONLY);
		gBackProjSpacing.getBuffer().put(backProjSpacing);
		gBackProjSpacing.getBuffer().rewind();
		
		CLBuffer<FloatBuffer> gBackProjOrigin = context.createFloatBuffer(backProjOrigin.length, Mem.READ_ONLY);
		gBackProjOrigin.getBuffer().put(backProjOrigin);
		gBackProjOrigin.getBuffer().rewind();
		
		CLBuffer<FloatBuffer> gBackProjSize = context.createFloatBuffer(backProjSize.length, Mem.READ_ONLY);
		gBackProjSize.getBuffer().put(backProjSize);
		gBackProjSize.getBuffer().rewind();

		CLImageFormat format = new CLImageFormat(ChannelOrder.INTENSITY, ChannelType.FLOAT);
		
		gSinogram.getDelegate().prepareForDeviceOperation();
		gSinogram.getDelegate().getCLBuffer().getBuffer().rewind();
		CLImage2d<FloatBuffer> gSinogramTex = null;
		gSinogramTex = context.createImage2d(gSinogram.getDelegate().getCLBuffer().getBuffer(), gSinogram.getWidth(), gSinogram.getHeight(), format, Mem.READ_ONLY);
		gSinogram.getDelegate().release();
				
//		__kernel void backProjectionKernel(
//				__read_only image2d_t sinogram, 
//				__constant float thetaIndex, 
//				__constant float* sinogramSpacing, 
//				__constant float* sinogramOrigin, 
//				__constant float* backProjSpacing, 
//				__constant float* backProjOrigin, 
//				__global float* backProjection, 
//				_constant float* backProjSize) {
		commandQueue
		.putWriteImage(gSinogramTex, true)
		.putWriteBuffer(gThetaIndexMax, true)
		.putWriteBuffer(gSinogramSpacing, true)
		.putWriteBuffer(gSinogramOrigin, true)
		.putWriteBuffer(gBackProjSpacing, true)
		.putWriteBuffer(gBackProjOrigin, true)
		.putWriteBuffer(gBackProjection.getDelegate().getCLBuffer(), true)
		.putWriteBuffer(gBackProjSize, true)
		.finish();
		
		kernelFunction.rewind();
		kernelFunction
			.putArg(gSinogramTex)
			.putArg(gThetaIndexMax)
			.putArg(gSinogramSpacing)
			.putArg(gSinogramOrigin)
			.putArg(gBackProjSpacing)
			.putArg(gBackProjOrigin)
			.putArg(gBackProjection.getDelegate().getCLBuffer())
			.putArg(gBackProjSize);

		int bpBlockSize[] = {32, 32};
		int maxWorkGroupSize = device.getMaxWorkGroupSize();
		int[] realLocalSize = new int[] {Math.min((int) Math.pow(maxWorkGroupSize, 1/2.0), bpBlockSize[0]), Math.min((int) Math.pow(maxWorkGroupSize, 1/2.0), bpBlockSize[1])};
		int[] globalWorkSize = new int[] {OpenCLUtil.roundUp(realLocalSize[0], (int) backProjSize[0]), OpenCLUtil.roundUp(realLocalSize[1], (int) backProjSize[1])};	

		commandQueue.put2DRangeKernel(kernelFunction, 0, 0, globalWorkSize[0], globalWorkSize[1], realLocalSize[0], realLocalSize[1]).finish();
			
		gBackProjection.getDelegate().notifyDeviceChange();
		
		Grid2D backProjection = new Grid2D(gBackProjection);
		backProjection.setSpacing(this.spacing);
		backProjection.setOrigin(new double[] {-(this.size[0] - 1) * this.spacing[0] / 2, -(this.size[1] - 1) * this.spacing[1] / 2});
		
		return backProjection;
		
	}
	
	public static Grid2D rampFilter(Grid2D sinogram) {
		int width = sinogram.getWidth();
		Grid1DComplex rampFilter = new Grid1DComplex(width, false);
		rampFilter.setSpacing(1 / (sinogram.getSpacing()[0] * (rampFilter.getSize()[0] - width)));
		for (int x = 0; x < width / 2; x++) {
			rampFilter.setRealAtIndex(x, x);
		}
		for (int x = width / 2; x < width; x++) {
			rampFilter.setRealAtIndex(x, (width - x));
		}
		if (debug) rampFilter.show();
		Grid2D sinoFiltered = new Grid2D(sinogram.getWidth(), sinogram.getHeight());
		sinoFiltered.setSpacing(sinogram.getSpacing());
		sinoFiltered.setOrigin(sinogram.getOrigin());
		for (int y = 0; y < sinogram.getHeight(); y++) {
			Grid1DComplex projection = new Grid1DComplex(sinogram.getSubGrid(y), false);
			projection.transformForward();
			for (int x = 0; x < projection.getSize()[0]; x++) {
				projection.multiplyAtIndex(x, rampFilter.getRealAtIndex(x), rampFilter.getImagAtIndex(x));
			}
			projection.transformInverse();
			for (int x = 0; x < sinoFiltered.getWidth(); x++) {
				sinoFiltered.setAtIndex(x, y, projection.getRealAtIndex(x));
			}
		}
		return sinoFiltered;
	}
	
	public static Grid2D ramLakFilter(Grid2D sinogram) {
		int width = sinogram.getWidth();
		Grid1DComplex ramLak = new Grid1DComplex(width, false);
		ramLak.setSpacing(1 / sinogram.getSpacing()[0] * (ramLak.getSize()[0] - width));
		ramLak.setAtIndex(0, 0.25f);
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
		if (debug) ramLak.show();
		ramLak.transformForward();
		if (debug) ramLak.show();
		Grid2D sinoFiltered = new Grid2D(sinogram.getWidth(), sinogram.getHeight());
		sinoFiltered.setSpacing(sinogram.getSpacing());
		sinoFiltered.setOrigin(sinogram.getOrigin());
		for (int y = 0; y < sinogram.getHeight(); y++) {
			Grid1DComplex projection = new Grid1DComplex(sinogram.getSubGrid(y), false);
			projection.transformForward();
			for (int x = 0; x < width; x++) {
				projection.multiplyAtIndex(x, ramLak.getRealAtIndex(x), ramLak.getImagAtIndex(x));
			}
			projection.transformInverse();
			for (int x = 0; x < sinoFiltered.getWidth(); x++) {
				sinoFiltered.setAtIndex(x, y, projection.getRealAtIndex(x));
			}
		}
		return sinoFiltered;
	}

	public static void main(String[] args) {
		
		int[] size = new int[] { 256, 256 };
		double[] spacing = new double[] {1.0, 1.0};
		int numProjections = 180;
		int numDetectorPixels = 2 * size[0];
		double detectorSpacing = 1.0;
		
		new ImageJ();
		
		Phantom phantom = new Phantom(size, spacing);
//		ImagePlus phan = VisualizationUtil.showGrid2D(phantom, "Test Phantom");
//		phan.show();
		
		RadonTransform radonTransform = new RadonTransform(numProjections, numDetectorPixels, detectorSpacing);
		Grid2D sinogram = radonTransform.transform(phantom);
//		ImagePlus sino = VisualizationUtil.showGrid2D(sinogram, "Radon Transform");
//		sino.show();
				
		BackProjection backProjection = new BackProjection(size, spacing);
//		Grid2D backProjected = backProjection.backProject(sinogram);
//		ImagePlus backProj = VisualizationUtil.showGrid2D(backProjected, "Unfiltered Backprojection");
//		backProj.show();
//		
//		Grid2D sinoRampFiltered = rampFilter(sinogram);
//		if (debug) {
//			ImagePlus sinoRampFil = VisualizationUtil.showGrid2D(sinoRampFiltered, "Ramp Filtered");
//			sinoRampFil.show();
//		}
//
//		Grid2D rampFilteredBackProjection = backProjection.backProject(sinoRampFiltered);
//		ImagePlus rampFilteredBackProj = VisualizationUtil.showGrid2D(rampFilteredBackProjection, "Ramp-filtered Backprojection");
//		rampFilteredBackProj.show();
		
		Grid2D sinoRamLakFiltered = ramLakFilter(sinogram);
		if (debug) {
			ImagePlus sinoRamLakFil = VisualizationUtil.showGrid2D(sinoRamLakFiltered, "RamLak Filtered");
			sinoRamLakFil.show();
		}

//		Grid2D ramLakFilteredBackProjection = backProjection.backProjectNoHelp(sinoRamLakFiltered);
		Grid2D ramLakFilteredBackProjection = backProjection.backProjectOnGPU(sinoRamLakFiltered);
		ImagePlus ramLakFilteredBackProj = VisualizationUtil.showGrid2D(ramLakFilteredBackProjection, "RamLak-filtered Backprojection");
		ramLakFilteredBackProj.show();

//		ParallelBackprojector2D backProjection2 = new ParallelBackprojector2D(256, 256, 1, 1);
//		backProjection2.backprojectPixelDriven(sinogram).show("The Reconstruction");
		

	}

}
