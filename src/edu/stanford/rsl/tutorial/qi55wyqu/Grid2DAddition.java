package edu.stanford.rsl.tutorial.qi55wyqu;

import java.io.IOException;
import java.nio.FloatBuffer;

import Ice.InputStream;

import com.jogamp.opencl.CLBuffer;
import com.jogamp.opencl.CLCommandQueue;
import com.jogamp.opencl.CLContext;
import com.jogamp.opencl.CLDevice;
import com.jogamp.opencl.CLImage;
import com.jogamp.opencl.CLImage2d;
import com.jogamp.opencl.CLImageFormat;
import com.jogamp.opencl.CLImageFormat.ChannelOrder;
import com.jogamp.opencl.CLImageFormat.ChannelType;
import com.jogamp.opencl.CLKernel;
import com.jogamp.opencl.CLProgram;
import com.jogamp.opencl.CLMemory.Mem;

import edu.stanford.rsl.conrad.data.numeric.Grid1D;
import edu.stanford.rsl.conrad.data.numeric.Grid2D;
import edu.stanford.rsl.conrad.data.numeric.opencl.OpenCLGrid1D;
import edu.stanford.rsl.conrad.data.numeric.opencl.OpenCLGrid2D;
import edu.stanford.rsl.conrad.data.numeric.opencl.OpenCLGridTest;
import edu.stanford.rsl.conrad.opencl.OpenCLUtil;

public class Grid2DAddition {
	
	public static void CPUAddition(Grid2D inputImage, int numIterations) {
		Grid2D image = new Grid2D(inputImage);
		long startTime = System.currentTimeMillis();
		for (int i = 0; i < numIterations; i++) {
			for (int y = 0; y < inputImage.getHeight(); y++) {
				for (int x = 0; x < inputImage.getWidth(); x++) {
					image.addAtIndex(x, y, inputImage.getAtIndex(x, y));
				}
			}
		}
		long endTime = System.currentTimeMillis() - startTime;
		System.out.println("Time taken for CPU addition: " + endTime + " ms");
	}
	
	public static void CPUAddition(Grid2D inputImage) {
		CPUAddition(inputImage, 1000000);
	}
	
	public static void GPUAddition(OpenCLGrid2D inputImage, int numIterations) {
//		OpenCLGrid2D imageCopy = inputImage.clone();
		float[] imgSize = new float[2];
		imgSize[0] = inputImage.getSize()[0];
		imgSize[1] = inputImage.getSize()[1];
		OpenCLGrid2D res = new OpenCLGrid2D(new Grid2D (inputImage.getSize()[0], inputImage.getSize()[1]));
		CLContext context = OpenCLUtil.getStaticContext();
		CLDevice device = context.getMaxFlopsDevice();
		CLCommandQueue commandQueue = device.createCommandQueue();
		java.io.InputStream is = Grid2DAddition.class.getResourceAsStream("gridAddKernel.cl");
		CLProgram program = null;
		try {
			program = context.createProgram(is).build();
		} catch (IOException e) {
			e.printStackTrace();
		}
		if (program == null) {
			System.err.println("Build failed");
			return;
		}
		
		CLKernel kernelFunction = program.createCLKernel("gridAddKernel");

		CLBuffer<FloatBuffer> gImgSize = context.createFloatBuffer(imgSize.length, Mem.READ_ONLY);
		gImgSize.getBuffer().put(imgSize);
		gImgSize.getBuffer().rewind();

		CLImageFormat format = new CLImageFormat(ChannelOrder.INTENSITY, ChannelType.FLOAT);
		
		inputImage.getDelegate().prepareForDeviceOperation();
		inputImage.getDelegate().getCLBuffer().getBuffer().rewind();
		CLImage2d<FloatBuffer> inputImageTex = null;
		inputImageTex = context.createImage2d(inputImage.getDelegate().getCLBuffer().getBuffer(), inputImage.getWidth(), inputImage.getHeight(), format, Mem.READ_WRITE);
		inputImage.getDelegate().release();

//		imageCopy.getDelegate().prepareForDeviceOperation();
//		imageCopy.getDelegate().getCLBuffer().getBuffer().rewind();
//		CLImage2d<FloatBuffer> imageCopyTex = null;
//		imageCopyTex = context.createImage2d(imageCopy.getDelegate().getCLBuffer().getBuffer(), imageCopy.getWidth(), imageCopy.getHeight(), format, Mem.READ_ONLY);
//		imageCopy.getDelegate().release();
		
//		res.getDelegate().prepareForDeviceOperation();
//		res.getDelegate().getCLBuffer().getBuffer().rewind();
//		CLImage2d<FloatBuffer> resTex = null;
//		resTex = context.createImage2d(res.getDelegate().getCLBuffer().getBuffer(), res.getWidth(), res.getHeight(), format, Mem.READ_ONLY);
//		res.getDelegate().release();
		
		commandQueue
			.putWriteImage(inputImageTex, true)
//			.putWriteImage(imageCopyTex, true)
			.putWriteBuffer(res.getDelegate().getCLBuffer(), true)
			.putWriteBuffer(gImgSize, true)
			.finish();
		kernelFunction.rewind();
		
		kernelFunction
			.putArg(inputImageTex)
//			.putArg(imageCopyTex)
			.putArg(res.getDelegate().getCLBuffer())
			.putArgs(gImgSize);
		
		int bpBlockSize[] = {32, 32};
		int maxWorkGroupSize = device.getMaxWorkGroupSize();
		int[] realLocalSize = new int[] {Math.min((int) Math.pow(maxWorkGroupSize, 1/2.0), bpBlockSize[0]), Math.min((int) Math.pow(maxWorkGroupSize, 1/2.0), bpBlockSize[1])};
		int[] globalWorkSize = new int[] {OpenCLUtil.roundUp(realLocalSize[0], (int) imgSize[0]), OpenCLUtil.roundUp(realLocalSize[1], (int) imgSize[1])};	
		
		long startTime = System.currentTimeMillis();
		
		for (int i = 0; i < numIterations; i++) {
			commandQueue.put2DRangeKernel(kernelFunction, 0, 0, globalWorkSize[0], globalWorkSize[1], realLocalSize[0], realLocalSize[1]).finish();
		}
		
		long endTime = System.currentTimeMillis() - startTime;
		System.out.println("Time taken for GPU addition: " + endTime + " ms");
		
		res.getDelegate().notifyDeviceChange();
		
	}
	
	public static void GPUAddition(OpenCLGrid2D inputImage) {
		GPUAddition(inputImage, 1000000);
	}
	
	public static void main(String[] args) {

		int[] size = {1024, 1024};
		double[] spacing = {0.5 ,0.75};
		Grid2D phantom = new Phantom(size, spacing);
		CPUAddition(phantom, 1000);
		
//		Grid1D phantom1D = new Grid1D(size[0]*size[1]);
//		for (int y = 0; y < phantom.getHeight(); y++) {
//			for (int x = 0; x < phantom.getWidth(); x++) {
//				phantom1D.setAtIndex(y * phantom.getWidth() + x, phantom.getAtIndex(x, y));
//			}
//		}
//		OpenCLGrid1D openCLPhantom = new OpenCLGrid1D(phantom1D);
		OpenCLGrid2D openCLPhantom = new OpenCLPhantom(size, spacing);
		GPUAddition(openCLPhantom, 1000);
		
	}
	
}
