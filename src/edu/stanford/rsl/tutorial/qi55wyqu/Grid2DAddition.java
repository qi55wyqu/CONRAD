package edu.stanford.rsl.tutorial.qi55wyqu;

import java.io.IOException;
import java.nio.FloatBuffer;

import Ice.InputStream;

import com.jogamp.opencl.CLBuffer;
import com.jogamp.opencl.CLCommandQueue;
import com.jogamp.opencl.CLContext;
import com.jogamp.opencl.CLDevice;
import com.jogamp.opencl.CLImage2d;
import com.jogamp.opencl.CLImageFormat;
import com.jogamp.opencl.CLImageFormat.ChannelOrder;
import com.jogamp.opencl.CLImageFormat.ChannelType;
import com.jogamp.opencl.CLKernel;
import com.jogamp.opencl.CLProgram;
import com.jogamp.opencl.CLMemory.Mem;

import edu.stanford.rsl.conrad.data.numeric.Grid2D;
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
		OpenCLGrid2D image = inputImage.clone();
		float[] imgSize = {inputImage.getWidth(), inputImage.getHeight()};
		CLContext context = OpenCLUtil.getStaticContext();
		CLDevice device = context.getMaxFlopsDevice();
		CLCommandQueue commandQueue = device.createCommandQueue();
		java.io.InputStream is = OpenCLGridTest.class.getResourceAsStream("openCLGriddAdd.cl");
		try {
			CLProgram program = context.createProgram(is).build();
		} catch (IOException e) {
			e.printStackTrace();
		}
		CLKernel kernelFunction = program.createCLKernel("griddAddKernel");
		CLBuffer<FloatBuffer> gImgSize = context.createFloatBuffer(imgSize.length, Mem.READ_ONLY);
		gImgSize.getBuffer().put(imgSize);
		gImgSize.getBuffer().rewind();
		CLImageFormat format = new CLImageFormat(ChannelOrder.INTENSITY, ChannelType.FLOAT);
		inputImage.getDelegate().prepareForDeviceOperation();
		inputImage.getDelegate().getCLBuffer().getBuffer().rewind();
		CLImage2d<FloatBuffer> inputImageTex = null;
		inputImageTex = context.createImage2d(inputImage.getDelegate().getCLBuffer().getBuffer(), inputImage.getWidth(), inputImage.getHeight(), format, Mem.READ_ONLY);
		inputImage.getDelegate().release();
		image.getDelegate().prepareForDeviceOperation();
		image.getDelegate().getCLBuffer().getBuffer().rewind();
		CLImage2d<FloatBuffer> imageTex = null;
		imageTex = context.createImage2d(image.getDelegate().getCLBuffer().getBuffer(), image.getWidth(), image.getHeight(), format, Mem.READ_ONLY);
		image.getDelegate().release();
		commandQueue
			.putWriteImage(inputImageTex, true)
			.putWriteImage(image, true);
			.putWriteBuffer(gImgSize, true)
			.finish();
		kernelFunction.rewind()
		kernelFunction
			.putArg(inputImageTex)
			.putArg(imageTex)
			.putArgs(gImgSize);
		int bpBlockSize[] = {32, 32};
		int maxWorkGroupSize = device.getMaxWorkGroupSize();
		int[] realLocalSize = new int[] {Math.min((int) Math.pow(maxWorkGroupSize, 1/2.0), bpBlockSize[0]), Math.min((int) Math.pow(maxWorkGroupSize, 1/2.0), bpBlockSize[1])};
		int[] globalWorkSize = new int[] {OpenCLUtil.roundUp(realLocalSize[0], (int) imgSize[0]), OpenCLUtil.roundUp(realLocalSize[1], (int) imgSize[1])};	
		commandQueue.put2DRangeKernel(kernelFunction, 0, 0, globalWorkSize[0], globalWorkSize[1], realLocalSize[0], realLocalSize[1]).finish();
		
	}
	
	public static void GPUAddition(OpenCLGrid2D inputImage) {
		GPUAddition(inputImage, 1000000);
	}
	
	public static void main(String[] args) {

		int[] size = {1024, 1024};
		double[] spacing = {0.5 ,0.75};
		Grid2D phantom = new Phantom(size, spacing);
		CPUAddition(phantom);
		
		OpenCLGrid2D openCLPhantom = new OpenCLPhantom(size, spacing);
		GPUAddition(openCLPhantom);
		
	}
	
}
