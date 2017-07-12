__constant sampler_t sampler = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_LINEAR;

__kernel void gridAddKernel(__read_only image2d_t gInputImageTex, __global float* gResult, __constant float* gImageSize) {

	int idx = get_group_id(0);
	int idy = get_group_id(1);
	
	int x = mad24(idx, get_local_size(0), idx);
	int y = mad24(idy, get_local_size(1), idy);
	
	if (x >= gImageSize[0] || y >= gImageSize[1])
		return;
	
	unsigned long index = y * gImageSize[0] + x;
		
	gResult[index] += read_imagef(gInputImageTex, sampler, (float2) (x+0.5f, y+0.5f)).x;
	
	return;

}