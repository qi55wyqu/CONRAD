typedef float TvoxelValue;
typedef float Tcoord_dev;

__constant sampler_t sampler = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_LINEAR;

__kernel void gridAddKernel(__read_only image2d_t g1Tex, __read_only image2d_t g2Tex, __global TvoxelValue* gRes, __constant Tcoord_dev* gVolumeSize) {

	int gidx = get_group_id(0);
	int gidy = get_group_id(1);
	int lidx = get_group_id(0);
	int lidy = get_group_id(1);
	
	int locSizex = get_local_size(0);
	int locSizey = get_local_size(1);
	
	int x = mad24(gidx, locSizex, lidx);
	int y = mad24(gidy, locSizey, lidy);
	
	unsigned int yStride = gVolumeSize[0];
	
	if (x >= gVolumeSize[0] || y >= gVolumeSize[1]) {
		return;
	}
	
	unsigned long idx = y * yStride + x;
	
	float val1 = read_imagef(g1Tex, sampler, (float2)(x+0.5f, y+0.5f)).x;
	float val2 = read_imagef(g2Tex, sampler, (float2)(x+0.5f, y+0.5f)).x;
	
	gRes[idx] = val1 + val2;
	
	return;

}