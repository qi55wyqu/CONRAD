__constant sampler_t sampler = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_LINEAR;

__kernel void backProjectionKernel(
	__read_only image2d_t sinogram, 
	__constant float* thetaIndexMax, 
	__constant float* sinogramSpacing, 
	__constant float* sinogramOrigin, 
	__constant float* backProjSpacing, 
	__constant float* backProjOrigin, 
	__global float* backProjection, 
	__constant float* backProjSize) {

	for (int thetaIndex = 0; thetaIndex < thetaIndexMax; thetaIndex++) {
	
		int IDx = get_group_id(0);
		int IDy = get_group_id(1);
		
		int x = mad24(IDx, get_local_size(0), IDx);
		int y = mad24(IDy, get_local_size(1), IDy);
		
		if (x >= backProjSize[0] || y >= backProjSize[1])
			return;
		
		float theta = thetaIndex * sinogramSpacing[1];
		float sinTheta = sin(radians(theta));
		float cosTheta = cos(radians(theta));
		float s = (x * backProjSpacing[0] + backProjOrigin[0]) * cosTheta + 
				(y * backProjSpacing[1] + backProjOrigin[1]) * sinTheta;
		float sinoIdxX = 1 / sinogramSpacing[0] * (s - sinogramOrigin[0]);
		float sinoIdxY = 1 / sinogramSpacing[1] * (theta - sinogramOrigin[1]);

		unsigned long index = y * backProjSize[0] + x;
			
		backProjection[index] += read_imagef(sinogram, sampler, (float2) (sinoIdxX, sinoIdxY)).x;
		
	}

	return;

}

