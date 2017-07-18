__constant sampler_t sampler = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_LINEAR;

__kernel void backProjectionKernel(
	__read_only image2d_t sinogram, 
	float sinogramSpacingX,
	float sinogramSpacingY, 
	float sinogramOriginX,
	float sinogramOriginY, 
	__global float* backProjection, 
	float backProjSpacingX,
	float backProjSpacingY, 
	float backProjOriginX,
	float backProjOriginY, 
	float backProjSizeX,
	float backProjSizeY,
	int thetaIndexMax) {

	
	int x = get_global_id(0);
	int y = get_global_id(1);
	
	if (x >= backProjSizeX || y >= backProjSizeY)
		return;

	float xWorld = x * backProjSpacingX + backProjOriginX;
	float yWorld = y * backProjSpacingY + backProjOriginY;

	float sum = 0.0f;
			
	for (int thetaIndex = 0; thetaIndex < thetaIndexMax; thetaIndex++) {

		float theta = thetaIndex * sinogramSpacingY;
		float sinTheta = sin(radians(theta));
		float cosTheta = cos(radians(theta));
		float s = xWorld * cosTheta + yWorld * sinTheta;
		float sinoIdxX = 1 / sinogramSpacingX * (s - sinogramOriginX);
		float sinoIdxY = 1 / sinogramSpacingY * (theta - sinogramOriginY);
			
		sum += read_imagef(sinogram, sampler, (float2) (sinoIdxX, sinoIdxY)).x;
	}

	unsigned long index = y * backProjSizeX + x;
	backProjection[index] = sum;
		
	return;

}

