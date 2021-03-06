/*
 * Copyright (C) 2014 Marcel Pohlmann
 * CONRAD is developed as an Open Source project under the GNU General Public License (GPL).
*/
package edu.stanford.rsl.tutorial.fourierConsistency.wedgefilter;

import java.util.ArrayList;

import edu.stanford.rsl.conrad.data.numeric.Grid2D;

/**
 * This class creates a binary representation of the double wedge region of the Fourier transform
 * of a fan beam sinogram measured with an equal spaced detector. The output Grid2D can be used as 
 * band-stop filter. The implementation to build the filter follows the mathematics described in 
 * the paper:
 * "Fourier properties of the fan-beam sinogram", 
 * Mazin, Samuel R., and Norbert J. Pelc, Medical physics 37.4 (2010): 1674-1680.
 *  
 * @author Marcel Pohlmann
 * 
 */

public class DoubleWedgeFilterFanES extends DoubleWedgeFilter{

	private double _L; // source-isocenter distance
	private double _D; // detector isocenter distance

	/**
	 * Constructor for the fan beam double wedge filter (equal spaced detector).
	 * 
	 * @param size Size of the 2D-Fourier transform
	 * @param spacing Spacing of the 2D-Fourier transform
	 * @param origin Origin of the 2D-Fourier transform
	 * @param rp Maximum object extend
	 * @param _L source-isocenter distance L
	 * @param _D detector-isocenter distance D
	 */
	public DoubleWedgeFilterFanES(int[] size, double[] spacing, double[] origin, double rp, double _L, double _D) {
		super(size, spacing, origin, rp);

		this._L = _L;
		this._D = _D;

		this.update();
	}

	public void setParameterL(double _L) {
		this._L = _L;

		this.update();
	}

	public void setParameterD(double _D) {
		this._D = _D;

		this.update();
	}

	public double getParameterL() {
		return this._L;
	}

	public double getParameterD() {
		return this._D;
	}

	protected void update(){
		double k, zeta;
		
		// computes the fan beam double wedge filter (equal spaced detector) following the parametrization
		// shown in Eq. (15) in the paper 
		for(int i = 0; i < size[1]; ++i) {
			for(int j = 0; j < size [0]; ++j) {
				k = this.indexToPhysical(j, i)[1];
				zeta = this.indexToPhysical(j, i)[0];

				// the factor 0.5 is added on the right hand side of the inequality, because the tutorial 
				// forwardprojector class locates the detector into the isocenter and we assume that _L = _D
				if (Math.abs(k/(k-(zeta*(this._L + this._D)))) > (0.5*this.rp/this._L)) {
					this.setAtIndex(j, i, 0.0f);
				}else{
					this.setAtIndex(j, i, 1.0f);
				}
			}
		}
	}

	/**
	 * Apply an erosion in the direction of q on the double wedge filter with a kernel 
	 * size that is defined by the parameter
	 * 
	 * @param diameter Kernel width
	 */
	public void erode(int diameter){
		Grid2D erodedVersion = new Grid2D(this.getWidth(), this.getHeight());

		for(int i = 0; i < erodedVersion.getHeight(); ++i){
			for(int j = 0; j < erodedVersion.getWidth(); ++j){
				erodedVersion.setAtIndex(j, i, 1.0f);
			}
		}

		for(int i = 0; i < this.getHeight(); ++i){
			for(int j = 0; j < this.getWidth(); ++j){
				ArrayList<Float> myList = new ArrayList<Float>();
				for(int k = 1; k < diameter/2; ++k){
					if(j-k < 0 || j+k >= this.getWidth()){

					}else{
						myList.add(this.getAtIndex(j-k, i));
						myList.add(this.getAtIndex(j+k, i));
					}
					float checker = 0;
					for(Float flag:myList){
						checker += (flag != null ? flag:Float.NaN);
					}
					if(checker == 0.0){
						erodedVersion.setAtIndex(j, i, 0.0f);
					}else{
						erodedVersion.setAtIndex(j, i, 1.0f);
					}
				}
			}
		}
		for(int i = 0; i < this.getHeight(); ++i){
			for(int j = 0; j < this.getWidth(); ++j){
				this.setAtIndex(j, i, erodedVersion.getAtIndex(j, i));
			}
		}
	}

}
