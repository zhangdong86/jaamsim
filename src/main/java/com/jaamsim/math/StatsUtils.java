/*
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package com.jaamsim.math;

// Apache-Commons-Math is required
// http://commons.apache.org/proper/commons-math/
import org.apache.commons.math3.distribution.TDistribution;


/**
 * Stores functions relating to statistics, e.g t-coefficients.
 * 
 *  @author Michael Bergman
 */
public class StatsUtils {

	/**
	 * To calculate the t-statistic for a given p-value and number of items
	 * For a 2 tail p value, t-test
	 */
	public static double tStat2Tail(double p, int n) {
		assert(0 <= p);
		assert(p < 1);
		assert(1 < n);
		
		double out = 0f;
		double degFree = n-1.0;
			
		p /= 2;
		// T distribution
		TDistribution T = new TDistribution(null, degFree, 1e-9);
		out = -T.inverseCumulativeProbability(p);	// as the distribution is symmetrical
		return out;
	}
		
	
	/**
	 * To calculate the t-statistic for a given p-value and number of items
	 * For a 1 tail value
	 */
	public static double tStat1Tail(double p, int n) {
		assert(0 < p);
		assert(p < 1);
		assert(1 < n);
		
		double out = 0f;
		double degFree = n-1.0;
			
		// T distribution
		TDistribution T = new TDistribution(null, degFree, 1e-9);
		out = T.inverseCumulativeProbability(p);	
		return out;
	}

	
	// TODO: maybe change parameter from array to (array)list
	/**
	 * Calculates the autocovariance for an array, with all elements having the
	 * same weight.
	 * 
	 * @param X the array of the sequence
	 * @param n the distance between elements of the sequence
	 * @return the autocovariance
	 */
	public static double AutoCovariance(double[] X, int n) {
		
		assert (n > 0);
		/*TODO: length check on X to make sure there are enough entries
		 *		for the given m
		*/
		
		double muX0 = 0d;
		double muXn = 0d;
		double autocov = 0.0d;
		
		// calculate the means
		for(int i = 0; i < X.length - n ; i++) {
			muX0 += X[i];
			muXn += X[i+n];
		}
		
		muX0 /= X.length - n;
		muXn /= X.length - n;
		
		
		for(int i = 0; i < X.length - n ; i++) {
			autocov += (X[i] - muX0)*(X[i+n] - muXn);
		}
		
		autocov /= X.length - n;
		
		return autocov;
	}
 }
