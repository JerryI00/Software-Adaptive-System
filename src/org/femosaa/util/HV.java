package org.femosaa.util;

public class HV {

	private boolean dominates(double point1[], double point2[], int noObjectives) {
	    int i;
	    int betterInAnyObjective;

	    betterInAnyObjective = 0;
	    for (i = 0; i < noObjectives && point1[i] >= point2[i]; i++) {
	      if (point1[i] > point2[i]) {
	        betterInAnyObjective = 1;
	      }
	    }

	    return ((i >= noObjectives) && (betterInAnyObjective > 0));
	  }

	  private void swap(double[][] front, int i, int j) {
	    double[] temp;

	    temp = front[i];
	    front[i] = front[j];
	    front[j] = temp;
	  }

	  /* all nondominated points regarding the first 'noObjectives' dimensions
	  are collected; the points referenced by 'front[0..noPoints-1]' are
	  considered; 'front' is resorted, such that 'front[0..n-1]' contains
	  the nondominated points; n is returned */
	  private int filterNondominatedSet(double[][] front, int noPoints, int noObjectives) {
	    int i, j;
	    int n;

	    n = noPoints;
	    i = 0;
	    while (i < n) {
	      j = i + 1;
	      while (j < n) {
	        if (dominates(front[i], front[j], noObjectives)) {
	  /* remove point 'j' */
	          n--;
	          swap(front, j, n);
	        } else if (dominates(front[j], front[i], noObjectives)) {
		/* remove point 'i'; ensure that the point copied to index 'i'
		   is considered in the next outer loop (thus, decrement i) */
	          n--;
	          swap(front, i, n);
	          i--;
	          break;
	        } else {
	          j++;
	        }
	      }
	      i++;
	    }
	    return n;
	  }

	  /* calculate next value regarding dimension 'objective'; consider
	     points referenced in 'front[0..noPoints-1]' */
	  private double surfaceUnchangedTo(double[][] front, int noPoints, int objective) {
	    int i;
	    double minValue, value;

	    if (noPoints < 1) {
	      new RuntimeException("run-time error");
	    }

	    minValue = front[0][objective];
	    for (i = 1; i < noPoints; i++) {
	      value = front[i][objective];
	      if (value < minValue) {
	        minValue = value;
	      }
	    }
	    return minValue;
	  }

	  /* remove all points which have a value <= 'threshold' regarding the
	     dimension 'objective'; the points referenced by
	     'front[0..noPoints-1]' are considered; 'front' is resorted, such that
	     'front[0..n-1]' contains the remaining points; 'n' is returned */
	  private int reduceNondominatedSet(double[][] front, int noPoints, int objective,
	                                    double threshold) {
	    int n;
	    int i;

	    n = noPoints;
	    for (i = 0; i < n; i++) {
	      if (front[i][objective] <= threshold) {
	        n--;
	        swap(front, i, n);
	      }
	    }

	    return n;
	  }

	  private double calculateHypervolume(double[][] front, int noPoints, int noObjectives) {
	    int n;
	    double volume, distance;

	    volume = 0;
	    distance = 0;
	    n = noPoints;
	    while (n > 0) {
	      int nonDominatedPoints;
	      double tempVolume, tempDistance;

	      nonDominatedPoints = filterNondominatedSet(front, n, noObjectives - 1);
	      if (noObjectives < 3) {
	        if (nonDominatedPoints < 1) {
	          new RuntimeException("run-time error");
	        }

	        tempVolume = front[0][0];
	      } else {
	        tempVolume = calculateHypervolume(front, nonDominatedPoints, noObjectives - 1);
	      }

	      tempDistance = surfaceUnchangedTo(front, n, noObjectives - 1);
	      volume += tempVolume * (tempDistance - distance);
	      distance = tempDistance;
	      n = reduceNondominatedSet(front, n, noObjectives - 1, distance);
	    }
	    return volume;
	  }

	  /**
	   * This implementation natraully handle non-dominated solution, thus
	   * no need to ensure the front is non-dominated.
	   * 
	   * Also, since the front is normalized, the reference point after invertion would be <0,0,0...>.
	   * 
	   * Returns the hypervolume value of a front of points
	   *
	   * @param front        The front
	   * @param referenceFront    The true pareto front
	   */
	  public double hypervolume(double[][] front) {

		 double[][] invertedFront;
	    invertedFront = getInvertedFront(front);

	    int numberOfObjectives = front[0].length;

	    // STEP4. The hypervolume (control is passed to the Java version of Zitzler code)
	    return this.calculateHypervolume(invertedFront,
	    		invertedFront.length, numberOfObjectives);
	  }
	  
	  public static double[][]  getInvertedFront(double[][] front) {
		 

		    int numberOfDimensions = front[0].length;
		    double[][] invertedFront = new  double[front.length][numberOfDimensions];

		    for (int i = 0; i < front.length; i++) {
		      for (int j = 0; j < numberOfDimensions; j++) {
		        if (front[i][j] <= 1.0
		            && front[i][j]  >= 0.0) {
		        	invertedFront[i][j] = 1.0 - front[i][j];
		        } else if (front[i][j] > 1.0) {
		        	invertedFront[i][j] = 0.0 ;
		        } else if (front[i][j] < 0.0) {
		        	invertedFront[i][j] = 1.0 ;
		        }
		      }
		    }
		    return invertedFront;
		  }
}
