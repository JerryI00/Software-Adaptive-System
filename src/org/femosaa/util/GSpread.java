package org.femosaa.util;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

public class GSpread extends Indicator {

	
	/**
	   *  Calculates the generalized spread metric. Given the 
	   *  pareto front, the true pareto front as <code>double []</code>
	   *  and the number of objectives, the method return the value for the
	   *  metric.
	   *  @param front The front.
	   *  @param referenceFront The reference pareto front. {1.0,0.0,0.0} {0.0,1.0,0.0} {0.0,0.0,1.0}
	   *  @return the value of the generalized spread metric
	   **/
	  @SuppressWarnings("unchecked")
	public double generalizedSpread(double[][] front, double[][] extremeValues) {
		  
		front = nondominating(front);
		  
	    int numberOfObjectives = front[0].length ;

	   // double[] extremeValues = new double[numberOfObjectives] ;
	    for (int i = 0; i < numberOfObjectives; i++) {
	      /*referenceFront.sort(new PointDimensionComparator(i));
	      Point newPoint = new ArrayPoint(numberOfObjectives) ;
	      for (int j = 0 ; j < numberOfObjectives; j++) {
	        newPoint.setValue(j,
	            referenceFront.getPoint(referenceFront.getNumberOfPoints()-1).getValue(j));
	      }
	      extremeValues[i] = newPoint ;*/
	    }

	    int numberOfPoints = front.length;

	    //front.sort(new LexicographicalPointComparator());
	    List<Double[]> list = new ArrayList<Double[]>();
	    for (double[] d : front) {
	    	Double[] newD = new Double[d.length];
	    	 for (int i = 0; i < d.length; i++) { 
	    		 newD[i] = d[i];
	    	 }
	    	list.add(newD);
	    }
	    
	    Collections.sort(list, new Comparator() {

			@Override
			public int compare(Object obj1, Object obj2) {
				
				Double[] pointOne = (Double[])obj1; 
				Double[] pointTwo = (Double[])obj2; 
				 // Determine the first i such as pointOne[i] != pointTwo[i];
			    int index = 0;
			    while ((index < pointOne.length)
			        && (index < pointTwo.length)
			        && pointOne[index] == pointTwo[index]) {
			      index++;
			    }

			    int result = 0 ;
			    if ((index >= pointOne.length) || (index >= pointTwo.length)) {
			      result = 0;
			    } else if (pointOne.length < pointTwo.length) {
			      result = -1;
			    } else if (pointOne[index] > pointTwo[index]) {
			      result = 1;
			    }
			    return result ;
			}
	    	
	    });
	    
	    
	    
	    for (int i = 0; i < list.size(); i++) {
	    	
	    	 for (int j = 0; j < list.get(i).length; j++) { 
	    		 front[i][j] = list.get(i)[j];
	    	 }
	    }
	    	
	    
	   
	    
	   
	    if ( distance(front[0], front[front.length-1]) == 0.0) {
	      return 1.0;
	    } else {
	      double dmean = 0.0;

	      for (int i = 0 ; i < numberOfPoints; i++) {
	        dmean += distanceToNearestPoint(front[i], front);
	      }

	      dmean = dmean / (numberOfPoints);

	      double dExtrems = 0.0;
	      for (int i = 0 ; i < extremeValues.length; i++) {
	        dExtrems += distanceToClosestPoint(extremeValues[i], front);
	      }

	      double mean = 0.0;
	      for (int i = 0; i < numberOfPoints; i++) {
	        mean += Math.abs(distanceToNearestPoint(front[i], front) -
	            dmean);
	      }

	      return (dExtrems + mean) / (dExtrems + (numberOfPoints*dmean));
	    }
	  }
	  

	  
	  
	  public double distanceToNearestPoint(double[] point, double[][] front) {
		
		    double minDistance = Double.MAX_VALUE;

		    for (int i = 0; i < front.length; i++) {
		      double aux = distance(point, front[i]);
		      if ((aux < minDistance) && (aux > 0.0)) {
		        minDistance = aux;
		      }
		    }

		    return minDistance;
		  }
	  
	  public double distanceToClosestPoint(double[] point, double[][] front) {
		  

		    double minDistance = distance(point, front[0]);

		    for (int i = 1; i < front.length; i++) {
		      double aux = distance(point, front[i]);
		      if (aux < minDistance) {
		        minDistance = aux;
		      }
		    }

		    return minDistance;
		  }
	  
		// EuclideanDistance
		private double distance(double[] a, double[] b) {
			double r = 0.0;
			for (int i = 0; i < a.length; i++) {
				r += Math.pow((a[i] + b[i]), 2);
			}

			return Math.sqrt(r);
		}
}
