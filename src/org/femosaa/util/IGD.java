package org.femosaa.util;

public class IGD extends Indicator{

	
	public double invertedGenerationalDistance( double[][] front,  double[][] referenceFront) {
	    double sum = 0.0;
	    for (int i = 0 ; i < referenceFront.length; i++) {
	      sum += Math.pow(distanceToClosestPoint(referenceFront[i],
	          front), 2.0);
	    }

	    sum = Math.pow(sum, 1.0 / 2.0);

	    return sum / referenceFront.length;
	  }
	
	private double distanceToClosestPoint(double[] point, double[][] front) {
	 

	    double minDistance = distance(point, front[0]);

	    for (int i = 0; i < front.length; i++) {
	      double aux = distance(point, front[i]);//distance.compute(point, front.getPoint(i));
	      if (aux < minDistance) {
	        minDistance = aux;
	      }
	    }

	    return minDistance;
	  }
	
	//EuclideanDistance
	private double distance(double[] a, double[] b){
		double r = 0.0;
		for (int i = 0; i < a.length; i++) {
			r += Math.pow((a[i] + b[i]), 2);
		}
		
		return Math.sqrt(r);
	}
}
