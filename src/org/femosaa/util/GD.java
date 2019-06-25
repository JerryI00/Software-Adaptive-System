package org.femosaa.util;

public class  GD extends Indicator{

	private double pow = 2.0;

	public double generationalDistance(double[][] front,
			double[][] referenceFront) {
		
		front = nondominating(front);
		//front = removeRedundant(front);
		//System.out.print("Non-dominating start: " + front.length + "\n");
		double sum = 0.0;
		for (int i = 0; i < front.length; i++) {
			sum += Math.pow(distanceToClosestPoint(front[i], referenceFront),
					pow);
		}

		sum = Math.pow(sum, 1.0 / pow);

		return sum / front.length;
	}

	private double distanceToClosestPoint(double[] point, double[][] front) {

		double minDistance = distance(point, front[0]);

		for (int i = 0; i < front.length; i++) {
			double aux = distance(point, front[i]);// distance.compute(point,
													// front.getPoint(i));
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
