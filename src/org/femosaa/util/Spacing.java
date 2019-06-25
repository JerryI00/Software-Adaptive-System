package org.femosaa.util;

public class Spacing extends Indicator{

	
	public static void main(String[] a) {
		HV hv = new HV();
		GD gd = new GD();
		Spacing sp = new Spacing();
		
		double[][] test  = new double[1][];
		test[0] = new double[]{0.0,0.0};
		
		double[][] d  = new double[4][];
		d[0] = new double[]{0.3,0.3,0.3};
		d[1] = new double[]{0.3,0.3,0.3};
		d[2] = new double[]{0.3,0.3,0.3};
		d[3] = new double[]{0.3,0.3,0.3};
		//d[2] = new double[]{0.5,0.5,0.5};
		//d[2] = new double[]{0.2,0.4,0.2};
		//d[3] = new double[]{0.6,0.6,0.6};
		//d[4] = new double[]{0.6,0.6,0.6};
		
		//d[4] = new double[]{0.9,0.6,0.9};
		//System.out.print(hv.hypervolume(d));
		//System.out.print(gd.generationalDistance(d,test));
		System.out.print(sp.spacing(d));
		//0.17320508075688773
	}
	
	
	public double spacing(double[][] approximationSet) {
		
		
		approximationSet = nondominating(approximationSet);
		//approximationSet = removeRedundant(approximationSet);
		
		if (approximationSet.length< 2) {
			return 1.0;
		}
		/*System.out.print("Non-dominating start: " + approximationSet.length + "\n");
		for (double[] dd : approximationSet) {
			System.out.print(dd[0] + " " + dd[1] + " " + dd[2] + "\n");
		}
		System.out.print("Non-dominating end: " + approximationSet.length + "\n");*/
		double[] d = new double[approximationSet.length];

		for (int i = 0; i < approximationSet.length; i++) {
			double min = Double.POSITIVE_INFINITY;
			double[] solutionI = approximationSet[i];
			
			
			
			
			for (int j = 0; j < approximationSet.length; j++) {
				if (i != j) {
					double[] solutionJ = approximationSet[j];
					
					
					
					min = Math.min(min, distance(solutionI, solutionJ, 1.0));
				}
			}

			d[i] = min;
		}

		double dbar = 0.0;
		
		for (double v : d) {
			dbar += v;
		}
		
		//dbar = StatUtils.sum(d) / approximationSet.length;
		dbar = dbar / approximationSet.length;
		double sum = 0.0;
		
		for (int i = 0; i < approximationSet.length; i++) {
			
			
			sum += Math.pow(d[i] - dbar, 2.0);
		}

		return Math.sqrt(sum / (approximationSet.length - 1));
	}
	
	private double distance(double[] a, double[] b,
			double power) {
		double distance = 0.0;

		for (int i = 0; i < a.length; i++) {
			distance += Math.pow(Math.abs(a[i] - 
					b[i]), power);
		}

		return Math.pow(distance, 1.0 / power);
	}
}
