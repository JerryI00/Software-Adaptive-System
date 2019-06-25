package org.femosaa.util;

import java.util.HashSet;
import java.util.Set;

import jmetal.core.Solution;
import jmetal.core.SolutionSet;
import jmetal.util.Ranking;

public abstract class Indicator {

	
	protected double[][] removeRedundant(double[][] front){
		
		Set<String> set = new HashSet<String>();
		
        for (int i = 0; i < front.length;i++){
			
			String s = "";
			for (int j = 0; j < front[i].length;j++){
				s += front[i][j] + " ";
			}
		
		
			set.add(s);
		}
		
		double[][] result = new double[set.size()][];
		int i = 0;
		for (String s : set) {
			double[] d = new double[3];
			d[0] = Double.valueOf(s.split(" ")[0]);
			d[1] = Double.valueOf(s.split(" ")[1]);
			d[2] = Double.valueOf(s.split(" ")[2]);
			result[i] = d;
			i++;
		}
		
		return result;
	}
	
	protected double[][] nondominating(double[][] front){
		SolutionSet set = new SolutionSet();
		for (int i = 0; i < front.length;i++){
			
			Solution s = new Solution(front[i].length);
			for (int j = 0; j < front[i].length;j++){
				s.setObjective(j, front[i][j]);
			}
		
			set.add(s);
		}
		
		Ranking ranking = new Ranking(set);
		set =  ranking.getSubfront(0);
		//System.out.print(set.size());
		double[][] result = new double[set.size()][];
		
		for (int i = 0; i < set.size();i++){
			
			double[] d = new double[set.get(i).numberOfObjectives()];
			for (int j = 0; j < set.get(i).numberOfObjectives();j++){
				d[j] = set.get(i).getObjective(j);
			}
			result[i] = d;
		}
		
		return result;
	}
}
