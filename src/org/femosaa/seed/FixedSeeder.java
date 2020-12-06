package org.femosaa.seed;

import java.util.List;

import jmetal.core.Operator;
import jmetal.core.Problem;
import jmetal.core.Solution;
import jmetal.core.SolutionSet;
import jmetal.util.JMException;

import org.femosaa.core.SASSolutionInstantiator;

public class FixedSeeder extends Seeder{
	
	private FixedSeeder(Operator mutationOperator) {
		super(mutationOperator);
		// TODO Auto-generated constructor stub
	}

	private static FixedSeeder seeder = null;
	private List<String> seeds = null;
	
	public static FixedSeeder getInstance(){
		if(seeder == null) {
			seeder = new FixedSeeder(null);
		}
		return seeder;
	}
	
	public void setSeeds(List<String> seeds){
		this.seeds = seeds;
	}
	
	public void seeding(SolutionSet population,
			SASSolutionInstantiator factory, Problem problem_,
			int populationSize) throws JMException {
		Solution newSolution;
		int l = populationSize < 20? 1 : populationSize/20;//20 is for s-vs-m; 1 is used for req-vs-mo
		for (int i = 0; i < l; i++) {
			newSolution = factory.getSolution(problem_);
			String[] d = seeds.get(i).split(":");
			for (int j = 0; j < newSolution.getDecisionVariables().length; j++) {
				newSolution.getDecisionVariables()[j].setValue(Double.parseDouble(d[j]));
			}
			problem_.evaluate(newSolution);
			problem_.evaluateConstraints(newSolution);
			population.add(newSolution);
		
			//System.out.print("Solution " + newSolution.getObjective(0) + "-" + newSolution.getObjective(1)+"\n");
		} 
		
		if(population.size() < populationSize) {
			int s = population.size();
			for (int i = 0; i < populationSize - s; i++) {
				newSolution = factory.getSolution(problem_);
				problem_.evaluate(newSolution);
				problem_.evaluateConstraints(newSolution);
				population.add(newSolution);
			}
		}
	}

}
