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
		for (int i = 0; i < populationSize; i++) {
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
	}

}
