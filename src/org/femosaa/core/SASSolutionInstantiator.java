package org.femosaa.core;


import jmetal.core.Problem;
import jmetal.core.Solution;
import jmetal.core.SolutionSet;
import jmetal.core.Variable;

/**
 * A facade model in order to lower code coupling
 * 
 * TODO: make sure that all the instances of Solution in 
 * MOEAD_STM_SAS are instantiated from here.
 * 
 * @author tao
 *
 */
public interface SASSolutionInstantiator {

	public Solution getSolution (Problem problem);
	
	public Solution getSolution (Solution solution);
	
	public Solution getSolution (Problem problem, Variable[] variables);
	
	public Solution getSolution();
	
	public Solution getSolution(int objective_number);
	
	public SolutionSet fuzzilize (SolutionSet set);
	
	public SolutionSet fuzzilize(SolutionSet set, double increment);
	
	public Solution defuzzilize(int i, SolutionSet newPopulation, SolutionSet oldPopulation);
	
	public Solution defuzzilize(Solution s, SolutionSet oldPopulation);
	
	public void defuzzilizeAndRemove(Solution s, SolutionSet oldPopulation);
	
	public double[][] getLambda();
	
	public double[] getWeights();
	
	public double[][] getFixedBounds();
	
	public int record(SolutionSet set);
	
	public int record(Solution s);
}
