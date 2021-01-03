package org.femosaa.core;

import java.io.IOException;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

import jmetal.core.Solution;
import jmetal.core.SolutionSet;
import jmetal.util.JMException;

public abstract class SASAlgorithmAdaptor {
	public static String valid1 = "", valid2 = "", invalid1 = "", invalid2 = "";
	//protected static double MUTATION_RATE = 0.1;
	//protected static double CROSSOVER_RATE = 0.9;
	private static final boolean PRINT_SOLUTIONS = true; /*important*/
	private static final boolean PRINT_INVALID_SOLUTIONS = false;
	private static final boolean LOG_SOLUTIONS = true; /*important*/
	private static final boolean LOG_NON_DOMINATED_SOLUTIONS = false;
	
	public static boolean isLogToD = false;
	public static double d = 14.5;
	private static final boolean LOG_SOLUTIONS_VALUES = false; /*important*/
	// This can be changed within SSASE
	public static boolean isPreserveInvalidSolution = false;
	// This can be changed within SSASE
	public static boolean isSeedSolution = false;
	// This can be changed within SSASE
	public static boolean isLogTheEvalNeededToRemiveNonSeed = false;
	// This is actually number of function evaluation
	public static int logGenerationOfObjectiveValue = 100;//500;//5000; // <=0 means disabled.
	// This is to control if fuzzy requirement is enable
	public static boolean isFuzzy = false; 
	public static boolean isWeightedSumNormalized = true; 
	// This is to control if all solutions that have ever been found needs to be recalled
	public static boolean isLogDiscardedSolutions = false; 
	
	public static boolean logMeasurementOfObjectiveValue = false;
	
	public static boolean logNumberOfNoDominated = false;

	// For SGA, this means the total time
	public static long seed_time = -1;//test only for seeding, -1 means disabled
	public Solution execute(SASSolutionInstantiator factory, int[][] vars,
			int numberOfObjectives_, int numberOfConstraints_)
			throws JMException, SecurityException, IOException,
			ClassNotFoundException {
		SolutionSet pareto_front = findParetoFront(factory, vars,
				numberOfObjectives_, numberOfConstraints_);
		
		//logDependencyAfterEvolution(getAllFoundSolutions());
		
		SolutionSet result = correctDependencyAfterEvolution(pareto_front);
		// Means no valid solutions.
		if(result.size() == 0) {
			for (int i = 0; i < pareto_front.size(); i++) {
				// Make sure the solution does not violate dependency.
				// The dependencyMap is temporarily reset within the function.
				 ((SASSolution)pareto_front.get(i)).correctDependency();
				 double[] f = ((SASSolution)pareto_front.get(i)).getObjectiveValuesFromIndexValue();
					
					for (int k = 0; k < f.length ; k ++) {
						((SASSolution)pareto_front.get(i)).setObjective(k, f[k]);
					}
			}
			
		} else {
			pareto_front = result;
		}
		//pareto_front = doRanking(pareto_front);
		pareto_front = filterRequirementsAfterEvolution(pareto_front);
		if(PRINT_SOLUTIONS) {
			printSolutions(pareto_front, numberOfObjectives_);
		}		
		if(LOG_SOLUTIONS) {
			org.femosaa.util.Logger.logSolutionSet(pareto_front, "SolutionSet.rtf");
			//org.femosaa.util.Logger.logPercentageOfMarkedSolution(pareto_front, "HowManyFromSeeds.rtf");
		}	
		if(LOG_SOLUTIONS_VALUES) {
			org.femosaa.util.Logger.logSolutionSetValues(pareto_front, "SolutionSetValue.rtf");
		}	
		pareto_front = doRanking(pareto_front);
		if(PRINT_SOLUTIONS) {
			printSolutions(pareto_front, numberOfObjectives_);
		}
		if(LOG_NON_DOMINATED_SOLUTIONS) {
			org.femosaa.util.Logger.logSolutionSet(pareto_front, "NonDominatedSolutionSet.rtf");
		}
		printParetoFront(pareto_front);
		if(PRINT_INVALID_SOLUTIONS) {
			printSolutionsIsFromInvalidity(pareto_front, numberOfObjectives_);
		}
		
		//System.out.print("Pareto front size after dependency and requirement check: " + pareto_front.size() + "\n");
		return findSoleSolutionAfterEvolution(pareto_front);
	}

	protected abstract SolutionSet findParetoFront(
			SASSolutionInstantiator factory, int[][] vars,
			int numberOfObjectives_, int numberOfConstraints_)
			throws JMException, SecurityException, IOException,
			ClassNotFoundException;

	protected abstract ApproachType getName();
	
	protected abstract SolutionSet doRanking(SolutionSet pf);
	
	

	/**
	 * This can return knee point or just randomly selected one.
	 * @param pareto_front
	 * @return
	 */
	protected abstract Solution findSoleSolutionAfterEvolution(
			SolutionSet pareto_front);

	
	/**
	 * Used mainly by the existing approach.
	 * @param pareto_front
	 * @return
	 */
	protected SolutionSet correctDependencyAfterEvolution(
			SolutionSet pareto_front) {
		return pareto_front;
	}
	
	protected void printParetoFront(SolutionSet pareto_front) {
    }
	
	
	protected void printSolutions(SolutionSet pareto_front, int numberOfObjectives_) {
		System.out.print("Total number : " + pareto_front.size() + "\n");
		System.out.print("=========== solutions' objective values ===========\n");
		Iterator itr = pareto_front.iterator();
		while(itr.hasNext()) {
			Solution s = (Solution)itr.next();
			String str = "(";
			for(int i = 0; i < numberOfObjectives_; i++) {
				str = (i +1 == numberOfObjectives_)? str + s.getObjective(i) : 
					str + s.getObjective(i) + ",";
			}
			System.out.print(str+ ")\n");
		}
    }
	
	protected void printSolutionsIsFromInvalidity(SolutionSet pareto_front, int numberOfObjectives_) {
		System.out.print("Check invalidity, total number : " + pareto_front.size() + "\n");
		int count = 0;
		Iterator itr = pareto_front.iterator();
		String valid1 = "", valid2 = "";
		String invalid1 = "", invalid2 = "";
		while(itr.hasNext()) {
			SASSolution s = (SASSolution)itr.next();
			if (s.isFromInValid) {
				count++;
				invalid1 += s.getObjective(0) + "\n";
				invalid2 += s.getObjective(1) + "\n";
			} else {
				valid1 += s.getObjective(0) + "\n";
				valid2 += s.getObjective(1) + "\n";
			}
		}
		
		SASAlgorithmAdaptor.invalid1 += invalid1;
		SASAlgorithmAdaptor.invalid2 += invalid2;
		SASAlgorithmAdaptor.valid1 += valid1;
		SASAlgorithmAdaptor.valid2 += valid2;
		
		System.out.print(count + " solutions are the descendants of invalid ones\n");
		System.out.print("Valid ones 1:\n");
		System.out.print(valid1);
		System.out.print("Valid ones 2:\n");
		System.out.print(valid2);
		System.out.print("Invalid ones 1:\n");
		System.out.print(invalid1);
		System.out.print("Invalid ones 2:\n");
		System.out.print(invalid2);
    }
	

//	protected void logDependencyAfterEvolution(
//			SolutionSet pareto_front_without_ranking) {
//	}
	/**
	 * If the requirements (constraints) are considered outside the evolution,
	 * then it should be processed here.
	 * @param pareto_front
	 * @return
	 */
	protected SolutionSet filterRequirementsAfterEvolution(
			SolutionSet pareto_front) {
		return pareto_front;
	}

	protected enum ApproachType {
		MOEAD_STM_D_K, NSGAII, GP, MIP, IBEA, SGA
	}
}
