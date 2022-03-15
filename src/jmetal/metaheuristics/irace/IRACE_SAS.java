//  NSGAII.java
//
//  Author:
//       Antonio J. Nebro <antonio@lcc.uma.es>
//       Juan J. Durillo <durillo@lcc.uma.es>
//
//  Copyright (c) 2011 Antonio J. Nebro, Juan J. Durillo
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU Lesser General Public License for more details.
// 
//  You should have received a copy of the GNU Lesser General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.

package jmetal.metaheuristics.irace;

import java.util.Arrays;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Random;
import java.util.Set;

import org.femosaa.core.EAConfigure;
import org.femosaa.core.SASAlgorithmAdaptor;
import org.femosaa.core.SASSolution;
import org.femosaa.core.SASSolutionInstantiator;
import org.femosaa.seed.Seeder;

import jmetal.core.*;
import jmetal.util.comparators.CrowdingComparator;
import jmetal.util.*;

/**
 * 
 * 
 *
 */
public class IRACE_SAS extends Algorithm {

	private SASSolutionInstantiator factory = null;
	
	private Set<String> full_set;

	boolean isNeigboring = true;
	private Seeder seeder = null;
	double[] weights = new double[0];
	
	// This to be used with single objective only and without any fuzzy setting
	double[][] fixed_bounds = null;
	// ideal point
	double[] z_;

	// nadir point
	double[] nz_;

	int populationSize_;

	SolutionSet population_;
	
	protected double t0 = 10.0;

	protected double tn = 1.0;


	/**
	 * Constructor
	 * 
	 * @param problem Problem to solve
	 */
	public IRACE_SAS(Problem problem) {
		super(problem);
	} // NSGAII

	/**
	 * Constructor
	 * 
	 * @param problem Problem to solve
	 */
	public IRACE_SAS(Problem problem, SASSolutionInstantiator factory) {
		super(problem);
		this.factory = factory;
	}

	/**
	 * Runs the NSGA-II algorithm.
	 * 
	 * @return a <code>SolutionSet</code> that is a set of non dominated solutions
	 *         as a result of the algorithm execution
	 * @throws JMException
	 */
	public SolutionSet execute() throws JMException, ClassNotFoundException {
		
		
		if(SASAlgorithmAdaptor.isLogSolutionsInFull) {
			full_set = new HashSet<String>();
		}
		

		if (factory == null) {
			throw new RuntimeException("No instance of SASSolutionInstantiator found!");
		}

		if (getInputParameter("seeder") != null) {
			seeder = (Seeder) getInputParameter("seeder");
		}
		if(getInputParameter("fixed_bounds") != null) {
			fixed_bounds = (double[][])getInputParameter("fixed_bounds");
		}

		int maxEvaluations = ((Integer) getInputParameter("maxEvaluations")).intValue();
		int evaluations;
		weights = (double[])getInputParameter("weights");


		z_ = new double[problem_.getNumberOfObjectives()];
		nz_ = new double[problem_.getNumberOfObjectives()];

		int populationSize =  ((Integer) getInputParameter("populationSize")).intValue();
		// Initialize the variables
		SolutionSet population = new SolutionSet();
		// SolutionSet nadir_population = new SolutionSet(populationSize);

		evaluations = 0;
		int measurement = 0;

		// Create the initial solutionSet
		Solution newSolution;
		if (seeder != null) {
			seeder.seeding(population, factory, problem_, populationSize);
			evaluations += populationSize;
			if(!SASAlgorithmAdaptor.isInvalidSolutionConsumeMeasurement) {
				for (int i = 0; i < populationSize; i++) {
					Solution s = population.get(i);
					if(s.getObjective(0) == Double.MAX_VALUE || s.getObjective(0) == Double.MAX_VALUE/100) {
						
					} else {
						measurement += factory.record(s);
					}
				}
				
				
			} else {
				measurement += factory.record(population);
			}
			fitnessAssignment(population.get(0));
		} else {
			// Create the initial solutionSet
			for (int i = 0; i < populationSize; i++) {
				newSolution = factory.getSolution(problem_);
				problem_.evaluate(newSolution);
				problem_.evaluateConstraints(newSolution);
				evaluations++;
				if (!SASAlgorithmAdaptor.isInvalidSolutionConsumeMeasurement) {

					if (newSolution.getObjective(0) == Double.MAX_VALUE
							|| newSolution.getObjective(0) == Double.MAX_VALUE / 100) {

					} else {
						measurement += factory.record(newSolution);
					}

				} else {
					measurement += factory.record(newSolution);
				}
				population.add(newSolution);
				fitnessAssignment(newSolution);
			} // for
		}

		if(SASAlgorithmAdaptor.isLogSolutionsInFull) {
			Iterator itr = population.iterator();
			while(itr.hasNext()) {
				full_set.add(convertFullInfo((Solution)itr.next()));
			}
			
		}
		
		initIdealPoint();
		initNadirPoint();

		SolutionSet old_population = new SolutionSet(populationSize);
		if (SASAlgorithmAdaptor.isFuzzy) {
			old_population = population;
			population = factory.fuzzilize(population);
		}

		for (int i = 0; i < population.size(); i++) {
			fitnessAssignment(population.get(i)); // assign fitness value to each solution
		}

		if (SASAlgorithmAdaptor.logGenerationOfObjectiveValue > 0) {
			if(SASAlgorithmAdaptor.isFuzzy) {
				org.femosaa.util.Logger.logSolutionSetWithGenerationAndFuzzyValue(population, old_population,
						"SolutionSetWithGen.rtf", 0);
			} else {
				org.femosaa.util.Logger.logSolutionSetWithGeneration(population, "SolutionSetWithGen.rtf", 
						0);
			}
		}

	
		
		int N_iter = 2 + (int)(Math.log(population.get(0).getDecisionVariables().length) / Math.log(2));
		
		double[][] para_p = new double[population.get(0).getDecisionVariables().length][];
		
		for (int i = 0; i < population.get(0).getDecisionVariables().length; i++) {
			int[] bounds = ((SASSolution) population.get(0)).getBounds(i);
			double[] temp = new double[bounds[1] + 1];
			
			for (int l = 0; l < temp.length; l++) {
				temp[l] = (double)(1.0 / temp.length);
				//System.out.print(temp[l] + "**\n");
			}
			para_p[i] = temp;
		}
	
		double te = 0.0;
		// Generations
		int j = 1;
		while (evaluations < maxEvaluations) {

			/**
			 * set a max limit to irace, in case it runs forerver.
			 */
			j++;
			if(j > N_iter) {
				j = N_iter;
			}
			

			if (EAConfigure.getInstance().measurement == measurement) {
				break;
			}
			
			
			int N_min = 5;
			int B_j = (EAConfigure.getInstance().measurement  - measurement) / (N_iter - j +1);
			int T_each = 1;
			int mu = 1;
			int N_j = B_j/(mu + T_each * Math.min(5, j));
			
			int N_new = N_j - population.size();
			
			N_new = N_new < N_min? N_min : N_new;
			
			System.out.print("evaluations: "+ evaluations + "\n");
			System.out.print("measurement: "+ measurement + "\n");
			System.out.print("B_j: " + B_j + "\n");
			System.out.print("N_j: " + N_j + "\n");
			System.out.print("N_new: " + N_new + "\n");
			SolutionSet offsprings = new SolutionSet();
			
			double[] parent_p = new double[population.size()];
			population.sort(new Comparator() {

				@Override
				public int compare(Object arg0, Object arg1) {
					Solution s1 = (Solution)arg0;
					Solution s2 = (Solution)arg1;
					
					if(s1.getFitness() < s2.getFitness()) {
						return -1;
					} else if(s1.getFitness() > s2.getFitness()) {
						return 1;
					}
					
					return 0;
				}
				
			});
			
			
			double N_elite = population.size();
			
			for (int i = 0; i < N_elite; i++) {
				parent_p [i] = (double) ((N_elite - (i+1) + 1) / ((N_elite * (N_elite + 1))/2));

				//System.out.print(parent_p[i] + "**\n");
			}
			
			for (int k = 0; k < N_new; k++) {
				if (EAConfigure.getInstance().measurement == measurement) {
					break;
				}
				// 1. Sample parent
				
				
				
				
				Solution parent = null;
				Solution child = null;
				double p = PseudoRandom.randDouble();
				double cumulativeProbability = 0.0;
				for (int i = parent_p.length - 1; i >= 0; i--) {
				    cumulativeProbability += parent_p[i];
				    if (p <= cumulativeProbability) {
				    	parent = population.get(i);
				    	break;
				    }
				}
				
				// 2. Sample parameter
				if(parent == null) {
					parent = population.get(0);
				}
				
				for (int i = 0; i < parent.getDecisionVariables().length; i++) {
					double[] temp = para_p[i];
					double v = parent.getDecisionVariables()[i].getValue();
					double total = 0.0;
					for (int l = 0; l < temp.length; l++) {
						
						if(l == v) {
							temp[l] = (double) (temp[l] * (1 - ((j-1)/N_iter)) + ((j-1)/N_iter));
						} else {
							temp[l] = (double) (temp[l] * (1 - ((j-1)/N_iter)));
						}
						//System.out.print(temp[l] + "**\n");
						total += temp[l];
					}
					
					for (int l = 0; l < temp.length; l++) {
						temp[l] = temp[l] / total;
					}
				}
				
				child = factory.getSolution(parent);
				
				for (int v = 0; v < child.getDecisionVariables().length; v++) {
					
					double[] prob = para_p[v];
					Arrays.sort(prob);
					
					p = PseudoRandom.randDouble();
					cumulativeProbability = 0.0;
					for (int i = 0; i < prob.length; i++) {
					    cumulativeProbability += prob[i];
					    if (p <= cumulativeProbability) {
					    	child.getDecisionVariables()[v].setValue(i);
					    	break;
					    }
					}
				}
				
				problem_.evaluate(child);
				problem_.evaluateConstraints(child);
				
				updateReference(child);
				updateNadirPoint(child);

				
				
				
				offsprings.add(child);
				measurement += factory.record(child);
				if(SASAlgorithmAdaptor.isLogSolutionsInFull) {
						full_set.add(convertFullInfo((Solution)child));
				}
				evaluations++;
			}
			
			
			
			
			
			
			// 3. Comnbination
			
			SolutionSet union = population.union(offsprings);
			double max = Double.NEGATIVE_INFINITY;
			for (int i = 0; i < union.size(); i++) {
				fitnessAssignment(union.get(i));
				if(union.get(i).getFitness() > max) {
					max = union.get(i).getFitness();
				}
			}
			
			// 4. Race

			population.clear();
			
			if(union.size() == 1) {
				population.add(union.get(0));
			} else {
				for (int i = 0; i < union.size(); i++) {
					if(union.get(i).getFitness() != max) {
						population.add(union.get(i));
					}
				}
			}
			
			if(population.size() == 0) {
				population = population.union(union);
			}
			
			
			
			

			
			
			if (SASAlgorithmAdaptor.logGenerationOfObjectiveValue > 0 && evaluations%SASAlgorithmAdaptor.logGenerationOfObjectiveValue == 0) {
				if(SASAlgorithmAdaptor.isFuzzy) {
					org.femosaa.util.Logger.logSolutionSetWithGenerationAndFuzzyValue(population, old_population,
							"SolutionSetWithGen.rtf", evaluations);
				} else {
					org.femosaa.util.Logger.logSolutionSetWithGeneration(population, "SolutionSetWithGen.rtf", 
							evaluations);
				}
			}
			
			//System.out.print("Measurement: " + measurement + "\n");

			if (SASAlgorithmAdaptor.logMeasurementOfObjectiveValue) {
				if (SASAlgorithmAdaptor.isFuzzy) {
					org.femosaa.util.Logger.logSolutionSetWithGeneration(old_population,
							"SolutionSetWithMeasurement.rtf", measurement);
				} else {
					org.femosaa.util.Logger.logSolutionSetWithGeneration(population, "SolutionSetWithMeasurement.rtf",
							measurement);
				}

			}
			
			

		} // while
		
		if (SASAlgorithmAdaptor.logGenerationOfObjectiveValue > 0 && evaluations%SASAlgorithmAdaptor.logGenerationOfObjectiveValue == 0) {
			if(SASAlgorithmAdaptor.isFuzzy) {
				org.femosaa.util.Logger.logSolutionSetWithGenerationAndFuzzyValue(population, old_population,
						"SolutionSetWithGen.rtf", evaluations);
			} else {
				org.femosaa.util.Logger.logSolutionSetWithGeneration(population, "SolutionSetWithGen.rtf", 
						evaluations);
			}
		}

		if (SASAlgorithmAdaptor.isFuzzy) {
			population = old_population;
		}
		// Return as output parameter the required evaluations
		// setOutputParameter("evaluations", requiredEvaluations);
		if (SASAlgorithmAdaptor.isLogToD) {
			System.out.print("Minimum evaluation " + te + "\n");
			org.femosaa.util.Logger.logFirstTod(te, "FirstToD.rtf");
		}
		
		if(SASAlgorithmAdaptor.isLogSolutionsInFull) {
			org.femosaa.util.Logger.logSolutionFull(full_set, "FullSolution.rtf");
		}

		return population;
	} // execute

	public SolutionSet doRanking(SolutionSet population) {
		SolutionSet set = new SolutionSet(1);
		set.add(population.get(0));

		return set;
	}

	/**
	 * This is used to assign fitness value to a solution, according to weighted sum
	 * strategy.
	 * 
	 * @param cur_solution
	 */
	public void fitnessAssignment(Solution cur_solution) {
		double cur_fitness = 0.0;
		// double weight = 1.0 / (double) problem_.getNumberOfObjectives();

		for (int i = 0; i < problem_.getNumberOfObjectives(); i++) {

			if (SASAlgorithmAdaptor.isWeightedSumNormalized) {

				if (fixed_bounds != null) {

					if (fixed_bounds[0][0] == 0 && fixed_bounds[0][1] == 0 && fixed_bounds[1][0] == 0
							&& fixed_bounds[1][1] == 0) {
						if(cur_solution.getObjective(i) < 0) {
							cur_fitness += weights[i] * (cur_solution.getObjective(i)/(cur_solution.getObjective(i)-1));

						} else {
							cur_fitness += weights[i] * (cur_solution.getObjective(i)/(cur_solution.getObjective(i)+1));
						}
					} else if (fixed_bounds[0][0] == -1 && fixed_bounds[0][1] == -1 && fixed_bounds[1][0] == -1
							&& fixed_bounds[1][1] == -1) {
						cur_fitness += weights[i] * cur_solution.getObjective(i);
					} else {
						if (cur_solution.getObjective(i) == Double.MAX_VALUE / 100) {
							cur_fitness += weights[i] * 1.0;
							// System.out.print(cur_fitness + " Find one fitness with MAX_VALUE!\n");
						} else {

							cur_fitness += weights[i]
									* (nz_[i] != z_[i] ? ((cur_solution.getObjective(i) - z_[i]) / (nz_[i] - z_[i]))
											: ((cur_solution.getObjective(i) - z_[i]) / (nz_[i])));
						}
					}

				} else {

					if (cur_solution.getObjective(i) == Double.MAX_VALUE / 100) {
						cur_fitness += weights[i] * 1.0;
						// System.out.print(cur_fitness + " Find one fitness with MAX_VALUE!\n");
					} else {

						cur_fitness += weights[i]
								* (nz_[i] != z_[i] ? ((cur_solution.getObjective(i) - z_[i]) / (nz_[i] - z_[i]))
										: ((cur_solution.getObjective(i) - z_[i]) / (nz_[i])));
					}
				}
			} else {
				cur_fitness += weights[i] * cur_solution.getObjective(i);
			}

		}

		if(Double.isNaN(cur_fitness)) {
			System.out.print("Find one fitness with NaN!\n");
			cur_fitness = 0;//-1.0e+30;
		}
		cur_solution.setFitness(cur_fitness);
	}

	/**
	 * Initialize the ideal point
	 * 
	 * @throws JMException
	 * @throws ClassNotFoundException
	 */
	void initIdealPoint() throws JMException, ClassNotFoundException {
		if (SASAlgorithmAdaptor.isFuzzy) {
			for (int i = 0; i < problem_.getNumberOfObjectives(); i++) {
				z_[i] = 0;
			}
			return;
		}
		
		if(fixed_bounds != null) {
			for (int i = 0; i < problem_.getNumberOfObjectives(); i++) {
				z_[i] = fixed_bounds[i][0];
			}
			return;
		}

		for (int i = 0; i < problem_.getNumberOfObjectives(); i++)
			z_[i] = 1.0e+30;

		for (int i = 0; i < populationSize_; i++)
			updateReference(population_.get(i));
	}

	/**
	 * Initialize the nadir point
	 * 
	 * @throws JMException
	 * @throws ClassNotFoundException
	 */
	void initNadirPoint() throws JMException, ClassNotFoundException {
		if (SASAlgorithmAdaptor.isFuzzy) {
			for (int i = 0; i < problem_.getNumberOfObjectives(); i++) {
				nz_[i] = 1;
			}
			return;
		}
		
		if(fixed_bounds != null) {
			for (int i = 0; i < problem_.getNumberOfObjectives(); i++) {
				nz_[i] = fixed_bounds[i][1];
			}
			return;
		}

		for (int i = 0; i < problem_.getNumberOfObjectives(); i++)
			nz_[i] = -1.0e+30;

		for (int i = 0; i < populationSize_; i++)
			updateNadirPoint(population_.get(i));
	}

	/**
	 * Update the ideal point, it is just an approximation with the best value for
	 * each objective
	 * 
	 * @param individual
	 */
	void updateReference(Solution individual) {
		if (SASAlgorithmAdaptor.isFuzzy) {
			return;
		}
		
		if(fixed_bounds != null) {
			return;
		}
		
		for (int i = 0; i < problem_.getNumberOfObjectives(); i++) {
			if (individual.getObjective(i) < z_[i])
				z_[i] = individual.getObjective(i);
		}
	}

	/**
	 * Update the nadir point, it is just an approximation with worst value for each
	 * objective
	 * 
	 * @param individual
	 */
	void updateNadirPoint(Solution individual) {
		if (SASAlgorithmAdaptor.isFuzzy) {
			return;
		}
		
		if(fixed_bounds != null) {
			return;
		}
		
		for (int i = 0; i < problem_.getNumberOfObjectives(); i++) {
			if (individual.getObjective(i) > nz_[i])
				nz_[i] = individual.getObjective(i);
		}
	}
	
	
	//CoolingScheduleHyperbolic from opt4j
	double getTemperature(int i, int n) {
		return t0 * Math.pow((tn / t0), (double) i / (double) n);
	}


	/**
	 * This is used to find the knee point from a set of solutions
	 * 
	 * @param population
	 * @return
	 */
//	public Solution kneeSelection(SolutionSet population_) {		
//		int[] max_idx    = new int[problem_.getNumberOfObjectives()];
//		double[] max_obj = new double[problem_.getNumberOfObjectives()];
//		int populationSize_ = population_.size();
//		// finding the extreme solution for f1
//		for (int i = 0; i < populationSize_; i++) {
//			for (int j = 0; j < problem_.getNumberOfObjectives(); j++) {
//				// search the extreme solution for f1
//				if (population_.get(i).getObjective(j) > max_obj[j]) {
//					max_idx[j] = i;
//					max_obj[j] = population_.get(i).getObjective(j);
//				}
//			}
//		}
//
//		if (max_idx[0] == max_idx[1])
//			System.out.println("Watch out! Two equal extreme solutions cannot happen!");
//		
//		int maxIdx;
//		double maxDist;
//		double temp1 = (population_.get(max_idx[1]).getObjective(0) - population_.get(max_idx[0]).getObjective(0)) * 
//				(population_.get(max_idx[0]).getObjective(1) - population_.get(0).getObjective(1)) - 
//				(population_.get(max_idx[0]).getObjective(0) - population_.get(0).getObjective(0)) * 
//				(population_.get(max_idx[1]).getObjective(1) - population_.get(max_idx[0]).getObjective(1));
//		double temp2 = Math.pow(population_.get(max_idx[1]).getObjective(0) - population_.get(max_idx[0]).getObjective(0), 2.0) + 
//				Math.pow(population_.get(max_idx[1]).getObjective(1) - population_.get(max_idx[0]).getObjective(1), 2.0);
//		double constant = Math.sqrt(temp2);
//		double tempDist = Math.abs(temp1) / constant;
//		maxIdx  = 0;
//		maxDist = tempDist;
//		for (int i = 1; i < populationSize_; i++) {
//			temp1 = (population_.get(max_idx[1]).getObjective(0) - population_.get(max_idx[0]).getObjective(0)) *
//					(population_.get(max_idx[0]).getObjective(1) - population_.get(i).getObjective(1)) - 
//					(population_.get(max_idx[0]).getObjective(0) - population_.get(i).getObjective(0)) * 
//					(population_.get(max_idx[1]).getObjective(1) - population_.get(max_idx[0]).getObjective(1));
//			tempDist = Math.abs(temp1) / constant;
//			if (tempDist > maxDist) {
//				maxIdx  = i;
//				maxDist = tempDist;
//			}
//		}
//		
//		return population_.get(maxIdx);
//	}
} // NSGA-II
