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

package jmetal.metaheuristics.hc;

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
 * @author keli, taochen
 *
 */
public class HC_SAS extends Algorithm {

	private SASSolutionInstantiator factory = null;

	boolean isNeigboring = true;
	private Seeder seeder = null;
	double[] weights = new double[0];
	
	private Set<String> full_set;
	// This to be used with single objective only and without any fuzzy setting
	double[][] fixed_bounds = null;
	// ideal point
	double[] z_;

	// nadir point
	double[] nz_;

	int populationSize_;

	SolutionSet population_;

	/**
	 * Constructor
	 * 
	 * @param problem Problem to solve
	 */
	public HC_SAS(Problem problem) {
		super(problem);
	} // NSGAII

	/**
	 * Constructor
	 * 
	 * @param problem Problem to solve
	 */
	public HC_SAS(Problem problem, SASSolutionInstantiator factory) {
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

		int populationSize = isNeigboring? 1 : 2;
		// Initialize the variables
		SolutionSet population = new SolutionSet(populationSize);
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

		int index = new Random().nextInt(population.get(0).numberOfVariables());
		double te = 0.0;
		Solution startingPoint = population.get(0);
		// Generations
		while (evaluations < maxEvaluations) {

			/**
			 * This is a hill climbing search, where the neighbour means the solution with only one gene change
			 * This is also a random restart hill climbing
			 */
			
			if (EAConfigure.getInstance().measurement == measurement) {
				break;
			}

			//double f = Double.MAX_VALUE;
			Solution bestS = null;
			// Only work for single objective
			if (isNeigboring) {
				
				for (int i = 0; i < startingPoint.getDecisionVariables().length; i++) {
					
					// Get the neibghour
					int[] bounds = ((SASSolution) startingPoint).getBounds(i);
					int new_v = (int) (PseudoRandom.randInt(
							// In the implementation of SASSolution, we can ensure the right boundary is 
							// always used even under variable dependency.
							bounds[0],
							bounds[1]));
					
					if (new_v == startingPoint.getDecisionVariables()[i].getValue()) {
						if(new_v + 1 > bounds[1]) {
							new_v = bounds[0];
						} else {
							new_v += 1;
						}
					}
					
					
					
					Solution nextSolution = factory.getSolution(startingPoint);
					nextSolution.getDecisionVariables()[i].setValue(new_v);

					problem_.evaluate(nextSolution);
					problem_.evaluateConstraints(nextSolution);
					
					if(SASAlgorithmAdaptor.isLogSolutionsInFull) {
						full_set.add(convertFullInfo(nextSolution));
					}
					
					updateReference(nextSolution);
					updateNadirPoint(nextSolution);

					fitnessAssignment(nextSolution);
					fitnessAssignment(population.get(0));
					if(bestS != null) {
					   fitnessAssignment(bestS);
					}
					if (bestS == null || bestS.getFitness() > nextSolution.getFitness()) {
						bestS = nextSolution;
					}

				
					
					
					if (SASAlgorithmAdaptor.logGenerationOfObjectiveValue > 0
							&& evaluations % SASAlgorithmAdaptor.logGenerationOfObjectiveValue == 0) {
						if(SASAlgorithmAdaptor.isFuzzy) {
							org.femosaa.util.Logger.logSolutionSetWithGenerationAndFuzzyValue(population, old_population,
									"SolutionSetWithGen.rtf", evaluations);
						} else {
							org.femosaa.util.Logger.logSolutionSetWithGeneration(population, "SolutionSetWithGen.rtf", 
									evaluations);
						}
					}
					
					evaluations++;
					if (!SASAlgorithmAdaptor.isInvalidSolutionConsumeMeasurement) {

						if (nextSolution.getObjective(0) == Double.MAX_VALUE
								|| nextSolution.getObjective(0) == Double.MAX_VALUE / 100) {

						} else {
							measurement += factory.record(nextSolution);
						}

					} else {
						measurement += factory.record(nextSolution);
					}
					
					if (EAConfigure.getInstance().measurement == measurement) {
						break;
					}

					if (evaluations >= maxEvaluations) {
						break;
					}
				}
				
				if (population.get(0).getFitness() > bestS.getFitness()) {
					population.clear();
					population.add(bestS);
					startingPoint = bestS;
					//System.out.print("has better one\n");
				} else {
					
					if (EAConfigure.getInstance().measurement == measurement) {
						break;
					}

					if (evaluations >= maxEvaluations) {
						break;
					}
					
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
					startingPoint = newSolution;
					fitnessAssignment(newSolution);
					
					if(SASAlgorithmAdaptor.isLogSolutionsInFull) {
						full_set.add(convertFullInfo(newSolution));
					}
				}
				

				//%%%%%%%%%%% older neigboring
				
				/*int next_index = index >= population.get(0).getDecisionVariables().length - 1 ? 0 : index + 1;
				int pre_index = index == 0 ? population.get(0).getDecisionVariables().length - 1 : index - 1;

				int[] bounds = ((SASSolution) population.get(0)).getBounds(next_index);
				for (int i = bounds[0]; i <= bounds[1]; i++) {
					Solution nextSolution = factory.getSolution(population.get(0));
					nextSolution.getDecisionVariables()[next_index].setValue(i);

					problem_.evaluate(nextSolution);
					problem_.evaluateConstraints(nextSolution);
					
					updateReference(nextSolution);
					updateNadirPoint(nextSolution);

					fitnessAssignment(nextSolution);

					if (f > nextSolution.getFitness()) {
						f = nextSolution.getFitness();
						bestS = nextSolution;
						index = next_index;
					}

					evaluations++;
					measurement += factory.record(nextSolution);
					if (EAConfigure.getInstance().measurement == measurement) {
						break;
					}

					if (evaluations >= maxEvaluations) {
						break;
					}
				}

				if (EAConfigure.getInstance().measurement > measurement && evaluations < maxEvaluations) {
					bounds = ((SASSolution) population.get(0)).getBounds(pre_index);
					for (int i = bounds[0]; i <= bounds[1]; i++) {
						Solution nextSolution = factory.getSolution(population.get(0));
						nextSolution.getDecisionVariables()[pre_index].setValue(i);

						problem_.evaluate(nextSolution);
						problem_.evaluateConstraints(nextSolution);
						
						updateReference(nextSolution);
						updateNadirPoint(nextSolution);

						fitnessAssignment(nextSolution);

						if (f > nextSolution.getFitness()) {
							f = nextSolution.getFitness();
							bestS = nextSolution;
							index = pre_index;
						}

						evaluations++;
						measurement += factory.record(nextSolution);
						if (EAConfigure.getInstance().measurement == measurement) {
							break;
						}

						if (evaluations >= maxEvaluations) {
							break;
						}
					}
				}

				if (population.get(0).getFitness() > f) {
					population.clear();
					population.add(bestS);
					System.out.print("has better one\n");
				} else {
					index = new Random().nextInt(population.get(0).numberOfVariables());
					//index = index >= population.get(0).getDecisionVariables().length - 1 ? 0 : index + 1;
				}*/
				
				//System.out.print("Measurement: " + measurement + "\n");
				
				if (SASAlgorithmAdaptor.logMeasurementOfObjectiveValue) {
					if(SASAlgorithmAdaptor.isFuzzy) {
						org.femosaa.util.Logger.logSolutionSetWithGeneration(old_population, "SolutionSetWithMeasurement.rtf", 
								measurement );
					} else {
						org.femosaa.util.Logger.logSolutionSetWithGeneration(population, "SolutionSetWithMeasurement.rtf", 
								measurement );
					}
					
				}
				
			
			
			} else {

				Solution nextSolution = factory.getSolution(population.get(0));

				((SASSolution) nextSolution).mutateWithDependency(index, true);
				index++;
				if (index >= nextSolution.getDecisionVariables().length) {
					index = 0;
				}

				problem_.evaluate(nextSolution);
				problem_.evaluateConstraints(nextSolution);

				updateReference(nextSolution);
				updateNadirPoint(nextSolution);

				if (SASAlgorithmAdaptor.isLogToD && nextSolution.getObjective(0) <= SASAlgorithmAdaptor.d
						&& te == 0.0) {
					System.out.print("Found one with " + evaluations + "\n");
					te = evaluations;
				}

				old_population.add(nextSolution);
				if (SASAlgorithmAdaptor.isFuzzy) {
					population = factory.fuzzilize(old_population);
				}

				fitnessAssignment(population.get(0));
				fitnessAssignment(population.get(1));

				if (population.get(0).getFitness() < population.get(1).getFitness()) {
					population.remove(1);
					old_population.remove(1);

				} else {
					population.remove(0);
					old_population.remove(0);
				}

				evaluations++;
				if (SASAlgorithmAdaptor.logGenerationOfObjectiveValue > 0
						&& evaluations % SASAlgorithmAdaptor.logGenerationOfObjectiveValue == 0) {
					if(SASAlgorithmAdaptor.isFuzzy) {
						org.femosaa.util.Logger.logSolutionSetWithGenerationAndFuzzyValue(population, old_population,
								"SolutionSetWithGen.rtf", evaluations);
					} else {
						org.femosaa.util.Logger.logSolutionSetWithGeneration(population, "SolutionSetWithGen.rtf", 
								evaluations);
					}
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
