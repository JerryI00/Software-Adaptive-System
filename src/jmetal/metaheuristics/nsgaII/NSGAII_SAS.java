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

package jmetal.metaheuristics.nsgaII;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import org.apache.commons.math3.distribution.BetaDistribution;
import org.femosaa.core.EAConfigure;
import org.femosaa.core.SASAlgorithmAdaptor;
import org.femosaa.core.SASSolution;
import org.femosaa.core.SASSolutionInstantiator;
import org.femosaa.invalid.SASValidityAndInvalidityCoEvolver;
import org.femosaa.seed.Seeder;

import jmetal.core.*;
import jmetal.util.comparators.CrowdingComparator;
import jmetal.util.*;

/**
 * 
 * @author keli, taochen
 *
 */
public class NSGAII_SAS extends Algorithm {

	private SASSolutionInstantiator factory = null;
	
	private SASValidityAndInvalidityCoEvolver vandInvCoEvolver = null;
	private Seeder seeder = null;
	SolutionSet population_;
	boolean logOnce = false;
	
	private int same_optimum_count = 0;
	private double mmo_weight = 1.0;//Double.MIN_VALUE;
	private double mmo_weight_step = 0.1;//0.00001;
	private double global_optimum = Double.MAX_VALUE;
	private HashMap<Double,Double> q_table = new HashMap<Double,Double> ();
	private double current_proportion = -1;
	private int target = 0;
	
	private Set<String> full_set;
	
	
	String weight_str = "";
	String proportion_str = "";
	/**
	 * Constructor
	 * @param problem Problem to solve
	 */
	public NSGAII_SAS(Problem problem) {
		super (problem) ;
	} // NSGAII


  	/**
  	 * Constructor
  	 * @param problem Problem to solve
  	 */
	public NSGAII_SAS(Problem problem, SASSolutionInstantiator factory) {
		super(problem);
        this.factory = factory;
	}

	/**   
	 * Runs the NSGA-II algorithm.
	 * @return a <code>SolutionSet</code> that is a set of non dominated solutions
	 * as a result of the algorithm execution
	 * @throws JMException 
	 */
	public SolutionSet execute() throws JMException, ClassNotFoundException {
		
		if(SASAlgorithmAdaptor.isLogSolutionsInFull) {
			full_set = new HashSet<String>();
		}
		
		if (factory == null) {
			throw new RuntimeException("No instance of SASSolutionInstantiator found!");
		}
		
		int type;
		
		int populationSize;
		int maxEvaluations;
		int evaluations;
		int measurement;

		int requiredEvaluations; // Use in the example of use of the
		// indicators object (see below)

		// knee point which might be used as the output
		Solution kneeIndividual = factory.getSolution(problem_);
		if(getInputParameter("vandInvCoEvolver") != null) {
		    vandInvCoEvolver = (SASValidityAndInvalidityCoEvolver)getInputParameter("vandInvCoEvolver");
		}
		if(getInputParameter("seeder") != null) {
			seeder = (Seeder)getInputParameter("seeder");
		}
		SolutionSet population;
		SolutionSet offspringPopulation;
		SolutionSet union;

		Operator mutationOperator;
		Operator crossoverOperator;
		Operator selectionOperator;

		Distance distance = new Distance();

		//Read the parameters
		populationSize = ((Integer) getInputParameter("populationSize")).intValue();
		maxEvaluations = ((Integer) getInputParameter("maxEvaluations")).intValue();

		//Initialize the variables
		population = new SolutionSet(populationSize);
		evaluations = 0;
		measurement = 0;

		requiredEvaluations = 0;
		//System.out.print("maxEvaluations " + maxEvaluations + "\n");
		//Read the operators
		mutationOperator  = operators_.get("mutation");
		crossoverOperator = operators_.get("crossover");
		selectionOperator = operators_.get("selection");

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
			} // for
		}

		if(SASAlgorithmAdaptor.isLogSolutionsInFull) {
			Iterator itr = population.iterator();
			while(itr.hasNext()) {
				full_set.add(convertFullInfo((Solution)itr.next()));
			}
			
		}
		
		SolutionSet old_population = new SolutionSet(populationSize);
		if(SASAlgorithmAdaptor.isFuzzy) {
			old_population = population;
			population = factory.fuzzilize(population);
		}
		
		if (SASAlgorithmAdaptor.logGenerationOfObjectiveValue > 0 || SASAlgorithmAdaptor.logMeasurementOfObjectiveValue) {			
			org.femosaa.util.Logger.logSolutionSetWithGeneration(population,
					"InitialSolutionSet.rtf", 0);
		}
		
		if (SASAlgorithmAdaptor.logMeasurementOfObjectiveValue) {
			if(SASAlgorithmAdaptor.isFuzzy) {
				org.femosaa.util.Logger.logSolutionSetWithGeneration(old_population, "SolutionSetWithMeasurement.rtf", 
						measurement );
			} else {
				org.femosaa.util.Logger.logSolutionSetWithGeneration(population, "SolutionSetWithMeasurement.rtf", 
						measurement );
			}
			
		}
		
		if (SASAlgorithmAdaptor.logMeasurementOfFuzzyObjectiveValue) {
			if(SASAlgorithmAdaptor.isFuzzy) {
				org.femosaa.util.Logger.logSolutionSetWithGeneration(population, "SolutionSetWithFuzzyMeasurement.rtf", 
						measurement );
			}
			
		}

		if(vandInvCoEvolver != null) {
			for (int i = 0; i < populationSize; i++) {
				newSolution = factory.getSolution(problem_);
				vandInvCoEvolver.createInitialSolution(newSolution, problem_);
//				if(vandInvCoEvolver.createInitialSolution(newSolution, problem_)){
//					evaluations++;
//					population.add(newSolution);
//				}
				
			} //for  
		}
		
		long time = Long.MAX_VALUE;
		
		// Generations 
		while (evaluations < maxEvaluations || (evaluations >= maxEvaluations && (System.currentTimeMillis() - time) < SASAlgorithmAdaptor.seed_time )) {
			//System.out.print("no" + evaluations + "***eval\n");
			
			if(EAConfigure.getInstance().measurement == measurement) {
				break;
			}
			
			
			
//
//			Iterator itr = population.iterator();
//			double no = 0.0;
//			while(itr.hasNext()) {
//				Solution s = (Solution)itr.next();
//				if(((SASSolution)s).isFromInValid) {
//					no++;
//				}
//			}
//			System.out.print("("+evaluations+","+(no/population.size()) + ")\n");
			// Create the offSpring solutionSet      
			offspringPopulation = new SolutionSet(populationSize);
			Solution[] parents = new Solution[2];		
			//for (int i = 0; i < (populationSize / 2); i++) {
			while (offspringPopulation.size() < populationSize) {
				int c = 2;
				if (evaluations < maxEvaluations || (evaluations >= maxEvaluations)) {
					Solution[] offSpring = null;
					//obtain parents			
					if(vandInvCoEvolver != null) {
						//parents[0] = (Solution) vandInvCoEvolver.doMatingSelection(population);
						//parents[1] = (Solution) vandInvCoEvolver.doMatingSelection(population);
						parents[0] = (Solution) vandInvCoEvolver.doMatingSelection(population, true);
						parents[1] = (Solution) vandInvCoEvolver.doMatingSelection(population, false);
                        offSpring = vandInvCoEvolver.doReproduction(parents, problem_);
						
						for(Solution s : offSpring) {
							if(offspringPopulation.size() >= populationSize) {
								break;
							}
							offspringPopulation.add(s);
							evaluations++;
							c--;
							if(((SASSolution)parents[0]).isFromInValid || ((SASSolution)parents[1]).isFromInValid) {
								((SASSolution)s).isFromInValid = true;
							}
						}
					} 
					//System.out.print("front: " + new Ranking(population).getNumberOfSubfronts() + "\n");
					if(SASAlgorithmAdaptor.isToFilterRedundantSolutionMating) {
						
						/*SolutionSet new_set = new SolutionSet();

						for (int i = 0; i < population.size(); i++) {
							boolean same = false;

							for (int j = 0; j < new_set.size(); j++) {
								same = isSame(population.get(i), new_set.get(j), true);
								if (same) {
									break;
								}
							}

							if (!same) {
								new_set.add(population.get(i));
							}

						}
						//System.out.print("size " + new_set.size() + "\n");
						
						parents[0] = (Solution) selectionOperator.execute(new_set);
						parents[1] = (Solution) selectionOperator.execute(new_set);*/
						
						parents[0] =  RemoveRedundantAdaptor.matingRemoval(population);
						parents[1] =  RemoveRedundantAdaptor.matingRemoval(population);
					} else {
						parents[0] = (Solution) selectionOperator.execute(population);
						parents[1] = (Solution) selectionOperator.execute(population);
					}
				
					
				
					//offSpring  = new Solution[2];
					//offSpring[0] = factory.getSolution(parents[0]);
					//offSpring[1] = factory.getSolution(parents[1]);
				
					offSpring = (Solution[]) crossoverOperator.execute(parents);
					
					
					mutationOperator.execute(offSpring[0]);
					mutationOperator.execute(offSpring[1]);
					
					if (SASAlgorithmAdaptor.logPreivousAndCurrentPopToBest) {
						if(SASAlgorithmAdaptor.isFuzzy) {
							((SASSolution)offSpring[0]).setParents(old_population.get(((SASSolution)parents[0]).index),old_population.get(((SASSolution)parents[1]).index));
							((SASSolution)offSpring[1]).setParents(old_population.get(((SASSolution)parents[0]).index),old_population.get(((SASSolution)parents[1]).index));
						} else {
							((SASSolution)offSpring[0]).setParents(parents[0], parents[1]);
							((SASSolution)offSpring[1]).setParents(parents[0], parents[1]);
						}
						
					}
					
					//long test_time = System.currentTimeMillis();
					problem_.evaluate(offSpring[0]);
					problem_.evaluateConstraints(offSpring[0]);
					if (!SASAlgorithmAdaptor.isInvalidSolutionConsumeMeasurement) {

						if (offSpring[0].getObjective(0) == Double.MAX_VALUE
								|| offSpring[0].getObjective(0) == Double.MAX_VALUE / 100) {

						} else {
							measurement += factory.record(offSpring[0]);
						}

					} else {
						measurement += factory.record(offSpring[0]);
					}
					if(EAConfigure.getInstance().measurement == measurement) {
						break;
					}
					problem_.evaluate(offSpring[1]);
					problem_.evaluateConstraints(offSpring[1]);
					if (!SASAlgorithmAdaptor.isInvalidSolutionConsumeMeasurement) {

						if (offSpring[1].getObjective(0) == Double.MAX_VALUE
								|| offSpring[1].getObjective(0) == Double.MAX_VALUE / 100) {

						} else {
							measurement += factory.record(offSpring[1]);
						}

					} else {
						measurement += factory.record(offSpring[1]);
					}
					if(EAConfigure.getInstance().measurement == measurement) {
						break;
					}
					//System.out.print("Evaluation time: " + (System.currentTimeMillis()-test_time) + "\n");
					if(c == 0) {
						continue;
					}
					if(offspringPopulation.size() >= populationSize) {
						break;
					}
					offspringPopulation.add(offSpring[0]);
					evaluations++;
					c--;
					if(c == 0) {
						continue;
					}
					if(offspringPopulation.size() >= populationSize) {
						break;
					}
					offspringPopulation.add(offSpring[1]);
					evaluations++;
					if(((SASSolution)parents[0]).isFromInValid || ((SASSolution)parents[1]).isFromInValid) {
						((SASSolution)offSpring[0]).isFromInValid = true;
						((SASSolution)offSpring[1]).isFromInValid = true;
					}
					
					if(SASAlgorithmAdaptor.isLogSolutionsInFull) {
						full_set.add(convertFullInfo((Solution)offSpring[0]));
						full_set.add(convertFullInfo((Solution)offSpring[1]));
					}
					
				} // if                            
			} // for
			
			SolutionSet old_union = null;
			
//			if (SASAlgorithmAdaptor.isToFilterRedundantSolution) {
//
//				int target = 0;
//				Set<String> record = new HashSet<String>();
//				List<Solution> removed = new ArrayList<Solution>();
//				if (SASAlgorithmAdaptor.isFuzzy) {
//					old_population = filter(old_population, record, removed, target);
//				} else {
//					population = filter(population, record, removed,  target);
//				}
//				
//				offspringPopulation = filter(offspringPopulation, record, removed,  target);
//				
//				int size = 0;
//				if (SASAlgorithmAdaptor.isFuzzy) {
//					size = old_population.size();
//				} else {
//					size = population.size();
//				}
//				
//				if (SASAlgorithmAdaptor.isFuzzy) {
//					printSolutions(old_population);
//				} else {
//					printSolutions(population);
//				}
//				printSolutions(offspringPopulation);
//			
//				
//				size = size + offspringPopulation.size();
//				
//				System.out.print("Unique: " + size + "\n");
//				
//				if (size < populationSize) {
//					int l = populationSize - size;
//					for (int i = 0; i < l; i++) {
//						offspringPopulation.add(removed.get(i));
//					}
//				} else {
//					
//				}
//											
//			}
//			
		
			
			// Create the solutionSet union of solutionSet and offSpring			
			if(SASAlgorithmAdaptor.isFuzzy) {			
				union = ((SolutionSet) old_population).union(offspringPopulation);
				old_union = union;
				if(SASAlgorithmAdaptor.isBoundNormalizationForTarget) {
					((SASSolution)old_population.get(0)).resetNormalizationBounds(0);
					((SASSolution)old_population.get(0)).resetNormalizationBounds(1);
					/*for(int i = 0; i < union.size(); i++) {
						((SASSolution)union.get(i)).updateNormalizationBounds(new double[] {union.get(i).getObjective(0),
								union.get(i).getObjective(1)});
					}*/
					
					for(int i = 0; i < old_population.size(); i++) {
						((SASSolution)old_population.get(i)).updateNormalizationBounds(new double[] {old_population.get(i).getObjective(0),
								old_population.get(i).getObjective(1)});
					}
				}
				
				if(SASAlgorithmAdaptor.isAdaptiveWeightinMMO) {

					boolean isChange = true;
					double pre = global_optimum;
					for (int i = 0; i < union.size(); i++) {
						// Assume 0 is the target objective
						if (union.get(i).getObjective(target) < global_optimum) {
							global_optimum = union.get(i).getObjective(target);
							isChange = false;
						}
					}
					
					if(q_table.size() == 0) {
						//q_table.put(0.1, 0.1);
						//q_table.put(0.2, 0.1);
						q_table.put(0.05, 0.1);
						q_table.put(0.1, 0.1);
						q_table.put(0.3, 0.1);
						q_table.put(0.5, 0.1);
						q_table.put(0.7, 0.1);
						q_table.put(0.9, 0.1);
						q_table.put(1.0, 0.1);
						/*q_table.put(0.6, 1.0 / 8);
						q_table.put(0.7, 1.0 / 8);
						q_table.put(0.8, 1.0 / 8);
						q_table.put(0.9, 1.0 / 8);*/
					}

					
					if (isChange) {
						if(current_proportion > 0 && q_table.size() != 0) {
							double p = q_table.get(current_proportion);
							//if (p > 1.0 / 8) {
								//q_table.put(current_proportion, p - 1.0 / 8);
							//q_table.put(current_proportion, p / 2);
							//if(p-0.1 < 0.1) {
							double d = (double) measurement / EAConfigure.getInstance().measurement;
                        	q_table.put(current_proportion, p / (1 + d));
							//q_table.put(current_proportion, p / 1.1);
                        	/*} else {
								q_table.put(current_proportion, p-0.1);
							}*/
							
							//}
						}
						same_optimum_count++;
					} else {
                        if(current_proportion > 0 && q_table.size() != 0) {
                        	double p = q_table.get(current_proportion);
							//q_table.put(current_proportion, p + 1.0 / 8);
                        	//q_table.put(current_proportion, p * 2);
                        	double d = (double) measurement / EAConfigure.getInstance().measurement;
                        	q_table.put(current_proportion, p * (1 + d));
                        	//q_table.put(current_proportion, p * 1.1);
						}
						same_optimum_count = 0;
					}

					if (MMOWeightAdaptor.isAdapt(same_optimum_count, measurement, EAConfigure.getInstance().measurement)) {
					//if(1==1) {
					
					    System.out.print("trigger weight change...\n");
					    current_proportion = 0.3;//0.3;
						//current_proportion = MMOWeightAdaptor.nondominatedProportion(q_table);
						System.out.print("current proportion: " + current_proportion +"\n");
						mmo_weight = MMOWeightAdaptor.smartGetWeight(problem_, factory, union, mmo_weight, -1, populationSize, mmo_weight_step, current_proportion);					
						System.out.print("weight = " + mmo_weight + "\n");

					}
					
					weight_str += "("+measurement + ","+mmo_weight+")\n";
					proportion_str += "("+measurement + ","+MMOWeightAdaptor.testProportion(factory, union, mmo_weight)+")\n";
					
					//mmo_weight = 1.0;
					union = factory.fuzzilize(union, mmo_weight);
				} else {
					weight_str += "("+measurement + ","+mmo_weight+")\n";
					proportion_str += "("+measurement + ","+MMOWeightAdaptor.testProportion(factory, union, mmo_weight)+")\n";
					
					union = factory.fuzzilize(union);//((double) measurement - 1) / (EAConfigure.getInstance().measurement - 1)
				}
				
			} else {
				union = ((SolutionSet) population).union(offspringPopulation);
				
				
			}
			// Create the solutionSet union of solutionSet and offSpring
			//union = ((SolutionSet) population).union(offspringPopulation);
			
			
	       
			
			if (SASAlgorithmAdaptor.logMeasurementOfObjectiveValueTwoPop) {
				if(SASAlgorithmAdaptor.isFuzzy) {
					org.femosaa.util.Logger.logSolutionSetWithGeneration(old_union, "SolutionSetWithMeasurementTwoPop.rtf", 
							measurement );
				} else {
					org.femosaa.util.Logger.logSolutionSetWithGeneration(union, "SolutionSetWithMeasurementTwoPop.rtf", 
							measurement );
				}
			}
			
			
			boolean isLog = false;
			
			if (SASAlgorithmAdaptor.logPreivousAndCurrentPopToBest) {
				
				if(!logOnce) {
					
					logOnce = org.femosaa.util.Logger.logSolutionSetWithGenerationOnBest(offspringPopulation, "SolutionSetWithMeasurementParentsOfBest.rtf", 
							measurement );
					
					if(logOnce) {
						if(SASAlgorithmAdaptor.isFuzzy) {
							org.femosaa.util.Logger.logSolutionSetWithGenerationAndBreaket(old_population, "SolutionSetWithMeasurementPreviousPop.rtf", 
									measurement );
						} else {
							org.femosaa.util.Logger.logSolutionSetWithGenerationAndBreaket(population, "SolutionSetWithMeasurementPreviousPop.rtf", 
									measurement );
						}
						
						isLog = true;
					}
					
					
				}
				
				
			}
			
			// for remove redundant ***********
			/*List<Solution> removedRedundant = new ArrayList<Solution>();
			if (SASAlgorithmAdaptor.isToFilterRedundantSolution) {
				//int target = 0;
				Set<String> record = new HashSet<String>();
				
		
				//System.out.print("old " + front.size() + "\n");
				union = filter(union, record, removedRedundant,  target);
				SolutionSet ss = new SolutionSet();
				for (int i = 0; i < removedRedundant.size(); i++) {
					ss.add(removedRedundant.get(i));
				}
				Ranking r = new Ranking(ss);
				removedRedundant.clear();
				for (int i = 0; i < r.getNumberOfSubfronts(); i++) {
					
					SolutionSet sub = r.getSubfront(i);
					distance.crowdingDistanceAssignment(sub, problem_.getNumberOfObjectives());
					sub.sort(new CrowdingComparator());
					for (int j = 0; j < sub.size(); j++) {
						removedRedundant.add(sub.get(j));
					}
				}
				
				System.out.print("unique " + union.size() + "\n");
				//System.out.print("removedRedundant " + removedRedundant.size() + "\n");
				//printSolutions(front);
			}*/
			// for remove redundant ***********
			// Ranking the union
			Ranking ranking = new Ranking(union);
			
			
			if (SASAlgorithmAdaptor.logNumberOfNoDominated) {
				/*if(SASAlgorithmAdaptor.isFuzzy) {
					org.femosaa.util.Logger.logNumberOfNoDominated(old_population, "NondominatedCount.rtf", 
							evaluations );
				} else {
					org.femosaa.util.Logger.logNumberOfNoDominated(population, "NondominatedCount.rtf", 
							evaluations );
				}*/
				org.femosaa.util.Logger.logNumberOfNoDominated(ranking.getSubfront(0), "NondominatedCount.rtf", 
						evaluations );
				if(EAConfigure.getInstance().measurement == measurement) {
					org.femosaa.util.Logger.logNumberOfNoDominated(ranking.getSubfront(0), "FinalNondominatedCount.rtf", 
							evaluations );
				}
			}

			int remain = populationSize;
			int index = 0;
			SolutionSet front = null;
			population.clear();
			old_population.clear();
			
			// Obtain the next front
			front = ranking.getSubfront(index);
			
		
		
			

		// for remove redundant ***********
			if(SASAlgorithmAdaptor.isToFilterRedundantSolution) {
				
				Set<String> record = new HashSet<String>();
				SolutionSet[] sets = new SolutionSet[ranking.getNumberOfSubfronts()];
				SolutionSet pre = new SolutionSet();
				
				for (int i = 0; i < ranking.getNumberOfSubfronts(); i++) {
					record.clear();
					
					sets[i] = new SolutionSet();
					SolutionSet front_set = new SolutionSet();
					
					for (int j = 0; j < ranking.getSubfront(i).size(); j++) {
						front_set.add(ranking.getSubfront(i).get(j));
					}
					
					for (int j = 0; j < pre.size(); j++) {
						front_set.add(pre.get(j));
					}
					
					pre.clear();
					
					
					for (int k = 0; k < front_set.size(); k++) {
						String v = "";
						for (int j = 0; j < front_set.get(k).getDecisionVariables().length; j++) {
							v += front_set.get(k).getDecisionVariables()[j].getValue()+",";
						}
						if(record.contains(v)) {
							pre.add(front_set.get(k));
						} else {
							sets[i].add(front_set.get(k));
						}
					}
					
					
				}
				
				front = sets[index];
				
				while ((remain > 0) && (remain >= front.size())) {
					//Assign crowding distance to individuals
					distance.crowdingDistanceAssignment(front, problem_.getNumberOfObjectives());
					//Add the individuals of this front
					for (int k = 0; k < front.size(); k++) {
						population.add(front.get(k));
						if(SASAlgorithmAdaptor.isFuzzy) {
							old_population.add(factory.defuzzilize(front.get(k), old_union));
							if (SASAlgorithmAdaptor.logPreivousAndCurrentPopToBest) {
								((SASSolution)front.get(k)).index = old_population.size()-1;
							}
						}
					} // for

					//Decrement remain
					remain = remain - front.size();

					//Obtain the next front
					index++;
					// for remove redundant ***********
					
					// for remove redundant ***********
					if (remain > 0) {
					     front = sets[index];
						
					} // if  
					
					
				} // while
				
				// Remain is less than front(index).size, insert only the best one
				if (remain > 0) {  // front contains individuals to insert                        
					distance.crowdingDistanceAssignment(front, problem_.getNumberOfObjectives());
					front.sort(new CrowdingComparator());
					for (int k = 0; k < remain; k++) {
						population.add(front.get(k));
						if(SASAlgorithmAdaptor.isFuzzy) {
							old_population.add(factory.defuzzilize(front.get(k), old_union));
							if (SASAlgorithmAdaptor.logPreivousAndCurrentPopToBest) {
								((SASSolution)front.get(k)).index = old_population.size()-1;
							}
						}
					} // for

					remain = 0;
				} // if 
				
				
				/*Set<String> record = new HashSet<String>();
				List<SolutionSet> removedRedundant = new ArrayList<SolutionSet>();
				
				for (int k = 0; k < ranking.getNumberOfSubfronts(); k++) {
					removedRedundant.add(new SolutionSet());
				}
				
				while ((remain > 0)) {
					//Assign crowding distance to individuals
					distance.crowdingDistanceAssignment(front, problem_.getNumberOfObjectives());
					//Add the individuals of this front
					
					SolutionSet new_set = new SolutionSet();
					for (int k = 0; k < front.size(); k++) {
						
						String v = "";
						for (int j = 0; j < front.get(k).getDecisionVariables().length; j++) {
							v += front.get(k).getDecisionVariables()[j].getValue()+",";
						}
						
						if(record.contains(v)) {
							removedRedundant.get(index).add(front.get(k));
						} else {
							new_set.add(front.get(k));
							
							record.add(v);
						}
						
					
					} // for
					
                    if(remain < new_set.size()) {
						break;
					}
					
					for (int k = 0; k < new_set.size(); k++) {
						population.add(new_set.get(k));
						if(SASAlgorithmAdaptor.isFuzzy) {
							old_population.add(factory.defuzzilize(new_set.get(k), old_union));
							if (SASAlgorithmAdaptor.logPreivousAndCurrentPopToBest) {
								((SASSolution)new_set.get(k)).index = old_population.size()-1;
							}
						}
					}
					
					

					//Decrement remain
					remain = remain - new_set.size();
					
					

					//Obtain the next front
					index++;
					// for remove redundant ***********
					if(index >= ranking.getNumberOfSubfronts()) {
						break;
					}
					// for remove redundant ***********
					if (remain > 0) {
					     front = ranking.getSubfront(index);
						
					} // if  
					
					
				} // while
				
				if(index >= ranking.getNumberOfSubfronts() && remain > 0) {
					
					for (int i = 0; i < removedRedundant.size(); i++) {
						
						for (int k = 0; k < removedRedundant.get(i).size(); k++) {
						
							removedRedundant.get(i).sort(new CrowdingComparator());
							population.add(removedRedundant.get(i).get(k));
							if(SASAlgorithmAdaptor.isFuzzy) {
								old_population.add(factory.defuzzilize(removedRedundant.get(i).get(k), old_union));
								if (SASAlgorithmAdaptor.logPreivousAndCurrentPopToBest) {
									((SASSolution)removedRedundant.get(i).get(k)).index = old_population.size()-1;
								}
							}
							remain--;
							if(remain == 0) {
								break;
							}
						}
						
						if(remain == 0) {
							break;
						}
						
					}
					
				} else {

					// Remain is less than front(index).size, insert only the best one
					if (remain > 0) {  // front contains individuals to insert   
						record.clear();
						distance.crowdingDistanceAssignment(front, problem_.getNumberOfObjectives());
						front.sort(new CrowdingComparator());
						for (int k = 0; k < front.size(); k++) {
							
							if(remain <= 0) {
								break;
							}
							
							String v = "";
							for (int j = 0; j < front.get(k).getDecisionVariables().length; j++) {
								v += front.get(k).getDecisionVariables()[j].getValue()+",";
							}
							
							if(record.contains(v)) {
								//removedRedundant.get(index).add(front.get(k));
							} else {
								remain--;
								population.add(front.get(k));
								if(SASAlgorithmAdaptor.isFuzzy) {
									old_population.add(factory.defuzzilize(front.get(k), old_union));
									if (SASAlgorithmAdaptor.logPreivousAndCurrentPopToBest) {
										((SASSolution)front.get(k)).index = old_population.size()-1;
									}
								}
								
								record.add(v);
							}
							
							
							
						} // for
						
						if(remain > 0) {
							for (int k = 0; k < remain; k++) {
								removedRedundant.get(index).sort(new CrowdingComparator());
								population.add(removedRedundant.get(index).get(k));
								if(SASAlgorithmAdaptor.isFuzzy) {
									old_population.add(factory.defuzzilize(removedRedundant.get(index).get(k), old_union));
									if (SASAlgorithmAdaptor.logPreivousAndCurrentPopToBest) {
										((SASSolution)removedRedundant.get(index).get(k)).index = old_population.size()-1;
									}
								}
							}
						}

						remain = 0;
					} // if 
				}*/
				
				
			} else {
				

				while ((remain > 0) && (remain >= front.size())) {
					//Assign crowding distance to individuals
					distance.crowdingDistanceAssignment(front, problem_.getNumberOfObjectives());
					//Add the individuals of this front
					for (int k = 0; k < front.size(); k++) {
						population.add(front.get(k));
						if(SASAlgorithmAdaptor.isFuzzy) {
							old_population.add(factory.defuzzilize(front.get(k), old_union));
							if (SASAlgorithmAdaptor.logPreivousAndCurrentPopToBest) {
								((SASSolution)front.get(k)).index = old_population.size()-1;
							}
						}
					} // for

					//Decrement remain
					remain = remain - front.size();

					//Obtain the next front
					index++;
					// for remove redundant ***********
					/*if(index >= ranking.getNumberOfSubfronts()) {
						break;
					}*/
					// for remove redundant ***********
					if (remain > 0) {
					     front = ranking.getSubfront(index);
						
					} // if  
					
					
				} // while
				
				// Remain is less than front(index).size, insert only the best one
				if (remain > 0) {  // front contains individuals to insert                        
					distance.crowdingDistanceAssignment(front, problem_.getNumberOfObjectives());
					front.sort(new CrowdingComparator());
					for (int k = 0; k < remain; k++) {
						population.add(front.get(k));
						if(SASAlgorithmAdaptor.isFuzzy) {
							old_population.add(factory.defuzzilize(front.get(k), old_union));
							if (SASAlgorithmAdaptor.logPreivousAndCurrentPopToBest) {
								((SASSolution)front.get(k)).index = old_population.size()-1;
							}
						}
					} // for

					remain = 0;
				} // if 
			}
			
		  
			
			
			if(vandInvCoEvolver != null) {
				vandInvCoEvolver.doEnvironmentalSelection(population);
			}
			if(SASAlgorithmAdaptor.isLogTheEvalNeededToRemiveNonSeed) {
				org.femosaa.util.Logger.printMarkedSolution(population, evaluations);
			}
		
			
			/*if(SASAlgorithmAdaptor.isFuzzy) {
				for (int i = 0; i < population.size(); i++) {
					System.out.print("***\n");
					System.out.print("fuzzy value = " + population.get(i).getObjective(0) + ":" + population.get(i).getObjective(1) + "\n");
					System.out.print("orignal value = " + factory.defuzzilize(population.get(i), old_union).getObjective(0) + ":" + factory.defuzzilize(population.get(i), old_union).getObjective(1) + "\n");
					System.out.print("***\n");
				}
			}*/
			
			if(SASAlgorithmAdaptor.logGenerationOfObjectiveValue > 0 && evaluations%SASAlgorithmAdaptor.logGenerationOfObjectiveValue == 0) {
				if(SASAlgorithmAdaptor.isFuzzy) {
					org.femosaa.util.Logger.logSolutionSetWithGeneration(old_population, "SolutionSetWithGen.rtf", 
							evaluations );
				} else {
					org.femosaa.util.Logger.logSolutionSetWithGeneration(population, "SolutionSetWithGen.rtf", 
							evaluations );
				}
				
				//org.femosaa.util.Logger.logSolutionSetValuesWithGen(population, "SolutionSetValuesWithGen.rtf", 
						//evaluations );
			}
			
			if (SASAlgorithmAdaptor.logMeasurementOfObjectiveValue) {
				if(SASAlgorithmAdaptor.isFuzzy) {
					org.femosaa.util.Logger.logSolutionSetWithGeneration(old_population, "SolutionSetWithMeasurement.rtf", 
							measurement );					
				} else {
					org.femosaa.util.Logger.logSolutionSetWithGeneration(population, "SolutionSetWithMeasurement.rtf", 
							measurement );
				}
				
			}
			
			if (SASAlgorithmAdaptor.logMeasurementOfFuzzyObjectiveValue) {
				if(SASAlgorithmAdaptor.isFuzzy) {
					org.femosaa.util.Logger.logSolutionSetWithGeneration(population, "SolutionSetWithFuzzyMeasurement.rtf", 
							measurement );
				}
				
			}
			
			if(SASAlgorithmAdaptor.logMeasurementOfRemovedObjectiveValue) {
				SolutionSet set = new SolutionSet();
				for (int k = 0; k < union.size(); k++) {
					
					boolean removed = true;
					for (int i = 0; i < population.size(); i++) {
						if (union.get(k) == population. get(i)) {
							removed = false;
							break;
						}
					}
					
					if(removed) {
						set.add(union.get(k));
					}
					
				}
				
				org.femosaa.util.Logger.logSolutionSetWithGeneration(set, "SolutionSetWithFuzzyMeasurementForRemoved.rtf", 
						measurement );					
		
				SolutionSet o_set = new SolutionSet();
				for (int k = 0; k < set.size(); k++) {
					o_set.add(factory.defuzzilize(set.get(k), old_union));
				}
				
				org.femosaa.util.Logger.logSolutionSetWithGeneration(o_set, "SolutionSetWithMeasurementForRemoved.rtf", 
						measurement );
				
			}
			
            if (SASAlgorithmAdaptor.logPreivousAndCurrentPopToBest) {
				
				if(isLog) {
					
						if(SASAlgorithmAdaptor.isFuzzy) {
							org.femosaa.util.Logger.logSolutionSetWithGenerationAndBreaket(old_population, "SolutionSetWithMeasurementCurrentPop.rtf", 
									measurement );
						} else {
							org.femosaa.util.Logger.logSolutionSetWithGenerationAndBreaket(population, "SolutionSetWithMeasurementCurrentPop.rtf", 
									measurement );
						}
						
					
					
					
				}
				
				
			}
			
			if(SASAlgorithmAdaptor.isLogDiscardedSolutions) {
				SolutionSet p = new SolutionSet(union.size() - population.size());
			
				for (int k = 0; k < union.size(); k++) {
					boolean has = false;
					for (int i = 0; i < population.size(); i++) {
						if(union.get(k).equals(population.get(i))){
							has = true;
							break;
						}
					}
					
					if(!has) {
						p.add(union.get(k));
					}
				}
				
				org.femosaa.util.Logger.logSolutionSetWithGeneration(p, "DiscardSolutionSetWithGen.rtf", 
						evaluations );
				org.femosaa.util.Logger.logSolutionSetValuesWithGen(p, "DiscardSolutionSetValuesWithGen.rtf", 
						evaluations );
				
			}
			
			if(evaluations >= maxEvaluations && time == Long.MAX_VALUE) {
				time = System.currentTimeMillis();
			}
			//System.out.print("ranking.getNumberOfSubfronts(): " + new Ranking(population).getNumberOfSubfronts() + "\n");
			
		} // while
//		Iterator itr = population.iterator();
//		double no = 0.0;
//		while(itr.hasNext()) {
//			Solution s = (Solution)itr.next();
//			if(((SASSolution)s).isFromInValid) {
//				no++;
//			}
//		}
//		System.out.print("("+evaluations+","+(no/population.size()) + ")\n");
		
		
		if(SASAlgorithmAdaptor.isFuzzy) {
			population = old_population;
			org.femosaa.util.Logger.logFinalEvaluation("FinalEvaluationCount.rtf", evaluations);
		}
		
		if(SASAlgorithmAdaptor.isLogSolutionsInFull) {
			org.femosaa.util.Logger.logSolutionFull(full_set, "FullSolution.rtf");
		}
		
		/*if (SASAlgorithmAdaptor.logMeasurementOfObjectiveValue) {
			org.femosaa.util.Logger.logSolutionSetWithGeneration(population, "SolutionSetWithMeasurement.rtf", 
					measurement );
		}*/
		// Return as output parameter the required evaluations
		setOutputParameter("evaluations", requiredEvaluations);
		population_ = population;
		
		System.out.print("-------final evalution: " + evaluations + "-------\n");
		
		
		System.out.print(weight_str + "\n\n\n");
		System.out.print(proportion_str + "\n\n\n");
		
		// Return the first non-dominated front
//		Ranking ranking = new Ranking(population);
//		return ranking.getSubfront(0);
		return population;
	} // execute
	
	
	public SolutionSet doRanking(SolutionSet population){
		Ranking ranking = new Ranking(population);
		return ranking.getSubfront(0);
	}
	
	private void printSolutions(SolutionSet pop) throws JMException {
		System.out.print("-----------------\n");
		for (int i = 0; i < pop.size(); i++) {
			
		   
			String v = "";
			for (int j = 0; j < pop.get(i).getDecisionVariables().length; j++) {
				v += pop.get(i).getDecisionVariables()[j].getValue()+",";
			}
			
			v += "=" + pop.get(i).getObjective(0) + ":" + pop.get(i).getObjective(1);
			
			System.out.print(v + "\n");
			
		}
		System.out.print("-----------------\n");
	}
	
	
	
	
	private SolutionSet filter(SolutionSet pop, Set<String> record, List<Solution> removed, final int target) throws JMException {
		SolutionSet newSet = new SolutionSet();
		
		
		List<Integer> index = new ArrayList<Integer>();
		
		
		for (int i = 0; i < pop.size(); i++) {
			
			String v = "";
			for (int j = 0; j < pop.get(i).getDecisionVariables().length; j++) {
				v += pop.get(i).getDecisionVariables()[j].getValue()+",";
			}
			
			if(record.contains(v) || pop.get(i).getObjective(0) == Double.MAX_VALUE / 100 ||
					pop.get(i).getObjective(1) == Double.MAX_VALUE / 100) {
			//if(record.contains(v)){
				index.add(i);
				
				/*if(removed.size() == 0 || pop.get(i).getObjective(target) > removed.get(removed.size() - 1).getObjective(target)) {
					removed.add(pop.get(i));
				} else {
					for (int k = 0; k < removed.size(); k++) {
						if (pop.get(i).getObjective(target) < removed.get(k).getObjective(target)) {
							removed.add(k, pop.get(i));
							break;
						}
					}
				}*/
				removed.add(pop.get(i));
				
				
				
			} else {
				record.add(v);
			}
			
		}
		
		/*Collections.sort(removed, new Comparator() {

			@Override
			public int compare(Object a, Object b) {
				Solution sa = (Solution)a;
				Solution sb = (Solution)b;
				
				if(sa.getObjective(target) < sb.getObjective(target)) {
					return -1;
				}
				if(sa.getObjective(target) > sb.getObjective(target)) {
					return 1;
				}
				
				return 0;
			}
			
			
		});*/
		
		for (int i = 0; i < pop.size(); i++) {
			

			boolean include = true;
			for (int term : index) {
				if(i == term) {
					include = false;
					break;
				}
						
			}
			
			if(include) {
				newSet.add(pop.get(i));
			}
		}
		
		// may mean all the one are invlaid
		if(newSet.size() == 0) {
			return pop;
		}
		
		return newSet;
	}
	
	/**
	 * This is used to find the knee point from a set of solutions
	 * 
	 * @param population
	 * @return
	 */
	public Solution kneeSelection(SolutionSet population_) {		
		int[] max_idx    = new int[problem_.getNumberOfObjectives()];
		double[] max_obj = new double[problem_.getNumberOfObjectives()];
		int populationSize_ = population_.size();
		// finding the extreme solution for f1
		for (int i = 0; i < populationSize_; i++) {
			for (int j = 0; j < problem_.getNumberOfObjectives(); j++) {
				// search the extreme solution for f1
				if (population_.get(i).getObjective(j) > max_obj[j]) {
					max_idx[j] = i;
					max_obj[j] = population_.get(i).getObjective(j);
				}
			}
		}

		if (max_idx[0] == max_idx[1])
			System.out.println("Watch out! Two equal extreme solutions cannot happen!");
		
		int maxIdx;
		double maxDist;
		double temp1 = (population_.get(max_idx[1]).getObjective(0) - population_.get(max_idx[0]).getObjective(0)) * 
				(population_.get(max_idx[0]).getObjective(1) - population_.get(0).getObjective(1)) - 
				(population_.get(max_idx[0]).getObjective(0) - population_.get(0).getObjective(0)) * 
				(population_.get(max_idx[1]).getObjective(1) - population_.get(max_idx[0]).getObjective(1));
		double temp2 = Math.pow(population_.get(max_idx[1]).getObjective(0) - population_.get(max_idx[0]).getObjective(0), 2.0) + 
				Math.pow(population_.get(max_idx[1]).getObjective(1) - population_.get(max_idx[0]).getObjective(1), 2.0);
		double constant = Math.sqrt(temp2);
		double tempDist = Math.abs(temp1) / constant;
		maxIdx  = 0;
		maxDist = tempDist;
		for (int i = 1; i < populationSize_; i++) {
			temp1 = (population_.get(max_idx[1]).getObjective(0) - population_.get(max_idx[0]).getObjective(0)) *
					(population_.get(max_idx[0]).getObjective(1) - population_.get(i).getObjective(1)) - 
					(population_.get(max_idx[0]).getObjective(0) - population_.get(i).getObjective(0)) * 
					(population_.get(max_idx[1]).getObjective(1) - population_.get(max_idx[0]).getObjective(1));
			tempDist = Math.abs(temp1) / constant;
			if (tempDist > maxDist) {
				maxIdx  = i;
				maxDist = tempDist;
			}
		}
		
		return population_.get(maxIdx);
	}
	
	private static boolean isSame(Solution a, Solution b, boolean removeInvalid) {
		
		
		if(removeInvalid) {
			if( a.getObjective(0) == Double.MAX_VALUE / 100 ||
					a.getObjective(1) == Double.MAX_VALUE / 100) {
				return true;
			}
		}
		
		
		
		try {
			
			//for(int i = 0; i < a.numberOfObjectives(); i++) {
			//	if(a.getObjective(i))
			//}
			
			
		String v1 = "";
		String v2 = "";
		for(int i = 0; i < a.getDecisionVariables().length; i++) {
			
				v1 += v1.equals("")? a.getDecisionVariables()[i].getValue() : ":" + a.getDecisionVariables()[i].getValue();
				v2 += v2.equals("")? b.getDecisionVariables()[i].getValue() : ":" + b.getDecisionVariables()[i].getValue();
				
		}
		
		return v1.equals(v2);
		} catch (JMException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		return false;
	}
	
	public static class MMOWeightAdaptor {
		
		public static double nondominatedProportion(HashMap<Double,Double> q_table) {
			
//			if(q_table.size() == 0) {
//				q_table.put(0.1, 1.0 / 8);
//				q_table.put(0.2, 1.0 / 8);
//				q_table.put(0.3, 1.0 / 8);
//				q_table.put(0.4, 1.0 / 8);
//				q_table.put(0.5, 1.0 / 8);
//				/*q_table.put(0.6, 1.0 / 8);
//				q_table.put(0.7, 1.0 / 8);
//				q_table.put(0.8, 1.0 / 8);
//				q_table.put(0.9, 1.0 / 8);*/
//			}
			
			
			double total = 0.0;
			for (Double d : q_table.keySet()) {
				total += q_table.get(d);
			}
			
			double[] array = new double[q_table.size()];
			double[] array_index = new double[] {0.3,0.4,0.5,0.6,0.7/*,0.6,0.7,0.8,0.9*/};
			for (int i = 0; i < array.length; i++) {
				array[i] = q_table.get(array_index[i]) / total;
				System.out.print("q table: " + array_index[i] + ", probability: " + (q_table.get(array_index[i]) / total) + "\n");
			}
			// +- 1.0/8
			double p = PseudoRandom.randDouble();
			double cumulativeProbability = 0.0;
			for (int i = 0; i < array.length; i++) {
			    cumulativeProbability += array[i];
			    if (p <= cumulativeProbability) {
			        return array_index[i];
			    }
			}
			
			return 0.1;
		}
		
		public static boolean isAdapt(int same_optimum_count, int current_meansurement, int budget) {
			double decay_point = 0.5;
			
			int offset = 1;//(int)Math.round(0.5 * budget);
			
			double slope = budget / current_meansurement;
			int d = same_optimum_count - offset;
			
			System.out.print("d : " + d + "\n");
			
			double ro = ((slope * slope) / (2 * Math.log(decay_point))) * -1.0;
			
			double t = 1 - Math.exp(((Math.pow(Math.max(0, d),2) / (2 * ro)) * -1.0));
			double p = PseudoRandom.randDouble();
			
			//System.out.print((p <= t ) + " : " + p +  " : " + t + "\n");
			
			return p <= t;
		}
		
	
		private static SolutionSet filterAndRank (Problem problem_, SolutionSet set, int populationSize)  throws JMException {
			int remain = populationSize;
			int index = 0;
			
			Ranking ranking = new Ranking(set);
			SolutionSet population = new SolutionSet();
			
			Set<String> record = new HashSet<String>();
			List<SolutionSet> removedRedundant = new ArrayList<SolutionSet>();
			
			SolutionSet front = ranking.getSubfront(index);
			Distance distance = new Distance();
			
			for (int k = 0; k < ranking.getNumberOfSubfronts(); k++) {
				removedRedundant.add(new SolutionSet());
			}
			
			while ((remain > 0)) {
				//Assign crowding distance to individuals
				distance.crowdingDistanceAssignment(front, problem_.getNumberOfObjectives());
				//Add the individuals of this front
				
				SolutionSet new_set = new SolutionSet();
				for (int k = 0; k < front.size(); k++) {
					
					String v = "";
					for (int j = 0; j < front.get(k).getDecisionVariables().length; j++) {
						v += front.get(k).getDecisionVariables()[j].getValue()+",";
					}
					
					if(record.contains(v)) {
						removedRedundant.get(index).add(front.get(k));
					} else {
						new_set.add(front.get(k));
						
						record.add(v);
					}
					
				
				} // for
				
	            if(remain < new_set.size()) {
					break;
				}
				
				for (int k = 0; k < new_set.size(); k++) {
					population.add(new_set.get(k));
				}
				
				

				//Decrement remain
				remain = remain - new_set.size();
				
				

				//Obtain the next front
				index++;
				// for remove redundant ***********
				if(index >= ranking.getNumberOfSubfronts()) {
					break;
				}
				// for remove redundant ***********
				if (remain > 0) {
				     front = ranking.getSubfront(index);
					
				} // if  
				
				
			} // while
			
			if(index >= ranking.getNumberOfSubfronts() && remain > 0) {
				
				for (int i = 0; i < removedRedundant.size(); i++) {
					
					for (int k = 0; k < removedRedundant.get(i).size(); k++) {
					
						removedRedundant.get(i).sort(new CrowdingComparator());
						population.add(removedRedundant.get(i).get(k));
						remain--;
						if(remain == 0) {
							break;
						}
					}
					
					if(remain == 0) {
						break;
					}
					
				}
				
			} else {

				// Remain is less than front(index).size, insert only the best one
				if (remain > 0) {  // front contains individuals to insert     
					record.clear();
					distance.crowdingDistanceAssignment(front, problem_.getNumberOfObjectives());
					front.sort(new CrowdingComparator());
					for (int k = 0; k < front.size(); k++) {
						
						if(remain <= 0) {
							break;
						}
						
						String v = "";
						for (int j = 0; j < front.get(k).getDecisionVariables().length; j++) {
							v += front.get(k).getDecisionVariables()[j].getValue()+",";
						}
						
						if(record.contains(v)) {
							//removedRedundant.get(index).add(front.get(k));
						} else {
							remain--;
							population.add(front.get(k));
							
							record.add(v);
						}
						
						
						
					} // for
					
					if(remain > 0) {
						for (int k = 0; k < remain; k++) {
							removedRedundant.get(index).sort(new CrowdingComparator());
							population.add(removedRedundant.get(index).get(k));
						}
					}

					remain = 0;
				} // if 
			}
			
			return population;
			 
		}
		
		public static double testProportion(SASSolutionInstantiator factory , SolutionSet set, 
				double w) {
			
			
			
			SolutionSet new_set = new SolutionSet();
			
			if(SASAlgorithmAdaptor.isToFilterRedundantSolution) {

				Set<String> record = new HashSet<String>();
				for (int i = 0; i < set.size(); i++) {
					
					
					
					String v = "";
					for (int j = 0; j < set.get(i).getDecisionVariables().length; j++) {
						try {
							v += set.get(i).getDecisionVariables()[j].getValue()+",";
						} catch (JMException e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
						}
					}
					
					
					if(record.contains(v)) {
						
					} else {
						new_set.add(set.get(i));
						record.add(v);
					}
					
					
				}
				/*try {
					new_set = filterAndRank(problem_, set, population_size);
				} catch (JMException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}*/
				
			} else {
				new_set = set;
			}
			
			SolutionSet newSet = factory.fuzzilize(new_set, w);
			
			Ranking ranking = new Ranking(newSet);
			double nondomainted_count = ranking.getSubfront(0).size();
			//int nondominated_threshold = (int) Math.round((int)new_set.size() * proportion);
			
			return nondomainted_count / (double)new_set.size();
			
		}
		
		public static double smartGetWeight(Problem problem_, SASSolutionInstantiator factory , SolutionSet set, 
				double w, int previous_nondomainted_count /*initialize with <=0*/, int population_size, double step, double proportion) {
			
			
			
			SolutionSet new_set = new SolutionSet();
			
			if(SASAlgorithmAdaptor.isToFilterRedundantSolution) {

				Set<String> record = new HashSet<String>();
				for (int i = 0; i < set.size(); i++) {
					
					
					
					String v = "";
					for (int j = 0; j < set.get(i).getDecisionVariables().length; j++) {
						try {
							v += set.get(i).getDecisionVariables()[j].getValue()+",";
						} catch (JMException e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
						}
					}
					
					
					if(record.contains(v)) {
						
					} else {
						new_set.add(set.get(i));
						record.add(v);
					}
					
					
				}
				/*try {
					new_set = filterAndRank(problem_, set, population_size);
				} catch (JMException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}*/
				
			} else {
				new_set = set;
			}
			
			/*try {
				MMOWeightAdaptor.printSolutions(new_set);
			} catch (JMException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}*/
			
			//Ranking ranking = new Ranking(new_set);
			//double nondomainted_count = ranking.getSubfront(0).size();
			//System.out.print("Percentage: " + (double)(nondomainted_count / new_set.size()) + "\n");
			//return 1.0;
			
			// for testing
			//w = 0;
			//proportion = 1.0;
			//step = 0.0001;
			// for testing
			return navielyLoopGetWeight(factory, new_set, w, previous_nondomainted_count, (int)new_set.size(),step, proportion);
		}
		
		public static double navielyLoopGetWeight(SASSolutionInstantiator factory , SolutionSet set, 
				double w_input, int previous_nondomainted_count /*initialize with <=0*/, int population_size, double step_input, double proportion) {
			
			boolean c = true;
			double final_w = 0.0;
			double w = w_input;
			double step = step_input;
			do {
				

				SolutionSet newSet = factory.fuzzilize(set, w);
				
				Ranking ranking = new Ranking(newSet);
				int nondomainted_count = ranking.getSubfront(0).size();
				int nondominated_threshold = (int) Math.round(population_size * proportion);
				
				System.out.print(newSet.size() +" nondomainted_count: " + nondomainted_count  + ", threshold: " + nondominated_threshold +  " = " + w +"\n");
				//if(w % 0.1 == 0)
				//double t = nondomainted_count;
				//double s = nondomainted_count / 20.0;
				//System.out.print("(" + w + "," + s + ")\n");
				
				
				
				if(nondomainted_count < nondominated_threshold) {
					
					//return 1.0;
					//step = 10 * step;
					//if(w + step > Double.MAX_VALUE) {
					//if(w + step > 1.0) {
					if(w /*+ step*/ > 1000) {
						//final_w = 1.0;
						//final_w = w;
						final_w = 1000;
						break;
					}
					
					if(w >= 0.1 && step == 0.0001) {
						w = 0.1;
						step = 0.1;
					}
				}
				
				if(nondomainted_count > nondominated_threshold) {
					
					if(w < 0.1 && step == 0.1) {
						step = 0.0001;
						w = 0.1 - step;
					}
					
				}
				
				/*if(w < 0) {
					final_w = Double.MIN_VALUE;
					break;
				}*/
				
				if (previous_nondomainted_count > 0) {
					
					if(nondomainted_count == nondominated_threshold) {
				
						final_w = w;
						break;
					}
					
					if(previous_nondomainted_count < nondominated_threshold && nondomainted_count > nondominated_threshold) {
						
						int a = nondominated_threshold - previous_nondomainted_count;
						int b = nondomainted_count - nondominated_threshold;
				
						final_w = a <= b? w - step : w;
						break;
					}
					
	                if(previous_nondomainted_count > nondominated_threshold && nondomainted_count < nondominated_threshold) {
						
	                	int a = previous_nondomainted_count - nondominated_threshold;
	                	int b = nondominated_threshold - nondomainted_count;
						
	                	final_w = a <= b? w + step : w;
	                	break;
					}
	                
	                
	                if(nondomainted_count < nondominated_threshold) {
	                	previous_nondomainted_count = nondomainted_count;
	                	if(w + step > Double.MAX_VALUE) {
	                		final_w = w;
	                		break;
	                	}
	                	w = w + step;
	                }
	                
	                if(nondomainted_count > nondominated_threshold) {
	                	previous_nondomainted_count = nondomainted_count;
	                	if(w - step < Double.MIN_VALUE) {
	                		final_w = w;
	                		break;
	                	}
	                	w = w - step;
	                }
				} else {
					
					
					
					if(nondomainted_count < nondominated_threshold) {
						previous_nondomainted_count = nondomainted_count;
						if(w + step > Double.MAX_VALUE) {
	                		final_w = w;
	                		break;
	                	}
						w = w + step;
	                } else if(nondomainted_count > nondominated_threshold) {
	                	previous_nondomainted_count = nondomainted_count;
	                	if(w - step < Double.MIN_VALUE) {
	                		final_w = w;
	                		break;
	                	}
	                	w = w - step;
	                	
	                } else {
						final_w = w;
						break;
	                }
				}
				
			} while(c);
			
			return final_w;
		}
		
		public static double timesLoopGetWeight(SASSolutionInstantiator factory , SolutionSet set, 
				double w_input, int previous_nondomainted_count /*initialize with <=0*/, int population_size, double step) {
			
			boolean c = true;
			double final_w = 0.0;
			double w = w_input;
			
			do {
				

				SolutionSet newSet = factory.fuzzilize(set, w);
				
				Ranking ranking = new Ranking(newSet);
				int nondomainted_count = ranking.getSubfront(0).size();
				int nondominated_threshold = population_size;
				
				System.out.print("nondomainted_count: " + nondomainted_count  + ", threshold: " + population_size +  " = " + w +"\n");
				
				
				
				
				if(w < 0) {
					final_w = Double.MIN_VALUE;
					break;
				}
				
				if (previous_nondomainted_count > 0) {
					
					if(nondomainted_count == nondominated_threshold) {
				
						final_w = w;
						break;
					}
					
					if(previous_nondomainted_count < nondominated_threshold && nondomainted_count > nondominated_threshold) {
						
						int a = nondominated_threshold - previous_nondomainted_count;
						int b = nondomainted_count - nondominated_threshold;
				
						final_w = a <= b? w / step : w;
						break;
					}
					
	                if(previous_nondomainted_count > nondominated_threshold && nondomainted_count < nondominated_threshold) {
						
	                	int a = previous_nondomainted_count - nondominated_threshold;
	                	int b = nondominated_threshold - nondomainted_count;
						
	                	final_w = a <= b? w * step : w;
	                	break;
					}
	                
	                
	                if(nondomainted_count < nondominated_threshold) {
	                	previous_nondomainted_count = nondomainted_count;
	                	if(w * step > Double.MAX_VALUE) {
	                		final_w = w;
	                		break;
	                	}
	                	w = w * step;
	                }
	                
	                if(nondomainted_count > nondominated_threshold) {
	                	previous_nondomainted_count = nondomainted_count;
	                	if(w / step < Double.MIN_VALUE) {
	                		final_w = w;
	                		break;
	                	}
	                	w = w / step;
	                }
				} else {
					
					
					
					if(nondomainted_count < nondominated_threshold) {
						previous_nondomainted_count = nondomainted_count;
						if(w * step > Double.MAX_VALUE) {
	                		final_w = w;
	                		break;
	                	}
						w = w * step;
	                } else if(nondomainted_count > nondominated_threshold) {
	                	previous_nondomainted_count = nondomainted_count;
	                	if(w / step < Double.MIN_VALUE) {
	                		final_w = w;
	                		break;
	                	}
	                	w = w / step;
	                	
	                } else {
						final_w = w;
						break;
	                }
				}
				
			} while(c);
			
			return final_w;
		}
		
		public static double navielyGetWeight(SASSolutionInstantiator factory , SolutionSet set, 
				double w, int previous_nondomainted_count /*initialize with <=0*/, int population_size, double step) {
			
			
			
			SolutionSet newSet = factory.fuzzilize(set, w);
			
			Ranking ranking = new Ranking(newSet);
			int nondomainted_count = ranking.getSubfront(0).size();
			int nondominated_threshold = population_size;
			
			System.out.print("nondomainted_count: " + nondomainted_count  + ", threshold: " + population_size +  " = " + w +"\n");
			
			
			
			if(w > 1.0 && nondomainted_count < nondominated_threshold) {
				
				//return 1.0;
				step = 10 * step;
				if(w + step > Double.MAX_VALUE) {
					return Double.MAX_VALUE;
				}
			}
			
			if(w < 0) {
				return Double.MIN_VALUE;
			}
			
		
			
			if (previous_nondomainted_count > 0) {
				
				if(nondomainted_count == nondominated_threshold) {
					return w;
				}
				
				if(previous_nondomainted_count < nondominated_threshold && nondomainted_count > nondominated_threshold) {
					
					int a = nondominated_threshold - previous_nondomainted_count;
					int b = nondomainted_count - nondominated_threshold;
					
					return a < b? w - step : w;
					
				}
				
                if(previous_nondomainted_count > nondominated_threshold && nondomainted_count < nondominated_threshold) {
					
                	int a = previous_nondomainted_count - nondominated_threshold;
                	int b = nondominated_threshold - nondomainted_count;
					
					
					return a < b? w + step : w;
				}
                
                
                if(nondomainted_count < nondominated_threshold) {
                	return navielyGetWeight(factory, set, w + step, nondomainted_count, population_size, step);
                }
                
                if(nondomainted_count > nondominated_threshold) {
                	return navielyGetWeight(factory, set, w - step, nondomainted_count, population_size, step);
                }
			} else {
				
				
				
				if(nondomainted_count < nondominated_threshold) {
                	return navielyGetWeight(factory, set, w + step, nondomainted_count, population_size, step);
                } else if(nondomainted_count > nondominated_threshold) {
                	
                	return navielyGetWeight(factory, set, w - step, nondomainted_count, population_size, step);
                } else {
                	return w;
                }
			}
			return w;
		}
		
		public static double getWeight(SolutionSet set, int population_size) {
			
			Ranking ranking = new Ranking(set);
			int nondomainted_count = ranking.getSubfront(0).size();
			
			
			
			if(nondomainted_count ==  population_size) {
				return Double.NaN;
			}
			
			boolean higher = nondomainted_count <  population_size;
			//int distance  = Math.abs(nondomainted_count - population_size);
			
			double alpha = 0.0;
			double beta = 0.0;
			
			double points = 1000.0;
			double gap = 1.0 / points;
			
			
			
			if (higher) {
				alpha = population_size / nondomainted_count;
				beta = 1;
				
				//alpha = alpha / population_size;
				//beta = beta / population_size;
				
				BetaDistribution bd = new BetaDistribution(alpha,beta); 
				
				double p = PseudoRandom.randDouble();
				double cumulativeProbability = 0.0;
				for (int i = 0; i < points; i++) {
					double w = gap * (i+1);
				    cumulativeProbability = bd.cumulativeProbability(w);
				    if (p <= cumulativeProbability) {
				    	//System.out.print(w + " : " + cumulativeProbability + "\n");
				        return w;
				    }
				}
			} else {
				alpha = 1;
				if (2*population_size - nondomainted_count == 0) { 
					beta = population_size;
				} else {
					beta = population_size / (2*population_size - nondomainted_count);
				}
				//alpha = alpha / population_size;
				//beta = beta / population_size;
				
				BetaDistribution bd = new BetaDistribution(alpha,beta); 
				
				double p = PseudoRandom.randDouble();
				double cumulativeProbability = 0.0;
				for (int i = (int)points - 1; i > 0; i--) {
					double w = gap * (i+1);
				    cumulativeProbability = 1.0 - bd.cumulativeProbability(w);
				    if (p <= cumulativeProbability) {
				    	//System.out.print(w + " : " + cumulativeProbability + "\n");
				        return w;
				    }
				}
			}
			return Double.NaN;
		}
		
        public static double getWeight(int nondomainted_count, int population_size) {
			
		
			
			
			
			if(nondomainted_count ==  population_size) {
				return Double.NaN;
			}
			
			boolean higher = nondomainted_count <  population_size;
			//int distance  = Math.abs(nondomainted_count - population_size);
			
			double alpha = 0.0;
			double beta = 0.0;
			
			double points = 1000.0;
			double gap = 1.0 / points;
			
			//System.out.print(higher);
			
			if (higher) {
				alpha = population_size / nondomainted_count;
				beta = 1;
				
				//alpha = alpha / population_size;
				//beta = beta / population_size;
				
				BetaDistribution bd = new BetaDistribution(alpha,beta); 
				
				double p = PseudoRandom.randDouble();
				double cumulativeProbability = 0.0;
				for (int i = 0; i < points; i++) {
					double w = gap * (i+1);
				    cumulativeProbability = bd.cumulativeProbability(w);
				    System.out.print(p + " : " + w + " : " + cumulativeProbability + "\n");
				    if (p <= cumulativeProbability) {
				    	//System.out.print(p + " : " + w + " : " + cumulativeProbability + "\n");
				        return w;
				    }
				}
			} else {
				alpha = 1;
				
				if (2*population_size - nondomainted_count == 0) { 
					beta = population_size;
				} else {
					beta = population_size / (2*population_size - nondomainted_count);
				}
				
				//alpha = alpha / population_size;
				//beta = beta / population_size;
				
				BetaDistribution bd = new BetaDistribution(alpha,beta); 
				
				double p = PseudoRandom.randDouble();
				double cumulativeProbability = 0.0;
				for (int i = (int)points - 1; i > 0; i--) {
					double w = gap * (i+1);
				    cumulativeProbability = 1.0 - bd.cumulativeProbability(w);
				  System.out.print(p + " : " + w + " : " + cumulativeProbability + "\n");
				    if (p <= cumulativeProbability) {
				    	//System.out.print(p + " : " + w + " : " + cumulativeProbability + "\n");
				        return w;
				    }
				}
			}
			return Double.NaN;
		}
        
        private static void printSolutions(SolutionSet pop) throws JMException {
    		System.out.print("-----------------\n");
    		for (int i = 0; i < pop.size(); i++) {
    			
    		   
    			String v = "";
    			for (int j = 0; j < pop.get(i).getDecisionVariables().length; j++) {
    				v += pop.get(i).getDecisionVariables()[j].getValue()+",";
    			}
    			
    			v += "=" + pop.get(i).getObjective(0) + ":" + pop.get(i).getObjective(1);
    			
    			System.out.print(v + "\n");
    			
    		}
    		System.out.print("-----------------\n");
    	}
	}
	
	
	public static class RemoveRedundantAdaptor {
		
		public static Solution matingRemoval(SolutionSet pop) {
			
			Ranking ranking = new Ranking(pop);
			
			double total = 0.0;
			int[] array = new int[ranking.getNumberOfSubfronts()];
			double[] prob = new double[ranking.getNumberOfSubfronts()];
			for (int i = 0; i < ranking.getNumberOfSubfronts(); i++) {
				array[i] = ranking.getSubfront(i).size();
				total += array[i];
			}
			for (int i = 0; i < ranking.getNumberOfSubfronts(); i++) {
				prob[i] = array[i] / total;
				
			}
			
			
			double p = PseudoRandom.randDouble();
			double cumulativeProbability = 0.0;
			
			SolutionSet front = null;
			
			for (int i = 0; i < prob.length; i++) {
			    cumulativeProbability += prob[i];
			    //System.out.print("cumulativeProbability: " + cumulativeProbability + "\n");
			    if (p <= cumulativeProbability) {
			    	front = ranking.getSubfront(i);
			    	break;
			    }
			}
			
			//System.out.print("front: " + ranking.getNumberOfSubfronts() + "\n");
            SolutionSet new_set = new SolutionSet();
			
		

			for (int i = 0; i < front.size(); i++) {
				boolean same = false;

				for (int j = 0; j < new_set.size(); j++) {
					same = isSame(front.get(i), new_set.get(j), false);
					if (same) {
						break;
					}
				}

				if (!same) {
					new_set.add(front.get(i));
				}
			}

			int r = PseudoRandom.randInt(0, new_set.size()-1);
			
			
			return new_set.get(r);
		}
	}
	
	public static void main (String[] arg) {
		
		int current_meansurement = 10; 
		int budget = 100;
		int population_size = 30;
		
		for (int i = 1; i <= 100; i++) {
			current_meansurement += i;
			
			//MMOWeightAdaptor.isAdapt(7, current_meansurement, budget);
			
		}
		//System.out.print(3%10);
		for (int i = 5; i <= 5; i++) {
			
			//MMOWeightAdaptor.getWeight(i, population_size);
			
		}
	}
} // NSGA-II
