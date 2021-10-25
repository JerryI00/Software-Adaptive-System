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

package jmetal.metaheuristics.msga;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

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
 * @author taochen
 *
 */
public class MNSGAII_SAS extends Algorithm {

	private SASSolutionInstantiator factory = null;
	
	private SASValidityAndInvalidityCoEvolver vandInvCoEvolver = null;
	private Seeder seeder = null;
	SolutionSet population_;
	boolean logOnce = false;
	private boolean isChange = true;

	private Set<String> full_set;
	/**
	 * Constructor
	 * @param problem Problem to solve
	 */
	public MNSGAII_SAS(Problem problem) {
		super (problem) ;
	} // NSGAII


  	/**
  	 * Constructor
  	 * @param problem Problem to solve
  	 */
	public MNSGAII_SAS(Problem problem, SASSolutionInstantiator factory) {
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
				
					parents[0] = (Solution) selectionOperator.execute(population);
					parents[1] = (Solution) selectionOperator.execute(population);
				
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
			
			if (SASAlgorithmAdaptor.isToFilterRedundantSolution) {

				int target = 0;
				Set<String> record = new HashSet<String>();
				List<Solution> removed = new ArrayList<Solution>();
				if (SASAlgorithmAdaptor.isFuzzy) {
					old_population = filter(old_population, record, removed, target);
				} else {
					population = filter(population, record, removed,  target);
				}
				
				offspringPopulation = filter(offspringPopulation, record, removed,  target);
				
				int size = 0;
				if (SASAlgorithmAdaptor.isFuzzy) {
					size = old_population.size();
				} else {
					size = population.size();
				}
				
				if (SASAlgorithmAdaptor.isFuzzy) {
					printSolutions(old_population);
				} else {
					printSolutions(population);
				}
				printSolutions(offspringPopulation);
			
				
				size = size + offspringPopulation.size();
				
				System.out.print("Unique: " + size + "\n");
				
				if (size < populationSize) {
					int l = populationSize - size;
					for (int i = 0; i < l; i++) {
						offspringPopulation.add(removed.get(i));
					}
				} else {
					
				}
											
			}
			
		
			
			// Create the solutionSet union of solutionSet and offSpring			
			if(SASAlgorithmAdaptor.isFuzzy) {			
				union = ((SolutionSet) old_population).union(offspringPopulation);
				//twoNeighbour(union);
				smallestDistance(union);
				// Max and min difference
				/*for(int i = 0; i < union.size(); i++) {
					Solution s = union.get(i);
					double max = Double.NEGATIVE_INFINITY;
					double min = Double.MAX_VALUE;
					if(s.getObjective(0) == Double.MAX_VALUE / 100) {
						s.setObjective(1,Double.MAX_VALUE / 100);
						continue;
					}
					
					
					for(int j = 0; j < union.size();j++) {
						
						
						
						if(i != j) {
							
							if(union.get(j).getObjective(0) == Double.MAX_VALUE / 100) {
								continue;
							}
							
							
							if(Math.abs(union.get(j).getObjective(0) - s.getObjective(0)) > max) {
								max = union.get(j).getObjective(0);
							}
							
							if(Math.abs(union.get(j).getObjective(0) - s.getObjective(0)) < min) {
								min = union.get(j).getObjective(0);
							}
							
							
						}
						
					}
					
					if( max == Double.NEGATIVE_INFINITY && min == Double.MAX_VALUE) {
						s.setObjective(1,s.getObjective(0));
					} else {
						s.setObjective(1, (max + min) / 2.0);
					}
					
					
					//System.out.print(total / (union.size()-1) + "****\n");
				}*/
				//System.out.print(total / (union.size()-1) + "****\n");
				
				// normal novelty score
				/*for(int i = 0; i < union.size(); i++) {
					Solution s = union.get(i);
					double total = 0;
					for(int j = 0; j < union.size();j++) {
						if(i != j) {
							Variable[] a = s.getDecisionVariables();
							Variable[] b = union.get(j).getDecisionVariables();
							double d = distance(a,b);
							total += d;
						}
						
					}
					s.setObjective(1, total / (union.size()-1));
					//System.out.print(total / (union.size()-1) + "****\n");
				}*/
				
				// $$$$$$$$$$$$$$$$$$$
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
				union = factory.fuzzilize(union);//((double) measurement - 1) / (EAConfigure.getInstance().measurement - 1)
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
			
			if(measurement >= EAConfigure.getInstance().measurement/2 && !isChange) {
				isChange = true;
				factory.fuzzilize(union, 0);
				
				((SASSolution)old_population.get(0)).resetNormalizationBounds(0);
				((SASSolution)old_population.get(0)).resetNormalizationBounds(1);
				/*for(int i = 0; i < union.size(); i++) {
					((SASSolution)union.get(i)).updateNormalizationBounds(new double[] {union.get(i).getObjective(0),
							union.get(i).getObjective(1)});
				}*/
				
				for(int i = 0; i < old_population.size(); i++) {
					problem_.evaluate(old_population.get(i));
					//((SASSolution)old_population.get(i)).updateNormalizationBounds(new double[] {old_population.get(i).getObjective(0),
					//		old_population.get(i).getObjective(1)});
				}
				
				//twoNeighbour(old_population);
				smallestDistance(old_population);
				
				//measurement += old_population.size();
				SolutionSet ss = factory.fuzzilize(old_population);
				population.clear();
				for(int i = 0; i < ss.size(); i++) {
					population.add(ss.get(i));
				}
				
				
				System.out.print("Finish changing...\n");
			}
			
			
			
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
		
		// Return the first non-dominated front
//		Ranking ranking = new Ranking(population);
//		return ranking.getSubfront(0);
		return population;
	} // execute
	
	
	private double distance(Variable[] a, Variable[] b) throws JMException{
	        double diff_square_sum = 0.0;
	        for (int i = 0; i < a.length; i++) {
	            diff_square_sum += (a[i].getValue() - b[i].getValue() ) * (a[i].getValue()  - b[i].getValue() );
	        }
	        return Math.sqrt(diff_square_sum);
	}
	
	void smallestDistance(SolutionSet union) throws JMException{
		for(int i = 0; i < union.size(); i++) {
			
			if(union.get(i).getObjective(0) == Double.MAX_VALUE / 100) {
				union.get(i).setObjective(1,Double.MAX_VALUE / 100);
				continue;
			}
			
			final Solution s = union.get(i);
			final Variable[] baseline = s.getDecisionVariables();
			List<Solution> list = new ArrayList<Solution>();
			for(int j = 0; j < union.size();j++) {
				
				if(union.get(j).getObjective(0) == Double.MAX_VALUE / 100) {
					continue;
				}
				
				if(i != j) {
					list.add(union.get(j));
				}
				
			}
			
			Collections.sort(list, new Comparator<Solution>() {

				@Override
				public int compare(Solution arg0, Solution arg1) {
					Variable[] a = arg0.getDecisionVariables();
					Variable[] b = arg1.getDecisionVariables();
					
					double d1 = 0.0;
					double d2 = 0.0;
					try {
						d1 = distance(baseline,a);
						d2 = distance(baseline,b);
					} catch (JMException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
					
					/*if(d1 == 0 && d2 != 0) {
						return 1;
					} else if(d1 != 0 && d2 == 0) {
						return -1;
					} else if(d1 == 0 && d2 == 0) {
						return 0;
					} else*/ 
						
					if(d1 < d2) {
						return -1;
					} else if(d2 > d1) {
						return 1;
					} else {
						if(Math.abs(s.getObjective(0)-arg0.getObjective(0)) > Math.abs(s.getObjective(0)-arg1.getObjective(0))) {
							return -1;
						} else if(Math.abs(s.getObjective(0)-arg0.getObjective(0)) < Math.abs(s.getObjective(0)-arg1.getObjective(0))) {
							return 1;
						}
					}
					return 0;
				}
				
			});
			
			
			
			s.setObjective(1, list.size() == 0? Double.MAX_VALUE / 100 : list.get(0).getObjective(0));
			//System.out.print(total / (union.size()-1) + "****\n");
		}
	}
	
	void twoNeighbour(SolutionSet union) {
		// assume only interested in objective 0
		// based on the average of closed objectives
		List<Solution> list = new ArrayList<Solution>();
		
		for(int i = 0; i < union.size(); i++) {
			list.add(union.get(i));
		}
		
		Collections.sort(list, new Comparator<Solution>() {

			@Override
			public int compare(Solution arg0, Solution arg1) {
				if(arg0.getObjective(0) < arg1.getObjective(0)) {
					return -1;
				} else if(arg0.getObjective(0) > arg1.getObjective(0)) {
					return 1;
				}
				return 0;
			}
			
		});
		
		for(int i = 0; i < list.size(); i++) {
			
			if(list.get(i).getObjective(0) == Double.MAX_VALUE / 100) {
				list.get(i).setObjective(1,Double.MAX_VALUE / 100);
				continue;
			}
			
			if(i == 0) {
				double a = list.get(list.size() - 1).getObjective(0) == Double.MAX_VALUE / 100? 0 : list.get(list.size() - 1).getObjective(0);
				double b = list.get(i+1).getObjective(0) == Double.MAX_VALUE / 100? 0 : list.get(i+1).getObjective(0);
				double n = a == 0 || b == 0? 1.0 : 2.0;
				list.get(i).setObjective(1, (a + b)/n);
			} else if(i == list.size() - 1) {
				double a = list.get(i - 1).getObjective(0) == Double.MAX_VALUE / 100? 0 : list.get(i - 1).getObjective(0);
				double b = list.get(0).getObjective(0) == Double.MAX_VALUE / 100? 0 : list.get(0).getObjective(0);
				double n = a == 0 || b == 0? 1.0 : 2.0;
				list.get(i).setObjective(1, (a + b)/n);		
			} else {
				double a = list.get(i - 1).getObjective(0) == Double.MAX_VALUE / 100? 0 : list.get(i - 1).getObjective(0);
				double b = list.get(i+1).getObjective(0) == Double.MAX_VALUE / 100? 0 : list.get(i+1).getObjective(0);
				double n = a == 0 || b == 0? 1.0 : 2.0;
				list.get(i).setObjective(1, (a + b)/n);
			}
		}
	}
	
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
	
	private SolutionSet filter(SolutionSet pop, Set<String> record, List<Solution> removed, int target) throws JMException {
		SolutionSet newSet = new SolutionSet();
		
		
		List<Integer> index = new ArrayList<Integer>();
		
		
		for (int i = 0; i < pop.size(); i++) {
			
			String v = "";
			for (int j = 0; j < pop.get(i).getDecisionVariables().length; j++) {
				v += pop.get(i).getDecisionVariables()[j].getValue()+",";
			}
			
			if(record.contains(v)) {
				index.add(i);
				
				for (int k = 0; k < removed.size(); k++) {
					if (pop.get(i).getObjective(target) < removed.get(k).getObjective(target)) {
						removed.add(k, pop.get(i));
						break;
					}
				}
				
			} else {
				record.add(v);
			}
			
		}
		
		
		
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
} // NSGA-II
