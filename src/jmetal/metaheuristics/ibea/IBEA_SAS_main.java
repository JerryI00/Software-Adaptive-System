//  IBEA_main.java
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

package jmetal.metaheuristics.ibea;

import jmetal.core.Algorithm;
import jmetal.core.Operator;
import jmetal.core.Problem;
import jmetal.core.Solution;
import jmetal.core.SolutionSet;
import jmetal.operators.crossover.CrossoverFactory;
import jmetal.operators.mutation.MutationFactory;
import jmetal.operators.selection.BinaryTournament;
import jmetal.operators.selection.SelectionFactory;
import jmetal.problems.ZDT.*;
import jmetal.problems.DTLZ.*;
import jmetal.problems.WFG.*;
import jmetal.problems.LZ09.*;
import jmetal.problems.cec2009Competition.*;
import jmetal.util.Configuration;
import jmetal.util.JMException;
import jmetal.util.PseudoRandom;
import jmetal.util.comparators.FitnessComparator;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.StringTokenizer;
import java.util.logging.FileHandler;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.femosaa.core.EAConfigure;
import org.femosaa.core.SAS;
import org.femosaa.core.SASAlgorithmAdaptor;
import org.femosaa.core.SASProblemFactory;
import org.femosaa.core.SASSolutionInstantiator;
import org.femosaa.invalid.SASValidityAndInvalidityCoEvolver;
import org.femosaa.seed.FixedSeeder;
import org.femosaa.seed.NewSeeder;

/**
 * Class for configuring and running the DENSEA algorithm
 */
public class IBEA_SAS_main extends SASAlgorithmAdaptor{
	public static Logger logger_; // Logger object
	public static FileHandler fileHandler_; // FileHandler object
	protected Problem problem; // The problem to solve
	protected Algorithm algorithm; // The algorithm to use
	Operator crossover; // Crossover operator
	Operator mutation; // Mutation operator
	Operator selection; // Selection operator

	/**
	 * @param args
	 *            Command line arguments.
	 * @throws JMException
	 * @throws IOException
	 * @throws SecurityException
	 *             Usage: three choices -
	 *             jmetal.metaheuristics.nsgaII.NSGAII_main -
	 *             jmetal.metaheuristics.nsgaII.NSGAII_main problemName -
	 *             jmetal.metaheuristics.nsgaII.NSGAII_main problemName
	 *             paretoFrontFile
	 */
	public static void main(String[] args) throws JMException, IOException,
			ClassNotFoundException {
		Problem problem; // The problem to solve
		Algorithm algorithm; // The algorithm to use
		Operator crossover; // Crossover operator
		Operator mutation; // Mutation operator
		Operator selection; // Selection operator

		HashMap parameters; // Operator parameters

		String problemName = "";
		String referenceParetoFront = "";
		int runIndex = 0; 
		int n        = 0;
		int iter     = 0;
		double[] ref_point = {};

		// Logger object and file to store log messages
		logger_ = Configuration.logger_;
		fileHandler_ = new FileHandler("IBEA.log");
		logger_.addHandler(fileHandler_);

		if (args.length == 1) {
			Object[] params = { "Real" };
			problem = (new SASProblemFactory()).getProblem(args[0], params);
		}
		else {
			if (args.length == 2) {
				Object[] params = {"Real"};
				problem = (new SASProblemFactory()).getProblem(args[0], params);
			}
			else { // Default problem
//				problem = new ZDT1("Real");
//				problem = new DTLZ2("Real");
				problem = new SAS("IntSolutionType", null, null, 2, 0);
			} // else
		}

		algorithm = new IBEA_SAS(problem, null);

		// Algorithm parameters
		int popsize = 331;
		int generations = 1000;
		algorithm.setInputParameter("populationSize", popsize);
		algorithm.setInputParameter("maxEvaluations", popsize * generations);

		algorithm.setInputParameter("dataDirectory", "weight");

		// Crossover operator
		int tag = 2;
		if (tag == 1) {
			parameters = new HashMap();
			parameters.put("CR", 0.5);
			parameters.put("F", 0.5);
			crossover = CrossoverFactory.getCrossoverOperator("DifferentialEvolutionCrossover", parameters);
		} else {
			parameters = new HashMap();
			parameters.put("probability", 0.9);
			parameters.put("distributionIndex", 20.0);
			// This needs to change in testing.
			parameters.put("jmetal.metaheuristics.moead.SASSolutionInstantiator", null);
			crossover = CrossoverFactory.getCrossoverOperator("UniformCrossoverSAS", parameters);
		}
		
		// Mutation operator
		parameters = new HashMap();
		parameters.put("probability", 1.0 / problem.getNumberOfVariables());
		parameters.put("distributionIndex", 20.0);
		mutation = MutationFactory.getMutationOperator("BitFlipMutation", parameters);

		/* Selection Operator (need check!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!) */
		//selection = SelectionFactory.getSelectionOperator("BinaryTournament2", parameters);
		parameters = new HashMap();
		parameters.put("comparator", new FitnessComparator());
		selection = new BinaryTournament(parameters);

		algorithm.addOperator("crossover", crossover);
		algorithm.addOperator("mutation", mutation);
		algorithm.addOperator("selection", selection);

		String curDir = System.getProperty("user.dir");
		String str 	  = curDir + "/" + problem.getName() + "M" + problem.getNumberOfObjectives();

		File dir = new File(str);
		if (Utils.deleteFolder(dir)) {
			System.out.println("Folders are deleted!");
		} else
			System.out.println("Folders can NOT be deleted!");
		Utils.createFolder(str);
		
		String str1 = "FUN";
		String str2;
		String str3 = "VAR";
		String str4;
		for (int i = 0; i < 1; i++) {
			str2 = str1 + Integer.toString(i);
			str4 = str3 + Integer.toString(i);
			// Execute the Algorithm
			long initTime = System.currentTimeMillis();
			System.out.println("The " + i + " run");
			SolutionSet population = algorithm.execute();
			long estimatedTime = System.currentTimeMillis() - initTime;

			// Result messages
			logger_.info("Total execution time: " + estimatedTime + "ms");
			logger_.info("Variables values have been writen to file VAR");
//			population.printVariablesToFile(str4);
			logger_.info("Objectives values have been writen to file FUN");
			population.printObjectivesToFile(curDir + "/" + problem.getName()
					+ "M" + problem.getNumberOfObjectives() + "/" + str2); 
		}
	} // main
		
	@SuppressWarnings({ "rawtypes", "unchecked" })
	protected SolutionSet findParetoFront(SASSolutionInstantiator factory, int[][] vars,  int numberOfObjectives_, int numberOfConstraints_) throws JMException,
			SecurityException, IOException, ClassNotFoundException {
	
		HashMap parameters; // Operator parameters
		Operator selection; // Selection operator

		// Logger object and file to store log messages
		if(SAS.isTest) { 
		logger_ = Configuration.logger_;
		fileHandler_ = new FileHandler("MOEAD.log");
		logger_.addHandler(fileHandler_);
		}
		problem = new SAS("SASSolutionType", factory, vars, numberOfObjectives_, numberOfConstraints_);
		
		algorithm = new IBEA_SAS(problem, factory);

		// Algorithm parameters
		int popsize = EAConfigure.getInstance().pop_size;
		int generations = EAConfigure.getInstance().generation;
		algorithm.setInputParameter("archiveSize", popsize*5);
		algorithm.setInputParameter("populationSize", popsize);
		algorithm.setInputParameter("maxEvaluations", popsize * generations);
		
		// Crossover operator
		parameters = new HashMap();
		parameters.put("probability", EAConfigure.getInstance().crossover_rate);
		parameters.put("distributionIndex", 20.0);
		// This needs to change in testing.
		parameters.put("jmetal.metaheuristics.moead.SASSolutionInstantiator",
				factory);
		crossover = CrossoverFactory.getCrossoverOperator(
				"UniformCrossoverSAS", parameters);
		//}
		
		
		if(SASAlgorithmAdaptor.isPreserveInvalidSolution) {
			algorithm.setInputParameter("vandInvCoEvolver", new SASValidityAndInvalidityCoEvolver(factory, 0.9, 0.1, 20));
		}

		// Mutation operator
		parameters = new HashMap();
		parameters.put("probability", EAConfigure.getInstance().mutation_rate);
		parameters.put("distributionIndex", 20.0);
		mutation = MutationFactory.getMutationOperator("BitFlipMutation",
				parameters);

		if(SASAlgorithmAdaptor.isSeedSolution) {
			//algorithm.setInputParameter("seeder", new Seeder(mutation));	
			//algorithm.setInputParameter("seeder", NewSeeder.getInstance(mutation));
			algorithm.setInputParameter("seeder", FixedSeeder.getInstance());	
		}
		
		selection = SelectionFactory.getSelectionOperator("BinaryTournament2", parameters);

		algorithm.addOperator("crossover", crossover);
		algorithm.addOperator("mutation", mutation);
		algorithm.addOperator("selection", selection);
		
		
		long initTime = System.currentTimeMillis();
		
		SolutionSet population = algorithm.execute();
		long estimatedTime = System.currentTimeMillis() - initTime;
		if(SAS.isTest) { 
		logger_.setLevel(Level.CONFIG);
		logger_.log(Level.CONFIG, "Total execution time: " + estimatedTime + "ms");
		
		String str = "data/IBEA/SAS";
		
		Utils.deleteFolder(new File(str+ "/results.dat"));
		Utils.createFolder(str);
		
		population.printObjectivesToFile(str + "/results.dat");
		}
		return population;

	} 

	protected ApproachType getName() {
		return ApproachType.IBEA;
	}

	protected Solution findSoleSolutionAfterEvolution(SolutionSet pareto_front) {
		// find the knee point
		Solution individual = pareto_front.get(PseudoRandom.randInt(0, pareto_front.size() - 1)); 
			
		
		for (int i = 0; i < problem.getNumberOfObjectives(); i++)
			System.out.print(individual.getObjective(i) + "\n");
		
		
		String str = "data/IBEA/SAS";
		if(SAS.isTest) 
		Utils.deleteFolder(new File(str+ "/knee_results.dat"));
		SolutionSet set = new SolutionSet(1);
		set.add(individual);
		if(SAS.isTest) 
		set.printObjectivesToFile(str + "/knee_results.dat");
		
		return individual;
	}

	protected SolutionSet doRanking(SolutionSet population) {
		// TODO Auto-generated method stub
		return ((IBEA_SAS)algorithm).doRanking(population);
	}
	
} // IBEA_main.java
