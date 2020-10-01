/**
 * MOEAD_SAS_main.java
 * 
 * @author Ke Li <keli.genius@gmail.com>
 * 
 * Copyright (c) 2016 Ke Li
 * 
 * Note: This is a free software developed based on the open source project 
 * jMetal<http://jmetal.sourceforge.net>. The copy right of jMetal belongs to 
 * its original authors, Antonio J. Nebro and Juan J. Durillo. Nevertheless, 
 * this current version can be redistributed and/or modified under the terms of 
 * the GNU Lesser General Public License as published by the Free Software 
 * Foundation, either version 3 of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, 
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *  
 */

package jmetal.metaheuristics.moead;

import jmetal.core.*;
import jmetal.metaheuristics.nsgaII.NSGAII_SAS;
import jmetal.operators.crossover.*;
import jmetal.operators.mutation.*;
import jmetal.operators.selection.*;
import jmetal.problems.*;
import jmetal.problems.DTLZ.*;
import jmetal.problems.ZDT.*;
import jmetal.problems.cec2009Competition.*;
import jmetal.problems.WFG.*;

import jmetal.util.Configuration;
import jmetal.util.JMException;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.logging.FileHandler;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.femosaa.core.EAConfigure;
import org.femosaa.core.SAS;
import org.femosaa.core.SASAlgorithmAdaptor;
import org.femosaa.core.SASProblemFactory;
import org.femosaa.core.SASSolutionInstantiator;
import org.femosaa.invalid.SASValidityAndInvalidityCoEvolver;
import org.femosaa.seed.NewSeeder;


public class MOEAD_STM_SAS_main extends SASAlgorithmAdaptor{
	public static Logger logger_; // Logger object
	public static FileHandler fileHandler_; // FileHandler object

	public static int popsize = EAConfigure.getInstance().pop_size;//100;
	public static int factor = EAConfigure.getInstance().generation;//10;
	
	protected Problem problem; // The problem to solve
	protected Algorithm algorithm; // The algorithm to use
	Operator crossover; // Crossover operator
	Operator mutation; // Mutation operator

	
	/**
	 * @param args
	 *            Command line arguments. The first (optional) argument
	 *            specifies the problem to solve.
	 * @throws JMException
	 * @throws IOException
	 * @throws SecurityException
	 *             Usage: three options - jmetal.metaheuristics.moead.MOEAD_main
	 *             - jmetal.metaheuristics.moead.MOEAD_main problemName -
	 *             jmetal.metaheuristics.moead.MOEAD_main problemName
	 *             ParetoFrontFile
	 * @throws ClassNotFoundException
	 */
	public static void main(String[] args) throws JMException,
			SecurityException, IOException, ClassNotFoundException {
		Problem problem; // The problem to solve
		Algorithm algorithm; // The algorithm to use
		Operator crossover; // Crossover operator
		Operator mutation; // Mutation operator

		HashMap parameters; // Operator parameters

	
 		
		// Logger object and file to store log messages
		logger_ = Configuration.logger_;
		fileHandler_ = new FileHandler("MOEAD.log");
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

		algorithm = new MOEAD_STM_SAS(problem, null);

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
		parameters.put("probability", 1.0 / problem.getNumberOfObjectives());
		parameters.put("distributionIndex", 20.0);
		mutation = MutationFactory.getMutationOperator("BitFlipMutation", parameters);

		algorithm.addOperator("crossover", crossover);
		algorithm.addOperator("mutation", mutation);

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

		// Logger object and file to store log messages
		if(SAS.isTest) { 
		logger_ = Configuration.logger_;
		fileHandler_ = new FileHandler("MOEAD.log");
		logger_.addHandler(fileHandler_);

		}
		problem = new SAS("SASSolutionType", factory, vars, numberOfObjectives_, numberOfConstraints_);
		
		//algorithm = new MOEAD_STM_SAS_STATIC(problem, factory);
		/**Important!!!! remember to change it***/
		algorithm = new MOEAD_STM_SAS(problem, factory);

		// Algorithm parameters
	
		algorithm.setInputParameter("populationSize", EAConfigure.getInstance().pop_size);
		algorithm.setInputParameter("maxEvaluations", EAConfigure.getInstance().pop_size * EAConfigure.getInstance().generation);
		
		algorithm.setInputParameter("dataDirectory", System.getProperty("os.name").startsWith("Mac")? "weight" : "/home/tao/weight");
		//algorithm.setInputParameter("lambda", factory.getLambda());
		// Crossover operator
//		int tag = 2;
//		if (tag == 1) {
//			parameters = new HashMap();
//			parameters.put("CR", 0.5);
//			parameters.put("F", 0.5);
//			crossover = CrossoverFactory.getCrossoverOperator(
//					"DifferentialEvolutionCrossover", parameters);
//		} else {
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
			algorithm.setInputParameter("seeder", NewSeeder.getInstance(mutation));			
		}
		
		algorithm.addOperator("crossover", crossover);
		algorithm.addOperator("mutation", mutation);


		//algorithm.addOperator("selection", SelectionFactory.getSelectionOperator("BinaryTournament2", parameters));
		
		long initTime = System.currentTimeMillis();
		
		SolutionSet population = algorithm.execute();
		long estimatedTime = System.currentTimeMillis() - initTime;
		if(SAS.isTest) { 
		logger_.setLevel(Level.CONFIG);
		logger_.log(Level.CONFIG, "Total execution time: " + estimatedTime + "ms");
		
		String str = "data/" + problem.getName()
		+ "M" + problem.getNumberOfObjectives() + "/SAS";
		
		Utils.deleteFolder(new File(str+ "/results.dat"));
		Utils.createFolder(str);
		
		population.printObjectivesToFile(str + "/results.dat");
		}
		return population;

	}


	@Override
	protected ApproachType getName() {
		return ApproachType.MOEAD_STM_D_K;
	}

	@Override
	protected Solution findSoleSolutionAfterEvolution(SolutionSet pareto_front) {
		// find the knee point
		Solution kneeIndividual = null;
		if(algorithm instanceof MOEAD_STM_SAS) {
		     kneeIndividual = ((MOEAD_STM_SAS)algorithm).kneeSelection(pareto_front);
		} else 	if(algorithm instanceof MOEAD) {
			 kneeIndividual = ((MOEAD)algorithm).kneeSelection(pareto_front);
		}
		
		
		
		for (int i = 0; i < problem.getNumberOfObjectives(); i++)
			System.out.print(kneeIndividual.getObjective(i) + "\n");
		
		
		String str = "data/" +problem.getName()
		+ "M" + problem.getNumberOfObjectives() + "/SAS";
		if(SAS.isTest) 
		Utils.deleteFolder(new File(str+ "/knee_results.dat"));
		
		SolutionSet set = new SolutionSet(1);
		set.add(kneeIndividual);
		
		if(SAS.isTest) 
		set.printObjectivesToFile(str + "/knee_results.dat");
		
		return kneeIndividual;
	}

	
	@Override
	protected SolutionSet doRanking(SolutionSet population) {
		if(algorithm instanceof MOEAD_STM_SAS) {
			return ((MOEAD_STM_SAS)algorithm).doRanking(population);
		} else 	if(algorithm instanceof MOEAD) {
			return ((MOEAD)algorithm).doRanking(population);
		}
		
		return null;
	}

} // MOEAD_main
