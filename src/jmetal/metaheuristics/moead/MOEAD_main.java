package jmetal.metaheuristics.moead;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Random;
import java.util.logging.FileHandler;
import java.util.logging.Level;

import org.femosaa.core.EAConfigure;
import org.femosaa.core.SAS;
import org.femosaa.core.SASAlgorithmAdaptor;
import org.femosaa.core.SASSolutionInstantiator;
import org.femosaa.invalid.SASValidityAndInvalidityCoEvolver;
import org.femosaa.seed.FixedSeeder;
import org.femosaa.seed.NewSeeder;

import jmetal.core.Solution;
import jmetal.core.SolutionSet;
import jmetal.operators.crossover.CrossoverFactory;
import jmetal.operators.mutation.MutationFactory;
import jmetal.util.Configuration;
import jmetal.util.JMException;
import jmetal.util.PseudoRandom;

public class MOEAD_main extends MOEAD_STM_SAS_main{
	
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
		algorithm = new MOEAD(problem, factory);

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
			//algorithm.setInputParameter("seeder", NewSeeder.getInstance(mutation));
			algorithm.setInputParameter("seeder", FixedSeeder.getInstance());	
		}

		algorithm.addOperator("crossover", crossover);
		algorithm.addOperator("mutation", mutation);

		
		
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
	protected Solution findSoleSolutionAfterEvolution(SolutionSet pareto_front) {
		// find the knee point
		Solution kneeIndividual = pareto_front.get(PseudoRandom.randInt(0, pareto_front.size() - 1)); 
		
		//for (int i = 0; i < problem.getNumberOfObjectives(); i++)
		//	System.out.print(kneeIndividual.getObjective(i) + "\n");
		
		
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
}
