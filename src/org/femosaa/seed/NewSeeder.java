package org.femosaa.seed;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import org.femosaa.core.SASSolutionInstantiator;

import jmetal.core.Operator;
import jmetal.core.Problem;
import jmetal.core.Solution;
import jmetal.core.SolutionSet;
import jmetal.util.JMException;
import jmetal.util.Ranking;


/**
 * Contain seeding strategies for WSDREAM
 * @author tao
 *
 */
public class NewSeeder extends Seeder{
	
	private SeedingStrategy strategy;
	public int no_of_seed = 10;
	// They should be all instance of SASSolution
	// size of history needs >= no_of_seed
	private SolutionSet history = new SolutionSet();

	private static NewSeeder seeder = null;
	private static long time = 0;
	private static int count = 0;
	
	
	public static NewSeeder getInstance(Operator mutationOperator){
		if(seeder == null) {
			seeder = new NewSeeder(mutationOperator);
		}
		return seeder;
	}
	
	public static NewSeeder getInstance(){
		if(seeder == null) {
			seeder = new NewSeeder(null);
		}
		return seeder;
	}
	
	private NewSeeder(Operator mutationOperator) {
		super(mutationOperator);
		// TODO Auto-generated constructor stub
	}
	
	public SeedingStrategy getSeedingStrategy(){
		return strategy;
	}
	
	public void setSeedingStrategy(SeedingStrategy strategy){
		this.strategy = strategy;
	}
	
	public void addHistory(SolutionSet set){
		history = history.union(set);
	}
	
	public void addHistory(Solution s){
		history.add(s);
	}
	
	public int getSeedsCount(){
		return history.size();
	}
	
	public void reset(){
		history.clear();
	}
	
	public void seeding(SolutionSet population,
			SASSolutionInstantiator factory, Problem problem_,
			int populationSize) throws JMException {
		
		// The solutions in history are often non-evaluated so need to get their fitness values.
		for (int i = 0; i < history.size(); i++) {
			problem_.evaluate(history.get(i));
			problem_.evaluateConstraints(history.get(i));
		}
		
		
		
		long t = System.currentTimeMillis();
		SolutionSet set = null;
		switch (strategy) {
		case AO_Seed:
			set = AO_Seed(no_of_seed);
			break;
		case SO_Seed:
			set = SO_Seed(no_of_seed);
			break;
		case H_Seed:
			set = H_Seed(no_of_seed);
			break;
		case R_Seed:
			set = R_Seed(no_of_seed);
			break;

		}
		time += System.currentTimeMillis() - t;
		count ++;
		System.out.print("Average seeding overhead " + (time/count) + " ms\n");
		
		for(int i = 0; i < set.size(); i++) {
			population.add(set.get(i));
		}
		System.out.print("Number of seeds: " +  set.size() + "\n");
		
        int remain = populationSize - set.size();
		
		if(remain < 0) {
			throw new JMException("no_of_seed should not be larger than the population size.");
		}
		
		Solution newSolution;
		for (int i = 0; i < remain; i++) {
			newSolution = factory.getSolution(problem_);
			problem_.evaluate(newSolution);
			problem_.evaluateConstraints(newSolution);
			population.add(newSolution);
		
		} 
		
	}
	
	private SolutionSet AO_Seed(int no){
		Solution s = history.get(0);
		SolutionSet front = new SolutionSet();
		for(int i = 0; i < no; i++) {
			front.add(s);
			printSeeds(s);
		}
		
	
		return front;
	}
	
	private SolutionSet SO_Seed(int no){
		
		SolutionSet front = new SolutionSet();
		// no needs to be about to divided by 3.
		for(int i = 0; i < no/3; i++) {
			front.add(history.get(0));
			front.add(history.get(1));
			front.add(history.get(2));
			printSeeds(history.get(0));
			printSeeds(history.get(1));
			printSeeds(history.get(2));
		}
		
		return front;
	}
	
	private SolutionSet H_Seed(int no){
		Ranking ranking = new Ranking(history);
		SolutionSet front = new SolutionSet();
		//front = front.union(ranking.getSubfront(0));
		
		Iterator itr = ranking.getSubfront(0).iterator();
		while(itr.hasNext()) {
			front.add((Solution)itr.next());
		}
		
		System.out.print("Nondominated front size " + front.size() + "\n");
		if (front.size() >= no) {
			for(int i = front.size() - 1; i > no - 1; i--) {
				// This is basically random
				front.remove(i);
			}
			//System.out.print("k " + front.size() + "*****\n" );
		} else {
			int k = 1;
			
			while (front.size() < no) {
				int count = front.size();
				if(ranking.getSubfront(k).size() > no - count) {
					for(int i = 0; i < no - count; i++) {
						front.add(ranking.getSubfront(k).get(i));
					}
				} else {
					
					itr = ranking.getSubfront(k).iterator();
					while(itr.hasNext()) {
						front.add((Solution)itr.next());
					}
				}
				
			
				k++;
			}
		}
		printSeeds(front);
		return front;
	}
	
	private SolutionSet R_Seed(int no){
		List<Solution> list = new ArrayList<Solution> ();
		for (int i = 0; i < history.size(); i++) {
			list.add(history.get(i));
		}
		Collections.shuffle(list);
		SolutionSet front = new SolutionSet(no);
		for (int i = 0; i < no; i++) {
			front.add(list.get(i));
		}
		printSeeds(front);
		return front;
	}
	

	public enum SeedingStrategy {
		AO_Seed, SO_Seed, H_Seed, R_Seed, NONE
	}
	
	private void printSeeds(SolutionSet front){
		for (int i = 0; i < front.size();i++) {
			printSeeds(front.get(i));
		}
	}

	
	private void printSeeds(Solution s){
		String s1 = "";
		for (int j = 0; j < s.numberOfObjectives(); j++) {
			s1 += s.getObjective(j) + " ";
		}
		System.out.print(s+"\n");
	}

}
