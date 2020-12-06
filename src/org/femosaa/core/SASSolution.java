package org.femosaa.core;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import jmetal.core.Problem;
import jmetal.core.Solution;
import jmetal.core.Variable;
import jmetal.util.JMException;
import jmetal.util.PseudoRandom;

/**
 * Need to ensure that the order of variables is the same as control primitives.
 * 
 * The variable here refer to the index value, whereas the one in xValue in QualityOfService is the actual value
 * 
 * @author tao
 *
 */
public abstract class SASSolution extends Solution {
	public boolean isFromInValid = false;
	
	// Key = variable index of dependent variable, the VarEntity has the same order as the original values array.
	protected final static Map<Integer, SASVarEntity[]> dependencyMap = new HashMap<Integer, SASVarEntity[]>();
	protected final static Map<Integer, SASVarEntity[]> validationDependencyMap = new HashMap<Integer, SASVarEntity[]>();
	// Key = variable index, Value = list of main/dependent variable index.
	protected final static Map<Integer, List<Integer>> crossoverMap = new HashMap<Integer, List<Integer>>();
	
	// Key = variable index, Value = list of dependent variable index.
	protected final static Map<Integer, List<Integer>> mutationMap = new HashMap<Integer, List<Integer>>();

	protected static double[][] optionalVariables;
	
	//public abstract double getVariableValueFromIndexValue(int indexValue);

	public SASSolution(Problem problem) throws ClassNotFoundException {
		super(problem);
		// TODO Auto-generated constructor stub
	}
	
	public SASSolution(Problem problem, Variable[] variables) throws ClassNotFoundException {
		super(problem, variables);
		// TODO Auto-generated constructor stub
	}

	public SASSolution(Solution solution) {
		super(solution);
		// TODO Auto-generated constructor stub
	}
	
	public SASSolution(int numberOfObjectives) {
		super(numberOfObjectives);
	}
	
	public SASSolution() {
		super();
		// TODO Auto-generated constructor stub
	}

	public abstract double[] getObjectiveValuesFromIndexValue();

	public abstract double getVariableValueFromIndex(int index);
	
	/**
	 * Used for fuzzy requirement.
	 * @param f
	 */
	public abstract void updateNormalizationBounds(double[] f);
	
	
	public static synchronized void init(double[][] optionalVariables) {
		SASSolution.optionalVariables = optionalVariables;
	}
	
	public static synchronized void clearAndStoreForValidationOnly(){
		if(validationDependencyMap.size() != 0) return;
		validationDependencyMap.putAll(dependencyMap);
		dependencyMap.clear();
	}
	
	/**
	 * For testing only
	 */
	public static synchronized void putDependencyChainBack(){
		if(dependencyMap.size() != 0) return;
		dependencyMap.putAll(validationDependencyMap);
		validationDependencyMap.clear();
	}
	
	
	
	public static Map<Integer, SASVarEntity[]> getDependencyMap(){
		return dependencyMap;
	}

	public static void main(String[] a) {
		System.out.print(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory() + "\n");
		//VarEntity[] vars = new VarEntity[1000000]; 
		//VarEntity v1 = new VarEntity(0, null, null);
		for (int i = 0; i < 100000; i ++) {
			SASVarEntity v1 = new SASVarEntity(0, null, null);
			//Object obj = new Object();
		//vars[i] = v1;
		}
		System.out.print(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory() + "\n");
	}
	
	
	private int getUpperBoundforVariable(int index) throws JMException {	
		if (dependencyMap.containsKey(index)) {
			SASVarEntity v = dependencyMap.get(index)[(int) super.getDecisionVariables()[dependencyMap.get(index)[0].getVarIndex()].getValue()];
			return v.getOptionalValues(super.getDecisionVariables()).length - 1;
		} else {
			return optionalVariables[index].length - 1;
		}
	}
	
	private Set<Integer> getOptionalValueSet(int index) throws JMException {
		if (dependencyMap.containsKey(index)) {
			SASVarEntity v = dependencyMap.get(index)[(int) super.getDecisionVariables()[dependencyMap.get(index)[0].getVarIndex()].getValue()];
			return v.getOptionalValuesSet(super.getDecisionVariables());
		} else {
			return null;
		}
	
	}

	
	private int getLowerBoundforVariable(int index) throws JMException {
			return 0;		
	}

	private int translateIntoIndexInMainVariable(int index, int subIndex) throws JMException {
		if (dependencyMap.containsKey(index)) {
			SASVarEntity v = dependencyMap.get(index)[(int) super.getDecisionVariables()[dependencyMap.get(index)[0].getVarIndex()].getValue()];
			return v.getOptionalValues(super.getDecisionVariables())[subIndex];
		} else {
			return subIndex;
		}
	}
	
	
	private int[] getMainVariablesByDependentVariable(int index) {
		if (dependencyMap.containsKey(index)) {
			Integer[] ints = dependencyMap.get(index)[0]
					.getMainVariablesByDependentVariable(new ArrayList<Integer>());
			int[] result = new int[ints.length];

			for (int i = 0; i < result.length; i++) {
				result[i] = ints[i];
			}

			return result;
		}
		return null;
	}

	private List<Integer> getMainVariables(int index) {
		
		if (crossoverMap.containsKey(index)) {
			return crossoverMap.get(index);
		}
		
		List<Integer> list = new ArrayList<Integer>();
		for (Map.Entry<Integer, SASVarEntity[]> entity : dependencyMap.entrySet()) {
			if (index == entity.getKey()) {
				entity.getValue()[0].getMainVariablesByDependentVariable(list);
			} 
		}
		

		crossoverMap.put(index, list);
		
		return list;
	}
	
	private List<Integer> getDependentVariables(int index) {
		
		if (mutationMap.containsKey(index)) {
			return mutationMap.get(index);
		}
		
		List<Integer> list = new ArrayList<Integer>();
		for (Map.Entry<Integer, SASVarEntity[]> entity : dependencyMap.entrySet()) {
			
				
				SASVarEntity v = entity.getValue()[0];
				
				do {
					if(index == v.varIndex) {
						if(!list.contains(entity.getKey())) {
							list.add(entity.getKey());
						}
						
						break;
					}
					v = v.next == null? null : v.next[0];
				} while(v != null);
				
			
		}
		

		mutationMap.put(index, list);
		
		return list;
	}
	
	public void mutateWithDependency() throws JMException{
		for (int i = 0; i < super.getDecisionVariables().length; i++) {		
			this.mutateWithDependency(i, false);
		}
	}
	
	public void correctDependency() throws JMException {
		// Only trigger when no dependency injection.
		if (dependencyMap.size() == 0) {
			dependencyMap.putAll(validationDependencyMap);
			for (int i = 0; i < super.getDecisionVariables().length; i++) {
				this.mutateWithDependency(i, false);
			}
			dependencyMap.clear();
		}
	}
	
	public void mutateWithDependency(int i, boolean isMutate /*This is can be only true for the initial entrance*/) throws JMException{
	
			
		int upper = this.getUpperBoundforVariable(i);
		int lower = this.getLowerBoundforVariable(i);
		
		//isMutate = !isValid(this, i);
		if (isMutate || !isValid(this, i)) {
		
		
			int v = (int) (PseudoRandom.randInt(
					// In the implementation of SASSolution, we can ensure the right boundary is 
					// always used even under variable dependency.
					lower,
					upper));
		
			v = this.translateIntoIndexInMainVariable(i, v);
			super.getDecisionVariables()[i].setValue(v);
			List<Integer> list = this.getDependentVariables(i);
			for (Integer j : list) {
				this.mutateWithDependency(j, false);
			}
		}
		
	
	}
	
	public int[] getBounds(int i)  throws JMException{
		return new int[] {this.getLowerBoundforVariable(i),this.getUpperBoundforVariable(i)};
	}
	
	public void crossoverWithDependency(Solution parent1, Solution parent2, Solution offSpring1, Solution offSpring2) throws JMException{
		for (int i = 0; i < parent1.numberOfVariables(); i++) {
			// Crossover has been completed
			this.crossoverWithDependency(i, parent1, parent2, offSpring1, offSpring2, false);
		}
	}
	
	public void crossoverWithDependency(int i, Solution parent1,
			Solution parent2, Solution offSpring1, Solution offSpring2, boolean isCrossover /*This is can be only true for the initial entrance*/)
			throws JMException {

//		// If has been swapped.
//		if (offSpring1.getDecisionVariables()[i].getValue() == parent2
//				.getDecisionVariables()[i].getValue()
//				&& parent1.getDecisionVariables()[i].getValue() != parent2
//						.getDecisionVariables()[i].getValue()) {
//			return;
//		}

//		
		if (isCrossover && offSpring1.getDecisionVariables()[i].getValue() == parent1
				.getDecisionVariables()[i].getValue()
				&& parent1.getDecisionVariables()[i].getValue() != parent2
						.getDecisionVariables()[i].getValue()) {
			
			int valueX1 = (int) parent1.getDecisionVariables()[i]
			                   							.getValue();
			int valueX2 = (int) parent2.getDecisionVariables()[i]
			                   							.getValue();
			offSpring1.getDecisionVariables()[i].setValue(valueX2);
			offSpring2.getDecisionVariables()[i].setValue(valueX1);
			
		}
		

		 
		
		if (i >= parent1.getDecisionVariables().length) {
			return;
		}

		// If it swap and they are originally unequal. This might not needed.
		if (offSpring1.getDecisionVariables()[i].getValue() == parent2
				.getDecisionVariables()[i].getValue()
				&& parent1.getDecisionVariables()[i].getValue() != parent2
						.getDecisionVariables()[i].getValue()) {
			// This should include both dependent and main variable of this variable.
			List<Integer> mainList = ((SASSolution) parent1)
					.getMainVariables(i);
			
			List<Integer> dependentList = ((SASSolution) parent1)
			.getDependentVariables(i);
			
			// Do mains
			if(!isValid((SASSolution)offSpring1, i) || !isValid((SASSolution)offSpring2, i) ) {
				for (Integer j : mainList) {
					// swap all main variable, if they have not been swapped. This should not occur
					if (offSpring1.getDecisionVariables()[j].getValue() == parent1
							.getDecisionVariables()[j].getValue()
							&& parent1.getDecisionVariables()[j].getValue() != parent2
									.getDecisionVariables()[j].getValue()) {
						
						
						
							int valueX1 = (int) parent1.getDecisionVariables()[j]
									.getValue();
							int valueX2 = (int) parent2.getDecisionVariables()[j]
									.getValue();
							offSpring1.getDecisionVariables()[j].setValue(valueX2);
							offSpring2.getDecisionVariables()[j].setValue(valueX1);
							// Ensure that the main/dependent variable of the newly swapped variable are also swapped.
							this.crossoverWithDependency(j, parent1, parent2,
									offSpring1, offSpring2, false);
							
						} 
//					else {
//							throw new RuntimeException("index " + j + " cause in valid but it has already been swapped!");
//						}
					} 
			}
			
			
			
			// Do dependents
			for (Integer j : dependentList) {
				// swap if it the prior swap causes any variables in the dependency becomes invalid.
				if(!isValid((SASSolution)offSpring1, j) || !isValid((SASSolution)offSpring2, j) ) {
				// swap all dependent variable, if they have not been swapped. 
				if (offSpring1.getDecisionVariables()[j].getValue() == parent1
						.getDecisionVariables()[j].getValue()
						&& parent1.getDecisionVariables()[j].getValue() != parent2
								.getDecisionVariables()[j].getValue()) {
					
					
					
						int valueX1 = (int) parent1.getDecisionVariables()[j]
								.getValue();
						int valueX2 = (int) parent2.getDecisionVariables()[j]
								.getValue();
						offSpring1.getDecisionVariables()[j].setValue(valueX2);
						offSpring2.getDecisionVariables()[j].setValue(valueX1);
						// Ensure that the main/dependent variable of the newly swapped variable are also swapped.
						this.crossoverWithDependency(j, parent1, parent2,
								offSpring1, offSpring2, false);
						
					} 
//				else {
//						throw new RuntimeException("index " + j + " cause in valid but it has already been swapped!");
//					}
				} 
			}
		}
//		
//		if((offSpring1.getDecisionVariables()[3].getValue() == offSpring1.getDecisionVariables()[9].getValue() && offSpring1.getDecisionVariables()[3].getValue() == 1) || 
//				(offSpring2.getDecisionVariables()[3].getValue() == offSpring2.getDecisionVariables()[9].getValue() && offSpring2.getDecisionVariables()[3].getValue() == 1)){
//			
//			@SuppressWarnings("unused")
//			boolean is1  = isValid((SASSolution)offSpring1, 9);
//			boolean is2  = isValid((SASSolution)offSpring2, 9);
//			System.currentTimeMillis();
//		
//		}

	}
	
	// This should only be used to check, regardless if dependency
	// has been injected or not.
	public boolean isSolutionValid(){
		for (int i = 0; i < super.getDecisionVariables().length; i ++) {
			try {
				if(!checkIsValidOnly(i)) {
					return false;
				}
			} catch (JMException e) {
				e.printStackTrace();
			}
		}
		
		return true;
	}
	
	
	public int countDegreeOfViolation(){
		int count = 0;
		for (int i = 0; i < super.getDecisionVariables().length; i ++) {
			try {
				if(!checkIsValidOnly(i)) {
					count++;
				}
			} catch (JMException e) {
				e.printStackTrace();
			}
		}
		
		return count;
	}
	
	public double getProbabilityToBeNaturallyRepaired(){
		double prob = 1;
		for (int i = 0; i < super.getDecisionVariables().length; i ++) {
			try {
				// Only check for invalid genes.
				prob = prob * getProbabilityToBeNaturallyRepaired(i);
				
			} catch (JMException e) {
				e.printStackTrace();
			}
		}
		
		return prob;
	}
	
	private double getProbabilityToBeNaturallyRepaired(int i) throws JMException {
		
		Map<Integer, SASVarEntity[]> map = validationDependencyMap.size() == 0? dependencyMap : validationDependencyMap;
		
		
		Set<Integer> set = null;
		if (map.containsKey(i)) {
			SASVarEntity v = map.get(i)[(int) super.getDecisionVariables()[map.get(i)[0].getVarIndex()].getValue()];
			set = v.getOptionalValuesSet(super.getDecisionVariables());
		}
		double p = 1.0;
		
		final int value = (int) super.getDecisionVariables()[i].getValue();	
		// Means no dependency or this gene is value, then do nothing.
		if(set == null || set.contains(value)) {
			return p;
		}
		

		// The number of valid option / the total number of options.
		p = (double)set.size() / (double)optionalVariables[i].length;
		//System.out.print(set.size() + "\n");
		// prob of mutation not crossover + prob of crossover not mutation + prob of mutation and crossover
		return EAConfigure.getInstance().mutation_rate * p * (1 - EAConfigure.getInstance().crossover_rate) +
		EAConfigure.getInstance().crossover_rate * p * (1 - EAConfigure.getInstance().mutation_rate) +
		EAConfigure.getInstance().crossover_rate * p * EAConfigure.getInstance().mutation_rate * p;
	}
	
	
	private boolean checkIsValidOnly(int i) throws JMException {

		Map<Integer, SASVarEntity[]> map = validationDependencyMap.size() == 0? dependencyMap : validationDependencyMap;
		
		
		int value = (int) super.getDecisionVariables()[i].getValue();
		Set<Integer> set = null;
		if (map.containsKey(i)) {
			SASVarEntity v = map.get(i)[(int) super.getDecisionVariables()[map.get(i)[0].getVarIndex()].getValue()];
			set = v.getOptionalValuesSet(super.getDecisionVariables());
		}

		// Means no dependency
		if(set == null) {
			return true;
		}
		

		//System.out.print(set.contains(value)+"\n");
		return set.contains(value);
	}
	
	private boolean isValid(SASSolution s, int i) throws JMException{
		
		int value = (int)s.getDecisionVariables()[i].getValue();	
		Set<Integer> set = s.getOptionalValueSet(i);
		
		// Means no dependency
		if(set == null) {
			return true;
		}
		
//		int upper = s.getUpperBoundforVariable(i);
//		int lower = s.getLowerBoundforVariable(i);
//		int traUpper = s.translateIntoIndexInMainVariable(i, upper);
//		int traLower = s.translateIntoIndexInMainVariable(i, lower);
//		
//		if  (value > traUpper || value < traLower) {
//			return false;
//		}
		
		return set.contains(value);
	}
	

}
