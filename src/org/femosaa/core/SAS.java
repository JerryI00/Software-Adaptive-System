package org.femosaa.core;

import jmetal.core.*;
import jmetal.encodings.solutionType.IntSolutionType;
import jmetal.util.JMException;
import jmetal.util.wrapper.XInt;

public class SAS extends Problem {

	public static final boolean isTest = false;
	public static final long TIME_THRESHOLD = 40000;
	
	public SASSolutionInstantiator factory;
	/**
	 * Constructor
	 */
//	public SAS(String solutionType) throws ClassNotFoundException {
//		this(solutionType, 10);
//	}

	/**
	 * vars [number of variables][upper and lower bounds]
	 * 
	 * Create a new instance of problem for SAS
	 */
	public SAS(String solutionType, SASSolutionInstantiator factory, int[][] vars, int numberOfObjectives_, int numberOfConstraints_) throws ClassNotFoundException {
		numberOfVariables_ = vars.length;
		this.factory = factory;
		this.numberOfObjectives_ = numberOfObjectives_;
		this.numberOfConstraints_ = numberOfConstraints_;
		problemName_ = "SAS";

		upperLimitSAS_ = new int[numberOfVariables_];
		lowerLimitSAS_ = new int[numberOfVariables_];
//
		// Establishes upper and lower limits for the variables
		for (int i = 0; i < numberOfVariables_; i++) {
			lowerLimitSAS_[i] = vars[i][0];
			upperLimitSAS_[i] = vars[i][1];
		}

		if (solutionType.compareTo("SASSolutionType") == 0)
			solutionType_ = new SASSolutionType(this);
		else {
			System.out.println("Error: solution type " + solutionType + " invalid");
			System.exit(-1);
		}
	}
	
	/**
	 * Evaluate a solution
	 */
	public void evaluate(Solution solution) throws JMException {
		//XInt x = new XInt(solution);
		
		//int[] var  = new int[numberOfVariables_];
//		double[] varValues  = new double[numberOfVariables_];
//		double[] f = new double[numberOfObjectives_];
		SASSolution sol = null;
		if (solution instanceof SASSolution) {
			sol = (SASSolution)solution;
		}
		
		// put "var[]" into Tao's System, note that var[] is only the index, we need to get the variable value
		// from a variable list
//		for (int i = 0; i < numberOfVariables_; i++) {
//			var[i] = x.getValue(i);
//			varValues[i] = sol.getVariableValueFromIndex(var[i]);
//		}
		
		//for (int i = 0; i < numberOfVariables_; i++) 
			//var[i] = x.getValue(i);
		
		// It is possible to directly retrieve objective' values from indices, and this
		// can be done from my end.
		double[] f = sol.getObjectiveValuesFromIndexValue();
		sol.updateNormalizationBounds(f);
		//System.out.print("Solution " + f[0] + " : " +  f[1] + "\n");
		for (int i = 0; i < f.length ; i ++) {
			solution.setObjective(i, f[i]);
		}
	}
}