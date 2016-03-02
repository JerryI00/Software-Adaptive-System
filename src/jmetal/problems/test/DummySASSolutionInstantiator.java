package jmetal.problems.test;

import jmetal.core.Problem;
import jmetal.core.Solution;
import jmetal.problems.SASSolutionInstantiator;

public class DummySASSolutionInstantiator implements SASSolutionInstantiator {

	@Override
	public Solution getSolution(Problem problem) {
		Solution sol = null;
		try {
			sol = new DummySASSolution(problem);
		} catch (ClassNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return sol;
	}

	@Override
	public Solution getSolution(Solution solution) {
		return new DummySASSolution(solution);
	}

}