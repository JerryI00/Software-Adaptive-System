package jmetal.problems.test;

import org.femosaa.core.SAS;
import org.femosaa.core.SASSolution;

import jmetal.core.Problem;
import jmetal.core.Solution;
import jmetal.core.SolutionSet;
import jmetal.metaheuristics.moead.MOEAD_STM_SAS_main;
import junit.framework.Assert;
import junit.framework.TestCase;

public class MainTest extends TestCase{

	public void testMain(){
		//Problem p;
		//DummySASSolution dss = null;
		try {
			//p = new SAS("SASSolutionType", DummySASSolution.vars, 2, 0);
			//dss = new DummySASSolution(p);
			
			DummySASSolutionInstantiator inst = new DummySASSolutionInstantiator();
			
			MOEAD_STM_SAS_main main = new MOEAD_STM_SAS_main();
			SASSolution.init(DummySASSolution.testOptionalVariables);
			Solution s = main.execute(inst, DummySASSolution.vars, 2, 0);
	
			
		//Assert.assertEquals(dss.getUpperBoundforVariable(1), 2);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	
	}
}
