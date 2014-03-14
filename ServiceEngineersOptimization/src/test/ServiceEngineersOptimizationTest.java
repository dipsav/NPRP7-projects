package test;

import static org.junit.Assert.*;
import ilog.concert.IloException;

import org.junit.Test;

import OptimizationProblem.ServiceEngineersOptimization;

public class ServiceEngineersOptimizationTest extends ServiceEngineersOptimization{

	public ServiceEngineersOptimizationTest() {
		super(0, null, null, 0, null, 0, false);
	}

//	@Test
//	public void testGetIndexFunctions() {
//		ServiceEngineersOptimization obj = new ServiceEngineersOptimization(3, 4, 0, null, null);
//				
//		int[] n = {1,1,1,1};
//		int k = obj.getIndex(n);
//		assertEquals(k, 85);
//				
//		int[] n1 = obj.getIndecies(85);
//		assertArrayEquals(n1, n);
//	
//	}

	//@Test
	public void testFormulation() throws IloException {
		int truncation_level=50;
		double lambda = 10.0;
		double[] mu={1.0,3.0};
		double[] alpha = {0.0, 1.0};
		double lostCost = 300;
		double[] engineerPartCost = {1.0, 1.0};
		
		ServiceEngineersOptimization obj = new ServiceEngineersOptimization(lambda, mu, alpha, lostCost, engineerPartCost, truncation_level, true);
				
		obj.formLP();
		
		int[] n={1,1};
		//obj.addIndicatorLimits(n);
		obj.exportModel();
		
		
		System.out.println(obj.Optimize());
		//obj.printIndicators();
		obj.printIndicatorSums();
		//obj.printMarginalProbabilities();
		//obj.printPvariables2D();
		
		//obj.printPvariables();
		//obj.printYvariables();
		
		obj.cleanupModel();
	}
	@Test
	public void testMultiDimentional() throws IloException {
		int truncation_level=8;
		double lambda = 10.0;
		double[] mu={1.0, 3.0, 2.0, 4.0};
		double[] alpha = {0.0, 0.2, 0.3, 0.5};
		double lostCost = 300;
		double[] engineerPartCost = {1.0, 1.0, 1.0, 1.0, 1.0};
		
		ServiceEngineersOptimization obj = new ServiceEngineersOptimization(lambda, mu, alpha, lostCost, engineerPartCost, truncation_level, true);
				
		obj.formLP();
		
		obj.setStartSolution();
		//obj.tuneModel();
		obj.setParameters();
		
		//int[] n={1,1};
		//obj.addIndicatorLimits(n);
		//obj.exportModel();
		
		System.out.println(obj.Optimize());
		//obj.printIndicators();
		obj.printIndicatorSums();
		//obj.printMarginalProbabilities();
		//obj.printPvariables2D();
		
		//obj.printPvariables();
		//obj.printYvariables();
		
		obj.cleanupModel();
	}
}
