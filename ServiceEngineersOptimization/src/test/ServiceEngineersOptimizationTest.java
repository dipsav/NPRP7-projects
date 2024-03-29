package test;

import static org.junit.Assert.*;

import ilog.concert.IloException;

import org.junit.Test;

import OptimizationProblem.ServiceEngineersOptimization;

public class ServiceEngineersOptimizationTest extends ServiceEngineersOptimization{

	public ServiceEngineersOptimizationTest() {
		super(0, null, null, null, null, null, null, false);
	}

	//@Test
	public void testGetIndexFunctions() {
		int[] truncation_levels={3,3,3,4};
		double[] mu={0,0,0,0};
		
		ServiceEngineersOptimization obj = new ServiceEngineersOptimization(0, null, mu, null, null, truncation_levels, truncation_levels, true);
				
		int[] n = {1,1,1,4};
		int k = obj.getIndex(n);
		assertEquals(k, 109);
				
		int[] n1 = obj.getIndices(109);
		assertArrayEquals(n1, n);
	
	}

	//@Test
	public void testFormulation() throws IloException {
		int[] truncation_levels={5,5};
		double[] lambda = {10.0, 10.0};
		double[] mu={1.0,3.0};
		//double[] alpha = {0.0, 1.0};
		double[] lostCost = {300, 300};
		double[] engineerPartCost = {1.0, 1.0};
		
		ServiceEngineersOptimization obj = new ServiceEngineersOptimization(1, lambda, mu, lostCost, engineerPartCost, truncation_levels, truncation_levels, true);
				
		obj.formLPbyRows();
		//obj.formLP();
		
		//int[] n={1,1};
		//obj.addIndicatorLimits(n);
		obj.exportModel();
		
		
		System.out.println(obj.optimize());
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
		int[] truncation_levels={30,8,8,8,8};
		//int[] truncation_levels={10,6,6,6};
		double[] lambda = {10.0, 2.0, 3.0, 2.0, 3.0};
		double[] mu={1.0, 3.0, 2.0, 3.0, 2.0};
		//double[] alpha = {0.0, 0.2, 0.3, 0.2, 0.3};
		double[] lostCost = {300, 300};
		double[] engineerPartCost = {1.0, 1.0, 1.0, 1.0, 1.0};
				
		ServiceEngineersOptimization obj = new ServiceEngineersOptimization(4, lambda, mu, lostCost, engineerPartCost, truncation_levels, truncation_levels, true);
		obj.StartTimer();
		obj.ElapsedTime("Starting");
		
		obj.formLPbyRows();
		obj.ElapsedTime("MIP is formed");
		
		obj.setStartSolution();
		//obj.tuneModel();
		obj.setParameters();
		
		int[] ll={8,3,3,3,3};
		int[] ul={10,5,5,5,5};
		obj.setIndicatorLimits(ll,ul);
		//obj.exportModel();
		//obj.setLoggingOff();
		
		obj.ElapsedTime("Start Solving");
		boolean ready = true;
		while(!ready){
			//obj.exportModel();
			System.out.println("Curent objective: " + obj.optimize());
			obj.printIndicatorShort();
			ready = obj.doIteration();
			obj.ElapsedTime("");
		}

		
		System.out.println(obj.optimize());
		obj.printIndicatorSums();
		//obj.printMarginalProbabilities();
		//obj.printPvariables2D();
		
		//obj.printPvariables();
		//obj.printYvariables();
		
		obj.StopTimer();
		
		obj.cleanupModel();
	}

	
}
