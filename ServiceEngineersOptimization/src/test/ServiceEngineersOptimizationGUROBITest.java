package test;

import static org.junit.Assert.*;
import gurobi.GRBException;

import org.junit.Test;

import OptimizationProblem.ServiceEngineersOptimizationGUROBI;

public class ServiceEngineersOptimizationGUROBITest extends ServiceEngineersOptimizationGUROBI{

	public ServiceEngineersOptimizationGUROBITest() {
		super(0, null, null, null, null, null, null, false);
	}

	//@Test
	public void testGetIndexFunctions() {
		int[] truncation_levels_lw={0,0,0,0};
		int[] truncation_levels_up={3,3,3,4};
		double[] mu={0,0,0,0};
		
		ServiceEngineersOptimizationGUROBI obj = new ServiceEngineersOptimizationGUROBI(0, mu, null, null, null, truncation_levels_lw, truncation_levels_up, true);
				
		int[] n = {1,1,1,4};
		int k = obj.getIndex(n);
		assertEquals(k, 109);
				
		int[] n1 = obj.getIndices(109);
		assertArrayEquals(n1, n);
	
	}

	//@Test
	public void testFormulation() throws GRBException {
		int[] truncation_levels_lw={0,0};
		int[] truncation_levels_up={50,50};
		double[] lambda = {10.0, 10.0};
		double[] mu={1.0,3.0};
		//double[] alpha = {0.0, 1.0};
		double[] lostCost = {0.0, 300};
		double[] engineerPartCost = {1.0, 1.0};
		
		ServiceEngineersOptimizationGUROBI obj = new ServiceEngineersOptimizationGUROBI(1, lambda, mu, lostCost, engineerPartCost, truncation_levels_lw, truncation_levels_up, true);
				
		obj.formLP(null);
		
		//int[] n={1,1};
		//obj.addIndicatorLimits(n);
		obj.exportModel();
		
		
		obj.optimize();
		//obj.printIndicators();
		obj.printIndicatorSums();
		//obj.printMarginalProbabilities();
		//obj.printPvariables2D();
		
		//obj.printPvariables();
		//obj.printYvariables();
		
		obj.cleanupModel();
	}
	
	@Test
	public void testMultiDimentional() throws GRBException {
		int[] truncation_levels_lw={0,0,0,0,0};
		int[] truncation_levels_up={30,8,8,8,8};
		double[] lambda = {10.0, 0, 2.0, 3.0, 2.0, 3.0};
		double[] mu={1.0, 3.0, 2.0, 3.0, 2.0};
		//double[] alpha = {0.0, 0.2, 0.3, 0.2, 0.3};
		double[] lostCost = {0.0, 300, 300, 300, 300};
		double[] engineerPartCost = {1.0, 1.0, 1.0, 1.0, 1.0};
		
		ServiceEngineersOptimizationGUROBI obj = new ServiceEngineersOptimizationGUROBI(4, lambda, mu, lostCost, engineerPartCost, truncation_levels_lw, truncation_levels_up, true);
		obj.StartTimer();		
		
		obj.formLP(null);
		
		//obj.setStartSolution(1);
		//obj.tuneModel();
		//obj.setParameters();
		
		int[] ll={0,0,0,0,0};
		int[] ul={2,2,2,2,2};
		obj.setIndicatorLimits(ll,ul);
		//obj.exportModel();
		obj.setLoggingOff();
		
		boolean ready = false;
		while(!ready){
			//obj.exportModel();
			System.out.println("Curent objective: " + obj.optimize());
			obj.printIndicatorShort();
			//obj.printIndicatorSums();
			ready = obj.doIteration();
			obj.ElapsedTime("");
		}

		//obj.addIndicatorLimits(ll,null);
		//obj.exportModel();
		
		//obj.optimize();
		//obj.printIndicators();
		obj.printIndicatorSums();
		obj.StopTimer();
		//obj.printMarginalProbabilities();
		//obj.printPvariables2D();
		
		//obj.printPvariables();
		//obj.printYvariables();
		
		obj.cleanupModel();
	}

	
}
