package test;

import static org.junit.Assert.*;

import java.util.Date;

import ilog.concert.IloException;

import org.junit.Test;

import OptimizationProblem.ServiceEngineersOptimization;

public class ServiceEngineersOptimizationTest extends ServiceEngineersOptimization{

	public ServiceEngineersOptimizationTest() {
		super(0, null, null, 0, null, null, false);
	}

	//@Test
	public void testGetIndexFunctions() {
		int[] truncation_levels={3,3,3,4};
		double[] mu={0,0,0,0};
		
		ServiceEngineersOptimization obj = new ServiceEngineersOptimization(0, mu, null, 0, null, truncation_levels, true);
				
		int[] n = {1,1,1,4};
		int k = obj.getIndex(n);
		assertEquals(k, 109);
				
		int[] n1 = obj.getIndices(109);
		assertArrayEquals(n1, n);
	
	}

	//@Test
	public void testFormulation() throws IloException {
		int[] truncation_levels={50,10};
		double lambda = 10.0;
		double[] mu={1.0,3.0};
		double[] alpha = {0.0, 1.0};
		double lostCost = 300;
		double[] engineerPartCost = {1.0, 1.0};
		
		ServiceEngineersOptimization obj = new ServiceEngineersOptimization(lambda, mu, alpha, lostCost, engineerPartCost, truncation_levels, true);
				
		obj.formLP();
		
		int[] n={1,1};
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
		int[] truncation_levels={30,10,10,10,10};
		double lambda = 10.0;
		double[] mu={1.0, 3.0, 2.0, 3.0, 2.0};
		double[] alpha = {0.0, 0.2, 0.3, 0.2, 0.3};
		double lostCost = 300;
		double[] engineerPartCost = {1.0, 1.0, 1.0, 1.0, 1.0};
		
		
		ServiceEngineersOptimization obj = new ServiceEngineersOptimization(lambda, mu, alpha, lostCost, engineerPartCost, truncation_levels, true);
		obj.StartTimer();
		
		obj.formLP();
		
		obj.setStartSolution();
		//obj.tuneModel();
		obj.setParameters();
		
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
			ready = obj.doIteration();
			obj.ElapsedTime();
		}

		//obj.optimize();

		//System.out.println(obj.optimize());
		//obj.printIndicators();
		obj.printIndicatorSums();
		//obj.printMarginalProbabilities();
		//obj.printPvariables2D();
		
		//obj.printPvariables();
		//obj.printYvariables();
		
		obj.StopTimer();
		
		obj.cleanupModel();
	}

	
}
