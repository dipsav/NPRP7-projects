package test;

import static org.junit.Assert.*;
import gurobi.GRBException;
import ilog.concert.IloException;

import org.junit.Test;

import OptimizationProblem.ServiceEngineersOptimizationGUROBI;

public class ServiceEngineersOptimizationGUROBITest extends ServiceEngineersOptimizationGUROBI{

	public ServiceEngineersOptimizationGUROBITest() {
		super(0, null, null, 0, null, null, false);
	}

	//@Test
	public void testGetIndexFunctions() {
		int[] truncation_levels={3,3,3,4};
		double[] mu={0,0,0,0};
		
		ServiceEngineersOptimizationGUROBI obj = new ServiceEngineersOptimizationGUROBI(0, mu, null, 0, null, truncation_levels, true);
				
		int[] n = {1,1,1,4};
		int k = obj.getIndex(n);
		assertEquals(k, 109);
				
		int[] n1 = obj.getIndices(109);
		assertArrayEquals(n1, n);
	
	}

	//@Test
	public void testFormulation() throws GRBException {
		int[] truncation_levels={50,50};
		double lambda = 10.0;
		double[] mu={1.0,3.0};
		double[] alpha = {0.0, 1.0};
		double lostCost = 300;
		double[] engineerPartCost = {1.0, 1.0};
		
		ServiceEngineersOptimizationGUROBI obj = new ServiceEngineersOptimizationGUROBI(lambda, mu, alpha, lostCost, engineerPartCost, truncation_levels, true);
				
		obj.formLP(null);
		
		int[] n={1,1};
		//obj.addIndicatorLimits(n);
		obj.exportModel();
		
		
		obj.Optimize();
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
		int[] truncation_levels={30,5,5,5};
		double lambda = 10.0;
		double[] mu={1.0, 3.0, 2.0, 4.0};
		double[] alpha = {0.0, 0.2, 0.3, 0.2};
		double lostCost = 300;
		double[] engineerPartCost = {1.0, 1.0, 1.0, 1.0, 1.0};
		
		ServiceEngineersOptimizationGUROBI obj = new ServiceEngineersOptimizationGUROBI(lambda, mu, alpha, lostCost, engineerPartCost, truncation_levels, true);
				
		obj.formLP(null);
		
		//obj.setStartSolution(1);
		//obj.tuneModel();
		//obj.setParameters();
		
		int[] ll={18,4,5,5};
		int[] ul={18,4,5,5};
		obj.addIndicatorLimits(ll,null);
		obj.exportModel();
		
		obj.Optimize();
		//obj.printIndicators();
		obj.printIndicatorSums();
		//obj.printMarginalProbabilities();
		//obj.printPvariables2D();
		
		//obj.printPvariables();
		//obj.printYvariables();
		
		obj.cleanupModel();
	}

	
}
