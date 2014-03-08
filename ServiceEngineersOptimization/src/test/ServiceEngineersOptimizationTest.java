package test;

import static org.junit.Assert.*;
import ilog.concert.IloException;

import org.junit.Test;

import OptimizationProblem.ServiceEngineersOptimization;

public class ServiceEngineersOptimizationTest extends ServiceEngineersOptimization{

	public ServiceEngineersOptimizationTest() {
		super(0, 0, 0, null, null);
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

	@Test
	public void testFormulation() throws IloException {
		double[] mu={2.0,3.0};
		double[] alpha = {0.0, 1.0};
		
		ServiceEngineersOptimization obj = new ServiceEngineersOptimization(1, 10, 1, mu, alpha);
				
		obj.formLP();
	}
}
