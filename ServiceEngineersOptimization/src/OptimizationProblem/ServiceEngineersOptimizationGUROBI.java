package OptimizationProblem;

import gurobi.*;

import ilog.concert.IloColumn;
import ilog.concert.IloException;
import ilog.concert.IloNumVar;

import java.util.Iterator;
import java.util.StringTokenizer;

public class ServiceEngineersOptimizationGUROBI {
	
	int N; //number of spare types
	int[] M; //state space truncation limit
	double lambda; //failure rate;
	double[] mu;   //service lead time service rates, mu[0] corresponds to engenders, the rest to spares
	double[] alpha; //probability that spare i is requested
	double lostCost;
	double engnCost;
	double[] engineerPartCost;
	
	int M_max;
	boolean logging;
	
	
	@SuppressWarnings("restriction")
	GRBVar[][] I_var;
	@SuppressWarnings("restriction")
	GRBVar[] p_var, y_var;
	
	double[][] p_mar;
	
    @SuppressWarnings("restriction")
	GRBModel  model;
	//IloCplex model;
    GRBEnv    env;

    
    public static String workpath = "/Users/andrei/Documents/Research/Service Engineers";
    static String defaultParamFile = workpath + "/gurobiParameters.prm";


	public ServiceEngineersOptimizationGUROBI(double lambda, double[] mu, double[] alpha, double lostCost, 
			double[] engineerPartCost, int[] trunc_level, boolean logging) {
		super();
		try {
			N = mu.length - 1;
		} catch (Exception e) {
			N=-1;
		}
		M = trunc_level;
		this.lambda = lambda;
		this.mu = mu;
		this.alpha = alpha;
		
		this.lostCost = lostCost;
		this.engineerPartCost = engineerPartCost;
		
		M_max = 1; 
		p_mar = new double[N+1][];
		for(int i=0; i<=N; i++){ 
			M_max *= (M[i]+1);
			p_mar[i] = new double[M[i]+1];
		}

		this.logging = logging;
		

	}

	
	public static void main(String[] args){
	
		
	
	
	};
	
	public void formLP(int[] indicators_limits) throws GRBException{
	    env   = new GRBEnv();//"mip1.log");
	    model = new GRBModel(env);

	    GRBLinExpr[] eq_constraint = new GRBLinExpr[M_max+1];
        
        p_var 	   = new GRBVar[M_max+1];	
        y_var 	   = new GRBVar[M_max+1];	
        
    	for(int i=0; i<M_max; i++){
        	int[] n = getIndices(i);
        	String varIndex = "";  for(int j=0; j<=N; j++) varIndex += n[j];
        	p_var[i] = model.addVar(0, 1, lostCost, GRB.CONTINUOUS, "p" + varIndex);
        	y_var[i] = model.addVar(0, 1, -lostCost, GRB.CONTINUOUS, "y" + varIndex);
        }

        I_var      = new GRBVar[N+1][];	
    	for(int i=0; i<=N; i++){
	        I_var[i]      = new GRBVar[M[i]+1];	
        	for(int j=0; j<=M[i]; j++){
            	int lb=0,ub=1;
            	if(indicators_limits!=null && indicators_limits.length>N){
            		if(j<=indicators_limits[i]) lb=1;
            		else 						ub=0;
            		
            	}
        		double cost = (j>0) ? engineerPartCost[i] : 0;
        		I_var[i][j] = model.addVar(lb, ub, cost, GRB.BINARY, "I"+i+j);
        	}
    	}
    	model.update();
        
        GRBLinExpr objective = new GRBLinExpr();
			
//    	for(int i=0; i<=N; i++)
//        	for(int j=0; j<=M; j++){
//        		int i1 = (j>0) ? 1 : 0;
//        		objective.addTerm(i1*engineerPartCost[i], I_var[i][j]);
//        	}

    	GRBLinExpr norm_constraint = new GRBLinExpr();
    	for(int i=0; i<M_max; i++){
        	int[] n = getIndices(i);
        	String varIndex = "";  for(int j=0; j<=N; j++) varIndex += n[j];
        	
        	//objective.addTerm(lostCost, p_var[i]);
           	//objective.addTerm(-lostCost, y_var[i]);

        	norm_constraint.addTerm(1.0, p_var[i]);
   
        	//equilibrium constraints
        	eq_constraint[i] = new GRBLinExpr();
        	
        	eq_constraint[i].addTerm(lambda, y_var[i]);
        	double coeff = 0;
        	for(int j=0; j<=N; j++){
        		coeff += n[j]*mu[j];
        		n[j]--;
            	int i1 = getIndex(n);
            	if(i1>=0)
            		eq_constraint[i1].addTerm(-(n[j]+1)*mu[j], p_var[i]);
            	n[j]++;
        	}
        	eq_constraint[i].addTerm(coeff, p_var[i]);
        	
        	n[0]--;
        	for(int j=1; j<=N; j++){
        		n[j]--;
            	int i1 = getIndex(n);
            	if(i1>=0)
            		eq_constraint[i].addTerm(-lambda*alpha[j], y_var[i1]);
            	n[j]++;
        	}
        	n[0]++;        	
        	
        	// y_{ne,ns} <= p_{ne,ns}
        	{
        		GRBLinExpr tmp_const = new GRBLinExpr();
            	tmp_const.addTerm(-1.0, p_var[i]);
	        	tmp_const.addTerm( 1.0, y_var[i]);
	        	model.addConstr(tmp_const, GRB.LESS_EQUAL, 0.0, "yLep_const" + varIndex);
        	}
        	
        	// p_{ne,ns} <= I_ne and p_{ne,ns} <= I_ns
        	for(int i1=0; i1<=N; i1++){
        		GRBLinExpr tmp_const = new GRBLinExpr();
            	tmp_const.addTerm(1.0, p_var[i]);
            	tmp_const.addTerm(-1.0, I_var[i1][n[i1]]);
            	model.addConstr(tmp_const, GRB.LESS_EQUAL, 0.0, "pIconst" + varIndex + i1);
        	}

        	// y_{ne,ns} <= I_{ne+1} and y_{ne,ns} <= I_{ns+1}
        	for(int i1=0; i1<=N; i1++){
        		GRBLinExpr tmp_const = new GRBLinExpr();
//            	model.addConstr(tmp_const, GRB.LESS_EQUAL, 0.0, "yIconst" + varIndex + i1);
            	tmp_const.addTerm(1.0, y_var[i]);
            	if(n[i1]<M[i1]) tmp_const.addTerm(-1.0, I_var[i1][n[i1]+1]);
            	model.addConstr(tmp_const, GRB.LESS_EQUAL, 0.0, "yIconst" + varIndex + i1);
            }

        	// y_{ne,ns} <= p_{ne,ns} - (N+1 - I_{ne+1} - sum I_{ns+1}) 
        	{
        		GRBLinExpr tmp_const = new GRBLinExpr();
            	tmp_const.addTerm(-1.0,  p_var[i]);
            	tmp_const.addTerm( 1.0,  y_var[i]);
	        	for(int i1=0; i1<=N; i1++)
	        		if(n[i1]<M[i1]) tmp_const.addTerm(-1.0,  I_var[i1][n[i1]+1]);
	        	model.addConstr(tmp_const, GRB.GREATER_EQUAL, -(N+1), "ypIconst" + varIndex);	
	        }

        	
        }
    	model.addConstr(norm_constraint, GRB.EQUAL, 1.0, "normConst");
    	for(int i=0; i<M_max; i++){
        	int[] n = getIndices(i);
        	String varIndex = "";  for(int j=0; j<=N; j++) varIndex += n[j];
        	model.addConstr(eq_constraint[i], GRB.EQUAL, 0.0, "eqConst"+varIndex);
    	}

    	
    	for(int i=0; i<=N; i++)
        	for(int j=0; j<=M[i]; j++){
        		GRBLinExpr tmp_const = new GRBLinExpr();
        		tmp_const.addTerm(1.0, I_var[i][j]);
        		if(j<M[i]) tmp_const.addTerm(-1.0, I_var[i][j+1]);
            	model.addConstr(tmp_const, GRB.GREATER_EQUAL, 0.0, "IIconst" + i+j);
        	}

        //model.setObjective(objective, GRB.MINIMIZE);
    	model.update();
        //exportModel();
	}
	
	public void exportModel() throws GRBException{
		model.write(workpath + "/model1.lp");
	}
	

	public void setStartSolution(int common_start) throws GRBException{
	     for (int i = 0, idx = 0; i <= N; ++i)
	         for (int j = 0; j <= M[i]; ++j) {
	             I_var[i][j].set(GRB.DoubleAttr.Start, common_start);
	         }
	}
	
	public void addIndicatorLimits(int[] ll, int[] ul) throws GRBException{
    	for(int i=0; i<=N; i++){
        	if(ll!=null)
	    		for(int j=0; j<=ll[i]; j++)
	        		I_var[i][j].set(GRB.DoubleAttr.LB, 1.0);
        	if(ul!=null)
	        	for(int j=ul[i]+1; j<=M[i]; j++)
	        		I_var[i][j].set(GRB.DoubleAttr.UB, 0.0);
    	}
    	model.update();
	}


	public void Optimize() throws GRBException {
        model.optimize();
    }
    
    
    public void printIndicators() throws GRBException{
    	for(int i=0; i<=N; i++){
        	for(int j=0; j<=M[i]; j++)
        		System.out.print(I_var[i][j].get(GRB.DoubleAttr.X) + "\t");
        	System.out.println();
    	}
    }
    public void printIndicatorSums() throws GRBException{
    	for(int i=0; i<=N; i++){
        	int sum = 0;
    		for(int j=1; j<=M[i]; j++)
        		sum += Math.round(I_var[i][j].get(GRB.DoubleAttr.X));
        	if(i==0) System.out.println("Number of Engineers: "  + sum);
        	else     System.out.println("Number of Parts "+i+":   " + sum);
    	}
    }
    public void printYvariables() throws GRBException{
    	for(int i=0; i<M_max; i++){
        	int[] n = getIndices(i);	
        	for(int j=0; j<=N; j++) System.out.print(n[j]);
        	System.out.println("\t" + y_var[i].get(GRB.DoubleAttr.X));
    	}
    }
    public void printPvariables() throws GRBException{
    	for(int i=0; i<M_max; i++){
        	int[] n = getIndices(i);	
        	for(int j=0; j<=N; j++) System.out.print(n[j]);
        	System.out.println("\t" + p_var[i].get(GRB.DoubleAttr.X));
    	}
    }

    public void printPvariables2D() throws GRBException{
    	if(N==1){
    		int[] n = new int[2];
    		for(int i=M[0]; i>=0; i--){
    			for(int j=M[1]; j>=0; j--){
    				n[0] = i; n[1] = j;
    				int i1 = getIndex(n);
    				System.out.print("\t" + p_var[i1].get(GRB.DoubleAttr.X));
    			}
    			System.out.println();
    		}
        }
    }

    public void printMarginalProbabilities() throws GRBException{
    	for(int i=0; i<M_max; i++){
        	int[] n = getIndices(i);	
        	for(int j=0; j<=N; j++) 
        		p_mar[j][n[j]] += p_var[i].get(GRB.DoubleAttr.X);
    	}
    	for(int j=0; j<=N; j++){
    		System.out.print("item " + j + ":");
    		for(int i=0; i<=M[j]; i++)
    			System.out.print("\t" + p_mar[j][i]);
    		System.out.println();
    	}
    }

//    public double getDual(int rowIndex) throws IloException{
//    	return model.getDual(constraint[rowIndex]);
//    }
//
//    public double getSlack(int rowIndex) throws IloException{
//    	return model.getSlack(constraint[rowIndex]);
//    }
//    
    public void cleanupModel() throws GRBException{
    	model.dispose();
    	env.dispose();
    }

    public void setParameters() throws GRBException{
    	setParametersFromFile(defaultParamFile);
    }

    
    public void setParametersFromFile(String fixedfile) throws GRBException{
        if ( fixedfile != null ) {
           env.readParams(fixedfile);
        }
    }
    public void tuneModel(String fixedfile) throws GRBException{
        if ( fixedfile != null && fixedfile.contains(".prm")) {
           	model.tune();
           	//env.writeParams(fixedfile);
           	model.write(fixedfile);
           	System.out.println("Tuned parameters written to file '" +
            		fixedfile + "'");
         }
    }
	
    public void tuneModel() throws GRBException{
    	tuneModel(defaultParamFile);
    }

    public int getIndex(int[] n){
		int k = 0, k1=M_max;
		for(int i=0; i<=N; i++){
			if(n[i]>=0){
				k1 /= (M[i]+1);
				k  += k1*n[i];
			}
			else
				return -1;
		}
		
		return k;
	}
	
	public int[] getIndices(int i){
		int[] tempk = new int[N+1];
		
		int k = i, k1=M_max;
		for(int j=0; j<=N; j++){
			k1 = k1 / (M[j]+1);
			tempk[j] = k/k1;
			k -= tempk[j]*k1;
		}
		
		return tempk;
	}
	
	
}
