package OptimizationProblem;

import gurobi.*;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.UnsupportedEncodingException;
import java.util.Date;
import java.util.StringTokenizer;

public class ServiceEngineersOptimizationGUROBI {
	
	int N; //number of spare types
	int[] M_lb, M_ub; //state space truncation limit
	double[] lambda; //failure rate;
	double[] mu;   //service lead time service rates, mu[0] corresponds to engenders, the rest to spares
	double[] lostCost;
	double engnCost;
	double[] engineerPartCost;
	
	int M_max;
	boolean logging;
	
	
	GRBVar[][] I_var;
	GRBVar[] p_var;
	GRBVar[][] y_var;
	//GRBVar[][][] z_var;
	int[] optIsum;
	int[] ll, ul;

	
	double[][] p_mar;
	
    GRBModel  model;
	//IloCplex model;
    GRBEnv    env;

    Date startTime, endTime;
    long computation_time;

    
    public static String workpath;// = "D:/Users/as14446/Documents/Service Engineers";
    public static String defaultParamFile;// = workpath + "/gurobiParameters.prm";


	public ServiceEngineersOptimizationGUROBI(int N, double lambda[], double[] mu, double[] lostCost, 
			double[] engineerPartCost, int[] trunc_level_lb, int[] trunc_level_ub, boolean logging) {
		super();
		this.N = N;
//		try {
//			N = mu.length - 1;
//		} catch (Exception e) {
//			N=-1;
//		}
		M_lb = trunc_level_lb;
		M_ub = trunc_level_ub;
		this.lambda = lambda;
		this.mu = mu;
		
		this.lostCost = lostCost;
		this.engineerPartCost = engineerPartCost;
		
		M_max = 1; 
		p_mar = new double[N+1][];
		for(int i=0; i<=N; i++){ 
			M_max *= (M_ub[i]+1);
			p_mar[i] = new double[M_ub[i]+1];
		}

		this.logging = logging;
		
		this.optIsum = null;
		
		this.ll=null;
		this.ul=null;
		

		

	}

	
	public static void main(String[] args){
	
		
		int N = 0; //number of spare types
		int[] M_lb = null, M_ub=null; //state space truncation limits
		//double lambdaTot = 0; //failure rate;
		double[] mu = null;   //service lead time service rates, mu[0] corresponds to engenders, the rest to spares
		double[] lambda = null; //probability that spare i is requested
		double[] engineerPartCost = null;
		double[] lostCost = null;
		boolean logging = true;
		
			
		String inputFileName = args[0]; //"D:/Users/as14446/Documents/Service Engineers/input.txt"; //args[0];
		

		
		BufferedReader input = null;
		try {
			File file = new File(inputFileName);
			String absolutePath = file.getAbsolutePath();
		    workpath = absolutePath.substring(0,absolutePath.lastIndexOf(File.separator));
		    defaultParamFile = workpath + "/gurobiParameters.prm";


		    input = new BufferedReader(new InputStreamReader(new FileInputStream(inputFileName), "utf-8"));
		} catch (UnsupportedEncodingException | FileNotFoundException e) {
			System.err.println("Problem with input file: openning file");
			e.printStackTrace();
	        System.exit(1);
		}
	   	
		//read lambda
	    StringTokenizer s1 = null;
		try {
			//read N
		    s1 = new StringTokenizer(input.readLine(), "\t");
		    s1.nextToken(); N = Integer.parseInt(s1.nextToken());
	        //read M[i]
		    M_lb = new int[N+1];
		    s1 = new StringTokenizer(input.readLine(), "\t");
		    s1.nextToken(); for(int i=0; i<=N; i++) M_lb[i] = Integer.parseInt(s1.nextToken());
		    M_ub = new int[N+1];
		    s1 = new StringTokenizer(input.readLine(), "\t");
		    s1.nextToken(); for(int i=0; i<=N; i++) M_ub[i] = Integer.parseInt(s1.nextToken());
			//read lambda
		    lambda = new double[N+1];
		    s1 = new StringTokenizer(input.readLine(), "\t");
		    s1.nextToken(); for(int i=0; i<=N; i++) lambda[i] = Double.parseDouble(s1.nextToken());
		    lambda[0] = 0; for(int i=1; i<=N; i++)  lambda[0] += lambda[i];
	        //read mu
		    mu = new double[N+1];
		    s1 = new StringTokenizer(input.readLine(), "\t");
		    s1.nextToken(); for(int i=0; i<=N; i++) mu[i] = Double.parseDouble(s1.nextToken());
			// //read lost cost
		    // s1 = new StringTokenizer(input.readLine(), "\t");
		    // s1.nextToken(); lostCost = Double.parseDouble(s1.nextToken());
	        //read lost costs
		    lostCost = new double[N+1];
		    s1 = new StringTokenizer(input.readLine(), "\t");
		    s1.nextToken(); for(int i=0; i<=N; i++) lostCost[i] = Double.parseDouble(s1.nextToken());
	        //read costs
		    engineerPartCost = new double[N+1];
		    s1 = new StringTokenizer(input.readLine(), "\t");
		    s1.nextToken(); for(int i=0; i<=N; i++) engineerPartCost[i] = Double.parseDouble(s1.nextToken());
		    //read logging flag
		    s1 = new StringTokenizer(input.readLine(), "\t");
		    s1.nextToken(); logging=Boolean.parseBoolean(s1.nextToken());
		} catch (IOException e) {
			System.err.println("Problem with input file: reading file, check the format");
			e.printStackTrace();
	        System.exit(1);		
	    }		
		

		
		for(int nparts=1; nparts<=6; nparts++ )
		{
			ServiceEngineersOptimizationGUROBI opt = new ServiceEngineersOptimizationGUROBI(nparts, lambda, mu, lostCost, engineerPartCost, M_lb, M_ub, logging);
		try {
			
			opt.StartTimer();
			opt.ElapsedTime("Starting");
			
			opt.formLP(null);
			opt.ElapsedTime("MIP is formed");
			
			opt.setStartSolution(1);
			opt.setParameters();
			
			int[] ll= new int[N+1];
			int[] ul= new int[N+1];
			for(int i=0; i<=N; i++){
				ll[i] = Math.max(0, M_lb[i]);
				ul[i] = Math.min(ll[i]+2, M_ub[i]);
			}
			opt.setIndicatorLimits(ll,ul);
			opt.exportModel();
			
			opt.setLoggingOff();
			
			opt.ElapsedTime("Start Solving");
			boolean ready = false;
			while(!ready){
				opt.exportModel();
				System.out.println("Curent objective: " + opt.optimize());
				opt.printIndicatorShort();
				ready = opt.doIteration();
				opt.ElapsedTime("");
			}

			opt.printIndicatorSums();
			opt.printLossProbability();

			opt.StopTimer();
			
			System.out.println("-- Optimal cost: " + opt.model.get(GRB.DoubleAttr.ObjVal));
			
			opt.cleanupModel();
		
		} catch (GRBException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}		
		}
	
	
	};
	
	public void formLP(int[] indicators_limits) throws GRBException{
	    env   = new GRBEnv();//"mip1.log");
	    model = new GRBModel(env);

	    GRBLinExpr[] eq_constraint = new GRBLinExpr[M_max+1];
        
        p_var 	   = new GRBVar[M_max+1];	
        y_var 	   = new GRBVar[N+1][M_max+1];	
        
        
        for(int i=0; i<M_max; i++){
        	int[] n = getIndices(i);
        	String varIndex = "";  for(int j=0; j<=N; j++) varIndex += n[j];
        	double tmpCost = 0;
        	for(int j=1; j<=N; j++){ 
        		tmpCost += lostCost[j]*lambda[j];
        		y_var[j][i] = model.addVar(0, 1, -lostCost[j]*lambda[j], GRB.CONTINUOUS, "y" + j + varIndex);
        	}
        	p_var[i] = model.addVar(0, 1, tmpCost, GRB.CONTINUOUS, "p" + varIndex);    
        	//y_var[0][i] = model.addVar(0, 1, -tmpCost, GRB.CONTINUOUS, "y0" + varIndex);
        }
        
        //z_var      = new GRBVar[N+1][][];
    	//for(int i=1; i<=N; i++){
    	//	z_var[i]      = new GRBVar[M_ub[0]+1][];	
        //	for(int j=0; j<=M_ub[0]; j++){
        //		z_var[i][j]      = new GRBVar[M_ub[i]+1];	
        //    	
        //		for(int k=0; k<=M_ub[i]; k++)
        //			z_var[i][j][k] = model.addVar(0, 1, -lostCost[i]*lambda[i], GRB.CONTINUOUS, "z" +i+j+k);
        //	}
    	//}        
    	

        I_var      = new GRBVar[N+1][];	
    	for(int i=0; i<=N; i++){
	        I_var[i]      = new GRBVar[M_ub[i]+1];	
        	for(int j=0; j<=M_ub[i]; j++){
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
        
//        GRBLinExpr objective = new GRBLinExpr();
			
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
        	
        	//eq_constraint[i].addTerm(lambda[0], y_var[0][i]);
        	for(int j=1; j<=N; j++){
            	eq_constraint[i].addTerm(lambda[j], y_var[j][i]);
            	//eq_constraint[i].addTerm(-lambda[j], z_var[j][n[0]][n[j]]);
        	}
        	
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
            	if(i1>=0){
            		//eq_constraint[i].addTerm(-lambda[j], y_var[0][i1]);
            		eq_constraint[i].addTerm(-lambda[j], y_var[j][i1]);
            		//eq_constraint[i].addTerm(lambda[j], z_var[j][n[0]][n[j]]);
            	}
            	n[j]++;
        	}
        	n[0]++;        	
        	
        	// y_{ne,ns} <= p_{ne,ns}
        	for(int j=1; j<=N; j++){
        		GRBLinExpr tmp_const = new GRBLinExpr();
            	tmp_const.addTerm(-1.0, p_var[i]);
	        	tmp_const.addTerm( 1.0, y_var[j][i]);
	        	model.addConstr(tmp_const, GRB.LESS_EQUAL, 0.0, "yLep_const" + j + varIndex);
        	}
        	
        	// p_{ne,ns} <= I_ne and p_{ne,ns} <= I_ns
        	for(int i1=0; i1<=N; i1++){
        		GRBLinExpr tmp_const = new GRBLinExpr();
            	tmp_const.addTerm(1.0, p_var[i]);
            	tmp_const.addTerm(-1.0, I_var[i1][n[i1]]);
            	model.addConstr(tmp_const, GRB.LESS_EQUAL, 0.0, "pIconst" + varIndex + i1);
        	}

        	// y_{ne,ns} <= I_{ne+1} and y_{ne,ns} <= I_{ns+1}
        	for(int i1=1; i1<=N; i1++){
        		GRBLinExpr tmp_const1 = new GRBLinExpr();
        		GRBLinExpr tmp_const2 = new GRBLinExpr();
//            	model.addConstr(tmp_const, GRB.LESS_EQUAL, 0.0, "yIconst" + varIndex + i1);
            	tmp_const1.addTerm(1.0, y_var[i1][i]);
            	tmp_const2.addTerm(1.0, y_var[i1][i]);
            	if(n[i1]<M_ub[i1]) tmp_const1.addTerm(-1.0, I_var[i1][n[i1]+1]);
            	if(n[0]<M_ub[0])   tmp_const2.addTerm(-1.0, I_var[0][n[0]+1]);
            	model.addConstr(tmp_const1, GRB.LESS_EQUAL, 0.0, "y0const" + i1+ varIndex);
            	model.addConstr(tmp_const2, GRB.LESS_EQUAL, 0.0, "yIconst" + i1+ varIndex);
            }

        	// y_{ne,ns} <= p_{ne,ns} - (N+1 - I_{ne+1} - sum I_{ns+1}) 
        	for(int j=1; j<=N; j++){
        		GRBLinExpr tmp_const = new GRBLinExpr();
            	tmp_const.addTerm(-1.0,  p_var[i]);
            	tmp_const.addTerm( 1.0,  y_var[j][i]);
    			if(n[0] < M_ub[0]) tmp_const.addTerm(-1.0,  I_var[0][n[0]+1]);
    			if(n[j] < M_ub[j]) tmp_const.addTerm(-1.0,  I_var[j][n[j]+1]);
	        	model.addConstr(tmp_const, GRB.GREATER_EQUAL, -2, "ypIconst" + j + varIndex);	
	        }

        	
        }
    	model.addConstr(norm_constraint, GRB.EQUAL, 1.0, "normConst");
    	for(int i=0; i<M_max; i++){
        	int[] n = getIndices(i);
        	String varIndex = "";  for(int j=0; j<=N; j++) varIndex += n[j];
        	model.addConstr(eq_constraint[i], GRB.EQUAL, 0.0, "eqConst"+varIndex);
    	}
    	
//    	for(int i=1; i<=N; i++){
//        	for(int j=0; j<=M_ub[0]; j++){
//        		for(int k=0; k<=M_ub[i]; k++){
//        			GRBLinExpr tmp_const1 = new GRBLinExpr();
//        			GRBLinExpr tmp_const2 = new GRBLinExpr();
//        			GRBLinExpr tmp_const3 = new GRBLinExpr();
//        			
//        			tmp_const1.addTerm(1.0,  z_var[i][j][k]);
//        			tmp_const2.addTerm(1.0,  z_var[i][j][k]);
//        			tmp_const3.addTerm(1.0,  z_var[i][j][k]);
//
//        			for(int i1=0; i1<M_max; i1++){
//        				int[] n = getIndices(i1);
//        				if(n[0]==j && n[i]==k) tmp_const1.addTerm(-1.0,  p_var[i1]);
//        				if(n[0]==j && n[i]==k) tmp_const2.addTerm(-1.0,  p_var[i1]);
//        			}
//        			if(j < M_ub[0]) tmp_const2.addTerm(-1.0,  I_var[0][j+1]);
//        			if(k < M_ub[i]) tmp_const2.addTerm(-1.0,  I_var[i][k+1]);
//        			
//        			if(j < M_ub[0]) tmp_const3.addTerm(-1.0,  I_var[0][j+1]);
//        			if(k < M_ub[i]) tmp_const3.addTerm(-1.0,  I_var[i][k+1]);
//        			
//        			model.addConstr(tmp_const1, GRB.LESS_EQUAL, 0.0, "Zconst1" + i+j+k);
//        			model.addConstr(tmp_const2, GRB.GREATER_EQUAL, -2.0, "Zconst2" + i+j+k);
//        			model.addConstr(tmp_const3, GRB.LESS_EQUAL, 0.0, "Zconst3" + i+j+k);
//        		}
//        	}
//    	}        


    	
    	for(int i=0; i<=N; i++)
        	for(int j=0; j<=M_ub[i]; j++){
        		GRBLinExpr tmp_const = new GRBLinExpr();
        		tmp_const.addTerm(1.0, I_var[i][j]);
        		if(j<M_ub[i]) tmp_const.addTerm(-1.0, I_var[i][j+1]);
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
	     for (int i = 0; i <= N; ++i)
	         for (int j = 0; j <= M_ub[i]; ++j) {
	             I_var[i][j].set(GRB.DoubleAttr.Start, common_start);
	         }
	}
	
	public void addIndicatorLimits(int[] ll, int[] ul) throws GRBException{
    	for(int i=0; i<=N; i++){
        	if(ll!=null)
	    		for(int j=0; j<=ll[i]; j++)
	        		I_var[i][j].set(GRB.DoubleAttr.LB, 1.0);
        	if(ul!=null)
	        	for(int j=ul[i]+1; j<=M_ub[i]; j++)
	        		I_var[i][j].set(GRB.DoubleAttr.UB, 0.0);
    	}
    	model.update();
	}
	
	public void setIndicatorLimits(int[] ll, int[] ul) throws GRBException{
    	this.ll = ll;
    	this.ul = ul;
		for(int i=0; i<=N; i++){
    		for(int j=0; j<=M_ub[i]; j++){
        		I_var[i][j].set(GRB.DoubleAttr.LB, 0.0);
        		I_var[i][j].set(GRB.DoubleAttr.UB, 1.0);
    		}
        	if(ll!=null)
	    		for(int j=0; j<=this.ll[i]; j++)
	        		I_var[i][j].set(GRB.DoubleAttr.LB, 1.0);
        	if(ul!=null)
	        	for(int j=this.ul[i]+1; j<=M_ub[i]; j++)
	        		I_var[i][j].set(GRB.DoubleAttr.UB, 0.0);
    	}
		model.update();
	}


	public double optimize() throws GRBException {
    	model.optimize();
		if(optIsum==null) optIsum = new int[N+1];
		for(int i=0; i<=N; i++){
        	optIsum[i] = 0;
    		for(int j=1; j<=M_ub[i]; j++) optIsum[i] += Math.round(I_var[i][j].get(GRB.DoubleAttr.X));
    	}
		return model.get(GRB.DoubleAttr.ObjVal);
    }
    
	public boolean doIteration() throws GRBException{
		boolean ready = true;
		for(int i=0; i<=N; i++){
			if(optIsum[i] == ul[i] && ul[i]<M_ub[i]){
				ul[i]++; ll[i]++;
				ready =  false;
			}else if(optIsum[i] == ll[i] && ll[i]>M_lb[i]){
				ul[i]--; ll[i]--;
				ready =  false;
			}
		}
		if(!ready) this.setIndicatorLimits(this.ll, this.ul);
		return ready;
	}
    
    public void printIndicators() throws GRBException{
    	for(int i=0; i<=N; i++){
        	for(int j=0; j<=M_ub[i]; j++)
        		System.out.print(I_var[i][j].get(GRB.DoubleAttr.X) + "\t");
        	System.out.println();
    	}
    }
    public void printIndicatorSums() throws GRBException{
    	for(int i=0; i<=N; i++){
        	int sum = 0;
    		for(int j=1; j<=M_ub[i]; j++)
        		sum += Math.round(I_var[i][j].get(GRB.DoubleAttr.X));
        	if(i==0) System.out.println("Number of Engineers: "  + sum);
        	else     System.out.println("Number of Parts "+i+":   " + sum);
    	}
    }
    public void printIndicatorShort(){
    	for(int i=0; i<=N; i++) System.out.print(optIsum[i] + "\t");
    	System.out.println();
    }
    public void printYvariables() throws GRBException{
    	for(int i=0; i<M_max; i++){
        	int[] n = getIndices(i);	
        	for(int j=0; j<=N; j++) System.out.print(n[j]);
        	for(int j=0; j<=N; j++) System.out.print("\t" + y_var[j][i].get(GRB.DoubleAttr.X));
        	System.out.println();
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
    		for(int i=M_ub[0]; i>=0; i--){
    			for(int j=M_ub[1]; j>=0; j--){
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
    		for(int i=0; i<=M_ub[j]; i++)
    			System.out.print("\t" + p_mar[j][i]);
    		System.out.println();
    	}
    }
    
    public void printLossProbability() throws GRBException{
    	double[] p_loss= new double[N+1];
    	for(int i=0; i<M_max; i++){
        	int[] n = getIndices(i);	
        	for(int j=1; j<=N; j++)
        		if(n[j]==optIsum[j] || n[0]==optIsum[0]){
        			p_loss[j] += p_var[i].get(GRB.DoubleAttr.X);
        		}
        	boolean loss = false;
        	for(int j=0; j<=N; j++)	if(n[j]==optIsum[j]) loss = true;
			if(loss) p_loss[0] += p_var[i].get(GRB.DoubleAttr.X);
    	}
    	System.out.println("Total Loss probability: "  + p_loss[0]);
    	for(int j=1; j<=N; j++)
        	System.out.println("Loss probability " + j +": "  + p_loss[j]);
    	
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
				k1 /= (M_ub[i]+1);
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
			k1 = k1 / (M_ub[j]+1);
			tempk[j] = k/k1;
			k -= tempk[j]*k1;
		}
		
		return tempk;
	}
	
    public void StartTimer(){
    	startTime = new Date();
    }

    public void ElapsedTime(String print_line){
   	   Date curTime = new Date();
   	   if(print_line!="") System.out.print(curTime + ": " + print_line + " ");
  	   System.out.println("Elapsed Time: " + (curTime.getTime() - startTime.getTime())/1000 + " sec.");
      }


    public void StopTimer(){
 	   endTime = new Date();
	   computation_time = (endTime.getTime() - startTime.getTime());
	   System.out.println("Optimization Time: " + computation_time/1000 + " sec.");
    }

	public void setLoggingOff() throws GRBException{
		if(model!=null){
			model.getEnv().set(GRB.IntParam.OutputFlag, 0);
		}
	}


	
}
