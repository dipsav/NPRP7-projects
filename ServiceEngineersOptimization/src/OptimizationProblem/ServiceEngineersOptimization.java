package OptimizationProblem;

import ilog.concert.IloColumn;
import ilog.concert.IloException;
import ilog.concert.IloLinearNumExpr;
import ilog.concert.IloNumVar;
import ilog.concert.IloNumVarType;
import ilog.concert.IloObjective;
import ilog.concert.IloRange;
import ilog.cplex.IloCplex;
import ilog.cplex.IloCplex.UnknownObjectException;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.UnsupportedEncodingException;
import java.util.Date;
import java.util.StringTokenizer;

public class ServiceEngineersOptimization {
	
	int N; //number of spare types
	int[] M_lb, M_ub; //state space truncation limits
	double lambda; //failure rate;
	double[] mu;   //service lead time service rates, mu[0] corresponds to engenders, the rest to spares
	double[] alpha; //probability that spare i is requested
	double lostCost;
	double engnCost;
	double[] engineerPartCost;
	
	int M_max;
	boolean logging;
	
	
	IloNumVar[][] I_var;
	int[] optIsum, prevOptIsum;
	int[] ll, ul;
		
	IloNumVar[] p_var, y_var;
	
	double[][] p_mar;
	
	IloCplex model;
    
    public static String workpath = "D:/Users/as14446/Documents/Service Engineers";
    static String defaultParamFile = workpath + "/cplexParameters";

    Date startTime, endTime;
    long computation_time;
    
    public ServiceEngineersOptimization(double lambda, double[] mu, double[] alpha, double lostCost, 

			double[] engineerPartCost, int[] trunc_level_lb, int[] trunc_level_ub, boolean logging) {
		super();
		try {
			N = trunc_level_ub.length - 1;
		} catch (Exception e) {
			N=-1;
		}
		this.M_ub = trunc_level_ub;
		this.M_lb = trunc_level_lb;
		this.lambda = lambda;
		this.mu = mu;
		this.alpha = alpha;
		
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
		this.prevOptIsum = null;
		
		this.ll=null;
		this.ul=null;
		

	}
	
	public static void main(String[] args){
	
		int N = 0; //number of spare types
		int[] M_lb = null, M_ub=null; //state space truncation limits
		double lambda = 0; //failure rate;
		double[] mu = null;   //service lead time service rates, mu[0] corresponds to engenders, the rest to spares
		double[] alpha = null; //probability that spare i is requested
		double[] engineerPartCost = null;
		double lostCost = 0;
		boolean logging = true;
		
			
		String inputFileName = "D:/Users/as14446/Documents/Service Engineers/input.txt"; //args[0];
		BufferedReader input = null;
		try {
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
	        //read M_lb[i]
		    M_lb = new int[N+1];
		    s1 = new StringTokenizer(input.readLine(), "\t");
		    s1.nextToken(); for(int i=0; i<=N; i++) M_lb[i] = Integer.parseInt(s1.nextToken());
		    //read M_ub[i]
		    M_ub = new int[N+1];
		    s1 = new StringTokenizer(input.readLine(), "\t");
		    s1.nextToken(); for(int i=0; i<=N; i++) M_ub[i] = Integer.parseInt(s1.nextToken());
			//read lambda
		    s1 = new StringTokenizer(input.readLine(), "\t");
		    s1.nextToken(); lambda = Double.parseDouble(s1.nextToken());
	        //read alpha
		    alpha = new double[N+1];
		    s1 = new StringTokenizer(input.readLine(), "\t");
		    s1.nextToken(); for(int i=0; i<=N; i++) alpha[i] = Double.parseDouble(s1.nextToken());
	        //read mu
		    mu = new double[N+1];
		    s1 = new StringTokenizer(input.readLine(), "\t");
		    s1.nextToken(); for(int i=0; i<=N; i++) mu[i] = Double.parseDouble(s1.nextToken());
			//read lost cost
		    s1 = new StringTokenizer(input.readLine(), "\t");
		    s1.nextToken(); lostCost = Double.parseDouble(s1.nextToken());
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
		
		ServiceEngineersOptimization opt = new ServiceEngineersOptimization(lambda, mu, alpha, lostCost, engineerPartCost, M_lb, M_ub, logging);

		try {
			opt.StartTimer();
			opt.ElapsedTime("Starting");
			
			opt.formLPbyRows();
			opt.ElapsedTime("MIP is formed");
			
			opt.setStartSolution();
			opt.setParameters();
			
			int[] ll= new int[N+1];
			int[] ul= new int[N+1];
			for(int i=0; i<=N; i++){
				ll[i] = Math.max(0, M_lb[i]);
				ul[i] = Math.min(ll[i]+2, M_ub[i]);
			}
			opt.setIndicatorLimits(ll,ul);
			opt.exportModel();
			
			opt.ElapsedTime("Start Solving");
			boolean ready = false;
			while(!ready){
				//obj.exportModel();
				System.out.println("Curent objective: " + opt.optimize());
				opt.printIndicatorShort();
				ready = opt.doIteration();
				opt.ElapsedTime("");
			}

			opt.printIndicatorSums();
			opt.printLossProbability();

			opt.StopTimer();
			
			opt.cleanupModel();
		
		} catch (IloException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	
	};
	
	public void formLP() throws IloException{
	    IloRange[] constraint;
	    IloColumn[] p_column, y_column;
	    IloColumn[][] I_column;

	    model = new IloCplex();
		if(!logging){
			model.setOut(null);
			model.setWarning(null);
		}
        constraint = new IloRange[M_max+1];
        p_column   = new IloColumn[M_max+1];
        y_column   = new IloColumn[M_max+1];
        
        p_var 	   = new IloNumVar[M_max+1];	
        y_var 	   = new IloNumVar[M_max+1];	

        I_column   = new IloColumn[N+1][];	
        I_var      = new IloNumVar[N+1][];	
		for(int i=0; i<=N; i++){ 
	        I_column[i]   = new IloColumn[M_ub[i]+1];	
	        I_var[i]      = new IloNumVar[M_ub[i]+1];	
		}
        
        IloObjective cost = model.addMinimize();
        //IloLinearNumExpr cost = model.linearNumExpr();
        
    	for(int i=0; i<=N; i++)
        	for(int j=0; j<=M_ub[i]; j++){
        		int i1 = (j>0) ? 1 : 0;
        		I_column[i][j] = model.column(cost, i1*engineerPartCost[i]);
        	}

    	IloRange norm_constraint = model.addEq(null, 1.0, "normConst");
    	//IloLinearNumExpr norm_constraint = model.linearNumExpr();
    	
    	for(int i=0; i<M_max; i++){
        	int[] n = getIndices(i);
        	String varIndex = "";  for(int j=0; j<=N; j++) varIndex += n[j];
        	
        	p_var[i] = model.numVar(0.0, 1.0, IloNumVarType.Float, "p" + varIndex);
        	y_var[i] = model.numVar(0.0, 1.0, IloNumVarType.Float, "y" + varIndex);
        	
    		p_column[i] = model.column(cost, lostCost); //TODO update coefficient
        	
        	y_column[i] = model.column(cost, -lostCost); //TODO update coefficient

        	//normalization constraint
        	p_column[i] = p_column[i].and( model.column(norm_constraint, 1.0));
        	//norm_constraint.addTerm(1.0, p_var[i]);
        	
        	//equilibrium constraints
        	constraint[i] = model.addEq(null, 0.0, "eqConst"+varIndex);
        	
        	y_column[i] = y_column[i].and(model.column(constraint[i], lambda));
        	double coeff = 0;
        	for(int j=0; j<=N; j++){
        		coeff += n[j]*mu[j];
        		n[j]--;
            	int i1 = getIndex(n);
            	if(i1>=0)
            		p_column[i]  = p_column[i].and( model.column(constraint[i1], -(n[j]+1)*mu[j]));
            	n[j]++;
        	}
        	
        	p_column[i] = p_column[i].and(model.column(constraint[i], coeff));
        	
        	n[0]--;
        	for(int j=1; j<=N; j++){
        		n[j]--;
            	int i1 = getIndex(n);
            	if(i1>=0)
            		y_column[i1] = y_column[i1].and(model.column(constraint[i], -lambda*alpha[j]));
            	n[j]++;
        	}
        	n[0]++;        	
        	
        	// y_{ne,ns} <= p_{ne,ns}
        	{
            	IloRange tmp_const = model.addLe(null, 0.0, "yLep_const" + varIndex);
            	p_column[i] = p_column[i].and(model.column(tmp_const, -1.0));
	        	y_column[i] = y_column[i].and(model.column(tmp_const, 1.0));
        	}
        	
        	// p_{ne,ns} <= I_ne and p_{ne,ns} <= I_ns
        	for(int i1=0; i1<=N; i1++){
            	IloRange tmp_const = model.addLe(null, 0.0, "pIconst" + varIndex + i1);
            	p_column[i] = p_column[i].and(model.column(tmp_const, 1.0));
            	I_column[i1][n[i1]] = I_column[i1][n[i1]].and(model.column(tmp_const, -1.0));
            }

        	// y_{ne,ns} <= I_{ne+1} and p_{ne,ns} <= I_{ns+1}
        	for(int i1=0; i1<=N; i1++){
            	IloRange tmp_const = model.addLe(null, 0.0, "yIconst" + varIndex + i1);
            	y_column[i] = y_column[i].and(model.column(tmp_const, 1.0));
            	if(n[i1]<M_ub[i1]) I_column[i1][n[i1]+1] = I_column[i1][n[i1]+1].and(model.column(tmp_const, -1.0));
            }

        	// y_{ne,ns} <= p_{ne,ns} - (N+1 - I_{ne+1} - sum I_{ns+1}) 
        	{
        		IloRange tmp_const = model.addGe(null, -(N+1),  "ypIconst" + varIndex);
	        	p_column[i] = p_column[i].and(model.column(tmp_const, -1.0));
	        	y_column[i] = y_column[i].and(model.column(tmp_const, 1.0));
	        	for(int i1=0; i1<=N; i1++)
		        	if(n[i1]<M_ub[i1]) I_column[i1][n[i1]+1] = I_column[i1][n[i1]+1].and(model.column(tmp_const, -1.0));
        	}

        	
        }
    	
    	for(int i=0; i<M_max; i++){
        	int[] n = getIndices(i);
        	String varIndex = "";  for(int j=0; j<=N; j++) varIndex += n[j];
        	p_var[i] = model.numVar(p_column[i], 0, 1, IloNumVarType.Float, "p" + varIndex);
        	y_var[i] = model.numVar(y_column[i], 0, 1, IloNumVarType.Float, "y" + varIndex);
        }

    	for(int i=0; i<=N; i++)
        	for(int j=0; j<=M_ub[i]; j++){
        		IloRange tmp_const = model.addGe(null, 0.0, "IIconst" + i+j);
        		I_column[i][j] = I_column[i][j].and(model.column(tmp_const, 1.0));
        		if(j<M_ub[i])
        			I_column[i][j+1] = I_column[i][j+1].and(model.column(tmp_const, -1.0));

            	I_var[i][j] = model.numVar(I_column[i][j], 0, 1, IloNumVarType.Bool, "I"+i+j);
        	}

    	//model.exportModel(workpath + "/model1.lp");
        	
	}
	
	public void formLPbyRows() throws IloException{
	    IloLinearNumExpr[] constraint;

	    model = new IloCplex();
		if(!logging){
			model.setOut(null);
			model.setWarning(null);
		}
        constraint = new IloLinearNumExpr[M_max+1];
        
        p_var 	   = new IloNumVar[M_max+1];	
        y_var 	   = new IloNumVar[M_max+1];	

        I_var      = new IloNumVar[N+1][];	
		for(int i=0; i<=N; i++){ 
	        I_var[i]      = new IloNumVar[M_ub[i]+1];	
        	for(int j=0; j<=M_ub[i]; j++)
            	I_var[i][j] = model.numVar(0, 1, IloNumVarType.Bool, "I"+i + "" + j);
		}
        
        IloLinearNumExpr cost = model.linearNumExpr();
        
    	for(int i=0; i<=N; i++)
        	for(int j=0; j<=M_ub[i]; j++){
        		int i1 = (j>0) ? 1 : 0;
        		cost.addTerm(i1*engineerPartCost[i], I_var[i][j]);
        	}

    	IloLinearNumExpr norm_constraint = model.linearNumExpr();
    	
    	for(int i=0; i<M_max; i++){
        	int[] n = getIndices(i);
        	String varIndex = "";  for(int j=0; j<=N; j++) varIndex += n[j];
        	
        	p_var[i] = model.numVar(0.0, 1.0, IloNumVarType.Float, "p" + varIndex);
        	y_var[i] = model.numVar(0.0, 1.0, IloNumVarType.Float, "y" + varIndex);
        	
        	cost.addTerm(lostCost, p_var[i]); //TODO update coefficient
        	cost.addTerm(-lostCost, y_var[i]); //TODO update coefficient
        	
        	//normalization constraint
        	norm_constraint.addTerm(1.0, p_var[i]);
        	
        	//equilibrium constraints
        	constraint[i] = model.linearNumExpr();
        	
        	constraint[i].addTerm(lambda, y_var[i]);
        	double coeff = 0;
        	for(int j=0; j<=N; j++){
        		coeff += n[j]*mu[j];
        		n[j]--;
            	int i1 = getIndex(n);
            	if(i1>=0)
            		constraint[i1].addTerm(-(n[j]+1)*mu[j], p_var[i]);
            	n[j]++;
        	}
        	constraint[i].addTerm(coeff, p_var[i]);
        	
        	n[0]--;
        	for(int j=1; j<=N; j++){
        		n[j]--;
            	int i1 = getIndex(n);
            	if(i1>=0)
            		constraint[i].addTerm(-lambda*alpha[j], y_var[i1]);
            	n[j]++;
        	}
        	n[0]++;        	
        	
        	// y_{ne,ns} <= p_{ne,ns}
        	{
            	IloLinearNumExpr tmp_const = model.linearNumExpr();
            	tmp_const.addTerm(-1.0, p_var[i]);
            	tmp_const.addTerm( 1.0, y_var[i]);
            	model.addLe(tmp_const, 0.0, "yLep_const" + varIndex);
        	}
        	
        	// p_{ne,ns} <= I_ne and p_{ne,ns} <= I_ns
        	for(int i1=0; i1<=N; i1++){
        		IloLinearNumExpr tmp_const = model.linearNumExpr();
        		tmp_const.addTerm( 1.0, p_var[i]);
            	tmp_const.addTerm(-1.0, I_var[i1][n[i1]]);
            	model.addLe(tmp_const, 0.0, "pIconst" + varIndex + i1);
            }

        	// y_{ne,ns} <= I_{ne+1} and p_{ne,ns} <= I_{ns+1}
        	for(int i1=0; i1<=N; i1++){
        		IloLinearNumExpr tmp_const = model.linearNumExpr();
        		tmp_const.addTerm( 1.0, y_var[i]);
        		if(n[i1]<M_ub[i1]) tmp_const.addTerm( -1.0, I_var[i1][n[i1]+1]); 
            	model.addLe(tmp_const, 0.0, "yIconst" + varIndex + i1);
            }

        	// y_{ne,ns} <= p_{ne,ns} - (N+1 - I_{ne+1} - sum I_{ns+1}) 
        	{
        		IloLinearNumExpr tmp_const = model.linearNumExpr();
            	tmp_const.addTerm(-1.0, p_var[i]);
            	tmp_const.addTerm( 1.0, y_var[i]);
	        	for(int i1=0; i1<=N; i1++)
		        	if(n[i1]<M_ub[i1]) tmp_const.addTerm( -1.0, I_var[i1][n[i1]+1]);
            	model.addGe(tmp_const, -(N+1), "ypIconst" + varIndex);
        	}
    	}

    	for(int i=0; i<M_max; i++){
        	int[] n = getIndices(i);
        	String varIndex = "";  for(int j=0; j<=N; j++) varIndex += n[j];
			model.addEq(constraint[i], 0.0, "eqConst"+varIndex);
        }
    	

    	for(int i=0; i<=N; i++)
        	for(int j=0; j<=M_ub[i]; j++){
        		IloLinearNumExpr tmp_const = model.linearNumExpr();
        		tmp_const.addTerm(1.0, I_var[i][j]);
        		if(j<M_ub[i])
        			tmp_const.addTerm(-1.0, I_var[i][j+1]);
        		model.addGe(tmp_const, 0.0, "IIconst" + i + "" + j);
        	}
    	
    	model.addMinimize(cost);
    	model.addEq(norm_constraint, 1.0, "normConst");

    	//model.exportModel(workpath + "/model1.lp");
        	
	}
	
	public void exportModel() throws IloException{
		model.exportModel(workpath + "/model1.lp");
	}

	public void setIndicatorLimits(int[] ll, int[] ul) throws IloException{
    	this.ll = ll;
    	this.ul = ul;
		for(int i=0; i<=N; i++){
    		for(int j=0; j<=M_ub[i]; j++){
        		I_var[i][j].setLB(0);
        		I_var[i][j].setUB(1);
    		}
        	if(ll!=null)
	    		for(int j=0; j<=this.ll[i]; j++)
	        		I_var[i][j].setLB(1);
        	if(ul!=null)
	        	for(int j=this.ul[i]+1; j<=M_ub[i]; j++)
	        		I_var[i][j].setUB(0);
    	}
	}
	
	public void setStartSolution() throws IloException{
	     int maxi=0;
		 for(int i=0; i<=N; i++) maxi += (M_ub[i]+1);
		 IloNumVar[] startVar = new IloNumVar[maxi];
	     double[] startVal = new double[maxi];
	     for (int i = 0, idx = 0; i <= N; ++i)
	         for (int j = 0; j <= M_ub[i]; ++j) {
	             startVar[idx] = I_var[i][j];
	             startVal[idx] = 1;
	             ++idx;
	         }
	     model.addMIPStart(startVar, startVal);
	     startVar = null;
	     startVal = null;	     
	}

	public double optimize() throws IloException {
        //model.addMIP
    	if ( model.solve() ) {
        	if(optIsum==null) optIsum = new int[N+1];
    		for(int i=0; i<=N; i++){
            	optIsum[i] = 0;
        		for(int j=1; j<=M_ub[i]; j++) optIsum[i] += Math.round(model.getValue(I_var[i][j]));
        	}
    		return model.getObjValue();
    	}
		return -100000000; //if not solved
    }
    
	public boolean doIteration() throws IloException{
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
    
    public void printIndicators() throws UnknownObjectException, IloException{
    	for(int i=0; i<=N; i++){
        	for(int j=0; j<=M_ub[i]; j++)
        		System.out.print(model.getValue(I_var[i][j]) + "\t");
        	System.out.println();
    	}
    }
    public void printIndicatorSums() throws UnknownObjectException, IloException{
    	for(int i=0; i<=N; i++){
        	int sum = 0;
    		for(int j=1; j<=M_ub[i]; j++)
        		sum += Math.round(model.getValue(I_var[i][j]));
        	if(i==0) System.out.println("Number of Engineers: "  + sum);
        	else     System.out.println("Number of Parts "+i+":   " + sum);
    	}
    }
    public void printIndicatorShort(){
    	for(int i=0; i<=N; i++) System.out.print(optIsum[i] + "\t");
    	System.out.println();
    }

    public void printYvariables() throws UnknownObjectException, IloException{
    	for(int i=0; i<M_max; i++){
        	int[] n = getIndices(i);	
        	for(int j=0; j<=N; j++) System.out.print(n[j]);
        	System.out.println("\t" + model.getValue(y_var[i]));
    	}
    }
    public void printPvariables() throws UnknownObjectException, IloException{
    	for(int i=0; i<M_max; i++){
        	int[] n = getIndices(i);	
        	for(int j=0; j<=N; j++) System.out.print(n[j]);
        	System.out.println("\t" + model.getValue(p_var[i]));
    	}
    }

    public void printPvariables2D() throws UnknownObjectException, IloException{
    	if(N==1){
    		int[] n = new int[2];
    		for(int i=M_ub[0]; i>=0; i--){
    			for(int j=M_ub[1]; j>=0; j--){
    				n[0] = i; n[1] = j;
    				int i1 = getIndex(n);
    				System.out.print("\t" + model.getValue(p_var[i1]));
    			}
    			System.out.println();
    		}
        }
    }

    public void printMarginalProbabilities() throws UnknownObjectException, IloException{
    	for(int i=0; i<M_max; i++){
        	int[] n = getIndices(i);	
        	for(int j=0; j<=N; j++) 
        		p_mar[j][n[j]] += model.getValue(p_var[i]);
    	}
    	for(int j=0; j<=N; j++){
    		System.out.print("item " + j + ":");
    		for(int i=0; i<=M_ub[j]; i++)
    			System.out.print("\t" + p_mar[j][i]);
    		System.out.println();
    	}
    }

    public void printLossProbability() throws UnknownObjectException, IloException{
    	double p_loss=0;
    	for(int i=0; i<M_max; i++){
        	int[] n = getIndices(i);	
        	boolean loss = false;
        	for(int j=0; j<=N; j++) if(n[j]==optIsum[j]) loss = true;
        	if(loss)
        		p_loss += model.getValue(p_var[i]);
    	}
    	System.out.println("Loss probability: "  + p_loss);
    	
    }

    public void cleanupModel() throws IloException{
    	model.clearModel();
    	model.end();
    }

    public void setParameters() throws IloException{
    	setParametersFromFile(defaultParamFile);
    }
    
    public void setParametersFromFile(String fixedfile){
        if ( fixedfile != "" ) {
           try {
        	   model.readParam(fixedfile);
           } catch (IloException e) {
				// don't do anything
           }
        }
    }
    
    public void tuneModel(String fixedfile) throws IloException{
        if ( fixedfile != null ) {
           	model.tuneParam();
           	model.writeParam(fixedfile);
            System.out.println("Tuned parameters written to file '" +
            		fixedfile + "'");
         }
    }
	
    public void tuneModel() throws IloException{
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
	
	public void setLoggingOff(){
		if(model!=null){
			model.setOut(null);
			model.setWarning(null);
		}
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
}
