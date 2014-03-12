package OptimizationProblem;

import ilog.concert.IloColumn;
import ilog.concert.IloException;
import ilog.concert.IloNumVar;
import ilog.concert.IloNumVarType;
import ilog.concert.IloObjective;
import ilog.concert.IloRange;
import ilog.cplex.IloCplex;
import ilog.cplex.IloCplex.UnknownObjectException;

import java.util.Iterator;

public class ServiceEngineersOptimization {
	
	int N; //number of spare types
	int M; //state space truncation limit
	double lambda; //failure rate;
	double[] mu;   //service lead time service rates, mu[0] corresponds to engenders, the rest to spares
	double[] alpha; //probability that spare i is requested
	
	int M_max;
	
	IloNumVar[][] I_var;
	IloNumVar[] p_var, y_var;

	IloCplex model;
    
    public static String workpath = "/Users/andrei/Documents/Research/Service Engineers";


	public ServiceEngineersOptimization(int n, int m, double lambda,
			double[] mu, double[] alpha) {
		super();
		N = n;
		M = m;
		this.lambda = lambda;
		this.mu = mu;
		this.alpha = alpha;
		
		M_max = 1; for(int i=0; i<=N; i++) M_max *= (M+1);
		

	}

	
	public static void main(String[] args){
	
		
	
	
	};
	
	public void formLP() throws IloException{
	    IloRange[] constraint;
	    IloColumn[] p_column, y_column;
	    IloColumn[][] I_column;

	    model = new IloCplex();
		//model.setOut(null);
		//model.setWarning(null);
        constraint = new IloRange[M_max+1];
        p_column   = new IloColumn[M_max+1];
        y_column   = new IloColumn[M_max+1];
        
        p_var 	   = new IloNumVar[M_max+1];	
        y_var 	   = new IloNumVar[M_max+1];	

        I_column   = new IloColumn[N+1][M+1];	
        I_var      = new IloNumVar[N+1][M+1];	
        
        IloObjective cost = model.addMinimize();
			
    	for(int i=0; i<=N; i++)
        	for(int j=0; j<=M; j++){
        		int i1 = (j>0) ? 1 : 0;
        		I_column[i][j] = model.column(cost, i1*1.0);
        	}

    	IloRange norm_constraint = model.addEq(null, 1.0, "normConst");
    	for(int i=0; i<M_max; i++){
        	int[] n = getIndices(i);
        	String varIndex = "";  for(int j=0; j<=N; j++) varIndex += n[j];
        	
    		double cost_factor = 10;
        	p_column[i] = model.column(cost, cost_factor); //TODO update coefficient
        	
        	y_column[i] = model.column(cost, -cost_factor); //TODO update coefficient

        	//normalization constraint
        	p_column[i] = p_column[i].and( model.column(norm_constraint, 1.0));

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
            	if(n[i1]<M) I_column[i1][n[i1]+1] = I_column[i1][n[i1]+1].and(model.column(tmp_const, -1.0));
            }

        	// y_{ne,ns} <= p_{ne,ns} - (1 - I_{ne+1}) and y_{ne,ns} <= p_{ne,ns} - (1 - I_{ns+1})
        	{
        		IloRange tmp_const = model.addGe(null, -(N+1),  "ypIconst" + varIndex);
	        	p_column[i] = p_column[i].and(model.column(tmp_const, -1.0));
	        	y_column[i] = y_column[i].and(model.column(tmp_const, 1.0));
	        	for(int i1=0; i1<=N; i1++)
		        	if(n[i1]<M) I_column[i1][n[i1]+1] = I_column[i1][n[i1]+1].and(model.column(tmp_const, -1.0));
        	}

        	
        }
    	
    	for(int i=0; i<M_max; i++){
        	int[] n = getIndices(i);
        	String varIndex = "";  for(int j=0; j<=N; j++) varIndex += n[j];
        	p_var[i] = model.numVar(p_column[i], 0, 1, IloNumVarType.Float, "p" + varIndex);
        	y_var[i] = model.numVar(y_column[i], 0, 1, IloNumVarType.Float, "y" + varIndex);
        }

    	for(int i=0; i<=N; i++)
        	for(int j=0; j<=M; j++){
        		IloRange tmp_const = model.addGe(null, 0.0, "IIconst" + i+j);
        		I_column[i][j] = I_column[i][j].and(model.column(tmp_const, 1.0));
        		if(j<M)
        			I_column[i][j+1] = I_column[i][j+1].and(model.column(tmp_const, -1.0));

            	I_var[i][j] = model.numVar(I_column[i][j], 0, 1, IloNumVarType.Bool, "I"+i+j);
        	}

    	model.exportModel(workpath + "/model1.lp");
        	
	}
	
	public void exportModel() throws IloException{
		model.exportModel(workpath + "/model1.lp");
	}
	

	public void addIndicatorLimits(int[] n) throws IloException{
    	for(int i=0; i<=N; i++)
        	for(int j=n[i]+1; j<=M; j++)
        		I_var[i][j].setUB(0);
	}
	
    public double Optimize() throws IloException {
        if ( model.solve() ) return model.getObjValue();
		return -100000000; //if not solved
    }
    
    
    public void printIndicators() throws UnknownObjectException, IloException{
    	for(int i=0; i<=N; i++){
        	for(int j=0; j<=M; j++)
        		System.out.print(model.getValue(I_var[i][j]) + "\t");
        	System.out.println();
    	}
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
//    public double getDual(int rowIndex) throws IloException{
//    	return model.getDual(constraint[rowIndex]);
//    }
//
//    public double getSlack(int rowIndex) throws IloException{
//    	return model.getSlack(constraint[rowIndex]);
//    }
//    
    public void cleanupModel() throws IloException{
    	model.clearModel();
    	model.end();
    }
	
	
	protected int getIndex(int[] n){
		int k = 0;
		for(int i=0; i<=N; i++){
			if(n[i]>=0)
				k = k*(M+1) + n[i];
			else
				return -1;
		}
		
		return k;
	}
	
	protected int[] getIndices(int i){
		int[] tempk = new int[N+1];
		
		int k = i;
		for(int j=N; j>=0; j--){
			tempk[j] = k - (k/(M+1))*(M+1);
			k /= (M+1);
		}
		
		return tempk;
	}
	
	
}
