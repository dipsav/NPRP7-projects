package OptimizationProblem;

import ilog.concert.IloColumn;
import ilog.concert.IloException;
import ilog.concert.IloNumVarType;
import ilog.concert.IloObjective;
import ilog.concert.IloRange;
import ilog.cplex.IloCplex;

import java.util.Iterator;

public class ServiceEngineersOptimization {
	
	int N; //number of spare types
	int M; //state space truncation limit
	double lambda; //failure rate;
	double[] mu;   //service lead time service rates, mu[0] corresponds to engenders, the rest to spares
	double[] alpha; //probability that spare i is requested
	
	int M_max;
	
    IloCplex model;
    


	public ServiceEngineersOptimization(int n, int m, double lambda,
			double[] mu, double[] alpha) {
		super();
		N = n;
		M = m;
		this.lambda = lambda;
		this.mu = mu;
		this.alpha = alpha;
		
		M_max = 1; for(int i=0; i<=N; i++) M_max *= M;

	}

	
	public static void main(String[] args){
	
		
	
	
	};
	
	public void formLP() throws IloException{
	    IloRange[] constraint;
	    IloColumn[] p_column, y_column;
	    IloColumn[][] I_column;

	    model = new IloCplex();
		model.setOut(null);
		model.setWarning(null);
        constraint = new IloRange[M_max];
        p_column   = new IloColumn[M_max];
        y_column   = new IloColumn[M_max];
        I_column   = new IloColumn[N+1][M+1];	
        
        IloObjective cost = model.addMinimize();
			
    	for(int i=0; i<=N; i++)
        	for(int j=0; j<=M; j++){
        		I_column[i][j] = model.column(cost, 0.0);
        		model.numVar(I_column[i][j], 0, 1, IloNumVarType.Int);
        	}

    	IloRange norm_constraint = model.addEq(null, 1.0);
    	for(int i=0; i<=M_max; i++){
        	p_column[i] = model.column(cost, 1.0); //TODO update coefficient
        	model.numVar(p_column[i], 0, 1, IloNumVarType.Float);
        	
        	y_column[i] = model.column(cost, -1.0); //TODO update coefficient
        	model.numVar(y_column[i], 0, 1, IloNumVarType.Float);

        	//normalization constraint
        	p_column[i] = p_column[i].and( model.column(norm_constraint, 1.0));

        	//equilibrium constraints
        	constraint[i] = model.addEq(null, 0.0);
        	int[] n = getIndecies(i);
        	
        	y_column[i] = y_column[i].and(model.column(constraint[i], lambda));
        	for(int j=0; j<=N; j++){
            	p_column[i] = p_column[i].and(model.column(constraint[i], n[j]*mu[j]));
        		
        		n[j]--;
            	int i1 = getIndex(n);
            	if(i1>=0)
            		p_column[i]  = p_column[i].and( model.column(constraint[i1], (n[j]+2)*mu[j]));
            	n[j]++;
        	}
        	n[0]--;
        	for(int j=1; j<=N; j++){
        		n[j]--;
            	int i1 = getIndex(n);
            	if(i1>=0)
            		y_column[i1] = y_column[i1].and(model.column(constraint[i], -lambda*alpha[j]));
            	n[j]++;
        	}
        	n[0]++;        	

        	
        	// p_{ne,ns} <= I_ne and p_{ne,ns} <= I_ns
        	IloRange tmp_const = model.addLe(null, 0.0);
        	p_column[i] = p_column[i].and(model.column(tmp_const, 1.0));
        	for(int i1=0; i1<=N; i1++)
            	for(int j=0; j<=M; j++){
            		I_column[i1][j] = I_column[i1][j].and(model.column(tmp_const, -1.0));
            	}

        	// y_{ne,ns} <= p_{ne,ns} - (1 - I_ne) and y_{ne,ns} <= p_{ne,ns} - (1 - I_ns)
        	tmp_const = model.addLe(null, -1.0);
        	p_column[i] = p_column[i].and(model.column(tmp_const, -1.0));
        	y_column[i] = y_column[i].and(model.column(tmp_const, 1.0));
        	for(int i1=0; i1<=N; i1++)
            	for(int j=0; j<=M; j++){
            		I_column[i1][j] = I_column[i1][j].and(model.column(tmp_const, -1.0));
            	}
        	// I_{ne+1} <= I_ne and I_{ns+1} <= I_ns
        	for(int i1=0; i1<=N; i1++)
            	for(int j=0; j<=M; j++){
                	tmp_const = model.addGe(null, 0.0);
            		I_column[i1][j] = I_column[i1][j].and(model.column(tmp_const, 1.0));
            		if(j<M)
            			I_column[i1][j+1] = I_column[i1][j+1].and(model.column(tmp_const, -1.0));
            	}
        }
    	
    	model.exportModel("model1.lp");
        	
	}
	
	
    static public int weightIndex(){return 0;}
    static public int volumeIndex(){return 1;}

    
    public double Optimize() throws IloException {
        if ( model.solve() ) return model.getObjValue();
		return -100000000; //if not solved
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
				k = k*M + n[i];
			else
				return -1;
		}
		
		return k;
	}
	
	protected int[] getIndecies(int i){
		int[] tempk = new int[N+1];
		
		int k = i;
		for(int j=N; j>=0; j--){
			tempk[j] = k - (k/M)*M;
			k /= M;
		}
		
		return tempk;
	}
	
	
}
