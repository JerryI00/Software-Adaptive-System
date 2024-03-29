package org.femosaa.core;

public class EAConfigure {
	
	static {
		//con = new EAConfigure(100, 10, 0.9, 0.1);
		//con = new EAConfigure(20, 50, 0.5, 0.05);
		con = new EAConfigure(100, 50, 0.9,  0.1);
		//con = new EAConfigure(100, 300, 0.8, 0.02);
	} 

	public int pop_size;
	public int generation;
	public double crossover_rate;
	public double mutation_rate;
	
	public int measurement = 400;//50;//400;//1000; //-
	
	
	public EAConfigure(int pop_size, int generation, double crossover_rate,
			double mutation_rate) {
		super();
		this.pop_size = pop_size;
		this.generation = generation;
		this.crossover_rate = crossover_rate;
		this.mutation_rate = mutation_rate;
	}


	private static EAConfigure con;
	
	
	public static EAConfigure getInstance(){
		return con;
	}
	
	public void setupWSConfiguration(double crossover){
		//con = new EAConfigure(100, 50, crossover,  0.1);
		//seeding 100 AS - 
		con = new EAConfigure(100, 300, crossover, 0.02);
	}
	
	public void setupWSConfiguration(){
		//con = new EAConfigure(100, 100, 0.9,  0.1);
		//con = new EAConfigure(100, 300, 0.9,  0.1);
		//seeding 100 AS - 
		con = new EAConfigure(100, 500, 0.8, 0.02);
	}
	
	public void setupFLASHConfiguration(){
		con = new EAConfigure(10/*4,10*/, 500*10, 0.9,  0.1);
		//con = new EAConfigure(10/*4,10*/, 100, 0.9,  0.1);// for wodup
		//con = new EAConfigure(30, 20, 0.9,  0.1);
	}
	
	public void setupFLASHConfigurationVariant(int i){
		int pop = 100;
		int gen = 300;
		if (i == 0) {
			con = new EAConfigure(pop, gen, 0.9,  0.1);	
		} else if (i == 1) {
			con = new EAConfigure(pop, gen, 0.9,  0.01);	
		} else if (i == 2) {
			con = new EAConfigure(pop, gen, 0.8,  0.1);	
		} else if (i == 3) {
			con = new EAConfigure(pop, gen, 0.8,  0.01);	
		} else if (i == 4) {
			con = new EAConfigure((int)pop/2, 2 * (int)gen, 0.9,  0.1);	
		} else if (i == 5) {
			con = new EAConfigure((int)pop/2, 2 * (int)gen, 0.9,  0.01);	
		} else if (i == 6) {
			con = new EAConfigure((int)pop/2, 2 * (int)gen, 0.8,  0.1);	
		} else if (i == 7) {
			con = new EAConfigure((int)pop/2, 2 *  (int)gen, 0.8,  0.01);	
		}
		
	}
	
	public void setupNRPConfigurationVariant(int i){
		int pop = 100;
		int gen = 200;
		if (i == 0) {
			con = new EAConfigure(pop, gen, 0.8,  0.01);	
		} else if (i == 1) {
			con = new EAConfigure(pop, gen, 0.9,  0.01);	
		} else if (i == 2) {
			con = new EAConfigure(pop, gen, 0.8,  0.1);	
		} else if (i == 3) {
			con = new EAConfigure(pop, gen, 0.9,  0.1);	
		} else if (i == 4) {
			con = new EAConfigure((int)pop/2, 2 * (int)gen, 0.9,  0.1);	
		} else if (i == 5) {
			con = new EAConfigure((int)pop/2, 2 * (int)gen, 0.9,  0.01);	
		} else if (i == 6) {
			con = new EAConfigure((int)pop/2, 2 * (int)gen, 0.8,  0.1);	
		} else if (i == 7) {
			con = new EAConfigure((int)pop/2, 2 *  (int)gen, 0.8,  0.01);	
		}
		
	}
	
	public void setupFLASHConfigurationInternal(){
		con = new EAConfigure(50, 20, 0.9,  0.1);
		//con = new EAConfigure(60, 500*10, 0.9,  0.2);
	}
	
	public void setupFLASHConfiguration(int pop){
		con = new EAConfigure(pop, 500*10, 0.9,  0.1);
	}
	
	public void setupNRPConfiguration(){
		con = new EAConfigure(100, 200, 0.8,  0.01);//0.1
		//con = new EAConfigure(100, 200, 0.8,  0.01);//0.1
	}
	
	public void setupRUBiSSimConfiguration(){
		con = new EAConfigure(100, 10, 0.9, 0.1);
	}
	
	public void setupWSConfigurationOnlyMutation(){
		//con = new EAConfigure(100, 50, 0.0,  0.1);
		con = new EAConfigure(100, 300, 0.0, 0.02);
	}
	
	public void setupNRPConfigurationOnRun(int gen){
		con = new EAConfigure(100, gen, 0.8,  0.01);//0.1
		//con = new EAConfigure(100, 200, 0.8,  0.01);//0.1
	}
	
	public void setupWSConfigurationOnRun(int gen){
		//con = new EAConfigure(100, gen, 0.9,  0.1);//0.1
		con = new EAConfigure(100, gen, 0.8,  0.02);//0.1
	}

}
