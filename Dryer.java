import java.io.File;
import java.util.*;
import java.io.FileNotFoundException;
import java.io.PrintWriter;

public class Dryer 
{
	/*variables that specific the index of different properties in the water sat table*/
	public static boolean did = false;
	public static final int T_indx = 0;
	public static final int Psat_indx = 1; 
	public static final int hliquid_indx = 2;
	public static final int hvapor_indx = 3;
	public static final int Cp_indx = 4;
	public static final int DynamicViscosity_indx = 5;
	public static final int density_indx = 6;
	public static final int K_indx = 7;
	
	//add prandtl, specific heat
	//jk air specific heat is .241 and barely changes                                 Cv
													//    T          Density       Specific Heat   Thermal Conductivity        Dynamic Viscosity   Prandtl Number
													//    F       1b/ft^2             Btu(lb*R)       Btu/(h*ft*R)               lb/fts                 NaN
	public static double [][] AirPropTable = {		 	 {70,		0.07489,           0.171,           .01457,                  .0000123,              .7306},
										 				 {80,		0.07350,            0.171,          .01481,                  .00001247,             .7290},
										 				 {90,       0.07217,            0.171,          .01505,                  .00001265,             .7275},
										 				 {100,      0.07088,             0.171,          0.01529,                .00001281,           .7260},
										 				 {110,      0.06963,            0.171,           0.01552,                .00001299,            .7245}};
	
	public static double [] [] AirSatTable={ //page 1047 A-21E of textbook
		/*              ENTHALPY       */
		/*T        h          Pr 		Cp			Dynamic Viscosity		Density */
		/*R    Btu/lbm       NaN  					lb/fts x 10^-6         lb/ft^3 */
		{520, 124.27, 1.2147},
		{537, 128.10, 1.3860},
		{540, 129.06, 1.3860},
		{560, 133.86, 1.5742},
		{580, 138.66, 1.78},
		{600, 143.47, 2.005},
		{620, 148.28, 2.249},
		{640, 153.09, 2.514},
		{660, 157.92, 2.801},
		{680, 162.73, 3.111},
		{700, 167.56, 3.446},
		{720, 172.39, 3.806},
		{740, 177.23, 4.193},
		{760, 182.08, 4.607},
		{780, 186.94, 5.051},
		{800, 191.81, 5.526},
		{820, 196.69, 6.033},
		{840, 201.56, 6.573},
		{860, 206.46, 7.149},
		{880, 211.35, 7.761},
		{900, 216.26, 8.411},
		{920, 221.18, 9.102},
		{940, 226.11, 9.834},
		{960, 231.06, 10.61},
		{980, 236.02, 11.43},
		{1000,240.98, 12.30}};
	
 //add prandtl, specific heat
	
	//public static final double Pr_air = 1.00;
	public static double [] [] WaterSatTable={
		/*              ENTHALPY       */
		/*T   Psat  h liquid  h vapor                  Dynamic Viscosity            Density       Thermal Conductivity*/
		/*F   psia  Btu/lbm    Btu/lbm        Cp             lb/fts x 10^-6         lb/ft^3        (Btu/hr-ft-F)*/
		{70, 0.3632, 38.09,      1092,        .45,     6.556/1000000,               .00115,          0.0106},
		{80, 0.5073, 48.09,      1096.4,      .451,    6.667/1000000,               .00158,            0.0108},
		{90, 0.6988, 58.07,      1100.7,      .453,    6.778/1000000,               .00214,            0.0110},
		{100,0.9503, 68.05,      1105.0,      .454,    6.889/1000000,               .00286,            0.0112},
		{110,1.2763, 78.02,      1109.3,      .456,    7.000/1000000,               .00377,            0.0115},
		{120,1.6945, 88.00,      1113.5,      .458,    7.111/1000000,               .00493,            0.0117},
		{130,2.225, 97.98,       1117.8,      .460,    7.222/1000000,               .00636,            0.0120},
		{140,2.892, 107.96,      1121.9,      .463,    7.333/1000000,               .00814,            0.0122},
		{150,3.722, 117.96,      1126.1,      .465,    7.472/1000000,               .0103,             0.0125},
		{160,4.745, 127.96,      1130.1,      .468,    7.583/1000000,               .0129,             0.0128},
		{170,5.996, 137.97,      1134.2,      .472,    7.722/1000000,               .0161,             0.0131},
		{180,7.515, 147.99,      1138.2,      .475,    7.833/1000000,               .0199,             0.0134},
		{190,9.343, 158.03,      1142.1,      .479,    7.972/1000000,               .0244,             0.0137},
		{200,11.529,168.07,      1145.9,      .483,    8.083/1000000,               .0297,             0.0141},
		{210,14.125,178.14,      1149.7,      .487,    8.222/1000000,               .0354,             0.0144},
		{212,14.698,180.16,      1150.5,      .488,    8.250/1000000,               .0373,             0.0145} };
	
	public static double F2R(double Temp){
		return 459.67 + Temp;
	}
	
	public static double btulbm_2_dimensionless(double val){
		return val * 778.76974716;
	}
	/*Given a temperature T, return the Psat and liquid and vapor enthalpies 
	 *whose values fall between respective values of T1 and T2 where T1<T<T2  */
	public double[] interpolate(double Temp, double[][] table)
	{
		int T1_indx = 0;
		int T2_indx = table.length-1;
		/*iterate row by row looking up given Temp. If we encounter a number larger than given temp, and the 
		preceding number was lower than given Temp we know to interpolate the two and return that result*/
		for(int row = 0; row < table.length; row++){
			
			//If the temp already has defined values in the table, return the row
			if (table[row][T_indx] == Temp){
				double [] temp = new double [table[row].length];
				for(int col = 0; col < table[row].length; col++){
					temp[col] = table[row][col];}
				return temp;
			}
			
			if (Temp > table[row][T_indx]){
				T1_indx = row;//new closest lower-bound definition
			}
			
			if (Temp < table[row][T_indx]){
				T2_indx = row; //encountered upper-bound definition
				break; //break out of the for statement, we interpolate with T1 and T2 now
			}
		}//end of for loop

		double [] temp = new double [table[0].length];
		for(int col = 0; col <table[0].length; col++){
			temp[col] = table[T1_indx][col] + (Temp - (table[T1_indx][T_indx]))
					* (table[T1_indx][col]-table[T2_indx][col])/((table[T1_indx][T_indx]-table[T2_indx][T_indx]));
		}
		return temp;
	}
	
	public double calc_Qin(double flow_air, double enthalpy_final_air, double absolute_humidity_s1, double enthalpy_final_vapor, double enthalpy_initial_air, double enthalpy_initial_vapor){
		return (flow_air * ((enthalpy_final_air + absolute_humidity_s1 * enthalpy_final_vapor)-(enthalpy_initial_air + absolute_humidity_s1 * enthalpy_initial_vapor)));
	}
	
	public double calc_flow_air(double flow_liquid, double enthalpy_vapor_s3, double enthalpy_liquid, double enthalpy_air_s2, double enthalpy_air_s3, double absolute_humidity_s1, double enthalpy_vapor_s2){
		return (flow_liquid * (enthalpy_vapor_s3 - enthalpy_liquid))/((enthalpy_air_s2 - enthalpy_air_s3) + absolute_humidity_s1 * (enthalpy_vapor_s2 - enthalpy_vapor_s3));
	}
	
	public double calc_LHS(double enthalpy_air_s2, double absolute_humidity_s2, double enthalpy_vapor_s2, double enthalpy_liquid, double Psat_s3, double P3, double relative_humidity_s3){
		return enthalpy_air_s2 + absolute_humidity_s2 * enthalpy_vapor_s2 + enthalpy_liquid * ((0.622 * relative_humidity_s3 * Psat_s3)/(P3 - (relative_humidity_s3 * Psat_s3) ));
	}
	
	public double calc_RHS(double enthalpy_air_s3, double enthalpy_vapor_s3, double Psat_s3, double P3, double relative_humidity_s3){
		return enthalpy_air_s3 + enthalpy_vapor_s3 * (0.622 * relative_humidity_s3 * Psat_s3)/(P3-relative_humidity_s3*Psat_s3);
	}
	
	public double calc_absolute_humidity(double partial_P, double P){
		return 0.622 * partial_P/(P-partial_P);
	}
	
	public double calc_absolute_humidity_s3(double absolute_humidity_s1, double flow_water, double flow_air){
		return absolute_humidity_s1 + (flow_water/flow_air);
	}
	public static void main(String [] args)throws FileNotFoundException
	{
		//These objects help eventually with the writing of the data to the excel spreadsheet, and the organization of the data
		PrintWriter writer = new PrintWriter(new File("/Users/mortum987789/Desktop//T3Speculations2016.csv"));
		//If you need to run the code yourself, change the filepath
		StringBuilder sb = new StringBuilder();
		/*
		sb.append("flow_liquid,"+"T1,"+"Psat_s1,"+"enthalpy_liquid_s1,"+"enthalpy_vapor_s1,"+"enthalpy_air_s1,"+"partial_pressure_H20_s1,"+"absolute_humidity_s1,"
		+"T2,"+"Psat_s2,"+"enthalpy_liquid_s2,"+"enthalpy_vapor_s2,"+"enthalpy_air_s2,"
				+"T3,"+"Psat_s3,"+"enthalpy_liquid_s3,"+"enthalpy_vapor_s3,"+"enthalpy_air_s3,"+"relative humidity," +"airflow,"+"Qin,"+"absolute humidity s3,"
		+"LHS," + "RHS");
		sb.append('\n');*/
		sb.append("T2" +","  + "T3" +"," + "relative humidity s3"+"," + "air flow"+"," + "qin" + "," + "absolute humidity s3");	
		sb.append('\n');
		sb.append("F,"+"F,"+"%," +"lb/s,"+"Btu/s,"+"%,");
		sb.append('\n');
		Dryer us = new Dryer();
		/* 	GIVENS	*/
		final double P1 = 14.7; //pressure in lb/in^2 @ sealevel
		final double relative_humidity_s1 = 0.50;
		final double flow_liquid = 5.0000000000 / (60 * 45); //lb per second  .45359 lbs to kg
		
		final double T1 = 70; //Room Temperature
			double[] temp = us.interpolate(T1, Dryer.WaterSatTable);
			double Psat_s1 = temp[Psat_indx];
			double enthalpy_liquid_s1 = temp[hliquid_indx];
			double enthalpy_vapor_s1 = temp[hvapor_indx];
			temp = us.interpolate(F2R(T1), Dryer.AirSatTable);
			double enthalpy_air_s1 = temp[1];
		double partial_pressure_H20_s1 = relative_humidity_s1 * Psat_s1;
		double absolute_humidity_s1 = us.calc_absolute_humidity(partial_pressure_H20_s1, P1);
		
		
		for(double T2 = 190; T2 <= 212; T2 +=1) {
		//final double T2 = 212; //Max dryer temp is 212F
			temp = us.interpolate(T2, Dryer.WaterSatTable);
			double Psat_s2 = temp[Psat_indx];
			double enthalpy_liquid_s2 = temp[hliquid_indx];
			double enthalpy_vapor_s2 = temp[hvapor_indx];
			temp = us.interpolate(F2R(T2), Dryer.AirSatTable);
			double enthalpy_air_s2 = temp[1];
		
			
			
		/*   EQ Flow   */
		/* Find exhaust temp T3 by making a guess between T1 & T2 and seeing if LHS and RHS equations match*/
		for(double temporaryTemperature = T1; temporaryTemperature <= T2; temporaryTemperature ++){
			temp = us.interpolate(temporaryTemperature, Dryer.WaterSatTable);
			double Psat_s3 = temp[Psat_indx];
			double enthalpy_liquid_s3 = temp[hliquid_indx];
			double enthalpy_vapor_s3 = temp[hvapor_indx];
			temp = us.interpolate(F2R(temporaryTemperature), Dryer.AirSatTable);
			double enthalpy_air_s3 = temp[1];
			for (double relative_humidity_s3 = 0.7; relative_humidity_s3 <= 0.90; relative_humidity_s3+= 0.01){
			//enthalpy of liquid, hl, is 38.09
			//P3 is comparable to P1, since water loaded clothes are cool even though the air is considerably hotter, cause there's also all this 
			//air being blown around [and recycled, you'll worry about that later Mr. Program]
			//Also relative humidity for stage3 can be assumed to be anywhere between 70-90
			double LHS = us.calc_LHS(enthalpy_air_s2, absolute_humidity_s1, enthalpy_vapor_s2, enthalpy_liquid_s1, Psat_s3, P1, relative_humidity_s3);
			double RHS = us.calc_RHS(enthalpy_air_s3, enthalpy_vapor_s3, Psat_s3, P1, relative_humidity_s3);
			
			if(Math.abs(LHS - RHS) <= 0.5){
				//     double calc_flow_air(double flow_liquid, double enthalpy_vapor_s3, double enthalpy_liquid, double enthalpy_air_s2, double enthalpy_air_s3, double absolute_humidity_s1, double enthalpy_vapor_s2){
				double air_flow = us.calc_flow_air(flow_liquid, enthalpy_vapor_s3, enthalpy_liquid_s1, enthalpy_air_s2, enthalpy_air_s3, absolute_humidity_s1, enthalpy_vapor_s2);
				double qin = us.calc_Qin(air_flow, enthalpy_air_s2, absolute_humidity_s1, enthalpy_vapor_s2, enthalpy_air_s1, enthalpy_vapor_s1);
				double absolute_humidity_s3 = us.calc_absolute_humidity_s3(absolute_humidity_s1, flow_liquid, air_flow);
				sb.append(T2+"," + temporaryTemperature+","  + relative_humidity_s3+"," + air_flow+","+ qin+"," + absolute_humidity_s3);
				sb.append('\n');
				/*if(temporaryTemperature == 103.0 && T2 == 212.0 && relative_humidity_s3 == 0.71){
					
				sb.append(flow_liquid+ "," +T1+ "," +Psat_s1+ "," +enthalpy_liquid_s1+ "," +enthalpy_vapor_s1+ "," +enthalpy_air_s1+ "," +partial_pressure_H20_s1+ "," +absolute_humidity_s1+ "," +
						+T2+ "," +Psat_s2+ "," +enthalpy_liquid_s2+ "," +enthalpy_vapor_s2+ "," +enthalpy_air_s2+ "," +
						+temporaryTemperature+ "," +Psat_s3+ "," +enthalpy_liquid_s3+ ","+enthalpy_vapor_s3+ ","+enthalpy_air_s3+ "," + relative_humidity_s3+ "," +air_flow+ "," + qin+ "," +absolute_humidity_s3+ "," +
				+LHS+ "," + RHS);
						sb.append('\n');
				
					sb.append(T2 +","  + temporaryTemperature +"," + relative_humidity_s3+"," + air_flow+"," + qin + "," + absolute_humidity_s3);		
					sb.append('\n');
				}*/
				double flow_vapor_3 = air_flow * absolute_humidity_s3;
				
				/*if( !did){
				System.out.println("mu wat 103" + us.interpolate(103, WaterSatTable)[5]);
				System.out.println("mu wat 103" + us.interpolate(70, WaterSatTable)[5]);
				System.out.println("mu air 103" + us.interpolate(103, AirPropTable)[4]);
				System.out.println("mu air 70" + us.interpolate(70, AirPropTable)[4]);
				did = true;}*/
				
				}//end of RHS=LHS if
			}//end of forloop through relative humidity guesses
			}//end of forloop through T3 guesses
		}//end of forloop through T2 [max drying temperature] experimental values
		 writer.write(sb.toString());
	     writer.close();
	     System.out.println("done!");
	     
	     /*Still remaining: 
	      * Pressure drop across lint filter
	      * T_final for the inlet air (I.E. the temp that needs to be heated to 'T2', which is dryer max, by the heating element
	      * T_final for the exhaust air (after having exchanged heat with fresh inlet air)
	      * Length of the heat exchanger
	      */
	     
	     /*
	      * Methodology: plug and chug guesstimate numbers for Tf_exhaust & Tf_inlet to solve for a heat exchange
	      */
	     
	     //there are literally too many pipes to think about all of the
	     //restricting down to some options related to 89895K777 from previous project
	     //to look for trends / prospect of heat exchanger
	     
	     int ID_indx = 0;
	     int OD_indx = 1;
	     double [] [] DimList = {
	    		 //IDp   ODp
	 //   		 {1.084, 1.25}, 
	 //   		 {1.277, 1.375},
	 //   		 {1.402, 1.5},
	 //   		 {1.62, 1.75},
	  //  		 {1.76, 2},
	  //  		 {1.995, 2.125},
	 //   		 {2.245, 2.375},
	 //   		 {2.37, 2.5},
	 //   		 {2.62, 2.75},
	  /*  		 {2.87, 3}*/
	    		 {6, 6.125},
	    		 {8, 8.125}}; //those silver bendy tubes you see sometimes
	     
	     double[] poop = us.interpolate(103, Dryer.WaterSatTable);
	     System.out.println(poop[Psat_indx]);
	     for(int IDp_indx = 0; IDp_indx < (DimList.length -1); IDp_indx ++){ //iterate through all but bottom row for potential inner pipe
	    	 int IDa_indx = IDp_indx + 1;
	    	 if (Math.abs((DimList[IDp_indx][OD_indx] - DimList[IDa_indx][ID_indx])) <= 0.1) //make sure the inside of outer pipe is resonably bigger than outside of inner pipe
	    	 {
	    		 IDa_indx += 1; //outer pipe goes a size up
	    	 }
	    	 for(int i = IDa_indx; i < DimList.length; i++){
	    		 	System.out.println("It's now a party. FUCKIN PARTY");
	    			 double IDa = DimList[i][ID_indx]/12;
	    			 double IDp = DimList[IDp_indx][ID_indx]/12;
	    			 double ODp = DimList[IDp_indx][OD_indx]/12;
	    			 //For debugging to know we're now in the correct piping configs
	    			 System.out.println();
	    			 System.out.println("IDp = " + IDp + "   ODp = " + ODp + "  IDa = " + IDa);
	    			 double Dh = IDa - ODp;
	    			 double De = ((IDa * IDa) - (ODp * ODp))/ODp;
	    			 System.out.println();
	    			 System.out.println("Dh = " + Dh + "   De = " + De);
	    			 
	    			 double Ap = Math.PI * (IDp/2)* (IDp/2);
	    			 double Aa = (Math.PI * (IDa/2)* (IDa/2)) - Ap;
	    			 
	    			 //the results of T3 = 103, relative_humidity_3 = .71
	    			 double air_flow = 0.074451013;
	    			 double absolute_humidity_s3 = 0.0327;
	    			 
	    			 double flow_vapor_3 = air_flow * absolute_humidity_s3;
	    			 double flow_total_exhaust = flow_vapor_3 + air_flow;
	    			 
	    			 double AreaFilter = (13*16)/144; //ft^2
	    			                                     //T = 103 in Rankine is 562.67
	    			 double density_exhaust = (((14.8989)/(562.67 * 1716))*4633.056)*(1 + .71)/(1+1.609*.71);
	    			 double velocity_s3 = flow_total_exhaust/(AreaFilter * density_exhaust);
	    			 double velocity_s4 = flow_total_exhaust/(.85*(AreaFilter * density_exhaust));
	    			 
	    			 System.out.println();
	    			 System.out.println("density_exhaust = " + density_exhaust + "   velocity_s3 = " + velocity_s3 + "  velocity_s4 = " + velocity_s4 + "  flow_vapor_3 = " + flow_vapor_3+ "  flow_total_exhaust = " + flow_total_exhaust);
	    			 
	    			 double Pdrop = ((velocity_s3 * velocity_s3) - (velocity_s4 * velocity_s4)) * (density_exhaust)/2;
	    			 
	    			 double density_inlet = (4633.056 * (14.8989) / (1716 * 562.67))*(1.5)/(1+1.609*.5);
	    			 
	    			 double flow_vapor_1 = air_flow * absolute_humidity_s1; 
	    			 double flow_total_inlet = air_flow + flow_vapor_1;
	    			 
	    			 System.out.println();
	    			 System.out.println("Pdrop = " + Pdrop + "   density_inlet = " + density_inlet + "  flow_vapor_1 = " + flow_vapor_1+ "  flow_total_inlet = " + flow_total_inlet);
	    			 
	    			 
	    			 double velocity_inlet = flow_total_inlet/(Aa * density_inlet);
	    			 double velocity_exhaust = flow_total_exhaust/(Ap * density_exhaust);
	    			 
	    			 //get thermal conductivity of water
	    			 double thermal_conductivity_water = 0.364;
	    			 //now assume that both inlet and outlet are just air in terms of thermal conductivity
	    			 
	    			 System.out.println();
	    			 System.out.println("velocity_inlet = " + velocity_inlet + "   velocity_exhaust = " + velocity_exhaust + "  thermal_conductivity_water = " + thermal_conductivity_water);
	    			 
	    			 
	    			 double T_avg = (103+70)/2;
	    			 double PHI_AV = Math.pow((1 + Math.pow(us.interpolate(T_avg, AirPropTable)[4]/us.interpolate(T_avg, WaterSatTable)[5], 0.5) * 0.887603),2)/4.57043;
	    			 double PHI_VA = Math.pow((1 + Math.pow(us.interpolate(T_avg, WaterSatTable)[5]/us.interpolate(T_avg, AirPropTable)[4], 0.5) * 1.12663),2)/3.60076;
	    			 

	    			 System.out.println();
	    			 System.out.println("T_avg = " + T_avg + "   PHI_AV = " + PHI_AV + "  PHI_VA = " + PHI_VA);
	    			 
	    			 double Dynamic_Viscosity_Outlet  = us.interpolate(T_avg, AirPropTable)[4]/(1+PHI_AV*1.61*absolute_humidity_s3) + us.interpolate(T_avg, WaterSatTable)[5]/(1+PHI_VA/(1.61*absolute_humidity_s3));
	    			 double Dynamic_Viscosity_Inlet  = us.interpolate(T_avg, AirPropTable)[4]/(1+PHI_AV*1.61*absolute_humidity_s1) + us.interpolate(T_avg, WaterSatTable)[5]/(1+PHI_VA/(1.61*absolute_humidity_s1));
	    			 
	    			 System.out.println();
	    			 System.out.println("Dynamic_Viscosity_Outlet = " + Dynamic_Viscosity_Outlet + "   Dynamic_Viscosity_Inlet = " + Dynamic_Viscosity_Inlet );
	    			 
	    			 double KinematicViscosity_Inlet = Dynamic_Viscosity_Inlet/density_inlet;
	    			 double KinematicViscosity_Outlet = Dynamic_Viscosity_Outlet/density_exhaust;
	    			 

	    			 System.out.println();
	    			 System.out.println("KinematicViscosity_Inlet = " + KinematicViscosity_Inlet + "   KinematicViscosity_Outlet = " + KinematicViscosity_Outlet );
	    			 
	    			 double Re_inlet = velocity_inlet * De/KinematicViscosity_Inlet;
	    			 double Re_outlet = velocity_exhaust * IDp/KinematicViscosity_Outlet;
	    			 

	    			 System.out.println();
	    			 System.out.println("Re_inlet = " + Re_inlet + "   Re_outlet = " + Re_outlet );
	    			 
	    			 double specific_heat_air_inlet = us.interpolate(70, AirPropTable)[2];
	    			 double specific_heat_air_outlet = us.interpolate(103, AirPropTable)[2];
	    			 
	    			 double specific_heat_water_inlet = us.interpolate(70, WaterSatTable)[4];
	    			 double specific_heat_water_outlet = us.interpolate(103, WaterSatTable)[4];
	    			 
	    			 double exhaust_cp = (specific_heat_air_outlet+absolute_humidity_s3*specific_heat_water_outlet)/(1+absolute_humidity_s3);
	    			 double inlet_cp = (specific_heat_air_inlet+absolute_humidity_s1*specific_heat_water_inlet)/(1+absolute_humidity_s1);
	    			 System.out.println();
	    			 System.out.println("exhaust_cp = " + exhaust_cp + "   inlet_cp = " + inlet_cp);
	    			 //exhaust [hot] is 0.4
	    			 //inlet [cold] is 0.3
	    			 
	    			 double thermal_conductivity_air_inlet = us.interpolate(70, AirPropTable)[3];
	    			 double thermal_conductivity_air_outlet = us.interpolate(103, AirPropTable)[3];
	    			 
	    			 double thermal_conductivity_water_inlet = us.interpolate(70, WaterSatTable)[7];
	    			 double thermal_conductivity_water_outlet = us.interpolate(103, WaterSatTable)[7];
	    			 
	    			 double exhaust_k = (thermal_conductivity_air_outlet+absolute_humidity_s3*thermal_conductivity_water_outlet)/(1+absolute_humidity_s3);
	    			 double inlet_k = (thermal_conductivity_air_inlet+absolute_humidity_s1*thermal_conductivity_water_inlet)/(1+absolute_humidity_s1);
	    			 System.out.println();
	    			 System.out.println("inlet_k = " + inlet_k + "   exhaust_k = " + exhaust_k);
	    			 
	    			 double exhaust_pr = 3600*(exhaust_cp * Dynamic_Viscosity_Outlet)/exhaust_k; //the 3600 is to cancel the time units
	    			 double inlet_pr = 3600*(inlet_cp * Dynamic_Viscosity_Inlet)/inlet_k;
	    			 
	    			 System.out.println();
	    			 System.out.println("inlet_pr = " + inlet_pr + "   exhaust_pr = " + exhaust_pr);
	    			
	    			 //0.023Re^(4/5)Pr^N where N = 0.4 for exhaust and 0.3 for inlet
	    			 
	    			 double inlet_Nu = 0.023 * Math.pow(Re_inlet, (4.0/5.0)) * Math.pow(inlet_pr, 0.3);
	    			 double exhaust_Nu = 0.023 * Math.pow(Re_outlet, (4.0/5.0)) * Math.pow(exhaust_pr, 0.4);
	    			 
	    			 System.out.println();
	    			 System.out.println("inlet_Nu = " + inlet_Nu + "   exhaust_Nu = " + exhaust_Nu);
	    			 
	    			 double inlet_h = inlet_k * inlet_Nu / De;
	    			 double exhaust_h = exhaust_k * exhaust_Nu / ODp;
	    			 double U = Math.pow((1/inlet_h + 1/exhaust_h), -1);
	    			 //Rair = .00227
	    			 double FouledU = Math.pow(1/U + 2 * .00227, -1);
	    			 
	    			 System.out.println();
	    			 System.out.println("inlet_h = " + inlet_h + "   exhaust_h = " + exhaust_h + "  U = " + U + "  FouledU = " + FouledU);
	    			 
	    			 double Cmax = exhaust_cp * flow_total_exhaust; double Cmin = inlet_cp * flow_total_inlet;
	    			 if (Cmax < Cmin){
	    				 double temporaryStationary = Cmin; 
	    				 Cmin = Cmax;
	    				 Cmax = temporaryStationary;
	    			 }
	    			 
	    			 System.out.println();
	    			 System.out.println("Cmin = " + Cmin + "   Cmax = " + Cmax);
	    			 
	    			 
	    			 double qmax = Cmin *(103-70);
	    			 double NTU = U * Ap / Cmin;
	    			 double effectiveness =  (1 - Math.exp(-NTU*(1 - Cmin/Cmax))) / ((1 - Cmin/Cmax) * Math.exp(-NTU*(1 - Cmin/Cmax)));
	    			 
	    			 System.out.println();
	    			 System.out.println("qmax = " + qmax + "   NTU = " + NTU + "  effectiveness = " + effectiveness);
	    			 
	    			 double q = effectiveness * qmax;
	    			 double exhaust_final = -q/(exhaust_cp * flow_total_exhaust) + 103;
	    			 double inlet_final = (q/(inlet_cp * flow_total_inlet) - 70)*-1;
	    			 
	    			 System.out.println();
	    			 System.out.println("q = " + q + "   exhaust_final = " + exhaust_final + "  inlet_final = " + inlet_final);
	    			 
	    			 
	    			 
	    			 
	    			 /*
	    			  * Frustrated, lets try solving from the other side since we've arbitrarily picked a length and pipe sizes
	    			  * L=Ao/(pi*ODeg)   
	    			  */
	    			 double L = 5; //ft          //converting ODp to feet
	    			 double Ao = L * Math.PI * ODp/12;  
	    			 double LMTD = Ao * q / FouledU;
	    			 
	    			 System.out.println();
	    			 System.out.println("Ao = " + Ao + "   LMTD = " + LMTD);
	    			 
	    			 
	    		 }
	    	 }
	     
			
//			}//end of RHS=LHS if
//		}//end of forloop through relative humidity guesses
//		}//end of forloop through T3 guesses
//	}//end of forloop through T2 [max drying temperature] experimental values
	     
	}//main statement end
}//end of dryer class
