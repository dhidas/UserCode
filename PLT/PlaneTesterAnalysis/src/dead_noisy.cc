#include "dead_noisy.h"


void Find_Dead_And_Noisy(TH2F * map_pixels, TH2F & map_maskedOff, TH2F & map_dead){
		//map_pixels is the main input, the other two can just be pointers, uninitialized
	
	/*
	 call with:
	 int nX = map_pixels->GetNbinsX();
	 int nY = map_pixels->GetNbinsY();
	 TH2F * map_maskedOff = new TH2F("map_maskedOff", "Noisy pixels to mask off", nX, 1, nX, nY, 1, nY);
	 TH2F * map_dead = new TH2F("map_dead", "dead and under-performing pixels", nX, 1, nX, nY, 1, nY);
	 Find_Dead_And_Noisy(map_pixels, *map_maskedOff, *map_dead);
	 */
	
	int ir_max = map_pixels->GetNbinsX();
	int ic_max = map_pixels->GetNbinsY();//these should be the bounds of the row and collumn for map_pixel
	
	TH2F * map_warned = new TH2F("map_warned", "warned pixels", ir_max, 1, ir_max, ic_max, 1, ic_max );
	
	Fill_loc(map_maskedOff,0);//0 is false = alive, 1 is true = noisy and masked off
	Fill_loc(map_dead,0);//0 is false = alive, 1 is true = noisy and masked off
	Fill_loc(*map_warned,0);//0 is false = alive, 1 is true = noisy and masked off
							//map_maskedOff and map_dead are already references, while map_warned is a pointer, needs the star.
	
	
		//#define MAX_NEIGHBORS 25   // size of neighbors[]
	/*
	 The warner looks at whether each pixel looks anomolously high or low compared to its neighborhood. If so, it gives it a warning, symbolized by a 1 in map_warned. 
	 This lets later averages be computed without skew from anomolous pixels. 
	 */
	
		//BEGIN WARNER     WWWW
		//mark pixels that look like they might be noisy with a warning. 
    for(int ir = 1; ir<=ir_max; ir++){
		for(int ic = 1; ic<=ic_max; ic++){//real
										  //	 if(map_dead.(ir, ic) ==1 ||
										  //	     map_maskedOff.GetBinContent(ir,ic)==1){
										  //	    	    map_warned->SetBinContent(ir, ic, 1);
										  //		    continue;
										  //	  }//endif, if already marked as bad, mark as warned
			
			float pixel = map_pixels->GetBinContent(ir, ic);
			float neighbors[]={-1,-1,-1,-1,-1,-1,
				-1,-1,-1,-1,-1,-1,
				-1,-1,-1,-1,-1,-1,
				-1,-1,-1,-1,-1,-1, -1};//25 -1's, only 24 may be filled with other things.
			
			int N_neighbors = get_neighbors(map_pixels,ir,ic, map_maskedOff, map_dead, neighbors);
			if(N_neighbors > 6){
				
				float local_mean = get_mean(neighbors);
				float local_stdev = get_stdev(neighbors, local_mean);
				float local_efficiency = get_efficiency(pixel, local_mean);
				
				if(   is_noisy(pixel,local_mean, local_stdev, local_efficiency)  )   map_warned->SetBinContent(ir, ic, 1);
				if(   is_dead(pixel,local_mean, local_stdev, local_efficiency)  ) map_warned->SetBinContent(ir, ic, 1);
			}//end if N_neighbors>6, if a pixel has fewer than 6 valid neighbors, rely 
			
		}//end for ic
    }//end for ir
	 //END WARNER
	
	
	
	
		//BEGIN GLOBAL FILTERING
		//This section uses global measures to declare a pixel dead or noisy
	
		//find global sum, not couting pixels that look anomolous. 
	float mean_sum=0;
	int n=0;
	for(int ir = 1; ir<ir_max; ir++){
		for(int ic = 1; ic<ic_max; ic++){//these are the real ones
			if(map_warned->GetBinContent(ir,ic) != 1){//xxxx
				mean_sum += map_pixels->GetBinContent(ir,ic);
				n++;
			}//endif
		}//end for ic
	}//end for ir
    float global_mean = mean_sum/n;
	
		//find global sample st.dev, not couting pixels that look anomolous
    n = 0;
	float stdev_sum=0;
    for(int ir = 1; ir<ir_max; ir++){
		for(int ic = 1; ic<ic_max; ic++){
			if(map_warned->GetBinContent(ir,ic) != 1){
				stdev_sum+=pow((map_pixels->GetBinContent(ir,ic) - global_mean),2);
				n++;
			}//endif
		}//end for ic
	}//end for ir
	
	float global_stdev = sqrt(stdev_sum/(n-1));
	
		//mark off those that are noisy or dead. 
    for(int ir = 1; ir<ir_max; ir++){
		for(int ic = 1; ic<ic_max; ic++){
			
			float pixel = map_pixels->GetBinContent(ir,ic);
			float efficiency=pixel/global_mean;
			if(is_very_noisy(pixel, global_mean, global_stdev, efficiency)) map_maskedOff.SetBinContent(ir, ic, 1);
			if(is_very_dead(pixel, global_mean, global_stdev, efficiency)) map_dead.SetBinContent(ir, ic, 1);	    
		}//end for ic
    }//end for ir
	 //END GLOBAL FILTERING 
	
	
	
	
		//BEGIN LOCAL MARKER      MMMMM
		//This section uses local measures to declare a pixel dead or noisy
	
		//mark bad pixels
	for(int ir = 1; ir<=ir_max; ir++){
		for(int ic = 1; ic<=ic_max; ic++){
			
			if(map_dead.GetBinContent(ir, ic) ==1 ||
			   map_maskedOff.GetBinContent(ir,ic)==1){
	            continue;
			}//endif, if already marked as bad, leave it
			float pixel = map_pixels->GetBinContent(ir, ic);
			
			float neighbors[]={-1,-1,-1,-1,-1,-1,
				-1,-1,-1,-1,-1,-1,
				-1,-1,-1,-1,-1,-1,
				-1,-1,-1,-1,-1,-1, -1};//25 -1's, only 24 may be filled with other things.
			int N_neighbors = get_neighbors(map_pixels,ir,ic, *map_warned, neighbors);
			if(N_neighbors > 6){
				float local_mean = get_mean(neighbors);
				float local_stdev = get_stdev(neighbors, local_mean);
				float local_efficiency = get_efficiency(pixel, local_mean);
				
				
				if(is_noisy(pixel,local_mean, local_stdev, local_efficiency))  map_maskedOff.SetBinContent(ir,ic,1);
				if(is_dead(pixel,local_mean, local_stdev, local_efficiency)) map_dead.SetBinContent(ir, ic, 1);
			}//end if N_neighbors > 6
			
				//else if a pixel has less than 6 acceptable neighbors, 
				//use global mean from prefiltering to determine its fate.
				//why 6: well, if it's a corner, there are at most 8 pixels
				//so 6 looked good. 
			
		}//end for ic
    }//end for ir
	return;
}//end Find_Dead_And_Noisy

void Find_Dead_And_Noisy(TH2F * map_pixels, TH2F & map_maskedOff, TH2F & map_dead, TH2F & map_abs_dead){
		//This is just like the other Find_Dead_And_Noisy, but we've added map_abs_dead, which is a map of pixels that are identically 0.
	
		//map_pixels is the main input, the other two can just be pointers, uninitialized
	
	/*
	 call with:
	 int nX = map_pixels->GetNbinsX();
	 int nY = map_pixels->GetNbinsY();
	 TH2F * map_maskedOff = new TH2F("map_maskedOff", "Noisy pixels to mask off", nX, 1, nX, nY, 1, nY);
	 TH2F * map_dead = new TH2F("map_dead", "dead and under-performing pixels", nX, 1, nX, nY, 1, nY);
	 Find_Dead_And_Noisy(map_pixels, *map_maskedOff, *map_dead);
	 */
	
	int ir_max = map_pixels->GetNbinsX();
	int ic_max = map_pixels->GetNbinsY();//these should be the bounds of the row and collumn for map_pixel
	
	TH2F * map_warned = new TH2F("map_warned", "warned pixels", ir_max, 1, ir_max, ic_max, 1, ic_max );
	
	Fill_loc(map_maskedOff,0);//0 is false = alive, 1 is true = noisy and masked off
	Fill_loc(map_dead,0);//0 is false = alive, 1 is true = noisy and masked off
	Fill_loc(map_abs_dead,0);//0 is false = alive, 1 is true = noisy and masked off
	Fill_loc(*map_warned,0);//0 is false = alive, 1 is true = noisy and masked off
							//map_maskedOff and map_dead are already references, while map_warned is a pointer, needs the star.
	
	
		//#define MAX_NEIGHBORS 25   // size of neighbors[]
	/*
	 The warner looks at whether each pixel looks anomolously high or low compared to its neighborhood. If so, it gives it a warning, symbolized by a 1 in map_warned. 
	 This lets later averages be computed without skew from anomolous pixels. 
	 */
	
		//BEGIN WARNER     WWWW
		//mark pixels that look like they might be noisy with a warning. 
    for(int ir = 1; ir<=ir_max; ir++){
		for(int ic = 1; ic<=ic_max; ic++){//real
										  //	 if(map_dead.(ir, ic) ==1 ||
										  //	     map_maskedOff.GetBinContent(ir,ic)==1){
										  //	    	    map_warned->SetBinContent(ir, ic, 1);
										  //		    continue;
										  //	  }//endif, if already marked as bad, mark as warned
			
			float pixel = map_pixels->GetBinContent(ir, ic);
			float neighbors[]={-1,-1,-1,-1,-1,-1,
				-1,-1,-1,-1,-1,-1,
				-1,-1,-1,-1,-1,-1,
				-1,-1,-1,-1,-1,-1, -1};//25 -1's, only 24 may be filled with other things.
			
			int N_neighbors = get_neighbors(map_pixels,ir,ic, map_maskedOff, map_dead, neighbors);
			if(N_neighbors > 6){
				
				float local_mean = get_mean(neighbors);
				float local_stdev = get_stdev(neighbors, local_mean);
				float local_efficiency = get_efficiency(pixel, local_mean);
				
				if(   is_noisy(pixel,local_mean, local_stdev, local_efficiency)  )   map_warned->SetBinContent(ir, ic, 1);
				if(   is_dead(pixel,local_mean, local_stdev, local_efficiency)  ) map_warned->SetBinContent(ir, ic, 1);
			}//end if N_neighbors>6, if a pixel has fewer than 6 valid neighbors, rely 
			
		}//end for ic
    }//end for ir
	 //END WARNER
	
	
	
	
		//BEGIN GLOBAL FILTERING
		//This section uses global measures to declare a pixel dead or noisy
	
		//find global sum, not couting pixels that look anomolous. 
	float mean_sum=0;
	int n=0;
	for(int ir = 1; ir<ir_max; ir++){
		for(int ic = 1; ic<ic_max; ic++){//these are the real ones
			if(map_warned->GetBinContent(ir,ic) != 1){
				mean_sum += map_pixels->GetBinContent(ir,ic);
				n++;
			}//endif
		}//end for ic
	}//end for ir
    float global_mean = mean_sum/n;
	
		//find global sample st.dev, not couting pixels that look anomolous
    n = 0;
	float stdev_sum=0;
    for(int ir = 1; ir<ir_max; ir++){
		for(int ic = 1; ic<ic_max; ic++){
			if(map_warned->GetBinContent(ir,ic) != 1){
				stdev_sum+=pow((map_pixels->GetBinContent(ir,ic) - global_mean),2);
				n++;
			}//endif
		}//end for ic
	}//end for ir
	
	float global_stdev = sqrt(stdev_sum/(n-1));
	
		//mark off those that are noisy or dead. 
    for(int ir = 1; ir<ir_max; ir++){
		for(int ic = 1; ic<ic_max; ic++){
			
			float pixel = map_pixels->GetBinContent(ir,ic);
			float efficiency=pixel/global_mean;
			if(is_very_noisy(pixel, global_mean, global_stdev, efficiency)) map_maskedOff.SetBinContent(ir, ic, 1);
			if(is_very_dead(pixel, global_mean, global_stdev, efficiency)) map_dead.SetBinContent(ir, ic, 1);
			if(is_abs_dead(pixel, global_mean, global_stdev, efficiency)) map_abs_dead.SetBinContent(ir, ic, 1);
		}//end for ic
    }//end for ir
	 //END GLOBAL FILTERING 
	
	
	
	
		//BEGIN LOCAL MARKER      MMMMM
		//This section uses local measures to declare a pixel dead or noisy
	
		//mark bad pixels
	for(int ir = 1; ir<=ir_max; ir++){
		for(int ic = 1; ic<=ic_max; ic++){
			
			if(map_dead.GetBinContent(ir, ic) ==1 ||
			   map_maskedOff.GetBinContent(ir,ic)==1){
	            continue;
			}//endif, if already marked as bad, leave it
			float pixel = map_pixels->GetBinContent(ir, ic);
			
			float neighbors[]={-1,-1,-1,-1,-1,-1,
				-1,-1,-1,-1,-1,-1,
				-1,-1,-1,-1,-1,-1,
				-1,-1,-1,-1,-1,-1, -1};//25 -1's, only 24 may be filled with other things.
			int N_neighbors = get_neighbors(map_pixels,ir,ic, *map_warned, neighbors);
			if(N_neighbors > 6){
				float local_mean = get_mean(neighbors);
				float local_stdev = get_stdev(neighbors, local_mean);
				float local_efficiency = get_efficiency(pixel, local_mean);
				
				
				if(is_noisy(pixel,local_mean, local_stdev, local_efficiency))  map_maskedOff.SetBinContent(ir,ic,1);
				if(is_dead(pixel,local_mean, local_stdev, local_efficiency)) map_dead.SetBinContent(ir, ic, 1);
			}//end if N_neighbors > 6
			
				//else if a pixel has less than 6 acceptable neighbors, 
				//use global mean from prefiltering to determine its fate.
				//why 6: well, if it's a corner, there are at most 8 pixels
				//so 6 looked good. 
			
		}//end for ic
    }//end for ir
	return;
}//end Find_Dead_And_Noisy



int get_neighbors(TH2F* map_pixels, int ir, int ic, TH2F & Bad1, TH2F & Bad2, float *neighbors){
		//returns an array of pixels in a 5x5 box around ic,ir (excluding ic,ir)
		//it excludes pixels that are out of range and which are marked as bad in Bad1 or Bad2
		//takes in a map of pixel values, and the coordinates of one pixel: ic, ir
		//it also takes in t
	
	
	int ir_max = map_pixels->GetNbinsX();
	int ic_max = map_pixels->GetNbinsY();//these should be the bounds of the row and collumn for map_pixel
	
	int i = 0;
	for(int R = ir-2; R<=ir+2;R++){
		for(int C = ic-2; C<=ic+2; C++){
			if(is_fiducial(R,C, ir_max, ic_max) && !(R == ir && C == ic) && Bad1.GetBinContent(R,C)!=1 && Bad2.GetBinContent(R,C)!=1 ){ 
				neighbors[i]=map_pixels->GetBinContent(R,C);
				i++;
			}//endif
		}//end for C
	}//end For R
	
	return i;
}//end get_neighbors

int get_neighbors(TH2F* map_pixels, int ir, int ic, TH2F & Bad1, float *neighbors){
		//returns an array of pixels in a 5x5 box around ic,ir (excluding ic,ir)
		//it excludes pixels that are out of range and which are marked as bad in Bad1
		//takes in a map of pixel values, and the coordinates of one pixel: ic, ir
		//it also takes in t
	
	
	int ir_max = map_pixels->GetNbinsX();
	int ic_max = map_pixels->GetNbinsY();//these should be the bounds of the row and collumn for map_pixel
	
	int i = 0;
	for(int R = ir-2; R<=ir+2;R++){
		for(int C = ic-2; C<=ic+2; C++){
			if(is_fiducial(R,C, ir_max, ic_max) && !(R == ir && C == ic) && Bad1.GetBinContent(R,C)!=1 ){ 
				neighbors[i]=map_pixels->GetBinContent(R,C);
				i++;
			}//endif
		}//end for C
	}//end For R
	
	return i;
}//end get_neighbors version 2

float get_mean(float *list){ 
		//computes mean from a float array, that terminates in -1
	int i=0;
	float sum = 0;
	while(list[i] != -1) sum += list[i++];
	return sum/i;
}//end get_mean

float get_stdev(float * list, float mean){ 
		//computes mean from a float array, assuming the array ends in -1
	int i=0;
	float sum = 0;
	
	while(list[i] != -1) sum += pow(list[i++] - mean,2);
	
	if(i==1) return sqrt(sum);
	return sqrt(sum/(i-1));
	
}//end get_stdev

inline float get_efficiency(float pixel, float local_mean){
		//computes efficiency for a pixel compared to a mean
	return pixel/local_mean;
}//end get_efficiency 

inline bool is_fiducial(int ir, int ic, int ir_max, int ic_max){
		//tells if a pixel is in the valid range. 
	return ir <= ir_max && ir > 0 && ic <= ic_max && ic > 0;
}//end is_fiducial

void Fill_loc(TH2F & hist, float z){
		//fills hist with value z
		// from inside Find dead and noisy, call with Fill_loc( hist, #); b/c the interesting hist's 
		//are already references, for one that isn't a reference, maybe call with
		//Fill_loc(*hist, ###);
	
	int nBinsX = hist.GetNbinsX();
	int nBinsY = hist.GetNbinsY();
	for(int i=1; i<=nBinsX; i++){
		for(int j=1; j<=nBinsY; j++){
			hist.SetBinContent(i,j,z);    
		}}//end fors
}//end Fill_loc 






	//in these, contain the definitions/procedures for testing whether a pixel is bad;
	//provide these with parameters about the pixel and it's surroundings. 

inline bool is_noisy(float pixel, float local_mean, float local_stdev, float local_efficiency){
		//Oleksiy reccomended 1.5
	return local_efficiency > 1.4;//Steve sugjested 1.5
}//end is_noisy

inline bool is_very_noisy(float pixel, float global_mean, float global_stdev, float global_efficiency){
		//loud pixels: those with 5 (maybe more) global-st.dev above the global mean.
	return pixel > (global_mean + 4*global_stdev);
}//end is_very_nois

inline bool is_dead(float pixel, float local_mean, float local_stdev, float local_efficiency){
		//Oleksiy reccomended 0.7, I was using 0.4
	return local_efficiency < 0.6 || pixel <= 0;//Steve sugjested 0.5
}//end is_dead

bool is_very_dead(float pixel, float global_mean, float global_stdev, float global_efficiency){
	if(pixel <= 0.25*global_stdev) return true;
	if((global_mean - 3.5*global_stdev) > 0){
		if(pixel < (global_mean - 5*global_stdev)) return true;
	}
	else if((global_mean - 3*global_stdev) > 0){
		if(pixel < (global_mean - 4*global_stdev)) return true;
	}
	else if((global_mean - 2.5*global_stdev) > 0){
		if(pixel < (global_mean - 4*global_stdev)) return true;
	}
	return false;
}//is very dead

inline bool is_abs_dead(float pixel, float global_mean, float global_stdev, float global_efficiency){
	if(pixel <= 0) return true;
	return false;
}//is very dead
