package main

import (
	
	"lo"
	"strings"
	"math"
) 

var pdb_path string
var chosen_chain_id string    //Chain ID
var probe_size float64       //vdw radius of the probe 
var qcv_cutoff float64      //qaulity control value as threshold for result
var asa_sample_size int    //the number of latest asa sample to use for quality control 
var keep_heteroatom bool  //heteroatom as default
var Atom_Info_DB = make(map[int][]string)
var Heteroatom_Info_DB = make(map[int][]string)
/**/
var VDW_Radius_List = make(map[string]float64)
var Max_RES_ASA     = make(map[string]float64)
var Res_R           = make(map[string]float64)
/**/
var Probing_TSS_Cutoff_List = make(map[string]float64)
var MAX_ATOM_ASA_List       = make(map[string]float64)


func main(){
	//Atom_Info_DB = map[int][]string{}
	//Heteroatom_Info_DB = map[int][]string{}
	checkArgs()
  fillUp()
  /**/
  

  


  check := Atom_Data_Extraction_Simple() //Atom info database
  if(check == false){
    lo.Print( "Error" )
    lo.Exit() 
  }


  Pos                := map[int][3]float64{}                // atom coordinates
  Alpha_Carbon_List  := map[int]int{}   // res_num to alpha carnbon atom number
  Atoms_in_Res       := map[int][]int{}         // atoms list under resnum
  AN_To_Element_List := map[int]string{}  // atom number to element list 
  ResNum_To_ResName  := map[int]string{}   // res_num to res_name

  //lo.Println(Atom_Info_DB) // *** test
  for _, Atom_Data := range Atom_Info_DB{
    atom_number := lo.Int(Atom_Data[1])   
    atom_type   := Atom_Data[2]   
    
    res_name := Atom_Data[4]   
    res_num  := lo.Int(Atom_Data[6])   
    
    x := lo.Float64(Atom_Data[8] )  
    y := lo.Float64(Atom_Data[9] ) 
    z := lo.Float64(Atom_Data[10])  
    
    element_name :=   Atom_Data[13]  
	
    Pos[atom_number] = [3]float64{x, y, z} 
    AN_To_Element_List[atom_number] = element_name 
	
    if (atom_type == "CA" ){
      Alpha_Carbon_List[res_num] = atom_number  
      ResNum_To_ResName[res_num] = res_name  
    }

    Atoms_in_Res[res_num] = append(Atoms_in_Res[res_num],atom_number) 
    
  } //for _, Atom_Data

 //-------------------------------------------------------------------------
  // ###Spherical grid preparation
	
  path := "s642-v9_20_5.pdb" // spherical grids in pdb form
  Spherical_Grid_Set := PDB_To_Shells(path) //extraction
  group_count := len(Spherical_Grid_Set) //get number of groups of grids
  _=group_count
  //lo.Println(path) // *** test
  Expanded_Spherical_Grid_Set :=  map[string]map[int]map[int][3]float64{} // pre-expanded spherical grids 
  for element_name, vdw_r := range VDW_Radius_List{   
    radius := vdw_r + probe_size 
		  Expanded_Spherical_Grid_Set[element_name] = map[int]map[int][3]float64{}

    for group_id, Spherical_Grid := range Spherical_Grid_Set{  
      Expanded_Spherical_Grid_Set[element_name][group_id] = map[int][3]float64{}
      
      //expand the grids to right size accoring to vdw radious of elements
      for dot_index, Dot := range Spherical_Grid {
        Expanded_Dot := [3]float64{}
        Expanded_Dot[0] = Dot[0] * radius 
        Expanded_Dot[1] = Dot[1] * radius 
        Expanded_Dot[2] = Dot[2] * radius 
        
        Expanded_Spherical_Grid_Set[element_name][group_id][dot_index] = Expanded_Dot
        
      } //for dot_index, Dot := range Spherical_Grid
        
    } //for group_id, Spherical_Grid := range Spherical_Grid_Set

  } //for element_name, vdw_r := range VDW_Radius_List


//-------------------------------------------------------------------------
	Neighboring_Res_List := map[int]map[int]bool{}
  /* 獨家冠名*/
	for res_num,_ := range Alpha_Carbon_List{
    Neighboring_Res_List[res_num] = map[int]bool{}
    
  }
	//lo.Println(len(Alpha_Carbon_List))
  /**/
	for res_num,atom_number := range Alpha_Carbon_List{

   
		
		for target_res_num,target_atom_number := range Alpha_Carbon_List{
      //lo.Println(target_res_num ," ", target_atom_number )
			if(Neighboring_Res_List[res_num]==nil){
				lo.Println(res_num)
			}
			Neighboring_Res_List[res_num][res_num] = true
      
			if (res_num >= target_res_num){
				continue
			}
      
			tss := TSS_3D(Pos[atom_number], Pos[target_atom_number])
      
			res_name        := ResNum_To_ResName[res_num]
			target_res_name := ResNum_To_ResName[target_res_num]
			
	
			Res_R = Res_R
  		res_name_r   := Res_R[res_name]
			target_res_r := Res_R[target_res_name]
      
			cutoff_dist_res :=  res_name_r  + 6.4 + target_res_r
			cutoff_tss_res  :=  cutoff_dist_res * cutoff_dist_res 
			
			if(tss <= cutoff_tss_res){
        
				Neighboring_Res_List[res_num][target_res_num] = true; 
				Neighboring_Res_List[target_res_num][res_num] = true; 
				
			} //save data if residue could be next to each others
      
		}
		
	
	}
  
	
 //-----------------------------------------------------------------------
	Proximity_Table := map[int]map[int]float64{}
	res_pair_log    := map[int]map[int]bool{}

	
	

	for res_num,Target_Res_List := range Neighboring_Res_List{

    if(res_pair_log[res_num] == nil){
			res_pair_log[res_num] = map[int]bool{}
		}
		
		for target_res_num,_ := range Target_Res_List{

			if(res_pair_log[target_res_num] == nil){
				res_pair_log[target_res_num] = map[int]bool{}
			}

			if res_pair_log[res_num][target_res_num]{
				continue
			}

			res_pair_log[res_num][target_res_num] = true
			res_pair_log[target_res_num][res_num] = true
			
			for _,atom_number := range Atoms_in_Res[res_num]{
				if(Proximity_Table[atom_number] == nil){
					Proximity_Table[atom_number] = map[int]float64{}
				}
				
				element_name := AN_To_Element_List[atom_number]
				vdw_r        := VDW_Radius_List[element_name]

				for _,target_atom_number := range Atoms_in_Res[target_res_num]{
					if(Proximity_Table[target_atom_number] == nil){
						Proximity_Table[target_atom_number] = map[int]float64{}
					}
					
					if(Proximity_Table[atom_number][target_atom_number] != 0){
						continue
					}

					//tss calculation
					tss := TSS_3D(Pos[atom_number], Pos[target_atom_number])
					
					//cutoff calculation

					target_element_name := AN_To_Element_List[target_atom_number]
					target_vdw_r        := VDW_Radius_List[target_element_name]

					cutoff_dist := vdw_r + 2.8 + target_vdw_r
					cutoff_tss := cutoff_dist * cutoff_dist

					if(tss <= cutoff_tss){
						Proximity_Table[atom_number][target_atom_number] = tss 
						Proximity_Table[target_atom_number][atom_number] = tss //pass info		
					} 

				}
			}
		}
	}



	Temp := map[int](map[int]float64){}

	for atom_number,Target_Data := range Proximity_Table {
	
		Temp[atom_number] = map[int]float64{}
		for target_atom_number,tss := range Target_Data{
			if (atom_number == target_atom_number){
				continue
			}

			Temp[atom_number][target_atom_number] = tss;
		} //for target_atom_number,tss := range Target_Data
		
		//!!!!!!
		//lo.SortFloat64Slice(Temp[atom_number])
	
	} //for atom_number,Target_Data := range Proximity_Table 

	Proximity_Table = Temp

  

	Log_Refined_ASA               := map[int]float64{}
  Resisue_ASA_Unnormalized_Sum  := map[int]float64{}
  

  for res_num, _ := range ResNum_To_ResName {
		Resisue_ASA_Unnormalized_Sum[res_num] = 0;
  } //for res_num, res_name := range ResNum_To_ResName

	protein_asa_unnormalized := 0.0 //unnormalized protein_asa
	used_dot_count           := 0.0  //number of dots of used spherical grid groups   

	
	for group_id := 0; group_id < ( asa_sample_size -1 ) ; group_id++{

		Expanded_Spherical_Grid := Expanded_Spherical_Grid_Set["C"][group_id]; 
		dot_count               := lo.Float64(len(Expanded_Spherical_Grid))
		
		
		
		for atom_number,Target_AN_List := range Proximity_Table{	
		
			x := Pos[atom_number][0]
			y := Pos[atom_number][1]
			z := Pos[atom_number][2]
			
			element_name := AN_To_Element_List[atom_number];
			
			Expanded_Spherical_Grid = Expanded_Spherical_Grid_Set[element_name][group_id] 
	
			//Shifting spherical grid around target atom		
		
	    Shifted_Dot_List := map[int]([3]float64){}
			for num,Expanded_Dot := range Expanded_Spherical_Grid{
				
				Shifted_Dot_List[num] = [3]float64{0 : x + Expanded_Dot[0], 1:y + Expanded_Dot[1], 2:z + Expanded_Dot[2]};
				
			}//foreach($Expanded_Spherical_Grid as $Expanded_Dot){
				
				
			//Probing solvent available spots  	

			for target_atom_number, tss := range Target_AN_List{
				_=tss
				target_element_name := AN_To_Element_List[target_atom_number]
				Target_Pos          := Pos[target_atom_number]
				
				probe_tss_cutoff    := Probing_TSS_Cutoff_List[target_element_name] ;
				
				for index, Shifted_Dot := range Shifted_Dot_List{
					if TSS_3D(Target_Pos, Shifted_Dot) < probe_tss_cutoff {
						
						delete(Shifted_Dot_List,index) 
             
					}
					
				}//foreach($Shifted_Dot as $Shifted_Dot){
			}//foreach($Target_AN_List as $target_atom_number => $tss){
			
			lefted_dot_count      := lo.Float64(len(Shifted_Dot_List))
			atom_asa_unnormalized := MAX_ATOM_ASA_List[element_name] * lefted_dot_count //unnormalized atom asa
			protein_asa_unnormalized += atom_asa_unnormalized //saving
			
			res_num := lo.Int(Atom_Info_DB[atom_number][6])
			Resisue_ASA_Unnormalized_Sum[res_num] += atom_asa_unnormalized //saving
		
		
			
		}//foreach($Proximity_Table as $atom_number => $Target_AN_List){
		
			
		used_dot_count += dot_count //number of dots of used spherical grid groups   
		
		refined_asa := protein_asa_unnormalized / lo.Float64(used_dot_count) //asa 
		Log_Refined_ASA[group_id] = refined_asa //saving asa for quality control 

	
	}



	var refined_asa float64 
	for group_id := asa_sample_size -1 ; group_id < group_count; group_id++ {

		Expanded_Spherical_Grid := Expanded_Spherical_Grid_Set["C"][group_id] 
		dot_count := len(Expanded_Spherical_Grid)
	
	
		for atom_number, Target_AN_List := range Proximity_Table{
			x := Pos[atom_number][0] 
			y := Pos[atom_number][1]
			z := Pos[atom_number][2]
			
			element_name := AN_To_Element_List[atom_number];
			
			Expanded_Spherical_Grid =  Expanded_Spherical_Grid_Set[element_name][group_id]; 
	
			
			//##Shifting spherical grid around target atom
			Shifted_Dot_List := map[int]([3]float64){}

			for num, Expanded_Dot := range Expanded_Spherical_Grid{

				Shifted_Dot_List[num] = [3]float64{0 : x + Expanded_Dot[0], 1:y + Expanded_Dot[1], 2:z + Expanded_Dot[2]};
			} //for num, Expanded_Dot := range Expanded_Spherical_Grid
				
			//Probing solvent available spots  		
			for target_atom_number, _ := range Target_AN_List{
				target_element_name := AN_To_Element_List[target_atom_number];
				Target_Pos          := Pos[target_atom_number];
				
				
				probe_tss_cutoff := Probing_TSS_Cutoff_List[target_element_name] ;
				
				for index, Shifted_Dot := range Shifted_Dot_List{

					if( TSS_3D(Target_Pos, Shifted_Dot) < probe_tss_cutoff ){

						delete(Shifted_Dot_List,index) 
					}
					
				}
			}//for target_atom_number, tss := range Target_AN_List{
			
			lefted_dot_count      := len(Shifted_Dot_List);
			atom_asa_unnormalized := MAX_ATOM_ASA_List[element_name] * lo.Float64(lefted_dot_count);
			protein_asa_unnormalized += atom_asa_unnormalized; //saving
			
			
			res_num := Atom_Info_DB[atom_number][6];
			Resisue_ASA_Unnormalized_Sum[lo.Int(res_num)] += atom_asa_unnormalized; //saving
			
			
		}//foreach($Proximity_Table as $atom_number => $Target_AN_List){

		
		used_dot_count += lo.Float64(dot_count) //number of dots of used spherical grid groups   
		refined_asa = protein_asa_unnormalized / lo.Float64(used_dot_count) //asa 

		Log_Refined_ASA[group_id] = refined_asa ; //saving asa for quality control 
		
		//$real_shell_num = $group_id + 1;
			
		//take latest asa results for quality check 
		
		Array_ASA_Samples := map[int]float64{}
		count := 0
		for sampe_group_id := group_id - asa_sample_size + 1 ; sampe_group_id < group_id + 1; sampe_group_id++ {
			Array_ASA_Samples[count] = Log_Refined_ASA[sampe_group_id] 
			count++
		}
		
		sd_asa_sample := Stand_Deviation(Array_ASA_Samples);
		cv := sd_asa_sample / refined_asa;
		if(cv < qcv_cutoff){break}  //pass quality checking
	
	}
		// Print out asa
	refined_asa = lo.Round(refined_asa, 3);

	
	// Get RSA and print out

	
  //lo.Println(Resisue_ASA_Unnormalized_Sum)  // *** test
	lo.Println("ASA: ",refined_asa) 
	lo.Println("RSA:") 
	lo.Println("Residue Number\tRSA") 

  keys := make([]int,0,len(Resisue_ASA_Unnormalized_Sum))
	for key,_ := range Resisue_ASA_Unnormalized_Sum{
		keys = append(keys,key)
	}
	lo.SortIntSlice(keys)


	for _, res_num := range keys{
	
		res_name := ResNum_To_ResName[res_num]
		
		
	
		residue_asa_unnormalized := Resisue_ASA_Unnormalized_Sum[res_num]
		if Max_RES_ASA[res_name]!=0{
			
			refined_residue_asa :=  residue_asa_unnormalized / lo.Float64(used_dot_count)
			
			refined_residue_rsa := refined_residue_asa / Max_RES_ASA[res_name];
			refined_residue_rsa = lo.Round(refined_residue_rsa, 3)
			lo.Println( lo.String(res_num) + "\t" + lo.String(refined_residue_rsa) )
			
		}else{
			lo.Println( lo.String(res_num) + "\tError" )
		}
	
	} //for res_num, residue_asa_unnormalized := range Resisue_ASA_Unnormalized_Sum{
	
	/**/
	//lo.Println(Atom_Data_Extraction_Simple("test.pdb", chosen_chain_id ) )
	//Atom_Data_Extraction_Simple(pdb_path,chosen_chain_id)
}




///---------------------------------------------------------------------------
///------------------------------------------------------------------------
func TSS_3D( Point_A [3]float64, Point_B [3]float64) float64 {

	d_x := Point_A[0] - Point_B[0] 
	d_y := Point_A[1] - Point_B[1] 
	d_z := Point_A[2] - Point_B[2] 	

	return d_x*d_x + d_y*d_y + d_z*d_z
} //func TSS_3D



func checkArgs() bool {
	size := len(lo.Args)
	argv := lo.Args

	if size < 2 {
		lo.Println("No PDB path info")
		return false

	} else if lo.File_Exists(argv[1]){
		pdb_path = argv[1]

	} else {
		lo.Println("PDB Unfound")
		return false
	}

	chosen_chain_id = "A"
	probe_size      = 1.4
	qcv_cutoff      = 0.005
	asa_sample_size = 3
	keep_heteroatom = false	

	if size > 2{
		chosen_chain_id = argv[2]
	} 
	if size > 3{
		probe_size = lo.Float64(argv[3])
	}
	if size > 4{
		qcv_cutoff = lo.Float64(argv[4])
	}
	if size > 5{
		asa_sample_size = lo.Int(argv[5])
	}

	return true

} //func checkInput() , true legal fasle not legal

func fillUp(){

	VDW_Radius_List["H"]  = 1.2  
	VDW_Radius_List["C"]  = 1.7  
	VDW_Radius_List["N"]  = 1.55 
	VDW_Radius_List["O"]  = 1.52 
	VDW_Radius_List["S"]  = 1.8  
	VDW_Radius_List["P"]  = 1.8  	
	VDW_Radius_List["SE"] = 1.9  

	Max_RES_ASA["ALA"] = 121 
	Max_RES_ASA["ARG"] = 265 
	Max_RES_ASA["ASN"] = 187 
	Max_RES_ASA["ASP"] = 187 
	Max_RES_ASA["CYS"] = 148 
	Max_RES_ASA["GLU"] = 214 
	Max_RES_ASA["GLN"] = 214 
	Max_RES_ASA["GLY"] = 97  
	Max_RES_ASA["HIS"] = 216 
	Max_RES_ASA["ILE"] = 195 
	Max_RES_ASA["LEU"] = 191 
	Max_RES_ASA["LYS"] = 230 
	Max_RES_ASA["MET"] = 203 
	Max_RES_ASA["PHE"] = 228 
	Max_RES_ASA["PRO"] = 154 
	Max_RES_ASA["SER"] = 143 
	Max_RES_ASA["THR"] = 163 
	Max_RES_ASA["TRP"] = 264 
	Max_RES_ASA["TYR"] = 255 
	Max_RES_ASA["VAL"] = 165 

	Res_R              = Residue_Radius()

	for element_name, vdw_r := range VDW_Radius_List{
		radious := vdw_r + lo.Float64(probe_size)
		Probing_TSS_Cutoff_List[element_name] = radious * radious
	} //for element_name, vdw_r
	//Pre-calculate the tss cutoff data used in progress of  probing for availibility ofspots
	for element_name, vdw_r := range VDW_Radius_List{
		MAX_ATOM_ASA_List[element_name] =  4 * math.Pi * (vdw_r + lo.Float64(probe_size)) * (vdw_r + lo.Float64(probe_size)) 
	} //for element_name, vdw_r
	//Pre-calculate the max asa data

}

func Stand_Deviation(arr map[int]float64) float64 {
	
    num_of_elements := float64(len(arr))
      
    variance := 0.0
    
    // calculating mean using array_sum() method
    sum := 0.0

    for _ , i := range arr {
      sum = sum + i
    }
    average := sum / num_of_elements
    
    for _ , i := range arr {
      variance += math.Pow((i - average), 2)
    }

    return lo.Float64(math.Sqrt(variance/num_of_elements))
} //func Stand_Deviation



func PDB_To_Shells(path string) [][][3]float64 {
	radius := 15.999983777


	get    := strings.Trim(lo.File_Get_Contents(path),"\t\n\n\r\x0B" )
  Chains := strings.Split( get,"\nTER\n")

	Chains[len(Chains)-1] = strings.TrimRight(Chains[len(Chains)-1],"\nTER")

	Shell_Set := [][][3]float64{}

	for _, chain_data := range Chains{
		Lines := strings.Split(chain_data,"\n")

		Shell := [][3]float64{}
		
		for _, line := range Lines {

			x := lo.Float64(strings.Trim(line[30:38],"\t\n\n\r\x0B" ))
			y := lo.Float64(strings.Trim(line[38:46],"\t\n\n\r\x0B" ))
			z := lo.Float64(strings.Trim(line[46:54],"\t\n\n\r\x0B" ))

			fixed_x := x / radius
			fixed_y := y / radius
			fixed_z := z / radius

			Shell = append(Shell,[3]float64{fixed_x,fixed_y,fixed_z});
		}

		Shell_Set = append(Shell_Set,Shell)
	}

	return Shell_Set
}






func Atom_Data_Extraction_Simple( keep_heter ... bool) bool {	
  keep_heteroatom := false
	if len(keep_heter) != 0 {
		keep_heteroatom = keep_heter[0]
	}
  
	//Atom_Info_DB       := map[int]([]string){}
  //Heteroatom_Info_DB := map[int]([]string){}
  

	
	Get := lo.File(pdb_path); 
	//for _, line := range Get
  for _, line := range Get{
		//Example: ATOM    802  CA  LYS A 105      30.356   2.148  10.394  1.00 29.41           C  
    record_type := lo.Trim(line[0:6])   //ATOM     #0
		atom_number := lo.Trim(line[6:11]) //802      #1
		atom_type   := lo.Trim(line[12:16]) //CA       #2
		altLoc      := lo.Trim(line[16:17]) //         #3
		res_name    := lo.Trim(line[17:20]) //LYS      #4
		chain_id    := lo.Trim(line[21:22]) //A        #5
		res_num     := lo.Trim(line[22:26]) //105      #6
		iCode       := lo.Trim(line[26:27]) //         #7
		
    x := lo.Trim(line[30:38]) //30.356   #8
		y := lo.Trim(line[38:46]) //2.148    #9
		z := lo.Trim(line[46:54]) //10.394   #10
    
    occupancy    := lo.Trim(line[54:60]) //1.00     #11
		temp_factor  := lo.Trim(line[60:66]) //29.41    #12
		element_name := lo.Trim(line[76:78]) //C        #13
		charge       := "?" //         #14

    /*test* /
    pri:=[]string{
        record_type, atom_number, atom_type, altLoc , res_name, 
        chain_id, res_num, iCode, x, y, z, occupancy, temp_factor, element_name, charge,
    }
		lo.Println(pri)
    /**/
    
    if (res_name == "HOH"){continue}                  //skip water data
		if (chain_id != chosen_chain_id){continue}       //skip unchosen chains 
		if(!(altLoc == "" || altLoc == "A")){continue}	//choose only single altLoc
    
    var check bool = false
    for keys, values := range  VDW_Radius_List{
        _=values
        if keys == element_name{
            check = true
            break
        }
    }
		if(!check){
			lo.Println("unknown element:"+lo.String(atom_number)+"\t"+element_name)
      return false
		} //abort and return false if uncommon atom is in the protein
		
		if(record_type == "ATOM"){
			//save infomation of atoms
			Atom_Info_DB[ lo.Int(atom_number) ] = []string{
        record_type, atom_number, atom_type, altLoc , res_name, 
        chain_id, res_num, iCode, x, y, z, occupancy, temp_factor, element_name, charge,
      }
		}else if((record_type == "HETATM") && (keep_heteroatom == true)){
			//save infomation of heteroatoms
			Heteroatom_Info_DB[ lo.Int(atom_number) ] = []string{
        record_type, atom_number, atom_type, altLoc , res_name, 
        chain_id, res_num, iCode, x, y, z, occupancy, temp_factor, element_name, charge,
      }
		}else{continue}
		//saving atom data 
		
	} //for _, line := range Get
	
	
	if ( len(Atom_Info_DB) == 0){ return false } //checking
	return true
  
}





func Residue_Radius(safty_fact ...int) map[string]float64 { 

	safty_factor := 1.0
	if len(safty_fact) == 1 {
		safty_factor = lo.Float64(safty_fact[0])
	}
	Res_R := map[string]float64{}  
	Res_R["CYS"] = 3.881 * safty_factor 
	Res_R["GLU"] = 5.192 * safty_factor
	Res_R["HIS"] = 5.816 * safty_factor
	Res_R["ASP"] = 4.369 * safty_factor
	Res_R["ARG"] = 8.241 * safty_factor
	Res_R["LEU"] = 4.873 * safty_factor
	Res_R["GLN"] = 5.862 * safty_factor
	Res_R["SER"] = 3.331 * safty_factor
	Res_R["GLY"] = 2.574 * safty_factor
	Res_R["ILE"] = 4.746 * safty_factor
	Res_R["THR"] = 3.854 * safty_factor
	Res_R["PHE"] = 6.342 * safty_factor
	Res_R["LYS"] = 7.245 * safty_factor
	Res_R["ALA"] = 2.712 * safty_factor
	Res_R["MET"] = 6.280 * safty_factor
	Res_R["TYR"] = 7.133 * safty_factor
	Res_R["VAL"] = 3.555 * safty_factor
	Res_R["TRP"] = 7.755 * safty_factor
	Res_R["PRO"] = 3.394 * safty_factor
	Res_R["ASN"] = 4.576 * safty_factor
	
	Res_R["SEC"] = 2.867 * safty_factor
	
	return Res_R
}
//Define possible distance from CA to farest atom in each residue for 










