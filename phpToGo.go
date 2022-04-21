package main

import (
	"fmt"
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
var Atom_Info_DB [][]string
var Heteroatom_Info_DB [][]string


func checkArgs() bool

func TSS_3D(Point_A []float64, Point_B []float64) float64
func Stand_Deviation(arr []float64) float64
func VDW_Radius_List() map[string]float64
func Max_RES_ASA() map[string]float64
func PDB_To_Shells(path string) [][][3]float64
func Atom_Data_Extraction_Simple(PDBpath string, chosen_chain_id string, keep_heter ...bool ) bool
func Residue_Radius(safty_fact ...int) map[string]float64
func Probing_TSS_Cutoff_List(VDW_Radius_List map[string]float64, probe_size int) map[string]float64
func MAX_ATOM_ASA_List(VDW_Radius_List map[string]float64, probe_size int ) map[string]float64




func main(){
	checkArgs()
	
  /** /
	VDW_Radius_List := VDW_Radius_List()
  Max_RES_ASA     := Max_RES_ASA()
  Res_R           := Residue_Radius()

  Probing_TSS_Cutoff_List := Probing_TSS_Cutoff_List(VDW_Radius_List, probe_size)
  MAX_ATOM_ASA_List := MAX_ATOM_ASA_List(VDW_Radius_List, probe_size)


  Atom_Info_DB := map[string]float64{}
  Heteroatom_Info_DB := map[string]float64{}
  check := Atom_Data_Extraction_Simple(pdb_path, chosen_chain_id) //Atom info database
  if(check == false){
    lo.Print( "Error" )
    lo.Exit() 
  }


  Pos := [][3]int{}                // atom coordinates
  Alpha_Carbon_List := []int{}   // res_num to alpha carnbon atom number
  Atoms_in_Res := [][]int{}         // atoms list under resnum
  AN_To_Element_List := []string{}  // atom number to element list 
  ResNum_To_ResName := []string{}   // res_num to res_name

  for _, Atom_Data := range Atom_Info_DB{
    atom_number := lo.Int(Atom_Data[1])   
    atom_type   := Atom_Data[2]   
    
    res_name := Atom_Data[4]   
    res_num  := lo.Int(Atom_Data[6])   
    
    x :=              lo.Float64(Atom_Data[8] )  
    y :=              lo.Float64(Atom_Data[9] ) 
    z :=              lo.Float64(Atom_Data[10])  
    
    element_name :=   Atom_Data[13]  
	
    Pos[atom_number] = []float64{x, y, z} 
    AN_To_Element_List[atom_number] = element_name 
	
    if (atom_type == "CA" ){
      Alpha_Carbon_List[res_num] = atom_number  
      ResNum_To_ResName[res_num] = res_name  
    }

    Atoms_in_Res[res_num] = atom_number 
    
  } //for _, Atom_Data

//-------------------------------------------------------------------------
  // ###Spherical grid preparation

  path := "s642-v9_20_5.pdb" // spherical grids in pdb form
  Spherical_Grid_Set := PDB_To_Shells(path) //extraction
  group_count := len(Spherical_Grid_Set) //get number of groups of grids

  Expanded_Spherical_Grid_Set :=  [][][][3]float64{} // pre-expanded spherical grids 
  for element_name, vdw_r := range VDW_Radius_List{   
    radius := vdw_r + probe_size 
    
    for group_id, Spherical_Grid := range Spherical_Grid_Set{  
      Expanded_Spherical_Grid_Set[element_name][group_id] := [][3]float64 
      
      //expand the grids to right size accoring to vdw radious of elements
      for dot_index, Dot := range Spherical_Grid {
        Expanded_Dot[0] = Dot[0] * radius 
        Expanded_Dot[1] = Dot[1] * radius 
        Expanded_Dot[2] = Dot[2] * radius 
        
        Expanded_Spherical_Grid_Set[element_name][group_id][dot_index] = Expanded_Dot
        
      } //for dot_index, Dot := range Spherical_Grid
        
    } //for group_id, Spherical_Grid := range Spherical_Grid_Set

  } //for element_name, vdw_r := range VDW_Radius_List

//-------------------------------------------------------------------------
	Neighboring_Res_List := map[int]map[int]bool{}

	for res_num,atom_number := range Alpha_Carbon_List{

	

		for target_res_num,target_atom_number := range Alpha_Carbon_List{
			Neighboring_Res_List[res_num][res_num] = true

			if res_num >= target_res_num{
				continue
			}

			tss = TSS_3D(Pos[atom_number], Pos[target_atom_number])

			res_name        = ResNum_To_ResName[res_num]
			target_res_name = ResNum_To_ResName[target_res_num]
			
			res_name_r   = Res_R[res_name]
			target_res_r = Res_R[target_res_name]

			cutoff_dist_res =  res_name_r  + 6.4 + target_res_r
			cutoff_tss_res  =  cutoff_dist_res * cutoff_dist_res 
			
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

		for target_res_num,boolen := range Target_Res_List{

			if res_pair_log[res_num][target_res_num]{
				continue
			}

			res_pair_log[res_num][target_res_num] = true
			res_pair_log[target_res_num][res_num] = true
			
			for _,atom_number := Atoms_in_Res[res_num]{

				element_name = AN_To_Element_List[atom_number]
				vdw_r        = VDW_Radius_List[element_name]

				for _,target_atom_number := range Atoms_in_Res[target_res_num]{
					
					if(Proximity_Table[atom_number][target_atom_number] != 0){
						continue
					}

					//tss calculation
					tss = TSS_3D(Pos[atom_number], Pos[target_atom_number])
					
					//cutoff calculation

					target_element_name := AN_To_Element_List[target_atom_number]
					target_vdw_r         = VDW_Radius_List[target_element_name]

					cutoff_dist = vdw_r + 2.8 + target_vdw_r
					cutoff_tss = cutoff_dist * cutoff_dist

					if(tss <= cutoff_tss){
						Proximity_Table[atom_number][target_atom_number] = tss 
						Proximity_Table[target_atom_number][atom_number] = tss //pass info		
					} 

				}
			}
		}
	}
	/**/
}



func TSS_3D(Point_A []float64, Point_B []float64) float64 {

	d_x := Point_A[0] - Point_B[0] 
	d_y := Point_A[1] - Point_B[1] 
	d_z := Point_A[2] - Point_B[2] 	

	return d_x*d_x + d_y*d_y + d_z*d_z
} //func TSS_3D
func Stand_Deviation(arr []float64) float64 {
	
    num_of_elements := len(arr)
      
    variance := 0.0
    
    // calculating mean using array_sum() method
    sum := 0
    for _ , i := range arr {
      sum=sum+i
    }
    average := sum/num_of_elements
    
    for _ , i := range arr {
      variance += math.pow((i - average), 2)
    }

    return (float64)math.sqrt(variance/num_of_elements)
} //func Stand_Deviation

 


func Max_RES_ASA() map[string]float64 { 
	
	Max_RES_ASA := map[string]float64{}
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
	//Max_RES_ASA["SEC"] =   ??? 
  
	return Max_RES_ASA
} //func   Max_RES_ASA
func checkInput() bool {
	size := len(lo.Args)
	argv := lo.Args

	if size < 2 {
		lo.Println("No PDB path info")
		return false

	} else if File_Exists(argv[1]){
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
		probe_size = argv[3]
	}
	if size > 4{
		qcv_cutoff = argv[4]
	}
	if size > 5{
		asa_sample_size = argv[5]
	}

	return true

} //func checkInput() , true legal fasle not legal



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



func Atom_Data_Extraction_Simple(PDBpath string, chosen_chain_id string, keep_heter ...bool ) bool{
	
	keep_heteroatom := false
	if len(keep_heter) != 0 {
		keep_heteroatom = keep_heter[0]
	}
	Atom_Info_DB       := [][]string{}
  Heteroatom_Info_DB := [][]string{}

	VDW_Radius_List := VDW_Radius_List()

	Get := strings.Split(strings.Trim(lo.File_Get_contents(PDBpath),"\t\n\n\r\x0B"),"\n")

	for _,line := range Get {

		record_type := strings.Trim(line[0:6],"\t\n\n\r\x0B");  //ATOM     #0
		atom_number := strings.Trim(line[6:11],"\t\n\n\r\x0B");  //802      #1
		atom_type   := strings.Trim(line[12:16],"\t\n\n\r\x0B"); //CA       #2
		altLoc      := strings.Trim(line[16:17],"\t\n\n\r\x0B"); //         #3
		res_name    := strings.Trim(line[17:20],"\t\n\n\r\x0B"); //LYS      #4
		chain_id    := strings.Trim(line[21:22],"\t\n\n\r\x0B"); //A        #5
		res_num     := strings.Trim(line[22:26],"\t\n\n\r\x0B"); //105      #6
		iCode       := strings.Trim(line[26:27],"\t\n\n\r\x0B"); //         #7
		
		x := strings.Trim(line[30:38],"\t\n\n\r\x0B"); //30.356   #8
		y := strings.Trim(line[38:46],"\t\n\n\r\x0B"); //2.148    #9
		z := strings.Trim(line[46:54],"\t\n\n\r\x0B"); //10.394   #10

		occupancy    := strings.Trim(line[54:60],"\t\n\n\r\x0B"); //1.00     #11
		temp_factor  := strings.Trim(line[60:66],"\t\n\n\r\x0B"); //29.41    #12
		element_name := strings.Trim(line[76:78],"\t\n\n\r\x0B"); //C        #13
		charge       := strings.Trim(line[78:80],"\t\n\n\r\x0B"); //         #14

		if res_name == "HOH" || chain_id != chosen_chain_id || !(altLoc == "" || altLoc == "A"){
			continue
		}         //skip water data

		if VDW_Radius_List[element_name]==""{
			lo.Println("unknown element:",atom_number,"\t",element_name)
		  return false
		}
		
		if record_type == "ATOM" {
			Atom_Info_DB[atom_number] = []string{record_type, atom_number, atom_type, altLoc , res_name, chain_id, res_num, iCode, x, y, z, occupancy, temp_factor, element_name, charge}
		} else if record_type == "HETATM" && keep_heteroatom == true {
			Heteroatom_Info_DB[atom_number] = []string{record_type, atom_number, atom_type, altLoc , res_name, chain_id, res_num, iCode, x, y, z, occupancy, temp_factor, element_name, charge}
		} else {
			continue
		}

		if len(Atom_Info_DB) == 0{
			return false 
		} //checking
	
		return true;
	}
}



func Residue_Radius(safty_fact ...int) map[string]float64 { 

	safty_factor := 1
	if len(safty_fact) == 1 {
		safty_factor = safty_fact[0]
	}
	Res_R = map[string]float64{}  
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



func Probing_TSS_Cutoff_List(VDW_Radius_List map[string]float64, probe_size int) map[string]float64 { 
	
	Probing_TSS_Cutoff_List := map[string]float64{}
	
	for element_name, vdw_r := range VDW_Radius_List{
		radious := vdw_r + probe_size
		Probing_TSS_Cutoff_List[element_name] = radious * radious
	} //for element_name, vdw_r
	
	return Probing_TSS_Cutoff_List
} //func Probing_TSS_Cutoff_List
//Pre-calculate the tss cutoff data used in progress of  probing for availibility ofspots



func MAX_ATOM_ASA_List(VDW_Radius_List map[string]float64, probe_size int ) map[string]float64 { 

	MAX_ATOM_ASA_List := map[string]float64{}

	for element_name, vdw_r := range VDW_Radius_List{
		MAX_ATOM_ASA_List[element_name] =  4 * math.Pi * (vdw_r + probe_size) * (vdw_r + probe_size) 
	} //for element_name, vdw_r

	return MAX_ATOM_ASA_List
}
//Pre-calculate the max asa data

func VDW_Radius_List() map[string]float64 { 
	
	VDW_Radius_List :=map[string]float64{}
	VDW_Radius_List["H"] = 1.2  
	VDW_Radius_List["C"] = 1.7  
	VDW_Radius_List["N"] = 1.55 
	VDW_Radius_List["O"] = 1.52 
	VDW_Radius_List["S"] = 1.8  
	VDW_Radius_List["P"] = 1.8  
	
	VDW_Radius_List["SE"] = 1.9  
	return VDW_Radius_List

} //func VDW_Radius_List
//define VDW_radius of atoms 