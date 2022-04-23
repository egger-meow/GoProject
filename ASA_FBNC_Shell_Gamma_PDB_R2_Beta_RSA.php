<?php


##Settings##

if(!Isset($argv[1])){echo "No PDB path info\n"; Exit;}else{$pdb_path = $argv[1];} 
if(!File_Exists($pdb_path)){echo "PDB Unfound\n"; Exit;} #checking pdb existance

$chosen_chain_id  = Isset($argv[2])? $argv[2] : "A" ;  #set Chain ID
$probe_size = Isset($argv[3])? $argv[3] : 1.4 ;        #set vdw radius of the probe 
$qcv_cutoff = Isset($argv[4])? $argv[4] : 0.005 ;      #set qaulity control value as threshold for result
$asa_sample_size = Isset($argv[5])? $argv[5] : 3 ;     #set the number of latest asa sample to use for quality control 
$keep_heteroatom = False; #skip heteroatom as default


##Functions##


function TSS_3D($Point_A , $Point_B){

	$d_x = $Point_A[0] - $Point_B[0];
	$d_y = $Point_A[1] - $Point_B[1];
	$d_z = $Point_A[2] - $Point_B[2];	

	return $d_x*$d_x + $d_y*$d_y + $d_z*$d_z;
}
#the function for total sum of square(TSS) for three-dimensional space


function Stand_Deviation($arr){

	
    $num_of_elements = count($arr);
      
    $variance = 0.0;
    
       // calculating mean using array_sum() method
    $average = array_sum($arr)/$num_of_elements;
      
    foreach($arr as $i)
    {
        // sum of squares of differences between 
         // all numbers and means.
        $variance += pow(($i - $average), 2);
    }
      
    return (float)sqrt($variance/$num_of_elements);
}
#the function for standard deviation
#credit: https://www.geeksforgeeks.org/php-program-find-standard-deviation-array/#:~:text=To%20calculate%20the%20standard%20deviation,%E2%88%9A(variance%2Fno_of_elements).



function PDB_To_Shells($path){
	
	$radius = 16; #15.999983777
	
	if (! function_exists("array_key_last")) {
		function array_key_last($array) {
			if (!is_array($array) || empty($array)) {
				return NULL;
			}

      $Keys = array_keys($array);
		
			return $Keys[ count($array) - 1 ];
		}
	}
	
	$get = File_Get_Contents($path);
	$get = Trim($get);
	
	$Chains = Explode("\nTER\n", $get);
	
	$key = array_key_last($Chains) ;
	
	
	$Chains[$key] = rTrim($Chains[$key], "\nTER");
	#echo $Chains[$key];
	
	$Shell_Set = Array();
	foreach ($Chains as $chain_index => $chain_data){
		
		$Lines = Explode("\n", $chain_data);
		
		
		$Shell = Array();
		foreach ($Lines as $line){
			
			$x = (float)Trim(SubStr($line, 30, 8)); 
			$y = (float)Trim(SubStr($line, 38, 8)); 
			$z = (float)Trim(SubStr($line, 46, 8)); 	
			
			$fixed_x = $x / $radius;
			$fixed_y = $y / $radius;
			$fixed_z = $z / $radius;
			
			$Shell[] = Array($fixed_x, $fixed_y, $fixed_z);
			
		}
		
		
		$Shell_Set[] = $Shell;
		
	}
	
	Return($Shell_Set);
}
#this extracts the Spherical Grids data in PDB format


function Atom_Data_Extraction_Simple($pdb_path, $chosen_chain_id, $keep_heteroatom = False){	
	
	Global $Atom_Info_DB;
	Global $Heteroatom_Info_DB;
	
	
	$VDW_Radius_List = VDW_Radius_List(); #set up to check if any uncommon atom is in the protein
	
	$Get = Explode("\n",Trim(File_Get_contents($pdb_path))); 
	print_r(strlen($Get[0]));
	foreach($Get as $line){
		
		//Example: ATOM    802  CA  LYS A 105      30.356   2.148  10.394  1.00 29.41           C  
		$record_type =         Trim(SubStr($line, 0, 6));  //ATOM     #0
		$atom_number =    (int)Trim(SubStr($line, 6, 5));  //802      #1
		$atom_type =           Trim(SubStr($line, 12, 4)); //CA       #2
		$altLoc =	           Trim(SubStr($line, 16, 1)); //         #3
		$res_name =            Trim(SubStr($line, 17, 3)); //LYS      #4
		$chain_id =            Trim(SubStr($line, 21, 1)); //A        #5
		$res_num =        (int)Trim(SubStr($line, 22, 4)); //105      #6
		$iCode =               Trim(SubStr($line, 26, 1)); //         #7
		$x =            (float)Trim(SubStr($line, 30, 8)); //30.356   #8
		$y =            (float)Trim(SubStr($line, 38, 8)); //2.148    #9
		$z =            (float)Trim(SubStr($line, 46, 8)); //10.394   #10
		$occupancy =    (float)Trim(SubStr($line, 54, 6)); //1.00     #11
		$temp_factor =  (float)Trim(SubStr($line, 60, 6)); //29.41    #12
		$element_name =        Trim(SubStr($line, 76, 2)); //C        #13
		$charge       =        Trim(SubStr($line, 78, 2)); //         #14
		
		
		if ($res_name == "HOH"){continue;}                  #skip water data
		if ($chain_id != $chosen_chain_id){continue;}       #skip unchosen chains 
		if(!($altLoc == "" || $altLoc == "A")){continue;}	#choose only single altLoc

		if(!Isset($VDW_Radius_List[$element_name])){
			echo "unknown element:$atom_number\t$element_name"; return False;
		} #abort and return false if uncommon atom is in the protein
		echo "$record_type\n";
		print_r(Array($record_type, $atom_number, $atom_type, $altLoc , $res_name, $chain_id, $res_num, $iCode, $x, $y, $z, $occupancy, $temp_factor, $element_name, $charge));
		if($record_type == "ATOM"){
			#save infomation of atoms
			$Atom_Info_DB[$atom_number] = Array($record_type, $atom_number, $atom_type, $altLoc , $res_name, $chain_id, $res_num, $iCode, $x, $y, $z, $occupancy, $temp_factor, $element_name, $charge);
		}
		elseif(($record_type == "HETATM") && ($keep_heteroatom === True)){
			#save infomation of heteroatoms
			$Heteroatom_Info_DB[$atom_number] = Array($record_type, $atom_number, $atom_type, $altLoc , $res_name, $chain_id, $res_num, $iCode, $x, $y, $z, $occupancy, $temp_factor, $element_name, $charge);
		}
		else{continue;}
		#saving atom data 
		
		
		
	}//foreach($Get as $line)
	
	
	if ($Atom_Info_DB === Array()){ Return False; } #checking
	
	return True;
}
#this function extract data of atoms in PDB file


function VDW_Radius_List(){ 
	
	$VDW_Radius_List = Array();
	$VDW_Radius_List["H"] = 1.2  ;
	//$VDW_Radius_List["H"] = 1  ;
	$VDW_Radius_List["C"] = 1.7  ;
	$VDW_Radius_List["N"] = 1.55 ;
	$VDW_Radius_List["O"] = 1.52 ;
	$VDW_Radius_List["S"] = 1.8  ;
	$VDW_Radius_List["P"] = 1.8  ;
	
	$VDW_Radius_List["SE"] = 1.9  ;
	Return $VDW_Radius_List;
}
#define VDW_radius of atoms 

function Max_RES_ASA(){ 
	
	$Max_RES_ASA = Array() ; 
	$Max_RES_ASA["ALA"] =   121 ;
	$Max_RES_ASA["ARG"] =   265 ;
	$Max_RES_ASA["ASN"] =   187 ;
	$Max_RES_ASA["ASP"] =   187 ;
	$Max_RES_ASA["CYS"] =   148 ;
	$Max_RES_ASA["GLU"] =   214 ;
	$Max_RES_ASA["GLN"] =   214 ;
	$Max_RES_ASA["GLY"] =   97  ;
	$Max_RES_ASA["HIS"] =   216 ;
	$Max_RES_ASA["ILE"] =   195 ;
	$Max_RES_ASA["LEU"] =   191 ;
	$Max_RES_ASA["LYS"] =   230 ;
	$Max_RES_ASA["MET"] =   203 ;
	$Max_RES_ASA["PHE"] =   228 ;
	$Max_RES_ASA["PRO"] =   154 ;
	$Max_RES_ASA["SER"] =   143 ;
	$Max_RES_ASA["THR"] =   163 ;
	$Max_RES_ASA["TRP"] =   264 ;
	$Max_RES_ASA["TYR"] =   255 ;
	$Max_RES_ASA["VAL"] =   165 ;
	//Max_RES_ASA["SEC"] =   ??? ;

	Return $Max_RES_ASA;
}
##Define MaxASA of residues

function Residue_Radius($safty_factor = 1){ 

	$Res_R = Array() ; 
	$Res_R["CYS"] = 3.881 * $safty_factor; 
	$Res_R["GLU"] = 5.192 * $safty_factor;
	$Res_R["HIS"] = 5.816 * $safty_factor;
	$Res_R["ASP"] = 4.369 * $safty_factor;
	$Res_R["ARG"] = 8.241 * $safty_factor;
	$Res_R["LEU"] = 4.873 * $safty_factor;
	$Res_R["GLN"] = 5.862 * $safty_factor;
	$Res_R["SER"] = 3.331 * $safty_factor;
	$Res_R["GLY"] = 2.574 * $safty_factor;
	$Res_R["ILE"] = 4.746 * $safty_factor;
	$Res_R["THR"] = 3.854 * $safty_factor;
	$Res_R["PHE"] = 6.342 * $safty_factor;
	$Res_R["LYS"] = 7.245 * $safty_factor;
	$Res_R["ALA"] = 2.712 * $safty_factor;
	$Res_R["MET"] = 6.280 * $safty_factor;
	$Res_R["TYR"] = 7.133 * $safty_factor;
	$Res_R["VAL"] = 3.555 * $safty_factor;
	$Res_R["TRP"] = 7.755 * $safty_factor;
	$Res_R["PRO"] = 3.394 * $safty_factor;
	$Res_R["ASN"] = 4.576 * $safty_factor;
	
	$Res_R["SEC"] = 2.867 * $safty_factor;
	
	return $Res_R;
}
#Define possible distance from CA to farest atom in each residue for 





function Probing_TSS_Cutoff_List($VDW_Radius_List, $probe_size){ 
	
	$Probing_TSS_Cutoff_List = Array();
	
	foreach ($VDW_Radius_List as $element_name => $vdw_r){
		$radious = $vdw_r + $probe_size;
		$Probing_TSS_Cutoff_List[$element_name] = $radious * $radious;
	}//foreach ($VDW_Radius_List as $element_name => $vdw_r){
	
	return $Probing_TSS_Cutoff_List;
}
#Pre-calculate the tss cutoff data used in progress of  probing for availibility ofspots

function MAX_ATOM_ASA_List($VDW_Radius_List, $probe_size){ 

	$MAX_ATOM_ASA_List = Array();
	
	foreach ($VDW_Radius_List as $element_name => $vdw_r){
		$MAX_ATOM_ASA_List[$element_name] =  4 * Pi() * ($vdw_r + $probe_size) * ($vdw_r + $probe_size) ;
	}//foreach ($VDW_Radius_List as $element_name => $vdw_r){
	
	return $MAX_ATOM_ASA_List;
}
#Pre-calculate the max asa data




##Bulding Dictionaries
/*
$VDW_Radius_List = VDW_Radius_List();
$Max_RES_ASA = Max_RES_ASA();
$Res_R = Residue_Radius();

$Probing_TSS_Cutoff_List = Probing_TSS_Cutoff_List($VDW_Radius_List, $probe_size);
Print_r($VDW_Radius_List);
echo $probe_size;


$MAX_ATOM_ASA_List = MAX_ATOM_ASA_List($VDW_Radius_List, $probe_size);
Print_r ($MAX_ATOM_ASA_List) ;

$Atom_Info_DB = Array();
$Heteroatom_Info_DB = Array();
$check = Atom_Data_Extraction_Simple($pdb_path, $chosen_chain_id); #Atom info database
if($check == False){echo "Error"; exit;}


$Pos = Array();                 # atom coordinates
$Alpha_Carbon_List = Array();   # res_num to alpha carnbon atom number
$Atoms_in_Res = Array();        # atoms list under resnum
$AN_To_Element_List = Array();  # atom number to element list 
$ResNum_To_ResName = Array();   # $res_num to res_name

foreach ($Atom_Info_DB as $Atom_Data){
	
	$atom_number =    $Atom_Data[1];   
	$atom_type =      $Atom_Data[2];   
	
	$res_name =       $Atom_Data[4];   
	$res_num =        $Atom_Data[6];   
	
	$x =              $Atom_Data[8];   
	$y =              $Atom_Data[9];   
	$z =              $Atom_Data[10];  
	
	$element_name =   $Atom_Data[13];  
	
	$Pos[$atom_number] = Array($x, $y, $z); 
	$AN_To_Element_List[$atom_number] = $element_name; 
	
	if ($atom_type == "CA" ){
		$Alpha_Carbon_List[$res_num] = $atom_number;  
		$ResNum_To_ResName[$res_num] = $res_name;  
	}

	$Atoms_in_Res[$res_num][] = $atom_number; 
	
}//foreach ($Atom_Info_DB as $Atom_Data){




##Spherical grid preparation

$path = "s642-v9_20_5.pdb"; # spherical grids in pdb form
$Spherical_Grid_Set = PDB_To_Shells($path); #extraction
$group_count = Count($Spherical_Grid_Set); #get number of groups of grids

$Expanded_Spherical_Grid_Set =  Array(); ## pre-expanded spherical grids 
foreach($VDW_Radius_List as $element_name => $vdw_r){
	
	$radius = $vdw_r + $probe_size; 
	
	foreach($Spherical_Grid_Set as $group_id => $Spherical_Grid){
		
		$Expanded_Spherical_Grid_Set[$element_name][$group_id] = Array(); 
		
		#expand the grids to right size accoring to vdw radious of elements
		foreach($Spherical_Grid as $dot_index => $Dot){
			
			
			$Expanded_Dot[0] = $Dot[0] * $radius; 
			$Expanded_Dot[1] = $Dot[1] * $radius; 
			$Expanded_Dot[2] = $Dot[2] * $radius; 
			
			$Expanded_Spherical_Grid_Set[$element_name][$group_id][$dot_index] = $Expanded_Dot;
			
		}//foreach($Spherical_Grid as $dot_index => $Dot){
			
	}//foreach($Spherical_Grid_Set as $group_id => $Spherical_Grid){


}//foreach($VDW_Radius_List as $element_name => $vdw_r){




##Get data of neighboring residues 

$Neighboring_Res_List = Array();
foreach($Alpha_Carbon_List as $res_num => $atom_number){	
	

	if(!Isset($Neighboring_Res_List[$res_num])){
		$Neighboring_Res_List[$res_num] = Array();
	}
	
	
	foreach($Alpha_Carbon_List as $target_res_num => $target_atom_number){	
		
		
		$Neighboring_Res_List[$res_num][$res_num] = TRUE; ##add residue self pairing first
		
		if ($res_num >= $target_res_num){continue;}	 ##only go downward
		
		
		$tss = TSS_3D($Pos[$atom_number], $Pos[$target_atom_number]); ##tss calculation between CA


		$res_name        =            $ResNum_To_ResName[$res_num] ;  
		$target_res_name =     $ResNum_To_ResName[$target_res_num] ;  

		$res_name_r      =         $Res_R[$res_name] ; 
		$target_res_r    =  $Res_R[$target_res_name] ; 
	
		$cutoff_dist_res =  $res_name_r  + 6.4 + $target_res_r ;  
		$cutoff_tss_res =  $cutoff_dist_res * $cutoff_dist_res ;
		
		if($tss <= $cutoff_tss_res){
			
			$Neighboring_Res_List[$res_num][$target_res_num] = TRUE; 
			$Neighboring_Res_List[$target_res_num][$res_num] = TRUE; 
			
		} #save data if residue could be next to each others
		
		
	} //foreach($Alpha_Carbon_List as $target_res_num => $target_atom_number){	

} // foreach($Alpha_Carbon_List as $res_num => $atom_number){	



##Get data of neighboring atoms 
$Proximity_Table  = Array();  ## database of neighboring atoms and their TSS 
$res_pair_log = Array(); ## log of calculated residue pairs
 
foreach($Neighboring_Res_List as $res_num => $Target_Res_List){		
	
	foreach($Target_Res_List as $target_res_num => $boolen){
		
		if(Isset($res_pair_log[$res_num][$target_res_num])){continue;} # avoid repeating calculation
		$res_pair_log[$res_num][$target_res_num] = True;
		$res_pair_log[$target_res_num][$res_num] = True;
		
		foreach($Atoms_in_Res[$res_num] as $atom_number ){
			
			$element_name = $AN_To_Element_List[$atom_number] ;
			$vdw_r = $VDW_Radius_List[$element_name];
			
			foreach($Atoms_in_Res[$target_res_num] as $target_atom_number ){
				
				if (Isset($Proximity_Table[$atom_number][$target_atom_number])){continue;}	# avoid repeating calculation
				
				
				##tss calculation
				$tss = TSS_3D($Pos[$atom_number], $Pos[$target_atom_number]);	
				
				##cutoff calculation
				
				$target_element_name = $AN_To_Element_List[$target_atom_number] ;
				$target_vdw_r = $VDW_Radius_List[$target_element_name];
				
				//$cutoff_dist = $vdw_r + $probe_size * 2 + $target_vdw_r;	
				$cutoff_dist = $vdw_r + 2.8 + $target_vdw_r;
				$cutoff_tss = $cutoff_dist * $cutoff_dist;
				
				##cehck and save neighboring atoms
				if($tss <= $cutoff_tss){
					$Proximity_Table[$atom_number][$target_atom_number] = $tss; 
					$Proximity_Table[$target_atom_number][$atom_number] = $tss; ##pass info
					
				} 
		
					
			} //foreach($Atoms_in_Res[$target_res_num]
			
		}//foreach($Atoms_in_Res[$res_num] 
		
	}//foreach($Target_Res_List
	
}//foreach($Neighboring_Res_List



##Delete self-paring data 
$Temp = Array();
foreach($Proximity_Table as $atom_number => $Target_Data ){
	
	foreach($Target_Data as $target_atom_number => $tss){
		if($atom_number == $target_atom_number){continue;}
		$Temp[$atom_number][$target_atom_number] =  $tss;
	}//foreach($Target_Data as $target_atom_number => $tss){
	
	aSort($Temp[$atom_number]);

}//foreach($Proximity_Table as $atom_number => $Target_Data ){
$Proximity_Table = $Temp;
$Temp = Array();




##Define array storages in asa calculation process
$Log_Refined_ASA = Array();
$Resisue_ASA_Unnormalized_Sum  = Array(); 
foreach($ResNum_To_ResName as $res_num => $res_name){
	$Resisue_ASA_Unnormalized_Sum[$res_num] = 0;
}//foreach($ResNum_To_ResName as $res_num => $res_name){

$protein_asa_unnormalized = 0; #unnormalized protein_asa
$used_dot_count = 0;  #number of dots of used spherical grid groups   



##ASA calculation
for($group_id = 0; $group_id <  ( $asa_sample_size -1 ) ; $group_id++){

	$Expanded_Spherical_Grid = $Expanded_Spherical_Grid_Set["C"][$group_id]; 
	$dot_count = Count($Expanded_Spherical_Grid);
	
	
	foreach($Proximity_Table as $atom_number => $Target_AN_List){	
	
		List($x, $y, $z) = $Pos[$atom_number] ;
		
		$element_name = $AN_To_Element_List[$atom_number];
		
		
		$Expanded_Spherical_Grid = $Expanded_Spherical_Grid_Set[$element_name][$group_id]; 

		##Shifting spherical grid around target atom		
		
		$Shifted_Dot_List = Array();
		foreach($Expanded_Spherical_Grid as $Expanded_Dot){
			
			$Shifted_Dot_List[] = Array($x + $Expanded_Dot[0],$y + $Expanded_Dot[1], $z + $Expanded_Dot[2]);
			
		}//foreach($Expanded_Spherical_Grid as $Expanded_Dot){
			
			
		##Probing solvent available spots  	
		foreach($Target_AN_List as $target_atom_number => $tss){
		
			$target_element_name = $AN_To_Element_List[$target_atom_number];
			$Target_Pos = $Pos[$target_atom_number];
			
			$probe_tss_cutoff = $Probing_TSS_Cutoff_List[$target_element_name] ;
			
			foreach($Shifted_Dot_List as $index => $Shifted_Dot){
				if( TSS_3D($Target_Pos, $Shifted_Dot) < $probe_tss_cutoff ){
					Unset($Shifted_Dot_List[$index]);
				}
				
			}//foreach($Shifted_Dot as $Shifted_Dot){
		}//foreach($Target_AN_List as $target_atom_number => $tss){
		

		$lefted_dot_count = Count($Shifted_Dot_List);
		$atom_asa_unnormalized = $MAX_ATOM_ASA_List[$element_name] * $lefted_dot_count; #unnormalized atom asa
		$protein_asa_unnormalized += $atom_asa_unnormalized; #saving
		
		$res_num = $Atom_Info_DB[$atom_number][6];
		$Resisue_ASA_Unnormalized_Sum[$res_num] += $atom_asa_unnormalized; #saving
		
	
		
	}//foreach($Proximity_Table as $atom_number => $Target_AN_List){
	

	$used_dot_count += $dot_count;  #number of dots of used spherical grid groups   
	
	$refined_asa = $protein_asa_unnormalized / $used_dot_count; #asa 
	$Log_Refined_ASA[$group_id] = $refined_asa ; #saving asa for quality control 
	
}



for($group_id = $asa_sample_size -1 ; $group_id < $group_count; $group_id++){


	$Expanded_Spherical_Grid = $Expanded_Spherical_Grid_Set["C"][$group_id]; 
	$dot_count = Count($Expanded_Spherical_Grid);

	foreach($Proximity_Table as $atom_number => $Target_AN_List){	
	
		List($x, $y, $z) = $Pos[$atom_number] ;
		
		$element_name = $AN_To_Element_List[$atom_number];
		
		$Expanded_Spherical_Grid =  $Expanded_Spherical_Grid_Set[$element_name][$group_id]; 

		
		##Shifting spherical grid around target atom
		$Shifted_Dot_List = Array();
		foreach($Expanded_Spherical_Grid as $Expanded_Dot){
			$Shifted_Dot_List[] = Array($x + $Expanded_Dot[0],$y + $Expanded_Dot[1], $z + $Expanded_Dot[2]);
		}//foreach($Expanded_Spherical_Grid as $Expanded_Dot){
			
		##Probing solvent available spots  		
		foreach($Target_AN_List as $target_atom_number => $tss){
		
			$target_element_name = $AN_To_Element_List[$target_atom_number];
			$Target_Pos = $Pos[$target_atom_number];
			
			
			$probe_tss_cutoff = $Probing_TSS_Cutoff_List[$target_element_name] ;
			
			foreach($Shifted_Dot_List as $index => $Shifted_Dot){
				if( TSS_3D($Target_Pos, $Shifted_Dot) < $probe_tss_cutoff ){
					Unset($Shifted_Dot_List[$index]);
				}
				
			}//foreach($Shifted_Dot as $Shifted_Dot){
		}//foreach($Target_AN_List as $target_atom_number => $tss){
		

		$lefted_dot_count = Count($Shifted_Dot_List);
		$atom_asa_unnormalized = $MAX_ATOM_ASA_List[$element_name] * $lefted_dot_count;
		$protein_asa_unnormalized += $atom_asa_unnormalized; #saving
		
		
		$res_num = $Atom_Info_DB[$atom_number][6];
		$Resisue_ASA_Unnormalized_Sum[$res_num] += $atom_asa_unnormalized; #saving
		
	
	}//foreach($Proximity_Table as $atom_number => $Target_AN_List){

	
	$used_dot_count += $dot_count; #number of dots of used spherical grid groups   
	$refined_asa = $protein_asa_unnormalized / $used_dot_count; #asa 
	
	$Log_Refined_ASA[$group_id] = $refined_asa ; #saving asa for quality control 
	
	//$real_shell_num = $group_id + 1;
		
	##take latest asa results for quality check 
	
	$Array_ASA_Samples = Array();
	for($sampe_group_id = $group_id - $asa_sample_size + 1 ; $sampe_group_id < $group_id + 1; $sampe_group_id++ ){
		$Array_ASA_Samples[] = $Log_Refined_ASA[$sampe_group_id] ;
	}
	
	$sd_asa_sample = Stand_Deviation($Array_ASA_Samples);
	$cv = $sd_asa_sample / $refined_asa;
	if($cv < $qcv_cutoff){break;}  #pass quality checking
	

}


## Print out asa
$refined_asa = Round($refined_asa, 3);
echo "ASA: $refined_asa (squared angstrom)\n";

## Get RSA and print out
echo "RSA:\n";
echo "Residue Number\tRSA\n";


foreach($Resisue_ASA_Unnormalized_Sum as $res_num => $residue_asa_unnormalized){
	
	$res_name = $ResNum_To_ResName[$res_num];
	if(Isset($Max_RES_ASA[$res_name])){
		
		$refined_residue_asa =  $residue_asa_unnormalized / $used_dot_count;
		
		$refined_residue_rsa = $refined_residue_asa / $Max_RES_ASA[$res_name];
		$refined_residue_rsa = Round($refined_residue_rsa, 3);
		echo "$res_num\t$refined_residue_rsa\n";
		
	}else{
		echo "$res_num\tError\n";
	}

	
}//foreach($Resisue_ASA_Unnormalized_Sum as $res_num => $residue_asa_unnormalized){





/**/


?>