for(group_id := asa_sample_size -1 ; group_id < group_count; group_id++){

	Expanded_Spherical_Grid = Expanded_Spherical_Grid_Set["C"][group_id] 
	dot_count = len(Expanded_Spherical_Grid)


	for atom_number, Target_AN_List := range Proximity_Table{
		x := Pos[atom_number][0] 
    y := Pos[atom_number][1]
    z := Pos[atom_number][2]
		
		element_name = AN_To_Element_List[atom_number];
		
		Expanded_Spherical_Grid =  Expanded_Spherical_Grid_Set[element_name][group_id]; 

		
		//##Shifting spherical grid around target atom
		Shifted_Dot_List = map[int](map[int]float64){}
    for num, Expanded_Dot := range Expanded_Spherical_Grid{
      Shifted_Dot_List[num] = map[int]float64(0:x + Expanded_Dot[0], 1:y + Expanded_Dot[1], 2:z + $Expanded_Dot[2];
		} //for num, Expanded_Dot := range Expanded_Spherical_Grid
			
		//Probing solvent available spots  		
		for target_atom_number, tss := range Target_AN_List{
			target_element_name = AN_To_Element_List[target_atom_number];
			Target_Pos = Pos[target_atom_number];
			
			
			probe_tss_cutoff = Probing_TSS_Cutoff_List[target_element_name] ;
			
      for index, Shifted_Dot := range Shifted_Dot_List{
				if( TSS_3D(Target_Pos, Shifted_Dot) < probe_tss_cutoff ){
					Shifted_Dot_List = append(Shifted_Dot_List[:index], Shifted_Dot_List[index+1:])
				}
				
			}
		}for target_atom_number, tss := range Target_AN_List{
		
		lefted_dot_count := len(Shifted_Dot_List);
		atom_asa_unnormalized := MAX_ATOM_ASA_List[element_name] * lefted_dot_count;
		protein_asa_unnormalized += atom_asa_unnormalized; //saving
		
		
		res_num := Atom_Info_DB[atom_number][6];
		Resisue_ASA_Unnormalized_Sum[res_num] += atom_asa_unnormalized; //saving
		
	
	}//foreach($Proximity_Table as $atom_number => $Target_AN_List){

	
	used_dot_count += dot_count; //number of dots of used spherical grid groups   
	refined_asa = protein_asa_unnormalized / used_dot_count; //asa 
	
	Log_Refined_ASA[group_id] = refined_asa ; //saving asa for quality control 
	
	//$real_shell_num = $group_id + 1;
		
	//take latest asa results for quality check 
	
	***Array_ASA_Samples = map[int]float64{}
  count := 0
	for(sampe_group_id := group_id - asa_sample_size + 1 ; sampe_group_id < group_id + 1; sampe_group_id++ ){
		Array_ASA_Samples[count] = Log_Refined_ASA[sampe_group_id] 
    count++
	}
	
	sd_asa_sample := Stand_Deviation(Array_ASA_Samples);
	cv := sd_asa_sample / refined_asa;
	if(cv < qcv_cutoff){break}  //pass quality checking
	
}

// Print out asa
refined_asa := lo.Round(refined_asa, 3);
lo.Println( "ASA", refined_asa, "(squared angstrom)\n" )

// Get RSA and print out
lo.Println( "RSA:\n" )
lo.Println( "Residue Number\tRSA\n" )

for res_num, residue_asa_unnormalized := range Resisue_ASA_Unnormalized_Sum{

	res_name = ResNum_To_ResName[res_num]
  
  var check bool = false
  for keys, values := range  Max_RES_ASA{
      _=values
      if keys == res_name{
        check = true
        break
      }
  }
  if(!check){
    
		refined_residue_asa =  residue_asa_unnormalized / used_dot_count;
		
		refined_residue_rsa = refined_residue_asa / Max_RES_ASA[res_name];
		refined_residue_rsa = lo.Round(refined_residue_rsa, 3)
		lo.Println( lo.String(res_num) + "\t" + lo.String(refined_residue_rsa) + "\n")
    
	}else{
    lo.Println( lo.String(res_num) + "\tError\n" )
	}

} //for res_num, residue_asa_unnormalized := range Resisue_ASA_Unnormalized_Sum{
